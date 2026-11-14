// Plotfile round trip (cell-centred variables and nodal z_phys_nd).
// Checkpoint/restart round trip, onto a different box layout.
#include "FileMagic.H"
#include "StateIO.H"
#include "State.H"
#include "ZLevels.H"
#include "TestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <gtest/gtest.h>

#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

using gust::test::count_bit_differences;
using gust::test::fill_pattern;

gust::MeshInputs test_mesh ()
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(20, 12, 10);
    // -100/3 is not a short decimal: a header that printed it with fewer
    // digits than a hex float would fail the restart's exact mesh check.
    m.prob_lo = {amrex::Real(-100.0 / 3.0), 50, 0};
    m.prob_hi = {300, 290, 700};
    m.is_periodic = {1, 0, 0};
    m.max_grid_size = {8, 8};   // 3 x 2 boxes, uneven in x
    return m;
}

// Ghost values that disagree with every valid value, so a restart that let a
// checkpoint ghost land on a valid point would show up.
void fill_state (gust::State& s)
{
    int c = 0;
    for (auto& f : s.fields()) {
        fill_pattern(*f.mf, ++c);
        f.mf->setBndry(amrex::Real(1.0e30) * c);
    }
    s.step = 7;
}

// Removes a directory written by a test, on the I/O rank.
struct ScopedDir
{
    std::string name;
    ~ScopedDir () {
        amrex::ParallelDescriptor::Barrier();
        if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(name); }
        amrex::ParallelDescriptor::Barrier();
    }
};

// Cell-centred values the plotfile should hold, computed independently of
// write_plotfile on the host from the same fields.
amrex::MultiFab expected_plot (gust::State& s)
{
    const auto& ba = s.p.boxArray();
    amrex::MultiFab out(ba, s.p.DistributionMap(), 7, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    for (amrex::MFIter mfi(out, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        auto copy_host = [&] (const amrex::MultiFab& mf) {
            amrex::FArrayBox h(mf[mfi].box(), 1, amrex::The_Pinned_Arena());
            h.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
            amrex::Gpu::streamSynchronize();
            return h;
        };
        const auto u = copy_host(s.u), v = copy_host(s.v), w = copy_host(s.w);
        const auto p = copy_host(s.p), th = copy_host(s.theta), k = copy_host(s.tke);
        const auto z = copy_host(s.z_phys_nd);
        auto& o = out[mfi];
        const amrex::Real h = 0.5;
        const amrex::Real e = 0.125;
        amrex::LoopOnCpu(bx, [&] (int i, int j, int kk) {
            using IV = amrex::IntVect;
            o(IV(i, j, kk), 0) = h * (u(IV(i, j, kk)) + u(IV(i + 1, j, kk)));
            o(IV(i, j, kk), 1) = h * (v(IV(i, j, kk)) + v(IV(i, j + 1, kk)));
            o(IV(i, j, kk), 2) = h * (w(IV(i, j, kk)) + w(IV(i, j, kk + 1)));
            o(IV(i, j, kk), 3) = p(IV(i, j, kk));
            o(IV(i, j, kk), 4) = th(IV(i, j, kk));
            o(IV(i, j, kk), 5) = k(IV(i, j, kk));
            o(IV(i, j, kk), 6) = e * (z(IV(i, j, kk)) + z(IV(i + 1, j, kk))
                                    + z(IV(i, j + 1, kk)) + z(IV(i + 1, j + 1, kk))
                                    + z(IV(i, j, kk + 1)) + z(IV(i + 1, j, kk + 1))
                                    + z(IV(i, j + 1, kk + 1)) + z(IV(i + 1, j + 1, kk + 1)));
        });
    }
    return out;
}

TEST (IO, PlotfileRoundTripIsBitIdentical)
{
    const auto mesh = test_mesh();
    const auto geom = gust::make_geometry(mesh);
    const auto ba = gust::make_box_array(mesh);
    ASSERT_EQ(ba.size(), 6);
    const amrex::DistributionMapping dm(ba);

    gust::State s;
    s.define(ba, dm);
    fill_state(s);

    const ScopedDir dir{gust::write_plotfile("ut_io_plt", s, geom)};
    EXPECT_EQ(dir.name, "ut_io_plt00007");

    amrex::PlotFileData pf(dir.name);
    ASSERT_EQ(pf.spaceDim(), 3);
    EXPECT_EQ(pf.levelStep(0), 7);
    EXPECT_EQ(pf.varNames(), std::vector<std::string>(gust::plot_var_names().begin(),
                                                      gust::plot_var_names().end()));
    EXPECT_EQ(pf.probDomain(0), geom.Domain());

    const amrex::MultiFab expect = expected_plot(s);
    const auto& names = gust::plot_var_names();
    for (int n = 0; n < static_cast<int>(names.size()); ++n) {
        const amrex::MultiFab& got = pf.get(0, names[n]);
        amrex::MultiFab a(ba, dm, 1, 0), b(ba, dm, 1, 0);
        a.ParallelCopy(got, 0, 0, 1);
        b.ParallelCopy(expect, n, 0, 1);
        EXPECT_EQ(count_bit_differences(a, b, 0), 0) << names[n];
    }

    // Nodal heights, ghost nodes included. Read into the same BoxArray and
    // DistributionMapping, so each box is compared with its own ghosts.
    amrex::MultiFab z_read(s.z_phys_nd.boxArray(), dm, 1, gust::ngrow_nodal);
    amrex::VisMF::Read(z_read, dir.name + "/Level_0/z_phys_nd");
    EXPECT_EQ(count_bit_differences(z_read, s.z_phys_nd, gust::ngrow_nodal), 0);
}

TEST (IO, CheckpointRestartIsBitIdenticalOnANewLayout)
{
    const auto mesh = test_mesh();
    const auto ba = gust::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);

    gust::State s;
    s.define(ba, dm);
    fill_state(s);
    const ScopedDir dir{gust::write_checkpoint("ut_io_chk", s, mesh)};
    EXPECT_EQ(dir.name, "ut_io_chk00007");

    // Restart onto another layout: 2 x 3 boxes instead of 3 x 2.
    auto mesh2 = mesh;
    mesh2.max_grid_size = {10, 4};
    const auto ba2 = gust::make_box_array(mesh2);
    ASSERT_EQ(ba2.size(), 6);
    ASSERT_NE(ba2, ba);
    const amrex::DistributionMapping dm2(ba2);
    gust::State r;
    r.define(ba2, dm2);
    gust::read_checkpoint(dir.name, r, mesh);
    EXPECT_EQ(r.step, 7);

    // Compare on the original layout, valid region.
    auto rf = r.fields();
    auto sf = s.fields();
    for (int n = 0; n < static_cast<int>(sf.size()); ++n) {
        const amrex::MultiFab& orig = *sf[n].mf;
        amrex::MultiFab back(orig.boxArray(), dm, 1, 0);
        back.ParallelCopy(*rf[n].mf, 0, 0, 1);
        amrex::MultiFab ref(orig.boxArray(), dm, 1, 0);
        amrex::MultiFab::Copy(ref, orig, 0, 0, 1, 0);
        EXPECT_EQ(count_bit_differences(back, ref, 0), 0) << sf[n].name;
    }
}

// The checkpoint's solution fields hold every ghost image as the
// neighbouring box's (or the periodic image's) valid value, even when the
// state's own ghosts are stale (fill_state sets them to 1e30); the state is
// not touched (audit A5, from Pressel's A3 review: theta and k images were
// written stale). On 3 x 2 uneven boxes, x periodic.
TEST (IO, CheckpointGhostImagesAreTheNeighboursValues)
{
    const auto mesh = test_mesh();
    const auto geom = gust::make_geometry(mesh);
    const auto ba = gust::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);
    gust::State s;
    s.define(ba, dm);
    fill_state(s);
    // One value per shared face, as the solver keeps them (the pattern gives
    // a periodic face and its image different values, and which copy an
    // image takes would then depend on the ranks).
    int c = 0;
    for (auto& f : s.fields()) {
        ++c;
        f.mf->OverrideSync(geom.periodicity());
        f.mf->setBndry(amrex::Real(1.0e30) * c);
    }
    const ScopedDir dir{gust::write_checkpoint("ut_io_chkghost", s, mesh)};
    c = 0;
    for (auto& f : s.fields()) {
        ++c;
        if (f.mf == &s.z_phys_nd) { continue; }
        amrex::MultiFab src;
        amrex::VisMF::Read(src, amrex::MultiFabFileFullPrefix(0, dir.name, "Level_", f.name));
        ASSERT_EQ(src.nGrowVect(), f.mf->nGrowVect()) << f.name;
        // The images from the file's own valid data: equal to what it holds.
        amrex::MultiFab ref(src.boxArray(), src.DistributionMap(), 1, src.nGrowVect());
        amrex::MultiFab::Copy(ref, src, 0, 0, 1, src.nGrowVect());
        ref.FillBoundary(geom.periodicity());
        EXPECT_EQ(count_bit_differences(src, ref, src.nGrow()), 0) << f.name;
        // The state keeps its own ghosts.
        EXPECT_EQ(f.mf->max(0, f.mf->nGrow()), amrex::Real(1.0e30) * c) << f.name;
    }
}

// A checkpoint written before the solution fields had two ghost cells (one,
// before the MUSCL convection operators) restarts bit for bit: only valid
// data is restored.
TEST (IO, CheckpointWithOneGhostCellRestarts)
{
    const auto mesh = test_mesh();
    const auto ba = gust::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);

    gust::State s;
    s.define(ba, dm);
    for (auto& f : s.fields()) {
        if (f.mf == &s.z_phys_nd) { continue; }
        const amrex::BoxArray fba = f.mf->boxArray();
        f.mf->clear();
        f.mf->define(fba, dm, 1, 1);
    }
    fill_state(s);
    const ScopedDir dir{gust::write_checkpoint("ut_io_chk1g", s, mesh)};

    gust::State r;
    r.define(ba, dm);
    ASSERT_EQ(r.u.nGrow(), 2);
    gust::read_checkpoint(dir.name, r, mesh);
    auto rf = r.fields();
    auto sf = s.fields();
    for (int n = 0; n < static_cast<int>(sf.size()); ++n) {
        EXPECT_EQ(count_bit_differences(*rf[n].mf, *sf[n].mf, 0), 0) << sf[n].name;
    }
}

// A checkpoint written before the rename to GUST starts its Header with
// "ablsteady-checkpoint" and restarts bit for bit.
TEST (IO, CheckpointWithTheOldMagicRestarts)
{
    const auto mesh = test_mesh();
    const auto ba = gust::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);

    gust::State s;
    s.define(ba, dm);
    fill_state(s);
    const ScopedDir dir{gust::write_checkpoint("ut_io_chkold", s, mesh)};
    // The I/O rank rewrites the Header; every rank learns whether it could,
    // so a failed check stops all ranks before the collective read.
    bool rewritten = true;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream in(dir.name + "/Header");
        std::stringstream ss;
        ss << in.rdbuf();
        std::string text = ss.str();
        const std::string now = "gust-checkpoint ";
        rewritten = text.compare(0, now.size(), now) == 0;
        if (rewritten) {
            text.replace(0, now.size(), "ablsteady-checkpoint ");
            std::ofstream(dir.name + "/Header") << text;
        }
    }
    amrex::ParallelDescriptor::ReduceBoolAnd(rewritten);
    ASSERT_TRUE(rewritten) << "the Header does not start with gust-checkpoint";

    gust::State r;
    r.define(ba, dm);
    gust::read_checkpoint(dir.name, r, mesh);
    EXPECT_EQ(r.step, 7);
    auto rf = r.fields();
    auto sf = s.fields();
    for (int n = 0; n < static_cast<int>(sf.size()); ++n) {
        EXPECT_EQ(count_bit_differences(*rf[n].mf, *sf[n].mf, 0), 0) << sf[n].name;
    }
}

// Format version 2 keeps the reference z-levels in the Header, bit for bit
// (hex floats), for a nest's check of its parent (NestParent.H); a count
// that does not match mesh.n_cell z aborts naming the key the caller gives.
TEST (IO, CheckpointHeaderKeepsTheZLevels)
{
    auto mesh = test_mesh();
    mesh.dz0 = amrex::Real(100.0 / 7.0);   // stretched: levels that are not short decimals
    const auto geom = gust::make_geometry(mesh);
    const auto zl = gust::make_zlevels(mesh, geom);
    const auto ba = gust::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);
    gust::State s;
    s.define(ba, dm);
    fill_state(s);
    const ScopedDir dir{gust::write_checkpoint("ut_io_chkzl", s, mesh, &zl)};
    const gust::CheckpointHeader h = gust::read_checkpoint_header(dir.name, "io.restart");
    EXPECT_EQ(h.version, gust::kCheckpointFormatVersion);
    EXPECT_EQ(h.version, 2);
    EXPECT_EQ(h.step, 7);
    EXPECT_EQ(h.n_cell, mesh.n_cell);
    ASSERT_EQ(h.zlevels.size(), zl.size());
    for (std::size_t k = 0; k < zl.size(); ++k) {
        EXPECT_EQ(std::memcmp(&h.zlevels[k], &zl[k], sizeof(amrex::Real)), 0) << "z-level " << k;
    }
    EXPECT_EQ(h.ba, ba);

    // One z-level fewer in the Header.
    bool rewritten = false;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream in(dir.name + "/Header");
        std::ostringstream os;
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("zlevels ", 0) == 0) {
                std::istringstream ls(line);
                std::string key;
                int n = 0;
                ls >> key >> n;
                std::ostringstream nl;
                nl << "zlevels " << n - 1;
                std::string z;
                for (int k = 0; k < n - 1 && (ls >> z); ++k) { nl << " " << z; }
                line = nl.str();
                rewritten = true;
            }
            os << line << "\n";
        }
        in.close();
        std::ofstream(dir.name + "/Header") << os.str();
    }
    amrex::ParallelDescriptor::ReduceBoolOr(rewritten);
    ASSERT_TRUE(rewritten);
    std::string msg;
    try { (void)gust::read_checkpoint_header(dir.name, "bc.nest_parent"); } catch (const std::runtime_error& e) { msg = e.what(); }
    EXPECT_NE(msg.find("bc.nest_parent"), std::string::npos) << msg;

    // No zlevels line at all: the message names the caller's key, not the
    // token found in its place.
    bool dropped = false;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream in(dir.name + "/Header");
        std::ostringstream os;
        std::string line;
        while (std::getline(in, line)) {
            if (line.rfind("zlevels ", 0) == 0) { dropped = true; continue; }
            os << line << "\n";
        }
        in.close();
        std::ofstream(dir.name + "/Header") << os.str();
    }
    amrex::ParallelDescriptor::ReduceBoolOr(dropped);
    ASSERT_TRUE(dropped);
    msg.clear();
    try { (void)gust::read_checkpoint_header(dir.name, "bc.nest_parent"); } catch (const std::runtime_error& e) { msg = e.what(); }
    EXPECT_NE(msg.find("gust: bc.nest_parent = "), std::string::npos) << msg;
    EXPECT_NE(msg.find("expected 'zlevels'"), std::string::npos) << msg;
}

// The readers take "gust-<kind>" and the older "ablsteady-<kind>" of the
// same kind, and nothing else.
TEST (IO, FileMagicTakesTheNewAndTheOldName)
{
    EXPECT_TRUE(gust::is_file_magic("gust-checkpoint", "gust-checkpoint"));
    EXPECT_TRUE(gust::is_file_magic("ablsteady-checkpoint", "gust-checkpoint"));
    EXPECT_TRUE(gust::is_file_magic("ablsteady-newton-stage", "gust-newton-stage"));
    EXPECT_FALSE(gust::is_file_magic("ablsteady-driver", "gust-checkpoint"));
    EXPECT_FALSE(gust::is_file_magic("gust-driver", "gust-checkpoint"));
    EXPECT_FALSE(gust::is_file_magic("ablsteady-", "gust-checkpoint"));
    EXPECT_FALSE(gust::is_file_magic("checkpoint", "gust-checkpoint"));
    EXPECT_FALSE(gust::is_file_magic("", "gust-checkpoint"));
    EXPECT_FALSE(gust::is_file_magic("GUST-checkpoint", "gust-checkpoint"));
}

TEST (IO, CheckpointRestartRestoresNodalGhostsOnTheSameLayout)
{
    const auto mesh = test_mesh();
    const auto ba = gust::make_box_array(mesh);
    const amrex::DistributionMapping dm(ba);

    gust::State s;
    s.define(ba, dm);
    gust::init_flat_z_phys_nd(s.z_phys_nd, gust::make_geometry(mesh));
    const ScopedDir dir{gust::write_checkpoint("ut_io_chkz", s, mesh)};

    gust::State r;
    r.define(ba, dm);
    gust::read_checkpoint(dir.name, r, mesh);
    EXPECT_EQ(count_bit_differences(r.z_phys_nd, s.z_phys_nd, gust::ngrow_nodal), 0);
}

} // namespace

// The checkpoint-overwrite guard asks the file system whether two paths name
// one directory (audit A10, P3): "./", trailing and doubled separators, an
// absolute path and a symbolic link all name the restart directory; another
// directory, or one that does not exist, does not. Before, the guard compared
// the strings after stripping trailing '/' only.
TEST (IO, SameDirectoryAsksTheFileSystem)
{
    namespace fs = std::filesystem;
    const std::string base = "ut_same_dir";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::remove_all(base);
        fs::create_directories(base + "/chk00040");
        fs::create_directories(base + "/chk00041");
        std::error_code ec;
        fs::create_directory_symlink("chk00040", base + "/link00040", ec);
    }
    amrex::ParallelDescriptor::Barrier();
    const std::string chk = base + "/chk00040";
    EXPECT_TRUE(gust::same_directory(chk, chk));
    EXPECT_TRUE(gust::same_directory(chk, "./" + chk));
    EXPECT_TRUE(gust::same_directory(chk, chk + "/"));
    EXPECT_TRUE(gust::same_directory(chk, base + "//chk00040//"));
    EXPECT_TRUE(gust::same_directory(chk, fs::absolute(chk).string()));
    EXPECT_TRUE(gust::same_directory(chk, base + "/../" + base + "/chk00040"));
    if (fs::is_symlink(base + "/link00040")) { EXPECT_TRUE(gust::same_directory(chk, base + "/link00040")); }
    EXPECT_FALSE(gust::same_directory(chk, base + "/chk00041"));
    EXPECT_FALSE(gust::same_directory(base + "/chk00042", chk));   // the target does not exist yet
    EXPECT_FALSE(gust::same_directory(chk, base + "/missing"));
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { fs::remove_all(base); }
}
