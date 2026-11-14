// The parent run of a nested patch (NestParent.H), plan rev. 6 sub-phase 9.1.
//   UT-9.1.1  a field a + b x + c y + d x y on each zeta-level (coefficients
//             that change with the level and the field) is reproduced to
//             round-off at every staggered point of the patch, its two
//             horizontal ghost layers included, on flat ground and on a hill
//             with stretched layers; the z ghost rows hold the nearest row
//   UT-9.1.2  each placement fault is refused with a message naming
//             bc.nest_parent; a checkpoint of format version 1 and a missing
//             checkpoint abort naming it
//   UT-9.1.3  the patch values do not depend on the parent's or the patch's
//             box layout or rank map (bit for bit), on 1, 2 and 4 ranks (the
//             unit_np* tests)
//   beyond the plan: across a periodic seam of the parent (a patch that
//             touches the seam, and a periodic patch spanning the parent) the
//             values are the bilinear interpolation of the parent's discrete
//             values with the periodic neighbours, computed here on the host
// The parent's ghost cells hold 1e30 in the checkpoint and the patch starts
// at 1e30, so a value taken from a parent ghost or a point left unfilled
// shows up.
#include "Inputs.H"
#include "Mesh.H"
#include "NestParent.H"
#include "State.H"
#include "StateIO.H"
#include "ZLevels.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double kTol = kDouble ? 1.0e-11 : 2.0e-5;
constexpr double kBig = 1.0e30;

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

// A run's level 0 as App.cpp lays it out: geometry, layout, mesh (node
// heights built from the terrain inputs) and state.
struct NestRun
{
    gust::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    gust::Mesh mesh;
    gust::State q;

    //! reverse: the boxes' ranks in reverse order (another rank map).
    void make (bool reverse = false)
    {
        geom = gust::make_geometry(in.mesh);
        ba = gust::make_box_array(in.mesh);
        dm = amrex::DistributionMapping(ba);
        if (reverse) {
            amrex::Vector<int> pmap = dm.ProcessorMap();
            const int np = amrex::ParallelDescriptor::NProcs();
            for (auto& p : pmap) { p = np - 1 - p; }
            dm = amrex::DistributionMapping(pmap);
        }
        q.define(ba, dm);
        mesh.define(ba, dm);
        mesh.build(in, geom, q.z_phys_nd);
    }
};

// The parent: 20 x 16 x 8 cells of 100 m, stretched layers (dz0 = 40), boxes
// of at most 7 x 6 (uneven). hill: a round Agnesi hill.
gust::Inputs parent_inputs (bool hill, int mgs = 7)
{
    gust::Inputs in;
    in.mesh.n_cell = amrex::IntVect(20, 16, 8);
    in.mesh.prob_lo = {-1000, 500, 0};
    in.mesh.prob_hi = {1000, 2100, 800};
    in.mesh.is_periodic = {0, 0, 0};
    in.mesh.max_grid_size = {mgs, mgs - 1};
    in.mesh.dz0 = 40;
    if (hill) {
        in.terrain.type = "agnesi";
        in.terrain.hill_shape = "round";
        in.terrain.hill_height = 150;
        in.terrain.hill_half_width = 300;
        in.terrain.smoothing = 1;
    }
    return in;
}

// The patch over parent cells 3 to 15 in x and 2 to 12 in y: 24 x 20 cells
// of 50 m, the parent's z-levels, boxes of at most 9 x 7 (uneven). Its own
// hill (centred in the patch) when hill: the node heights are the patch's.
gust::Inputs patch_inputs (const gust::Inputs& parent, int mgs = 9)
{
    gust::Inputs in = parent;
    in.mesh.n_cell = amrex::IntVect(24, 20, 8);
    in.mesh.prob_lo = {-700, 700, 0};
    in.mesh.prob_hi = {500, 1700, 800};
    in.mesh.max_grid_size = {mgs, mgs - 2};
    return in;
}

// The closed form: bilinear in x and y on each level k, other coefficients
// for each field c.
double bilinear (double x, double y, int k, int c)
{
    return (1.0 + c) + 0.3 * k + (0.002 + 1.0e-4 * c) * x - (0.0015 - 1.0e-4 * k) * y
         + 1.0e-6 * (1.0 + 0.1 * k + 0.05 * c) * x * y;
}

// The horizontal position of point (i, j) of mf's type on geometry g.
double xpos (const amrex::Geometry& g, const amrex::MultiFab& mf, int i, int d)
{
    const double off = mf.ixType().nodeCentered(d) ? 0.0 : 0.5;
    return static_cast<double>(g.ProbLo(d)) + (i + off) * static_cast<double>(g.CellSize(d));
}

// Valid points of mf = fn(i, j, k); ghost points kBig.
void fill_valid (amrex::MultiFab& mf, const std::function<double(int, int, int)>& fn)
{
    mf.setVal(amrex::Real(kBig));
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(i, j, k)); });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

// The six fields of a state with their field numbers.
std::vector<std::pair<amrex::MultiFab*, int>> fields6 (gust::State& q)
{
    return {{&q.u, 0}, {&q.v, 1}, {&q.w, 2}, {&q.p, 3}, {&q.theta, 4}, {&q.tke, 5}};
}

// Writes the parent's checkpoint with field values fn(field, mf, i, j, k).
std::string write_parent (NestRun& p, const std::string& prefix,
                          const std::function<double(const amrex::MultiFab&, int, int, int, int)>& fn)
{
    for (auto [mf, c] : fields6(p.q)) {
        const amrex::MultiFab* m = mf;
        const int cc = c;
        fill_valid(*mf, [&] (int i, int j, int k) { return fn(*m, cc, i, j, k); });
    }
    return gust::write_checkpoint(prefix, p.q, p.in.mesh, &p.mesh.zlevels());
}

// Largest |patch - expected| over the valid points grown by ngrow_state in x
// and y (all valid levels), and the number of z ghost points that differ from
// the nearest valid row. expected(mf, c, i, j, k) by patch index.
struct Check { double err{0}; long zghost_bad{0}; long points{0}; };
Check check_patch (gust::State& q, const std::function<double(const amrex::MultiFab&, int, int, int, int)>& expected)
{
    Check r;
    for (auto [mf, c] : fields6(q)) {
        const amrex::IntVect ngh(gust::ngrow_state, gust::ngrow_state, 0);
        const int ngz = mf->nGrowVect()[2];
        for (amrex::MFIter mfi(*mf, false); mfi.isValid(); ++mfi) {
            const amrex::Box vb = amrex::grow(mfi.validbox(), ngh);
            amrex::FArrayBox h((*mf)[mfi].box(), 1, amrex::The_Pinned_Arena());
            h.copy<amrex::RunOn::Device>((*mf)[mfi], (*mf)[mfi].box());
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(vb, [&] (int i, int j, int k) {
                const double d = std::abs(static_cast<double>(h(amrex::IntVect(i, j, k))) - expected(*mf, c, i, j, k));
                r.err = std::max(r.err, std::isfinite(d) ? d : kBig);
                ++r.points;
            });
            const int klo = vb.smallEnd(2), khi = vb.bigEnd(2);
            for (int g = 1; g <= ngz; ++g) {
                amrex::LoopOnCpu(amrex::Box(vb).setSmall(2, klo).setBig(2, klo), [&] (int i, int j, int) {
                    if (h(amrex::IntVect(i, j, klo - g)) != h(amrex::IntVect(i, j, klo))) { ++r.zghost_bad; }
                    if (h(amrex::IntVect(i, j, khi + g)) != h(amrex::IntVect(i, j, khi))) { ++r.zghost_bad; }
                });
            }
        }
    }
    const auto comm = amrex::ParallelDescriptor::Communicator();
    amrex::ParallelAllReduce::Max(r.err, comm);
    amrex::ParallelAllReduce::Sum(r.zghost_bad, comm);
    amrex::ParallelAllReduce::Sum(r.points, comm);
    return r;
}

// The parent read for the patch, then interpolated into it (starting at kBig).
void nest (const std::string& dir, NestRun& patch)
{
    for (auto [mf, c] : fields6(patch.q)) { mf->setVal(amrex::Real(kBig)); }
    const gust::NestParent parent(dir, patch.in.mesh, patch.mesh.zlevels(), patch.ba, patch.dm);
    parent.interpolate(patch.q);
}

// UT-9.1.1 on the given ground.
void bilinear_is_reproduced (bool hill, const std::string& prefix)
{
    NestRun p;
    p.in = parent_inputs(hill);
    p.make();
    const ScopedDir dir{write_parent(p, prefix, [&] (const amrex::MultiFab& mf, int c, int i, int j, int k) {
        return bilinear(xpos(p.geom, mf, i, 0), xpos(p.geom, mf, j, 1), k, c);
    })};
    NestRun f;
    f.in = patch_inputs(p.in);
    f.make();
    nest(dir.name, f);
    const Check r = check_patch(f.q, [&] (const amrex::MultiFab& mf, int c, int i, int j, int k) {
        return bilinear(xpos(f.geom, mf, i, 0), xpos(f.geom, mf, j, 1), k, c);
    });
    EXPECT_GT(r.points, 0);
    EXPECT_LT(r.err, kTol) << (hill ? "hill" : "flat");
    EXPECT_EQ(r.zghost_bad, 0) << (hill ? "hill" : "flat");
}

TEST (NestParent, BilinearFieldIsReproducedOnFlatGround) { bilinear_is_reproduced(false, "ut_nest_flat"); }
TEST (NestParent, BilinearFieldIsReproducedOnAHill) { bilinear_is_reproduced(true, "ut_nest_hill"); }

// The Header of the parent as check_nest_placement sees it.
gust::CheckpointHeader header_of (const gust::Inputs& in, int version = 2)
{
    gust::CheckpointHeader h;
    h.version = version;
    h.n_cell = in.mesh.n_cell;
    h.prob_lo = {in.mesh.prob_lo[0], in.mesh.prob_lo[1], in.mesh.prob_lo[2]};
    h.prob_hi = {in.mesh.prob_hi[0], in.mesh.prob_hi[1], in.mesh.prob_hi[2]};
    h.is_periodic = in.mesh.is_periodic;
    if (version >= 2) { h.zlevels = gust::make_zlevels(in.mesh, gust::make_geometry(in.mesh)); }
    return h;
}

std::vector<std::string> placement (const gust::CheckpointHeader& h, const gust::MeshInputs& m, amrex::IntVect* off = nullptr)
{
    return gust::check_nest_placement(h, m, gust::make_zlevels(m, gust::make_geometry(m)), off);
}

// UT-9.1.2: one fault at a time; each is refused, every message names the key.
TEST (NestParent, EachPlacementFaultIsRefusedNamingTheKey)
{
    const gust::Inputs pin = parent_inputs(false);
    const gust::CheckpointHeader h = header_of(pin);
    const gust::MeshInputs good = patch_inputs(pin).mesh;
    amrex::IntVect off(-1);
    EXPECT_TRUE(placement(h, good, &off).empty());
    EXPECT_EQ(off, amrex::IntVect(3, 2, 0));

    struct Fault { const char* what; std::function<void(gust::MeshInputs&)> make; gust::CheckpointHeader hdr; };
    gust::CheckpointHeader h_per = h;
    h_per.is_periodic = {0, 0, 1};
    const std::vector<Fault> faults = {
        {"outside", [] (gust::MeshInputs& m) { m.prob_lo[0] = -1100; m.n_cell[0] = 32; }, h},
        {"inside the margin (lo)", [] (gust::MeshInputs& m) { m.prob_lo[0] = -900; m.n_cell[0] = 28; }, h},
        {"inside the margin (hi)", [] (gust::MeshInputs& m) { m.prob_hi[1] = 2000; m.n_cell[1] = 26; }, h},
        {"misaligned", [] (gust::MeshInputs& m) { m.prob_lo[0] += 30; m.prob_hi[0] += 30; }, h},
        {"misaligned by a quarter cell", [] (gust::MeshInputs& m) { m.prob_lo[1] += 25; m.prob_hi[1] += 25; }, h},
        {"ratio 4", [] (gust::MeshInputs& m) { m.n_cell[0] = 48; }, h},
        {"ratio 1", [] (gust::MeshInputs& m) { m.n_cell[1] = 10; }, h},
        {"other n_cell z", [] (gust::MeshInputs& m) { m.n_cell[2] = 10; }, h},
        {"other prob_hi z", [] (gust::MeshInputs& m) { m.prob_hi[2] = 900; }, h},
        {"other z-levels", [] (gust::MeshInputs& m) { m.dz0 = 30; }, h},
        {"periodic patch, parent not", [] (gust::MeshInputs& m) { m.is_periodic[0] = 1; }, h},
        {"periodic z", [] (gust::MeshInputs&) {}, h_per},
        {"version 1", [] (gust::MeshInputs&) {}, header_of(pin, 1)},
    };
    for (const auto& f : faults) {
        gust::MeshInputs m = good;
        f.make(m);
        const auto err = placement(f.hdr, m);
        EXPECT_FALSE(err.empty()) << f.what;
        for (const auto& e : err) { EXPECT_EQ(e.rfind("bc.nest_parent: ", 0), 0U) << f.what << ": " << e; }
    }

    // A periodic parent: a patch touching its x faces is accepted (the
    // periodic images are its neighbours); a periodic patch must span it.
    gust::CheckpointHeader hp = header_of(pin);
    hp.is_periodic = {1, 0, 0};
    gust::MeshInputs edge = good;
    edge.prob_lo[0] = -1000;
    edge.prob_hi[0] = 0;
    edge.n_cell[0] = 20;
    EXPECT_TRUE(placement(hp, edge, &off).empty());
    EXPECT_EQ(off[0], 0);
    gust::MeshInputs span = edge;
    span.prob_hi[0] = 1000;
    span.n_cell[0] = 40;
    span.is_periodic[0] = 1;
    EXPECT_TRUE(placement(hp, span).empty());
    edge.is_periodic[0] = 1;
    EXPECT_FALSE(placement(hp, edge).empty()) << "a periodic patch narrower than the parent";
}

// The message of the runtime_error fn throws (amrex.throw_exception = 1).
std::string abort_message (const std::function<void()>& fn)
{
    try { fn(); } catch (const std::runtime_error& e) { return e.what(); }
    return "";
}

// UT-9.1.2: a version-1 checkpoint and a missing one abort naming the key.
TEST (NestParent, OldOrMissingCheckpointAbortsNamingTheKey)
{
    NestRun p;
    p.in = parent_inputs(false);
    p.make();
    const ScopedDir dir{write_parent(p, "ut_nest_v1", [] (const amrex::MultiFab&, int c, int i, int j, int k) {
        return 1.0 + c + 0.01 * i - 0.02 * j + 0.03 * k;
    })};
    // The I/O rank turns the Header into version 1 (no zlevels line).
    bool rewritten = false;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ifstream is(dir.name + "/Header");
        std::ostringstream os;
        std::string line;
        int n = 0;
        while (std::getline(is, line)) {
            if (n++ == 0) {
                const auto sp = line.rfind(' ');
                rewritten = line.substr(sp + 1) == "2";
                line = line.substr(0, sp) + " 1";
            }
            if (line.rfind("zlevels ", 0) == 0) { continue; }
            os << line << "\n";
        }
        is.close();
        std::ofstream(dir.name + "/Header") << os.str();
    }
    amrex::ParallelDescriptor::ReduceBoolOr(rewritten);
    ASSERT_TRUE(rewritten) << "the Header's first line does not end in format version 2";
    // A version-1 checkpoint still restarts (io.restart), bit for bit.
    {
        gust::State r;
        r.define(p.ba, p.dm);
        gust::read_checkpoint(dir.name, r, p.in.mesh);
        auto rf = fields6(r);
        auto pf = fields6(p.q);
        for (std::size_t n = 0; n < rf.size(); ++n) {
            const amrex::MultiFab& a = *rf[n].first;
            const amrex::MultiFab& b = *pf[n].first;
            amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(d, a, 0, 0, 1, 0);
            amrex::MultiFab::Subtract(d, b, 0, 0, 1, 0);
            EXPECT_EQ(d.norminf(0, 0), amrex::Real(0)) << "field " << rf[n].second;   // norminf reduces over ranks
        }
    }
    NestRun f;
    f.in = patch_inputs(p.in);
    f.make();
    const std::string m1 = abort_message([&] { nest(dir.name, f); });
    EXPECT_NE(m1.find("bc.nest_parent"), std::string::npos) << m1;
    EXPECT_NE(m1.find("format version 1"), std::string::npos) << m1;
    const std::string m2 = abort_message([&] { nest("ut_nest_no_such_dir", f); });
    EXPECT_NE(m2.find("bc.nest_parent"), std::string::npos) << m2;
}

// The patch values gathered on every rank, by global index (valid points;
// a nodal point shared by two boxes holds their mean, which is its value).
std::vector<double> gather_valid (const amrex::MultiFab& mf)
{
    const amrex::Box dom = mf.boxArray().minimalBox();
    const auto npts = static_cast<std::size_t>(dom.numPts());
    std::vector<double> sum(npts, 0.0), cnt(npts, 0.0);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox h(bx, 1, amrex::The_Pinned_Arena());
        h.copy<amrex::RunOn::Device>(mf[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const auto m = static_cast<std::size_t>(dom.index(amrex::IntVect(i, j, k)));
            sum[m] += static_cast<double>(h(amrex::IntVect(i, j, k)));
            cnt[m] += 1;
        });
    }
    const auto comm = amrex::ParallelDescriptor::Communicator();
    amrex::ParallelAllReduce::Sum(sum.data(), static_cast<int>(npts), comm);
    amrex::ParallelAllReduce::Sum(cnt.data(), static_cast<int>(npts), comm);
    for (std::size_t m = 0; m < npts; ++m) { sum[m] /= cnt[m]; }
    return sum;
}

// A smooth parent field that is not bilinear (so the layouts' results are not
// all the same closed form by accident).
double wavy (double x, double y, int k, int c)
{
    return 5.0 + c + 0.2 * k + std::sin(0.003 * x + 0.4 * c) * std::cos(0.002 * y - 0.1 * k) + 1.0e-3 * x;
}

// UT-9.1.3: the parent's and the patch's layouts and rank maps change; the
// patch values are the same bit for bit.
TEST (NestParent, ValuesDoNotDependOnLayoutOrRanks)
{
    struct Layout { int parent_mgs; int patch_mgs; bool reverse; };
    const Layout layouts[] = {{7, 9, false}, {64, 64, false}, {5, 6, true}, {11, 4, true}};
    std::vector<std::vector<double>> ref;
    int n = 0;
    for (const auto& L : layouts) {
        NestRun p;
        p.in = parent_inputs(true, L.parent_mgs);
        p.make(L.reverse);
        const ScopedDir dir{write_parent(p, "ut_nest_layout" + std::to_string(n), [&] (const amrex::MultiFab& mf, int c, int i, int j, int k) {
            return wavy(xpos(p.geom, mf, i, 0), xpos(p.geom, mf, j, 1), k, c);
        })};
        NestRun f;
        f.in = patch_inputs(p.in, L.patch_mgs);
        f.make(!L.reverse);
        nest(dir.name, f);
        std::vector<std::vector<double>> got;
        for (auto [mf, c] : fields6(f.q)) { got.push_back(gather_valid(*mf)); }
        if (n == 0) {
            ref = got;
        } else {
            for (std::size_t c = 0; c < got.size(); ++c) {
                long ndiff = 0;
                for (std::size_t m = 0; m < got[c].size(); ++m) { ndiff += (got[c][m] != ref[c][m]) ? 1 : 0; }
                EXPECT_EQ(ndiff, 0) << "layout " << n << ", field " << c;
            }
        }
        ++n;
    }
}

// Beyond the plan: across a periodic seam of the parent in x. The expected
// value is the bilinear interpolation of the parent's discrete values,
// wrapped, computed here from the parent's indices.
TEST (NestParent, PeriodicSeamUsesThePeriodicNeighbours)
{
    NestRun p;
    p.in = parent_inputs(false);
    p.in.mesh.is_periodic = {1, 0, 0};
    p.make();
    const int N = p.in.mesh.n_cell[0];
    auto g = [] (int I, int J, int k, int c) {   // the parent's value at its index (I wrapped)
        return 3.0 + c + 0.1 * k + std::sin(0.7 * I + 0.3 * c) + 0.05 * J * J - 0.02 * I * J;
    };
    auto wrap = [N] (int I) { return ((I % N) + N) % N; };
    const ScopedDir dir{write_parent(p, "ut_nest_seam", [&] (const amrex::MultiFab&, int c, int i, int j, int k) {
        return g(wrap(i), j, k, c);
    })};
    // Patch 1: parent cells 0 to 10 in x, not periodic (its ghosts cross the
    // seam); patch 2: periodic, spanning the parent.
    for (int span = 0; span < 2; ++span) {
        NestRun f;
        f.in = patch_inputs(p.in);
        f.in.mesh.prob_lo[0] = -1000;
        f.in.mesh.prob_hi[0] = (span != 0) ? 1000 : 0;
        f.in.mesh.n_cell[0] = (span != 0) ? 40 : 20;
        f.in.mesh.is_periodic[0] = span;
        f.make();
        nest(dir.name, f);
        const int oy = 2;   // the patch's y offset in parent cells
        const Check r = check_patch(f.q, [&] (const amrex::MultiFab& mf, int c, int i, int j, int k) {
            // The parent's (fractional) index of patch point (i, j): nodal
            // i / 2, cell-centred (i + 1/2) / 2 - 1/2, plus the offset.
            auto pos = [&] (int ii, int d, int o) {
                return mf.ixType().nodeCentered(d) ? 0.5 * ii + o : 0.5 * (ii + 0.5) - 0.5 + o;
            };
            const double X = pos(i, 0, 0), Y = pos(j, 1, oy);
            const int I = static_cast<int>(std::floor(X)), J = static_cast<int>(std::floor(Y));
            const double tx = X - I, ty = Y - J;
            return (1 - tx) * (1 - ty) * g(wrap(I), J, k, c) + tx * (1 - ty) * g(wrap(I + 1), J, k, c)
                 + (1 - tx) * ty * g(wrap(I), J + 1, k, c) + tx * ty * g(wrap(I + 1), J + 1, k, c);
        });
        EXPECT_LT(r.err, kTol) << "span " << span;
        EXPECT_EQ(r.zghost_bad, 0) << "span " << span;
    }
}

} // namespace
