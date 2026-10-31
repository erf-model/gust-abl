#include "StateIO.H"

#include "FileMagic.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Utility.H>
#include <AMReX_VisMF.H>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

namespace gust {

namespace {

constexpr const char* kCheckpointMagic = "gust-checkpoint";
constexpr int kCheckpointVersion = 1;

// Hex-float text, so a restart reads back the exact bits.
std::string hexfloat (amrex::Real x)
{
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%a", static_cast<double>(x));
    return buf;
}

void restart_error (const std::string& dir, const std::string& what)
{
    amrex::Abort("gust: io.restart = " + dir + ": " + what);
}

} // namespace

const amrex::Vector<std::string>& plot_var_names ()
{
    static const amrex::Vector<std::string> names{"u", "v", "w", "p", "theta", "tke", "z"};
    return names;
}

std::string write_plotfile (const std::string& prefix, const State& state,
                            const amrex::Geometry& geom, const amrex::MultiFab* detJ,
                            const amrex::MultiFab* canopy_drag, const amrex::MultiFab* canopy_lad,
                            const amrex::MultiFab* building_beta, const amrex::MultiFab* building_wall)
{
    const std::string name = amrex::Concatenate(prefix, state.step, 5);
    const amrex::BoxArray& ba = state.p.boxArray();
    const int ncomp = 7 + (detJ ? 1 : 0) + (canopy_drag ? 1 : 0) + (canopy_lad ? 1 : 0)
                      + (building_beta ? 1 : 0) + (building_wall ? 1 : 0);
    amrex::MultiFab plt(ba, state.p.DistributionMap(), ncomp, 0);

    for (amrex::MFIter mfi(plt, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto out = plt.array(mfi);
        const auto u = state.u.const_array(mfi);
        const auto v = state.v.const_array(mfi);
        const auto w = state.w.const_array(mfi);
        const auto p = state.p.const_array(mfi);
        const auto th = state.theta.const_array(mfi);
        const auto k = state.tke.const_array(mfi);
        const auto z = state.z_phys_nd.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int kk) noexcept
        {
            constexpr amrex::Real half = amrex::Real(0.5);
            constexpr amrex::Real eighth = amrex::Real(0.125);
            out(i, j, kk, 0) = half * (u(i, j, kk) + u(i + 1, j, kk));
            out(i, j, kk, 1) = half * (v(i, j, kk) + v(i, j + 1, kk));
            out(i, j, kk, 2) = half * (w(i, j, kk) + w(i, j, kk + 1));
            out(i, j, kk, 3) = p(i, j, kk);
            out(i, j, kk, 4) = th(i, j, kk);
            out(i, j, kk, 5) = k(i, j, kk);
            out(i, j, kk, 6) = eighth * (z(i, j, kk)     + z(i + 1, j, kk)
                                       + z(i, j + 1, kk) + z(i + 1, j + 1, kk)
                                       + z(i, j, kk + 1) + z(i + 1, j, kk + 1)
                                       + z(i, j + 1, kk + 1) + z(i + 1, j + 1, kk + 1));
        });
    }

    amrex::Vector<std::string> names = plot_var_names();
    if (detJ != nullptr) {
        amrex::MultiFab::Copy(plt, *detJ, 0, 7, 1, 0);
        names.emplace_back("detJ");
    }
    if (canopy_drag != nullptr) {
        amrex::MultiFab::Copy(plt, *canopy_drag, 0, static_cast<int>(names.size()), 1, 0);
        names.emplace_back("canopy_drag");
    }
    if (canopy_lad != nullptr) {
        amrex::MultiFab::Copy(plt, *canopy_lad, 0, static_cast<int>(names.size()), 1, 0);
        names.emplace_back("canopy_lad");
    }
    if (building_beta != nullptr) {
        amrex::MultiFab::Copy(plt, *building_beta, 0, static_cast<int>(names.size()), 1, 0);
        names.emplace_back("building_beta");
    }
    if (building_wall != nullptr) {
        amrex::MultiFab::Copy(plt, *building_wall, 0, static_cast<int>(names.size()), 1, 0);
        names.emplace_back("building_wall");
    }
    amrex::WriteSingleLevelPlotfile(name, plt, names, geom, 0.0, state.step);
    amrex::VisMF::Write(state.z_phys_nd, name + "/Level_0/z_phys_nd");
    return name;
}

std::string write_checkpoint (const std::string& prefix, State& state, const MeshInputs& mesh)
{
    const std::string name = amrex::Concatenate(prefix, state.step, 5);
    amrex::PreBuildDirectorHierarchy(name, "Level_", 1, true);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(name + "/Header");
        if (!os.good()) { amrex::FileOpenFailed(name + "/Header"); }
        os << kCheckpointMagic << " " << kCheckpointVersion << "\n";
        os << "real_bytes " << sizeof(amrex::Real) << "\n";
        os << "step " << state.step << "\n";
        os << "n_cell " << mesh.n_cell[0] << " " << mesh.n_cell[1] << " " << mesh.n_cell[2] << "\n";
        os << "prob_lo";
        for (int d = 0; d < 3; ++d) { os << " " << hexfloat(mesh.prob_lo[d]); }
        os << "\nprob_hi";
        for (int d = 0; d < 3; ++d) { os << " " << hexfloat(mesh.prob_hi[d]); }
        os << "\nis_periodic " << mesh.is_periodic[0] << " " << mesh.is_periodic[1] << " "
           << mesh.is_periodic[2] << "\n";
        state.p.boxArray().writeOn(os);
        os << "\n";
    }

    // The solution fields with their ghost images refreshed: the solver
    // leaves some images stale between fills (theta and k after their
    // sweeps), and a reader of the ghosts would take them. The restart reads
    // valid data only. A copy, so the running state is not touched (a run
    // that writes a checkpoint and goes on stays bit for bit the run without
    // it). The node heights are written as they are (their ghost nodes are
    // the terrain's, read back on restart).
    const amrex::Periodicity period = make_geometry(mesh).periodicity();
    for (auto& f : state.fields()) {
        const std::string file = amrex::MultiFabFileFullPrefix(0, name, "Level_", f.name);
        if (f.mf == &state.z_phys_nd) {
            amrex::VisMF::Write(*f.mf, file);
            continue;
        }
        amrex::MultiFab out(f.mf->boxArray(), f.mf->DistributionMap(), f.mf->nComp(), f.mf->nGrowVect());
        amrex::MultiFab::Copy(out, *f.mf, 0, 0, f.mf->nComp(), f.mf->nGrowVect());
        out.FillBoundary(period);
        amrex::VisMF::Write(out, file);
    }
    return name;
}

void read_checkpoint (const std::string& dir, State& state, const MeshInputs& mesh)
{
    amrex::Vector<char> buf;
    amrex::ParallelDescriptor::ReadAndBcastFile(dir + "/Header", buf);
    std::istringstream is(std::string(buf.dataPtr()));

    std::string magic;
    int version = 0;
    is >> magic >> version;
    if (!is || !is_file_magic(magic, kCheckpointMagic) || version != kCheckpointVersion) {
        restart_error(dir, "not a GUST checkpoint of format version "
                           + std::to_string(kCheckpointVersion));
    }

    // Each Header line is "key value(s)"; a missing key or value is a corrupt Header.
    auto expect_key = [&] (const char* want) {
        std::string key;
        is >> key;
        if (!is || key != want) {
            restart_error(dir, std::string("corrupt Header: expected '") + want + "'");
        }
    };
    auto check_stream = [&] (const char* what) {
        if (!is) { restart_error(dir, std::string("corrupt Header: bad value for '") + what + "'"); }
    };
    auto read_real = [&] (const char* what) {
        std::string tok;
        is >> tok;
        char* end = nullptr;
        const double x = std::strtod(tok.c_str(), &end);
        if (!is || tok.empty() || *end != '\0') {
            restart_error(dir, std::string("corrupt Header: bad value for '") + what + "'");
        }
        return static_cast<amrex::Real>(x);
    };

    int real_bytes = 0;
    expect_key("real_bytes");
    is >> real_bytes;
    check_stream("real_bytes");
    if (real_bytes != static_cast<int>(sizeof(amrex::Real))) {
        restart_error(dir, "checkpoint was written with " + std::to_string(real_bytes)
                           + "-byte reals; this build uses " + std::to_string(sizeof(amrex::Real)));
    }

    expect_key("step");
    is >> state.step;
    check_stream("step");

    amrex::IntVect n_cell;
    expect_key("n_cell");
    is >> n_cell[0] >> n_cell[1] >> n_cell[2];
    check_stream("n_cell");
    std::array<amrex::Real, 3> lo{}, hi{};
    expect_key("prob_lo");
    for (int d = 0; d < 3; ++d) { lo[d] = read_real("prob_lo"); }
    expect_key("prob_hi");
    for (int d = 0; d < 3; ++d) { hi[d] = read_real("prob_hi"); }
    amrex::Array<int, 3> per{};
    expect_key("is_periodic");
    is >> per[0] >> per[1] >> per[2];
    check_stream("is_periodic");

    if (n_cell != mesh.n_cell) {
        std::ostringstream os;
        os << "checkpoint mesh.n_cell = " << n_cell << " differs from the inputs' " << mesh.n_cell;
        restart_error(dir, os.str());
    }
    for (int d = 0; d < 3; ++d) {
        if (lo[d] != mesh.prob_lo[d] || hi[d] != mesh.prob_hi[d]) {
            restart_error(dir, "checkpoint mesh.prob_lo/prob_hi differ from the inputs");
        }
        if (per[d] != mesh.is_periodic[d]) {
            restart_error(dir, "checkpoint mesh.is_periodic differs from the inputs");
        }
    }

    amrex::BoxArray ba_chk;
    ba_chk.readFrom(is);
    if (!is) { restart_error(dir, "could not read the BoxArray from the Header"); }
    const amrex::DistributionMapping dm_chk(ba_chk);

    for (auto& f : state.fields()) {
        amrex::MultiFab& dst = *f.mf;
        amrex::MultiFab src(amrex::convert(ba_chk, dst.ixType()), dm_chk,
                            dst.nComp(), dst.nGrowVect());
        amrex::VisMF::Read(src, amrex::MultiFabFileFullPrefix(0, dir, "Level_", f.name));
        // Valid data exactly; nodal heights also with their ghost nodes,
        // since those come from the terrain and are not recomputed. A copy
        // from grown source boxes can land a checkpoint ghost on a valid
        // destination point, so ghosts go first and valid data overwrites.
        if (&dst == &state.z_phys_nd) {
            dst.ParallelCopy(src, 0, 0, dst.nComp(), ngrow_nodal, ngrow_nodal);
        }
        dst.ParallelCopy(src, 0, 0, dst.nComp(), 0, 0);
    }
}

bool same_directory (const std::string& target, const std::string& restart)
{
    std::error_code ec;
    const bool same = std::filesystem::equivalent(target, restart, ec);
    return !ec && same;
}

} // namespace gust
