#include "SurfaceTemperature.H"

#include "SurfaceRoughness.H"
#include "TerrainSurface.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

namespace gust {

namespace {

// Index of cell i along a direction: periodic wraps, otherwise clamped to the domain (the column a ghost column
// images), as SurfaceRoughness.cpp.
int cell_index (int i, int lo, int hi, bool periodic)
{
    if (!periodic) { return std::min(std::max(i, lo), hi); }
    const int n = hi - lo + 1;
    return lo + (((i - lo) % n) + n) % n;
}

} // namespace

void check_surface_temperature_map (const char* op, const amrex::MultiFab& map, const amrex::MultiFab& like)
{
    const amrex::IntVect ng = map.nGrowVect();
    if (map.nComp() != 1 || !map.ixType().cellCentered() || ng[0] < 1 || ng[1] < 1
        || map.boxArray() != ground_slab(amrex::convert(like.boxArray(), amrex::IntVect(0)))
        || map.DistributionMap() != like.DistributionMap()) {
        amrex::Abort(std::string("gust: ") + op + ": the surface temperature map is not on this level's ground slab");
    }
}

void read_surface_temperature_map (const std::string& file, const std::string& key, int k, double power,
                                   const amrex::Geometry& geom, const amrex::BoxArray& cell_ba,
                                   const amrex::DistributionMapping& dm, amrex::MultiFab& theta_s)
{
    const TerrainPoints pts = read_xyz_points(file, key);
    const auto n = static_cast<long>(pts.z.size());
    std::vector<std::string> err;
    const std::string range = surface_temperature_range();
    if (pts.nonfinite > 0) {
        err.push_back(std::to_string(pts.nonfinite)
                      + " lines hold a value that is not finite; every theta_s must be finite and in " + range);
    }
    long bad = 0;
    for (const double v : pts.z) {
        if (!(v >= kSurfaceTemperatureMin && v <= kSurfaceTemperatureMax)) { ++bad; }
    }
    if (bad > 0) { err.push_back(std::to_string(bad) + " points have theta_s outside " + range); }
    if (k < 1 || k > n) {
        err.push_back("terrain.idw_neighbors = " + std::to_string(k) + " must be between 1 and the number of points ("
                      + std::to_string(n) + ")");
    }
    if (!err.empty()) {
        std::string all = "gust: " + key + " = " + file + ":";
        for (const auto& e : err) { all += "\n  " + e; }
        amrex::Abort(all);
    }
    const IdwInterpolator idw(pts, k, power);

    theta_s.clear();
    theta_s.define(ground_slab(cell_ba), dm, 1, amrex::IntVect(1, 1, 0));
    const amrex::Box& dom = geom.Domain();
    const bool perx = geom.isPeriodic(0), pery = geom.isPeriodic(1);
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    for (amrex::MFIter mfi(theta_s, false); mfi.isValid(); ++mfi) {
        // Host work (the search), copied to the device; ghost columns included.
        const amrex::Box fb = mfi.fabbox();
        amrex::FArrayBox h(fb, 1, amrex::The_Pinned_Arena());
        const auto ha = h.array();
        amrex::LoopOnCpu(fb, [&] (int i, int j, int kk) {
            const int ii = cell_index(i, dom.smallEnd(0), dom.bigEnd(0), perx);
            const int jj = cell_index(j, dom.smallEnd(1), dom.bigEnd(1), pery);
            const double x = static_cast<double>(plo[0]) + (ii + 0.5) * static_cast<double>(dx[0]);
            const double y = static_cast<double>(plo[1]) + (jj + 0.5) * static_cast<double>(dx[1]);
            ha(i, j, kk) = static_cast<amrex::Real>(idw(x, y));
        });
        theta_s[mfi].copy<amrex::RunOn::Device>(h, fb, 0, fb, 0, 1);
        // h is freed at the end of this iteration: the copy must be done.
        amrex::Gpu::streamSynchronize();
    }
}

std::string surface_temperature_summary (const amrex::MultiFab& theta_s)
{
    amrex::ReduceOps<amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpSum> ops;
    amrex::ReduceData<amrex::Real, amrex::Real, double> data(ops);   // the sum in double for the mean
    for (amrex::MFIter mfi(theta_s, false); mfi.isValid(); ++mfi) {
        const auto a = theta_s.const_array(mfi);
        ops.eval(mfi.validbox(), data,
                 [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Real, amrex::Real, double> {
            return {a(i, j, k), a(i, j, k), static_cast<double>(a(i, j, k))};
        });
    }
    const auto r = data.value();
    amrex::Real lo = amrex::get<0>(r), hi = amrex::get<1>(r);
    double sum = amrex::get<2>(r);
    amrex::ParallelAllReduce::Min(lo, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Max(hi, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(sum, amrex::ParallelDescriptor::Communicator());
    const double mean = sum / static_cast<double>(theta_s.boxArray().numPts());
    char b[128];
    std::snprintf(b, sizeof(b), "min %.6g, max %.6g, mean %.6g K", static_cast<double>(lo), static_cast<double>(hi), mean);
    return b;
}

} // namespace gust
