// The buildings' terms of the closure (BuildingTke.H), beyond the plan, sub-phase 5.4 (UT-5.4.x):
//   a wall cell's k row is a (k - u*^2 / Cmu0^2) with the law's u*, zero at the wall value, in the
//   residual and in the implicit rows, and it is added to the ground's row in the first cells
//   (UT-5.4.1); the closure's distance in wall cells is the fluid centroid's, (1 - beta) Delta_n / 2,
//   at beta 0, 0.5 and 0.9, floored at buildings.z0 (UT-5.4.2); solid cells hold k at its floor
//   (UT-5.4.3); where the ground is nearer than any building the wall distance keeps its bits
//   (UT-5.4.4); the distance to a block equals the closed form of a box, on every box layout and
//   across a periodic seam, and that to an L with a step the smallest over its boxes; the diffusive
//   flux of k and its couplings are closed on blocked faces.
// After 5.4 (tracker row N-88): the closure's d + z0 is the smaller of the ground's and the nearest
//   wall's, each with its own roughness, with a uniform ground z0 and with a map; the distance beside a
//   step between two roofs whose lower one ends inside a layer, and to a block whose sides and top lie
//   inside cells, are the boxes'; the norm of R_k leaves the solid cells out; the rows through
//   residual() without a surface layer.
#include "BuildingTke.H"
#include "Closure.H"
#include "ClosureWake.H"
#include "Reductions.H"
#include "BoundaryConditions.H"
#include "PoissonTestHelpers.H"
#include "PseudoTime.H"
#include "Residual.H"
#include "Simplec.H"
#include "SurfaceRoughness.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#ifdef _WIN32
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-12 : 2e-4;
constexpr double kKappa = 0.41, kCmu0 = 0.5562, kFloor = 1e-6;

std::string write_text (const std::string& name, const std::function<void(std::ostream&)>& fill)
{
    static const long id = [] () {
        long pid = static_cast<long>(::getpid());
        amrex::ParallelDescriptor::Bcast(&pid, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        return pid;
    }();
    const fs::path dir = fs::temp_directory_path() / ("gust_building_tke_" + std::to_string(id));
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        // Binary: the same bytes, and so the same hash of the file, on every platform.
        std::ofstream os(dir / name, std::ios::binary);
        os.precision(17);
        fill(os);
    }
    amrex::ParallelDescriptor::Barrier();
    return (dir / name).string();
}

// 8 x 8 x 24 cells of 10 m (dz: the layers' height in place of 10 m).
PCase mesh_case (int mgs, bool periodic = false, double dz = 10)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(8, 8, 24);
    c.in.mesh.prob_hi = {80, 80, static_cast<amrex::Real>(24 * dz)};
    c.in.mesh.is_periodic = {periodic ? 1 : 0, periodic ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

// A map in ERF's format on n x n points 20 m apart from 0: tile (ix, iy) covers 20 ix - 10 .. 20 ix + 10
// along x, so tile 2 is the cells 3 and 4. height(ix, iy) in m.
// spacing: the points' distance in place of 20 m (15 m: tile 2 covers 22.5 .. 37.5 m, three quarters of the
// cells 2 and 3).
gust::BuildingsOptions tile_options (const std::string& name, const std::function<double(int, int)>& height, int n = 5, double z0 = 0.1,
                                     double spacing = 20)
{
    gust::BuildingsOptions o;
    o.type = gust::BuildingsType::immersed_forcing;
    o.type_name = "immersed_forcing";
    o.small_volfrac = 0;
    o.z0 = z0;
    o.file = write_text(name, [&] (std::ostream& os) {
        os << n << "\n" << n << "\n";
        for (int r = 0; r < 2; ++r) { for (int p = 0; p < n; ++p) { os << spacing * p << "\n"; } }
        for (int i = 0; i < n; ++i) { for (int j = 0; j < n; ++j) { os << height(i, j) << "\n"; } }
    });
    return o;
}

// The value at one point, on every rank, or -1e30 outside the field's points.
double value_at (const amrex::MultiFab& mf, int i, int j, int k, int n = 0)
{
    double v = -1e30;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(amrex::IntVect(i, j, k))) { v = static_cast<double>(mf.const_array(mfi)(i, j, k, n)); }
    }
    amrex::ParallelDescriptor::ReduceRealMax(v);
    return v;
}

// A case's buildings with their wall law and closure terms.
struct Set
{
    PCase c;
    gust::BuildingsOptions o;
    gust::Buildings b;
    gust::BuildingForcing law;
    gust::BuildingTke tke;
    amrex::MultiFab ox, oy, oz;
    double dz;
    Set (int mgs, const gust::BuildingsOptions& opt, bool periodic = false, double layer = 10)
        : c(mesh_case(mgs, periodic, layer)), o(opt), dz(layer)
    {
        b.build(o, c.geom, c.z, c.ba, c.dm);
        law.build(b, o, c.geom, c.mesh.detJ);
        tke.build(b, law, o, c.geom, c.z, amrex::Real(kFloor), amrex::Real(kCmu0));
        amrex::MultiFab* op[3] = {&ox, &oy, &oz};
        for (int d = 0; d < 3; ++d) {
            op[d]->define(b.beta_face(d).boxArray(), c.dm, 1, 1);
            op[d]->setVal(1);
            amrex::MultiFab::LinComb(*op[d], 1, *op[d], 0, -1, b.beta_face(d), 0, 0, 1, 0);
            op[d]->FillBoundary(c.geom.periodicity());
        }
    }
    [[nodiscard]] amrex::Array<const amrex::MultiFab*, 3> open () const { return {&ox, &oy, &oz}; }
};

// The distance from a point to the nearest of some boxes it lies outside of ({xlo, xhi, ylo, yhi, zhi}, from z = 0),
// with the nearest image along x and y when periodic (period 80 m).
double box_distance (const std::vector<std::array<double, 5>>& boxes, double x, double y, double z, bool periodic = false)
{
    double best = 1e30;
    for (const auto& bx : boxes) {
        for (int sx = (periodic ? -1 : 0); sx <= (periodic ? 1 : 0); ++sx) {
            for (int sy = (periodic ? -1 : 0); sy <= (periodic ? 1 : 0); ++sy) {
                const double ex = std::max({bx[0] + 80 * sx - x, 0.0, x - bx[1] - 80 * sx});
                const double ey = std::max({bx[2] + 80 * sy - y, 0.0, y - bx[3] - 80 * sy});
                const double ez = std::max(0.0, z - bx[4]);
                best = std::min(best, std::sqrt(ex * ex + ey * ey + ez * ez));
            }
        }
    }
    return best;
}

// The largest |d_b - closed form| / closed form over the fluid cells that are no wall cells, and how many there are.
std::pair<double, amrex::Long> distance_error (const Set& s, const std::vector<std::array<double, 5>>& boxes, bool periodic = false)
{
    double err = 0;
    amrex::Long n = 0;
    for (amrex::MFIter mfi(s.tke.distance(), false); mfi.isValid(); ++mfi) {
        const auto d = s.tke.distance().const_array(mfi);
        const auto beta = s.b.beta().const_array(mfi);
        const auto wl = s.b.wall().const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            if (beta(i, j, k) > 0 || wl(i, j, k, gust::kWallMask) != 0) { return; }
            const double want = box_distance(boxes, 10 * i + 5, 10 * j + 5, s.dz * (k + 0.5), periodic);
            err = std::max(err, std::abs(static_cast<double>(d(i, j, k)) - want) / want);
            ++n;
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(err);
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return {err, n};
}

// UT-5.4.2: the fluid centroid's distance in wall cells, and the floor.
TEST (BuildingTke, WallCellDistanceIsTheFluidCentroids)
{
    // Roofs at 20, 25 and 29 m: the cell above 20 m is fluid, the others are half and nine tenths solid.
    const auto height = [] (int i, int j) { return (j != 2) ? 0 : (i == 1) ? 20 : (i == 2) ? 25 : (i == 3) ? 29 : 0; };
    for (const int mgs : {8, 3}) {
        Set s(mgs, tile_options("roofs.txt", height));
        EXPECT_NEAR(value_at(s.tke.distance(), 1, 3, 2), 5.0, 100 * rtol) << "beta 0";
        EXPECT_NEAR(value_at(s.tke.distance(), 3, 3, 2), 2.5, 100 * rtol) << "beta 0.5";
        EXPECT_NEAR(value_at(s.tke.distance(), 5, 3, 2), 0.5, 100 * rtol) << "beta 0.9";
        EXPECT_EQ(value_at(s.tke.distance(), 3, 3, 1), static_cast<double>(gust::kBuildingFar)) << "a solid cell has none";
        // A side's wall cell: half a cell from the side.
        ASSERT_EQ(static_cast<int>(value_at(s.b.wall(), 3, 5, 1, gust::kWallMask)), gust::kWallYlo);
        EXPECT_NEAR(value_at(s.tke.distance(), 3, 5, 1), 5.0, 100 * rtol);
    }
    // Ghost cells outside a non-periodic side take the nearest domain column's value. A block touching the
    // low x side (the tile -10..10 m; the log warns of it) has a 5 m wall cell at j = 5 in the first column.
    Set edge(4, tile_options("edge.txt", [] (int i, int j) { return (i == 0 && j == 2) ? 45 : 0; }));
    // The values are reduced over the ranks before they are compared, so every rank sees a failure.
    double in5 = -1e30, out5 = -1e30, in7 = -1e30, out7 = -1e30;
    for (amrex::MFIter mfi(edge.tke.distance(), false); mfi.isValid(); ++mfi) {
        if (!mfi.validbox().contains(amrex::IntVect(0, 5, 1)) || !mfi.validbox().contains(amrex::IntVect(0, 7, 1))) { continue; }
        const auto d = edge.tke.distance().const_array(mfi);
        in5 = d(0, 5, 1);
        out5 = d(-1, 5, 1);
        in7 = d(0, 7, 1);
        out7 = d(-1, 7, 1);
    }
    for (double* x : {&in5, &out5, &in7, &out7}) { amrex::ParallelDescriptor::ReduceRealMax(*x); }
    EXPECT_NEAR(in5, 5.0, 100 * rtol);
    EXPECT_EQ(out5, in5) << "the ghost outside x = 0";
    EXPECT_GT(in7, 5.0);
    EXPECT_EQ(out7, in7);
    Set rough(4, tile_options("roofs.txt", height, 5, 1.0));
    EXPECT_NEAR(value_at(rough.tke.distance(), 5, 3, 2), 1.0, 100 * rtol) << "floored at buildings.z0 = 1 m";
    EXPECT_NEAR(value_at(rough.tke.distance(), 3, 3, 2), 2.5, 100 * rtol);
}

// The distance to one block (the cells 3, 4 in x and y, 45 m high) is that to its box.
TEST (BuildingTke, DistanceToABlockIsTheBoxs)
{
    const std::vector<std::array<double, 5>> box = {{30, 50, 30, 50, 45}};
    for (const int mgs : {8, 5, 3}) {
        Set s(mgs, tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }));
        const auto [err, n] = distance_error(s, box);
        EXPECT_GT(n, 1000);
        EXPECT_LT(err, 100 * rtol) << "boxes " << mgs;
        EXPECT_EQ(s.tke.distance_passes(), 19) << "the same on every layout";
        // By hand: above the roof's corner, and off the top edge of a side (the side ends at 45 m, mid-layer).
        EXPECT_NEAR(value_at(s.tke.distance(), 3, 3, 7), 30.0, 100 * rtol * 30);
        EXPECT_NEAR(value_at(s.tke.distance(), 5, 3, 5), std::hypot(5.0, 10.0), 100 * rtol * 12);
        EXPECT_NEAR(value_at(s.tke.distance(), 6, 6, 5), std::sqrt(15.0 * 15 + 15 * 15 + 10 * 10), 100 * rtol * 24);
    }
}

// The same block across the periodic seam in x (the tile -10..10 m).
TEST (BuildingTke, DistanceAcrossAPeriodicSeam)
{
    const std::vector<std::array<double, 5>> box = {{-10, 10, 30, 50, 45}};
    for (const int mgs : {8, 3}) {
        Set s(mgs, tile_options("seam.txt", [] (int i, int j) { return (i == 0 && j == 2) ? 45 : 0; }, 4), true);
        const auto [err, n] = distance_error(s, box, true);
        EXPECT_GT(n, 1000);
        EXPECT_LT(err, 100 * rtol) << "boxes " << mgs;
        EXPECT_NEAR(value_at(s.tke.distance(), 5, 3, 1), 15.0, 100 * rtol * 15) << "nearer through the seam: x = 55 m to the image at 70 m, not back to 10 m";
    }
}

// An L of three 45 m tiles and a 20 m step: not convex. The smallest distance over the boxes, to round-off, and
// the number of passes: walls passed through face neighbours alone would give the same distances here in more passes.
TEST (BuildingTke, DistanceToAnLWithAStep)
{
    const std::vector<std::array<double, 5>> boxes = {{30, 50, 30, 50, 45}, {10, 30, 30, 50, 45}, {30, 50, 10, 30, 45}, {50, 70, 30, 50, 20}};
    Set s(3, tile_options("ell.txt", [] (int i, int j) {
        if ((i == 2 && j == 2) || (i == 1 && j == 2) || (i == 2 && j == 1)) { return 45; }
        return (i == 3 && j == 2) ? 20 : 0;
    }));
    const auto [err, n] = distance_error(s, boxes);
    EXPECT_GT(n, 1000);
    EXPECT_LT(err, 100 * rtol) << "largest relative difference from the boxes' distance";
    EXPECT_EQ(s.tke.distance_passes(), 19);
}

// UT-5.4.4: the limit leaves the cells nearer to the ground than to a building as they were.
TEST (BuildingTke, LimitKeepsTheGroundsDistanceWhereItIsNearer)
{
    for (const int mgs : {8, 3}) {
        Set s(mgs, tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }));
        amrex::MultiFab dist, before;
        gust::height_above_ground(s.c.z, dist);
        before.define(dist.boxArray(), dist.DistributionMap(), 1, 1);
        amrex::MultiFab::Copy(before, dist, 0, 0, 1, 1);
        const amrex::Long n = s.tke.limit_wall_distance(dist, gust::Roughness(amrex::Real(0.1)));   // buildings.z0, so no shift
        amrex::Long want = 0, wrong = 0;
        for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
            const auto a = dist.const_array(mfi);
            const auto o = before.const_array(mfi);
            const auto db = s.tke.distance().const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (db(i, j, k) < o(i, j, k)) {
                    ++want;
                    if (a(i, j, k) != db(i, j, k)) { ++wrong; }
                } else if (a(i, j, k) != o(i, j, k)) {   // the ground's distance, bit for bit
                    ++wrong;
                }
            });
        }
        amrex::ParallelDescriptor::ReduceLongSum(want);
        amrex::ParallelDescriptor::ReduceLongSum(wrong);
        EXPECT_EQ(wrong, 0) << "cells with neither the building's distance where it is nearer nor the ground's bits";
        EXPECT_EQ(n, want);
        EXPECT_GT(n, 0);
        EXPECT_EQ(value_at(dist, 0, 0, 0), 5.0) << "first cell far from the block: 5 m above the ground, 42 m from the block";
        EXPECT_EQ(value_at(dist, 3, 3, 1), 15.0) << "a solid cell keeps the ground's distance";
        EXPECT_NEAR(value_at(dist, 5, 3, 3), 5.0, 100 * rtol) << "a side's wall cell at 35 m: 5 m from the side";
        EXPECT_NEAR(value_at(dist, 3, 3, 7), 30.0, 100 * rtol * 30) << "above the roof: 30 m, not 75 m";
    }
}

// The closure adds the ground's z0 to the distance it is given. Where a building's wall is nearer by d + z0, with
// its own z0, the sum is the wall's: 5.01 m in a side's wall cell with buildings.z0 = 0.01 m on ground of z0 = 1 m,
// not 6 m. Elsewhere the ground's distance keeps its bits. With a uniform z0 and with a map (1 m for x < 40 m,
// 0.05 m beyond).
TEST (BuildingTke, LengthTakesTheNearerWallsOwnRoughness)
{
    const double z0b = 0.01;
    for (const int mgs : {8, 3}) {
        for (const bool mapped : {false, true}) {
            Set s(mgs, tile_options("block45r.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }, 5, z0b));
            const auto z0g = [mapped] (int i) { return (mapped && i >= 4) ? 0.05 : 1.0; };
            amrex::MultiFab map;
            if (mapped) {
                const std::string file = write_text("z0map.xyz", [&] (std::ostream& os) {
                    for (int j = 0; j < 8; ++j) { for (int i = 0; i < 8; ++i) { os << 10 * i + 5 << " " << 10 * j + 5 << " " << z0g(i) << "\n"; } }
                });
                gust::read_roughness_map(file, "bc.surface_z0_file", 1, 2.0, s.c.geom, s.c.ba, s.c.dm, map);
            }
            const gust::Roughness rough = mapped ? gust::Roughness(amrex::Real(0), &map) : gust::Roughness(amrex::Real(1));
            amrex::MultiFab dist, before;
            gust::height_above_ground(s.c.z, dist);
            before.define(dist.boxArray(), dist.DistributionMap(), 1, 1);
            amrex::MultiFab::Copy(before, dist, 0, 0, 1, 1);
            const amrex::Long n = s.tke.limit_wall_distance(dist, rough);
            amrex::Long want = 0, kept_wrong = 0;
            double err = 0;
            for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
                const auto a = dist.const_array(mfi);
                const auto o = before.const_array(mfi);
                const auto db = s.tke.distance().const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    const double ground = static_cast<double>(o(i, j, k)) + z0g(i), wall = static_cast<double>(db(i, j, k)) + z0b;
                    if (wall < ground) {
                        ++want;
                        err = std::max(err, std::abs(static_cast<double>(a(i, j, k)) + z0g(i) - wall) / wall);
                    } else if (a(i, j, k) != o(i, j, k)) {
                        ++kept_wrong;
                    }
                });
            }
            amrex::ParallelDescriptor::ReduceLongSum(want);
            amrex::ParallelDescriptor::ReduceLongSum(kept_wrong);
            amrex::ParallelDescriptor::ReduceRealMax(err);
            EXPECT_EQ(kept_wrong, 0) << "cells nearer to the ground whose distance changed; boxes " << mgs << " map " << mapped;
            EXPECT_LT(err, 100 * rtol) << "d + z0 against the wall's d_b + z0_b";
            EXPECT_EQ(n, want);
            EXPECT_GT(n, 100);
            // A side's wall cell at x = 55 m, 35 m up, 5 m from the side: 5.01 m, the log law's kappa u* (d + z0_b).
            EXPECT_NEAR(value_at(dist, 5, 3, 3) + z0g(5), 5.01, 100 * rtol * 5);
            // The same on the low side, on the 1 m ground in both cases.
            EXPECT_NEAR(value_at(dist, 2, 3, 3) + 1.0, 5.01, 100 * rtol * 5);
            // The first cell beside the side at x = 25 m: the ground is 5 m below and the wall 5 m away, the wall
            // smoother, so the wall's 5.01 m, not the ground's 6 m. The stored distance is 4.01 m.
            EXPECT_NEAR(value_at(dist, 2, 3, 0), 4.01, 100 * rtol * 5);
            EXPECT_EQ(value_at(dist, 0, 0, 0), 5.0) << "far from the block: the ground's";
        }
    }
    // A ground rougher than the wall's distance and roughness together: the stored distance is below zero and the
    // sum the closure forms is still the wall's.
    Set s(4, tile_options("block45r.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }, 5, z0b));
    amrex::MultiFab dist;
    gust::height_above_ground(s.c.z, dist);
    s.tke.limit_wall_distance(dist, gust::Roughness(amrex::Real(8)));
    EXPECT_NEAR(value_at(dist, 5, 3, 3), 5.01 - 8, 100 * rtol * 8);
}

// A tower (90 m) on a podium (35 m, ending inside the second of the 20 m layers). The podium's cells beside the
// tower are three quarters solid and carry the tower's side above their roof: the side lies on their face. Drawn at
// the cell's solid fraction from that face, as before, it stood 7.5 m out over the podium up to 40 m, and the cell
// two columns from the tower in the layer above (50 m up, 15 m above the roof) was 12.5 m from it.
TEST (BuildingTke, DistanceBesideAStepBetweenTwoRoofsIsTheBoxes)
{
    const std::vector<std::array<double, 5>> boxes = {{30, 50, 30, 50, 90}, {50, 70, 30, 50, 35}, {30, 50, 10, 30, 35}};
    for (const int mgs : {8, 3}) {
        Set s(mgs, tile_options("podium.txt", [] (int i, int j) {
            if (i == 2 && j == 2) { return 90; }
            return ((i == 3 && j == 2) || (i == 2 && j == 1)) ? 35 : 0;
        }), false, 20);
        ASSERT_EQ(static_cast<int>(value_at(s.b.wall(), 5, 3, 1, gust::kWallMask)), gust::kWallRoof | gust::kWallXlo);
        ASSERT_EQ(value_at(s.b.beta(), 5, 3, 1), 0.75);
        const auto [err, n] = distance_error(s, boxes);
        EXPECT_GT(n, 1000);
        EXPECT_LT(err, 100 * rtol) << "largest relative difference from the boxes' distance; boxes " << mgs;
        // Above the podium: one column from the tower 5 m to its side, two columns out 15 m to it and to the roof.
        EXPECT_NEAR(value_at(s.tke.distance(), 5, 3, 2), 5.0, 100 * rtol * 5);
        EXPECT_NEAR(value_at(s.tke.distance(), 6, 3, 2), 15.0, 100 * rtol * 15);
        // Beside the podium's side at y = 50 m, 5 m below its roof, in the column next to the tower: 15 m to that
        // side. The cell across the face holds the podium and the tower's side, not one block, so the wall keeps
        // the face's blocked heights over the whole width; cut to that cell's walls' depths it would be a line at
        // the tower, 15.8 m away, as the roof's edge is.
        EXPECT_NEAR(value_at(s.tke.distance(), 5, 6, 1), 15.0, 100 * rtol * 15);
    }
}

// A block whose sides and top lie inside cells: 15 m square (22.5 .. 37.5 m, three quarters of the cells 2 and 3
// each way) and 45 m high. Its corner cells hold a quarter column, with walls along both horizontal directions,
// and in the top layer along all three. Each wall is drawn over the block's own face.
TEST (BuildingTke, DistanceToABlockWithEdgesInsideCellsIsTheBoxs)
{
    const std::vector<std::array<double, 5>> box = {{22.5, 37.5, 22.5, 37.5, 45}};
    for (const int mgs : {8, 3}) {
        Set s(mgs, tile_options("inside.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }, 6, 0.1, 15));
        ASSERT_EQ(value_at(s.b.beta(), 2, 2, 1), 0.5625) << "a corner cell: three quarters each way";
        ASSERT_EQ(value_at(s.b.beta(), 2, 2, 4), 0.28125) << "and half of that in the top layer";
        ASSERT_EQ(static_cast<int>(value_at(s.b.wall(), 2, 2, 4, gust::kWallMask)), gust::kWallRoof | gust::kWallXhi | gust::kWallYhi);
        const auto [err, n] = distance_error(s, box);
        EXPECT_GT(n, 1000);
        EXPECT_LT(err, 100 * rtol) << "largest relative difference from the box's distance; boxes " << mgs;
        // Off the corner in the first layers, and above the roof's corner.
        EXPECT_NEAR(value_at(s.tke.distance(), 1, 1, 1), std::hypot(7.5, 7.5), 100 * rtol * 11);
        EXPECT_NEAR(value_at(s.tke.distance(), 1, 1, 5), std::sqrt(2 * 7.5 * 7.5 + 100), 100 * rtol * 15);
        EXPECT_NEAR(value_at(s.tke.distance(), 4, 3, 5), std::hypot(7.5, 10.0), 100 * rtol * 13);
        // Two tiles along x (22.5 .. 52.5 m): the cells 3 and 4 are no corners. They hold a slice of the block along y,
        // with one side wall, and their roof covers the slice, not the cell: above and beside it, at (35, 15, 55) m,
        // the roof's edge is 7.5 m off in y and 10 m below.
        Set wide(mgs, tile_options("inside_wide.txt", [] (int i, int j) { return ((i == 2 || i == 3) && j == 2) ? 45 : 0; }, 6, 0.1, 15));
        ASSERT_EQ(static_cast<int>(value_at(wide.b.wall(), 3, 2, 4, gust::kWallMask)), gust::kWallRoof | gust::kWallYhi);
        const auto [werr, wn] = distance_error(wide, {{22.5, 52.5, 22.5, 37.5, 45}});
        EXPECT_GT(wn, 1000);
        EXPECT_LT(werr, 100 * rtol) << "the wide block; boxes " << mgs;
        EXPECT_NEAR(value_at(wide.tke.distance(), 3, 1, 5), std::hypot(7.5, 10.0), 100 * rtol * 13);
    }
}

// A block with two faces on cell faces (x = 30 and 50 m) and two inside cells (y = 22.5 and 37.5 m). The fluid
// cells beside the aligned faces see three quarters of their face blocked, along y, and those above the roof three
// quarters of their floor: each wall is the block's cross-section, read from the partly solid cell across the face.
// Drawn over the height's three quarters and the whole width, or over the whole floor, they stood 2.5 m out past
// the block's edges. With the top on a layer's face (40 m) and inside a layer (45 m).
TEST (BuildingTke, DistanceToABlockWithTwoFacesInsideCellsIsTheBoxs)
{
    for (const double top : {40.0, 45.0}) {
        const std::vector<std::array<double, 5>> box = {{30, 50, 22.5, 37.5, top}};
        gust::BuildingsOptions o;
        o.type = gust::BuildingsType::immersed_forcing;
        o.type_name = "immersed_forcing";
        o.small_volfrac = 0;
        o.z0 = 0.1;
        o.file = write_text(top > 40 ? "two_inside45.txt" : "two_inside40.txt", [&] (std::ostream& os) {
            os << 5 << "\n" << 6 << "\n";
            for (int p = 0; p < 5; ++p) { os << 20 * p << "\n"; }
            for (int p = 0; p < 6; ++p) { os << 15 * p << "\n"; }
            for (int i = 0; i < 5; ++i) { for (int j = 0; j < 6; ++j) { os << ((i == 2 && j == 2) ? top : 0.0) << "\n"; } }
        });
        for (const int mgs : {8, 3}) {
            Set s(mgs, o);
            ASSERT_EQ(value_at(s.b.beta(), 3, 2, 1), 0.75);
            ASSERT_EQ(value_at(s.b.beta_face(0), 5, 2, 1), 0.75) << "the face at x = 50 m, blocked along y";
            ASSERT_EQ(static_cast<int>(value_at(s.b.wall(), 5, 2, 1, gust::kWallMask)), gust::kWallXlo);
            const auto [err, n] = distance_error(s, box);
            EXPECT_GT(n, 1000);
            EXPECT_LT(err, 100 * rtol) << "largest relative difference from the box's distance; top " << top << " boxes " << mgs;
            // Beside the corner at (50, 22.5) m: 5 m out in x and 7.5 m in y. And above it.
            EXPECT_NEAR(value_at(s.tke.distance(), 5, 1, 1), std::hypot(5.0, 7.5), 100 * rtol * 9);
            EXPECT_NEAR(value_at(s.tke.distance(), 5, 1, 6), std::sqrt(25 + 7.5 * 7.5 + (65 - top) * (65 - top)), 100 * rtol * 30);
            EXPECT_NEAR(value_at(s.tke.distance(), 3, 1, 6), std::hypot(7.5, 65 - top), 100 * rtol * 30);
        }
    }
}

// The scaled norm of R_k counts the solid cells as zeros: their rows k - k_floor are no part of the flow. With
// R_k = 3 everywhere the norm is 3 sqrt(1 - share of solid cells) over the scale, not 3; the other norms and the
// fields of R are as they were.
TEST (BuildingTke, NormOfKLeavesTheSolidCellsOut)
{
    for (const int mgs : {8, 3}) {
        Set s(mgs, tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }));
        gust::ResidualFields R;
        R.define(s.c.ba, s.c.dm);
        for (amrex::MultiFab* f : {&R.u, &R.v, &R.w, &R.p, &R.theta}) { f->setVal(2); }
        R.tke.setVal(3);
        gust::State q;
        q.define(s.c.ba, s.c.dm);
        for (amrex::MultiFab* f : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { f->setVal(1.0); }
        const gust::JvpScales sc;
        const double sk = static_cast<double>(sc.velocity * sc.tke);   // U k_s / H with H = 1
        std::vector<amrex::Real> all, fluid, mon;
        amrex::MultiFab work;
        gust::stage_norms(R, q, nullptr, sc, 1, s.c.geom, all, mon);
        gust::stage_norms(R, q, nullptr, sc, 1, s.c.geom, fluid, mon, &s.tke, &work);
        // The block: 2 x 2 columns, 4 solid layers (the fifth is half solid).
        const double solid = 16.0 / (8 * 8 * 24);
        EXPECT_NEAR(static_cast<double>(all[5]), 3 / sk, 100 * rtol * 3 / sk);
        EXPECT_NEAR(static_cast<double>(fluid[5]), 3 * std::sqrt(1 - solid) / sk, 100 * rtol * 3 / sk) << "boxes " << mgs;
        for (int n = 0; n < 5; ++n) { EXPECT_EQ(fluid[static_cast<std::size_t>(n)], all[static_cast<std::size_t>(n)]) << n; }
        EXPECT_EQ(value_at(work, 3, 3, 1), 0.0) << "a solid cell";
        EXPECT_EQ(value_at(work, 3, 3, 4), 3.0) << "the half solid cell under the roof stays";
        EXPECT_EQ(value_at(R.tke, 3, 3, 1), 3.0) << "R_k itself is not changed: the sweep reads it";
    }
}

// The same in the count of steady.confidence (SimplecStage::confidence): after one iteration from k = 1 m2/s2 in
// every cell, the residual the stage holds has k - k_floor, about 1, in the solid cells and at most the wall
// cells' relaxation elsewhere. With a threshold between the two every cell meets the criterion.
TEST (BuildingTke, ConfidenceLeavesTheSolidCellsOut)
{
    Set s(4, tile_options("block45p.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }, 4), true);
    s.c.mesh.set_open_faces(s.b.beta_face(0), s.b.beta_face(1), s.b.beta_face(2), s.b.beta(), s.c.geom);
    gust::State q;
    q.define(s.c.ba, s.c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, s.c.z, 0, 0, 1, s.c.z.nGrow());
    q.u.setVal(3);
    q.v.setVal(0);
    q.w.setVal(0);
    q.p.setVal(0);
    q.theta.setVal(300);
    q.tke.setVal(1);
    amrex::MultiFab nu(s.c.ba, s.c.dm, 1, 1);
    nu.setVal(0.5);
    using gust::FaceType;
    gust::BoundaryInputs bin;
    bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall, FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
    gust::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.theta_ref.z = {-1000, 5000};
    in.theta_ref.theta = {300, 300};
    in.buoyancy.gravity = 0;
    in.buildings = &s.law;
    in.building_tke = &s.tke;
    const gust::BoundaryConditions bc(bin, s.c.geom, s.c.mesh, in.coriolis);
    in.bc = &bc;
    gust::PoissonOptions popt;
    popt.reltol = amrex::Real(0.1);
    gust::SimplecOptions opt;
    gust::SimplecStage stage(opt, 2, in, nullptr, popt, s.c.mesh, s.c.geom, s.c.z, amrex::Real(kFloor), s.c.ba, s.c.dm);
    stage.init(q);
    bc.fill_ghosts(q, nullptr, in.buoyancy);
    amrex::MultiFab dtau;
    gust::local_dtau(q, &stage.nu(), s.c.mesh, s.c.geom, amrex::Real(1), dtau);
    stage.iterate(q, dtau);
    // The residual of the first state: the solid cell's row, and the largest of the others.
    const amrex::MultiFab& Rk = stage.residual_fields().tke;
    EXPECT_NEAR(value_at(Rk, 3, 3, 1), 1.0 - kFloor, 100 * rtol);
    amrex::MultiFab fluid(s.c.ba, s.c.dm, 1, 0);
    amrex::MultiFab::Copy(fluid, Rk, 0, 0, 1, 0);
    s.tke.zero_solid_cells(fluid);
    const double largest = static_cast<double>(fluid.norm0());
    ASSERT_LT(largest, 0.4) << "the wall cells' rows are well below the solid cells'";
    // |R_k| <= 0.5 m2/s3 as a scaled threshold; the other fields always pass.
    const gust::JvpScales sc;
    const double sk = static_cast<double>(sc.velocity * sc.tke) / static_cast<double>(stage.norm_height());
    const std::vector<amrex::Real> threshold = {amrex::Real(1e20), amrex::Real(1e20), amrex::Real(1e20), amrex::Real(1e20), amrex::Real(1e20),
                                                static_cast<amrex::Real>(0.5 / sk)};
    amrex::Real share = 0, height = -1;
    ASSERT_TRUE(stage.confidence(threshold, share, height));
    EXPECT_EQ(share, amrex::Real(1)) << "the 16 solid cells, with |R_k| near 1, do not count against the criterion";
    EXPECT_EQ(height, amrex::Real(0));
}

// The rows through residual() without a surface layer: the wall's row takes the place of the cell's equation in
// the first layer too (with a surface layer it is added to the ground's row there, or to the equation).
TEST (BuildingTke, ResidualWithoutASurfaceLayerPutsTheWallsRowInPlace)
{
    Set s(3, tile_options("block25.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 25 : 0; }));
    s.c.mesh.set_open_faces(s.b.beta_face(0), s.b.beta_face(1), s.b.beta_face(2), s.b.beta(), s.c.geom);
    gust::State q;
    q.define(s.c.ba, s.c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, s.c.z, 0, 0, 1, s.c.z.nGrow());
    const double U = 3;
    q.u.setVal(amrex::Real(U));
    q.v.setVal(0);
    q.w.setVal(0);
    q.p.setVal(0);
    q.theta.setVal(300);
    const auto kval = [] (int i, int j, int kk) { return 0.7 + 0.01 * i + 0.02 * j + 0.03 * kk; };
    for (amrex::MFIter mfi(q.tke, false); mfi.isValid(); ++mfi) {
        const auto a4 = q.tke.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int kk) { a4(i, j, kk) = static_cast<amrex::Real>(kval(i, j, kk)); });
    }
    amrex::MultiFab nu(s.c.ba, s.c.dm, 1, 1);
    nu.setVal(0.5);
    gust::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.theta_ref.z = {0};
    in.theta_ref.theta = {300};
    in.buildings = &s.law;
    in.building_tke = &s.tke;
    gust::ResidualFields R;
    gust::residual(q, in, s.c.mesh, s.c.geom, R);
    // The side at y = 50 m with u along it: its own u*, 5 m away; the face is wholly blocked in the first two layers.
    const double ss = kKappa * U / std::log(50.0), ks = ss * ss / (kCmu0 * kCmu0), as = ss / 5.0;
    EXPECT_NEAR(value_at(R.tke, 3, 5, 0), as * (kval(3, 5, 0) - ks), 1e3 * rtol * as) << "first layer: in the equation's place";
    EXPECT_NEAR(value_at(R.tke, 3, 5, 1), as * (kval(3, 5, 1) - ks), 1e3 * rtol * as);
    EXPECT_NEAR(value_at(R.tke, 3, 3, 0), kval(3, 3, 0) - kFloor, 100 * rtol) << "a solid cell: k - floor";
    EXPECT_NEAR(value_at(R.tke, 3, 3, 1), kval(3, 3, 1) - kFloor, 100 * rtol);
}

// The zone of the lee options (buildings.limiter_zone = c): the weight is 0 within c H_b / 2 of a building, 1 beyond
// c H_b, the smooth step t^2 (3 - 2 t) between, t = 2 d_b / (c H_b) - 1, with H_b the height of the building whose
// wall is the nearest. Two blocks, 45 m and 25 m high, 20 m apart (tiles 1 and 3: cells 1-2 and 5-6 along x).
TEST (BuildingTke, LimiterZoneFollowsTheDistanceAndTheBuildingsHeight)
{
    for (const int mgs : {8, 3}) {
        gust::BuildingsOptions o = tile_options("two_blocks.txt", [] (int i, int j) {
            if (j != 2) { return 0; }
            return i == 1 ? 45 : i == 3 ? 25 : 0;
        });
        o.limiter_zone = 2;
        o.limiter_zone_acts = true;   // as the application sets it with turb.realizability > 0
        Set s(mgs, o);
        ASSERT_TRUE(s.tke.zone_active());
        // H_b: beside the tall block 45 m, beside the short one 25 m, above both the tall one's.
        EXPECT_NEAR(value_at(s.tke.building_height(), 1, 1, 1), 45.0, 100 * rtol * 45) << "boxes " << mgs;
        EXPECT_NEAR(value_at(s.tke.building_height(), 6, 1, 1), 25.0, 100 * rtol * 25);
        EXPECT_NEAR(value_at(s.tke.building_height(), 7, 3, 1), 25.0, 100 * rtol * 25);
        EXPECT_NEAR(value_at(s.tke.building_height(), 1, 3, 10), 45.0, 100 * rtol * 45) << "above the tall roof";
        // The weight against the closed form at every fluid cell with a wall found.
        double err = 0;
        amrex::Long zero = 0, part = 0, one = 0;
        for (amrex::MFIter mfi(s.tke.limiter_zone(), false); mfi.isValid(); ++mfi) {
            const auto w = s.tke.limiter_zone().const_array(mfi);
            const auto d = s.tke.distance().const_array(mfi);
            const auto h = s.tke.building_height().const_array(mfi);
            const auto beta = s.b.beta().const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (beta(i, j, k) >= 1) { return; }
                double want = 1;
                if (d(i, j, k) < gust::kBuildingFar) {
                    const double t = 2 * static_cast<double>(d(i, j, k)) / (2 * static_cast<double>(h(i, j, k))) - 1;
                    want = t <= 0 ? 0 : t >= 1 ? 1 : t * t * (3 - 2 * t);
                }
                err = std::max(err, std::abs(static_cast<double>(w(i, j, k)) - want));
                // The classes from the field itself (the closed form is checked through err).
                const double got = static_cast<double>(w(i, j, k));
                if (got == 0) { ++zero; } else if (got < 1) { ++part; } else { ++one; }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(err);
        amrex::ParallelDescriptor::ReduceLongSum(zero);
        amrex::ParallelDescriptor::ReduceLongSum(part);
        amrex::ParallelDescriptor::ReduceLongSum(one);
        EXPECT_LT(err, 100 * rtol);
        EXPECT_GT(zero, 100);
        EXPECT_GT(part, 100);
        EXPECT_GT(one, 100);
        EXPECT_EQ(s.tke.zone_off_cells(), zero);
        EXPECT_EQ(s.tke.zone_cells(), zero + part);
        // By hand: 5 m from the tall block's side (cell 0, beside x = 10 m) and 5 m from the short one's (cell 7,
        // beside x = 70 m): inside both zones.
        EXPECT_EQ(value_at(s.tke.limiter_zone(), 0, 3, 1), 0.0);
        EXPECT_EQ(value_at(s.tke.limiter_zone(), 7, 3, 1), 0.0) << "5 m from the short block: within 25 m";
        // 30 m above the short roof at x = 65 m (cell 6, layer 5): the tall block's top edge is 36 m away, so the
        // short roof is the nearest wall, H_b = 25 m, t = 2 * 30 / 50 - 1. At x = 55 m the tall block's edge, 27 m
        // away, is the nearer, with H_b = 45 m: inside.
        const double t = 2 * 30.0 / 50 - 1;
        EXPECT_NEAR(value_at(s.tke.limiter_zone(), 6, 3, 5), t * t * (3 - 2 * t), 100 * rtol);
        EXPECT_NEAR(value_at(s.tke.building_height(), 5, 3, 5), 45.0, 100 * rtol * 45);
        EXPECT_EQ(value_at(s.tke.limiter_zone(), 5, 3, 5), 0.0);
        // Far above both: 1.
        EXPECT_EQ(value_at(s.tke.limiter_zone(), 6, 3, 20), 1.0);
    }
    // c = 0, or a zone that does not act (no turb.realizability): none.
    gust::BuildingsOptions off = tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; });
    off.limiter_zone = 0;
    off.limiter_zone_acts = true;
    Set s0(4, off);
    EXPECT_FALSE(s0.tke.zone_active());
    off.limiter_zone = 2;
    off.limiter_zone_acts = false;
    Set s1(4, off);
    EXPECT_FALSE(s1.tke.zone_active());
}

// The restart hash gains a line for the zone only where it acts: without it the hash is the one from before the zone,
// whatever buildings.limiter_zone holds (OptionDefaultsAndRefusals pins those hashes for this file), and with it two
// sizes give two hashes.
TEST (BuildingTke, LimiterZoneEntersTheRestartHashWhereItActs)
{
    gust::BuildingsOptions o = tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; });
    o.limiter_zone_acts = false;
    const std::string before = gust::Buildings::hash_inputs(o);
    o.limiter_zone = 7;
    EXPECT_EQ(gust::Buildings::hash_inputs(o), before) << "a zone that does not act is not in the hash";
    o.limiter_zone_acts = true;
    const std::string seven = gust::Buildings::hash_inputs(o);
    EXPECT_NE(seven, before);
    o.limiter_zone = 3;
    EXPECT_NE(gust::Buildings::hash_inputs(o), seven);
    EXPECT_NE(gust::Buildings::hash_inputs(o), before);
}

// The zone's lateral ghost ring, which the cap reads, equals its neighbours' valid cells: across a periodic seam (the
// block straddles it) and between boxes; outside a non-periodic face the ring keeps the nearest column's value, so
// the comparison after a periodic fill is exact there too.
TEST (BuildingTke, LimiterZoneGhostRingMatchesTheNeighbours)
{
    for (const bool periodic : {true, false}) {
        gust::BuildingsOptions o = tile_options(periodic ? "seam_zone.txt" : "edge_zone.txt",
                                                [] (int i, int j) { return (i == 0 && j == 2) ? 45 : 0; }, periodic ? 4 : 5);
        o.limiter_zone = 2;
        o.limiter_zone_acts = true;
        Set s(3, o, periodic);
        const amrex::MultiFab& z = s.tke.limiter_zone();
        amrex::MultiFab tmp(z.boxArray(), z.DistributionMap(), 1, 1);
        amrex::MultiFab::Copy(tmp, z, 0, 0, 1, 1);
        tmp.FillBoundary(s.c.geom.periodicity());
        EXPECT_EQ(gust::count_bit_differences(tmp, z, amrex::IntVect(1, 1, 0)), 0) << (periodic ? "periodic" : "not periodic");
        EXPECT_GT(s.tke.zone_off_cells(), 0);
        if (periodic) {
            // The block (tile -10 .. 10 m) fills the cells 7 and 0 across the seam: the cells beside it on either side, and
            // two cells out, hold the same weight.
            for (const int k : {1, 4, 6}) {
                EXPECT_EQ(value_at(z, 6, 3, k), value_at(z, 1, 3, k)) << "beside the block, layer " << k;
                EXPECT_EQ(value_at(z, 5, 3, k), value_at(z, 2, 3, k)) << "two cells out, layer " << k;
            }
        }
    }
}

// The cap's weight (Closure::set_cap_weight): the cap divided by w, the largest Real at w = 0, bit for bit the same
// at w = 1 and without a weight.
TEST (BuildingTke, CapWeightLiftsTheCap)
{
    Set s(3, tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }));
    gust::State q;
    q.define(s.c.ba, s.c.dm);
    for (amrex::MFIter mfi(q.tke, false); mfi.isValid(); ++mfi) {
        const auto a4 = q.tke.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int kk) { a4(i, j, kk) = static_cast<amrex::Real>(0.5 + 0.01 * (i + j + kk)); });
    }
    amrex::MultiFab S2(s.c.ba, s.c.dm, 1, 1);
    S2.setVal(amrex::Real(0.04));
    gust::ClosureOptions copt;
    copt.realizability = amrex::Real(0.309);
    const gust::ClosureConstants cc = gust::closure_constants(copt);
    amrex::MultiFab plain, ones, weighted, w(s.c.ba, s.c.dm, 1, 1);
    gust::nut_cap_field(q, S2, copt, cc, plain);
    w.setVal(1);
    gust::nut_cap_field(q, S2, copt, cc, ones, &w);
    EXPECT_EQ(gust::count_bit_differences(plain, ones, amrex::IntVect(1, 1, 0)), 0) << "w = 1 is the cap bit for bit";
    for (amrex::MFIter mfi(w, false); mfi.isValid(); ++mfi) {
        const auto a4 = w.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int kk) { a4(i, j, kk) = (i < 2) ? amrex::Real(0) : (i < 4 ? amrex::Real(0.25) : amrex::Real(1)); });
    }
    gust::nut_cap_field(q, S2, copt, cc, weighted, &w);
    const double big = static_cast<double>(std::numeric_limits<amrex::Real>::max());
    EXPECT_EQ(value_at(weighted, 1, 3, 2), big) << "w = 0: no cap";
    EXPECT_NEAR(value_at(weighted, 3, 3, 2), 4 * value_at(plain, 3, 3, 2), 100 * rtol * value_at(plain, 3, 3, 2)) << "w = 1/4";
    EXPECT_EQ(value_at(weighted, 5, 3, 2), value_at(plain, 5, 3, 2)) << "w = 1";
}

// UT-5.4.1, UT-5.4.3: the rows.
TEST (BuildingTke, WallCellsRelaxToTheLawsValueAndSolidCellsToTheFloor)
{
    for (const int mgs : {8, 4, 3}) {
        // A roof at 25 m: the cell (3, 3, 2) is half solid, 2.5 m from the roof, its reference 10 m.
        Set s(mgs, tile_options("block25.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 25 : 0; }));
        amrex::MultiFab u(s.c.mesh.ax.boxArray(), s.c.dm, 1, 1), v(s.c.mesh.ay.boxArray(), s.c.dm, 1, 1), w(s.c.mesh.az.boxArray(), s.c.dm, 1, 1);
        const double U = 3;
        u.setVal(amrex::Real(U));
        v.setVal(0);
        w.setVal(0);
        s.law.update(u, v, w);
        const double so = kKappa * U / std::log(25.0), sr = kKappa * U / std::log(100.0);
        const double us2 = 0.5 * so * so + 0.5 * sr * sr, kw = us2 / (kCmu0 * kCmu0), a = std::sqrt(us2) / 2.5;
        EXPECT_NEAR(value_at(s.law.rate(), 3, 3, 2, gust::kBuildingUstar2), us2, 100 * rtol * us2);

        amrex::MultiFab k(s.c.ba, s.c.dm, 1, 1), Rk(s.c.ba, s.c.dm, 1, 0);
        // k differs from cell to cell, so a row reading a neighbour's k fails.
        const auto kval = [] (int i, int j, int kk) { return 0.7 + 0.01 * i + 0.02 * j + 0.03 * kk; };
        const auto fill_k = [&] () {
            for (amrex::MFIter mfi(k, false); mfi.isValid(); ++mfi) {
                const auto a4 = k.array(mfi);
                amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int kk) { a4(i, j, kk) = static_cast<amrex::Real>(kval(i, j, kk)); });
            }
        };
        fill_k();
        Rk.setVal(9);
        s.tke.apply_rows(k, Rk, gust::BuildingKRow::replace);
        EXPECT_NEAR(value_at(Rk, 3, 3, 2), a * (kval(3, 3, 2) - kw), 1e3 * rtol * a) << "boxes " << mgs;
        EXPECT_NEAR(value_at(Rk, 3, 3, 1), kval(3, 3, 1) - kFloor, 100 * rtol) << "a solid cell: k - floor";
        EXPECT_EQ(value_at(Rk, 3, 3, 4), 9.0) << "no wall cell: its equation stays";
        EXPECT_EQ(value_at(Rk, 0, 0, 0), 9.0);
        // A side's wall cell with u along the side at y = 50 m: its own u* (beta 0), 5 m from the side.
        const double ss = kKappa * U / std::log(50.0), ks = ss * ss / (kCmu0 * kCmu0), as = ss / 5.0;
        EXPECT_NEAR(value_at(Rk, 3, 5, 1), as * (kval(3, 5, 1) - ks), 1e3 * rtol * as);
        // The same side in the third layer, where the block's top at 25 m leaves half the face blocked: half
        // the cell's own row (9) stays.
        ASSERT_EQ(value_at(s.b.beta_face(1), 3, 5, 2), 0.5);
        EXPECT_NEAR(value_at(Rk, 3, 5, 2), 0.5 * 9 + 0.5 * as * (kval(3, 5, 2) - ks), 1e3 * rtol * 9);
        // A side's wall cell with the flow normal to its side (x = 50 m): u* = 0, so the floor's value at the lowest rate.
        EXPECT_NEAR(value_at(Rk, 5, 3, 1), 0.01 / 5.0 * (kval(5, 3, 1) - kFloor), 1e3 * rtol * 0.002);
        // At the wall value the row is zero.
        k.setVal(amrex::Real(kw));
        s.tke.apply_rows(k, Rk, gust::BuildingKRow::replace);
        EXPECT_NEAR(value_at(Rk, 3, 3, 2), 0.0, 1e3 * rtol * a * kw);

        // The ground's first cell beside the side: added to the ground's row, or in its place.
        fill_k();
        Rk.setVal(9);
        s.tke.apply_rows(k, Rk, gust::BuildingKRow::ground_row);
        EXPECT_NEAR(value_at(Rk, 3, 5, 0), 9 + as * (kval(3, 5, 0) - ks), 1e3 * rtol * 9);
        EXPECT_NEAR(value_at(Rk, 3, 5, 1), as * (kval(3, 5, 1) - ks), 1e3 * rtol * as) << "above the first cell: in its place";
        EXPECT_NEAR(value_at(Rk, 3, 3, 0), kval(3, 3, 0) - kFloor, 100 * rtol) << "a solid first cell: the floor's row alone";
        // bc.surface_wall_k = source: added to the cell's row everywhere, by the blocked part.
        Rk.setVal(9);
        s.tke.apply_rows(k, Rk, gust::BuildingKRow::source);
        EXPECT_NEAR(value_at(Rk, 3, 5, 1), 9 + as * (kval(3, 5, 1) - ks), 1e3 * rtol * 9);
        EXPECT_NEAR(value_at(Rk, 3, 5, 2), 9 + 0.5 * as * (kval(3, 5, 2) - ks), 1e3 * rtol * 9);
        EXPECT_NEAR(value_at(Rk, 3, 3, 1), kval(3, 3, 1) - kFloor, 100 * rtol) << "solid cells hold the floor in every form";
        EXPECT_EQ(value_at(Rk, 0, 0, 0), 9.0);

        amrex::MultiFab A(s.c.ba, s.c.dm, gust::n_implicit_coefs, 0);
        A.setVal(2);
        s.tke.add_to_implicit(s.open(), A, gust::BuildingKRow::ground_row);
        EXPECT_NEAR(value_at(A, 3, 3, 2, gust::coef_P), a, 1e3 * rtol * a);
        EXPECT_NEAR(value_at(A, 3, 5, 2, gust::coef_P), 0.5 * 2 + 0.5 * as, 1e3 * rtol) << "half the face blocked: half of each";
        EXPECT_EQ(value_at(A, 3, 5, 2, gust::coef_ZP), 1.0) << "half of the cell's own coupling stays";
        EXPECT_NEAR(value_at(A, 3, 5, 1, gust::coef_P), as, 1e3 * rtol * as);
        EXPECT_NEAR(value_at(A, 3, 5, 0, gust::coef_P) - 2, as, 1e3 * rtol * as) << "added to the ground's row";
        EXPECT_EQ(value_at(A, 3, 5, 0, gust::coef_ZP), 2.0) << "the ground's row keeps what it had";
        EXPECT_EQ(value_at(A, 3, 3, 1, gust::coef_P), 1.0) << "a solid cell: the identity";
        for (const int n : {gust::coef_XM, gust::coef_XP, gust::coef_YM, gust::coef_YP, gust::coef_ZM, gust::coef_ZP}) {
            EXPECT_EQ(value_at(A, 3, 3, 2, n), 0.0);
            EXPECT_EQ(value_at(A, 3, 3, 1, n), 0.0);
        }
        // Next to no building: the couplings through open faces stay, the one through a blocked face goes.
        EXPECT_EQ(value_at(A, 0, 0, 0, gust::coef_P), 2.0);
        EXPECT_EQ(value_at(A, 3, 5, 0, gust::coef_YM), 0.0) << "the face on the side is blocked";
        EXPECT_EQ(value_at(A, 3, 5, 0, gust::coef_YP), 2.0);
        A.setVal(2);
        s.tke.add_to_implicit(s.open(), A, gust::BuildingKRow::replace);
        EXPECT_NEAR(value_at(A, 3, 5, 0, gust::coef_P), as, 1e3 * rtol * as) << "no ground row: in its place";
        EXPECT_EQ(value_at(A, 3, 5, 0, gust::coef_ZP), 0.0);
        A.setVal(2);
        s.tke.add_to_implicit(s.open(), A, gust::BuildingKRow::source);
        EXPECT_NEAR(value_at(A, 3, 5, 1, gust::coef_P) - 2, as, 1e3 * rtol * as) << "source: added to the cell's row";
        EXPECT_EQ(value_at(A, 3, 5, 1, gust::coef_ZP), 2.0);
        EXPECT_EQ(value_at(A, 3, 3, 1, gust::coef_P), 1.0) << "a solid cell: the identity in every form";
    }
}

// k does not diffuse through a blocked face; half a face lets half through.
TEST (BuildingTke, DiffusiveFluxOfKIsClosedOnBlockedFaces)
{
    for (const int mgs : {8, 3}) {
        Set s(mgs, tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; }));
        amrex::Array<amrex::MultiFab, 3> G;
        gust::define_cv_fluxes(gust::Loc::cell, s.c.ba, s.c.dm, G);
        for (auto& g : G) { g.setVal(2); }
        s.tke.close_fluxes(s.open(), G);
        EXPECT_EQ(value_at(G[0], 5, 3, 1), 0.0) << "the side at x = 50 m";
        EXPECT_EQ(value_at(G[0], 3, 3, 1), 0.0) << "the side at x = 30 m";
        EXPECT_EQ(value_at(G[0], 5, 3, 4), 1.0) << "half blocked in the fifth layer";
        EXPECT_EQ(value_at(G[0], 6, 3, 1), 2.0);
        EXPECT_EQ(value_at(G[1], 3, 5, 1), 0.0) << "the side at y = 50 m";
        EXPECT_EQ(value_at(G[1], 3, 6, 1), 2.0);
        EXPECT_EQ(value_at(G[2], 3, 3, 4), 0.0) << "the face under the roof's cell";
        EXPECT_EQ(value_at(G[2], 3, 3, 5), 2.0);
        amrex::MultiFab A(s.c.ba, s.c.dm, gust::n_implicit_coefs, 0);
        A.setVal(2);
        s.tke.add_to_implicit(s.open(), A, gust::BuildingKRow::replace);
        // A fluid cell that is no wall cell, beside the half-blocked face's cell: its couplings stay.
        EXPECT_EQ(value_at(A, 6, 3, 4, gust::coef_XM), 2.0);
        EXPECT_EQ(value_at(A, 6, 3, 4, gust::coef_P), 2.0);
        // Above the half-solid roof cell: the face between them is open.
        EXPECT_EQ(value_at(A, 3, 3, 5, gust::coef_ZM), 2.0);
    }
}

// An edge cell (a roof below and a side at low x) takes the larger u*^2 of its two walls, and a corner cell
// with a side 5 m away along x and one 2.5 m away along y (cells of 10 by 5 m) the smaller distance.
TEST (BuildingTke, SeveralWallsGiveTheLargestUstarAndTheNearestDistance)
{
    {
        Set s(3, tile_options("ell.txt", [] (int i, int j) {
            if ((i == 2 && j == 2) || (i == 1 && j == 2) || (i == 2 && j == 1)) { return 45; }
            return (i == 3 && j == 2) ? 20 : 0;
        }));
        ASSERT_EQ(static_cast<int>(value_at(s.b.wall(), 5, 3, 2, gust::kWallMask)), gust::kWallRoof | gust::kWallXlo);
        amrex::MultiFab u(s.c.mesh.ax.boxArray(), s.c.dm, 1, 1), v(s.c.mesh.ay.boxArray(), s.c.dm, 1, 1), w(s.c.mesh.az.boxArray(), s.c.dm, 1, 1);
        const auto us2_of = [] (double sp) { const double us = kKappa * sp / std::log(50.0); return us * us; };
        // The roof sees (u, v), the side (v, w): first the roof's speed is the larger, then the side's.
        for (const auto& uvw : {std::array<double, 3>{3, 2, 1}, std::array<double, 3>{1, 2, 3}}) {
            u.setVal(amrex::Real(uvw[0]));
            v.setVal(amrex::Real(uvw[1]));
            w.setVal(amrex::Real(uvw[2]));
            s.law.update(u, v, w);
            const double roof = us2_of(std::hypot(uvw[0], uvw[1])), side = us2_of(std::hypot(uvw[1], uvw[2]));
            EXPECT_NEAR(value_at(s.law.rate(), 5, 3, 2, gust::kBuildingUstar2), std::max(roof, side), 100 * rtol * roof);
            EXPECT_GT(std::abs(roof - side), 0.1 * roof);
        }
    }
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(8, 16, 24);
    c.in.mesh.prob_hi = {80, 80, 240};
    c.in.mesh.is_periodic = {0, 0, 0};
    c.in.mesh.max_grid_size = {4, 4};
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    const gust::BuildingsOptions o = tile_options("ell2.txt", [] (int i, int j) {
        return ((i == 2 && j == 2) || (i == 1 && j == 2) || (i == 2 && j == 1)) ? 45 : 0;
    });
    gust::Buildings b;
    b.build(o, c.geom, c.z, c.ba, c.dm);
    gust::BuildingForcing law;
    law.build(b, o, c.geom, c.mesh.detJ);
    gust::BuildingTke tke;
    tke.build(b, law, o, c.geom, c.z, amrex::Real(kFloor), amrex::Real(kCmu0));
    // The inner corner of the L: x 20..30 m (cell 2), y 25..30 m (cell 5).
    ASSERT_EQ(static_cast<int>(value_at(b.wall(), 2, 5, 1, gust::kWallMask)), gust::kWallXhi | gust::kWallYhi);
    EXPECT_NEAR(value_at(tke.distance(), 2, 5, 1), 2.5, 100 * rtol);
}

// The options: on by default with the law on every wall, off with the other forms; on needs that law.
TEST (BuildingTke, OptionDefaultsAndRefusals)
{
    gust::BuildingsOptions o = tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; });
    EXPECT_TRUE(gust::check_buildings_options(o).empty());
    // The pins are from before the open faces' mean of the wall law's speed: buildings.wall_speed = all_faces.
    o.wall_speed = gust::BuildingWallSpeed::all_faces;
    const std::string on = gust::Buildings::hash_inputs(o);
    o.wall_tke = false;
    const std::string off = gust::Buildings::hash_inputs(o);
    EXPECT_NE(on, off);
    // The hash without the closure's terms is the one from before them (pinned from main at 8f1b605 for this file),
    // and the one with them that of main at c57ad06.
    EXPECT_EQ(off, "e9cc23b2f4423a16596f835888fe3b16683c322fda07bbaae2be64ad4ad1f8a7");
    EXPECT_EQ(on, "a2c54457eb7d18fe03ee0c3ed4aad70a4fbcf06c6c9dcdec9cee245a0feb230a");
    o.wall_speed = gust::BuildingWallSpeed::open_faces;
    EXPECT_NE(gust::Buildings::hash_inputs(o), off);
    o.wall_tke = true;
    EXPECT_NE(gust::Buildings::hash_inputs(o), on);
    o.wall_tke = true;
    o.wall_form = gust::BuildingWallForm::roof_stress;
    o.wall_form_name = "roof_stress";
    const auto err = gust::check_buildings_options(o);
    ASSERT_EQ(err.size(), 1U);
    EXPECT_NE(err[0].find("buildings.wall_tke = on needs buildings.wall_form = fraction_stress"), std::string::npos);
    o.wall_tke_name = "maybe";
    const auto err2 = gust::check_buildings_options(o);
    ASSERT_EQ(err2.size(), 1U);
    EXPECT_NE(err2[0].find("buildings.wall_tke must be on or off (got maybe)"), std::string::npos);
}

} // namespace
