// The buildings' terms of the closure (BuildingTke.H), beyond the plan, sub-phase 5.4 (UT-5.4.x):
//   a wall cell's k row is a (k - u*^2 / Cmu0^2) with the law's u*, zero at the wall value, in the
//   residual and in the implicit rows, and it is added to the ground's row in the first cells
//   (UT-5.4.1); the closure's distance in wall cells is the fluid centroid's, (1 - beta) Delta_n / 2,
//   at beta 0, 0.5 and 0.9, floored at buildings.z0 (UT-5.4.2); solid cells hold k at its floor
//   (UT-5.4.3); where the ground is nearer than any building the wall distance keeps its bits
//   (UT-5.4.4); the distance to a block equals the closed form of a box, on every box layout and
//   across a periodic seam, and that to an L with a step the smallest over its boxes; the diffusive
//   flux of k and its couplings are closed on blocked faces.
#include "BuildingTke.H"
#include "Closure.H"
#include "PoissonTestHelpers.H"

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

// 8 x 8 x 24 cells of 10 m.
PCase mesh_case (int mgs, bool periodic = false)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(8, 8, 24);
    c.in.mesh.prob_hi = {80, 80, 240};
    c.in.mesh.is_periodic = {periodic ? 1 : 0, periodic ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

// A map in ERF's format on n x n points 20 m apart from 0: tile (ix, iy) covers 20 ix - 10 .. 20 ix + 10
// along x, so tile 2 is the cells 3 and 4. height(ix, iy) in m.
gust::BuildingsOptions tile_options (const std::string& name, const std::function<double(int, int)>& height, int n = 5, double z0 = 0.1)
{
    gust::BuildingsOptions o;
    o.type = gust::BuildingsType::immersed_forcing;
    o.type_name = "immersed_forcing";
    o.small_volfrac = 0;
    o.z0 = z0;
    o.file = write_text(name, [&] (std::ostream& os) {
        os << n << "\n" << n << "\n";
        for (int r = 0; r < 2; ++r) { for (int p = 0; p < n; ++p) { os << 20 * p << "\n"; } }
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
    Set (int mgs, const gust::BuildingsOptions& opt, bool periodic = false) : c(mesh_case(mgs, periodic)), o(opt)
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
            const double want = box_distance(boxes, 10 * i + 5, 10 * j + 5, 10 * k + 5, periodic);
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
        const amrex::Long n = s.tke.limit_wall_distance(dist);
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
    const std::string on = gust::Buildings::hash_inputs(o);
    o.wall_tke = false;
    const std::string off = gust::Buildings::hash_inputs(o);
    EXPECT_NE(on, off);
    // The hash without the closure's terms is the one from before them (pinned from main at 8f1b605 for this file).
    EXPECT_EQ(off, "e9cc23b2f4423a16596f835888fe3b16683c322fda07bbaae2be64ad4ad1f8a7");
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
