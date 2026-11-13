// PR-49 (plan rev. 6, sub-phase 6.4), UT-6.4.2 and UT-6.4.3, beyond the plan: the force on buildings
// (BuildingLoads.H), summed from the momentum residual on the blocked faces.
//   UT-6.4.2: air at rest in a linear pressure field p = G.x around a block of nx ny nz cells: the
//   pressure gradient telescopes through the block to -G_x (nx + 1) ny nz V on x, likewise on y,
//   and -G_z nz nx ny V on z (the ground's w rows left out), V the cell volume; with the
//   buildings' fixed rows the residual on those faces is u = 0 and the sum 0. One box, four boxes
//   (a face shared by two boxes) and any number of ranks give the same sum. R = 1 on every face gives the
//   count of blocked faces, the ground's w left out.
//   Beyond the plan: R = 1 with detJ linear in the cell indices: each face weighs the mean of its two
//   cells' J, the value of the linear J at the face (a J from one cell, or from the wrong neighbour, fails).
//   UT-6.4.3: a linear shear u = S (z - dz/2) over a solid bottom layer: the viscous flux into the
//   layer's fixed faces is the plate drag nu S Lx Ly.
#include "BuildingForcing.H"
#include "BuildingLoads.H"
#include "LaminarTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#ifdef _WIN32
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <type_traits>

namespace fs = std::filesystem;

namespace {

using gust::FaceType;
using gust::ptest::PCase;

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-12 : 1e-5;

// A map in ERF's format on 4 x 4 points 20 m apart from 0 (periodic: the tiles wrap, tile p covers
// 20 p - 10 .. 20 p + 10), height(ix, iy) in m.
std::string write_map (const std::string& name, const std::function<double(int, int)>& height)
{
    static const long id = [] () {
        long pid = static_cast<long>(::getpid());
        amrex::ParallelDescriptor::Bcast(&pid, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        return pid;
    }();
    const fs::path dir = fs::temp_directory_path() / ("gust_building_loads_" + std::to_string(id));
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(dir / name, std::ios::binary);
        os << "4\n4\n";
        for (int r = 0; r < 2; ++r) { for (int p = 0; p < 4; ++p) { os << 20 * p << "\n"; } }
        for (int i = 0; i < 4; ++i) { for (int j = 0; j < 4; ++j) { os << height(i, j) << "\n"; } }
    }
    amrex::ParallelDescriptor::Barrier();
    return (dir / name).string();
}

gust::BuildingsOptions drag_options (const std::string& file)
{
    gust::BuildingsOptions o;
    o.type = gust::BuildingsType::immersed_forcing;
    o.type_name = "immersed_forcing";
    o.small_volfrac = 0;
    o.wall_form_name = "drag";
    o.wall_form = gust::BuildingWallForm::drag;
    o.wall_speed = gust::default_wall_speed(o.wall_form);
    o.wall_tke = false;
    o.file = file;
    return o;
}

void fill_by_index (amrex::MultiFab& mf, const std::function<double(int, int, int)>& f)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { a(i, j, k) = static_cast<amrex::Real>(f(i, j, k)); });
    }
}

// 8 x 8 x 6 cells of 10 m, periodic in x and y, slip ground and top, the buildings of the map on it.
struct LoadCase
{
    PCase c;
    gust::lamtest::Laminar lam;
    gust::Buildings b;
    gust::BuildingForcing forcing;
    gust::State q;
    LoadCase (const std::string& map, int mgs, double nu)
        : c(gust::lamtest::flat_case(8, 8, 6, 80, 80, 60, true, true, mgs))
    {
        lam.make(c, nu, FaceType::slip_wall, FaceType::slip_wall);
        const gust::BuildingsOptions o = drag_options(map);
        b.build(o, c.geom, c.z, c.ba, c.dm);
        forcing.build(b, o, c.geom, c.mesh.detJ);
        c.mesh.set_open_faces(b.beta_face(0), b.beta_face(1), b.beta_face(2), b.beta(), c.geom);
        lam.in.buildings = &forcing;
        gust::lamtest::define_state(q, c);
        for (amrex::MultiFab* mf : {&q.u, &q.v, &q.w, &q.p, &q.tke}) { mf->setVal(0); }
        q.theta.setVal(300);
    }
    std::array<amrex::Real, 3> force (bool rows)
    {
        gust::ResidualInputs in = lam.in;
        in.building_rows = rows;
        gust::ResidualFields R;
        gust::residual(q, in, c.mesh, c.geom, R);
        return gust::building_force(b, c.mesh, c.geom, R);
    }
};

// UT-6.4.2
TEST (BuildingLoads, PressureGradientTelescopesThroughABlock)
{
    // Tile (2, 2): cells 3..4 in x and y; 30 m high: cells 0..2.
    const std::string map = write_map("block.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 30.0 : 0.0; });
    const double Gx = 0.7, Gy = -0.3, Gz = 0.2, V = 1000;
    for (int mgs : {8, 4}) {
        LoadCase lc(map, mgs, 0.5);
        // p at the cell centres, ghost cells included (no fill: the field is not periodic, and only the
        // block's faces are read).
        fill_by_index(lc.q.p, [&] (int i, int j, int k) { return Gx * 10 * (i + 0.5) + Gy * 10 * (j + 0.5) + Gz * 10 * (k + 0.5); });
        const auto F = lc.force(false);
        EXPECT_NEAR(F[0], -Gx * 3 * 2 * 3 * V, rtol * 1e4) << "max_grid_size " << mgs;
        EXPECT_NEAR(F[1], -Gy * 3 * 2 * 3 * V, rtol * 1e4) << "max_grid_size " << mgs;
        EXPECT_NEAR(F[2], -Gz * 3 * 2 * 2 * V, rtol * 1e4) << "max_grid_size " << mgs << ": the roof and the two faces "
                                                              "inside, not the ground's w";
        const auto F0 = lc.force(true);
        for (int n = 0; n < 3; ++n) { EXPECT_EQ(F0[n], 0) << "the fixed rows: R = u = 0 on the blocked faces"; }
        // R = 1 on every face: the count of blocked faces times V, 3 x 2 x 3 on x and y, 2 x 2 x 3 on z (the
        // ground's w, a boundary row, left out: the residual's own row there is w = 0, so only this pins it).
        gust::ResidualFields R;
        R.define(lc.c.ba, lc.c.dm);
        for (amrex::MultiFab* mf : {&R.u, &R.v, &R.w}) { mf->setVal(1); }
        const auto F1 = gust::building_force(lc.b, lc.c.mesh, lc.c.geom, R);
        EXPECT_NEAR(F1[0], -18 * V, rtol * 1e4);
        EXPECT_NEAR(F1[1], -18 * V, rtol * 1e4);
        EXPECT_NEAR(F1[2], -12 * V, rtol * 1e4);
    }
}

// Beyond the plan: the face volume J_f dx dy dzeta, J_f the mean of the face's two cells.
TEST (BuildingLoads, FaceVolumeIsTheMeanOfItsTwoCells)
{
    const std::string map = write_map("block_j.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 30.0 : 0.0; });
    const double V = 1000;
    // J at a (possibly half-integer) index point: linear, so a face's two-cell mean is J at the face.
    auto Jat = [] (double i, double j, double k) { return 1.0 + 0.1 * i + 0.03 * j + 0.05 * k; };
    for (int mgs : {8, 4}) {
        LoadCase lc(map, mgs, 0.5);
        fill_by_index(lc.c.mesh.detJ, [&] (int i, int j, int k) { return Jat(i, j, k); });
        gust::ResidualFields R;
        R.define(lc.c.ba, lc.c.dm);
        for (amrex::MultiFab* mf : {&R.u, &R.v, &R.w}) { mf->setVal(1); }
        const auto F = gust::building_force(lc.b, lc.c.mesh, lc.c.geom, R);
        // The block: cells 3..4 in x and y, 0..2 in z. Its x faces i = 3..5, y faces j = 3..5, z faces k = 1..3
        // (the ground's w left out).
        double ex = 0, ey = 0, ez = 0;
        for (int a = 3; a <= 5; ++a) {
            for (int b = 3; b <= 4; ++b) {
                for (int k = 0; k <= 2; ++k) { ex += Jat(a - 0.5, b, k); ey += Jat(b, a - 0.5, k); }
            }
        }
        for (int i = 3; i <= 4; ++i) {
            for (int j = 3; j <= 4; ++j) {
                for (int k = 1; k <= 3; ++k) { ez += Jat(i, j, k - 0.5); }
            }
        }
        EXPECT_NEAR(F[0], -ex * V, rtol * 1e4) << "max_grid_size " << mgs;
        EXPECT_NEAR(F[1], -ey * V, rtol * 1e4) << "max_grid_size " << mgs;
        EXPECT_NEAR(F[2], -ez * V, rtol * 1e4) << "max_grid_size " << mgs;
    }
}

// UT-6.4.3
TEST (BuildingLoads, ShearOverASolidLayerGivesThePlateDrag)
{
    // Every tile 10 m high: the bottom layer of cells is solid; its u faces are fixed rows at z = 5 m.
    const std::string map = write_map("layer.txt", [] (int, int) { return 10.0; });
    const double nu = 0.5, S = 0.03, Lx = 80, Ly = 80;
    for (int mgs : {8, 4}) {
        LoadCase lc(map, mgs, nu);
        // u = S (z - 5) above the layer, 0 in it and in the ground's ghost (the slip ground's mirror of 0).
        fill_by_index(lc.q.u, [&] (int, int, int k) { return k >= 1 ? S * (10 * (k + 0.5) - 5) : 0.0; });
        const auto F = lc.force(false);
        EXPECT_NEAR(F[0], nu * S * Lx * Ly, rtol * 100) << "max_grid_size " << mgs;
        EXPECT_EQ(F[1], 0);
        EXPECT_EQ(F[2], 0);
    }
}

} // namespace
