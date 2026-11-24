// The immersed ground (terrain.representation = immersed_forcing; ImmersedTerrain.H, ImmersedStrain.H,
// BuildingForcing.H):
//   the wall law at the exact distances returns u* of an exact log profile for any beta, below the
//   buildings' 10 z0 floor too, and its heat velocities; the vertical couplings of u and v across the
//   cut cell's top face span the distance between the fluid centroids, in the implicit rows (a_P =
//   sum a_nb kept) and in the viscous fluxes; the strain correction is the identity for a ground on a
//   cell face and the log-law gradient times the flat ground's error above a cut cell; the body
//   forces act on the open part of each face; theta is not conducted through the ground; the refusals of check_immersed_terrain; the restart record;
//   with a roughness map (bc.surface_z0_file) each cut cell takes its column's z0, in the law and in the strain.
// turb.shear_stencil = log_law (LogShear.H): the strain of an exact log profile is the log-law gradient
// in every cell above the wall cell, over a body-fitted ground and above the cut cell of an immersed
// one, where the central difference is 20 % high in the second cell; the body-fitted wall cell and
// the top cell keep the closure's strain, the cut cell takes the wall law's gradient.
#include "BoundaryConditions.H"
#include "BuildingForcing.H"
#include "Closure.H"
#include "Heights_K.H"
#include "ImmersedStrain.H"
#include "ImmersedTerrain.H"
#include "ImplicitOperator.H"
#include "LogShear.H"
#include "MassFlux.H"
#include "LaminarTestHelpers.H"
#include "Residual.H"
#include "SurfaceRoughness.H"

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
#include <stdexcept>
#include <array>
#include <string>
#include <type_traits>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-12 : 2e-4;
// The log-law gradients: tight enough in single precision to tell them from the central difference (3 % at k = 3).
constexpr double ltol = std::is_same_v<amrex::Real, double> ? 1e-9 : 2e-3;
constexpr double kKappa = 0.41;
constexpr double kZ0 = 0.1;
constexpr double kDz = 10;   // the layers of ground_case

// 4 x 4 x 24 cells of 10 m, periodic in x and y, boxes of mgs cells.
PCase ground_case (int mgs)
{
    return gust::lamtest::flat_case(4, 4, 24, 40, 40, 240, true, true, mgs);
}

gust::TerrainInputs terrain (double flat_height)
{
    gust::TerrainInputs t;
    t.representation = "immersed_forcing";
    t.flat_height = static_cast<amrex::Real>(flat_height);
    return t;
}

// The immersed ground at flat_height and its wall law, as App.cpp builds them.
struct Ground
{
    gust::BuildingsOptions opt;
    gust::Buildings b;
    gust::BuildingForcing law;
    Ground (const PCase& c, double flat_height, bool z0_map = false)
    {
        const gust::TerrainInputs t = terrain(flat_height);
        opt = gust::immersed_terrain_options(t, kZ0, z0_map);
        b.build_from_map(opt, gust::immersed_terrain_map(t, c.geom), "terrain.flat_height", c.geom, c.z, c.ba, c.dm);
        law.build(b, opt, c.geom, c.mesh.detJ);
    }
};

std::string abort_message (const std::function<void()>& f)
{
    try {
        f();
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return "";
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

// f(k) at every point of mf, ghost points included.
void fill_by_level (amrex::MultiFab& mf, const std::function<double(int)>& f)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { a(i, j, k) = static_cast<amrex::Real>(f(k)); });
    }
}

// 0 below the ground (the ghost layer under the mesh).
double log_law (double us, double d) { return d > 0 ? us / kKappa * std::log((d + kZ0) / kZ0) : 0.0; }

std::string scratch_dir (const std::string& name)
{
    static const long id = [] () {
        long pid = static_cast<long>(::getpid());
        amrex::ParallelDescriptor::Bcast(&pid, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        return pid;
    }();
    const fs::path dir = fs::temp_directory_path() / ("gust_immersed_terrain_" + std::to_string(id)) / name;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::remove_all(dir);
        fs::create_directories(dir);
    }
    amrex::ParallelDescriptor::Barrier();
    return dir.string();
}

}   // namespace

// The law at the exact distances: an exact log profile gives its u*, whatever beta, also below 10 z0
// (where the buildings' law reads the speed as at 10 z0).
TEST (ImmersedTerrain, ExactLawReturnsTheLogProfilesUstar)
{
    const double us = 0.35;
    for (const double beta : {0.0, 0.5, 0.9}) {
        for (const double d1 : {0.05, 0.5, 5.0}) {
            const double d2 = d1 + 5;   // any d2 > d1
            const double got = gust::form7_ustar2_exact(amrex::Real(log_law(us, d1)), amrex::Real(log_law(us, d2)), amrex::Real(beta),
                                                        amrex::Real(d1), amrex::Real(d2), amrex::Real(kZ0), amrex::Real(kKappa));
            EXPECT_NEAR(got, us * us, 10 * rtol * us * us) << "beta " << beta << " d1 " << d1;
        }
    }
    // Below 10 z0 the floored law of the buildings differs: the floor is gone.
    const double d1 = 0.05, d2 = 5;
    const double floored = gust::form7_ustar2(amrex::Real(log_law(us, d1)), amrex::Real(log_law(us, d2)), amrex::Real(0.5), amrex::Real(d1),
                                              amrex::Real(d2), amrex::Real(kZ0), amrex::Real(kKappa));
    EXPECT_GT(std::abs(floored - us * us), 0.1 * us * us);
    // No reference cell (d2 = 0): the own value alone.
    EXPECT_NEAR(gust::form7_ustar2_exact(amrex::Real(log_law(us, d1)), 0, amrex::Real(0.5), amrex::Real(d1), 0, amrex::Real(kZ0), amrex::Real(kKappa)),
                us * us, 10 * rtol * us * us);
    // The heat velocities: (1 - beta) kappa u* / ln((d1 + z0) / z0) and beta kappa u* / ln((d2 + z0) / z0).
    amrex::Real own = 0, ref = 0;
    gust::form7_heat_velocities_exact(amrex::Real(us * us), amrex::Real(0.3), amrex::Real(d1), amrex::Real(d2), amrex::Real(kZ0),
                                      amrex::Real(kKappa), own, ref);
    EXPECT_NEAR(own, 0.7 * kKappa * us / std::log((d1 + kZ0) / kZ0), 10 * rtol);
    EXPECT_NEAR(ref, 0.3 * kKappa * us / std::log((d2 + kZ0) / kZ0), 10 * rtol);
}

// The cut cell (k = 2, beta = f) and its wall distances; the centroid spacing s = dz / (d2 - d1) in the
// couplings of u and v across its top face, in the implicit rows and in the viscous fluxes.
TEST (ImmersedTerrain, CentroidSpacingAcrossTheCutCellsTopFace)
{
    for (const int mgs : {4, 2}) {
        const PCase c = ground_case(mgs);
        const double f = 0.7;
        Ground g(c, (2 + f) * kDz);
        EXPECT_NEAR(value_at(g.b.beta(), 1, 2, 2), f, 100 * rtol);
        EXPECT_EQ(value_at(g.b.beta(), 1, 2, 1), 1.0);
        EXPECT_EQ(value_at(g.b.beta(), 1, 2, 3), 0.0);
        const double d1 = (1 - f) * kDz / 2, d2 = (1 - f) * kDz + kDz / 2;
        EXPECT_NEAR(value_at(g.b.wall(), 1, 2, 2, gust::kWallD1z), d1, 100 * rtol);
        EXPECT_NEAR(value_at(g.b.wall(), 1, 2, 2, gust::kWallD2z), d2, 100 * rtol);
        ASSERT_TRUE(g.law.immersed_ground());
        const double s = kDz / (d2 - d1);
        EXPECT_NEAR(value_at(g.law.centroid_spacing(), 1, 2, 2), s, 100 * rtol * s);
        EXPECT_EQ(value_at(g.law.centroid_spacing(), 1, 2, 3), 1.0);
        EXPECT_EQ(value_at(g.law.centroid_spacing(), 1, 2, 1), 1.0);

        amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1), w(c.mesh.az.boxArray(), c.dm, 1, 1);
        u.setVal(3);
        v.setVal(0);
        w.setVal(0);
        g.law.update(u, v, w);
        amrex::MultiFab ox, oy, oz;
        amrex::MultiFab* o[3] = {&ox, &oy, &oz};
        for (int d = 0; d < 3; ++d) {
            o[d]->define(g.b.beta_face(d).boxArray(), c.dm, 1, 1);
            o[d]->setVal(1);
            amrex::MultiFab::LinComb(*o[d], 1, *o[d], 0, -1, g.b.beta_face(d), 0, 0, 1, 0);
            o[d]->FillBoundary(c.geom.periodicity());
        }
        for (int comp = 0; comp < 2; ++comp) {
            amrex::MultiFab A(comp == 0 ? c.mesh.ax.boxArray() : c.mesh.ay.boxArray(), c.dm, gust::n_implicit_coefs, 0);
            A.setVal(1);
            g.law.add_to_implicit(comp, {&ox, &oy, &oz}, A);
            EXPECT_NEAR(value_at(A, 1, 2, 2, gust::coef_ZP), s, 100 * rtol * s) << comp;
            EXPECT_NEAR(value_at(A, 1, 2, 3, gust::coef_ZM), s, 100 * rtol * s) << comp;
            EXPECT_EQ(value_at(A, 1, 2, 3, gust::coef_ZP), 1.0) << comp;
            // a_P = sum a_nb: the row above the cut cell gains s - 1 with its lower coupling.
            EXPECT_NEAR(value_at(A, 1, 2, 3, gust::coef_P), s, 100 * rtol * s) << comp;
            // On the periodic seam (0) and, with boxes of 2 cells, on a box seam (2): a face there reads
            // the spacing's ghost layer.
            for (const int n : {0, 2}) {
                const int i = comp == 0 ? n : 1, j = comp == 0 ? 2 : n;
                EXPECT_NEAR(value_at(A, i, j, 2, gust::coef_ZP), s, 100 * rtol * s) << comp << " seam " << n;
                EXPECT_NEAR(value_at(A, i, j, 3, gust::coef_ZM), s, 100 * rtol * s) << comp << " seam " << n;
            }
        }
        amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, AMREX_SPACEDIM> G;
        for (int n = 0; n < 3; ++n) {
            gust::define_cv_fluxes(gust::velocity_loc(n), c.ba, c.dm, G[n]);
            for (auto& gg : G[n]) { gg.setVal(2); }
        }
        g.law.close_viscous_fluxes({&ox, &oy, &oz}, G);
        EXPECT_NEAR(value_at(G[0][2], 1, 2, 3), 2 * s, 100 * rtol * s) << "u across the cut cell's top face";
        EXPECT_NEAR(value_at(G[1][2], 1, 2, 3), 2 * s, 100 * rtol * s) << "v across the cut cell's top face";
        for (const int n : {0, 2}) {
            EXPECT_NEAR(value_at(G[0][2], n, 2, 3), 2 * s, 100 * rtol * s) << "u on seam " << n;
            EXPECT_NEAR(value_at(G[1][2], 1, n, 3), 2 * s, 100 * rtol * s) << "v on seam " << n;
        }
        EXPECT_EQ(value_at(G[0][2], 1, 2, 4), 2.0);
        EXPECT_EQ(value_at(G[2][2], 1, 2, 3), 2.0) << "w keeps its flux";
    }
}

// The strain in the first full cell above the cut cell: unchanged for a ground on a cell face (beta = 0),
// whatever the profile; for an exact log profile above a cut cell, R u* / (kappa (d2 + z0)), R the flat
// ground's error of the central difference at the same cell sizes.
TEST (ImmersedTerrain, StrainCorrectionFollowsTheFlatGround)
{
    const PCase c = ground_case(2);
    gust::BuoyancyOptions buoy;
    buoy.gravity = 0;
    auto strain = [&] (double flat_height, const std::function<double(int)>& uk, bool fix) {
        Ground g(c, flat_height);
        gust::State q;
        gust::lamtest::define_state(q, c);
        fill_by_level(q.u, uk);
        for (amrex::MultiFab* mf : {&q.v, &q.w, &q.p, &q.tke}) { mf->setVal(0); }
        q.theta.setVal(300);
        amrex::MultiFab N2, S2;
        gust::stratification_and_strain(q, buoy, c.geom, N2, S2);
        if (fix) { gust::ImmersedStrain(g.law, amrex::Real(kZ0)).apply(q, c.geom, S2); }
        return std::array<double, 2>{value_at(S2, 1, 2, 3), value_at(S2, 1, 2, 2)};   // above the cut cell; the cut cell
    };
    // beta = 0: d1 = dz / 2 and d2 = 3 dz / 2 are the flat ground's centres.
    auto any = [] (int k) { return 0.3 * k * k + std::sin(double(k)); };
    const double plain = strain(2 * kDz, any, false)[0];
    ASSERT_GT(plain, 0);
    EXPECT_NEAR(strain(2 * kDz, any, true)[0], plain, 1e3 * rtol * plain);

    const double us = 0.35;
    for (const double f : {0.3, 0.7, 0.9}) {
        const double d1 = (1 - f) * kDz / 2, d2 = (1 - f) * kDz + kDz / 2;
        auto ulog = [&] (int k) {
            if (k < 2) { return 0.0; }
            return log_law(us, k == 2 ? d1 : d2 + (k - 3) * kDz);
        };
        const double f2 = kDz / 2, f3 = 1.5 * kDz, f4 = 2.5 * kDz;
        const double R = (f3 + kZ0) * std::log((f4 + kZ0) / (f2 + kZ0)) / (f4 - f2);
        const double dudz = R * us / (kKappa * (d2 + kZ0));
        const auto with = strain((2 + f) * kDz, ulog, true);
        EXPECT_NEAR(std::sqrt(with[0]), dudz, 1e3 * rtol * dudz) << "f " << f;
        // The cut cell: the wall law's gradient at its fluid centroid; the closure's difference (which reads the
        // solid cell below) differs from it.
        const double wall = us / (kKappa * (d1 + kZ0));
        EXPECT_NEAR(std::sqrt(with[1]), wall, ltol * wall) << "cut cell, f " << f;
        EXPECT_GT(std::abs(std::sqrt(strain((2 + f) * kDz, ulog, false)[1]) - wall), 0.05 * wall) << "cut cell, f " << f;
    }
}

// A point centred in the layer the ground cuts stands at the centroid of the layer's open part, half the open
// height; the other points keep their height above the ground (below it: negative), and w's faces theirs.
TEST (ImmersedTerrain, HeightInTheCutLayerIsItsOpenCentroid)
{
    const PCase c = ground_case(2);
    for (const double f : {0.3, 0.9}) {
        gust::State q;
        q.define(c.ba, c.dm, true);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        gust::fill_immersed_ground(q.z_phys_nd, c.geom, terrain((2 + f) * kDz));
        const double ground = (2 + f) * kDz;
        // The largest error over all ranks: an EXPECT inside the loop would not fail on a non-root rank.
        amrex::Real err = 0;
        for (amrex::MFIter mfi(q.z_phys_nd, false); mfi.isValid(); ++mfi) {
            const amrex::Box zb = mfi.fabbox();
            const gust::NodeZ nz{q.z_phys_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
            const amrex::Box bx = amrex::enclosedCells(mfi.validbox());
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                const double centre = (k + 0.5) * kDz - ground;
                const double want = k == 2 ? 0.5 * (1 - f) * kDz : centre;
                for (const double h : {double(gust::height_above_ground<gust::stencil::CC>(nz, i, j, k, 0)),
                                       double(gust::height_above_ground<gust::stencil::XF>(nz, i, j, k, 0)),
                                       double(gust::height_above_ground<gust::stencil::YF>(nz, i, j, k, 0))}) {
                    err = std::max(err, amrex::Real(std::abs(h - want)));
                }
                const double hw = gust::height_above_ground<gust::stencil::ZF>(nz, i, j, k, 0);
                err = std::max(err, amrex::Real(std::abs(hw - (k * kDz - ground))));
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(err);
        EXPECT_LT(err, 1e3 * rtol * kDz) << "f " << f;
    }
}

// turb.shear_stencil = log_law over a body-fitted ground: for an exact log profile at 30 degrees to x,
// sqrt(S2) = u* / (kappa (d + z0)) in every cell between the wall cell and the top cell, across the box
// seams; the closure's central difference is higher by (d + z0) ln((d_{k+1} + z0)/(d_{k-1} + z0)) / (2 dz),
// 1.20 in the second cell; the wall cell and the top cell keep the closure's value.
TEST (LogShear, LogProfileGivesTheLogLawGradient)
{
    const PCase c = ground_case(2);
    gust::BuoyancyOptions buoy;
    buoy.gravity = 0;
    const double us = 0.35, ca = std::cos(M_PI / 6), sa = std::sin(M_PI / 6);
    auto dk = [] (int k) { return (k + 0.5) * kDz; };
    gust::State q;
    gust::lamtest::define_state(q, c);
    fill_by_level(q.u, [&] (int k) { return ca * log_law(us, dk(k)); });
    fill_by_level(q.v, [&] (int k) { return sa * log_law(us, dk(k)); });
    for (amrex::MultiFab* mf : {&q.w, &q.p, &q.tke}) { mf->setVal(0); }
    q.theta.setVal(300);
    amrex::MultiFab dist(c.ba, c.dm, 1, 1);
    fill_by_level(dist, dk);
    amrex::MultiFab N2, plain, S2;
    gust::stratification_and_strain(q, buoy, c.geom, N2, plain);
    gust::stratification_and_strain(q, buoy, c.geom, N2, S2);
    gust::LogShear(dist, gust::Roughness(amrex::Real(kZ0)), nullptr).apply(q, c.geom, S2);
    for (int i = 0; i < 4; ++i) {
        for (const int j : {0, 1, 2}) {
            EXPECT_EQ(value_at(S2, i, j, 0), value_at(plain, i, j, 0)) << "the wall cell";
            EXPECT_EQ(value_at(S2, i, j, 23), value_at(plain, i, j, 23)) << "the top cell";
            for (int k = 1; k < 23; ++k) {
                const double exact = us / (kKappa * (dk(k) + kZ0));
                EXPECT_NEAR(std::sqrt(value_at(S2, i, j, k)), exact, ltol * exact) << i << " " << j << " " << k;
                const double R = (dk(k) + kZ0) * std::log((dk(k + 1) + kZ0) / (dk(k - 1) + kZ0)) / (2 * kDz);
                EXPECT_NEAR(std::sqrt(value_at(plain, i, j, k)), R * exact, ltol * exact) << "central " << k;
            }
        }
    }
    const double R1 = (dk(1) + kZ0) * std::log((dk(2) + kZ0) / (dk(0) + kZ0)) / (2 * kDz);
    EXPECT_NEAR(R1, 1.2032, 1e-4) << "the central difference in the second cell (z0 = dz / 100)";
}

// turb.shear_stencil = log_law above an immersed ground with solid fraction f in layer 2: the solid
// layers keep the closure's strain; in the cut cell (f = 0: the first fluid cell) and every cell above,
// the exact log profile at the fluid centroids' distances gives u* / (kappa (d + z0)). The cut cell has the
// wall law's gradient (apply_cut_cell_strain), which differs from the closure's difference there. The solid
// layers hold a distance that grows upward (1 and 2 cm), so that only the solid fraction keeps them out.
TEST (LogShear, ImmersedGroundIsExactFromTheCutCellUp)
{
    const PCase c = ground_case(2);
    gust::BuoyancyOptions buoy;
    buoy.gravity = 0;
    const double us = 0.35;
    for (const double f : {0.0, 0.3, 0.7, 0.9}) {
        Ground g(c, (2 + f) * kDz);
        const double d1 = (1 - f) * kDz / 2, d2 = (1 - f) * kDz + kDz / 2;
        auto dk = [&] (int k) { return k < 2 ? 0.01 * (k + 1) : (k == 2 ? d1 : d2 + (k - 3) * kDz); };
        gust::State q;
        gust::lamtest::define_state(q, c);
        fill_by_level(q.u, [&] (int k) { return k < 2 ? 0.0 : log_law(us, dk(k)); });
        for (amrex::MultiFab* mf : {&q.v, &q.w, &q.p, &q.tke}) { mf->setVal(0); }
        q.theta.setVal(300);
        amrex::MultiFab dist(c.ba, c.dm, 1, 1);
        fill_by_level(dist, dk);
        amrex::MultiFab N2, plain, S2;
        gust::stratification_and_strain(q, buoy, c.geom, N2, plain);
        gust::stratification_and_strain(q, buoy, c.geom, N2, S2);
        gust::LogShear(dist, gust::Roughness(amrex::Real(kZ0)), &g.b).apply(q, c.geom, S2);
        for (int i = 0; i < 4; ++i) {
            for (int k = 0; k < 2; ++k) { EXPECT_EQ(value_at(S2, i, 2, k), value_at(plain, i, 2, k)) << "f " << f << " layer " << k; }
            // The cut cell: the wall law's gradient at its fluid centroid, not the closure's difference (which
            // reads the solid cell below).
            for (int k = 2; k < 23; ++k) {
                const double exact = us / (kKappa * (dk(k) + kZ0));
                EXPECT_NEAR(std::sqrt(value_at(S2, i, 2, k)), exact, ltol * exact) << "f " << f << " " << i << " " << k;
            }
            EXPECT_GT(std::abs(value_at(S2, i, 2, 2) - value_at(plain, i, 2, 2)), 0.05 * value_at(plain, i, 2, 2)) << "f " << f;
        }
    }
}

// turb.shear_stencil = log_law with a roughness map: each column's own z0 enters its gradient. With one
// profile u(k) in every column and z0 = 0.05 (1 + i) m, sqrt(S2) in column i is
// (u_{k+1} - u_{k-1}) / ((d_k + z0_i) ln((d_{k+1} + z0_i) / (d_{k-1} + z0_i))).
TEST (LogShear, RoughnessMapEntersByColumn)
{
    const PCase c = ground_case(2);
    gust::BuoyancyOptions buoy;
    buoy.gravity = 0;
    auto dk = [] (int k) { return (k + 0.5) * kDz; };
    auto uk = [&] (int k) { return log_law(0.35, dk(k)); };
    auto z0i = [] (int i) { return 0.05 * (1 + ((i % 4) + 4) % 4); };   // periodic images included
    gust::State q;
    gust::lamtest::define_state(q, c);
    fill_by_level(q.u, uk);
    for (amrex::MultiFab* mf : {&q.v, &q.w, &q.p, &q.tke}) { mf->setVal(0); }
    q.theta.setVal(300);
    amrex::MultiFab dist(c.ba, c.dm, 1, 1);
    fill_by_level(dist, dk);
    amrex::MultiFab map(gust::ground_slab(c.ba), c.dm, 1, amrex::IntVect(1, 1, 0));
    for (amrex::MFIter mfi(map, false); mfi.isValid(); ++mfi) {
        const auto a = map.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { a(i, j, k) = static_cast<amrex::Real>(z0i(i)); });
    }
    amrex::MultiFab N2, S2;
    gust::stratification_and_strain(q, buoy, c.geom, N2, S2);
    gust::LogShear(dist, gust::Roughness(amrex::Real(0), &map), nullptr).apply(q, c.geom, S2);
    for (int i = 0; i < 4; ++i) {
        for (const int k : {1, 2, 5}) {
            const double z0 = z0i(i);
            const double want = (uk(k + 1) - uk(k - 1)) / ((dk(k) + z0) * std::log((dk(k + 1) + z0) / (dk(k - 1) + z0)));
            EXPECT_NEAR(std::sqrt(value_at(S2, i, 1, k)), want, ltol * want) << i << " " << k;
        }
    }
    const double a = (uk(2) - uk(0)) / ((dk(1) + 0.05) * std::log((dk(2) + 0.05) / (dk(0) + 0.05)));
    const double b = (uk(2) - uk(0)) / ((dk(1) + 0.2) * std::log((dk(2) + 0.2) / (dk(0) + 0.2)));
    EXPECT_GT(std::abs(a / b - 1), 4e-3) << "the columns' gradients differ";
}

// bc.surface_z0_file at an immersed ground: each cut cell takes its own column's z0, in the wall law's u*^2, in
// the wall law's gradient of the cut cell's strain, in the strain of the cell above it and with
// turb.shear_stencil = log_law. One wind profile in every column and z0 = 0.05 (1 + i)(1 + j / 2) m: in column (i, j)
//   u*^2 = (1 - f) (kappa S_1 / ln((d1 + z0_ij) / z0_ij))^2 + f (kappa S_2 / ln((d2 + z0_ij) / z0_ij))^2,
// S the speed in the cut cell (1) and in the cell above (2). A uniform map gives the uniform z0's rates
// bit for bit. The law refuses to run without the map it was told to expect, and a map on another layout.
TEST (ImmersedTerrain, RoughnessMapEntersTheCutCellByColumn)
{
    gust::BuoyancyOptions buoy;
    buoy.gravity = 0;
    // By column in x and in y; periodic images included.
    auto z0i = [] (int i, int j) { return 0.05 * (1 + ((i % 4) + 4) % 4) * (1 + 0.5 * (((j % 4) + 4) % 4)); };
    auto uk = [] (int k) { return k < 2 ? 0.0 : (k == 2 ? 3.0 : (k == 3 ? 5.0 : (k == 4 ? 6.5 : 7.0))); };
    const double vfac = 0.4, sp = std::sqrt(1 + vfac * vfac);   // v = 0.4 u: the speed is sp u
    for (const int mgs : {2, 4}) {
        const PCase c = ground_case(mgs);
        amrex::MultiFab map(gust::ground_slab(c.ba), c.dm, 1, amrex::IntVect(1, 1, 0)), flat(gust::ground_slab(c.ba), c.dm, 1, amrex::IntVect(1, 1, 0));
        for (amrex::MFIter mfi(map, false); mfi.isValid(); ++mfi) {
            const auto a = map.array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { a(i, j, k) = static_cast<amrex::Real>(z0i(i, j)); });
        }
        flat.setVal(amrex::Real(kZ0));
        gust::State q;
        gust::lamtest::define_state(q, c);
        fill_by_level(q.u, uk);
        fill_by_level(q.v, [&] (int k) { return vfac * uk(k); });
        for (amrex::MultiFab* mf : {&q.w, &q.p, &q.tke}) { mf->setVal(0); }
        q.theta.setVal(300);
        amrex::MultiFab dist(c.ba, c.dm, 1, 1);
        for (const double f : {0.3, 0.7}) {
            const double d1 = (1 - f) * kDz / 2, d2 = (1 - f) * kDz + kDz / 2;
            Ground g(c, (2 + f) * kDz, true);
            EXPECT_NE(abort_message([&] { g.law.update(q.u, q.v, q.w); }).find("bc.surface_z0_file"), std::string::npos)
                << "the law does not run on the uniform z0 when a map is expected";
            amrex::MultiFab bare(gust::ground_slab(c.ba), c.dm, 1, 0);
            bare.setVal(amrex::Real(kZ0));
            EXPECT_NE(abort_message([&] { g.law.set_ground_roughness(bare); }).find("BuildingForcing::set_ground_roughness"), std::string::npos);
            g.law.set_ground_roughness(map);
            g.law.update(q.u, q.v, q.w);
            auto us2 = [&] (double z0, double kappa) {
                const double own = kappa * sp * uk(2) / std::log((d1 + z0) / z0), ref = kappa * sp * uk(3) / std::log((d2 + z0) / z0);
                return (1 - f) * own * own + f * ref * ref;
            };
            amrex::MultiFab N2, S2, L2;
            gust::stratification_and_strain(q, buoy, c.geom, N2, S2);
            gust::stratification_and_strain(q, buoy, c.geom, N2, L2);
            gust::ImmersedStrain(g.law, gust::Roughness(amrex::Real(0), &map)).apply(q, c.geom, S2);
            fill_by_level(dist, [&] (int k) { return k < 2 ? 0.01 * (k + 1) : (k == 2 ? d1 : d2 + (k - 3) * kDz); });
            gust::LogShear(dist, gust::Roughness(amrex::Real(0), &map), &g.b).apply(q, c.geom, L2);
            for (int i = 0; i < 4; ++i) {
                for (const int j : {0, 3}) {
                    const double z0 = z0i(i, j);
                    const double want = us2(z0, kKappa);
                    EXPECT_NEAR(value_at(g.law.rate(), i, j, 2, gust::kBuildingUstar2), want, 100 * rtol * want) << "u*^2, f " << f << " column " << i << " " << j;
                    // The cut cell's strain: the wall law's gradient u* / (kappa (d1 + z0_i)).
                    const double wall = std::sqrt(us2(z0, 1.0)) / (d1 + z0);
                    EXPECT_NEAR(std::sqrt(value_at(S2, i, j, 2)), wall, ltol * wall) << "cut cell, f " << f << " column " << i << " " << j;
                    EXPECT_NEAR(std::sqrt(value_at(L2, i, j, 2)), wall, ltol * wall) << "log_law cut cell, f " << f << " column " << i << " " << j;
                    // The cell above: the difference over (d2 + z0_i) ln((d2 + dz + z0_i) / (d1 + z0_i)) times the
                    // flat ground's ratio (ImmersedStrain.H).
                    const double f2 = kDz / 2, f3 = 1.5 * kDz, f4 = 2.5 * kDz;
                    const double R = (f3 + z0) * std::log((f4 + z0) / (f2 + z0)) / (f4 - f2);
                    const double above = R * sp * (uk(4) - uk(2)) / ((d2 + z0) * std::log((d2 + kDz + z0) / (d1 + z0)));
                    EXPECT_NEAR(std::sqrt(value_at(S2, i, j, 3)), above, ltol * above) << "above, f " << f << " column " << i << " " << j;
                }
            }
            EXPECT_GT(us2(z0i(3, 0), kKappa), 1.3 * us2(z0i(0, 0), kKappa)) << "the columns' u*^2 differ along x";
            EXPECT_GT(us2(z0i(0, 3), kKappa), 1.3 * us2(z0i(0, 0), kKappa)) << "the columns' u*^2 differ along y";
            // A uniform map: the uniform z0's rates, bit for bit.
            Ground gu(c, (2 + f) * kDz), gm(c, (2 + f) * kDz, true);
            gu.law.update(q.u, q.v, q.w);
            gm.law.set_ground_roughness(flat);
            gm.law.update(q.u, q.v, q.w);
            for (int i = 0; i < 4; ++i) {
                for (int n = 0; n <= gust::kBuildingUstar2; ++n) {
                    EXPECT_EQ(value_at(gm.law.rate(), i, 1, 2, n), value_at(gu.law.rate(), i, 1, 2, n)) << "uniform map, component " << n;
                }
            }
        }
    }
    // A body-fitted ground's buildings take no ground map.
    gust::BuildingForcing none;
    amrex::MultiFab any;
    EXPECT_NE(abort_message([&] { none.set_ground_roughness(any); }).find("needs an immersed ground"), std::string::npos);
}

// The refusals of LogShear: a roughness length that is not positive, a wall distance without ghost cells.
TEST (LogShear, RefusesWhatItCannotUse)
{
    const PCase c = ground_case(4);
    amrex::MultiFab dist(c.ba, c.dm, 1, 1), bare(c.ba, c.dm, 1, 0);
    dist.setVal(1);
    bare.setVal(1);
    EXPECT_NE(abort_message([&] { gust::LogShear(dist, gust::Roughness(amrex::Real(0)), nullptr); })
                  .find("turb.shear_stencil = log_law needs a roughness length z0 > 0"),
              std::string::npos);
    EXPECT_NE(abort_message([&] { gust::LogShear(bare, gust::Roughness(amrex::Real(kZ0)), nullptr); })
                  .find("the wall distance needs one lateral ghost cell"),
              std::string::npos);
}

// The body forces on the open part of each face: with q = 0 the residual of u is the pressure gradient
// times the x-face's open fraction, so the cut cell's face carries 1 - f of a fluid face's.
TEST (ImmersedTerrain, BodyForcesActOnTheOpenPart)
{
    for (const int mgs : {4, 2}) {
        PCase c = ground_case(mgs);
        const double f = 0.7, F = 1.225e-4;
        Ground g(c, (2 + f) * kDz);
        gust::lamtest::Laminar lam;
        lam.make(c, 1.0, gust::FaceType::slip_wall, gust::FaceType::slip_wall);
        c.mesh.set_open_faces(g.b.beta_face(0), g.b.beta_face(1), g.b.beta_face(2), g.b.beta(), c.geom);
        lam.in.buildings = &g.law;
        lam.in.pressure_gradient = {static_cast<amrex::Real>(F), 0, 0};
        gust::State q;
        gust::lamtest::define_state(q, c);
        for (amrex::MultiFab* mf : {&q.u, &q.v, &q.w, &q.p, &q.tke}) { mf->setVal(0); }
        q.theta.setVal(300);
        gust::ResidualFields R;
        gust::residual(q, lam.in, c.mesh, c.geom, R);
        const double fluid = value_at(R.u, 1, 2, 5), cut = value_at(R.u, 1, 2, 2);
        ASSERT_NEAR(std::abs(fluid), F, 100 * rtol * F);
        EXPECT_NEAR(cut, (1 - f) * fluid, 100 * rtol * F) << "max_grid_size " << mgs;
        EXPECT_EQ(value_at(R.u, 1, 2, 1), 0.0) << "nothing in the solid";
    }
}

// The immersed ground takes no heat: theta is not conducted through the solid part of a face. The cut cell's
// residual of theta and its implicit row do not read theta in the solid cell below; the first fluid cell's
// row is unchanged. Without the closure the solid cells would carry heat in and out of the air.
TEST (ImmersedTerrain, ThetaIsNotConductedThroughTheGround)
{
    for (const int mgs : {4, 2}) {
        PCase c = ground_case(mgs);
        const double f = 0.7;
        Ground g(c, (2 + f) * kDz);
        gust::lamtest::Laminar lam;
        lam.make(c, 1.0, gust::FaceType::slip_wall, gust::FaceType::slip_wall);
        c.mesh.set_open_faces(g.b.beta_face(0), g.b.beta_face(1), g.b.beta_face(2), g.b.beta(), c.geom);
        lam.in.buildings = &g.law;
        auto R_theta = [&] (double solid, int k) {
            gust::State q;
            gust::lamtest::define_state(q, c);
            for (amrex::MultiFab* mf : {&q.u, &q.v, &q.w, &q.p, &q.tke}) { mf->setVal(0); }
            fill_by_level(q.theta, [&] (int kk) { return kk <= 1 ? solid : 300 + 0.01 * (kk + 0.5) * kDz; });
            gust::ResidualFields R;
            gust::residual(q, lam.in, c.mesh, c.geom, R);
            return value_at(R.theta, 1, 2, k);
        };
        const double cut = R_theta(300, 2), hot = R_theta(400, 2);
        ASSERT_GT(std::abs(cut), 1e-6) << "no flux through the cut cell's top face: the case tests nothing otherwise";
        EXPECT_EQ(hot, cut) << "the cut cell's residual reads the solid's theta, max_grid_size " << mgs;
        EXPECT_EQ(R_theta(400, 3), R_theta(300, 3)) << "the first fluid cell, max_grid_size " << mgs;

        gust::State q;
        gust::lamtest::define_state(q, c);
        for (amrex::MultiFab* mf : {&q.u, &q.v, &q.w}) { mf->setVal(0); }
        amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;
        gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
        amrex::MultiFab A;
        gust::implicit_rows({&F[0], &F[1], &F[2]}, q.z_phys_nd, lam.in, gust::Equation::theta, c.mesh, c.geom, A);
        EXPECT_EQ(value_at(A, 1, 2, 2, gust::coef_ZM), 0.0) << "the coupling to the solid below, max_grid_size " << mgs;
        EXPECT_GT(value_at(A, 1, 2, 2, gust::coef_ZP), 0.0) << "the coupling to the air above stays";
        EXPECT_GT(value_at(A, 1, 2, 2, gust::coef_XM), 0.0) << "the open part of the side face stays";
    }
}

// An inflow face over an immersed ground: the normal velocity is 0 on the solid faces. A layer whose open part
// is below terrain.small_volfrac is solid but still has a height above the ground, so the profile alone would
// give its inflow face a wind, which the solid rows of the residual (they ask for 0) could never remove.
TEST (ImmersedTerrain, InflowIsZeroOnSolidFaces)
{
    for (const int mgs : {4, 2}) {
        PCase c = gust::lamtest::flat_case(4, 4, 24, 40, 40, 240, false, true, mgs);
        gust::TerrainInputs t = terrain(2.95 * kDz);   // layer 2 is 5 % open
        t.small_volfrac = amrex::Real(0.1);
        const gust::BuildingsOptions opt = gust::immersed_terrain_options(t, kZ0);
        gust::Buildings b;
        b.build_from_map(opt, gust::immersed_terrain_map(t, c.geom), "terrain.flat_height", c.geom, c.z, c.ba, c.dm);
        gust::State q;
        q.define(c.ba, c.dm, true);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        gust::fill_immersed_ground(q.z_phys_nd, c.geom, t);
        for (amrex::MultiFab* mf : {&q.u, &q.v, &q.w, &q.p, &q.tke}) { mf->setVal(1); }
        q.theta.setVal(300);
        gust::BoundaryInputs in;
        in.face = {gust::FaceType::inflow, gust::FaceType::periodic, gust::FaceType::slip_wall,
                   gust::FaceType::outflow, gust::FaceType::periodic, gust::FaceType::slip_wall};
        gust::BoundaryConditions bc(in, c.geom, c.mesh, gust::CoriolisOptions{});
        const gust::BuoyancyOptions buoy{};
        bc.fill_ghosts(q, nullptr, buoy);
        const double without = value_at(q.u, 0, 1, 2);
        ASSERT_GT(without, 0.1) << "the profile gives the solid layer's face a wind: the case tests nothing otherwise";
        bc.set_solid_faces(&b.beta_face(0), &b.beta_face(1));
        bc.fill_ghosts(q, nullptr, buoy);
        ASSERT_EQ(value_at(b.beta_face(0), 0, 1, 2), 1.0) << "the 5 % layer is solid";
        for (int k = 0; k < 3; ++k) { EXPECT_EQ(value_at(q.u, 0, 1, k), 0.0) << "solid inflow face, k = " << k << ", max_grid_size " << mgs; }
        EXPECT_GT(value_at(q.u, 0, 1, 3), 0.1) << "the first open layer keeps the profile";
        EXPECT_EQ(value_at(q.u, 4, 1, 2), 1.0) << "the outflow face is not touched";
    }
}

// bc.rayleigh_height = absolute over an immersed ground: the layer starts bc.rayleigh_zdamp above the
// terrain's zero, as on a mesh that follows the terrain, not above the mesh's lower face (terrain.flat_height
// below it), and its ramp ends at the top level. w alone is damped (towards 0), so R.w is the coefficient.
TEST (ImmersedTerrain, AbsoluteRayleighHeightIsFromTheTerrainZero)
{
    for (const int mgs : {4, 2}) {
        PCase c = gust::lamtest::flat_case(4, 4, 24, 40, 40, 240, false, true, mgs);
        const double hf = 2.5 * kDz;   // the terrain's zero, 25 m above the mesh's lower face
        gust::State q;
        q.define(c.ba, c.dm, true);
        amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
        gust::fill_immersed_ground(q.z_phys_nd, c.geom, terrain(hf));
        for (amrex::MultiFab* mf : {&q.u, &q.v, &q.p, &q.tke}) { mf->setVal(1); }
        q.w.setVal(1);
        q.theta.setVal(300);
        gust::BoundaryInputs in;
        in.face = {gust::FaceType::inflow, gust::FaceType::periodic, gust::FaceType::slip_wall,
                   gust::FaceType::outflow, gust::FaceType::periodic, gust::FaceType::slip_wall};
        in.rayleigh_zdamp = 150;
        in.rayleigh_dampcoef = amrex::Real(0.05);
        in.rayleigh_fields = {false, false, true, false, false};
        in.height_origin = static_cast<amrex::Real>(hf);   // as App.cpp sets it for terrain.flat_height
        gust::BoundaryConditions bc(in, c.geom, c.mesh, gust::CoriolisOptions{});
        gust::ResidualFields R;
        R.define(c.ba, c.dm);
        for (auto& f : R.fields()) { f.mf->setVal(0.0); }
        bc.add_damping(q, R);
        double err = 0;
        int damped = 0;
        for (int k = 0; k <= 24; ++k) {
            const double za = k * kDz - hf;   // the w point's height above the terrain's zero
            const double s = std::sin(0.5 * M_PI * (za - 150) / (240 - hf - 150));
            const double expect = za > 150 ? 0.05 * s * s : 0.0;
            damped += expect > 0 ? 1 : 0;
            err = std::max(err, std::abs(value_at(R.w, 1, 1, k) - expect));
        }
        ASSERT_EQ(damped, 7) << "the w points 155 to 215 m above the terrain's zero are in the layer";
        // single precision: the round-off of 0.05 times sin^2, 6e-9
        EXPECT_LT(err, (std::is_same_v<amrex::Real, double>) ? 1e-12 : 1e-7) << "max_grid_size " << mgs;
    }
}

// The refusals that need other namespaces; nothing without a solve.
TEST (ImmersedTerrain, RefusesWhatItDoesNotSupport)
{
    gust::ImmersedTerrainFacts ok;
    ok.surface_layer = true;
    EXPECT_TRUE(gust::check_immersed_terrain(ok).empty());
    struct Case { std::function<void(gust::ImmersedTerrainFacts&)> set; std::string key; };
    const std::vector<Case> cases{
        {[] (auto& x) { x.surface_layer = false; }, "bc.zlo = surface_layer"},
        {[] (auto& x) { x.neutral = false; }, "bc.surface_stability = neutral"},
        {[] (auto& x) { x.buildings = true; }, "buildings.type"},
        {[] (auto& x) { x.nests = true; }, "nests"},
        {[] (auto& x) { x.newton = true; }, "steady.stage = simplec_newton"},
        {[] (auto& x) { x.poisson_mlmg = true; }, "solver.poisson_backend = mlmg"},
        {[] (auto& x) { x.unsteady = true; }, "steady.stage = unsteady"},
        {[] (auto& x) { x.canopy_heat = true; }, "canopy.heat = fixed"},
        {[] (auto& x) { x.stratified = true; }, "bc.theta_profile"},
        {[] (auto& x) { x.outflow_from_velocity = true; }, "bc.outflow_transport"},
        {[] (auto& x) { x.column_guess = true; }, "init.first_guess = column"},
        {[] (auto& x) { x.probes = true; }, "steady.probes_file"},
        {[] (auto& x) { x.zref_given = true; }, "bc.surface_zref"},
        {[] (auto& x) { x.floor_given = true; }, "bc.surface_denominator_floor"},
        {[] (auto& x) { x.wall_k_source = true; }, "bc.surface_wall_k"},
        {[] (auto& x) { x.laminar = true; }, "turb.closure = kl"}};
    for (const auto& cs : cases) {
        gust::ImmersedTerrainFacts x = ok;
        cs.set(x);
        const auto err = gust::check_immersed_terrain(x);
        ASSERT_EQ(err.size(), 1u) << cs.key;
        EXPECT_NE(err[0].find(cs.key), std::string::npos) << err[0];
        EXPECT_NE(err[0].find("terrain.representation = immersed_forcing"), std::string::npos) << err[0];
    }
    // The ground inside the first layer or the top layer is refused; on their faces and between it is accepted.
    const double top = 23 * kDz;
    EXPECT_NE(gust::check_immersed_terrain_height(terrain(0.5 * kDz), kDz, top).find("terrain.flat_height must lie between"), std::string::npos);
    EXPECT_NE(gust::check_immersed_terrain_height(terrain(23.5 * kDz), kDz, top).find("terrain.flat_height must lie between"), std::string::npos);
    EXPECT_TRUE(gust::check_immersed_terrain_height(terrain(kDz), kDz, top).empty());
    EXPECT_TRUE(gust::check_immersed_terrain_height(terrain(top), kDz, top).empty());
    EXPECT_TRUE(gust::check_immersed_terrain_height(terrain(1.5 * kDz), kDz, top).empty());
    // A hybrid mesh: from 0 (the lower face is the ground where the envelope meets the terrain).
    auto hyb = [] (double h) { auto t = terrain(h); t.representation = "hybrid"; return t; };
    EXPECT_TRUE(gust::check_immersed_terrain_height(hyb(0), kDz, top).empty());
    EXPECT_TRUE(gust::check_immersed_terrain_height(hyb(0.5 * kDz), kDz, top).empty());
    EXPECT_NE(gust::check_immersed_terrain_height(hyb(-0.1), kDz, top).find("between 0 and"), std::string::npos);
    EXPECT_NE(gust::check_immersed_terrain_height(hyb(23.5 * kDz), kDz, top).find("between 0 and"), std::string::npos);
    // The advice: immersed_forcing over a hill or a file recommends hybrid; flat terrain, hybrid and body_fitted get none.
    auto rep = [] (const std::string& r, const std::string& type) { auto t = terrain(kDz); t.representation = r; t.type = type; return t; };
    for (const std::string type : {"agnesi", "xyz", "erf_grid"}) {
        EXPECT_NE(gust::immersed_terrain_advice(rep("immersed_forcing", type)).find("terrain.representation = hybrid is recommended"),
                  std::string::npos) << type;
        EXPECT_TRUE(gust::immersed_terrain_advice(rep("hybrid", type)).empty()) << type;
        EXPECT_TRUE(gust::immersed_terrain_advice(rep("body_fitted", type)).empty()) << type;
    }
    EXPECT_TRUE(gust::immersed_terrain_advice(rep("immersed_forcing", "flat")).empty());
}

// The restart record: the hash follows flat_height, small_volfrac and z0; a checkpoint is refused for a
// changed ground and across representations; the ground writes no Buildings record.
TEST (ImmersedTerrain, RestartRecordRefusesAChangedGround)
{
    const gust::TerrainInputs t = terrain(27);
    const std::string h = gust::immersed_terrain_hash(t, kZ0);
    EXPECT_EQ(h.size(), 64u);
    EXPECT_EQ(h, gust::immersed_terrain_hash(t, kZ0));
    gust::TerrainInputs t2 = t;
    t2.flat_height = 27.5;
    EXPECT_NE(h, gust::immersed_terrain_hash(t2, kZ0));
    t2 = t;
    t2.small_volfrac = amrex::Real(0.01);
    EXPECT_NE(h, gust::immersed_terrain_hash(t2, kZ0));
    EXPECT_NE(h, gust::immersed_terrain_hash(t, 0.2));

    // bc.surface_z0_file: the file's contents and the IDW choices stand for z0, which is not read.
    {
        const std::string dir = scratch_dir("z0");
        const std::string a = dir + "/a.xyz", b = dir + "/b.xyz", c = dir + "/c.xyz";
        if (amrex::ParallelDescriptor::IOProcessor()) {
            std::ofstream(a) << "0 0 0.03\n40 40 0.3\n";
            std::ofstream(b) << "0 0 0.03\n40 40 0.31\n";
            std::ofstream(c) << "0 0 0.03\n40 40 0.3\n";
        }
        amrex::ParallelDescriptor::Barrier();
        const std::string ha = gust::immersed_terrain_hash(t, kZ0, a);
        EXPECT_EQ(ha.size(), 64u);
        EXPECT_NE(ha, h);
        EXPECT_EQ(ha, gust::immersed_terrain_hash(t, 0.2, a)) << "z0 is not read with a map";
        EXPECT_NE(ha, gust::immersed_terrain_hash(t, kZ0, b)) << "the file's contents";
        EXPECT_EQ(ha, gust::immersed_terrain_hash(t, kZ0, c)) << "the file's path is not in the hash";
        EXPECT_NE(abort_message([&] { gust::immersed_terrain_hash(t, kZ0, dir + "/none.xyz"); }).find("bc.surface_z0_file"), std::string::npos)
            << "a map that cannot be read";
        t2 = t;
        t2.idw_neighbors = t.idw_neighbors + 1;
        EXPECT_NE(ha, gust::immersed_terrain_hash(t2, kZ0, a)) << "terrain.idw_neighbors";
        t2 = t;
        t2.idw_power = t.idw_power + amrex::Real(1);
        EXPECT_NE(ha, gust::immersed_terrain_hash(t2, kZ0, a)) << "terrain.idw_power";
    }
    t2 = t;
    t2.small_volfrac = amrex::Real(0.01);

    const std::string imm = scratch_dir("imm"), bf = scratch_dir("bf");
    gust::write_immersed_terrain_checkpoint(imm, h);
    gust::write_immersed_terrain_checkpoint(bf, "");
    amrex::ParallelDescriptor::Barrier();
    EXPECT_TRUE(fs::exists(gust::immersed_terrain_checkpoint_file(imm)));
    EXPECT_FALSE(fs::exists(gust::immersed_terrain_checkpoint_file(bf))) << "a body-fitted run writes none";
    EXPECT_EQ(abort_message([&] () { gust::check_immersed_terrain_restart(imm, h); }), "");
    EXPECT_EQ(abort_message([&] () { gust::check_immersed_terrain_restart(bf, ""); }), "");
    EXPECT_NE(abort_message([&] () { gust::check_immersed_terrain_restart(imm, gust::immersed_terrain_hash(t2, kZ0)); }).find("differ"),
              std::string::npos);
    EXPECT_NE(abort_message([&] () { gust::check_immersed_terrain_restart(imm, ""); }).find("immersed"), std::string::npos);
    EXPECT_NE(abort_message([&] () { gust::check_immersed_terrain_restart(bf, h); }).find("body-fitted"), std::string::npos);

    // A record of another version (the flat ground's of the first release) is refused by its version.
    const std::string v1 = scratch_dir("v1");
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(gust::immersed_terrain_checkpoint_file(v1));
        os << "gust_immersed_terrain 1\nsha256 " << h << "\n";
    }
    amrex::ParallelDescriptor::Barrier();
    EXPECT_NE(abort_message([&] () { gust::check_immersed_terrain_restart(v1, h); }).find("record is version 1"), std::string::npos);

    // The Buildings record stays empty: the immersed ground owns its own.
    const PCase c = ground_case(4);
    Ground g(c, 27);
    EXPECT_TRUE(g.b.input_hash().empty());
    const std::string dir = scratch_dir("buildings");
    g.b.write_checkpoint(dir);
    amrex::ParallelDescriptor::Barrier();
    EXPECT_FALSE(fs::exists(gust::buildings_checkpoint_file(dir)));
}
