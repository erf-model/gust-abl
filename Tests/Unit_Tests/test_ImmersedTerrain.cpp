// The immersed ground (terrain.representation = immersed_forcing; ImmersedTerrain.H, ImmersedStrain.H,
// BuildingForcing.H), plan sub-phase PR-87:
//   the wall law at the exact distances returns u* of an exact log profile for any beta, below the
//   buildings' 10 z0 floor too, and its heat velocities; the vertical couplings of u and v across the
//   cut cell's top face span the distance between the fluid centroids, in the implicit rows (a_P =
//   sum a_nb kept) and in the viscous fluxes; the strain correction is the identity for a ground on a
//   cell face and the log-law gradient times the flat ground's error above a cut cell; the body
//   forces act on the open part of each face; the refusals of check_immersed_terrain; the restart record.
#include "BuildingForcing.H"
#include "Closure.H"
#include "ImmersedStrain.H"
#include "ImmersedTerrain.H"
#include "ImplicitOperator.H"
#include "LaminarTestHelpers.H"
#include "Residual.H"

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
#include <string>
#include <type_traits>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-12 : 2e-4;
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
    Ground (const PCase& c, double flat_height)
    {
        const gust::TerrainInputs t = terrain(flat_height);
        opt = gust::immersed_terrain_options(t, kZ0);
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

double log_law (double us, double d) { return us / kKappa * std::log((d + kZ0) / kZ0); }

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
        return value_at(S2, 1, 2, 3);
    };
    // beta = 0: d1 = dz / 2 and d2 = 3 dz / 2 are the flat ground's centres.
    auto any = [] (int k) { return 0.3 * k * k + std::sin(double(k)); };
    const double plain = strain(2 * kDz, any, false);
    ASSERT_GT(plain, 0);
    EXPECT_NEAR(strain(2 * kDz, any, true), plain, 1e3 * rtol * plain);

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
        EXPECT_NEAR(std::sqrt(strain((2 + f) * kDz, ulog, true)), dudz, 1e3 * rtol * dudz) << "f " << f;
    }
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
        {[] (auto& x) { x.z0_map = true; }, "bc.surface_z0_file"},
        {[] (auto& x) { x.buildings = true; }, "buildings.type"},
        {[] (auto& x) { x.nests = true; }, "nests"},
        {[] (auto& x) { x.newton = true; }, "steady.stage = simplec_newton"},
        {[] (auto& x) { x.zline_sweeps = false; }, "steady.sweep = zjacobi"},
        {[] (auto& x) { x.poisson_mlmg = true; }, "solver.poisson_backend = mlmg"},
        {[] (auto& x) { x.unsteady = true; }, "steady.stage = unsteady"},
        {[] (auto& x) { x.canopy = true; }, "canopy.type"},
        {[] (auto& x) { x.stratified = true; }, "bc.theta_profile"},
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
