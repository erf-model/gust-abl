// The time-marching column (init.column_time > 0, ColumnGuess.H): backward Euler in physical time with a fixed step
// to a fixed end time, the form of ERF's near-1-D column. The force-driven column of 64 x 1 m layers over flat ground
// (z0 = 0.1 m), buoyancy on (g = 9.81, theta0 = 300 K), a neutral sounding of 300 K, a slip top with zero theta
// gradient and a Rayleigh layer on u, v and k (not theta) above 48 m: the column's only heat source is the surface.

#include "BoundaryConditions.H"
#include "Closure_K.H"
#include "ColumnGuess.H"
#include "ERF_MOSTUtils.H"
#include "Initializer.H"
#include "InflowProfile.H"
#include "SurfaceLayer.H"
#include "TestHelpers.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

using gust::BoundaryInputs;
using gust::FaceType;
using gust::InitOptions;
using CS = gust::InitOptions::ColumnStability;
using gust::trtest::tol;

struct Column
{
    gust::Inputs in;
    InitOptions io;
    gust::SteadyOptions so;
    BoundaryInputs bin;
    gust::SurfaceLayerOptions surf;
    gust::ClosureOptions copt;
    gust::ForcingOptions fopt;
    Column ()
    {
        in.mesh.n_cell = amrex::IntVect(4, 4, 64);
        in.mesh.prob_lo = {0, 0, 0};
        in.mesh.prob_hi = {256, 256, 64};
        in.mesh.is_periodic = {1, 1, 0};
        in.mesh.max_grid_size = {4, 4};
        in.terrain.type = "flat";
        io.first_guess = InitOptions::FirstGuess::column;
        so.stage = gust::StageKind::simplec;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer, FaceType::periodic,
                    FaceType::periodic, FaceType::slip_wall};
        bin.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
        bin.inflow.ustar = amrex::Real(0.4);
        bin.inflow.z0 = amrex::Real(0.1);
        bin.theta_ref.z = {0};
        bin.theta_ref.theta = {300};
        bin.rayleigh_zdamp = 48;
        bin.rayleigh_dampcoef = amrex::Real(0.05);
        bin.rayleigh_fields = {true, true, true, false, true};
        surf.z0 = amrex::Real(0.1);
        fopt.pressure_gradient = {amrex::Real(2.5e-3), 0, 0};
        fopt.buoyancy.gravity = amrex::Real(9.81);
        fopt.buoyancy.theta0 = 300;
    }
    void heat_flux (double q, double t, double dt)
    {
        io.column_stability = CS::heat_flux;
        surf.stability = gust::SurfaceLayerOptions::Stability::heat_flux;
        surf.heat_flux = static_cast<amrex::Real>(q);
        io.column_time = static_cast<amrex::Real>(t);
        io.column_dt = static_cast<amrex::Real>(dt);
    }
    gust::ColumnGuessResult run () const { return gust::run_column_guess(in, io, so, bin, surf, true, copt, fopt); }
};

bool mentions (const std::vector<std::string>& err, const std::string& what)
{
    for (const auto& e : err) { if (e.find(what) != std::string::npos) { return true; } }
    return false;
}

constexpr double kappa = 0.41, g = 9.81, theta0 = 300;

// The column's heat content above the sounding, sum((theta - 300) dz) over the 1 m layers (the table's centres).
double heat_content (const gust::ColumnGuessResult& r)
{
    double h = 0;
    for (std::size_t i = 1; i < r.theta.size(); ++i) { h += static_cast<double>(r.theta[i]) - 300.0; }
    return h;
}

} // namespace

// The surface-flux budget: with the fixed flux q at the ground, a slip top with zero theta gradient and no Rayleigh
// layer on theta, d/dt sum(theta dz) = q, so the heat content at t is q t to round-off (the implicit diffusion
// conserves sum(theta dz) in exact arithmetic). t = 3605 s with 10 s steps: 361 steps, the last 5 s. Heated (q = 0.02 K m/s, the
// first cell warmer than the second) and cooled (q = -0.002 K m/s, the first cell colder).
TEST(ColumnMarching, SurfaceFluxBudget)
{
    for (const double q : {0.02, -0.002}) {
        Column c;
        c.heat_flux(q, 3605, 10);
        const auto r = c.run();
        EXPECT_FALSE(r.steady);
        EXPECT_EQ(r.time, 3605);
        EXPECT_EQ(r.iterations, 361);
        EXPECT_EQ(r.columns, 1);
        EXPECT_EQ(r.flux_iterations, 0);
        EXPECT_NEAR(static_cast<double>(r.heat_flux), q, tol(1e-15, 1e-9));
        // round-off: the implicit solves of theta near 300 K at a diffusion number near 100 (about 2e-10 of q t
        // here), and the table's Real (64 sums of theta near 300 K in single precision)
        EXPECT_NEAR(heat_content(r), q * 3605, tol(1e-8 * std::abs(q * 3605), 2e-3)) << q;
        if (q > 0) {
            EXPECT_GT(r.theta[1], r.theta[2]);
            EXPECT_LT(r.inv_L, 0);
        } else {
            EXPECT_LT(r.theta[1], r.theta[2]);
            EXPECT_GT(r.inv_L, 0);
        }
    }
}

// surface_temperature: the surface 2 K below the sounding cools the column; the reported surface is the input's, the
// flux downward and from the final state: theta_1 - theta_s = theta* (ln - psi_h) / kappa with theta* = -q/u* at the
// table's first cell (MOST at the end, not one step old).
TEST(ColumnMarching, SurfaceTemperature)
{
    Column c;
    c.io.column_stability = CS::surface_temperature;
    c.io.column_surface_theta = 298;
    c.io.column_time = 1800;
    c.io.column_dt = 5;
    const auto r = c.run();
    EXPECT_FALSE(r.steady);
    EXPECT_EQ(r.iterations, 360);
    EXPECT_EQ(r.theta_surface, amrex::Real(298));
    EXPECT_EQ(r.flux_iterations, 0);
    EXPECT_LT(r.heat_flux, 0);
    EXPECT_LT(r.theta[1], r.theta[2]);
    EXPECT_LT(r.theta[2], amrex::Real(300));
    const double z1 = 0.5, z0 = 0.1, lz = std::log((z1 + z0) / z0), il = static_cast<double>(r.inv_L);
    const double us = static_cast<double>(r.ustar);
    const double ph = static_cast<double>(similarity_funs{}.calc_psi_h(static_cast<amrex::Real>(z1 * il)));
    const double pm = static_cast<double>(similarity_funs{}.calc_psi_m(static_cast<amrex::Real>(z1 * il)));
    EXPECT_NEAR(static_cast<double>(r.theta[1]) - 298, -static_cast<double>(r.heat_flux) / us / kappa * (lz - ph),
                tol(1e-6, 1e-3));
    EXPECT_NEAR(us, kappa * std::hypot(static_cast<double>(r.u[1]), static_cast<double>(r.v[1])) / (lz - pm),
                tol(1e-6, 1e-3) * us);
}

// The step count: T / dt rounded in the inputs' precision (9 s in 0.9 s steps is 10 steps, not 11 with a sliver;
// 0.1 s in 0.01 s steps is 10), and a time below one step is one step of that time (the budget q t holds).
TEST(ColumnMarching, StepCount)
{
    const std::tuple<double, double, int> cases[] = {{9, 0.9, 10}, {0.1, 0.01, 10}, {1e-6, 1000, 1}};
    for (const auto& [t, dt, n] : cases) {
        Column c;
        c.heat_flux(0.02, t, dt);
        const auto r = c.run();
        EXPECT_EQ(r.iterations, n) << t << " " << dt;
        EXPECT_NEAR(heat_content(r), 0.02 * static_cast<double>(static_cast<amrex::Real>(t)), tol(1e-8 * 0.02 * t + 1e-10, 2e-3))
            << t << " " << dt;
    }
}

// Marched long enough, the neutral column (theta held) reaches the steady column's state: the fixed point of the
// backward Euler step does not depend on the step (the steady column relaxes nu_t between steps, which leaves the
// fixed point alone).
TEST(ColumnMarching, NeutralReachesTheSteadyColumn)
{
    Column s;
    s.io.column_tol = amrex::Real(1e-10);
    const auto rs = s.run();
    ASSERT_TRUE(rs.steady);
    Column m;
    m.io.column_time = amrex::Real(2.0e6);
    m.io.column_dt = amrex::Real(200);
    const auto rm = m.run();
    EXPECT_FALSE(rm.steady);
    ASSERT_EQ(rm.u.size(), rs.u.size());
    for (std::size_t i = 0; i < rs.u.size(); ++i) {
        EXPECT_NEAR(static_cast<double>(rm.u[i]), static_cast<double>(rs.u[i]), tol(1e-6, 1e-3)) << i;
        EXPECT_NEAR(static_cast<double>(rm.tke[i]), static_cast<double>(rs.tke[i]), tol(1e-7, 1e-4)) << i;
        EXPECT_EQ(rm.theta[i], amrex::Real(300));
    }
    EXPECT_NEAR(static_cast<double>(rm.ustar), static_cast<double>(rs.ustar), tol(1e-7, 1e-4));
}

// init.column_wall_k = solved: the first cell's k is solved with no k flux through the ground. The neutral column
// marched to its steady state then satisfies the first cell's discrete k balance with no ground term,
//   gamma_f (k_2 - k_1) / (dz dz_f) + P_1 - eps_1 = 0,  gamma_f = (gamma_k,1 + gamma_k,2) / 2,
// rebuilt here from the table with the closure's point kernels (S^2 from the wall law's gradient u*/(kappa (z_1 + z0))
// at the ground face and the layers' differences above, as the column's). The fixed wall value u*^2/Cmu0^2 is not
// a solution of it (k_1 is 0.76 of it here). With a heat flux the solved k_1 is below the wall value
// (u*^3 + kappa B z_1)^(2/3)/Cmu0^2, and the heat budget is q t as with the fixed one.
TEST(ColumnMarching, SolvedWallK)
{
    Column m;
    m.io.column_time = amrex::Real(2.0e6);
    m.io.column_dt = amrex::Real(200);
    m.io.column_wall_k = InitOptions::ColumnWallK::solved;
    const auto r = m.run();
    const gust::ClosureConstants cc = gust::closure_constants(m.copt);
    const double z0 = 0.1, us = static_cast<double>(r.ustar);
    auto at = [&] (std::size_t i) { return std::complex<double>(static_cast<double>(r.u[i]), static_cast<double>(r.v[i])); };
    // the closure at cells 1 and 2 of the table (centres 0.5 and 1.5 m): S^2 the mean of the cell's two faces
    const double s0 = std::pow(us / (kappa * (0.5 + z0)), 2), s1 = std::norm(at(2) - at(1)), s2 = std::norm(at(3) - at(2));
    auto point = [&] (std::size_t i, double S2) {
        const auto k = r.tke[i];
        const amrex::Real lg = gust::closure::geom_length(static_cast<amrex::Real>(kappa * (static_cast<double>(r.z[i]) + z0)),
                                                          static_cast<amrex::Real>(m.copt.max_geom_lscale));
        const auto sw = gust::closure::switches_for(k, 0, lg, cc);
        return gust::closure::closure_point<amrex::Real>(k, 0, static_cast<amrex::Real>(S2), lg, cc, sw);
    };
    const auto p1 = point(1, 0.5 * (s0 + s1)), p2 = point(2, 0.5 * (s1 + s2));
    const double g1 = static_cast<double>(m.copt.nu_molecular + p1.nut / m.copt.sigma_k);
    const double g2 = static_cast<double>(m.copt.nu_molecular + p2.nut / m.copt.sigma_k);
    const double k1 = static_cast<double>(r.tke[1]), k2 = static_cast<double>(r.tke[2]);
    const double P1 = static_cast<double>(p1.nut) * 0.5 * (s0 + s1), eps1 = static_cast<double>(p1.eps);
    const double balance = 0.5 * (g1 + g2) * (k2 - k1) + P1 - eps1;
    EXPECT_NEAR(balance / P1, 0, tol(1e-6, 1e-3)) << P1 << " " << eps1;
    const double kw = us * us / static_cast<double>(cc.Cmu0 * cc.Cmu0);
    EXPECT_NEAR(k1 / kw, 0.76, 0.02);
    EXPECT_EQ(r.tke[0], r.tke[1]);   // the table's ground row repeats the first cell

    Column f;   // the default holds the wall value
    f.io.column_time = amrex::Real(2.0e6);
    f.io.column_dt = amrex::Real(200);
    const auto rf = f.run();
    const double usf = static_cast<double>(rf.ustar);
    EXPECT_NEAR(static_cast<double>(rf.tke[1]), usf * usf / static_cast<double>(cc.Cmu0 * cc.Cmu0), tol(1e-12, 1e-6));

    Column h;
    h.heat_flux(0.02, 3605, 10);
    h.io.column_wall_k = InitOptions::ColumnWallK::solved;
    const auto rh = h.run();
    const double uh = static_cast<double>(rh.ustar);
    const double kwh = std::pow(uh * uh * uh + kappa * g / theta0 * 0.02 * 0.5, 2.0 / 3.0) / static_cast<double>(cc.Cmu0 * cc.Cmu0);
    EXPECT_LT(static_cast<double>(rh.tke[1]), 0.85 * kwh);
    EXPECT_NEAR(heat_content(rh), 0.02 * 3605, tol(1e-8 * 0.02 * 3605, 2e-3));
}

// The surface for a given flux (column_surface_most_flux). With Businger-Dyer's stable psi = -5 zeta, u* = kappa U1 /
// (ln + 5 zeta) and zeta = -kappa g q z1 / (theta0 u*^3); on the unstable side the same with ERF's psi_m. The surface
// temperature that gives this flux (theta_1 - theta_s = theta* (ln - psi_h) / kappa) gives the same flux back through
// the surface-temperature solve (column_surface_most). A downward flux too large for the wind has no root.
TEST(ColumnMarching, FluxSurfaceMatchesTheClosedForm)
{
    const double z1 = 7.8, z0 = 0.1, U1 = 3, lz = std::log((z1 + z0) / z0);
    for (const double q : {-0.005, 0.1}) {
        double us = 0, ts = 0, il = 0;
        ASSERT_TRUE(gust::column_surface_most_flux(kappa, z0, z1, U1, q, g, theta0, 0, false, {}, us, ts, il)) << q;
        const double zeta = z1 * il;
        const double pm = static_cast<double>(similarity_funs{}.calc_psi_m(static_cast<amrex::Real>(zeta)));
        const double ph = static_cast<double>(similarity_funs{}.calc_psi_h(static_cast<amrex::Real>(zeta)));
        if (q < 0) { EXPECT_NEAR(pm, -5 * zeta, 1e-12); }
        EXPECT_NEAR(us, kappa * U1 / (lz - pm), tol(1e-9, 1e-5) * us) << q;
        EXPECT_NEAR(zeta, -kappa * g * q * z1 / (theta0 * us * us * us), tol(1e-8, 1e-5) * std::abs(zeta)) << q;
        EXPECT_NEAR(ts, -q / us, 1e-15);
        // back through the surface-temperature solve
        const double dth = ts * (lz - ph) / kappa;
        double us2 = 0, ts2 = 0, q2 = 0, il2 = 0, ch2 = 0;
        ASSERT_TRUE(gust::column_surface_most(kappa, z0, z1, U1, 300 + dth, 300, g, theta0, 0, false, {}, -4, 4, us2, ts2,
                                              q2, il2, ch2)) << q;
        EXPECT_NEAR(q2, q, tol(1e-7, 1e-4) * std::abs(q)) << q;
    }
    double us = 0, ts = 0, il = 0;
    // a seed above the upper stable root (zeta = 20 at 3 m/s, q = -0.004 K m/s, Pressel's case) still finds the root
    il = 20 / z1;
    EXPECT_TRUE(gust::column_surface_most_flux(kappa, z0, z1, U1, -0.004, g, theta0, 0, false, {}, us, ts, il));
    EXPECT_LT(z1 * il, 1);
    il = 0;
    EXPECT_FALSE(gust::column_surface_most_flux(kappa, z0, z1, 1.0, -0.5, g, theta0, 0, false, {}, us, ts, il));
    // in the column: the run stops at the step where the first cell's wind has fallen below the flux's root (60 s
    // here, from the log-law start), naming the time and bc.surface_heat_flux
    Column c;
    c.heat_flux(-0.5, 600, 10);
    std::string m = "no abort";
    try { (void)c.run(); } catch (const std::exception& e) { m = e.what(); }
    EXPECT_NE(m.find("at t = 60 s the surface flux"), std::string::npos) << m;
    EXPECT_NE(m.find("bc.surface_heat_flux"), std::string::npos) << m;
    // ending at 60 s, every step has a root but the final state does not: the run stops there too (it returned
    // u* = 1e-3 and the clamp's 1/L before)
    Column e;
    e.heat_flux(-0.5, 60, 10);
    m = "no abort";
    try { (void)e.run(); } catch (const std::exception& x) { m = x.what(); }
    EXPECT_NE(m.find("at t = 60 s (the end)"), std::string::npos) << m;
    e.io.column_time = 50;
    EXPECT_NO_THROW((void)e.run());
}

// The inputs: init.column_time >= 0 and init.column_dt > 0, finite, at most 1e7 steps; init.column_dt without a time,
// init.column_tol and init.column_max_iterations with one, a target wind and init.column_mol with one are refused.
// The marching column needs no Rayleigh layer on theta (the steady stratified column does).
TEST(ColumnMarching, InputChecks)
{
    InitOptions io;
    io.first_guess = InitOptions::FirstGuess::column;
    io.column_time = -1;
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_time"));
    {
        const gust::test::ScopedNoInvalidTrap no_trap;   // the checks test isnan first; a compiler may still compare
        io.column_time = std::numeric_limits<amrex::Real>::quiet_NaN();
        EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_time"));
    }
    io.column_time = 3600;
    EXPECT_TRUE(gust::check_init_options(io).empty());
    io.column_dt = 0;
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_dt"));
    io.column_dt = std::numeric_limits<amrex::Real>::infinity();
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_dt"));
    io.column_dt = amrex::Real(1e-4);
    EXPECT_TRUE(mentions(gust::check_init_options(io), "at most 1e7 steps"));
    io.column_dt = 1;
    io.given = {"column_tol", "column_max_iterations"};
    const auto e1 = gust::check_init_options(io);
    EXPECT_TRUE(mentions(e1, "init.column_tol has no effect with init.column_time"));
    EXPECT_TRUE(mentions(e1, "init.column_max_iterations has no effect with init.column_time"));
    io.given = {};
    io.column_target_u = 5;
    io.column_target_height = 10;
    EXPECT_TRUE(mentions(gust::check_init_options(io), "need the steady column"));
    io.column_target_u = 0;
    io.column_target_height = 0;
    io.column_mol = 100;
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_mol with init.column_time"));
    io.column_mol = 0;
    io.column_time = 0;
    io.given = {"column_dt"};
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_dt has no effect"));
    io.given = {"column_wall_k"};
    EXPECT_TRUE(mentions(gust::check_init_options(io), "init.column_wall_k has no effect without init.column_time > 0"));
    io.column_time = 3600;
    EXPECT_TRUE(gust::check_init_options(io).empty());
    io.column_time = 0;
    InitOptions none;
    none.given = {"column_wall_k"};
    EXPECT_TRUE(mentions(gust::check_init_options(none), "init.column_wall_k has no effect"));
    InitOptions prof;
    prof.first_guess = InitOptions::FirstGuess::profile;
    prof.given = {"column_wall_k"};
    EXPECT_TRUE(mentions(gust::check_init_options(prof), "init.column_wall_k has no effect without init.first_guess = column"));
    none.given = {"column_time"};
    EXPECT_TRUE(mentions(gust::check_init_options(none), "init.column_time has no effect"));

    Column r;   // the bulk Richardson bounds have no effect with the given flux
    r.heat_flux(0.01, 3600, 1);
    r.io.given = {"column_rib_max"};
    EXPECT_TRUE(mentions(gust::check_column_guess(r.io, r.so, r.fopt, r.bin, r.surf, true), "init.column_rib_max has no effect"));
    r.io.column_time = 0;
    EXPECT_FALSE(mentions(gust::check_column_guess(r.io, r.so, r.fopt, r.bin, r.surf, true), "init.column_rib_max"));

    Column c;
    c.heat_flux(0.01, 3600, 1);
    c.bin.rayleigh_zdamp = 0;
    EXPECT_TRUE(gust::check_column_guess(c.io, c.so, c.fopt, c.bin, c.surf, true).empty());
    c.io.column_time = 0;
    EXPECT_TRUE(mentions(gust::check_column_guess(c.io, c.so, c.fopt, c.bin, c.surf, true), "a Rayleigh layer on theta"));
}
