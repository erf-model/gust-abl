// bc.surface_wstar = beljaars (SurfaceWstar_K.H): w* and U_eff match their closed forms; with a ground temperature
// the column's fixed point U_eff = sqrt(U^2 + (1.2 w*(q(U_eff)))^2) matches an independent bisection on U_eff; a
// stable or neutral column keeps U_eff = U and the solve without w*, bit for bit. In the wall model
// (compute_surface_state) a column with w* has the u*, q and drag of the same column without w* at the wind U_eff,
// for a given flux and a given ground temperature, Businger-Dyer and Jimenez; the stress keeps the wind's own
// direction and size, tau = drag U_t; z_i comes from the value or from the map; the inputs are checked.
#include "SurfaceLayer.H"
#include "SurfaceTemperature.H"
#include "SurfaceWstar_K.H"
#include "TestHelpers.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>

namespace {

using amrex::Real;

constexpr double kKappa = 0.41;
constexpr double kGoth = 9.81 / 300.0;
constexpr double kD = 5, kZ0 = 0.1;   // the 5 m first centre over z0 = 0.1 m (test_SurfaceTemperature.cpp)

double rtol () { return 1000 * static_cast<double>(std::numeric_limits<Real>::epsilon()); }

gust::SurfaceTemperatureMost solve (double U, double dtheta, double kb = 0)
{
    const double lz = std::log((kD + kZ0) / kZ0);
    return gust::surface_temperature_most(Real(kKappa), Real(U), Real(kD), Real(lz), Real(1), Real(kGoth),
                                          Real(dtheta), gust::kSurfaceTemperatureRibMin,
                                          gust::kSurfaceTemperatureRibMax, Real(kb));
}

gust::SurfaceTemperatureMost solve_wstar (double U, double dtheta, double zi, double& ueff, double kb = 0)
{
    const double lz = std::log((kD + kZ0) / kZ0);
    Real u_eff = 0;
    const auto r = gust::surface_temperature_most_wstar(Real(kKappa), Real(U), Real(kD), Real(lz), Real(1),
                                                        Real(kGoth), Real(dtheta), gust::kSurfaceTemperatureRibMin,
                                                        gust::kSurfaceTemperatureRibMax, Real(kb), Real(zi), u_eff);
    ueff = static_cast<double>(u_eff);
    return r;
}

// The independent root: bisection on h(V) = V - sqrt(U^2 + (1.2 (g/theta_0 z_i q(V))^(1/3))^2), q(V) the
// ground-temperature solve at V (tested against its own bisection in test_SurfaceTemperature.cpp). h(U) <= 0
// (the square root is at least U) and h(V) > 0 for V large (q grows at most linearly, so its cube root falls
// behind V).
double bisect_ueff (double U, double dtheta, double zi, double kb = 0)
{
    auto h = [&] (double V) {
        const double q = static_cast<double>(solve(V, dtheta, kb).q);
        const double w = (q > 0) ? std::cbrt(kGoth * zi * q) : 0.0;
        return V - std::sqrt(U * U + 1.44 * w * w);
    };
    double lo = U, hi = U + 1;
    while (h(hi) <= 0) { hi *= 2; }
    for (int it = 0; it < 200 && hi - lo > 4 * std::numeric_limits<double>::epsilon() * hi; ++it) {
        const double mid = 0.5 * (lo + hi);
        if (h(mid) <= 0) { lo = mid; } else { hi = mid; }
    }
    return 0.5 * (lo + hi);
}

TEST(SurfaceWstar, ClosedForms)
{
    // w* = (9.81 / 300 * 1000 * 0.2)^(1/3) = 6.54^(1/3) = 1.87008...
    const double w = static_cast<double>(gust::convective_velocity(Real(kGoth), Real(1000), Real(0.2)));
    EXPECT_NEAR(w, 1.8700760, 1e-6);
    EXPECT_NEAR(w * w * w, 6.54, 6.54 * rtol());
    EXPECT_EQ(gust::convective_velocity(Real(kGoth), Real(1000), Real(0)), Real(0));
    EXPECT_EQ(gust::convective_velocity(Real(kGoth), Real(1000), Real(-0.05)), Real(0));
    // U_eff = sqrt(3^2 + (1.2 w*)^2) = sqrt(9 + 1.44 * 6.54^(2/3)) = 3.74646...
    const double u_eff = static_cast<double>(gust::beljaars_speed(Real(3), Real(kGoth), Real(1000), Real(0.2)));
    EXPECT_NEAR(u_eff, std::sqrt(9 + 1.44 * w * w), 4 * rtol());
    EXPECT_NEAR(u_eff, 3.7464577, 1e-6);
    EXPECT_EQ(gust::beljaars_speed(Real(3), Real(kGoth), Real(1000), Real(-0.1)), Real(3));
}

TEST(SurfaceWstar, FixedPointMatchesTheBisection)
{
    struct Case { double U, dtheta, zi, kb; };
    // Unstable columns (theta_r < theta_s): calm to windy, shallow to deep boundary layers, with and without kB.
    for (const Case c : {Case{0.3, -4, 1500, 0}, Case{1, -3, 800, 0}, Case{3, -5, 1900, 0}, Case{8, -2, 1200, 0},
                         Case{2, -4, 1000, 2}, Case{0.05, -1, 500, 0}}) {
        SCOPED_TRACE("U " + std::to_string(c.U) + ", dtheta " + std::to_string(c.dtheta) + ", z_i "
                     + std::to_string(c.zi) + ", kB " + std::to_string(c.kb));
        double u_eff = 0;
        const auto r = solve_wstar(c.U, c.dtheta, c.zi, u_eff, c.kb);
        EXPECT_TRUE(r.converged);
        const double ref = bisect_ueff(c.U, c.dtheta, c.zi, c.kb);
        EXPECT_NEAR(u_eff, ref, 10 * rtol() * ref);
        EXPECT_GT(u_eff, c.U);
        // The flux and u* are the ground-temperature solve's at U_eff.
        const auto at = solve(u_eff, c.dtheta, c.kb);
        EXPECT_NEAR(static_cast<double>(r.q), static_cast<double>(at.q), rtol() * std::abs(at.q));
        EXPECT_NEAR(static_cast<double>(r.ustar), static_cast<double>(at.ustar), rtol() * at.ustar);
        // w* raises u*. (Not always q: with Businger-Dyer's psi_h and its floor ch >= 1, q(U) falls between
        // about 0.7 and 1.6 m/s under a 3 K excess, so U = 1 m/s ends with less flux than without w*.)
        EXPECT_GT(r.ustar, solve(c.U, c.dtheta, c.kb).ustar);
    }
    // Pinned: U = 1 m/s, dtheta = -3 K, z_i = 800 m: U_eff = 2.0396978 m/s and q = 0.12429 K m/s from an
    // independent Python fixed point (Businger-Dyer psi written out, bisection on zeta).
    double u_eff = 0;
    const auto r = solve_wstar(1, -3, 800, u_eff);
    EXPECT_NEAR(u_eff, 2.0396978, 1e-6);
    EXPECT_NEAR(static_cast<double>(r.q), 0.12429, 1e-5);
}

TEST(SurfaceWstar, StableAndNeutralColumnsKeepTheWind)
{
    for (const double dtheta : {0.0, 0.5, 3.0}) {
        double u_eff = -1;
        const auto r = solve_wstar(2, dtheta, 1000, u_eff);
        const auto s = solve(2, dtheta);
        EXPECT_EQ(u_eff, 2.0);
        EXPECT_EQ(r.q, s.q);
        EXPECT_EQ(r.ustar, s.ustar);
        EXPECT_EQ(r.zeta, s.zeta);
    }
}

// ---------------------------------------------------------------------------------------------- the wall model

using gust::SurfaceLayerOptions;
using gust::ptest::PCase;
using gust::trtest::tol;
namespace fs = std::filesystem;

const gust::BuoyancyOptions kBuoyancy{};   // g = 9.81 m/s^2, theta_0 = 300 K

// A flat mesh, nx x nx x nz cells over 400 m x 400 m, geometric layers (dz0 = 10 m: the first centre at 5 m), boxes
// of mgs cells (test_SurfaceTemperature.cpp's).
PCase flat_column (int nz, int nx = 4, int mgs = 4, bool periodic = true)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, nx, nz);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {periodic ? 1 : 0, periodic ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.mesh.dz0 = amrex::Real(10);
    c.in.mesh.stretch_ratio = amrex::Real(1.1);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

void uniform_state (gust::State& q, const PCase& c, double U)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
    q.u.setVal(static_cast<Real>(U));
    q.v.setVal(0.0);
    q.w.setVal(0.0);
    q.theta.setVal(300.0);
    q.tke.setVal(0.1);
}

double at00 (const amrex::MultiFab& mf, int comp = 0)
{
    double v = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(amrex::IntVect(0, 0, 0))) { v = mf[mfi](amrex::IntVect(0, 0, 0), comp); }
    }
    amrex::ParallelDescriptor::ReduceRealSum(v);
    return v;
}

SurfaceLayerOptions temperature (double theta_s)
{
    SurfaceLayerOptions o;
    o.z0 = Real(0.1);
    o.stability = SurfaceLayerOptions::Stability::surface_temperature;
    o.surface_temperature = static_cast<Real>(theta_s);
    return o;
}

SurfaceLayerOptions flux (double qs)
{
    SurfaceLayerOptions o;
    o.z0 = Real(0.1);
    o.stability = SurfaceLayerOptions::Stability::heat_flux;
    o.heat_flux = static_cast<Real>(qs);
    return o;
}

SurfaceLayerOptions with_wstar (SurfaceLayerOptions o, double zi)
{
    o.wstar = SurfaceLayerOptions::Wstar::beljaars;
    o.pblh = static_cast<Real>(zi);
    return o;
}

// The column with w* at wind U against the same column without w* at the wind U_eff (its speed): u*, q, 1/L, k_w and
// the drag coefficient agree; the stress is drag U along x (the wind's own), so tau(w*) = tau(no w*) U / U_eff.
void expect_wind_of_ueff (const PCase& c, const SurfaceLayerOptions& with, double U, const std::string& tag)
{
    gust::State q;
    uniform_state(q, c, U);
    gust::SurfaceState a;
    gust::compute_surface_state(q, with, kBuoyancy, c.mesh, c.geom, a);
    const double u_eff = at00(a.speed);
    EXPECT_GT(u_eff, U) << tag;
    SurfaceLayerOptions without = with;
    without.wstar = SurfaceLayerOptions::Wstar::none;
    without.pblh = 0;
    gust::State qe;
    uniform_state(qe, c, u_eff);
    gust::SurfaceState b;
    gust::compute_surface_state(qe, without, kBuoyancy, c.mesh, c.geom, b);
    const double r = tol(1e-12, 1e-5);
    EXPECT_NEAR(at00(b.speed), u_eff, r * u_eff) << tag;
    EXPECT_NEAR(at00(a.ustar), at00(b.ustar), r * at00(b.ustar)) << tag;
    EXPECT_NEAR(at00(a.qflux), at00(b.qflux), r * std::abs(at00(b.qflux))) << tag;
    EXPECT_NEAR(at00(a.inv_L), at00(b.inv_L), r * std::abs(at00(b.inv_L))) << tag;
    EXPECT_NEAR(at00(a.kwall), at00(b.kwall), r * at00(b.kwall)) << tag;
    EXPECT_NEAR(at00(a.drag), at00(b.drag), r * at00(b.drag)) << tag;
    EXPECT_NEAR(at00(a.tau, 0), at00(a.drag) * U, r * at00(a.drag) * U) << tag;
    EXPECT_NEAR(at00(a.tau, 0), at00(b.tau, 0) * U / u_eff, r * at00(a.tau, 0)) << tag;
    EXPECT_EQ(at00(a.pblh), static_cast<double>(static_cast<Real>(with.pblh))) << tag;
    EXPECT_EQ(at00(b.pblh), 0.0) << tag;
}

TEST(SurfaceWstar, WallModelTakesTheWindOfUeff)
{
    PCase c = flat_column(24);
    // A given flux: U_eff = sqrt(3^2 + (1.2 w*)^2) = 3.7464577 m/s for q = 0.2 K m/s, z_i = 1000 m (ClosedForms).
    {
        gust::State q;
        uniform_state(q, c, 3.0);
        gust::SurfaceState s;
        gust::compute_surface_state(q, with_wstar(flux(0.2), 1000), kBuoyancy, c.mesh, c.geom, s);
        EXPECT_NEAR(at00(s.speed), 3.7464577, 1e-6);
        EXPECT_EQ(at00(s.qflux), static_cast<double>(Real(0.2)));
    }
    expect_wind_of_ueff(c, with_wstar(flux(0.2), 1000), 3.0, "heat_flux 0.2, z_i 1000, U 3");
    expect_wind_of_ueff(c, with_wstar(flux(0.05), 400), 0.5, "heat_flux 0.05, z_i 400, U 0.5");
    // A given ground temperature, 3 K above the air: U = 1 m/s, z_i = 800 m gives U_eff = 2.0396978 m/s
    // (FixedPointMatchesTheBisection's pinned case, the same 5 m first centre over z0 = 0.1 m).
    {
        gust::State q;
        uniform_state(q, c, 1.0);
        gust::SurfaceState s;
        gust::compute_surface_state(q, with_wstar(temperature(303), 800), kBuoyancy, c.mesh, c.geom, s);
        EXPECT_NEAR(at00(s.speed), 2.0396978, tol(1e-6, 1e-4));
        EXPECT_NEAR(at00(s.qflux), 0.12429, tol(1e-5, 1e-4));
    }
    expect_wind_of_ueff(c, with_wstar(temperature(303), 800), 1.0, "theta_s 303, z_i 800, U 1");
    expect_wind_of_ueff(c, with_wstar(temperature(305), 1900), 3.0, "theta_s 305, z_i 1900, U 3");
    SurfaceLayerOptions kb = with_wstar(temperature(304), 1200);
    kb.kb_inverse = 2;
    expect_wind_of_ueff(c, kb, 2.0, "theta_s 304, kB 2, z_i 1200, U 2");
    // Jimenez's unstable side (a convective blend, unlike Businger-Dyer's): the kernel's own fixed point.
    SurfaceLayerOptions jz = with_wstar(temperature(303), 800);
    jz.stable_functions = SurfaceLayerOptions::StableFunctions::jimenez;
    expect_wind_of_ueff(c, jz, 1.0, "jimenez, theta_s 303, z_i 800, U 1");
    {
        gust::State q;
        uniform_state(q, c, 1.0);
        gust::SurfaceState s;
        gust::compute_surface_state(q, jz, kBuoyancy, c.mesh, c.geom, s);
        const double lz = std::log((kD + kZ0) / kZ0);
        Real u_eff = 0;
        const auto r = gust::surface_temperature_most_wstar<true>(
            Real(kKappa), Real(1), Real(kD), Real(lz), Real(1), Real(kGoth), Real(-3), gust::kSurfaceTemperatureRibMin,
            gust::kSurfaceTemperatureRibMaxJimenez, Real(0), Real(800), u_eff);
        EXPECT_NEAR(at00(s.speed), static_cast<double>(u_eff), tol(1e-12, 1e-5) * u_eff);
        EXPECT_NEAR(at00(s.qflux), static_cast<double>(r.q), tol(1e-12, 1e-5) * r.q);
    }
    // A stable column: w* = 0, the state of the run without w*, bit for bit.
    {
        gust::State q;
        uniform_state(q, c, 2.0);
        gust::SurfaceState a, b;
        gust::compute_surface_state(q, with_wstar(temperature(298), 800), kBuoyancy, c.mesh, c.geom, a);
        gust::compute_surface_state(q, temperature(298), kBuoyancy, c.mesh, c.geom, b);
        EXPECT_EQ(at00(a.speed), 2.0);
        EXPECT_EQ(at00(a.ustar), at00(b.ustar));
        EXPECT_EQ(at00(a.qflux), at00(b.qflux));
    }
}

// bc.surface_pblh_file: each column's z_i, read as the temperature map on any decomposition, reaches the state.
TEST(SurfaceWstar, PblhMapReachesTheState)
{
    PCase c = flat_column(12, 8, 4, false);
    const fs::path dir = fs::temp_directory_path()
                         / ("gust_pblh_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / "pblh.xyz").string();
    auto value = [] (double x, double y) { return 600 + 2 * x + y; };
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream out(file);
        out.precision(17);
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i < 8; ++i) {
                const int x = 25 + 50 * i, y = 25 + 50 * j;
                out << x << " " << y << " " << value(x, y) << "\n";
            }
        }
    }
    amrex::ParallelDescriptor::Barrier();
    amrex::MultiFab map;
    gust::read_pblh_map(file, "bc.surface_pblh_file", 4, 2.0, c.geom, c.ba, c.dm, map);
    EXPECT_GT(c.ba.size(), 3);
    EXPECT_NE(gust::surface_temperature_summary(map, "m").find("min 675, max 1725"), std::string::npos)
        << gust::surface_temperature_summary(map, "m");
    gust::State q;
    uniform_state(q, c, 1.0);
    SurfaceLayerOptions o = with_wstar(temperature(303), 0);
    o.pblh_file = file;
    o.pblh_map = &map;
    EXPECT_TRUE(gust::check_surface_layer_options(o).empty());
    gust::SurfaceState s;
    gust::compute_surface_state(q, o, kBuoyancy, c.mesh, c.geom, s);
    amrex::Long bad = 0;
    for (amrex::MFIter mfi(s.pblh, false); mfi.isValid(); ++mfi) {
        const auto zi = s.pblh.const_array(mfi);
        const auto sp = s.speed.const_array(mfi);
        amrex::LoopOnCpu(amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0)), [&] (int i, int j, int k) {
            const int ii = std::min(std::max(i, 0), 7), jj = std::min(std::max(j, 0), 7);
            const double want = value(25 + 50 * ii, 25 + 50 * jj);
            if (std::abs(zi(i, j, k) - want) > tol(1e-12, 1e-5) * want) { ++bad; }
            double u_eff = 0;
            solve_wstar(1, -3, want, u_eff);
            if (std::abs(sp(i, j, k) - u_eff) > tol(1e-12, 1e-5) * u_eff) { ++bad; }
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(bad);
    EXPECT_EQ(bad, 0);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { fs::remove_all(dir); }
}

// The options: beljaars needs heat_flux > 0 or surface_temperature and exactly one of z_i's value and map, the value
// in [10, 10000] m; neither z_i key without beljaars.
TEST(SurfaceWstar, OptionsAreChecked)
{
    auto first = [] (const SurfaceLayerOptions& o) {
        const auto e = gust::check_surface_layer_options(o);
        return e.empty() ? std::string() : e.front();
    };
    EXPECT_EQ(first(with_wstar(temperature(303), 800)), "");
    EXPECT_EQ(first(with_wstar(flux(0.1), 800)), "");
    EXPECT_EQ(first(with_wstar(temperature(303), 10)), "");
    EXPECT_EQ(first(with_wstar(temperature(303), 10000)), "");
    EXPECT_NE(first(with_wstar(temperature(303), 9)).find("[10, 10000]"), std::string::npos);
    EXPECT_NE(first(with_wstar(temperature(303), 10001)).find("[10, 10000]"), std::string::npos);
    {
        const gust::test::ScopedNoInvalidTrap no_trap;   // a NaN in the check's comparisons
        EXPECT_NE(first(with_wstar(temperature(303), std::numeric_limits<double>::quiet_NaN())).find("[10, 10000]"),
                  std::string::npos);
    }
    EXPECT_NE(first(with_wstar(temperature(303), 0)).find("exactly one"), std::string::npos);
    SurfaceLayerOptions both = with_wstar(temperature(303), 800);
    both.pblh_file = "z.xyz";
    EXPECT_NE(first(both).find("exactly one"), std::string::npos);
    EXPECT_NE(first(with_wstar(flux(0), 800)).find("<= 0"), std::string::npos);
    EXPECT_NE(first(with_wstar(flux(-0.02), 800)).find("<= 0"), std::string::npos);
    SurfaceLayerOptions neutral = with_wstar(SurfaceLayerOptions{}, 800);
    EXPECT_NE(first(neutral).find("wstar = beljaars needs stability"), std::string::npos);
    SurfaceLayerOptions zi_only = temperature(303);
    zi_only.pblh = 800;
    EXPECT_NE(first(zi_only).find("pblh is set but wstar is none"), std::string::npos);
    zi_only.pblh = 0;
    zi_only.pblh_file = "z.xyz";
    EXPECT_NE(first(zi_only).find("pblh_file is set but wstar is none"), std::string::npos);
}

} // namespace
