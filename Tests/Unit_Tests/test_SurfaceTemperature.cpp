// bc.surface_stability = surface_temperature (SurfaceLayer.H, SurfaceTemperature_K.H): the wall model's u* and q
// from a ground temperature match an independent bisection of the bulk-Richardson relation for stable and unstable
// air; a calm column keeps a finite flux at the bounded Ri_b; the theta rows take the flux's Picard coefficient;
// the map of bc.surface_temperature_file is read the same on any decomposition; the inputs are checked.
#include "ImplicitOperator.H"
#include "ERF_MOSTUtils.H"
#include "SurfaceLayer.H"
#include "SurfaceTemperature.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cfenv>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

using gust::Loc;
using gust::SurfaceLayerOptions;
using gust::ptest::PCase;
using gust::trtest::tol;

const gust::BuoyancyOptions kBuoyancy{};   // g = 9.81 m/s^2, theta_0 = 300 K

// A flat mesh, nx x nx x nz cells over 400 m x 400 m, geometric layers (dz0 = 10 m, ratio 1.1), boxes of mgs
// cells, periodic sides unless told otherwise.
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

void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

// A surface-state field at column (0, 0).
double at00 (const amrex::MultiFab& mf)
{
    double v = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(amrex::IntVect(0, 0, 0))) { v = mf[mfi](amrex::IntVect(0, 0, 0)); }
    }
    amrex::ParallelDescriptor::ReduceRealSum(v);
    return v;
}

struct Reference { double ustar, q, zeta, c; };

// The independent solve: bisection to round-off on zeta = Ri_b cm^2 / ch over [-100, 0] or [0, 100] with ERF's
// Businger-Dyer psi (extern_erf/ERF_MOSTUtils.H), Ri_b = g d dtheta / (theta_0 U^2) bounded to [-4, 0.19] (a
// bracket without a root gives its far end), cm = max(lz - psi_m, 1), ch = max(lz - psi_h, 1) (lz = 3.93 > 1 at the
// 5 m first centre over z0 = 0.1 m: both floors are 1). Then u* = kappa U / cm, q = -kappa u* dtheta / ch.
Reference reference_most (double U, double d, double z0, double dtheta)
{
    const auto prev_fpe = amrex::disableFPExcept(amrex::FPExcept::invalid);   // ERF's branches (test_SurfaceLayer.cpp)
    const similarity_funs sf;
    const double kappa = 0.41, lz = std::log((d + z0) / z0);
    const double rib = std::max(-4.0, std::min(0.19, 9.81 * d * dtheta / (300 * U * U)));
    auto cm = [&] (double z) { return std::max(lz - double(sf.calc_psi_m(static_cast<amrex::Real>(z))), 1.0); };
    auto ch = [&] (double z) { return std::max(lz - double(sf.calc_psi_h(static_cast<amrex::Real>(z))), 1.0); };
    auto G = [&] (double z) { return z - rib * cm(z) * cm(z) / ch(z); };
    double zeta = 0;
    if (rib != 0) {
        double a = rib < 0 ? -100 : 0, b = rib < 0 ? 0 : 100;
        if (rib > 0 && G(b) <= 0) { zeta = b; }
        else if (rib < 0 && G(a) >= 0) { zeta = a; }
        else {
            for (int it = 0; it < 300; ++it) {
                const double mid = 0.5 * (a + b);
                ((G(mid) < 0) ? a : b) = mid;
            }
            zeta = 0.5 * (a + b);
        }
    }
    const double us = kappa * U / cm(zeta);
    const double c = kappa * us / ch(zeta);
    std::feclearexcept(FE_ALL_EXCEPT);
    amrex::setFPExcept(prev_fpe);
    return {us, -c * dtheta, zeta, c};
}

SurfaceLayerOptions temperature_options (double theta_s)
{
    SurfaceLayerOptions opt;
    opt.z0 = amrex::Real(0.1);
    opt.stability = SurfaceLayerOptions::Stability::surface_temperature;
    opt.surface_temperature = static_cast<amrex::Real>(theta_s);
    return opt;
}

// A uniform wind U along x over theta = 300 K, every column.
void uniform_state (gust::State& q, double U)
{
    q.u.setVal(static_cast<amrex::Real>(U));
    q.v.setVal(0.0);
    q.w.setVal(0.0);
    q.theta.setVal(300.0);
    q.tke.setVal(0.1);
}

// u*, q, 1/L and the Picard coefficient against the bisection, at the 5 m first centre of a flat column: a warm
// ground (theta_s = 305 K under 300 K air, unstable) and a cool one (298 K, stable), at 8 and 3 m/s. The stable case
// at 3 m/s (Ri_b = 0.073, zeta = 0.5) is inside Businger-Dyer's root range. The flux is upward over a warm ground and
// downward over a cool one, and k_w carries Axell's buoyancy term only when q > 0.
TEST(SurfaceTemperature, MatchesTheBisectionStableAndUnstable)
{
    PCase c = flat_column(24);
    gust::State q;
    define_state(q, c);
    const double z0 = 0.1, d = 5.0;
    struct Case { double U, theta_s; };
    for (const Case cs : {Case{8.0, 305.0}, Case{8.0, 298.0}, Case{3.0, 305.0}, Case{3.0, 298.0}}) {
        uniform_state(q, cs.U);
        const SurfaceLayerOptions opt = temperature_options(cs.theta_s);
        gust::SurfaceState s;
        gust::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
        const Reference want = reference_most(cs.U, d, z0, 300.0 - cs.theta_s);
        const std::string tag = "U " + std::to_string(cs.U) + " theta_s " + std::to_string(cs.theta_s);
        EXPECT_NEAR(at00(s.ustar), want.ustar, tol(1e-11, 1e-4) * want.ustar) << tag;
        EXPECT_NEAR(at00(s.qflux), want.q, tol(1e-11, 1e-4) * std::abs(want.q)) << tag;
        EXPECT_NEAR(at00(s.inv_L) * d, want.zeta, tol(1e-10, 1e-4) * (1 + std::abs(want.zeta))) << tag;
        EXPECT_NEAR(at00(s.heat_coef), want.c, tol(1e-11, 1e-4) * want.c) << tag;   // flat: |A| = 1
        EXPECT_EQ(at00(s.theta_s), cs.theta_s) << tag;
        EXPECT_EQ(at00(s.qflux) > 0, cs.theta_s > 300) << tag;
        const double kw = std::pow(std::pow(want.ustar, 3) + 0.41 * std::max(9.81 / 300 * want.q, 0.0) * d, 2.0 / 3.0)
                          / (0.5562 * 0.5562);
        EXPECT_NEAR(at00(s.kwall), kw, tol(1e-11, 1e-4) * kw) << tag;
    }
    // The pinned numbers of the unstable 8 m/s case (from the bisection above, written out): a change in the
    // stability functions or in the formula shows here too.
    uniform_state(q, 8.0);
    gust::SurfaceState s;
    gust::compute_surface_state(q, temperature_options(305.0), kBuoyancy, c.mesh, c.geom, s);
    const Reference r = reference_most(8.0, d, z0, -5.0);
    EXPECT_LT(r.zeta, -0.05);
    EXPECT_GT(r.q, 0.05);
    EXPECT_NEAR(at00(s.qflux), r.q, tol(1e-11, 1e-4) * r.q);
}

// Zero wind: U is floored at 0.1 m/s, so Ri_b = g d dtheta / (theta_0 0.01) is -81.75 over a 5 K warmer ground and
// 32.7 over a 2 K cooler one. Bounded to [-4, 0.19], u* and q stay finite and equal the bisection's at the bound
// (without the bound, the stable column has no root below zeta = 100 and the unstable one a different root).
TEST(SurfaceTemperature, ZeroWindGivesTheBoundedFlux)
{
    PCase c = flat_column(24);
    gust::State q;
    define_state(q, c);
    for (const double theta_s : {305.0, 298.0}) {
        uniform_state(q, 0.0);
        gust::SurfaceState s;
        gust::compute_surface_state(q, temperature_options(theta_s), kBuoyancy, c.mesh, c.geom, s);
        const Reference want = reference_most(0.1, 5.0, 0.1, 300.0 - theta_s);
        EXPECT_TRUE(std::isfinite(at00(s.qflux)) && std::isfinite(at00(s.ustar))) << theta_s;
        EXPECT_NEAR(at00(s.ustar), want.ustar, tol(1e-11, 1e-4) * want.ustar) << theta_s;
        EXPECT_NEAR(at00(s.qflux), want.q, tol(1e-11, 1e-4) * std::abs(want.q)) << theta_s;
        // The bound acted: the root is inside the bracket (unbounded, the stable column ends at zeta = 100).
        EXPECT_NEAR(at00(s.inv_L) * 5.0, want.zeta, tol(1e-10, 1e-4) * (1 + std::abs(want.zeta))) << theta_s;
        EXPECT_GT(want.zeta, -99.0) << theta_s;
        EXPECT_LT(want.zeta, 99.0) << theta_s;
    }
}

// The theta rows of the first cells: apply_surface_temperature_to_implicit adds heat_coef / dz (dzeta^-1 / J = 1/dz
// on a flat mesh; dz = 10 m) to the diagonal and nothing else; for the other stabilities it adds nothing.
TEST(SurfaceTemperature, ImplicitThetaRowTakesThePicardCoefficient)
{
    PCase c = flat_column(24);
    gust::State q;
    define_state(q, c);
    uniform_state(q, 8.0);
    gust::SurfaceState s;
    const SurfaceLayerOptions opt = temperature_options(305.0);
    gust::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
    amrex::MultiFab A(c.ba, c.dm, gust::n_implicit_coefs, 0);
    A.setVal(1.0);
    gust::apply_surface_temperature_to_implicit(s, opt, c.mesh, c.geom, A);
    const double hc = at00(s.heat_coef);
    for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
        const auto a = A.const_array(mfi);
        const amrex::Box& b = mfi.validbox();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            for (int n = 0; n < gust::n_implicit_coefs; ++n) {
                const double want = (k == 0 && n == gust::coef_P) ? 1 + hc / 10.0 : 1.0;
                EXPECT_NEAR(a(i, j, k, n), want, tol(1e-13, 1e-6) * want) << i << " " << j << " " << k << " " << n;
            }
        });
    }
    SurfaceLayerOptions flux;
    flux.stability = SurfaceLayerOptions::Stability::heat_flux;
    flux.heat_flux = amrex::Real(0.1);
    A.setVal(1.0);
    gust::apply_surface_temperature_to_implicit(s, flux, c.mesh, c.geom, A);
    EXPECT_EQ(A.min(gust::coef_P), 1.0);
    EXPECT_EQ(A.max(gust::coef_P), 1.0);
}

// The map of bc.surface_temperature_file on an 8 x 8 mesh of four 4 x 4 boxes (seams between them), not periodic:
// the file holds one point at every cell centre with theta_s = 290 + 0.01 x + 0.02 y, so the map takes each centre's
// own value (an input point closer than 1e-6 m gives its value), and a ghost column the nearest domain column's.
// The unit tests run on 1, 2 and 4 ranks (unit_np1, unit_np2, unit_np4): every rank's columns match. The wall model
// then reads each column's theta_s.
TEST(SurfaceTemperature, MapIsReadTheSameOnEveryDecomposition)
{
    PCase c = flat_column(12, 8, 4, false);
    const fs::path dir = fs::temp_directory_path() / ("gust_tsurf_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / "theta_s.xyz").string();
    auto value = [] (double x, double y) { return 290 + 0.01 * x + 0.02 * y; };
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream out(file);
        out.precision(17);
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i < 8; ++i) {
                const double x = 25 + 50 * i, y = 25 + 50 * j;
                out << x << " " << y << " " << value(x, y) << "\n";
            }
        }
    }
    amrex::ParallelDescriptor::Barrier();
    amrex::MultiFab map;
    gust::read_surface_temperature_map(file, "bc.surface_temperature_file", 4, 2.0, c.geom, c.ba, c.dm, map);
    EXPECT_GT(c.ba.size(), 3);
    long bad = 0;
    for (amrex::MFIter mfi(map, false); mfi.isValid(); ++mfi) {
        const auto a = map.const_array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
            const int ii = std::min(std::max(i, 0), 7), jj = std::min(std::max(j, 0), 7);
            const double want = value(25 + 50 * ii, 25 + 50 * jj);
            if (std::abs(a(i, j, k) - want) > tol(1e-12, 1e-5) * want) { ++bad; }
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(bad);
    EXPECT_EQ(bad, 0);
    EXPECT_NE(gust::surface_temperature_summary(map).find("min 290.75"), std::string::npos)
        << gust::surface_temperature_summary(map);

    gust::State q;
    define_state(q, c);
    uniform_state(q, 8.0);
    SurfaceLayerOptions opt = temperature_options(0);
    opt.surface_temperature = 0;
    opt.surface_temperature_file = file;
    opt.theta_s_map = &map;
    gust::SurfaceState s;
    gust::compute_surface_state(q, opt, kBuoyancy, c.mesh, c.geom, s);
    bad = 0;
    for (amrex::MFIter mfi(s.theta_s, false); mfi.isValid(); ++mfi) {
        const auto ts = s.theta_s.const_array(mfi);
        const auto qf = s.qflux.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const double t = value(25 + 50 * i, 25 + 50 * j);
            const Reference r = reference_most(8.0, 5.0, 0.1, 300.0 - t);
            if (ts(i, j, k) != static_cast<amrex::Real>(t) && std::abs(ts(i, j, k) - t) > tol(1e-12, 1e-5) * t) { ++bad; }
            if (std::abs(qf(i, j, k) - r.q) > tol(1e-10, 1e-4) * std::abs(r.q)) { ++bad; }
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(bad);
    EXPECT_EQ(bad, 0);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { fs::remove_all(dir); }
}

// The options: exactly one of the value and the file; the value in [150, 350] K; neither without the stability;
// no qnse; no heat flux beside it.
TEST(SurfaceTemperature, OptionsAreChecked)
{
    auto first = [] (const SurfaceLayerOptions& o) {
        const auto e = gust::check_surface_layer_options(o);
        return e.empty() ? std::string() : e.front();
    };
    EXPECT_EQ(first(temperature_options(290)), "");
    EXPECT_NE(first(temperature_options(0)).find("exactly one"), std::string::npos);
    SurfaceLayerOptions both = temperature_options(290);
    both.surface_temperature_file = "t.xyz";
    EXPECT_NE(first(both).find("exactly one"), std::string::npos);
    EXPECT_NE(first(temperature_options(100)).find("[150, 350]"), std::string::npos);
    EXPECT_NE(first(temperature_options(400)).find("[150, 350]"), std::string::npos);
    EXPECT_NE(first(temperature_options(std::numeric_limits<double>::quiet_NaN())).find("[150, 350]"), std::string::npos);
    SurfaceLayerOptions qn = temperature_options(290);
    qn.stable_functions = SurfaceLayerOptions::StableFunctions::qnse;
    EXPECT_NE(first(qn).find("qnse"), std::string::npos);
    SurfaceLayerOptions hf = temperature_options(290);
    hf.heat_flux = amrex::Real(0.1);
    EXPECT_NE(first(hf).find("heat_flux"), std::string::npos);
    SurfaceLayerOptions neutral;
    neutral.surface_temperature = 290;
    EXPECT_NE(first(neutral).find("surface_temperature is set"), std::string::npos);
    neutral.surface_temperature = 0;
    neutral.surface_temperature_file = "t.xyz";
    EXPECT_NE(first(neutral).find("surface_temperature_file is set"), std::string::npos);
}

} // namespace
