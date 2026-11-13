// The time driver of dual time stepping (TimeDriver.H), plan rev. 6
// sub-phase 6.2 (beyond the plan): Stokes' first problem against the erfc
// profile with first order in dt (UT-6.2.1), a steady problem reached and
// then held (UT-6.2.2), the step from a target CFL, the last step clipped to
// the stop time, steps that reach unsteady.max_inner, and the option checks.
#include "Driver.H"
#include "LaminarTestHelpers.H"
#include "Simplec.H"
#include "TimeDriver.H"
#include "TimeTerm.H"
#include "TransportTestHelpers.H"

#include <AMReX_Math.H>
#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <vector>

namespace {

using gust::FaceType;
using gust::lamtest::define_state;
using gust::lamtest::flat_case;
using gust::lamtest::Laminar;
using gust::ptest::PCase;
using gust::trtest::tol;

constexpr double kNu = 0.01;
constexpr double kLid = 1.0;

// The plane average of u at each cell layer k (u is uniform in x and y here).
std::vector<double> u_profile (const gust::State& q, int nz)
{
    std::vector<double> sum(nz, 0), cnt(nz, 0);
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        amrex::FArrayBox f(q.u[mfi].box(), 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(q.u[mfi], q.u[mfi].box());
        amrex::Gpu::streamSynchronize();
        const auto a = f.const_array();
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            sum[k] += a(i, j, k);
            cnt[k] += 1;
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(sum.data(), nz);
    amrex::ParallelDescriptor::ReduceRealSum(cnt.data(), nz);
    for (int k = 0; k < nz; ++k) { sum[k] /= cnt[k]; }
    return sum;
}

double max_diff (const std::vector<double>& a, const std::vector<double>& b)
{
    double m = 0;
    for (std::size_t k = 0; k < a.size(); ++k) { m = std::max(m, std::abs(a[k] - b[k])); }
    return m;
}

// u on the Couette profile U (k + 1/2) / nz. A free function: nvcc does not
// take an extended device lambda inside a constructor.
void set_couette (gust::State& q, int nz)
{
    for (amrex::MFIter mfi(q.u, false); mfi.isValid(); ++mfi) {
        const auto u = q.u.array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) {
            u(i, j, k) = amrex::Real(kLid) * (amrex::Real(k) + amrex::Real(0.5)) / amrex::Real(nz);
        });
    }
}

// A channel 1 m deep, periodic in x and y, ground at rest, top lid at kLid
// (laminar, nu = kNu), its state at rest or on the Couette profile, and the
// time driver over it.
struct Channel
{
    static constexpr int nz = 64;
    PCase c;
    Laminar lam;
    gust::State q;
    std::unique_ptr<gust::TimeTerm> term;
    std::unique_ptr<gust::Driver> driver;
    std::unique_ptr<gust::TimeDriver> td;

    // bump: a periodic ground h = 0.05 sin(2 pi x) m under the layers instead of flat ground.
    // tweak: changes to the driver's options after the defaults here.
    explicit Channel (gust::UnsteadyOptions uo, double residual_tol, bool couette = false, bool bump = false,
                      const std::function<void (gust::SteadyOptions&)>& tweak = {})
    {
        c = bump ? flat_case(4, 4, nz, 1.0, 1.0, 1.0, true, true, 2,
                             [] (double x, double) { return 0.05 * std::sin(2 * amrex::Math::pi<double>() * x); })
                 : flat_case(4, 4, nz, 1.0, 1.0, 1.0, true, true, 2);
        lam.bin.wall_velocity[5] = {amrex::Real(kLid), 0, 0};
        lam.make(c, kNu, FaceType::no_slip_wall, FaceType::no_slip_wall);
        lam.in.jvp_scales.velocity = 1;
        lam.in.jvp_scales.pressure = 1;
        define_state(q, c);
        q.theta.setVal(300.0);
        q.tke.setVal(1.0);
        if (couette) { set_couette(q, nz); }
        term = std::make_unique<gust::TimeTerm>(gust::TimeTermOptions{uo.cfl > 0 ? amrex::Real(1) : uo.dt}, c.ba, c.dm);
        lam.in.time = term.get();
        gust::SteadyOptions so;
        so.stage = gust::StageKind::unsteady;
        so.cfl0 = 1;
        so.cfl_max = 1000;
        so.zline_sweeps = 4;
        so.criterion.residual_tol = static_cast<amrex::Real>(residual_tol);
        if (tweak) { tweak(so); }
        gust::SimplecOptions opt;
        auto stage = std::make_unique<gust::SimplecStage>(opt, so.zline_sweeps, lam.in, nullptr, lam.popt, c.mesh, c.geom,
                                                              c.z, 0, c.ba, c.dm);
        const gust::SimplecStage& simplec = *stage;
        std::vector<std::unique_ptr<gust::Stage>> stages;
        stages.push_back(std::move(stage));
        driver = std::make_unique<gust::Driver>(so, std::move(stages), c.mesh, c.geom, &simplec.nu());
        gust::BoundaryConditions* bc = lam.bc.get();
        td = std::make_unique<gust::TimeDriver>(uo, *term, *driver, c.mesh, c.geom, amrex::Real(residual_tol),
                                                [bc] (gust::State& s) {
                                                    gust::BuoyancyOptions buoy;
                                                    buoy.gravity = 0;
                                                    bc->fill_ghosts(s, nullptr, buoy);
                                                });
    }

    void run ()
    {
        while (!td->finished()) { td->advance(q); }
    }
};

gust::UnsteadyOptions fixed_dt (double stop, double dt, double drop, int max_inner)
{
    gust::UnsteadyOptions uo;
    uo.stop_time = amrex::Real(stop);
    uo.dt = amrex::Real(dt);
    uo.inner_drop_orders = amrex::Real(drop);
    uo.max_inner = max_inner;
    uo.given = {"stop_time", "dt", "inner_drop_orders", "max_inner"};
    return uo;
}

} // namespace

// UT-6.2.1: the impulsively started lid over fluid at rest is Stokes' first
// problem while the layer it sets moving is thin next to the depth: u =
// U erfc((H - z) / (2 sqrt(nu t))). Backward Euler is first order in dt: the
// differences between runs at dt, dt/2 and dt/4 halve.
TEST(TimeDriver, StokesFirstProblemFirstOrder)
{
    if (!std::is_same_v<amrex::Real, double>) { GTEST_SKIP() << "the order needs a tight inner solve (double)"; }
    const double T = 1.0;
    std::vector<std::vector<double>> u;
    std::vector<long> inner;
    for (const double dt : {0.1, 0.05, 0.025}) {
        Channel ch(fixed_dt(T, dt, 8, 400), 1e-11);
        ch.run();
        ASSERT_FALSE(ch.td->diverged());
        EXPECT_EQ(ch.td->unconverged_steps(), 0) << "dt " << dt;
        EXPECT_EQ(ch.td->records().size(), static_cast<std::size_t>(std::lround(T / dt))) << "dt " << dt;
        EXPECT_NEAR(ch.td->time(), T, 1e-12);
        u.push_back(u_profile(ch.q, Channel::nz));
        inner.push_back(ch.td->inner_iterations());
    }
    const double d1 = max_diff(u[0], u[1]);
    const double d2 = max_diff(u[1], u[2]);
    const double order = std::log2(d1 / d2);
    EXPECT_GT(order, 0.85) << d1 << " " << d2;
    EXPECT_LT(order, 1.15) << d1 << " " << d2;
    // Against erfc at dt/4: the time error (about d2) plus the second-order space error.
    std::vector<double> ex(Channel::nz);
    for (int k = 0; k < Channel::nz; ++k) {
        const double z = (k + 0.5) / Channel::nz;
        ex[k] = kLid * std::erfc((1.0 - z) / (2 * std::sqrt(kNu * T)));
    }
    const double err = max_diff(u[2], ex);
    std::vector<double> rich(Channel::nz);
    for (int k = 0; k < Channel::nz; ++k) { rich[k] = 2 * u[2][k] - u[1][k]; }
    const double err_rich = max_diff(rich, ex);
    amrex::Print() << "Stokes first problem: differences " << d1 << " " << d2 << ", order " << order
                   << ", error against erfc at dt/4 " << err << ", Richardson " << err_rich << ", inner iterations " << inner[0] << " " << inner[1]
                   << " " << inner[2] << "\n";
    // Backward Euler's error at dt/4 is about C dt/4 = d2; with the space error (measured below) err is about d2 +
    // 8e-4 (0.0042 against d2 = 0.0035). Richardson's 2 u(dt/4) - u(dt/2) removes the first-order time error and
    // leaves the space error and O(dt^2) (measured 8.2e-4; 4.2e-3 without the extrapolation).
    EXPECT_LT(err, 1.5 * d2) << "d2 " << d2;
    EXPECT_LT(err_rich, 1.5e-3);
    EXPECT_LT(err_rich, 0.5 * err);
    EXPECT_GT(d2, 1e-3);   // the time error is resolved (not round-off)
}

// UT-6.2.2: a long run reaches the steady Couette profile u = U z / H within
// the inner tolerance, and a run started on it holds it, one inner iteration
// a step (the residual at q^n is already below the floor).
TEST(TimeDriver, SteadyProblemReachedAndHeld)
{
    {
        Channel ch(fixed_dt(300, 10, 6, 400), tol(1e-10, 1e-5));
        ch.run();
        ASSERT_FALSE(ch.td->diverged());
        EXPECT_EQ(ch.td->unconverged_steps(), 0);
        const auto u = u_profile(ch.q, Channel::nz);
        double worst = 0;
        for (int k = 0; k < Channel::nz; ++k) {
            worst = std::max(worst, std::abs(u[k] - kLid * (k + 0.5) / Channel::nz));
        }
        EXPECT_LT(worst, tol(1e-8, 1e-3));
        amrex::Print() << "Couette start-up: worst |u - U z/H| " << worst << ", " << ch.td->inner_iterations()
                       << " inner iterations over 30 steps\n";
    }
    {
        Channel ch(fixed_dt(30, 10, 6, 400), tol(1e-10, 1e-5), true);
        const auto u0 = u_profile(ch.q, Channel::nz);
        ch.run();
        ASSERT_EQ(ch.td->records().size(), 3U);
        for (const auto& s : ch.td->records()) {
            EXPECT_EQ(s.inner, 1) << "step " << s.step;
            EXPECT_TRUE(s.converged);
        }
        EXPECT_LT(max_diff(u_profile(ch.q, Channel::nz), u0), tol(1e-10, 1e-5));
    }
}

// The last step ends at stop_time: 0.25 in steps of 0.1 is 0.1, 0.1, 0.05.
TEST(TimeDriver, LastStepEndsAtStopTime)
{
    Channel ch(fixed_dt(0.25, 0.1, 2, 50), tol(1e-10, 1e-5));
    ch.run();
    const auto& r = ch.td->records();
    ASSERT_EQ(r.size(), 3U);
    EXPECT_NEAR(r[0].dt, 0.1, tol(1e-12, 1e-6));
    EXPECT_NEAR(r[2].dt, 0.05, tol(1e-12, 1e-6));
    EXPECT_EQ(ch.td->time(), amrex::Real(0.25));
    EXPECT_EQ(ch.term->dt(), r[2].dt);
    EXPECT_EQ(ch.term->inverse_dt(), amrex::Real(1) / r[2].dt);
    EXPECT_TRUE(ch.td->finished());
}

// With a target CFL the step is cfl over the largest convective rate; on the
// Couette profile (w = 0) that is |u| / dx in the top layer: dt = cfl dx / u_top.
TEST(TimeDriver, StepFromTargetCfl)
{
    gust::UnsteadyOptions uo;
    uo.stop_time = 1000;
    uo.cfl = amrex::Real(0.5);
    uo.given = {"stop_time", "cfl"};
    Channel ch(uo, tol(1e-10, 1e-5), true);
    const double dx = 1.0 / 4;
    const double u_top = kLid * (Channel::nz - 0.5) / Channel::nz;
    // advance fills the ghost cells itself before it takes dt (the state's ghosts are not filled here)
    const auto& s = ch.td->advance(ch.q);
    EXPECT_NEAR(s.dt, 0.5 * dx / u_top, tol(1e-12, 1e-6));
    const double dt = static_cast<double>(gust::convective_dt(ch.q, ch.c.mesh, ch.c.geom, uo.cfl));
    EXPECT_NEAR(dt, 0.5 * dx / u_top, tol(1e-12, 1e-6));
    EXPECT_EQ(ch.term->dt(), s.dt);
    EXPECT_EQ(ch.term->inverse_dt(), amrex::Real(1) / s.dt);
}

// Over terrain the contravariant velocity of the top and bottom faces reads u
// and v in the ghost cells, so the step from a target CFL needs them filled:
// advance fills them before it takes dt. Ghost values of 10 m/s planted in u
// would otherwise set the step.
TEST(TimeDriver, StepFromTargetCflFillsTheGhostsFirst)
{
    gust::UnsteadyOptions uo;
    uo.stop_time = 1000;
    uo.cfl = amrex::Real(0.5);
    uo.given = {"stop_time", "cfl"};
    Channel ch(uo, tol(1e-10, 1e-5), true, true);
    gust::State filled;
    define_state(filled, ch.c);
    for (auto* mf : {&filled.u, &filled.v, &filled.w, &filled.theta, &filled.tke}) { mf->setVal(0); }
    amrex::MultiFab::Copy(filled.u, ch.q.u, 0, 0, 1, 0);
    filled.theta.setVal(300.0);
    filled.tke.setVal(1.0);
    gust::BuoyancyOptions buoy;
    buoy.gravity = 0;
    ch.lam.bc->fill_ghosts(filled, nullptr, buoy);
    const amrex::Real want = gust::convective_dt(filled, ch.c.mesh, ch.c.geom, uo.cfl);
    ch.q.u.setBndry(10);
    ch.q.v.setBndry(10);
    const amrex::Real stale = gust::convective_dt(ch.q, ch.c.mesh, ch.c.geom, uo.cfl);
    EXPECT_LT(stale, amrex::Real(0.9) * want);   // the planted ghosts do set the step when they are not refilled
    const auto& s = ch.td->advance(ch.q);
    EXPECT_NEAR(s.dt, want, tol(1e-12, 1e-6) * want);
}

// The inner stop is relative to the residual at the step's start when 10^-d of
// it is above steady.residual_tol: with d = 2 and a floor far below, every step
// ends with its residual at most 1/100 of its first value, in fewer inner
// iterations than with d = 4.
TEST(TimeDriver, InnerStopIsRelativeToTheStepStart)
{
    Channel two(fixed_dt(0.3, 0.1, 2, 200), 1e-30);
    two.run();
    Channel four(fixed_dt(0.3, 0.1, 4, 200), 1e-30);
    four.run();
    ASSERT_EQ(two.td->records().size(), 3U);
    for (const auto& s : two.td->records()) {
        EXPECT_TRUE(s.converged) << "step " << s.step;
        EXPECT_LE(s.inner_residual, 1e-2) << "step " << s.step;
        EXPECT_GT(s.inner, 1) << "step " << s.step;
    }
    EXPECT_EQ(four.td->unconverged_steps(), 0);
    EXPECT_LT(two.td->inner_iterations(), four.td->inner_iterations());
}

// A step that reaches unsteady.max_inner is reported and the run goes on.
TEST(TimeDriver, MaxInnerIsReported)
{
    Channel ch(fixed_dt(0.2, 0.1, 12, 2), 1e-14);
    ch.run();
    ASSERT_EQ(ch.td->records().size(), 2U);
    EXPECT_EQ(ch.td->unconverged_steps(), 2);
    for (const auto& s : ch.td->records()) {
        EXPECT_EQ(s.inner, 2);
        EXPECT_FALSE(s.converged);
        EXPECT_NE(gust::format_step(s).find("unsteady.max_inner reached"), std::string::npos);
    }
    EXPECT_EQ(ch.td->inner_iterations(), 4);
}

TEST(TimeDriver, OptionsNameTheInput)
{
    const auto expect_one = [] (const gust::UnsteadyOptions& o, const std::string& what) {
        const auto err = gust::check_unsteady_options(o);
        ASSERT_FALSE(err.empty()) << what;
        bool found = false;
        for (const auto& e : err) { found = found || e.find(what) != std::string::npos; }
        EXPECT_TRUE(found) << what << ": " << err.front();
    };
    const gust::UnsteadyOptions good = fixed_dt(1, 0.1, 2, 50);
    EXPECT_TRUE(gust::check_unsteady_options(good).empty());
    {
        auto o = good;
        o.given = {"dt"};
        expect_one(o, "unsteady.stop_time is required");
    }
    for (const double t : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        auto o = good;
        o.stop_time = amrex::Real(t);
        expect_one(o, "unsteady.stop_time");
    }
    {
        auto o = good;
        o.given = {"stop_time"};
        expect_one(o, "unsteady.dt is required");
    }
    {
        auto o = good;
        o.dt = 0;
        expect_one(o, "unsteady.dt");
    }
    {
        auto o = good;
        o.cfl = amrex::Real(0.5);
        o.given.emplace_back("cfl");
        expect_one(o, "exclusive");
    }
    for (const double c : {-1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        auto o = good;
        o.cfl = amrex::Real(c);
        expect_one(o, "unsteady.cfl");
    }
    for (const double d : {0.0, 16.0, std::numeric_limits<double>::quiet_NaN()}) {
        auto o = good;
        o.inner_drop_orders = amrex::Real(d);
        expect_one(o, "unsteady.inner_drop_orders");
    }
    {
        auto o = good;
        o.max_inner = 0;
        expect_one(o, "unsteady.max_inner");
    }
    {
        gust::UnsteadyOptions o;
        o.stop_time = 1;
        o.cfl = amrex::Real(0.5);
        o.given = {"stop_time", "cfl"};
        EXPECT_TRUE(gust::check_unsteady_options(o).empty());
    }
}

// The time is n dt for a fixed step (not a sum of rounded steps: in double
// nine sums of 0.1 are 0.8999999999999999, 9 x 0.1 is 0.9), and a last step
// within 0.1 % of the time left takes all of it: 1.00005 in steps of 0.1 is
// ten steps, the last 0.10005, not eleven with a sliver of 5e-5.
TEST(TimeDriver, TimeIsStepCountTimesDtAndNoSliver)
{
    Channel ch(fixed_dt(1.00005, 0.1, 2, 1), tol(1e-10, 1e-5));
    ch.run();
    const auto& r = ch.td->records();
    ASSERT_EQ(r.size(), 10U);
    for (std::size_t n = 0; n + 1 < r.size(); ++n) {
        EXPECT_EQ(r[n].time, amrex::Real(double(n + 1) * double(amrex::Real(0.1)))) << "step " << n + 1;
    }
    EXPECT_NEAR(r.back().dt, 0.10005, tol(1e-12, 1e-6));
    EXPECT_EQ(ch.td->time(), amrex::Real(1.00005));
    EXPECT_TRUE(gust::is_last_step(0.1001, 0.1));
    EXPECT_FALSE(gust::is_last_step(0.10011, 0.1));
}

// A step that is not finite, and a target CFL over a field at rest (the
// largest rate at the floor kMinRate: dt = cfl 1e12), are refused by name
// instead of running the whole time in one step.
TEST(TimeDriver, StepDtChecks)
{
    EXPECT_TRUE(gust::check_step_dt(amrex::Real(0.1), 0).empty());
    EXPECT_TRUE(gust::check_step_dt(amrex::Real(0.1), amrex::Real(0.5)).empty());
    for (const double d : {0.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        EXPECT_NE(gust::check_step_dt(amrex::Real(d), 0).find("unsteady.dt"), std::string::npos) << d;
        EXPECT_NE(gust::check_step_dt(amrex::Real(d), amrex::Real(0.5)).find("unsteady.cfl"), std::string::npos) << d;
    }
    gust::UnsteadyOptions uo;
    uo.stop_time = 1;
    uo.cfl = amrex::Real(0.5);
    uo.given = {"stop_time", "cfl"};
    Channel ch(uo, tol(1e-10, 1e-5));   // at rest
    gust::BuoyancyOptions buoy;
    buoy.gravity = 0;
    ch.lam.bc->fill_ghosts(ch.q, nullptr, buoy);
    const amrex::Real dt = gust::convective_dt(ch.q, ch.c.mesh, ch.c.geom, uo.cfl);
    EXPECT_NE(gust::check_step_dt(dt, uo.cfl).find("at rest"), std::string::npos) << dt;
}

// With steady.stage = unsteady the driver ramps the pseudo-time CFL whatever
// its steady verdict says. The criterion here is met at once (a tolerance of
// 1e30, no drop, a window of 2), as in a run that settles or a restart from
// a steady checkpoint. With the ramp tied to a
// NotSteady verdict the CFL stays at cfl0 = 1 after the hold.
TEST(TimeDriver, CflRampsOnWhenTheVerdictIsSteady)
{
    Channel ch(fixed_dt(1, 0.1, 6, 400), tol(1e-10, 1e-5), false, false, [] (gust::SteadyOptions& so) {
        so.criterion.residual_tol = amrex::Real(1.0e30);
        so.criterion.residual_drop_orders = 0;
        so.criterion.monitor_rel_tol = amrex::Real(1.0e30);
        so.criterion.monitor_window = 2;
        so.ser_max_growth = amrex::Real(1.1);
    });
    ch.run();
    ASSERT_FALSE(ch.td->diverged());
    EXPECT_EQ(ch.driver->verdict().status, gust::SteadyStatus::Steady);
    EXPECT_GT(ch.driver->cfl(), amrex::Real(2)) << ch.td->inner_iterations() << " inner iterations";
}
