// simplec.theta = relaxed_buoyant (ThetaRelax.H): the reader and its refusals,
// the source (theta - theta_p)/tau and the diagonal 1/tau against closed form
// (the kernels alone, and through residual() and implicit_rows()), and the two
// limits of the SIMPLEC stage on 4 boxes: tau -> 0 holds theta at theta_p and
// moves u, v, w, p as frozen_buoyant does; tau -> infinity moves them as
// theta solved does.
#include "ImplicitOperator.H"
#include "LaminarTestHelpers.H"
#include "MassFlux.H"
#include "PseudoTime.H"
#include "Residual.H"
#include "Simplec.H"
#include "ThetaRelax.H"

#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using gust::FaceType;
using gust::lamtest::define_state;
using gust::lamtest::flat_case;
using gust::lamtest::Laminar;
using gust::ptest::PCase;

// Random values on the valid points of every component, then the periodic images (seeded per box index, so the
// field does not depend on the rank count).
void fill_random (amrex::MultiFab& mf, double lo, double hi, unsigned seed, const PCase& c)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        std::mt19937 gen(seed + 7919u * static_cast<unsigned>(mfi.index()));
        std::uniform_real_distribution<double> dist(lo, hi);
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox f(b, mf.nComp(), amrex::The_Pinned_Arena());
        for (int n = 0; n < mf.nComp(); ++n) {
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k), n) = static_cast<amrex::Real>(dist(gen)); });
        }
        mf[mfi].copy<amrex::RunOn::Device>(f, b, 0, b, 0, mf.nComp());
        amrex::Gpu::streamSynchronize();
    }
    if (!mf.ixType().cellCentered()) { mf.OverrideSync(c.geom.periodicity()); }
    mf.FillBoundary(c.geom.periodicity());
}

double max_abs_diff (const amrex::MultiFab& a, const amrex::MultiFab& b, int comp = 0)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(d, a, comp, 0, 1, 0);
    amrex::MultiFab::Subtract(d, b, comp, 0, 1, 0);
    return static_cast<double>(d.norm0());
}

void copy_state (gust::State& to, const gust::State& from)
{
    for (auto [t, f] : {std::pair{&to.u, &from.u}, std::pair{&to.v, &from.v}, std::pair{&to.w, &from.w},
                        std::pair{&to.p, &from.p}, std::pair{&to.theta, &from.theta}, std::pair{&to.tke, &from.tke}}) {
        amrex::MultiFab::Copy(*t, *f, 0, 0, 1, 2);
    }
}

// A random state on c with theta far from the reference (its buoyancy force is of order g/theta_0 * 5 K).
void random_state (gust::State& q, const PCase& c)
{
    define_state(q, c);
    fill_random(q.u, 4, 8, 11, c);
    fill_random(q.v, -2, 2, 12, c);
    fill_random(q.w, -0.5, 0.5, 13, c);
    q.p.setVal(0.0);
    fill_random(q.theta, 295, 305, 14, c);
    q.tke.setVal(0.5);
}

// The laminar inputs with gravity on, periodic sides and slip walls at the ground and top.
void make_laminar (Laminar& lam, const PCase& c)
{
    lam.make(c, 5.0, FaceType::slip_wall, FaceType::slip_wall);
    lam.in.buoyancy.gravity = static_cast<amrex::Real>(9.81);
    lam.popt.reltol = amrex::Real(1e-10);
}

// n SIMPLEC iterations of the laminar stage with opt from q (in place).
void run_stage (gust::State& q, const gust::SimplecOptions& opt, const PCase& c, int n)
{
    Laminar lam;
    make_laminar(lam, c);
    lam.in.buoyancy_in_momentum = gust::buoyancy_in_momentum(opt);
    lam.in.buoyancy_point_implicit = gust::buoyancy_point_implicit(opt);
    gust::SimplecStage stage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm);
    stage.init(q);
    amrex::MultiFab dtau;
    for (int it = 0; it < n; ++it) {
        lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy);
        gust::local_dtau(q, &lam.nu, c.mesh, c.geom, amrex::Real(1), dtau);
        stage.iterate(q, dtau);
    }
}

gust::SimplecOptions relaxed (double tau)
{
    gust::SimplecOptions o;
    o.theta_relaxed = true;
    o.theta_relax_time = static_cast<amrex::Real>(tau);
    return o;
}

void clear_simplec ()
{
    amrex::ParmParse pp("simplec");
    while (pp.remove("theta") > 0) {}
    while (pp.remove("theta_relax_time") > 0) {}
}

void expect_refused (const std::string& key)
{
    try {
        (void)gust::read_simplec_options();
        ADD_FAILURE() << "accepted; expected a refusal naming " << key;
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find(key), std::string::npos) << e.what();
    }
}

} // namespace

// ---------------------------------------------------------------------------
// The reader: relaxed_buoyant with a time scale solves theta (not frozen) with the buoyancy in w and the w rows'
// point-implicit buoyancy on; without the time scale, with a time scale of 0 or < 0 (NaN and infinity: the range
// check alone), and with the
// time scale under any other simplec.theta, it aborts naming simplec.theta_relax_time.
TEST (ThetaRelax, ReaderAndRefusals)
{
    amrex::ParmParse pp("simplec");
    clear_simplec();
    pp.add("theta", std::string("relaxed_buoyant"));
    pp.add("theta_relax_time", amrex::Real(1000));
    const gust::SimplecOptions o = gust::read_simplec_options();
    EXPECT_TRUE(o.theta_relaxed);
    EXPECT_FALSE(o.theta_frozen);
    EXPECT_FALSE(o.theta_buoyant);
    EXPECT_EQ(o.theta_relax_time, amrex::Real(1000));
    EXPECT_EQ(gust::theta_value(o), "relaxed_buoyant");
    EXPECT_TRUE(gust::buoyancy_in_momentum(o));
    EXPECT_TRUE(gust::buoyancy_point_implicit(o));
    clear_simplec();

    pp.add("theta", std::string("relaxed_buoyant"));
    expect_refused("simplec.theta_relax_time");
    clear_simplec();
    for (const char* bad : {"0", "-5", "-1e-30"}) {
        pp.add("theta", std::string("relaxed_buoyant"));
        pp.add("theta_relax_time", std::string(bad));
        expect_refused(std::string("simplec.theta_relax_time = ") + bad + ":");   // the value as given (%g)
        clear_simplec();
    }
    for (const char* other : {"solve", "frozen", "frozen_buoyant"}) {
        pp.add("theta", std::string(other));
        pp.add("theta_relax_time", amrex::Real(1000));
        expect_refused("simplec.theta_relax_time");
        clear_simplec();
    }
    pp.add("theta_relax_time", amrex::Real(1000));   // with the default simplec.theta = solve
    expect_refused("simplec.theta_relax_time");
    clear_simplec();
    // The range check alone (the stage repeats it).
    gust::SimplecOptions r = relaxed(0);
    ASSERT_FALSE(gust::check_simplec_options(r).empty());
    EXPECT_NE(gust::check_simplec_options(r)[0].find("simplec.theta_relax_time"), std::string::npos);
    for (const amrex::Real bad : {std::numeric_limits<amrex::Real>::quiet_NaN(), std::numeric_limits<amrex::Real>::infinity()}) {
        const auto err = gust::check_simplec_options(relaxed(static_cast<double>(bad)));
        ASSERT_EQ(err.size(), 1u);
        EXPECT_NE(err[0].find("simplec.theta_relax_time"), std::string::npos) << err[0];
    }
    EXPECT_TRUE(gust::check_simplec_options(relaxed(300)).empty());
}

// ---------------------------------------------------------------------------
// The kernels against closed form on 4 boxes: R += (theta - theta_p)/tau (the expected value formed with MultiFab
// arithmetic, another path), and a_P += 1/tau with the six neighbour coefficients unchanged to the last bit.
TEST (ThetaRelax, SourceAndDiagonalClosedForm)
{
    PCase c = flat_case(8, 8, 8, 800, 800, 400, true, true, 4);
    const amrex::Real tau = 300;
    amrex::MultiFab theta(c.ba, c.dm, 1, 1), tp(c.ba, c.dm, 1, 0), R(c.ba, c.dm, 1, 0), R0(c.ba, c.dm, 1, 0);
    fill_random(theta, 295, 305, 21, c);
    fill_random(tp, 295, 305, 22, c);
    fill_random(R, -1, 1, 23, c);
    amrex::MultiFab::Copy(R0, R, 0, 0, 1, 0);
    gust::add_theta_relaxation(theta, tp, 1 / tau, R);
    amrex::MultiFab e(c.ba, c.dm, 1, 0);
    amrex::MultiFab::Copy(e, theta, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(e, tp, 0, 0, 1, 0);
    e.mult(1 / tau);
    amrex::MultiFab::Add(e, R0, 0, 0, 1, 0);
    const double eps = static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
    EXPECT_LE(max_abs_diff(R, e), 8 * eps);   // |R| <= 1.04: a few ulps (contracted multiply-adds)
    // Not vacuous: the source is of order 10 K / 300 s.
    EXPECT_GT(max_abs_diff(R, R0), 0.02);

    amrex::MultiFab A(c.ba, c.dm, gust::n_implicit_coefs, 0), A0(c.ba, c.dm, gust::n_implicit_coefs, 0);
    fill_random(A, 0, 2, 24, c);
    amrex::MultiFab::Copy(A0, A, 0, 0, gust::n_implicit_coefs, 0);
    gust::add_theta_relaxation_diagonal(1 / tau, A);
    amrex::MultiFab d(c.ba, c.dm, 1, 0);
    amrex::MultiFab::Copy(d, A, gust::coef_P, 0, 1, 0);
    amrex::MultiFab::Subtract(d, A0, gust::coef_P, 0, 1, 0);
    d.plus(-1 / tau, 0, 1, 0);
    EXPECT_LE(static_cast<double>(d.norm0()), 4 * eps);
    for (int n = 0; n < gust::n_implicit_coefs; ++n) {
        if (n != gust::coef_P) { EXPECT_EQ(max_abs_diff(A, A0, n), 0.0) << "coefficient " << n; }
    }
}

// ---------------------------------------------------------------------------
// The wiring on 4 boxes: with ResidualInputs::theta_relax set, residual()'s R_theta gains (theta - theta_p)/tau
// on every cell and R_u, R_v, R_w, R_p, R_k are unchanged to the last bit; implicit_rows() gains 1/tau on theta's
// a_P only (u's and k's rows unchanged to the last bit).
TEST (ThetaRelax, ResidualAndImplicitRows)
{
    PCase c = flat_case(8, 8, 8, 800, 800, 400, true, true, 4);
    Laminar lam;
    make_laminar(lam, c);
    gust::State q;
    random_state(q, c);
    lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy);
    const amrex::Real tau = 1000;
    gust::ThetaRelax tr(tau);
    amrex::MultiFab tp(c.ba, c.dm, 1, 0);
    fill_random(tp, 295, 305, 31, c);
    tr.set_target(tp);

    gust::ResidualFields off, on;
    off.define(c.ba, c.dm);
    on.define(c.ba, c.dm);
    gust::residual(q, lam.in, c.mesh, c.geom, off);
    gust::ResidualInputs relax_in = lam.in;
    relax_in.theta_relax = &tr;
    gust::residual(q, relax_in, c.mesh, c.geom, on);
    amrex::MultiFab e(c.ba, c.dm, 1, 0), d(c.ba, c.dm, 1, 0);
    amrex::MultiFab::Copy(e, q.theta, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(e, tp, 0, 0, 1, 0);
    e.mult(1 / tau);
    amrex::MultiFab::Copy(d, on.theta, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, off.theta, 0, 0, 1, 0);
    const double scale = static_cast<double>(amrex::max(off.theta.norm0(), e.norm0()));
    const double eps = static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
    EXPECT_LE(max_abs_diff(d, e), 16 * eps * scale);
    EXPECT_GT(static_cast<double>(e.norm0()), 1e-3);
    for (auto [a, b, name] : {std::tuple{&on.u, &off.u, "u"}, std::tuple{&on.v, &off.v, "v"}, std::tuple{&on.w, &off.w, "w"},
                              std::tuple{&on.p, &off.p, "p"}, std::tuple{&on.tke, &off.tke, "tke"}}) {
        EXPECT_EQ(max_abs_diff(*a, *b), 0.0) << name;
    }

    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;
    gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM> Fp{&F[0], &F[1], &F[2]};
    for (const gust::Equation eq : {gust::Equation::theta, gust::Equation::u, gust::Equation::tke}) {
        amrex::MultiFab A0, A1;
        gust::implicit_rows(Fp, q.z_phys_nd, lam.in, eq, c.mesh, c.geom, A0);
        gust::implicit_rows(Fp, q.z_phys_nd, relax_in, eq, c.mesh, c.geom, A1);
        const amrex::Real add = eq == gust::Equation::theta ? 1 / tau : amrex::Real(0);
        amrex::MultiFab dd(A0.boxArray(), A0.DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(dd, A1, gust::coef_P, 0, 1, 0);
        amrex::MultiFab::Subtract(dd, A0, gust::coef_P, 0, 1, 0);
        dd.plus(-add, 0, 1, 0);
        const double ap = static_cast<double>(A0.norm0(gust::coef_P));
        if (eq == gust::Equation::theta) {
            EXPECT_LE(static_cast<double>(dd.norm0()), 4 * eps * ap) << "theta's a_P";
        } else {
            EXPECT_EQ(static_cast<double>(dd.norm0()), 0.0) << "eq " << static_cast<int>(eq);
        }
        for (int n = 0; n < gust::n_implicit_coefs; ++n) {
            if (n != gust::coef_P) { EXPECT_EQ(max_abs_diff(A1, A0, n), 0.0) << "eq " << static_cast<int>(eq) << " coef " << n; }
        }
    }
}

// ---------------------------------------------------------------------------
// The limits of the SIMPLEC stage, 20 laminar iterations with gravity on from a random state whose theta is far from
// the reference, on 4 boxes:
//   tau = 1e-6 s: theta stays within 1e-4 K of its start (theta_p) and u, v, w, p are within 1e-5 of
//                 frozen_buoyant's (theta held to the last bit, buoyancy in w, no point-implicit buoyancy);
//   tau = 1e12 s: theta, u, v, w, p are within 1e-9 of theta solved's.
// Not vacuous: frozen_buoyant and solve differ by more than 0.01 in theta and w after these iterations, and
// tau = 100 s lies away from both. The mutants (the source's sign flipped, the diagonal left out) diverge at
// tau = 1e-6 s.
TEST (ThetaRelax, LimitsMatchFrozenBuoyantAndSolve)
{
    PCase c = flat_case(8, 8, 8, 800, 800, 400, true, true, 4);
    gust::State q0;
    random_state(q0, c);
    const int n = 20;
    auto run = [&] (const gust::SimplecOptions& opt, gust::State& q) {
        define_state(q, c);
        copy_state(q, q0);
        run_stage(q, opt, c, n);
    };
    gust::SimplecOptions fb;
    fb.theta_frozen = true;
    fb.theta_buoyant = true;
    gust::State frozen, solved, small, large, mid;
    run(fb, frozen);
    run(gust::SimplecOptions{}, solved);
    run(relaxed(1e-6), small);
    run(relaxed(1e12), large);
    run(relaxed(100), mid);

    EXPECT_EQ(max_abs_diff(frozen.theta, q0.theta), 0.0);
    EXPECT_GT(max_abs_diff(frozen.theta, solved.theta), 1e-2);
    EXPECT_GT(max_abs_diff(frozen.w, solved.w), 1e-2);

    EXPECT_LT(max_abs_diff(small.theta, q0.theta), 1e-4);
    for (auto [a, b, name] : {std::tuple{&small.u, &frozen.u, "u"}, std::tuple{&small.v, &frozen.v, "v"},
                              std::tuple{&small.w, &frozen.w, "w"}, std::tuple{&small.p, &frozen.p, "p"}}) {
        EXPECT_LT(max_abs_diff(*a, *b), 1e-5) << name << " (tau -> 0 against frozen_buoyant)";
    }
    for (auto [a, b, name] : {std::tuple{&large.theta, &solved.theta, "theta"}, std::tuple{&large.u, &solved.u, "u"},
                              std::tuple{&large.v, &solved.v, "v"}, std::tuple{&large.w, &solved.w, "w"},
                              std::tuple{&large.p, &solved.p, "p"}}) {
        EXPECT_LT(max_abs_diff(*a, *b), 1e-9) << name << " (tau -> infinity against solve)";
    }
    EXPECT_GT(max_abs_diff(mid.theta, frozen.theta), 1e-4);
    EXPECT_GT(max_abs_diff(mid.theta, solved.theta), 1e-4);
}
