// simplec.theta = frozen: the reader and its refusals, theta held to the last
// bit by the SIMPLEC stage, no buoyancy in the momentum residual or the
// predictor's w rows (the frozen stage with gravity on moves u, v, w and p
// bit for bit as the one with gravity off), the buoyancy kept in the closure
// (a stable theta lowers nu_t), the stage's refusal of residual inputs
// that do not match the option, and the ground ghosts of a frozen theta
// (the linear extension of the first two cells). simplec.theta =
// frozen_buoyant: the reader, theta held with the force kept in w, and the
// hydrostatic pressure of a held stratification at rest (closed form).
#include "BoundaryConditions.H"
#include "Closure.H"
#include "Poisson.H"
#include "PseudoTime.H"
#include "Residual.H"
#include "Simplec.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using gust::FaceType;
using gust::ptest::PCase;

// A flat periodic mesh with uniform layers, several boxes.
PCase flat_case (int nx, int ny, int nz, double Lx, double Ly, double H, int max_grid)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, ny, nz);
    c.in.mesh.prob_hi = {static_cast<amrex::Real>(Lx), static_cast<amrex::Real>(Ly), static_cast<amrex::Real>(H)};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {max_grid, max_grid};
    c.make([] (double, double) { return 0.0; });
    c.z.FillBoundary(c.geom.periodicity());
    c.mesh.build_from_ground(c.in, c.geom, c.z);
    gust::trtest::smooth_vertical_ghost_nodes(c);
    return c;
}

void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
}

// Random values on the valid points, then the periodic images (seeded per box
// index, so the field does not depend on the rank count).
void fill_random (amrex::MultiFab& mf, double lo, double hi, unsigned seed, const PCase& c)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        std::mt19937 gen(seed + 7919u * static_cast<unsigned>(mfi.index()));
        std::uniform_real_distribution<double> dist(lo, hi);
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(dist(gen)); });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
    if (!mf.ixType().cellCentered()) { mf.OverrideSync(c.geom.periodicity()); }
    mf.FillBoundary(c.geom.periodicity());
}

// theta = 300 + gamma z at the cell centres (a uniform stratification).
void fill_stable (amrex::MultiFab& theta, double gamma, const PCase& c)
{
    const amrex::Real dz = c.geom.CellSize(2);
    const auto g = static_cast<amrex::Real>(gamma);
    for (amrex::MFIter mfi(theta, false); mfi.isValid(); ++mfi) {
        const auto a = theta.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::ignore_unused(i, j);
            a(i, j, k) = amrex::Real(300) + g * (static_cast<amrex::Real>(k) + amrex::Real(0.5)) * dz;
        });
    }
    amrex::Gpu::streamSynchronize();
}

// d = p - a k (k+1) / 2 on the valid cells (a device lambda may not sit in a TEST body under nvcc).
void subtract_hydrostatic (amrex::MultiFab& d, const amrex::MultiFab& p, amrex::Real a)
{
    for (amrex::MFIter mfi(d, false); mfi.isValid(); ++mfi) {
        const auto o = d.array(mfi);
        const auto pp = p.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            o(i, j, k) = pp(i, j, k) - a * amrex::Real(0.5) * static_cast<amrex::Real>(k) * static_cast<amrex::Real>(k + 1);
        });
    }
    amrex::Gpu::streamSynchronize();
}

double max_abs_diff (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(d, a, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, b, 0, 0, 1, 0);
    return static_cast<double>(d.norm0());
}

void copy_state (gust::State& to, const gust::State& from)
{
    for (auto [t, f] : {std::pair{&to.u, &from.u}, std::pair{&to.v, &from.v}, std::pair{&to.w, &from.w},
                        std::pair{&to.p, &from.p}, std::pair{&to.theta, &from.theta}, std::pair{&to.tke, &from.tke}}) {
        amrex::MultiFab::Copy(*t, *f, 0, 0, 1, 2);
    }
}

// Fixed uniform diffusivities and slip walls (no closure, no wall model).
struct Laminar
{
    amrex::MultiFab nu, gt, gk;
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::ResidualInputs in;
    gust::PoissonOptions popt;
    void make (const PCase& c, double visc)
    {
        for (auto* mf : {&nu, &gt, &gk}) {
            mf->define(c.ba, c.dm, 1, 1);
            mf->setVal(static_cast<amrex::Real>(visc));
        }
        in.nu = &nu;
        in.gamma_theta = &gt;
        in.gamma_tke = &gk;
        in.theta_ref.z = {-1000, 5000};
        in.theta_ref.theta = {300, 300};
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall,
                    FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        in.bc = bc.get();
        popt.reltol = amrex::Real(1e-10);
    }
};

// A random state with a theta far from the reference (its buoyancy force is
// of order g/theta_0 * 5 K, 0.16 m/s^2, everywhere).
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

// n SIMPLEC iterations of the laminar stage from q (in place); with theta
// frozen, theta's residual norm is 0 at every iteration.
void run_laminar (gust::State& q, Laminar& lam, const gust::SimplecOptions& opt, const PCase& c, int n)
{
    gust::SimplecStage stage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm);
    stage.init(q);
    amrex::MultiFab dtau;
    for (int it = 0; it < n; ++it) {
        lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy);
        gust::local_dtau(q, &lam.nu, c.mesh, c.geom, amrex::Real(1), dtau);
        stage.iterate(q, dtau);
        const auto& names = stage.residual_names();
        for (std::size_t r = 0; r < names.size(); ++r) {
            if (opt.theta_frozen && names[r].find("theta") != std::string::npos) {
                EXPECT_EQ(stage.residual_norms()[r], 0) << "it " << it;
            }
        }
    }
}

} // namespace

// ---------------------------------------------------------------------------
// The reader: solve by default, frozen when given, any other value, a
// repeated key or two values refused with the key named.
TEST (FrozenTheta, ReaderAndRefusals)
{
    amrex::ParmParse pp("simplec");
    auto clear = [&] { while (pp.remove("theta") > 0) {} };
    clear();
    EXPECT_FALSE(gust::read_simplec_options().theta_frozen);
    pp.add("theta", std::string("solve"));
    EXPECT_FALSE(gust::read_simplec_options().theta_frozen);
    clear();
    pp.add("theta", std::string("frozen"));
    gust::SimplecOptions o = gust::read_simplec_options();
    EXPECT_TRUE(o.theta_frozen);
    ASSERT_EQ(o.given.size(), 1u);
    EXPECT_EQ(o.given[0], "theta");
    clear();
    for (const char* bad : {"Frozen", "fixed", "", "1"}) {
        pp.add("theta", std::string(bad));
        try {
            (void)gust::read_simplec_options();
            ADD_FAILURE() << "simplec.theta = '" << bad << "' was accepted";
        } catch (const std::runtime_error& e) {
            EXPECT_NE(std::string(e.what()).find("simplec.theta"), std::string::npos) << e.what();
        }
        clear();
    }
    pp.add("theta", std::string("frozen"));
    pp.add("theta", std::string("solve"));   // repeated
    EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
    clear();
    pp.addarr("theta", std::vector<std::string>{"frozen", "solve"});   // two values
    EXPECT_THROW(gust::read_simplec_options(), std::runtime_error);
    clear();
}

// ---------------------------------------------------------------------------
// Three iterations of the laminar stage from a random state whose theta is
// far from the reference, on 4 boxes: with theta frozen, theta is unchanged to
// the last bit and its residual is 0 at every iteration, and u, v, w and p are
// bit for bit those of the same frozen stage with gravity off (no buoyancy in
// the momentum residual nor in the predictor's w rows). The solved stage with
// gravity on moves theta and differs from it in w (the buoyancy acts there).
TEST (FrozenTheta, ThetaHeldAndNoBuoyancyInMomentum)
{
    PCase c = flat_case(8, 8, 8, 800, 800, 400, 4);
    gust::State q0;
    random_state(q0, c);

    gust::SimplecOptions frozen;
    frozen.theta_frozen = true;
    auto run = [&] (bool freeze, double gravity, gust::State& q) {
        Laminar lam;
        lam.make(c, 5.0);
        lam.in.buoyancy.gravity = static_cast<amrex::Real>(gravity);
        lam.in.buoyancy_in_momentum = !freeze;
        define_state(q, c);
        copy_state(q, q0);
        run_laminar(q, lam, freeze ? frozen : gust::SimplecOptions{}, c, 3);
    };
    gust::State a, b, s;
    run(true, 9.81, a);
    run(true, 0.0, b);
    run(false, 9.81, s);

    EXPECT_EQ(max_abs_diff(a.theta, q0.theta), 0);
    for (auto [f, g, name] : {std::tuple{&a.u, &b.u, "u"}, std::tuple{&a.v, &b.v, "v"}, std::tuple{&a.w, &b.w, "w"},
                              std::tuple{&a.p, &b.p, "p"}, std::tuple{&a.tke, &b.tke, "k"}}) {
        EXPECT_EQ(max_abs_diff(*f, *g), 0) << name;
    }
    // The old path: theta solved, buoyancy in w.
    EXPECT_GT(max_abs_diff(s.theta, q0.theta), 1e-6);
    EXPECT_GT(max_abs_diff(s.w, a.w), 1e-6);
}

// ---------------------------------------------------------------------------
// The closure keeps the buoyancy: with the k-L closure and theta frozen at a
// stable 10 K/km, one iteration gives a sum of nu_t below 0.9 of the one
// with gravity off (the stability functions and P_b see N^2), with 4- and 8-cell boxes.
TEST (FrozenTheta, ClosureKeepsTheBuoyancy)
{
    for (const int max_grid : {4, 8}) {
        PCase c = flat_case(8, 8, 8, 800, 800, 400, max_grid);
        gust::ClosureOptions copt;
        gust::Closure closure(copt);
        amrex::MultiFab dist;
        gust::height_above_ground(c.z, dist);
        gust::BoundaryInputs bin;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall, FaceType::periodic, FaceType::periodic,
                    FaceType::slip_wall};
        gust::PoissonOptions popt;
        gust::SimplecOptions opt;
        opt.theta_frozen = true;
        auto nu_sum = [&] (double gravity) {
            gust::ResidualInputs in;
            in.closure = &closure;
            in.wall_distance = &dist;
            in.theta_ref.z = {-1000, 5000};
            in.theta_ref.theta = {300, 300};
            in.buoyancy.gravity = static_cast<amrex::Real>(gravity);
            in.buoyancy_in_momentum = false;
            const gust::BoundaryConditions bc(bin, c.geom, c.mesh, in.coriolis);
            in.bc = &bc;
            gust::State q;
            define_state(q, c);
            q.u.setVal(6.0);
            q.v.setVal(0.0);
            q.w.setVal(0.0);
            q.p.setVal(0.0);
            fill_stable(q.theta, 0.01, c);
            q.tke.setVal(0.5);
            gust::SimplecStage stage(opt, 2, in, &closure, popt, c.mesh, c.geom, c.z, amrex::Real(1e-6), c.ba, c.dm);
            stage.init(q);
            amrex::MultiFab dtau;
            bc.fill_ghosts(q, nullptr, in.buoyancy);
            gust::local_dtau(q, &stage.nu(), c.mesh, c.geom, amrex::Real(1), dtau);
            stage.iterate(q, dtau);
            return static_cast<double>(stage.nu().sum(0));
        };
        const double with = nu_sum(9.81), without = nu_sum(0.0);
        EXPECT_LT(with, 0.9 * without) << with << " " << without << ", boxes of " << max_grid;
    }
}

// ---------------------------------------------------------------------------
// The stage refuses residual inputs that do not match the option: frozen with
// the momentum buoyancy on, and the momentum buoyancy off without frozen.
TEST (FrozenTheta, StageRefusesMismatchedResidualInputs)
{
    PCase c = flat_case(4, 4, 4, 400, 400, 200, 4);
    Laminar lam;
    lam.make(c, 1.0);
    gust::SimplecOptions opt;
    opt.theta_frozen = true;
    lam.in.buoyancy_in_momentum = true;
    EXPECT_THROW(gust::SimplecStage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm), std::runtime_error);
    opt.theta_frozen = false;
    lam.in.buoyancy_in_momentum = false;
    EXPECT_THROW(gust::SimplecStage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm), std::runtime_error);
    lam.in.buoyancy_in_momentum = true;
    EXPECT_NO_THROW(gust::SimplecStage(opt, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm));
}

// ---------------------------------------------------------------------------
// frozen_buoyant: theta held and the buoyancy kept in momentum; the switch
// that decides the residual inputs' buoyancy_in_momentum for each value.
TEST (FrozenTheta, BuoyantReaderAndSwitch)
{
    amrex::ParmParse pp("simplec");
    auto clear = [&] { while (pp.remove("theta") > 0) {} };
    clear();
    const gust::SimplecOptions solve = gust::read_simplec_options();
    EXPECT_TRUE(gust::buoyancy_in_momentum(solve));
    EXPECT_TRUE(gust::buoyancy_point_implicit(solve));
    pp.add("theta", std::string("frozen"));
    const gust::SimplecOptions frozen = gust::read_simplec_options();
    EXPECT_FALSE(frozen.theta_buoyant);
    EXPECT_FALSE(gust::buoyancy_in_momentum(frozen));
    EXPECT_FALSE(gust::buoyancy_point_implicit(frozen));
    clear();
    pp.add("theta", std::string("frozen_buoyant"));
    const gust::SimplecOptions fb = gust::read_simplec_options();
    EXPECT_TRUE(fb.theta_frozen);
    EXPECT_TRUE(fb.theta_buoyant);
    EXPECT_TRUE(gust::buoyancy_in_momentum(fb));
    EXPECT_FALSE(gust::buoyancy_point_implicit(fb));
    EXPECT_EQ(gust::theta_value(solve), "solve");
    EXPECT_EQ(gust::theta_value(frozen), "frozen");
    EXPECT_EQ(gust::theta_value(fb), "frozen_buoyant");
    clear();
    pp.add("theta", std::string("Frozen_Buoyant"));
    try {
        (void)gust::read_simplec_options();
        ADD_FAILURE() << "simplec.theta = Frozen_Buoyant was accepted";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("frozen_buoyant"), std::string::npos) << e.what();
    }
    clear();
    // The stage takes frozen_buoyant only with the momentum buoyancy on and its
    // point-implicit part off.
    PCase c = flat_case(4, 4, 4, 400, 400, 200, 4);
    Laminar lam;
    lam.make(c, 1.0);
    lam.in.buoyancy_in_momentum = false;
    lam.in.buoyancy_point_implicit = false;
    EXPECT_THROW(gust::SimplecStage(fb, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm), std::runtime_error);
    lam.in.buoyancy_in_momentum = true;
    lam.in.buoyancy_point_implicit = true;
    EXPECT_THROW(gust::SimplecStage(fb, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm), std::runtime_error);
    lam.in.buoyancy_point_implicit = false;
    EXPECT_NO_THROW(gust::SimplecStage(fb, 2, lam.in, nullptr, lam.popt, c.mesh, c.geom, c.z, 0, c.ba, c.dm));
}

// ---------------------------------------------------------------------------
// Three iterations from the random state of ThetaHeldAndNoBuoyancyInMomentum
// on 4 boxes: frozen_buoyant holds theta to the last bit, and its w differs
// from frozen's (the force acts); with gravity off the two are bit for bit
// the same (only the force differs).
TEST (FrozenTheta, BuoyantHoldsThetaAndKeepsTheForce)
{
    PCase c = flat_case(8, 8, 8, 800, 800, 400, 4);
    gust::State q0;
    random_state(q0, c);
    gust::SimplecOptions frozen, fb;
    frozen.theta_frozen = true;
    fb.theta_frozen = true;
    fb.theta_buoyant = true;
    auto run = [&] (const gust::SimplecOptions& opt, double gravity, gust::State& q) {
        Laminar lam;
        lam.make(c, 5.0);
        lam.in.buoyancy.gravity = static_cast<amrex::Real>(gravity);
        lam.in.buoyancy_in_momentum = gust::buoyancy_in_momentum(opt);
        lam.in.buoyancy_point_implicit = gust::buoyancy_point_implicit(opt);
        define_state(q, c);
        copy_state(q, q0);
        run_laminar(q, lam, opt, c, 3);
    };
    gust::State a, b, a0, b0;
    run(fb, 9.81, a);
    run(frozen, 9.81, b);
    run(fb, 0.0, a0);
    run(frozen, 0.0, b0);
    EXPECT_EQ(max_abs_diff(a.theta, q0.theta), 0);
    EXPECT_GT(max_abs_diff(a.w, b.w), 1e-3);
    for (auto [f, g, name] : {std::tuple{&a0.u, &b0.u, "u"}, std::tuple{&a0.v, &b0.v, "v"}, std::tuple{&a0.w, &b0.w, "w"},
                              std::tuple{&a0.p, &b0.p, "p"}}) {
        EXPECT_EQ(max_abs_diff(*f, *g), 0) << name;
    }
}

// ---------------------------------------------------------------------------
// The predictor's point-implicit buoyancy follows buoyancy_point_implicit:
// with theta solved (where it is on) three iterations from the random state
// give another w with it off. frozen_buoyant turns it off (the switch above).
TEST (FrozenTheta, PointImplicitSwitchActs)
{
    PCase c = flat_case(8, 8, 8, 800, 800, 400, 4);
    gust::State q0;
    random_state(q0, c);
    gust::SimplecOptions solve;
    auto run = [&] (bool point_implicit, gust::State& q) {
        Laminar lam;
        lam.make(c, 5.0);
        lam.in.buoyancy.gravity = static_cast<amrex::Real>(9.81);
        lam.in.buoyancy_point_implicit = point_implicit;
        define_state(q, c);
        copy_state(q, q0);
        run_laminar(q, lam, solve, c, 3);
    };
    gust::State on, off;
    run(true, on);
    run(false, off);
    EXPECT_GT(max_abs_diff(on.w, off.w), 1e-6);
}

// ---------------------------------------------------------------------------
// A held stratification theta = 300 + gamma z at rest, the reference
// theta_ref = 300 + gamma_ref z and theta_0 = 290 K, on 4 boxes: the steady
// state is rest with the hydrostatic pressure dp/dz = g (theta - theta_ref) /
// theta_0. On uniform layers the w face k+1/2 sits at z = (k+1) dz and its
// theta - theta_ref is the mean of its two cells', so p(k+1) - p(k) =
// (g (gamma - gamma_ref) / theta_0) (k+1) dz^2 and p(k) = p(0) + a k (k+1) / 2
// with a = g (gamma - gamma_ref) dz^2 / theta_0. A swap of theta_0 and the
// reference, or the reference at another height, fails. frozen_buoyant reaches it (u, v, w at rest to
// round-off); frozen keeps p at 0 (no force).
TEST (FrozenTheta, BuoyantHydrostaticPressure)
{
    const double gamma = 0.01;
    const double gamma_ref = 0.004;
    PCase c = flat_case(8, 8, 16, 800, 800, 800, 4);
    const amrex::Real dz = c.geom.CellSize(2);
    for (const bool buoyant : {true, false}) {
        gust::SimplecOptions opt;
        opt.theta_frozen = true;
        opt.theta_buoyant = buoyant;
        Laminar lam;
        lam.make(c, 5.0);
        lam.in.theta_ref.theta = {static_cast<amrex::Real>(300 - 1000 * gamma_ref),
                                  static_cast<amrex::Real>(300 + 5000 * gamma_ref)};
        lam.in.buoyancy.theta0 = 290;
        lam.in.buoyancy_in_momentum = gust::buoyancy_in_momentum(opt);
        lam.in.buoyancy_point_implicit = false;
        gust::State q;
        define_state(q, c);
        q.u.setVal(0.0);
        q.v.setVal(0.0);
        q.w.setVal(0.0);
        q.p.setVal(0.0);
        fill_stable(q.theta, gamma, c);
        q.tke.setVal(0.5);
        run_laminar(q, lam, opt, c, 200);
        const amrex::Real a = buoyant ? static_cast<amrex::Real>(lam.in.buoyancy.gravity * (gamma - gamma_ref) / 290.0) * dz * dz
                                      : amrex::Real(0);
        amrex::MultiFab d(c.ba, c.dm, 1, 0);
        subtract_hydrostatic(d, q.p, a);
        const double spread = static_cast<double>(d.max(0) - d.min(0));
        const double prange = static_cast<double>(a) * 0.5 * 15 * 16;   // p(15) - p(0)
        EXPECT_LT(spread, buoyant ? 1e-6 * prange : 1e-12) << "buoyant " << buoyant << ", p range " << prange;
        if (buoyant) { EXPECT_GT(prange, 1.0); }
        // At rest to round-off (the force is of order 1e-3 m/s^2, p of order 100 m^2/s^2).
        const double rest = 1e4 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
        for (auto [f, name] : {std::pair{&q.u, "u"}, std::pair{&q.v, "v"}, std::pair{&q.w, "w"}}) {
            EXPECT_LT(static_cast<double>(f->norm0()), rest) << name << ", buoyant " << buoyant;
        }
    }
}

// ---------------------------------------------------------------------------
// With the ground extension (simplec.theta = frozen), theta's two ground ghost
// layers continue the first two cells linearly: on uniform layers
// theta(-1) = 2 theta(0) - theta(1) and theta(-2) = 3 theta(0) - 2 theta(1),
// on every column of 4 boxes (lateral ghost columns included), so the first
// cell's centred d theta/dz is (theta(1) - theta(0))/dz. Off (the default),
// and in a homogeneous fill with it on, the slip ground's mirror holds:
// theta(-1) = theta(0), theta(-2) = theta(1).
TEST (FrozenTheta, GroundGhostsExtendTheFirstTwoCells)
{
    PCase c = flat_case(8, 8, 6, 800, 800, 300, 4);
    gust::State q0;
    random_state(q0, c);
    for (const int mode : {0, 1, 2}) {   // 0: off, 1: on, 2: on but homogeneous
        Laminar lam;
        lam.make(c, 1.0);
        lam.bc->set_theta_ground_extension(mode != 0);
        gust::State q;
        define_state(q, c);
        copy_state(q, q0);
        lam.bc->fill_ghosts(q, nullptr, lam.in.buoyancy, mode == 2);
        double worst = 0;
        for (amrex::MFIter mfi(q.theta, false); mfi.isValid(); ++mfi) {
            amrex::FArrayBox f(q.theta[mfi].box(), 1, amrex::The_Pinned_Arena());
            f.copy<amrex::RunOn::Device>(q.theta[mfi], q.theta[mfi].box());
            amrex::Gpu::streamSynchronize();
            amrex::Box b = mfi.growntilebox(amrex::IntVect(2, 2, 0));
            b.setRange(2, 0, 1);
            amrex::LoopOnCpu(b, [&] (int i, int j, int) {
                const double t0 = f(amrex::IntVect(i, j, 0)), t1 = f(amrex::IntVect(i, j, 1));
                const double g1 = f(amrex::IntVect(i, j, -1)), g2 = f(amrex::IntVect(i, j, -2));
                const double e1 = (mode == 1) ? 2 * t0 - t1 : t0, e2 = (mode == 1) ? 3 * t0 - 2 * t1 : t1;
                worst = std::max({worst, std::abs(g1 - e1), std::abs(g2 - e2)});
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        EXPECT_LT(worst, 1e-10) << "mode " << mode;
    }
}
