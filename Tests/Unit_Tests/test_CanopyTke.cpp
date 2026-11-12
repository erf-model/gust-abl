// The forest canopy's terms of the k–L closure (CanopyTke.H, the length cap
// in Closure_K.H), plan rev. 6 beyond the plan, sub-phase 8.3 (UT-8.3.x):
//   the wake production and short-circuit sink have their closed forms and
//   balance at beta_p |U|^2 / beta_d; the sink is on the Patankar and exact
//   diagonals (UT-8.3.1);
//   the length cap acts only where f > 0 and alpha < L f, and eps, nu_t and
//   Rt follow the capped length (UT-8.3.2);
//   the exact derivative of S_k with the canopy terms matches central
//   differences with the switches frozen, in stably stratified cells with and
//   without canopy (UT-8.3.3; it found the missing d(cmu'/cmu)/dk of P_b);
//   a canopy whose drag is zero everywhere, or no canopy, gives every closure
//   field bit for bit, and the residual of k gains D_c - P_c exactly
//   (UT-8.3.4);
//   the cap is a switch: frozen, it holds where the state has moved past it,
//   and changed_switches counts those cells (UT-8.3.5);
//   the ground-roughness warning counts the canopy columns.
#include "BoundaryConditions.H"
#include "Canopy.H"
#include "CanopyTke.H"
#include "Closure.H"
#include "PoissonTestHelpers.H"
#include "Reductions.H"
#include "Residual.H"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>

namespace fs = std::filesystem;

namespace {

using gust::FaceType;
using gust::ptest::PCase;

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;

// The patch: centre (200, 200), diameter 300, height 17.3 m, Cd 0.25, LAI 4,
// uniform profile: f = Cd LAI / h in the cells wholly below the top.
constexpr double kCd = 0.25, kLai = 4.0, kH = 17.3, kR = 150.0;
constexpr double kF = kCd * kLai / kH;

// 16 x 16 x 24 cells over 400 m x 400 m, 1 m first layer stretched by 1.1,
// boxes of mgs cells, periodic in x and y, flat.
PCase make_case (int mgs)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 16, 24);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.mesh.dz0 = amrex::Real(1);
    c.in.mesh.stretch_ratio = amrex::Real(1.1);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

gust::Canopy make_canopy (const PCase& c)
{
    const fs::path dir = fs::temp_directory_path() / ("gust_canopy_tke_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / "patch.txt").string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        os << "1 200 200 " << kH << " " << 2 * kR << " " << kCd << " " << kLai << " 0.8\n";
    }
    amrex::ParallelDescriptor::Barrier();
    gust::CanopyOptions o;
    o.type = gust::CanopyType::patches;
    o.type_name = "patches";
    o.file = file;
    o.given = {"type", "file"};
    gust::Canopy cn;
    cn.build(o, c.geom, c.z, c.ba, c.dm);
    return cn;
}

// Fill a field's points, ghosts included, from g at its index-space position.
void fill_component (amrex::MultiFab& mf, int c, const std::function<double(double, double, double)>& g)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const double X = i + (c == 0 ? 0.0 : 0.5), Y = j + (c == 1 ? 0.0 : 0.5), Z = k + (c == 2 ? 0.0 : 0.5);
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(g(X, Y, Z));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

constexpr double kw = 2 * M_PI / 16;

// A periodic, sheared, stratified state (stable: N^2 > 0) with k varying
// around k0.
void fill_state (gust::State& q, const PCase& c, double k0, double dtheta)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
    fill_component(q.u, 0, [] (double X, double Y, double Z) { return 2 + 0.3 * Z + 0.5 * std::sin(kw * Y); });
    fill_component(q.v, 1, [] (double X, double, double Z) { return 0.4 * std::cos(kw * X) + 0.05 * Z; });
    fill_component(q.w, 2, [] (double X, double Y, double) { return 0.05 * std::sin(kw * (X + Y)); });
    q.p.setVal(0.0);
    fill_component(q.theta, 3, [dtheta] (double, double, double Z) { return 300 + dtheta * Z; });
    fill_component(q.tke, 3, [k0] (double X, double Y, double) { return k0 * (1 + 0.3 * std::sin(kw * X) * std::cos(kw * Y)); });
}

amrex::FArrayBox host (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    const amrex::Box b = mf[mfi].box();
    amrex::FArrayBox h(b, mf.nComp(), amrex::The_Pinned_Arena());
    h.copy<amrex::RunOn::Device>(mf[mfi], b, 0, b, 0, mf.nComp());
    amrex::Gpu::streamSynchronize();
    return h;
}

struct TkeSetup
{
    PCase c;
    gust::Canopy canopy;
    gust::State q;
    amrex::MultiFab dist;
    explicit TkeSetup (double k0 = 0.5, double dtheta = 0.0) : c(make_case(8))
    {
        canopy = make_canopy(c);
        fill_state(q, c, k0, dtheta);
        gust::height_above_ground(q.z_phys_nd, dist);
    }
};

const gust::Roughness kZ0{amrex::Real(0.05)};

} // namespace

// UT-8.3.1: with a uniform wind of 5 m/s, P_c = beta_p f |U|^3 and
// D_c = beta_d f |U| k on every cell (f the drag field), and where
// k = beta_p |U|^2 / beta_d they balance; the Patankar diagonal gains
// beta_d f |U| and the exact one loses it (no cap: alpha = 0, so the rest of
// the closure is the run without a canopy).
TEST(CanopyTke, SourcesHaveTheirClosedForms)
{
    TkeSetup s(0.5);
    s.q.u.setVal(3.0);
    s.q.v.setVal(4.0);
    s.q.w.setVal(0.0);
    const double U = 5, bp = 0.8, bd = 3.5, kbal = bp * U * U / bd;
    // Half the cells at the balance k, the others at 0.5.
    fill_component(s.q.tke, 3, [kbal] (double X, double, double) { return (static_cast<int>(X) % 2 == 0) ? kbal : 0.5; });
    gust::Closure plain(gust::ClosureOptions{}), canopy(gust::ClosureOptions{});
    canopy.set_canopy(&s.canopy.drag(), static_cast<amrex::Real>(bp), static_cast<amrex::Real>(bd), 0);
    gust::ClosureFields f0, f1;
    plain.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f0);
    canopy.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f1);
    ASSERT_TRUE(f1.canopy);
    EXPECT_FALSE(f0.canopy);
    const double t = kDouble ? 1e-13 : 1e-5;
    double worst = 0, balance = 0;
    amrex::Long in_canopy = 0, balanced = 0;
    for (amrex::MFIter mfi(f1.Pc, false); mfi.isValid(); ++mfi) {
        const auto pc = host(f1.Pc, mfi), dc = host(f1.Dc, mfi), d1 = host(f1.dSdk, mfi), d0 = host(f0.dSdk, mfi),
                   e1 = host(f1.dSdk_exact, mfi), e0 = host(f0.dSdk_exact, mfi), fd = host(s.canopy.drag(), mfi),
                   kk = host(s.q.tke, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double f = fd(iv), kv = kk(iv);
            const double scale = bp * kF * U * U * U;
            worst = std::max({worst, std::abs(pc(iv) - bp * f * U * U * U) / scale, std::abs(dc(iv) - bd * f * U * kv) / scale,
                              std::abs((d1(iv) - d0(iv)) - bd * f * U) / (bd * kF * U),
                              std::abs((e0(iv) - e1(iv)) - bd * f * U) / (bd * kF * U)});
            if (f > 0) { ++in_canopy; }
            if (f > 0 && kv == static_cast<amrex::Real>(kbal)) {
                balance = std::max(balance, static_cast<double>(std::abs(pc(iv) - dc(iv)) / pc(iv)));
                ++balanced;
            }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    amrex::ParallelDescriptor::ReduceRealMax(balance);
    amrex::ParallelDescriptor::ReduceLongSum(in_canopy);
    amrex::ParallelDescriptor::ReduceLongSum(balanced);
    EXPECT_LE(worst, t);
    EXPECT_LE(balance, t);
    EXPECT_GT(in_canopy, 100L);
    EXPECT_GT(balanced, 40L);
}

// UT-8.3.1, with the cap: where the canopy terms and the capped dissipation
// alone act (a uniform wind: no shear production; neutral), P_c = D_c + eps at
// k = kappa_c |U|^2 with beta_d kappa_c + Cmu0^3 kappa_c^1.5 / alpha = beta_p,
// kappa_c = 0.1728 for ERF's defaults (docs/canopy.rst), in every capped cell.
TEST(CanopyTke, CappedBalanceIsKappaU2)
{
    const double bp = 1, bd = 4, alpha = 0.04, C3 = 0.5562 * 0.5562 * 0.5562, U = 5;
    double lo = 0, hi = 1;
    for (int it = 0; it < 200; ++it) {
        const double m = 0.5 * (lo + hi);
        (bd * m + C3 * std::pow(m, 1.5) / alpha < bp ? lo : hi) = m;
    }
    const double kappa = 0.5 * (lo + hi);
    EXPECT_NEAR(kappa, 0.1728, 1e-4);
    TkeSetup s(kappa * U * U);
    s.q.u.setVal(3.0);
    s.q.v.setVal(4.0);
    s.q.w.setVal(0.0);
    s.q.tke.setVal(static_cast<amrex::Real>(kappa * U * U));
    gust::Closure cl(gust::ClosureOptions{});
    cl.set_canopy(&s.canopy.drag(), static_cast<amrex::Real>(bp), static_cast<amrex::Real>(bd), static_cast<amrex::Real>(alpha));
    gust::ClosureFields f;
    cl.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f);
    double worst = 0;
    amrex::Long capped = 0;
    for (amrex::MFIter mfi(f.Pc, false); mfi.isValid(); ++mfi) {
        const auto pc = host(f.Pc, mfi), dc = host(f.Dc, mfi), ep = host(f.eps, mfi), ps = host(f.Ps, mfi),
                   L = host(f.length, mfi), fd = host(s.canopy.drag(), mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (!(fd(iv) > 0) || L(iv) != static_cast<amrex::Real>(alpha) / fd(iv)) { return; }
            worst = std::max({worst, static_cast<double>(std::abs(pc(iv) - dc(iv) - ep(iv)) / pc(iv)), std::abs(static_cast<double>(ps(iv)))});
            ++capped;
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    amrex::ParallelDescriptor::ReduceLongSum(capped);
    EXPECT_LE(worst, kDouble ? 1e-12 : 1e-5);
    EXPECT_GT(capped, 100L);
}

// UT-8.3.2: L = min(L_free, alpha / f) where f > 0 (L_free the run without
// a canopy), the same bits where f = 0 or with alpha = 0; eps and nu_t follow
// the capped L. Neutral: L_free = l_g, so the first cells (l_g < alpha / f)
// keep it and the cells above are capped.
TEST(CanopyTke, LengthCapActsOnlyInTheCanopy)
{
    TkeSetup s(0.5);
    const double alpha = 0.04;
    gust::Closure plain(gust::ClosureOptions{}), capped(gust::ClosureOptions{}), uncapped(gust::ClosureOptions{});
    capped.set_canopy(&s.canopy.drag(), 1, 4, static_cast<amrex::Real>(alpha));
    uncapped.set_canopy(&s.canopy.drag(), 1, 4, 0);
    gust::ClosureFields f0, f1, f2;
    plain.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f0);
    capped.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f1);
    uncapped.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f2);
    const gust::ClosureConstants cst = plain.constants();
    const double t = kDouble ? 1e-13 : 1e-5;
    double worst = 0;
    amrex::Long changed_outside = 0, uncapped_diff = 0, n_capped = 0, n_kept = 0;
    for (amrex::MFIter mfi(f0.length, false); mfi.isValid(); ++mfi) {
        const auto L0 = host(f0.length, mfi), L1 = host(f1.length, mfi), L2 = host(f2.length, mfi),
                   e1 = host(f1.eps, mfi), n1 = host(f1.nut, mfi), n0 = host(f0.nut, mfi), fd = host(s.canopy.drag(), mfi),
                   kk = host(s.q.tke, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double f = fd(iv);
            if (L2(iv) != L0(iv)) { ++uncapped_diff; }
            if (!(f > 0)) {
                if (L1(iv) != L0(iv) || n1(iv) != n0(iv)) { ++changed_outside; }
                return;
            }
            const double want = std::min(static_cast<double>(L0(iv)), alpha / f);
            worst = std::max(worst, std::abs(L1(iv) - want) / want);
            const double kv = kk(iv);
            worst = std::max(worst, std::abs(e1(iv) - cst.Cmu0_pow3 * std::pow(kv, 1.5) / want) / e1(iv));
            if (want < L0(iv)) { ++n_capped; } else { ++n_kept; }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    for (amrex::Long* n : {&changed_outside, &uncapped_diff, &n_capped, &n_kept}) { amrex::ParallelDescriptor::ReduceLongSum(*n); }
    EXPECT_LE(worst, t);
    EXPECT_EQ(changed_outside, 0L);
    EXPECT_EQ(uncapped_diff, 0L);
    EXPECT_GT(n_capped, 50L);
    EXPECT_GT(n_kept, 20L);
}

// UT-8.3.3: with the switches frozen, dS_k/dk of the closure with the canopy
// (P_s + P_b - eps + P_c - D_c) by central differences equals dSdk_exact;
// stably stratified, so the stable length depends on k, and capped cells
// (whose L does not) are among them.
TEST(CanopyTke, ExactDerivativeMatchesDifferences)
{
    TkeSetup s(0.2, 0.5);
    gust::Closure cl(gust::ClosureOptions{});
    cl.set_canopy(&s.canopy.drag(), 1, 4, amrex::Real(0.04));
    gust::ClosureFields f;
    cl.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f);
    cl.freeze_switches(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom);
    const double h = 1e-6;
    const double alpha_r = static_cast<double>(amrex::Real(0.04));
    const double cap_tol = kDouble ? 1e-12 : 4 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());
    auto source_at = [&] (double rel, amrex::MultiFab& S) {
        gust::State qp;
        qp.define(s.c.ba, s.c.dm);
        for (const auto& p : {std::make_pair(&qp.u, &s.q.u), std::make_pair(&qp.v, &s.q.v), std::make_pair(&qp.w, &s.q.w),
                        std::make_pair(&qp.theta, &s.q.theta), std::make_pair(&qp.tke, &s.q.tke),
                        std::make_pair(&qp.p, &s.q.p), std::make_pair(&qp.z_phys_nd, &s.q.z_phys_nd)}) {
            amrex::MultiFab::Copy(*p.first, *p.second, 0, 0, 1, p.second->nGrow());
        }
        qp.tke.mult(static_cast<amrex::Real>(1 + rel), 0, 1, qp.tke.nGrow());
        gust::ClosureFields g;
        cl.evaluate(qp, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, g);
        S.define(s.c.ba, s.c.dm, 1, 0);
        amrex::MultiFab::Copy(S, g.Ps, 0, 0, 1, 0);
        amrex::MultiFab::Add(S, g.Pb, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(S, g.eps, 0, 0, 1, 0);
        amrex::MultiFab::Add(S, g.Pc, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(S, g.Dc, 0, 0, 1, 0);
    };
    amrex::MultiFab Sp, Sm;
    source_at(h, Sp);
    source_at(-h, Sm);
    double worst = 0;
    amrex::Long n_capped_cells = 0;
    for (amrex::MFIter mfi(Sp, false); mfi.isValid(); ++mfi) {
        const auto sp = host(Sp, mfi), sm = host(Sm, mfi), de = host(f.dSdk_exact, mfi), kk = host(s.q.tke, mfi),
                   fd = host(s.canopy.drag(), mfi), L = host(f.length, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double dk = 2 * h * kk(iv);
            const double fdiff = (sp(iv) - sm(iv)) / dk;
            worst = std::max(worst, std::abs(fdiff - de(iv)) / (std::abs(de(iv)) + 1e-3));
            // capped: L is the cap amrex::Real(alpha) / f, to its last rounding (1e-12 in double)
            if (fd(iv) > 0 && std::abs(L(iv) - alpha_r / fd(iv)) < cap_tol * L(iv)) { ++n_capped_cells; }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    amrex::ParallelDescriptor::ReduceLongSum(n_capped_cells);
    cl.thaw();
    if (kDouble) { EXPECT_LE(worst, 1e-6); }
    EXPECT_GT(n_capped_cells, 20L);
}

// UT-8.3.4: a canopy whose drag field is zero everywhere leaves every closure
// field bit for bit, with P_c = D_c = 0; and with a canopy and no cap the k
// residual gains exactly D_c - P_c.
TEST(CanopyTke, ZeroCanopyIsBitIdenticalAndResidualCarriesTheTerms)
{
    TkeSetup s(0.5, 0.2);
    amrex::MultiFab zero(s.c.ba, s.c.dm, 1, amrex::IntVect(1, 1, 0));
    zero.setVal(0.0);
    gust::Closure plain(gust::ClosureOptions{}), empty(gust::ClosureOptions{});
    empty.set_canopy(&zero, 1, 4, amrex::Real(0.04));
    gust::ClosureFields f0, f1;
    plain.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f0);
    empty.evaluate(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, f1);
    const amrex::MultiFab* a[] = {&f0.nu, &f0.gamma_theta, &f0.gamma_tke, &f0.nut, &f0.length, &f0.Rt, &f0.Ps, &f0.Pb,
                                  &f0.eps, &f0.dSdk, &f0.dSdk_exact, &f0.dnutdk};
    const amrex::MultiFab* b[] = {&f1.nu, &f1.gamma_theta, &f1.gamma_tke, &f1.nut, &f1.length, &f1.Rt, &f1.Ps, &f1.Pb,
                                  &f1.eps, &f1.dSdk, &f1.dSdk_exact, &f1.dnutdk};
    for (int n = 0; n < 12; ++n) {
        EXPECT_EQ(gust::count_bit_differences(*a[n], *b[n], amrex::IntVect(0)), 0) << "field " << n;
    }
    EXPECT_EQ(static_cast<double>(f1.Pc.norm0()), 0.0);
    EXPECT_EQ(static_cast<double>(f1.Dc.norm0()), 0.0);

    // The residual: diffusivities from the closure (no boundary conditions, periodic).
    gust::Closure with(gust::ClosureOptions{});
    with.set_canopy(&s.canopy.drag(), 1, 4, 0);
    gust::ResidualInputs in0;
    in0.closure = &plain;
    in0.wall_distance = &s.dist;
    in0.theta_ref.z = {0};
    in0.theta_ref.theta = {300};
    gust::ResidualInputs in1 = in0;
    in1.closure = &with;
    gust::ResidualFields R0, R1;
    gust::ClosureFields c1;
    gust::residual(s.q, in0, s.c.mesh, s.c.geom, R0);
    gust::residual(s.q, in1, s.c.mesh, s.c.geom, R1, &c1);
    amrex::MultiFab want(s.c.ba, s.c.dm, 1, 0);
    amrex::MultiFab::Copy(want, c1.Dc, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(want, c1.Pc, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(R1.tke, R0.tke, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(R1.tke, want, 0, 0, 1, 0);
    EXPECT_GT(static_cast<double>(want.norm0()), 0.0);
    EXPECT_LE(static_cast<double>(R1.tke.norm0()), (kDouble ? 1e-12 : 1e-4) * static_cast<double>(want.norm0()) + 1e-300);
}

// UT-8.3.5: the cap is a switch. At the point level, a capped decision held
// for a smaller k keeps L = alpha / f while the live decision does not cap;
// on the field, the cells whose live decision moved are counted, and the
// frozen evaluation keeps their capped length.
TEST(CanopyTke, CapSwitchIsFrozen)
{
    const gust::ClosureConstants c{};
    const amrex::Real f = static_cast<amrex::Real>(kF), alpha = amrex::Real(0.04), l_g = 5, N2 = amrex::Real(0.01);
    const amrex::Real L_cap = alpha / f;
    const gust::Switches s_hi = gust::closure::switches_for(amrex::Real(0.1), N2, l_g, c, f, alpha);
    const gust::Switches s_lo = gust::closure::switches_for(amrex::Real(0.005), N2, l_g, c, f, alpha);
    EXPECT_TRUE(s_hi.canopy_capped);
    EXPECT_FALSE(s_lo.canopy_capped);
    const auto held = gust::closure::closure_point<amrex::Real>(amrex::Real(0.005), N2, 0, l_g, c, s_hi, L_cap);
    const auto live_pt = gust::closure::closure_point<amrex::Real>(amrex::Real(0.005), N2, 0, l_g, c, s_lo, L_cap);
    EXPECT_EQ(held.L, L_cap);
    EXPECT_LT(live_pt.L, L_cap);
    const gust::Switches none = gust::closure::switches_for(amrex::Real(0.1), N2, l_g, c, f, 0);
    EXPECT_FALSE(none.canopy_capped);

    TkeSetup s(0.2, 0.5);
    gust::Closure cl(gust::ClosureOptions{});
    cl.set_canopy(&s.canopy.drag(), 1, 4, alpha);
    cl.freeze_switches(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom);
    EXPECT_EQ(cl.changed_switches(s.q, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom), 0);
    gust::State q2;
    q2.define(s.c.ba, s.c.dm);
    for (const auto& p : {std::make_pair(&q2.u, &s.q.u), std::make_pair(&q2.v, &s.q.v), std::make_pair(&q2.w, &s.q.w),
                    std::make_pair(&q2.theta, &s.q.theta), std::make_pair(&q2.tke, &s.q.tke), std::make_pair(&q2.p, &s.q.p),
                    std::make_pair(&q2.z_phys_nd, &s.q.z_phys_nd)}) {
        amrex::MultiFab::Copy(*p.first, *p.second, 0, 0, 1, p.second->nGrow());
    }
    q2.tke.mult(amrex::Real(0.02), 0, 1, q2.tke.nGrow());
    const amrex::Long moved = cl.changed_switches(q2, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom);
    EXPECT_GT(moved, 0);
    gust::ClosureFields frozen, live;
    cl.evaluate(q2, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, frozen);
    cl.thaw();
    cl.evaluate(q2, s.dist, kZ0, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, live);
    amrex::Long held_cells = 0;
    for (amrex::MFIter mfi(frozen.length, false); mfi.isValid(); ++mfi) {
        const auto Lf = host(frozen.length, mfi), Ll = host(live.length, mfi), fd = host(s.canopy.drag(), mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (fd(iv) > 0 && Lf(iv) == static_cast<amrex::Real>(0.04) / fd(iv) && Ll(iv) < Lf(iv)) { ++held_cells; }
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(held_cells);
    EXPECT_GT(held_cells, 0L);
}

// The ground-roughness warning: every canopy column is counted when the
// ground z0 is forest-sized, none when it is the litter's.
TEST(CanopyTke, RoughGroundUnderTheCanopyIsCounted)
{
    PCase c = make_case(4);
    const gust::Canopy cn = make_canopy(c);
    amrex::Long want = 0;
    for (int j = 0; j < 16; ++j) {
        for (int i = 0; i < 16; ++i) {
            const double x = (i + 0.5) * 25 - 200, y = (j + 0.5) * 25 - 200;
            if (x * x + y * y <= kR * kR) { ++want; }
        }
    }
    EXPECT_EQ(gust::canopy_columns_on_rough_ground(cn.drag(), gust::Roughness{amrex::Real(1.0)}, gust::kForestSizedZ0), want);
    EXPECT_EQ(gust::canopy_columns_on_rough_ground(cn.drag(), gust::Roughness{amrex::Real(0.05)}, gust::kForestSizedZ0), 0);
}

// The canopy cells the inflow-face and damping-layer warnings count
// (BoundaryConditions::canopy_in_boundary_zones), on a domain without periodic
// sides: four narrow patches each cover only the first cell column (or row)
// inside one face, one patch sits inside. With one face inflow at a time, its
// patch's cells are counted, and only those; an outflow face counts nothing.
// With a 100 m sponge on the x faces and a Rayleigh layer above 15 m, the
// canopy cells whose centre lies in either are counted; the Rayleigh layer on
// w alone (it relaxes w to zero, not towards the profile) is not counted.
TEST(CanopyTke, BoundaryZonesCountTheCanopyCells)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 16, 24);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {0, 0, 0};
    c.in.mesh.max_grid_size = {4, 4};
    c.in.mesh.dz0 = amrex::Real(1);
    c.in.mesh.stretch_ratio = amrex::Real(1.1);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    const fs::path dir = fs::temp_directory_path() / ("gust_canopy_zones_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / "patches.txt").string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        for (const char* xy : {"12.5 200", "387.5 100", "200 12.5", "300 387.5"}) {   // one column or row each
            os << "1 " << xy << " " << kH << " 30 " << kCd << " " << kLai << " 0.8\n";
        }
        os << "1 260 200 " << kH << " 120 " << kCd << " " << kLai << " 0.8\n";   // inside
    }
    amrex::ParallelDescriptor::Barrier();
    gust::CanopyOptions o;
    o.type = gust::CanopyType::patches;
    o.type_name = "patches";
    o.file = file;
    o.given = {"type", "file"};
    gust::Canopy cn;
    cn.build(o, c.geom, c.z, c.ba, c.dm);

    // The expected counts from the drag field on the host: face 0..3 = xlo, xhi, ylo, yhi.
    const auto& zl = c.mesh.zlevels();
    auto expect = [&] (int face, double sponge, double zdamp) {
        amrex::Long n_in = 0, n_damp = 0;
        for (amrex::MFIter mfi(cn.drag(), false); mfi.isValid(); ++mfi) {
            const auto fd = host(cn.drag(), mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (!(fd(amrex::IntVect(i, j, k)) > 0)) { return; }
                const double x = (i + 0.5) * 25, zc = 0.5 * (zl[static_cast<std::size_t>(k)] + zl[static_cast<std::size_t>(k + 1)]);
                const bool on = (face == 0 && i == 0) || (face == 1 && i == 15) || (face == 2 && j == 0) || (face == 3 && j == 15);
                if (on) { ++n_in; }
                const bool sp = sponge > 0 && (x < sponge || 400 - x < sponge);
                const bool ry = zdamp > 0 && zc > zdamp;
                if (sp || ry) { ++n_damp; }
            });
        }
        amrex::ParallelDescriptor::ReduceLongSum(n_in);
        amrex::ParallelDescriptor::ReduceLongSum(n_damp);
        return std::make_pair(n_in, n_damp);
    };
    gust::CoriolisOptions cor;
    auto counts = [&] (int inflow_face, double sponge, double zdamp, bool w_only = false) {
        gust::BoundaryInputs bin;
        bin.face = {FaceType::outflow, FaceType::slip_wall, FaceType::surface_layer, FaceType::outflow, FaceType::slip_wall,
                    FaceType::slip_wall};
        const int slot[4] = {0, 3, 1, 4};   // xlo, xhi, ylo, yhi in AMReX's orientation order
        if (inflow_face >= 0) { bin.face[static_cast<std::size_t>(slot[inflow_face])] = FaceType::inflow; }
        bin.sponge_width = static_cast<amrex::Real>(sponge);
        bin.sponge_strength = sponge > 0 ? amrex::Real(0.01) : amrex::Real(0);
        bin.rayleigh_zdamp = static_cast<amrex::Real>(zdamp);
        bin.rayleigh_dampcoef = zdamp > 0 ? amrex::Real(0.2) : amrex::Real(0);
        if (w_only) { bin.rayleigh_fields = {false, false, true, false, false}; }
        gust::BoundaryConditions bc(bin, c.geom, c.mesh, cor);
        return bc.canopy_in_boundary_zones(cn.drag(), c.z);
    };
    const auto none = counts(-1, 0, 0);
    EXPECT_EQ(none.inflow, 0);
    EXPECT_EQ(none.damping, 0);
    for (int face = 0; face < 4; ++face) {
        const auto in = counts(face, 0, 0);
        EXPECT_EQ(in.inflow, expect(face, 0, 0).first) << "face " << face;
        EXPECT_GT(in.inflow, 10) << "face " << face;
        EXPECT_EQ(in.damping, 0) << "face " << face;
    }
    const auto damp = counts(0, 100, 15);
    const auto want = expect(0, 100, 15);
    EXPECT_EQ(damp.inflow, want.first);
    EXPECT_EQ(damp.damping, want.second);
    EXPECT_GT(damp.damping, damp.inflow);
    // The Rayleigh layer on w alone: only the sponge counts.
    EXPECT_EQ(counts(0, 100, 15, true).damping, expect(0, 100, 0).second);
}
