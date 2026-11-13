// The physical time term of dual time stepping (TimeTerm.H), plan rev. 6
// beyond the plan, sub-phase 6.1 (UT-6.1.x), on the 2-D hill of
// HillTestHelpers.H as the application sets it up (the k–L closure, the wall
// law with the Dirichlet wall-k rows, the slip top, Rayleigh damping of w and
// theta, the ground and top w faces as Dirichlet rows):
//   the residual gains exactly (q - q^n)/dt on every row of u, v, w, theta and
//   k that keeps its own residual, against the closed form on a non-uniform
//   field, nothing on p, and nothing on the rows that replace theirs (the
//   ground and top w faces, the first cells' wall-k rows) (UT-6.1.1);
//   implicit_rows gains 1/dt on the diagonal of exactly those rows, no other
//   coefficient changes, and the change of the residual along a direction v
//   is that diagonal times v, by finite differences (UT-6.1.2);
//   at q = q^n the residual is unchanged and without a time term nothing
//   reads it (UT-6.1.3; the steady path's bits are the regression golds);
//   the option checks name unsteady.dt.
#include "Closure.H"
#include "HillTestHelpers.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"
#include "Residual.H"
#include "SurfaceLayer.H"
#include "TestHelpers.H"
#include "TimeTerm.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>

namespace {

using gust::Equation;
using gust::Loc;
using gust::hilltest::HillRun;

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double kDt = 2.0;   // s
constexpr double kPi = 3.14159265358979323846;

// The hill at nx x 2 x nz on four boxes in x (max_grid 4), stably stratified.
HillRun make_hill () { return HillRun(0.003, 16, 12, 4); }

amrex::MultiFab& field (gust::State& q, int f)
{
    amrex::MultiFab* p[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return *p[f];
}
const amrex::MultiFab& rfield (const gust::ResidualFields& R, int f)
{
    const amrex::MultiFab* p[6] = {&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke};
    return *p[f];
}

// q = q0 + s dq with dq a smooth non-uniform field, nonzero everywhere; ghosts refilled.
void perturb (HillRun& h, const gust::State& q0, gust::State& q, double s)
{
    q.define(h.c.ba, h.c.dm);
    auto t = q.fields();
    // State::fields() is non-const; copy through a mutable alias of q0's fields.
    const amrex::MultiFab* f0[7] = {&q0.u, &q0.v, &q0.w, &q0.p, &q0.theta, &q0.tke, &q0.z_phys_nd};
    for (std::size_t n = 0; n < t.size(); ++n) {
        amrex::MultiFab::Copy(*t[n].mf, *f0[n], 0, 0, 1, t[n].mf->nGrowVect());
    }
    const double L = HillRun::Lx, H = HillRun::H;
    const Loc locs[6] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell, Loc::cell};
    const double amp[6] = {0.3, 0.2, 0.05, 0.4, 0.5, 0.02};
    for (int f = 0; f < 6; ++f) {
        amrex::MultiFab d(field(q, f).boxArray(), field(q, f).DistributionMap(), 1, field(q, f).nGrow());
        const double a = amp[f] * s;
        gust::trtest::fill_loc(d, locs[f], h.c, [=] (double x, double, double z) {
            return a * (1.2 + std::sin(2 * kPi * x / L + f) * std::cos(kPi * z / H));
        });
        amrex::MultiFab::Add(field(q, f), d, 0, 0, 1, field(q, f).nGrow());
    }
    h.ls->bc.fill_ghosts(q, h.ls->surface(), h.s.fopt.buoyancy);
}

struct Compare
{
    double max_err{0};      // max |dR - expect| / scale over the rows that keep their residual
    long replaced{0};       // rows where dR == 0 exactly (the perturbation is nonzero except where the
                            // ghost fill sets a boundary row, w = 0 on the slip top)
    long kept{0};
};

// Over the valid points of field f: dR = Ron - Roff against expect = c (q - qn)/dt.
Compare compare (const amrex::MultiFab& Ron, const amrex::MultiFab& Roff, const amrex::MultiFab& q,
                 const amrex::MultiFab& qn, double c)
{
    Compare out;
    for (amrex::MFIter mfi(Ron, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox a(b, 1, amrex::The_Pinned_Arena()), o(b, 1, amrex::The_Pinned_Arena()),
            x(b, 1, amrex::The_Pinned_Arena()), xn(b, 1, amrex::The_Pinned_Arena());
        a.copy<amrex::RunOn::Device>(Ron[mfi], b);
        o.copy<amrex::RunOn::Device>(Roff[mfi], b);
        x.copy<amrex::RunOn::Device>(q[mfi], b);
        xn.copy<amrex::RunOn::Device>(qn[mfi], b);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double dR = static_cast<double>(a(iv)) - static_cast<double>(o(iv));
            const double expect = c * (static_cast<double>(x(iv)) - static_cast<double>(xn(iv))) / kDt;
            // A replaced row has dR == 0 exactly. The smallest time term here is about 2e-3, far above the
            // round-off of any R, so a kept row cannot reach 0 by cancellation; keep that if the amplitudes change.
            if (dR == 0) { ++out.replaced; return; }
            ++out.kept;
            const double scale = std::abs(static_cast<double>(a(iv))) + std::abs(static_cast<double>(o(iv))) + std::abs(expect);
            out.max_err = std::max(out.max_err, std::abs(dR - expect) / scale);
        });
    }
    amrex::ParallelAllReduce::Max(out.max_err, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(out.replaced, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(out.kept, amrex::ParallelDescriptor::Communicator());
    return out;
}

void eval (HillRun& h, const gust::State& q, const gust::ResidualInputs& in, gust::ResidualFields& R)
{
    gust::residual(q, in, h.c.mesh, h.c.geom, R);
}

// The implicit rows of eq at q, with in's time term or none.
void rows (HillRun& h, const gust::State& q, const gust::ResidualInputs& in, Equation eq, amrex::MultiFab& A)
{
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;
    gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, h.c.mesh, h.c.geom, F);
    in.bc->outflow_transport(F);
    gust::SurfaceState ss;
    gust::compute_surface_state(q, *in.surface, in.buoyancy, h.c.mesh, h.c.geom, ss);
    gust::ClosureFields cf;
    in.closure->evaluate(q, *in.wall_distance, gust::ground_roughness(in.surface), in.buoyancy, h.c.mesh, h.c.geom, cf);
    gust::implicit_rows({&F[0], &F[1], &F[2]}, q.z_phys_nd, in, eq, h.c.mesh, h.c.geom, A,
                        eq == Equation::w ? nullptr : &ss, &cf);
}

// The (row count of each field) of the ground and top w faces and of the first cells of k.
long columns (const HillRun& h)
{
    const amrex::Box& d = h.c.geom.Domain();
    return static_cast<long>(d.length(0)) * d.length(1);
}

} // namespace

TEST(TimeTerm, ResidualGainsTheClosedFormOnItsOwnRows)   // UT-6.1.1
{
    HillRun h = make_hill();
    ASSERT_EQ(h.ls->res_in.surface->wall_k, gust::SurfaceLayerOptions::WallK::dirichlet);
    gust::TimeTerm tt(gust::TimeTermOptions{amrex::Real(kDt)}, h.c.ba, h.c.dm);
    tt.store(h.q);
    gust::State q;
    perturb(h, h.q, q, 1.0);
    gust::ResidualInputs off = h.ls->res_in;
    gust::ResidualInputs on = off;
    on.time = &tt;
    gust::ResidualFields Ron, Roff;
    eval(h, q, on, Ron);
    eval(h, q, off, Roff);
    const double tol = kDouble ? 1e-12 : 1e-4;
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    const int fi[5] = {0, 1, 2, 4, 5};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    const long expect_replaced[5] = {0, 0, 2 * columns(h), 0, columns(h)};
    for (int n = 0; n < 5; ++n) {
        const Compare c = compare(rfield(Ron, fi[n]), rfield(Roff, fi[n]), field(q, fi[n]), tt.previous(eqs[n]), 1.0);
        EXPECT_LT(c.max_err, tol) << names[n];
        EXPECT_EQ(c.replaced, expect_replaced[n]) << names[n] << ": rows without the time term";
        EXPECT_GT(c.kept, 0) << names[n];
    }
    // Continuity has no time term.
    EXPECT_EQ(gust::test::count_bit_differences(Ron.p, Roff.p, 0), 0);
}

TEST(TimeTerm, DiagonalIsTheDerivativeOfTheTimeTerm)   // UT-6.1.2
{
    HillRun h = make_hill();
    gust::TimeTerm tt(gust::TimeTermOptions{amrex::Real(kDt)}, h.c.ba, h.c.dm);
    tt.store(h.q);
    gust::State q, q2;
    perturb(h, h.q, q, 1.0);
    // The term is linear, so the step size only sets the round-off: under float the change of R must stay
    // well above ulp(R), hence the larger step.
    const double eps = kDouble ? 1e-3 : 1e-1;
    perturb(h, q, q2, eps);
    gust::ResidualInputs off = h.ls->res_in;
    gust::ResidualInputs on = off;
    on.time = &tt;
    // dR(x) = R_on(x) - R_off(x), the time term's contribution at x.
    gust::ResidualFields a1, b1, a2, b2;
    eval(h, q, on, a1);
    eval(h, q, off, b1);
    eval(h, q2, on, a2);
    eval(h, q2, off, b2);
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    const int fi[5] = {0, 1, 2, 4, 5};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    const double tol_fd = kDouble ? 1e-6 : 1e-2;
    for (int n = 0; n < 5; ++n) {
        amrex::MultiFab Aon, Aoff;
        rows(h, q, on, eqs[n], Aon);
        rows(h, q, off, eqs[n], Aoff);
        // Every coefficient but a_P bit for bit; a_P + 1/dt or a_P (an identity row).
        long others = 0, plus = 0, same = 0, wrong = 0;
        double fd_err = 0;
        const double inv = 1 / kDt;
        for (amrex::MFIter mfi(Aon, false); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox();
            auto host = [&] (const amrex::MultiFab& mf, int nc) {
                amrex::FArrayBox f(b, nc, amrex::The_Pinned_Arena());
                f.copy<amrex::RunOn::Device>(mf[mfi], b, 0, b, 0, nc);
                return f;
            };
            const auto ha = host(Aon, gust::n_implicit_coefs), hb = host(Aoff, gust::n_implicit_coefs);
            const auto r1a = host(rfield(a1, fi[n]), 1), r1b = host(rfield(b1, fi[n]), 1);
            const auto r2a = host(rfield(a2, fi[n]), 1), r2b = host(rfield(b2, fi[n]), 1);
            const auto x1 = host(field(q, fi[n]), 1), x2 = host(field(q2, fi[n]), 1);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                for (int c = 1; c < gust::n_implicit_coefs; ++c) {
                    if (ha(iv, c) != hb(iv, c)) { ++others; }
                }
                const double dP = static_cast<double>(ha(iv, gust::coef_P)) - static_cast<double>(hb(iv, gust::coef_P));
                if (dP == 0) { ++same; }
                else if (std::abs(dP - inv) <= (kDouble ? 1e-12 : 1e-4) * (std::abs(static_cast<double>(ha(iv, gust::coef_P))) + inv)) { ++plus; }
                else { ++wrong; }
                // The change of dR along q2 - q, against dP times that change.
                const double ddR = (static_cast<double>(r2a(iv)) - static_cast<double>(r2b(iv)))
                                   - (static_cast<double>(r1a(iv)) - static_cast<double>(r1b(iv)));
                const double dx = static_cast<double>(x2(iv)) - static_cast<double>(x1(iv));
                // Only rows that keep the term: an identity row (dP == 0) can have dx == 0 (the w faces the
                // ghost fill pins), where the ratio below would be 0/0 (FE_INVALID under the trap).
                if (dP == 0) { return; }
                if (dx == 0) { ++wrong; return; }
                fd_err = std::max(fd_err, std::abs(ddR - dP * dx) / (std::abs(dx) * inv));
            });
        }
        auto comm = amrex::ParallelDescriptor::Communicator();
        amrex::ParallelAllReduce::Sum(others, comm);
        amrex::ParallelAllReduce::Sum(plus, comm);
        amrex::ParallelAllReduce::Sum(same, comm);
        amrex::ParallelAllReduce::Sum(wrong, comm);
        amrex::ParallelAllReduce::Max(fd_err, comm);
        const long expect_same = (n == 2) ? 2 * columns(h) : (n == 4) ? columns(h) : 0;
        EXPECT_EQ(others, 0) << names[n] << ": off-diagonal coefficients changed";
        EXPECT_EQ(wrong, 0) << names[n] << ": a_P changed by neither 0 nor 1/dt";
        EXPECT_EQ(same, expect_same) << names[n] << ": identity rows";
        EXPECT_GT(plus, 0) << names[n];
        EXPECT_LT(fd_err, tol_fd) << names[n] << ": the diagonal is not the derivative of the time term";
    }
}

TEST(TimeTerm, NoChangeAtTheStoredState)   // UT-6.1.3
{
    HillRun h = make_hill();
    gust::TimeTerm tt(gust::TimeTermOptions{amrex::Real(kDt)}, h.c.ba, h.c.dm);
    tt.store(h.q);
    gust::ResidualInputs off = h.ls->res_in;
    gust::ResidualInputs on = off;
    on.time = &tt;
    gust::ResidualFields Ron, Roff;
    eval(h, h.q, on, Ron);
    eval(h, h.q, off, Roff);
    // (q - q^n)/dt = 0 adds nothing (by value: R + 0 turns -0 into +0). This shows the term vanishes at the
    // stored state; that the steady path (time == nullptr) is bit-identical is shown by the regression golds.
    const auto a = Ron.fields();
    const auto b = Roff.fields();
    for (std::size_t n = 0; n < a.size(); ++n) {
        amrex::MultiFab d(a[n].mf->boxArray(), a[n].mf->DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(d, *a[n].mf, 0, 0, 1, 0);
        amrex::MultiFab::Subtract(d, *b[n].mf, 0, 0, 1, 0);
        EXPECT_EQ(d.norm0(), 0) << a[n].name;
    }
}

TEST(TimeTerm, OptionsNameTheInput)
{
    EXPECT_TRUE(gust::check_time_term_options(gust::TimeTermOptions{amrex::Real(0.5)}).empty());
    const amrex::Real bad[] = {0, amrex::Real(-1), std::numeric_limits<amrex::Real>::infinity(),
                               std::numeric_limits<amrex::Real>::quiet_NaN(), std::numeric_limits<amrex::Real>::denorm_min()};
    for (const amrex::Real dt : bad) {
        const auto err = gust::check_time_term_options(gust::TimeTermOptions{dt});
        ASSERT_EQ(err.size(), 1u) << dt;
        EXPECT_NE(err.front().find("unsteady.dt"), std::string::npos);
    }
}
