// The lee turn of an immersed ground (terrain.lee_turn, ImmersedLeeTurn.H):
//   the cell's rate against the closed form: off on the windward side, the ramp, the layer count n;
//   on a plane falling along x and along y, the sink on the w rows and their a_P against the closed form
//     at every z face, with a reference wind that changes with height (it pins the reference cell), a
//     non-uniform w and a non-uniform J (it pins the cut cell's height), on one and several boxes;
//   the same plane with the wind reversed (windward everywhere) gives nothing; off when terrain.lee_turn = 0.
#include "BuildingForcing.H"
#include "ImmersedLeeTurn.H"
#include "ImmersedTerrain.H"
#include "ImplicitOperator.H"
#include "LaminarTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <string>
#include <type_traits>

namespace {

using gust::ptest::PCase;

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double rtol = kDouble ? 1e-11 : 2e-5;
constexpr double kZ0 = 0.1;

// The closed form of the cell's rate, in double: lambda and lambda w_t.
void closed_form (double sx, double sy, double ur, double vr, double dz, double dx, double dy, double c0, double eps, double& lam,
                  double& lam_wt)
{
    lam = 0;
    lam_wt = 0;
    const double wt = sx * ur + sy * vr;
    if (wt >= 0) { return; }
    const double U = std::max(std::hypot(ur, vr), 1e-3);
    const double f = std::min(1.0, -wt / (eps * U));
    const double n = std::max(1.0, std::hypot(sx * dx, sy * dy) / dz);
    lam = f * c0 * n * n * U / dz;
    lam_wt = lam * wt;
}

// The node heights with the ground h(x, y) as their second component.
amrex::MultiFab with_ground (const PCase& c, const std::function<double(double, double)>& h)
{
    amrex::MultiFab z(c.z.boxArray(), c.dm, 2, c.z.nGrow());
    amrex::MultiFab::Copy(z, c.z, 0, 0, 1, c.z.nGrow());
    const auto dx = c.geom.CellSizeArray();
    for (amrex::MFIter mfi(z, false); mfi.isValid(); ++mfi) {
        const auto a = z.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
            a(i, j, k, 1) = static_cast<amrex::Real>(h(i * static_cast<double>(dx[0]), j * static_cast<double>(dx[1])));
        });
    }
    return z;
}

// J = 1 + 0.05 i + 0.03 j at the cells (one ghost layer): the cut cell's height J dzeta varies.
double jac (int i, int j) { return 1.0 + 0.05 * i + 0.03 * j; }

// The w residual and a_P of the immersed ground's law on a plane, with the lee turn at c0 (0: off).
struct Rows
{
    amrex::MultiFab Rw, A;
};

struct Plane
{
    gust::Buildings b;
    amrex::MultiFab detJ, u, v, w;
    amrex::MultiFab ox, oy, oz;
    const PCase& c;
    Plane (const PCase& cc, const std::function<double(double, double)>& h, int dir, double along, double cross) : c(cc)
    {
        gust::TerrainInputs t;
        t.representation = "immersed_forcing";
        t.flat_height = 25;
        const gust::BuildingsOptions opt = gust::immersed_terrain_options(t, kZ0);
        const amrex::MultiFab z = with_ground(c, h);
        b.build_from_surface(opt, "test", c.geom, z, c.ba, c.dm);
        detJ.define(c.ba, c.dm, 1, 1);
        for (amrex::MFIter mfi(detJ, false); mfi.isValid(); ++mfi) {
            const auto a = detJ.array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { a(i, j, k) = static_cast<amrex::Real>(jac(i, j)); });
        }
        // The wind along the slope grows with height and varies along x and y, so that a cell's two faces
        // differ (A (1 + 0.05 k + 0.03 i + 0.02 j)); across it, uniform.
        u.define(amrex::convert(c.ba, amrex::IntVect(1, 0, 0)), c.dm, 1, 1);
        v.define(amrex::convert(c.ba, amrex::IntVect(0, 1, 0)), c.dm, 1, 1);
        w.define(amrex::convert(c.ba, amrex::IntVect(0, 0, 1)), c.dm, 1, 1);
        for (int n = 0; n < 3; ++n) {
            amrex::MultiFab& q = (n == 0) ? u : (n == 1) ? v : w;
            for (amrex::MFIter mfi(q, false); mfi.isValid(); ++mfi) {
                const auto a = q.array(mfi);
                amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
                    double val = 0;
                    if (n == 2) {
                        val = 0.2 * std::sin(0.7 * i + 0.3 * j + 0.5 * k);
                    } else {
                        val = (n == dir) ? along * (1 + 0.05 * k + 0.03 * i + 0.02 * j) : cross;
                    }
                    a(i, j, k) = static_cast<amrex::Real>(val);
                });
            }
        }
        amrex::MultiFab* o[3] = {&ox, &oy, &oz};
        for (int d = 0; d < 3; ++d) {
            o[d]->define(b.beta_face(d).boxArray(), c.dm, 1, 1);
            o[d]->setVal(1);
            amrex::MultiFab::LinComb(*o[d], 1, *o[d], 0, -1, b.beta_face(d), 0, 0, 1, 0);
            o[d]->FillBoundary(c.geom.periodicity());
        }
    }
    Rows rows (double c0, double eps)
    {
        gust::TerrainInputs t;
        t.representation = "immersed_forcing";
        t.flat_height = 25;
        gust::BuildingsOptions opt = gust::immersed_terrain_options(t, kZ0);
        opt.ground_lee_turn = c0;
        opt.ground_lee_turn_ramp = eps;
        gust::BuildingForcing law;
        law.build(b, opt, c.geom, detJ);
        EXPECT_EQ(law.lee_turn().active(), c0 > 0);
        law.update(u, v, w);
        Rows r;
        amrex::MultiFab Ru(u.boxArray(), c.dm, 1, 0), Rv(v.boxArray(), c.dm, 1, 0);
        r.Rw.define(w.boxArray(), c.dm, 1, 0);
        Ru.setVal(0);
        Rv.setVal(0);
        r.Rw.setVal(0);
        law.add_sink(u, v, w, Ru, Rv, r.Rw);
        r.A.define(w.boxArray(), c.dm, gust::n_implicit_coefs, 0);
        r.A.setVal(1);
        law.add_to_implicit(2, {&ox, &oy, &oz}, r.A);
        return r;
    }
};

struct Errors
{
    double rw{0}, ap{0}, scale{0}, faces{0};
};

// The lee turn's terms (on minus off) against the closed form at every z face.
Errors compare (Plane& p, int dir, double m, double c0, double eps)
{
    const Rows on = p.rows(c0, eps);
    const Rows off = p.rows(0, eps);
    const auto dxa = p.c.geom.CellSizeArray();
    const double dx = dxa[0], dy = dxa[1], dzeta = dxa[2];
    const int khi = p.c.geom.Domain().bigEnd(2);
    const double sx = dir == 0 ? -m : 0.0, sy = dir == 1 ? -m : 0.0;
    Errors e;
    for (amrex::MFIter mfi(on.Rw, false); mfi.isValid(); ++mfi) {
        const auto ron = on.Rw.const_array(mfi), roff = off.Rw.const_array(mfi);
        const auto aon = on.A.const_array(mfi), aoff = off.A.const_array(mfi);
        const auto beta = p.b.beta().const_array(mfi);
        const auto fz = p.b.beta_face(2).const_array(mfi);
        const auto uu = p.u.const_array(mfi), vv = p.v.const_array(mfi), ww = p.w.const_array(mfi);
        // The cell's lambda and lambda w_t; 0 outside a cut cell and above the top.
        auto cell = [&] (int i, int j, int k, double& lam, double& lw) {
            lam = 0;
            lw = 0;
            if (k > khi) { return; }
            const double be = beta(i, j, k);
            if (be <= 0 || be >= 1) { return; }
            int kr = k + 1;
            while (kr < khi && beta(i, j, kr) > 0) { ++kr; }
            const double ur = 0.5 * (uu(i, j, kr) + uu(i + 1, j, kr));
            const double vr = 0.5 * (vv(i, j, kr) + vv(i, j + 1, kr));
            closed_form(sx, sy, ur, vr, jac(i, j) * dzeta, dx, dy, c0, eps, lam, lw);
        };
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            double want_r = 0, want_a = 0;
            if (k > p.c.geom.Domain().smallEnd(2) && fz(i, j, k) < 1) {
                double l0 = 0, w0 = 0, l1 = 0, w1 = 0;
                cell(i, j, k - 1, l0, w0);
                cell(i, j, k, l1, w1);
                want_a = 0.5 * (l0 + l1);
                want_r = want_a * ww(i, j, k) - 0.5 * (w0 + w1);
            }
            const double got_r = static_cast<double>(ron(i, j, k)) - roff(i, j, k);
            const double got_a = static_cast<double>(aon(i, j, k, gust::coef_P)) - aoff(i, j, k, gust::coef_P);
            e.rw = std::max(e.rw, std::abs(got_r - want_r));
            e.ap = std::max(e.ap, std::abs(got_a - want_a));
            e.scale = std::max(e.scale, std::abs(want_a));
            e.faces += (want_a > 0) ? 1 : 0;
            // The other couplings do not change.
            for (int n = 1; n < gust::n_implicit_coefs; ++n) { e.ap = std::max(e.ap, std::abs(static_cast<double>(aon(i, j, k, n)) - aoff(i, j, k, n))); }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(e.rw);
    amrex::ParallelDescriptor::ReduceRealMax(e.ap);
    amrex::ParallelDescriptor::ReduceRealMax(e.scale);
    amrex::ParallelDescriptor::ReduceRealSum(e.faces);
    return e;
}

}   // namespace

// The cell's rate: windward (w_t >= 0) nothing; lee with the ramp's weight and the layer count n.
TEST (ImmersedLeeTurn, RateMatchesTheClosedForm)
{
    struct Case { double sx, sy, ur, vr, dz, c0, eps; };
    const Case cases[] = {
        {-1.5, 0.0, 8.0, 3.0, 10.0, 0.25, 0.05},    // n = 1.5, f = 1
        {-0.02, 0.0, 8.0, 3.0, 10.0, 0.25, 0.05},   // n = 1, f = 0.02 * 8 / (0.05 |U|) < 1
        {0.3, -0.4, 2.0, 5.0, 4.0, 1.0, 0.2},        // both slopes: w_t = -1.4, n = 5 / 4
        {0.0, -0.1, -1.0, 0.5, 2.0, 0.25, 0.05},     // the wind's other component
        {-1.5, 0.0, 0.0, 0.0, 10.0, 0.25, 0.05},     // at rest: w_t = 0, nothing
        {1.5, 0.0, 8.0, 3.0, 10.0, 0.25, 0.05},      // windward: nothing
        {-0.5, 0.0, 1e-5, 0.0, 10.0, 0.25, 0.05},    // below the smallest speed: |U| = 1e-3 in the rate
    };
    for (const auto& k : cases) {
        double lam = 0, lw = 0;
        closed_form(k.sx, k.sy, k.ur, k.vr, k.dz, 10.0, 10.0, k.c0, k.eps, lam, lw);
        amrex::Real gl = -1, gw = -1;
        gust::lee_turn_rate(static_cast<amrex::Real>(k.sx), static_cast<amrex::Real>(k.sy), static_cast<amrex::Real>(k.ur),
                            static_cast<amrex::Real>(k.vr), static_cast<amrex::Real>(k.dz), amrex::Real(10), amrex::Real(10),
                            static_cast<amrex::Real>(k.c0), static_cast<amrex::Real>(k.eps), gust::kBuildingMinSpeed, gl, gw);
        EXPECT_NEAR(gl, lam, rtol * lam) << k.sx << " " << k.sy << " " << k.ur << " " << k.vr;
        EXPECT_NEAR(gw, lw, rtol * std::abs(lw)) << k.sx << " " << k.sy << " " << k.ur << " " << k.vr;
    }
    // Hand values: n = 1.5, |U| = sqrt(73), f = 1: lambda = 0.25 * 2.25 * sqrt(73) / 10.
    amrex::Real gl = 0, gw = 0;
    gust::lee_turn_rate(amrex::Real(-1.5), 0, 8, 3, 10, 10, 10, amrex::Real(0.25), amrex::Real(0.05), gust::kBuildingMinSpeed, gl, gw);
    EXPECT_NEAR(gl, 0.25 * 2.25 * std::sqrt(73.0) / 10, rtol);
    EXPECT_NEAR(gw, -12.0 * 0.25 * 2.25 * std::sqrt(73.0) / 10, 10 * rtol);
    // The ramp: f = 0.02 * 8 / (0.05 sqrt(73)), n = 1.
    gust::lee_turn_rate(amrex::Real(-0.02), 0, 8, 3, 10, 10, 10, amrex::Real(0.25), amrex::Real(0.05), gust::kBuildingMinSpeed, gl, gw);
    EXPECT_NEAR(gl, 0.16 / (0.05 * std::sqrt(73.0)) * 0.25 * std::sqrt(73.0) / 10, rtol);
}

// The sink on the w rows and their a_P on planes falling along x and along y, steep (n = 1.5 / J, f = 1)
// and gentle (n = 1, the ramp), on one box and on boxes of 4 x 4 columns.
TEST (ImmersedLeeTurn, PlaneMatchesTheClosedForm)
{
    for (const int mgs : {8, 4}) {
        // 8 x 8 x 24 cells of 10 m x 10 m x 10 m.
        const PCase c = gust::lamtest::flat_case(8, 8, 24, 80, 80, 240, false, false, mgs);
        for (const double m : {1.5, 0.02}) {
            for (const int dir : {0, 1}) {
                SCOPED_TRACE("boxes of " + std::to_string(mgs) + ", slope " + std::to_string(m) + " along " + (dir == 0 ? "x" : "y"));
                const double h0 = m > 1 ? 150.0 : 25.0;
                Plane p(c, [=] (double x, double y) { return h0 - m * (dir == 0 ? x : y); }, dir, 8.0, 3.0);
                const Errors e = compare(p, dir, m, 0.25, 0.05);
                EXPECT_GE(e.faces, 64.0) << "a lee face in every column at least";
                EXPECT_GT(e.scale, 0.0);
                EXPECT_LT(e.rw, rtol * std::max(1.0, 10 * e.scale));
                EXPECT_LT(e.ap, rtol * std::max(1.0, e.scale));
            }
        }
    }
}

// The wind reversed: every cut cell is windward and the terms are 0; and terrain.lee_turn = 0 is off.
TEST (ImmersedLeeTurn, WindwardAndOffGiveNothing)
{
    const PCase c = gust::lamtest::flat_case(8, 8, 24, 80, 80, 240, false, false, 4);
    Plane p(c, [] (double x, double) { return 150.0 - 1.5 * x; }, 0, -8.0, 3.0);
    const Errors e = compare(p, 0, 1.5, 0.25, 0.05);
    EXPECT_EQ(e.faces, 0.0);
    EXPECT_EQ(e.rw, 0.0);
    EXPECT_EQ(e.ap, 0.0);
    const Rows off = p.rows(0, 0.05);   // checks lee_turn().active() is false
    (void)off;
}
