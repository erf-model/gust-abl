// The patches' feedback by injection (NestInject.H; plan rev. 6 tracker, row
// N-114 item 2, nest.feedback = inject, beyond the plan):
//   UT-I.1: on a hill parent with uneven boxes, the nest restricted into a
//     copy q* (restrict_patch): inject sets u, v, w, p, theta and k to q*'s
//     values bit for bit on the inner box's points (cells; faces from L to
//     H + 1 along their normal), with the point counts written out here, and
//     leaves every other point bit for bit; an inactive injection changes
//     nothing. On the inner cells J div equals that of q*, the mean of the
//     nest's J div over each cell's 2 x 2 nest cells (flux conservation); of
//     the cells outside, only the ring cells that share a face with the inner
//     box change their divergence.
//   UT-I.2: on a flat parent with uneven boxes and dx != dy (two footprints,
//     one spanning a periodic direction): zero_rows sets R to 0 on the inner
//     points of u, v, w, p, theta and k and leaves every other point bit for
//     bit; differences() is the max and RMS of |q - q*| over the inner
//     points (each point once), computed here from a gathered field; the same
//     values bit for bit on reversed ranks and to round-off on one box (and
//     on 1, 2 and 4 ranks, the unit_np* tests). Inactive: nothing.
//   UT-I.3: on the hill of HillTestHelpers.H as the application sets it up
//     (wall law, wall-k rows, Dirichlet ground and top w faces): an inactive
//     injection leaves the residual and the implicit rows bit for bit; an
//     active one zeroes the inner rows of u, v, w, p, theta and k, including
//     the rows the wall law and the boundaries replace, keeps the other rows
//     bit for bit, and leaves the implicit rows bit for bit.
//   UT-I.4 (the stop test): the hill's SIMPLEC driver with a target that is
//     not the parent's solution (the converged state with its velocity times
//     1.1 on the inner box), injected before every iteration: the residual
//     falls by six orders from its first value (to 2e-4 of it in single
//     precision, whose round-off floor is 1e-4), and the state outside the
//     inner box moves away from the converged one (the injection acts on it).
//   UT-I.7: the target's change that the joint run's stop window reads:
//     omega (q* - q) on the first inject, omega (1 - omega) (q* - q) on the
//     second, and the sizes, from gathered fields on uneven boxes.
#include "HillTestHelpers.H"
#include "ImplicitOperator.H"
#include "Inputs.H"
#include "MassFlux.H"
#include "Mesh.H"
#include "NestDefect.H"
#include "NestInject.H"
#include "Poisson.H"
#include "PoissonTestHelpers.H"
#include "Residual.H"
#include "State.H"
#include "SurfaceLayer.H"
#include "TestHelpers.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using gust::Equation;
using gust::Loc;

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double kTol = kDouble ? 1.0e-11 : 2.0e-4;

struct Case
{
    gust::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    gust::Mesh mesh;
    gust::State q;

    void make (bool reverse = false)
    {
        geom = gust::make_geometry(in.mesh);
        ba = gust::make_box_array(in.mesh);
        dm = amrex::DistributionMapping(ba);
        if (reverse) {
            amrex::Vector<int> pmap = dm.ProcessorMap();
            const int np = amrex::ParallelDescriptor::NProcs();
            for (auto& p : pmap) { p = np - 1 - p; }
            dm = amrex::DistributionMapping(pmap);
        }
        q.define(ba, dm);
        mesh.define(ba, dm);
        mesh.build(in, geom, q.z_phys_nd);
    }
};

// The parent: 20 x 16 x 8 cells of 100 m x 100 m (hill) or 100 m x 80 m (flat, periodic_y optional), stretched
// layers, uneven boxes of at most mgs x (mgs - 1) cells.
gust::Inputs parent_inputs (bool hill, bool periodic_y = false, int mgs = 7)
{
    gust::Inputs in;
    in.mesh.n_cell = amrex::IntVect(20, 16, 8);
    in.mesh.prob_lo = {-1000, 500, 0};
    in.mesh.prob_hi = {1000, hill ? amrex::Real(2100) : amrex::Real(1780), 800};
    in.mesh.is_periodic = {0, periodic_y ? 1 : 0, 0};
    in.mesh.max_grid_size = {mgs, mgs - 1};
    in.mesh.dz0 = 40;
    if (hill) {
        in.terrain.type = "agnesi";
        in.terrain.hill_shape = "round";
        in.terrain.hill_height = 150;
        in.terrain.hill_half_width = 300;
        in.terrain.smoothing = 0;
    }
    return in;
}

// The hill's nest over parent cells 3 to 14 in x and 2 to 11 in y: 24 x 20 cells of 50 m, boxes of at most 8 x 6.
gust::Inputs nest_inputs (const gust::Inputs& parent)
{
    gust::Inputs in = parent;
    in.mesh.n_cell = amrex::IntVect(24, 20, 8);
    in.mesh.prob_lo = {-700, 700, 0};
    in.mesh.prob_hi = {500, 1700, 800};
    in.mesh.max_grid_size = {8, 6};
    in.terrain.hill_centre = {0, 1300};
    return in;
}

const amrex::IntVect kOffset(3, 2, 0);
constexpr int kR = 2;

void fill_all (amrex::MultiFab& mf, const std::function<double(int, int, int)>& fn)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(i, j, k)); });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

void each_valid (const amrex::MultiFab& mf, const std::function<void(int, int, int, double)>& fn)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox h(b, 1, amrex::The_Pinned_Arena());
        h.copy<amrex::RunOn::Device>(mf[mfi], b);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { fn(i, j, k, static_cast<double>(h(amrex::IntVect(i, j, k)))); });
    }
}

// The whole of mf (valid points) on every rank, indexed by its domain box; a face shared by two boxes averaged.
std::vector<double> gather (const amrex::MultiFab& mf, const amrex::Box& dom)
{
    const amrex::Box b = amrex::convert(dom, mf.ixType());
    std::vector<double> out(static_cast<std::size_t>(b.numPts()), 0.0), cnt(out.size(), 0.0);
    each_valid(mf, [&] (int i, int j, int k, double v) {
        const auto n = static_cast<std::size_t>(b.index(amrex::IntVect(i, j, k)));
        out[n] += v;
        cnt[n] += 1.0;
    });
    amrex::ParallelDescriptor::ReduceRealSum(out.data(), static_cast<int>(out.size()));
    amrex::ParallelDescriptor::ReduceRealSum(cnt.data(), static_cast<int>(cnt.size()));
    for (std::size_t n = 0; n < out.size(); ++n) { if (cnt[n] > 0.0) { out[n] /= cnt[n]; } }
    return out;
}

double pos (const amrex::Geometry& g, amrex::IndexType t, int i, int d)
{
    const double off = t.nodeCentered(d) ? 0.0 : 0.5;
    return static_cast<double>(g.ProbLo(d)) + (i + off) * static_cast<double>(g.CellSize(d));
}

// A smooth, non-solenoidal field, different for each component c and for s (the state 0, the target 1).
double field (double x, double y, int k, int c, int s)
{
    return (1.0 + 0.5 * c - 0.3 * s) + 0.2 * k + 0.003 * x * (1 + c) - 0.002 * y * (1 + s)
         + 1.0e-6 * (1 + 0.1 * k) * x * y + 0.4 * std::sin(0.004 * x + 0.3 * c + s) * std::cos(0.003 * y);
}

amrex::MultiFab* fld (gust::State& q, int c)
{
    amrex::MultiFab* f[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return f[c];
}

// u, v, w, p, theta, k of q (ghosts included) = field(., s) (offset makes the parent distinct from the nest).
void fill_state (const amrex::Geometry& g, gust::State& q, int s, double offset = 0)
{
    for (int c = 0; c < 6; ++c) {
        amrex::MultiFab* m = fld(q, c);
        const amrex::IndexType t = m->ixType();
        // A face on a periodic high side is its low side's image (the same value, as in a solver state).
        const amrex::Box d = g.Domain();
        auto wrap = [&] (int i, int dir) {
            return (g.isPeriodic(dir) && t.nodeCentered(dir) && i == d.bigEnd(dir) + 1) ? d.smallEnd(dir) : i;
        };
        fill_all(*m, [&] (int i, int j, int k) {
            return offset + field(pos(g, t, wrap(i, 0), 0), pos(g, t, wrap(j, 1), 1), k, c, s);
        });
    }
}

void copy_state (const gust::State& a, gust::State& b)
{
    const amrex::MultiFab* s[7] = {&a.u, &a.v, &a.w, &a.p, &a.theta, &a.tke, &a.z_phys_nd};
    amrex::MultiFab* d[7] = {&b.u, &b.v, &b.w, &b.p, &b.theta, &b.tke, &b.z_phys_nd};
    for (int n = 0; n < 7; ++n) { amrex::MultiFab::Copy(*d[n], *s[n], 0, 0, 1, s[n]->nGrow()); }
}

std::function<void(gust::State&)> ghosts (const amrex::Geometry& g)
{
    return [g] (gust::State& s) {
        for (amrex::MultiFab* m : {&s.u, &s.v, &s.w, &s.p, &s.theta, &s.tke}) { m->FillBoundary(g.periodicity()); }
    };
}

// w = 0 on the ground and top z faces (ghosts included), as a Dirichlet condition there sets it: the injection does
// not write those faces (NestInject.cpp, data_box), so the divergence checks need fields that obey the condition.
void wall_w (gust::State& q, const amrex::Geometry& g)
{
    const int lo = g.Domain().smallEnd(2), hi = g.Domain().bigEnd(2) + 1;
    for (amrex::MFIter mfi(q.w); mfi.isValid(); ++mfi) {
        const auto w = q.w.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k == lo || k == hi) { w(i, j, k) = amrex::Real(0); }
        });
    }
}

amrex::MultiFab jdiv (const gust::State& q, Case& r)
{
    amrex::MultiFab div(r.ba, r.dm, 1, 0);
    gust::velocity_divergence(q.u, q.v, q.w, q.z_phys_nd, r.mesh, r.geom, div);
    amrex::MultiFab::Multiply(div, r.mesh.detJ, 0, 0, 1, 0);
    return div;
}

} // namespace

TEST(NestInject, InjectsTheRestrictedNestOnTheInnerBoxOnAHill)   // UT-I.1
{
    ASSERT_EQ(gust::kDefectMargin, 2);
    Case p{parent_inputs(true)};
    p.make();
    Case n{nest_inputs(p.in)};
    n.make(true);
    fill_state(p.geom, p.q, 0, 10.0);
    fill_state(n.geom, n.q, 1);
    wall_w(p.q, p.geom);
    wall_w(n.q, n.geom);
    const amrex::Box dom = p.geom.Domain();
    const gust::DefectRegion g = gust::defect_region("nest.a", n.geom.Domain().size(), kR, kOffset, dom);
    // The inner box: footprint cells 3-14 x 2-11, two cells in from each face: 5-12 x 4-9, all z.
    ASSERT_EQ(g.inner, amrex::Box(amrex::IntVect(5, 4, 0), amrex::IntVect(12, 9, 7)));
    gust::State qstar, q;
    qstar.define(p.ba, p.dm);
    q.define(p.ba, p.dm);
    copy_state(p.q, qstar);
    copy_state(p.q, q);
    gust::restrict_patch(n.q, n.mesh, n.geom, kR, kOffset, p.mesh, p.geom, ghosts(p.geom), qstar);

    gust::NestInject I(p.ba, p.dm, p.geom);
    EXPECT_FALSE(I.active());
    I.inject(qstar, q);   // inactive: nothing
    for (int c = 0; c < 6; ++c) {
        EXPECT_EQ(gust::test::count_bit_differences(*fld(q, c), *fld(p.q, c), 0), 0) << c << ": inactive injection";
    }
    I.set({g});
    EXPECT_TRUE(I.active());
    I.inject(qstar, q, true);   // full: the final injection after the joint loop
    ghosts(p.geom)(q);

    // Points on the inner box: 8 x 6 x 8 cells; x-faces 9 x 6 x 8, y-faces 8 x 7 x 8, z-faces 8 x 6 x 9.
    const long expect_on[6] = {9 * 6 * 8, 8 * 7 * 8, 8 * 6 * 9, 8 * 6 * 8, 8 * 6 * 8, 8 * 6 * 8};
    const char* names[6] = {"u", "v", "w", "p", "theta", "tke"};
    for (int c = 0; c < 6; ++c) {
        const amrex::Box b = amrex::convert(dom, fld(q, c)->ixType());
        const amrex::Box inner = amrex::convert(g.inner, fld(q, c)->ixType());
        const std::vector<double> a = gather(*fld(q, c), dom), s = gather(*fld(qstar, c), dom), o = gather(*fld(p.q, c), dom);
        long on = 0, wrong_on = 0, wrong_off = 0, differs = 0, walls = 0;
        for (amrex::BoxIterator it(b); it.ok(); ++it) {
            const auto m = static_cast<std::size_t>(b.index(it()));
            if (inner.contains(it())) {
                ++on;
                // w's ground and top faces are the boundary condition's, not data: the injection keeps q's (on the
                // terrain the divergence reads no w there, so the J div checks below hold).
                if (c == 2 && (it()[2] == b.smallEnd(2) || it()[2] == b.bigEnd(2))) {
                    ++walls;
                    if (a[m] != o[m]) { ++wrong_off; }
                    continue;
                }
                if (a[m] != s[m]) { ++wrong_on; }
                if (s[m] != o[m]) { ++differs; }
            } else if (a[m] != o[m]) {
                ++wrong_off;
            }
        }
        EXPECT_EQ(on, expect_on[c]) << names[c];
        EXPECT_EQ(wrong_on, 0) << names[c] << ": an inner point is not q*'s";
        EXPECT_EQ(wrong_off, 0) << names[c] << ": a point outside the inner box changed";
        EXPECT_EQ(walls, c == 2 ? 2 * 8 * 6 : 0) << names[c];
        EXPECT_EQ(differs, on - walls) << names[c] << ": q* equals the parent somewhere on the inner box: the test is vacuous";
    }

    // The exchange's injection: the default nest.feedback_factor (0.5) moves the inner box halfway to q*, nothing else.
    {
        gust::State qh;
        qh.define(p.ba, p.dm);
        copy_state(p.q, qh);
        gust::NestInject H(p.ba, p.dm, p.geom);
        ASSERT_EQ(H.factor(), gust::kInjectFactorDefault);
        H.set({g});
        H.inject(qstar, qh);
        for (int c = 0; c < 6; ++c) {
            const amrex::Box b = amrex::convert(dom, fld(qh, c)->ixType());
            const amrex::Box inner = amrex::convert(g.inner, fld(qh, c)->ixType());
            const std::vector<double> a = gather(*fld(qh, c), dom), s = gather(*fld(qstar, c), dom),
                                      o = gather(*fld(p.q, c), dom);
            long wrong_on = 0, wrong_off = 0;
            for (amrex::BoxIterator it(b); it.ok(); ++it) {
                const auto m = static_cast<std::size_t>(b.index(it()));
                if (c == 2 && inner.contains(it()) && (it()[2] == b.smallEnd(2) || it()[2] == b.bigEnd(2))) {
                    if (a[m] != o[m]) { ++wrong_off; }   // w's ground and top faces: q's (above)
                } else if (inner.contains(it())) {
                    const double half = o[m] + 0.5 * (s[m] - o[m]);
                    if (std::abs(a[m] - half) > 1e-6 * (std::abs(half) + 1)) { ++wrong_on; }
                } else if (a[m] != o[m]) {
                    ++wrong_off;
                }
            }
            EXPECT_EQ(wrong_on, 0) << names[c] << ": an inner point is not halfway to q*";
            EXPECT_EQ(wrong_off, 0) << names[c] << ": a point outside the inner box changed (blend)";
        }
    }

    // Divergence: the inner cells carry q*'s J div (the nest's 2 x 2 mean, restrict_patch's conservation); outside,
    // only the cells sharing a face with the inner box change.
    const std::vector<double> dq = gather(jdiv(q, p), dom), ds = gather(jdiv(qstar, p), dom), dp = gather(jdiv(p.q, p), dom);
    const std::vector<double> df = gather(jdiv(n.q, n), n.geom.Domain());
    const amrex::Box dom_f = n.geom.Domain();
    double err_star = 0, err_nest = 0, scale = 0;
    long ring_changed = 0, far_changed = 0, ground_cells = 0;
    for (amrex::BoxIterator it(dom); it.ok(); ++it) {
        const amrex::IntVect c = it();
        const auto m = static_cast<std::size_t>(dom.index(c));
        if (g.inner.contains(c)) {
            double s = 0;
            for (int b = 0; b < kR; ++b) {
                for (int a = 0; a < kR; ++a) {
                    s += df[static_cast<std::size_t>(dom_f.index(amrex::IntVect(kR * (c[0] - kOffset[0]) + a,
                                                                               kR * (c[1] - kOffset[1]) + b, c[2])))];
                }
            }
            s /= kR * kR;
            // The ground cells: their ground w face is the boundary condition's (on the terrain it follows u, set
            // by the application's ghost fill, not by this test's), so their J div is not compared here; UT-I.4
            // checks the column fluxes with the application's fill.
            if (c[2] == dom.smallEnd(2)) { ++ground_cells; continue; }
            err_star = std::max(err_star, std::abs(dq[m] - ds[m]));
            err_nest = std::max(err_nest, std::abs(dq[m] - s));
            scale = std::max(scale, std::abs(s));
            continue;
        }
        amrex::Box touch = g.inner;
        touch.grow(0, 1);
        const bool xn = touch.contains(c);
        touch = g.inner;
        touch.grow(1, 1);
        const bool yn = touch.contains(c);
        if (xn || yn) {
            if (dq[m] != dp[m]) { ++ring_changed; }
        } else if (dq[m] != dp[m]) {
            ++far_changed;
        }
    }
    EXPECT_EQ(ground_cells, 8 * 6) << "the inner box's ground cells";
    EXPECT_GT(scale, 1.0e-4) << "the nest field is not divergent enough to test the fluxes";
    EXPECT_LT(err_star, kTol * std::max(1.0, scale)) << "the inner cells' J div differs from q*'s";
    EXPECT_LT(err_nest, kTol * std::max(1.0, scale)) << "the inner cells' J div is not the nest's 2 x 2 mean";
    EXPECT_EQ(far_changed, 0) << "a cell not touching the inner box changed its divergence";
    // The ring cells beside the four lateral faces: 2 x 6 x 8 + 2 x 8 x 8 = 224.
    EXPECT_EQ(ring_changed, 224) << "the ring cells beside the inner box";
}

TEST(NestInject, ZeroRowsAndDifferencesOnUnevenBoxes)   // UT-I.2
{
    std::vector<std::vector<amrex::Real>> diffs;
    for (int layout = 0; layout < 3; ++layout) {
        // 0: uneven boxes; 1: the same on reversed ranks; 2: one box. Periodic in y.
        Case p{parent_inputs(false, true, layout == 2 ? 64 : 7)};
        p.make(layout == 1);
        const amrex::Box dom = p.geom.Domain();
        // A: parent cells 3-14 x 2-11; B: cells 15-19 spanning the periodic y (no y faces; inner x cell 17 only).
        const gust::DefectRegion a = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, amrex::IntVect(3, 2, 0), dom);
        const gust::DefectRegion b = gust::defect_region("nest.b", amrex::IntVect(10, 32, 8), 2, amrex::IntVect(15, 0, 0), dom,
                                                         {false, true});
        ASSERT_EQ(b.inner.length(1), dom.length(1));
        const std::vector<gust::DefectRegion> regions{a, b};
        gust::State qstar;
        qstar.define(p.ba, p.dm);
        fill_state(p.geom, p.q, 0);
        fill_state(p.geom, qstar, 1);
        gust::NestInject I(p.ba, p.dm, p.geom);

        gust::ResidualFields R, R0;
        R.define(p.ba, p.dm);
        R0.define(p.ba, p.dm);
        const auto rf = R.fields();
        const auto r0 = R0.fields();
        for (std::size_t c = 0; c < rf.size(); ++c) {
            const amrex::IndexType t = rf[c].mf->ixType();
            fill_all(*rf[c].mf, [&] (int i, int j, int k) { return 3.0 + field(pos(p.geom, t, i, 0), pos(p.geom, t, j, 1), k, int(c), 0); });
            amrex::MultiFab::Copy(*r0[c].mf, *rf[c].mf, 0, 0, 1, 0);
        }
        I.zero_rows(R);   // inactive
        for (std::size_t c = 0; c < rf.size(); ++c) {
            EXPECT_EQ(gust::test::count_bit_differences(*rf[c].mf, *r0[c].mf, 0), 0) << rf[c].name << ": inactive";
        }
        for (const amrex::Real d : I.differences(p.q, qstar)) { EXPECT_EQ(d, 0); }
        I.set(regions);
        I.zero_rows(R);
        for (std::size_t c = 0; c < rf.size(); ++c) {
            const amrex::Box bx = amrex::convert(dom, rf[c].mf->ixType());
            const std::vector<double> got = gather(*rf[c].mf, dom), was = gather(*r0[c].mf, dom);
            long on = 0, wrong = 0, kept_on = 0;
            for (amrex::BoxIterator it(bx); it.ok(); ++it) {
                const auto m = static_cast<std::size_t>(bx.index(it()));
                bool in = false;
                for (const auto& g : regions) { in = in || amrex::convert(g.inner, rf[c].mf->ixType()).contains(it()); }
                if (in) {
                    ++on;
                    if (got[m] != 0.0) { ++wrong; }
                    if (was[m] != 0.0) { ++kept_on; }
                } else if (got[m] != was[m]) {
                    ++wrong;
                }
            }
            EXPECT_EQ(wrong, 0) << rf[c].name << ": R is not 0 on the inner points and unchanged elsewhere";
            EXPECT_GT(on, 0) << rf[c].name;
            EXPECT_EQ(kept_on, on) << rf[c].name << ": R was 0 before: the test is vacuous";
        }
        // differences(): max and RMS of |q - q*| over the inner points, each once.
        const std::vector<amrex::Real> d = I.differences(p.q, qstar);
        ASSERT_EQ(d.size(), 10u);
        const int ci[5] = {0, 1, 2, 4, 5};
        for (int e = 0; e < 5; ++e) {
            const amrex::MultiFab& mf = *fld(p.q, ci[e]);
            const amrex::IndexType t = mf.ixType();
            const amrex::Box bx = amrex::convert(dom, t);
            const std::vector<double> x = gather(mf, dom), y = gather(*fld(qstar, ci[e]), dom);
            double mx = 0, ss = 0;
            long cnt = 0;
            for (amrex::BoxIterator it(bx); it.ok(); ++it) {
                // A periodic face is counted once (the high copy is the low one's image).
                if (t.nodeCentered(1) && it()[1] == bx.bigEnd(1)) { continue; }
                bool in = false;
                for (const auto& g : regions) { in = in || amrex::convert(g.inner, t).contains(it()); }
                if (!in) { continue; }
                const auto m = static_cast<std::size_t>(bx.index(it()));
                const double v = std::abs(x[m] - y[m]);
                mx = std::max(mx, v);
                ss += v * v;
                ++cnt;
            }
            ASSERT_GT(cnt, 0);
            EXPECT_GT(mx, 0.1) << e << ": q - q* is too small for the test";
            EXPECT_NEAR(static_cast<double>(d[static_cast<std::size_t>(e)]), mx, kTol * mx) << e << ": max";
            const double rms = std::sqrt(ss / static_cast<double>(cnt));
            EXPECT_NEAR(static_cast<double>(d[static_cast<std::size_t>(5 + e)]), rms, kTol * rms) << e << ": RMS";
        }
        diffs.push_back(d);
    }
    for (std::size_t m = 0; m < diffs[0].size(); ++m) {
        EXPECT_EQ(diffs[0][m], diffs[1][m]) << m << ": differences() depends on the rank map";
        EXPECT_NEAR(static_cast<double>(diffs[2][m]), static_cast<double>(diffs[0][m]),
                    kTol * static_cast<double>(diffs[0][m])) << m << ": differences() depends on the box layout";
    }
}

namespace {

using gust::hilltest::HillRun;

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

// The hill (16 x 2 x 12 cells, four boxes in x, periodic in x and y), a footprint over cells 3 to 12 in x spanning
// y (no y faces), its inner box 5 to 10.
gust::DefectRegion hill_region ()
{
    return gust::DefectRegion{"nest.a", amrex::Box(amrex::IntVect(3, 0, 0), amrex::IntVect(12, 1, 11)),
                              amrex::Box(amrex::IntVect(5, 0, 0), amrex::IntVect(10, 1, 11))};
}

} // namespace

TEST(NestInject, ResidualRowsOnTheHill)   // UT-I.3
{
    HillRun h(0.003, 16, 12, 4);
    gust::NestInject I(h.c.ba, h.c.dm, h.c.geom);
    gust::State& q = h.q;
    gust::ResidualInputs off = h.ls->res_in;
    ASSERT_EQ(off.inject, nullptr);
    gust::ResidualInputs on = off;
    on.inject = &I;
    gust::ResidualFields Ron, Roff;
    gust::residual(q, on, h.c.mesh, h.c.geom, Ron);
    gust::residual(q, off, h.c.mesh, h.c.geom, Roff);
    {
        const auto a = Ron.fields();
        const auto b = Roff.fields();
        for (std::size_t n = 0; n < a.size(); ++n) {
            EXPECT_EQ(gust::test::count_bit_differences(*a[n].mf, *b[n].mf, 0), 0) << a[n].name << ": inactive";
        }
    }
    I.set({hill_region()});
    gust::residual(q, on, h.c.mesh, h.c.geom, Ron);
    const amrex::Box inner = hill_region().inner;
    const amrex::Box dom = h.c.geom.Domain();
    const amrex::MultiFab* ra[6] = {&Ron.u, &Ron.v, &Ron.w, &Ron.p, &Ron.theta, &Ron.tke};
    const amrex::MultiFab* rb[6] = {&Roff.u, &Roff.v, &Roff.w, &Roff.p, &Roff.theta, &Roff.tke};
    const char* names[6] = {"u", "v", "w", "p", "theta", "tke"};
    for (int n = 0; n < 6; ++n) {
        const amrex::IndexType t = ra[n]->ixType();
        const amrex::Box bx = amrex::convert(dom, t);
        const amrex::Box ib = amrex::convert(inner, t);
        const std::vector<double> a = gather(*ra[n], dom), b = gather(*rb[n], dom);
        long on_pts = 0, wrong = 0, nonzero_before = 0, ground = 0;
        for (amrex::BoxIterator it(bx); it.ok(); ++it) {
            const auto m = static_cast<std::size_t>(bx.index(it()));
            if (ib.contains(it())) {
                ++on_pts;
                if (a[m] != 0.0) { ++wrong; }
                if (b[m] != 0.0) { ++nonzero_before; }
                // The rows the wall law and the boundaries replace: the first cells (k) and the ground faces (w).
                if (it()[2] == 0 && b[m] != 0.0) { ++ground; }
            } else if (a[m] != b[m]) {
                ++wrong;
            }
        }
        EXPECT_EQ(wrong, 0) << names[n] << ": R is not 0 on the inner points and unchanged elsewhere";
        EXPECT_GT(on_pts, 0) << names[n];
        // v is 0 on this 2-D hill (its rows are UT-I.2's); the others must have a residual to zero. The
        // continuity residual of this field is round-off, and how many of its cells are exactly 0 depends on
        // the compiler (60 of 144 nonzero under clang with -ffp-exception-behavior=maytrap), so p needs one.
        if (n == 3) {
            EXPECT_GT(nonzero_before, 0) << names[n] << ": the residual was 0 there: vacuous";
        } else if (n != 1) {
            EXPECT_GT(nonzero_before, on_pts / 2) << names[n] << ": the residual was mostly 0 there: vacuous";
        }
        if (n == 5) { EXPECT_GT(ground, 0) << "no wall-k row with a residual on the inner box: vacuous"; }
    }
    // The implicit rows do not change.
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    for (const Equation eq : eqs) {
        amrex::MultiFab Aon, Aoff;
        rows(h, q, on, eq, Aon);
        rows(h, q, off, eq, Aoff);
        EXPECT_EQ(gust::test::count_bit_differences(Aon, Aoff, 0), 0) << static_cast<int>(eq);
    }
}

namespace {

// The stop test in single precision, where r and the default steady criterion cannot be met on this hill: v is 0
// (2-D) and its residual is round-off from the start (first 4e-5, then 6e-7 to 9e-7, so r stays at 0.02), and k's
// floor is 3e-4, above steady.residual_tol. So each field has its own bound above its measured floor: u, w, p and
// theta to 2e-4 of their reference (the largest of the first five; their floors are at most 2e-5 of it in the
// injected run), k to k_drop of its, and v to 5e-6 itself.
struct SingleStop
{
    double k_drop;
    std::array<double, 6> ref{};    // u, v, w, p, theta, tke (stage_residual_names)
    int count{0};

    bool operator() (const std::vector<amrex::Real>& res)
    {
        if (res.size() != ref.size()) { return false; }
        if (count < 5) {
            for (std::size_t n = 0; n < ref.size(); ++n) { ref[n] = std::max(ref[n], static_cast<double>(res[n])); }
            ++count;
            return false;
        }
        for (const std::size_t n : {0u, 2u, 3u, 4u}) {
            if (!(res[n] < 2e-4 * ref[n])) { return false; }
        }
        return res[1] < 5e-6 && res[5] < k_drop * ref[5];
    }
};

} // namespace

TEST(NestInject, InjectedRunReachesSteady)   // UT-I.4, the stop test
{
    // The hill's state near steady.
    const double lapse = 0.003;
    HillRun h0(lapse, 16, 12, 4);
    if (kDouble) {
        h0.simplec_to(1e-4, 4000);
        ASSERT_LT(h0.r(), 1e-4) << "the hill did not converge";
    } else {
        SingleStop stop{5e-4};  // about 1700 iterations; k's floor is 1.5e-4 of its reference here
        bool done = false;
        for (int n = 0; n < 4000 && !done; ++n) {
            h0.ls->driver().step(h0.q);
            done = stop(h0.ls->stage().residual_norms());
        }
        ASSERT_TRUE(done) << "the hill did not converge";
    }
    HillRun h(lapse, 16, 12, 4, 0, gust::Limiter::van_leer, gust::StageKind::simplec, nullptr, true);
    ASSERT_NE(h.inject, nullptr);
    // The two runs build their own distribution maps: the converged state is copied onto h's layout.
    const amrex::MultiFab* s0[6] = {&h0.q.u, &h0.q.v, &h0.q.w, &h0.q.p, &h0.q.theta, &h0.q.tke};
    amrex::MultiFab* s1[6] = {&h.q.u, &h.q.v, &h.q.w, &h.q.p, &h.q.theta, &h.q.tke};
    for (int n = 0; n < 6; ++n) { s1[n]->ParallelCopy(*s0[n], 0, 0, 1, 0, 0); }
    h.ls->bc.fill_ghosts(h.q, h.ls->surface(), h.s.fopt.buoyancy);
    // A target that is not the parent's solution: the converged velocity times 1.1 on the inner box (its
    // divergence there is 1.1 times the converged one's, which is small).
    gust::State qstar;
    qstar.define(h.c.ba, h.c.dm);
    copy_state(h.q, qstar);
    constexpr double kFactor = 1.1;
    for (amrex::MultiFab* m : {&qstar.u, &qstar.v, &qstar.w}) { m->mult(amrex::Real(kFactor), 0, 1, m->nGrow()); }
    h.inject->set({hill_region()});
    const int limit = 2000;
    const double drop = 1e-6;   // double precision
    SingleStop stop{4e-3};      // single precision: k's floor is 1.8e-3 of its reference in this run
    bool done = false;
    double first = -1, last = 0;
    int it = 0;
    // The held state: the iteration leaves every field on the inner box at the injected value (checked on the
    // first iteration, all six fields, exactly; w's ground and top faces are boundary values).
    const amrex::Box dom = h.c.geom.Domain();
    double held_change = -1;
    for (; it < limit; ++it) {
        h.inject->inject(qstar, h.q);
        h.ls->bc.fill_ghosts(h.q, h.ls->surface(), h.s.fopt.buoyancy);
        std::vector<std::vector<double>> before;
        if (it == 0) {
            for (const amrex::MultiFab* m : s1) { before.push_back(gather(*m, dom)); }
        }
        h.ls->driver().step(h.q);
        if (it == 0) {
            held_change = 0;
            for (int n = 0; n < 6; ++n) {
                const std::vector<double> after = gather(*s1[n], dom);
                const amrex::Box fx = amrex::convert(dom, s1[n]->ixType());
                const amrex::Box in = amrex::convert(hill_region().inner, s1[n]->ixType());
                for (amrex::BoxIterator p(in); p.ok(); ++p) {
                    // w on the ground and top faces is the boundary condition's (on the terrain it follows u).
                    if (n == 2 && (p()[2] == fx.smallEnd(2) || p()[2] == fx.bigEnd(2))) { continue; }
                    const auto m = static_cast<std::size_t>(fx.index(p()));
                    held_change = std::max(held_change, std::abs(after[m] - before[static_cast<std::size_t>(n)][m]));
                }
            }
        }
        const std::vector<amrex::Real> res = h.ls->stage().residual_norms();
        last = static_cast<double>(*std::max_element(res.begin(), res.end()));
        if (first < 0) { first = last; }
        done = kDouble ? last < drop * first : stop(res);
        if (done) { break; }
    }
    EXPECT_EQ(held_change, 0.0) << "the iteration moved the held state on the inner box";
    EXPECT_GT(first, 0.0);
    EXPECT_TRUE(done) << "the injected run did not reach steady in " << limit << " iterations (residual "
                                  << last << " from " << first << ")";
    // The state outside the inner box moved from the converged one: the injection reaches it.
    const amrex::Box ux = amrex::convert(dom, h.q.u.ixType());
    const std::vector<double> a = gather(h.q.u, dom), b = gather(h0.q.u, dom);
    double moved = 0;
    for (amrex::BoxIterator p(ux); p.ok(); ++p) {
        if (amrex::convert(hill_region().inner, h.q.u.ixType()).contains(p())) { continue; }
        moved = std::max(moved, std::abs(a[static_cast<std::size_t>(ux.index(p()))] - b[static_cast<std::size_t>(ux.index(p()))]));
    }
    EXPECT_GT(moved, 0.05) << "the state outside the inner box did not respond to the injection";
    // Closed form: the channel is periodic in x and y and closed at the ground and top, so at steady every
    // x-face column carries the same volume flux, the injected one: kFactor times the converged hill's (the
    // inner box holds kFactor u0 on its faces). Pinned on every column, the inner box's included.
    const std::vector<double> ax = gather(h.c.mesh.ax, dom);
    double worst = 0, q_ref = 0;
    for (int i = ux.smallEnd(0); i <= ux.bigEnd(0); ++i) {
        double qa = 0, q0 = 0;
        for (int j = ux.smallEnd(1); j <= ux.bigEnd(1); ++j) {
            for (int k = ux.smallEnd(2); k <= ux.bigEnd(2); ++k) {
                const auto m = static_cast<std::size_t>(ux.index(amrex::IntVect(i, j, k)));
                qa += a[m] * ax[m];
                q0 += b[m] * ax[m];
            }
        }
        if (i == ux.smallEnd(0)) { q_ref = kFactor * q0; }
        worst = std::max(worst, std::abs(qa - kFactor * q0) / std::abs(kFactor * q0));
    }
    EXPECT_GT(std::abs(q_ref), 0.0);
    EXPECT_LT(worst, 1e-6) << "a column's volume flux is not the injected one";
    amrex::Print() << "NestInject stop test: " << it + 1 << " iterations, residual " << last << " from " << first
                   << ", largest change of u outside the inner box " << moved << " m/s, worst column flux error "
                   << worst << "\n";
}

//! UT-I.5: nest.feedback_factor with inject lies in (0, 1]; the message names the key and the feedback.
TEST(NestInject, FactorCheckNamesTheKey)
{
    EXPECT_TRUE(gust::check_inject_factor(amrex::Real(1)).empty());
    EXPECT_TRUE(gust::check_inject_factor(amrex::Real(0.2)).empty());
    EXPECT_TRUE(gust::check_inject_factor(gust::kInjectFactorDefault).empty());
    const amrex::Real bad[] = {0, amrex::Real(-0.5), amrex::Real(1.01), amrex::Real(2),
                               std::numeric_limits<amrex::Real>::infinity(),
                               std::numeric_limits<amrex::Real>::quiet_NaN()};
    for (const amrex::Real f : bad) {
        const std::string err = gust::check_inject_factor(f);
        EXPECT_NE(err.find("nest.feedback_factor"), std::string::npos) << f;
        EXPECT_NE(err.find("inject"), std::string::npos) << f;
    }
}

//! UT-I.6: the Poisson solve with NestInject's masks (FaceCoeffs::held_fluid, held_open). The held cells are rows
//! phi = 0, no flux crosses a held face, the fluid cells satisfy A phi = rhs minus its mean over the fluid cells
//! alone, and the mean removed is that mean. On a bump, 4 boxes the held box crosses, both GMRES back ends.
TEST(NestInject, PoissonHeldCellsAreData)
{
    for (const char* backend : {"gmres_fft", "gmres_mlmg"}) {
        gust::ptest::PCase c;
        c.in.mesh.n_cell = amrex::IntVect(16, 16, 8);
        c.in.mesh.prob_hi = {amrex::Real(160), amrex::Real(160), amrex::Real(80)};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.in.mesh.max_grid_size = {8, 8};
        c.in.terrain.smoothing = 0;
        c.make([] (double x, double y) { return 6.0 * std::exp(-((x - 80) * (x - 80) + (y - 80) * (y - 80)) / 1600.0); });
        const amrex::Box inner(amrex::IntVect(5, 4, 0), amrex::IntVect(10, 11, 7));
        gust::NestInject I(c.ba, c.dm, c.geom);
        I.set({{"nest.a", amrex::grow(amrex::grow(inner, 0, gust::kDefectMargin), 1, gust::kDefectMargin), inner}});
        gust::PoissonOptions opt;
        opt.backend = backend;
        gust::PoissonSolver s(c.mesh, c.geom, c.z, opt);
        amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1), keep(c.ba, c.dm, 1, 0);
        // Non-zero on the held cells too, with a mean: the solve masks and removes it.
        gust::ptest::fill_cells(rhs, c, [] (double x, double y, double z) {
            return std::cos(3.14159265358979 * x / 80) * std::sin(3.14159265358979 * y / 40) + 0.01 * z + 0.3; });
        amrex::MultiFab::Copy(keep, rhs, 0, 0, 1, 0);
        phi.setVal(0);
        amrex::Array<amrex::MultiFab, 3> F;
        F[0].define(c.mesh.ax.boxArray(), c.dm, 1, 0);
        F[1].define(c.mesh.ay.boxArray(), c.dm, 1, 0);
        F[2].define(c.mesh.az.boxArray(), c.dm, 1, 0);
        gust::FaceCoeffs fc;
        fc.held_fluid = &I.fluid_mask();
        fc.held_open = {&I.open_mask(0), &I.open_mask(1), &I.open_mask(2)};
        gust::PoissonBC bc;
        for (const int d : {0, 1}) {
            bc.face[amrex::Orientation(d, amrex::Orientation::low)] = gust::PoissonFace::periodic;
            bc.face[amrex::Orientation(d, amrex::Orientation::high)] = gust::PoissonFace::periodic;
        }
        const gust::PoissonResult res = s.solve(rhs, phi, fc, bc, &F);
        amrex::MultiFab div(c.ba, c.dm, 1, 0);
        gust::face_divergence(c.mesh, c.geom, {&F[0], &F[1], &F[2]}, div);
        // The mean over the fluid cells (J-weighted, as the solve removes it).
        double sum = 0, vol = 0;
        for (amrex::MFIter mfi(keep, false); mfi.isValid(); ++mfi) {
            const auto a = keep.const_array(mfi);
            const auto J = c.mesh.detJ.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (!inner.contains(amrex::IntVect(i, j, k))) { sum += a(i, j, k) * J(i, j, k); vol += J(i, j, k); }
            });
        }
        amrex::ParallelDescriptor::ReduceRealSum(sum);
        amrex::ParallelDescriptor::ReduceRealSum(vol);
        const double mean = sum / vol;
        EXPECT_NEAR(res.removed_mean, mean, kDouble ? 1e-9 : 1e-6) << backend << ": the mean removed is not the fluid cells' mean";
        amrex::Long held_phi = 0, fluid_bad = 0, held_flux = 0;   // ReduceLongSum takes amrex::Long (MSVC)
        double worst = 0;
        for (amrex::MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
            const auto p = phi.const_array(mfi);
            const auto dv = div.const_array(mfi);
            const auto a = keep.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (inner.contains(amrex::IntVect(i, j, k))) {
                    if (p(i, j, k) != 0) { ++held_phi; }
                } else {
                    // A phi = -div F (face_divergence of the solve's fluxes).
                    const double e = std::abs(dv(i, j, k) + (a(i, j, k) - mean));
                    worst = std::max(worst, e);
                    if (e > (kDouble ? 1e-6 : 1e-4)) { ++fluid_bad; }   // single precision: at most 2.5e-5
                }
            });
        }
        for (int d = 0; d < 3; ++d) {
            const amrex::Box fb = amrex::convert(inner, amrex::IntVect::TheDimensionVector(d));
            for (amrex::MFIter mfi(F[d], false); mfi.isValid(); ++mfi) {
                const auto f = F[d].const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox() & fb, [&] (int i, int j, int k) { if (f(i, j, k) != 0) { ++held_flux; } });
            }
        }
        amrex::ParallelDescriptor::ReduceLongSum(held_phi);
        amrex::ParallelDescriptor::ReduceLongSum(fluid_bad);
        amrex::ParallelDescriptor::ReduceLongSum(held_flux);
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        EXPECT_EQ(held_phi, 0) << backend << ": a held cell's phi is not 0";
        EXPECT_EQ(held_flux, 0) << backend << ": a flux crosses a held face";
        EXPECT_EQ(fluid_bad, 0) << backend << ": a fluid cell does not satisfy A phi = rhs - mean (worst " << worst << ")";
    }
}

//! UT-I.7: the change of the target that the joint run's stop window reads (change(), FeedbackChange). With the
//! state q and the restricted nest q* fixed, the first inject moves the target omega (q* - q) from q and the
//! second omega (1 - omega) (q* - q) from the first: max over the inner data points, computed here from gathered
//! fields. The size is the largest velocity component of the new target for u, v and w, 0 for theta and the
//! largest k. Never the first exchange (a restart measures it the same way). Uneven boxes, reversed ranks.
TEST(NestInject, TargetChangeIsTheWindowMeasure)   // UT-I.7
{
    constexpr double kOmega = 0.4;
    for (int layout = 0; layout < 2; ++layout) {
        Case p{parent_inputs(false, true, 7)};
        p.make(layout == 1);
        const amrex::Box dom = p.geom.Domain();
        const gust::DefectRegion a = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, amrex::IntVect(3, 2, 0), dom);
        gust::State qstar;
        qstar.define(p.ba, p.dm);
        fill_state(p.geom, p.q, 0);
        fill_state(p.geom, qstar, 1);
        gust::NestInject I(p.ba, p.dm, p.geom, amrex::Real(kOmega));
        I.set({a});
        EXPECT_FALSE(I.change().first);
        const int ci[5] = {0, 1, 2, 4, 5};
        // d[e]: max |q* - q| and qmax[e]: max |q + s (q* - q)| over the inner data points (w less its ground and
        // top faces), for the step s of each inject.
        auto measure = [&] (double s, std::array<double, 5>& d, std::array<double, 5>& t) {
            for (int e = 0; e < 5; ++e) {
                const amrex::MultiFab& mf = *fld(p.q, ci[e]);
                const amrex::IndexType ty = mf.ixType();
                const amrex::Box bx = amrex::convert(dom, ty);
                const std::vector<double> x = gather(mf, dom), y = gather(*fld(qstar, ci[e]), dom);
                d[e] = 0;
                t[e] = 0;
                for (amrex::BoxIterator it(amrex::convert(a.inner, ty)); it.ok(); ++it) {
                    if (ty.nodeCentered(2) && (it()[2] == bx.smallEnd(2) || it()[2] == bx.bigEnd(2))) { continue; }
                    const auto m = static_cast<std::size_t>(bx.index(it()));
                    d[e] = std::max(d[e], std::abs(y[m] - x[m]));
                    t[e] = std::max(t[e], std::abs(x[m] + s * (y[m] - x[m])));
                }
            }
        };
        for (int pass = 0; pass < 2; ++pass) {
            std::array<double, 5> d{}, t{};
            measure(kOmega, d, t);
            I.inject(qstar, p.q);
            const gust::FeedbackChange& c = I.change();
            EXPECT_FALSE(c.first);
            for (int e = 0; e < 5; ++e) {
                ASSERT_GT(d[e], 0.1) << e << ": q* - q is too small for the test";
                const double want = kOmega * d[e];
                EXPECT_NEAR(static_cast<double>(c.change[static_cast<std::size_t>(e)]), want, kTol * want)
                    << "pass " << pass << ", field " << e << ": the target's change";
            }
            const double U = std::max(t[0], std::max(t[1], t[2]));
            for (int e = 0; e < 3; ++e) {
                EXPECT_NEAR(static_cast<double>(c.size[static_cast<std::size_t>(e)]), U, kTol * U) << "pass " << pass;
            }
            EXPECT_EQ(c.size[3], 0) << "theta has no size (its floor is absolute)";
            EXPECT_NEAR(static_cast<double>(c.size[4]), t[4], kTol * t[4]) << "pass " << pass << ": k's size";
        }
    }
}
