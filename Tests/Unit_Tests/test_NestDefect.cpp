// The restriction and the defect source of the joint two-level run
// (NestDefect.H; plan rev. 6 tracker, PR-70, nest.feedback = defect).
//   restrict_patch, on a hill parent with an uneven patch layout:
//     - flux conservation: on every footprint cell the parent's J div after
//       the restriction is the mean of J div over its 2 x 2 nest cells, for
//       an arbitrary nest field (u, v and Omega restricted as fluxes); a
//       divergence-free nest field therefore gives a divergence-free parent
//       footprint;
//     - theta linear in x and y on flat ground: the parent cell's centre
//       value (the mean of the four nest cells), exactly;
//     - k constant on the hill: unchanged by the detJ weights;
//     - p: the detJ-weighted mean of the nest cells shifted by one constant, so that
//       its footprint mean is the parent's footprint mean before;
//     - outside the footprint every field keeps its values bit for bit;
//     - a footprint on a periodic seam sets both copies of the seam face;
//     - u, v, theta, k do not depend on the nest's box layout or rank map
//       (bit for bit), on 1, 2 and 4 ranks (the unit_np* tests);
//   defect_region: the footprint and the inner box kDefectMargin cells in;
//   NestDefect: S = R on the inner box's points (faces from L to H + 1 along
//     their normal), 0 elsewhere; add_to_residual subtracts it; nothing
//     while inactive.
#include "Inputs.H"
#include "MassFlux.H"
#include "Mesh.H"
#include "NestDefect.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>

namespace {

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double kTol = kDouble ? 1.0e-11 : 2.0e-4;

struct DefectCase
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

// The parent: 20 x 16 x 8 cells of 100 m, stretched layers, uneven boxes;
// hill: a round Agnesi hill at the parent's centre (BTF, as nests need).
gust::Inputs parent_inputs (bool hill)
{
    gust::Inputs in;
    in.mesh.n_cell = amrex::IntVect(20, 16, 8);
    in.mesh.prob_lo = {-1000, 500, 0};
    in.mesh.prob_hi = {1000, 2100, 800};
    in.mesh.is_periodic = {0, 0, 0};
    in.mesh.max_grid_size = {7, 6};
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

// The nest over parent cells 3 to 14 in x and 2 to 11 in y: 24 x 20 cells of
// 50 m, the parent's z-levels, boxes of at most mgs x (mgs - 2) cells (even:
// they must coarsen by 2).
gust::Inputs nest_inputs (const gust::Inputs& parent, int mgs = 8)
{
    gust::Inputs in = parent;
    in.mesh.n_cell = amrex::IntVect(24, 20, 8);
    in.mesh.prob_lo = {-700, 700, 0};
    in.mesh.prob_hi = {500, 1700, 800};
    in.mesh.max_grid_size = {mgs, mgs - 2};
    in.terrain.hill_centre = {0, 1300};  // the parent's centre, as App.cpp sets it
    return in;
}

const amrex::IntVect kOffset(3, 2, 0);
constexpr int kR = 2;

// Every point of mf (ghosts included) = fn(i, j, k), through a host copy.
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

// fn(i, j, k, value) on every valid point of mf, through a host copy.
void each_valid (const amrex::MultiFab& mf, const std::function<void(int, int, int, double)>& fn)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox h(mf[mfi].box(), 1, amrex::The_Pinned_Arena());
        h.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { fn(i, j, k, static_cast<double>(h(amrex::IntVect(i, j, k)))); });
    }
}

// The whole of mf (valid points) on every rank, indexed by its domain box.
std::vector<double> gather (const amrex::MultiFab& mf, const amrex::Box& dom)
{
    const amrex::Box b = amrex::convert(dom, mf.ixType());
    std::vector<double> out(static_cast<std::size_t>(b.numPts()), 0.0);
    std::vector<double> cnt(out.size(), 0.0);
    each_valid(mf, [&] (int i, int j, int k, double v) {
        const auto n = static_cast<std::size_t>(b.index(amrex::IntVect(i, j, k)));
        out[n] += v;
        cnt[n] += 1.0;
    });
    // A face shared by two boxes is a valid point of both: average the copies.
    amrex::ParallelDescriptor::ReduceRealSum(out.data(), static_cast<int>(out.size()));
    amrex::ParallelDescriptor::ReduceRealSum(cnt.data(), static_cast<int>(cnt.size()));
    for (std::size_t n = 0; n < out.size(); ++n) { if (cnt[n] > 0.0) { out[n] /= cnt[n]; } }
    return out;
}

double pos (const amrex::Geometry& g, const amrex::MultiFab& mf, int i, int d)
{
    const double off = mf.ixType().nodeCentered(d) ? 0.0 : 0.5;
    return static_cast<double>(g.ProbLo(d)) + (i + off) * static_cast<double>(g.CellSize(d));
}

// A smooth, non-solenoidal field: different for each component c.
double field (double x, double y, int k, int c)
{
    return (1.0 + 0.5 * c) + 0.2 * k + 0.003 * x * (1 + c) - 0.002 * y + 1.0e-6 * (1 + 0.1 * k) * x * y
         + 0.4 * std::sin(0.004 * x + 0.3 * c) * std::cos(0.003 * y);
}

// The parent's state: its own values everywhere (ghosts included), distinct from the nest's.
void fill_parent (DefectCase& p)
{
    gust::State& q = p.q;
    amrex::MultiFab* f[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    for (int c = 0; c < 6; ++c) {
        const amrex::MultiFab* m = f[c];
        fill_all(*f[c], [&] (int i, int j, int k) { return 10.0 - field(pos(p.geom, *m, i, 0), pos(p.geom, *m, j, 1), k, c); });
    }
}

void fill_nest (DefectCase& n, const std::function<double(const amrex::MultiFab&, int, int, int, int)>& fn)
{
    gust::State& q = n.q;
    amrex::MultiFab* f[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    for (int c = 0; c < 6; ++c) {
        const amrex::MultiFab* m = f[c];
        fill_all(*f[c], [&] (int i, int j, int k) { return fn(*m, c, i, j, k); });
    }
}

// The parent's ghost fill for restrict_patch: the interior ghosts; the footprint keeps two cells from the domain.
std::function<void(gust::State&)> ghosts (const amrex::Geometry& g)
{
    return [g] (gust::State& s) {
        for (amrex::MultiFab* m : {&s.u, &s.v, &s.w, &s.p, &s.theta, &s.tke}) { m->FillBoundary(g.periodicity()); }
    };
}

amrex::Box footprint (const DefectCase& p)
{
    return amrex::Box(amrex::IntVect(kOffset[0], kOffset[1], 0),
                      amrex::IntVect(kOffset[0] + 11, kOffset[1] + 9, p.geom.Domain().bigEnd(2)));
}

// J div of state q on its cells (u, v ghosts filled).
amrex::MultiFab jdiv (DefectCase& r)
{
    amrex::MultiFab div(r.ba, r.dm, 1, 0);
    gust::velocity_divergence(r.q.u, r.q.v, r.q.w, r.q.z_phys_nd, r.mesh, r.geom, div);
    amrex::MultiFab::Multiply(div, r.mesh.detJ, 0, 0, 1, 0);
    return div;
}

} // namespace

TEST(NestDefect, RestrictionConservesTheFluxesOnAHill)
{
    DefectCase p{parent_inputs(true)};
    p.make();
    DefectCase n{nest_inputs(p.in)};
    n.make(true);
    fill_parent(p);
    fill_nest(n, [&] (const amrex::MultiFab& m, int c, int i, int j, int k) {
        return field(pos(n.geom, m, i, 0), pos(n.geom, m, j, 1), k, c);
    });
    const amrex::Box dom_c = p.geom.Domain();
    const std::vector<double> before_u = gather(p.q.u, dom_c), before_p = gather(p.q.p, dom_c);
    gust::restrict_patch(n.q, n.mesh, n.geom, kR, kOffset, p.mesh, p.geom, ghosts(p.geom), p.q);

    const std::vector<double> dc = gather(jdiv(p), dom_c);
    const std::vector<double> df = gather(jdiv(n), n.geom.Domain());
    const amrex::Box fp = footprint(p), dom_f = n.geom.Domain();
    double err = 0, scale = 0;
    for (amrex::BoxIterator it(fp); it.ok(); ++it) {
        const amrex::IntVect c = it();
        double s = 0;
        for (int b = 0; b < kR; ++b) {
            for (int a = 0; a < kR; ++a) {
                const amrex::IntVect f(kR * (c[0] - kOffset[0]) + a, kR * (c[1] - kOffset[1]) + b, c[2]);
                s += df[static_cast<std::size_t>(dom_f.index(f))];
            }
        }
        s /= kR * kR;
        err = std::max(err, std::abs(dc[static_cast<std::size_t>(dom_c.index(c))] - s));
        scale = std::max(scale, std::abs(s));
    }
    EXPECT_GT(scale, 1.0e-4) << "the nest field is not divergent enough to test the fluxes";
    EXPECT_LT(err, kTol * std::max(1.0, scale)) << "the parent's J div differs from the nest's mean on the footprint";

    // Outside the footprint u and p keep their values bit for bit.
    const std::vector<double> after_u = gather(p.q.u, dom_c), after_p = gather(p.q.p, dom_c);
    const amrex::Box box_u = amrex::convert(dom_c, p.q.u.ixType());
    const amrex::Box fpu = amrex::convert(fp, p.q.u.ixType());
    amrex::Long changed = 0;
    for (amrex::BoxIterator it(box_u); it.ok(); ++it) {
        if (fpu.contains(it())) { continue; }
        const auto m = static_cast<std::size_t>(box_u.index(it()));
        if (after_u[m] != before_u[m]) { ++changed; }
    }
    for (amrex::BoxIterator it(dom_c); it.ok(); ++it) {
        if (fp.contains(it())) { continue; }
        const auto m = static_cast<std::size_t>(dom_c.index(it()));
        if (after_p[m] != before_p[m]) { ++changed; }
    }
    EXPECT_EQ(changed, 0) << "points outside the footprint changed";
}

// A parent periodic in x and a nest whose footprint starts at its face 0: face 20 is the same face, and the
// restriction sets both (a ParallelCopy without the periodicity sets face 0 only).
TEST(NestDefect, RestrictionSetsBothCopiesOfAPeriodicSeamFace)
{
    gust::Inputs parent_in = parent_inputs(false);
    parent_in.mesh.is_periodic = {1, 0, 0};
    DefectCase p{parent_in};
    p.make();
    gust::Inputs nest_in = nest_inputs(p.in);
    nest_in.mesh.is_periodic = {0, 0, 0};
    nest_in.mesh.prob_lo = {-1000, 700, 0};
    nest_in.mesh.prob_hi = {200, 1700, 800};
    DefectCase n{nest_in};
    n.make();
    fill_parent(p);
    fill_nest(n, [&] (const amrex::MultiFab& m, int c, int i, int j, int k) {
        return field(pos(n.geom, m, i, 0), pos(n.geom, m, j, 1), k, c);
    });
    const amrex::IntVect off(0, 2, 0);
    gust::restrict_patch(n.q, n.mesh, n.geom, kR, off, p.mesh, p.geom, ghosts(p.geom), p.q);
    const amrex::Box dom_c = p.geom.Domain();
    const std::vector<double> u = gather(p.q.u, dom_c);
    const amrex::Box box_u = amrex::convert(dom_c, p.q.u.ixType());
    const int nx = dom_c.length(0);
    double diff = 0, scale = 0;
    for (int k = 0; k <= dom_c.bigEnd(2); ++k) {
        for (int j = off[1]; j <= off[1] + 9; ++j) {
            const double a = u[static_cast<std::size_t>(box_u.index(amrex::IntVect(0, j, k)))];
            const double b = u[static_cast<std::size_t>(box_u.index(amrex::IntVect(nx, j, k)))];
            diff = std::max(diff, std::abs(a - b));
            scale = std::max(scale, std::abs(a));
        }
    }
    EXPECT_GT(scale, 0.1);
    EXPECT_LT(diff, kTol * std::max(1.0, scale)) << "u differs between face 0 and its periodic image, face " << nx;
}

TEST(NestDefect, RestrictionOfCellFieldsAndThePressureLevel)
{
    // Flat ground: theta linear in x and y gives the parent's cell-centre value exactly.
    {
        DefectCase p{parent_inputs(false)};
        p.make();
        DefectCase n{nest_inputs(p.in)};
        n.make();
        fill_parent(p);
        fill_nest(n, [&] (const amrex::MultiFab& m, int c, int i, int j, int k) {
            return c == 4 ? 300.0 + 0.01 * pos(n.geom, m, i, 0) - 0.02 * pos(n.geom, m, j, 1) + 0.5 * k
                          : field(pos(n.geom, m, i, 0), pos(n.geom, m, j, 1), k, c);
        });
        gust::restrict_patch(n.q, n.mesh, n.geom, kR, kOffset, p.mesh, p.geom, ghosts(p.geom), p.q);
        const amrex::Box fp = footprint(p);
        double err = 0;
        each_valid(p.q.theta, [&] (int i, int j, int k, double v) {
            if (!fp.contains(amrex::IntVect(i, j, k))) { return; }
            const double e = 300.0 + 0.01 * pos(p.geom, p.q.theta, i, 0) - 0.02 * pos(p.geom, p.q.theta, j, 1) + 0.5 * k;
            err = std::max(err, std::abs(v - e));
        });
        amrex::ParallelDescriptor::ReduceRealMax(err);
        EXPECT_LT(err, 300.0 * kTol) << "theta is not the mean of the four nest cells";
    }
    // Hill: k constant stays constant; p is the nest's mean shifted to the parent's footprint mean.
    DefectCase p{parent_inputs(true)};
    p.make();
    DefectCase n{nest_inputs(p.in)};
    n.make();
    fill_parent(p);
    fill_nest(n, [&] (const amrex::MultiFab& m, int c, int i, int j, int k) {
        if (c == 5) { return 0.37; }
        if (c == 3) { return 50.0 + field(pos(n.geom, m, i, 0), pos(n.geom, m, j, 1), k, c); }
        return field(pos(n.geom, m, i, 0), pos(n.geom, m, j, 1), k, c);
    });
    const amrex::Box fp = footprint(p), dom_c = p.geom.Domain();
    double mean_before = 0;
    each_valid(p.q.p, [&] (int i, int j, int k, double v) { if (fp.contains(amrex::IntVect(i, j, k))) { mean_before += v; } });
    amrex::ParallelDescriptor::ReduceRealSum(mean_before);
    mean_before /= static_cast<double>(fp.numPts());
    const std::vector<double> pf = gather(n.q.p, n.geom.Domain());
    gust::restrict_patch(n.q, n.mesh, n.geom, kR, kOffset, p.mesh, p.geom, ghosts(p.geom), p.q);

    double kerr = 0;
    each_valid(p.q.tke, [&] (int i, int j, int k, double v) {
        if (fp.contains(amrex::IntVect(i, j, k))) { kerr = std::max(kerr, std::abs(v - 0.37)); }
    });
    amrex::ParallelDescriptor::ReduceRealMax(kerr);
    EXPECT_LT(kerr, kTol) << "a constant k changed under the detJ weights";

    const std::vector<double> pc = gather(p.q.p, dom_c);
    const amrex::Box dom_f = n.geom.Domain();
    const std::vector<double> jf = gather(n.mesh.detJ, dom_f);
    double mean_after = 0, dmin = 1e300, dmax = -1e300, vmin = 1e300, vmax = -1e300;
    for (amrex::BoxIterator it(fp); it.ok(); ++it) {
        const amrex::IntVect c = it();
        const double v = pc[static_cast<std::size_t>(dom_c.index(c))];
        mean_after += v;
        // The shift: the restricted p less the detJ-weighted mean of its nest cells is one constant.
        double sj = 0, sjp = 0;
        for (int b = 0; b < kR; ++b) {
            for (int a = 0; a < kR; ++a) {
                const auto m = static_cast<std::size_t>(
                    dom_f.index(amrex::IntVect(kR * (c[0] - kOffset[0]) + a, kR * (c[1] - kOffset[1]) + b, c[2])));
                sj += jf[m];
                sjp += jf[m] * pf[m];
            }
        }
        dmin = std::min(dmin, v - sjp / sj);
        dmax = std::max(dmax, v - sjp / sj);
        vmin = std::min(vmin, v);
        vmax = std::max(vmax, v);
    }
    mean_after /= static_cast<double>(fp.numPts());
    EXPECT_NEAR(mean_after, mean_before, 1.0e3 * kTol) << "the restricted p's footprint mean is not the parent's";
    EXPECT_LT(dmax - dmin, 1.0e3 * kTol) << "the restricted p is not the nest's mean plus one constant";
    EXPECT_GT(vmax - vmin, 0.1) << "the restricted p lost the nest's variation";
}

// The pressure level of the restriction (one constant per nest, RestrictionOfCellFieldsAndThePressureLevel) does not
// reach S: a constant added to p over the whole footprint leaves the parent's residual unchanged on the source's
// region (kDefectMargin cells in), because the momentum rows there see p only through differences inside the
// footprint. On the footprint's outer cells it does change (the test is not vacuous).
TEST(NestDefect, ThePressureLevelDoesNotReachTheSource)
{
    DefectCase p{parent_inputs(true)};
    p.make();
    fill_parent(p);
    p.q.theta.setVal(300.0);
    p.q.tke.setVal(0.3);
    amrex::MultiFab nu(p.ba, p.dm, 1, 1);
    nu.setVal(2.0);
    gust::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.theta_ref.z = {0};
    in.theta_ref.theta = {300};
    const amrex::Box dom_c = p.geom.Domain();
    const gust::DefectRegion g = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), kR, kOffset, dom_c);
    gust::ResidualFields R0, R1;
    gust::residual(p.q, in, p.mesh, p.geom, R0);
    // p + C on the footprint's cells and on the ghosts below and above them (as restrict_patch then the ghost fill
    // leave it: the ground and top ghosts of p copy the cells next to them).
    constexpr double kC = 1234.5;
    const amrex::Box fp_z = amrex::grow(g.footprint, 2, p.q.p.nGrow());
    for (amrex::MFIter mfi(p.q.p, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox() & fp_z;
        if (b.ok()) { p.q.p[mfi].plus<amrex::RunOn::Device>(amrex::Real(kC), b, 0, 1); }
    }
    amrex::Gpu::streamSynchronize();
    gust::residual(p.q, in, p.mesh, p.geom, R1);
    const amrex::MultiFab* r0[3] = {&R0.u, &R0.v, &R0.w};
    const amrex::MultiFab* r1[3] = {&R1.u, &R1.v, &R1.w};
    double inner = 0, ring = 0, scale = 0;
    for (int n = 0; n < 3; ++n) {
        const std::vector<double> a = gather(*r0[n], dom_c), b = gather(*r1[n], dom_c);
        const amrex::Box all = amrex::convert(dom_c, r0[n]->ixType());
        const amrex::Box in_b = amrex::convert(g.inner, r0[n]->ixType());
        const amrex::Box fp_b = amrex::convert(g.footprint, r0[n]->ixType());
        for (amrex::BoxIterator it(all); it.ok(); ++it) {
            const auto m = static_cast<std::size_t>(all.index(it()));
            const double d = std::abs(b[m] - a[m]);
            scale = std::max(scale, std::abs(a[m]));
            if (in_b.contains(it())) { inner = std::max(inner, d); }
            else if (fp_b.contains(it())) { ring = std::max(ring, d); }
        }
    }
    EXPECT_GT(scale, 0.0);
    EXPECT_LT(inner, kTol * (scale + kC)) << "the pressure level changed the residual on the source's region";
    EXPECT_GT(ring, 1.0) << "the constant did not reach the footprint's edge rows: the test is vacuous";
}

TEST(NestDefect, RestrictionDoesNotDependOnTheNestLayout)
{
    // Each layout builds its own node heights (Mesh::build); they do not depend on the layout
    // (Mesh.NodeHeightsDoNotDependOnTheLayout), so the restriction must not either.
    std::vector<std::vector<double>> res;
    for (int layout = 0; layout < 2; ++layout) {
        DefectCase p{parent_inputs(true)};
        p.make();
        DefectCase n{nest_inputs(p.in, layout == 0 ? 8 : 24)};
        n.make(layout == 1);
        fill_parent(p);
        fill_nest(n, [&] (const amrex::MultiFab& m, int c, int i, int j, int k) {
            return field(pos(n.geom, m, i, 0), pos(n.geom, m, j, 1), k, c);
        });
        gust::restrict_patch(n.q, n.mesh, n.geom, kR, kOffset, p.mesh, p.geom, ghosts(p.geom), p.q);
        std::vector<double> all;
        for (const amrex::MultiFab* m : {&p.q.u, &p.q.v, &p.q.w, &p.q.theta, &p.q.tke}) {
            const std::vector<double> g = gather(*m, p.geom.Domain());
            all.insert(all.end(), g.begin(), g.end());
        }
        res.push_back(all);
    }
    ASSERT_EQ(res[0].size(), res[1].size());
    long differ = 0;
    for (std::size_t m = 0; m < res[0].size(); ++m) { if (res[0][m] != res[1][m]) { ++differ; } }
    EXPECT_EQ(differ, 0) << "the restriction depends on the nest's boxes or ranks";
}

TEST(NestDefect, RegionAndSource)
{
    const amrex::Box dom(amrex::IntVect(0, 0, 0), amrex::IntVect(19, 15, 7));
    const gust::DefectRegion g = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, kOffset, dom);
    EXPECT_EQ(g.footprint, amrex::Box(amrex::IntVect(3, 2, 0), amrex::IntVect(14, 11, 7)));
    EXPECT_EQ(g.inner, amrex::Box(amrex::IntVect(3 + gust::kDefectMargin, 2 + gust::kDefectMargin, 0),
                                  amrex::IntVect(14 - gust::kDefectMargin, 11 - gust::kDefectMargin, 7)));
    // A nest periodic in y (it spans a periodic parent direction) has no y faces: the region keeps the whole span.
    const gust::DefectRegion gp = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, kOffset, dom, {false, true});
    EXPECT_EQ(gp.inner, amrex::Box(amrex::IntVect(3 + gust::kDefectMargin, 2, 0), amrex::IntVect(14 - gust::kDefectMargin, 11, 7)));

    DefectCase p{parent_inputs(false)};
    p.make();
    gust::ResidualFields R, Rsum;
    R.define(p.ba, p.dm);
    Rsum.define(p.ba, p.dm);
    auto rv = [] (int i, int j, int k, int c) { return 1.0 + i + 100.0 * j + 1.0e4 * k + 0.25 * c; };
    int c = 0;
    for (auto& f : R.fields()) {
        const int cc = c++;
        fill_all(*f.mf, [&] (int i, int j, int k) { return rv(i, j, k, cc); });
    }
    gust::NestDefect S(p.ba, p.dm);
    EXPECT_FALSE(S.active());
    c = 0;
    for (auto& f : Rsum.fields()) {
        const int cc = c++;
        fill_all(*f.mf, [&] (int i, int j, int k) { return rv(i, j, k, cc); });
    }
    S.add_to_residual(Rsum);   // inactive: nothing
    S.set(R, {g});
    EXPECT_TRUE(S.active());
    S.add_to_residual(Rsum);   // R - S
    const amrex::MultiFab* src[5] = {&S.source().u, &S.source().v, &S.source().w, &S.source().theta, &S.source().tke};
    const amrex::MultiFab* rs[5] = {&Rsum.u, &Rsum.v, &Rsum.w, &Rsum.theta, &Rsum.tke};
    const int comp[5] = {0, 1, 2, 4, 5};
    amrex::Long bad = 0;
    for (int n = 0; n < 5; ++n) {
        const amrex::Box in = amrex::convert(g.inner, src[n]->ixType());
        each_valid(*src[n], [&] (int i, int j, int k, double v) {
            const bool inside = in.contains(amrex::IntVect(i, j, k));
            if (v != (inside ? rv(i, j, k, comp[n]) : 0.0)) { ++bad; }
        });
        each_valid(*rs[n], [&] (int i, int j, int k, double v) {
            const bool inside = in.contains(amrex::IntVect(i, j, k));
            if (v != (inside ? 0.0 : rv(i, j, k, comp[n]))) { ++bad; }
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(bad);
    EXPECT_EQ(bad, 0) << "S is not R on the inner box's points and 0 elsewhere, or R - S is wrong";
}
