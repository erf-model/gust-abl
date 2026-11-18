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
//   restrict_patches: two nests that touch, with different u on the shared
//     face, keep the restricted nest Omega on both footprints;
//   partly open faces (buildings' face fractions, made up): u, v and w are
//     the nest's flux through the open parts over the parent's open area,
//     0 on a parent face closed in full, against values by hand;
//   defect_region: the footprint and the inner box kDefectMargin cells in;
//   NestDefect: S = R on the inner box's points (faces from L to H + 1 along
//     their normal), 0 elsewhere; add_to_residual subtracts it; nothing
//     while inactive;
//   the change of S between two sets (the joint run's stop test, JointRun.H):
//     max |S_new - S_old| over the inner box's points and the largest
//     |S_new|, against values computed here; first on the first set only;
//     relative_feedback_change's closed form (floor, 0/0, x/0, NaN).
#include "Inputs.H"
#include "MassFlux.H"
#include "Mesh.H"
#include "NestDefect.H"
#include "Residual.H"
#include "State.H"
#include "TestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
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

// Two nests that touch along parent face i = 9 (footprints x 3..8 and 9..14), with different u on the shared face:
// restrict_patches restricts u, v of both before w of either, so on both footprints the parent's Omega is the
// restricted nest Omega (the az-weighted mean of the nest's Omega over its 2 x 2 faces). Restricting the nests
// one after the other (the order before) leaves the first nest's w from its own u on the shared face, which the
// second nest then replaces: on the hill the parent's Omega there is off by z_x times the u difference.
TEST(NestDefect, TouchingNestsKeepTheRestrictedOmegaOnBothFootprints)
{
    DefectCase p{parent_inputs(true)};
    p.make();
    std::vector<std::unique_ptr<DefectCase>> nests;
    const amrex::IntVect offs[2] = {amrex::IntVect(3, 2, 0), amrex::IntVect(9, 2, 0)};
    for (int m = 0; m < 2; ++m) {
        gust::Inputs in = nest_inputs(p.in);
        in.mesh.n_cell = amrex::IntVect(12, 20, 8);
        in.mesh.prob_lo = {amrex::Real(-700.0 + 600.0 * m), 700, 0};
        in.mesh.prob_hi = {amrex::Real(-100.0 + 600.0 * m), 1700, 800};
        nests.push_back(std::make_unique<DefectCase>(DefectCase{in}));
        nests.back()->make(m == 1);
        DefectCase& n = *nests.back();
        fill_nest(n, [&] (const amrex::MultiFab& f, int c, int i, int j, int k) {
            // The second nest's u is larger: the two restrictions of the shared face differ.
            return field(pos(n.geom, f, i, 0), pos(n.geom, f, j, 1), k, c) + ((m == 1 && c == 0) ? 0.3 + 0.1 * k : 0.0);
        });
    }
    fill_parent(p);
    std::vector<gust::NestImage> images;
    for (int m = 0; m < 2; ++m) { images.push_back({&nests[m]->q, &nests[m]->mesh, &nests[m]->geom, kR, offs[m]}); }
    gust::restrict_patches(images, p.mesh, p.geom, ghosts(p.geom), p.q);

    const amrex::IndexType ZF(amrex::IntVect(0, 0, 1));
    amrex::MultiFab om_c(p.q.w.boxArray(), p.dm, 1, 0);
    gust::omega_from_w(p.q.u, p.q.v, p.q.w, p.q.z_phys_nd, p.mesh, p.geom, om_c, 0);
    const amrex::Box dom_c = p.geom.Domain();
    const std::vector<double> oc = gather(om_c, dom_c), ac = gather(p.mesh.az, dom_c);
    const amrex::Box zc = amrex::convert(dom_c, ZF);
    double err = 0, scale = 0;
    amrex::Long bad = 0;
    for (int m = 0; m < 2; ++m) {
        DefectCase& n = *nests[m];
        amrex::MultiFab om_f(n.q.w.boxArray(), n.dm, 1, 0);
        gust::omega_from_w(n.q.u, n.q.v, n.q.w, n.q.z_phys_nd, n.mesh, n.geom, om_f, 0);
        const amrex::Box dom_f = n.geom.Domain();
        const std::vector<double> of = gather(om_f, dom_f), af = gather(n.mesh.az, dom_f);
        const amrex::Box zf = amrex::convert(dom_f, ZF);
        const amrex::Box fp(offs[m], offs[m] + amrex::IntVect(5, 9, dom_c.bigEnd(2) + 1));
        for (amrex::BoxIterator it(fp); it.ok(); ++it) {
            const amrex::IntVect c = it();
            double s = 0;
            for (int b = 0; b < kR; ++b) {
                for (int a = 0; a < kR; ++a) {
                    const amrex::IntVect f(kR * (c[0] - offs[m][0]) + a, kR * (c[1] - offs[m][1]) + b, c[2]);
                    const auto q = static_cast<std::size_t>(zf.index(f));
                    s += af[q] * of[q];
                }
            }
            const auto q = static_cast<std::size_t>(zc.index(c));
            s /= kR * kR * ac[q];
            const double d = std::abs(oc[q] - s);
            if (!std::isfinite(d)) { ++bad; }
            err = std::max(err, d);
            scale = std::max(scale, std::abs(s));
        }
    }
    EXPECT_EQ(bad, 0) << "non-finite Omega on a footprint";
    EXPECT_GT(scale, 1.0e-3) << "the nest Omega is too small to test";
    EXPECT_LT(err, kTol * std::max(1.0, scale)) << "the parent's Omega differs from the restricted nest Omega";
}

// Partly open faces (buildings, Mesh::set_open_faces) on flat ground, with face fractions made up for the test.
// The restriction divides the nest's flux through the open parts, sum_f a_f (1 - beta_f) q_f (1/r or 1/r^2), by
// the parent face's open area a_c (1 - beta_c); a parent face closed in full gets 0. Flat ground: the two levels
// share their layers, so a_f = a_c and Omega = w, and the restricted values below follow by hand.
//   x-faces (u), parent face (I, J), nest faces j = 2 (J - 2) + b with u_f = 2 + 2 b + 0.1 k:
//     I - 3 even, J - 2 even: beta_f = 0, 0.5; beta_c = 0.25:  (2 + 0.1 k + 0.5 (4 + 0.1 k)) / 2 / 0.75 = 8/3 + 0.1 k
//     I - 3 even, J - 2 odd:  beta_f = 1, 1;   beta_c = 1:     0
//     I - 3 odd,  J - 2 even: beta_f = 1, 0;   beta_c = 0:     (4 + 0.1 k) / 2                         = 2 + 0.05 k
//     I - 3 odd,  J - 2 odd:  beta_f = 0, 0;   beta_c = 0:     (6 + 0.2 k) / 2                         = 3 + 0.1 k
//   y-faces (v): the same with x and y swapped.
//   z-faces (w), nest faces (2 (I - 3) + a, 2 (J - 2) + b) with w_f = 1 + a + 2 b + 0.1 k:
//     I + J odd:  beta_f = 0, 0.5, 1, 0 (a + 2 b = 0..3); beta_c = 0.375: (6 + 0.25 k) / 4 / 0.625 = 2.4 + 0.1 k
//     I + J even: beta_f = 1 (all four); beta_c = 1: 0
// With the full areas (the restriction before N-114 item 7) these are 3 + 0.1 k, 3 + 0.1 k, 3 + 0.1 k, 3 + 0.1 k and
// 2.5 + 0.1 k.
TEST(NestDefect, RestrictionWeightsTheOpenFaceAreas)
{
    DefectCase p{parent_inputs(false)};
    p.make();
    DefectCase n{nest_inputs(p.in)};
    n.make(true);
    // The face fractions: beta of the nest's faces by their indices, of the parent's by the parent face's case.
    const auto beta_fine = [] (int d, int i, int j) -> double {
        if (d == 2) {
            const double b[4] = {0.0, 0.5, 1.0, 0.0};
            const int a = i / 2, b2 = j / 2;  // the parent face's offsets in the footprint, I - 3 and J - 2
            return ((a + b2) % 2 == 0) ? b[i % 2 + 2 * (j % 2)] : 1.0;
        }
        const int nrm = (d == 0) ? i : j, tan = (d == 0) ? j : i;
        if (nrm % 2 != 0) { return 0.0; }  // between the parent faces: not used
        const double even[4] = {0.0, 0.5, 1.0, 1.0}, odd[4] = {1.0, 0.0, 0.0, 0.0};
        return ((nrm / 2) % 2 == 0) ? even[tan % 4] : odd[tan % 4];
    };
    const amrex::Box fp = footprint(p);
    const auto beta_coarse = [] (int d, int I, int J) -> double {
        const int a = I - kOffset[0], b = J - kOffset[1];
        const int na = (d == 0) ? 12 : 11, nb = (d == 1) ? 10 : 9;  // the footprint's faces of type d
        if (a < 0 || b < 0 || a > na || b > nb) { return 0.0; }
        if (d == 2) { return ((a + b) % 2 == 0) ? 0.375 : 1.0; }
        const int nrm = (d == 0) ? a : b, tan = (d == 0) ? b : a;
        if (nrm % 2 == 0) { return (tan % 2 == 0) ? 0.25 : 1.0; }
        return 0.0;
    };
    for (DefectCase* c : {&p, &n}) {
        const amrex::MultiFab* area[3] = {&c->mesh.ax, &c->mesh.ay, &c->mesh.az};
        std::array<amrex::MultiFab, 3> beta;
        for (int d = 0; d < 3; ++d) {
            beta[static_cast<std::size_t>(d)].define(area[d]->boxArray(), c->dm, 1, 0);
            fill_all(beta[static_cast<std::size_t>(d)], [&] (int i, int j, int) {
                return (c == &n) ? beta_fine(d, i, j) : beta_coarse(d, i, j);
            });
        }
        amrex::MultiFab beta_cell(c->ba, c->dm, 1, 0);
        beta_cell.setVal(0.0);
        c->mesh.set_open_faces(beta[0], beta[1], beta[2], beta_cell, c->geom);
    }
    fill_parent(p);
    fill_nest(n, [&] (const amrex::MultiFab&, int c, int i, int j, int k) {
        if (c == 0) { return 2.0 + 2.0 * (j % 2) + 0.1 * k; }
        if (c == 1) { return 2.0 + 2.0 * (i % 2) + 0.1 * k; }
        if (c == 2) { return 1.0 + (i % 2) + 2.0 * (j % 2) + 0.1 * k; }
        return 300.0;
    });
    gust::restrict_patch(n.q, n.mesh, n.geom, kR, kOffset, p.mesh, p.geom, ghosts(p.geom), p.q);

    const amrex::Box dom_c = p.geom.Domain(), dom_f = n.geom.Domain();
    const amrex::MultiFab* qc[3] = {&p.q.u, &p.q.v, &p.q.w};
    const amrex::MultiFab* ac[3] = {&p.mesh.ax, &p.mesh.ay, &p.mesh.az};
    const amrex::MultiFab* af[3] = {&n.mesh.ax, &n.mesh.ay, &n.mesh.az};
    for (int d = 0; d < 3; ++d) {
        const amrex::IndexType t(amrex::IntVect::TheDimensionVector(d));
        const amrex::Box bc = amrex::convert(dom_c, t), bf = amrex::convert(dom_f, t);
        const std::vector<double> q = gather(*qc[d], dom_c), a_c = gather(*ac[d], dom_c), a_f = gather(*af[d], dom_f);
        const amrex::Box faces = amrex::convert(fp, t);
        double err = 0, aerr = 0;
        amrex::Long bad = 0, closed = 0, cases[4] = {0, 0, 0, 0};
        for (amrex::BoxIterator it(faces); it.ok(); ++it) {
            const amrex::IntVect c = it();
            const int a = c[0] - kOffset[0], b = c[1] - kOffset[1], k = c[2];
            // The precondition of the values by hand: the nest's areas equal the parent's.
            const amrex::IntVect f(kR * a, kR * b, k);
            aerr = std::max(aerr, std::abs(a_f[static_cast<std::size_t>(bf.index(f))]
                                           - a_c[static_cast<std::size_t>(bc.index(c))]));
            double e = 0;
            if (d == 2) {
                e = ((a + b) % 2 == 0) ? 2.4 + 0.1 * k : 0.0;
                ++cases[(a + b) % 2];
            } else {
                const int nrm = (d == 0) ? a : b, tan = (d == 0) ? b : a;
                const int m = 2 * (nrm % 2) + (tan % 2);
                const double v[4] = {8.0 / 3.0 + 0.1 * k, 0.0, 2.0 + 0.05 * k, 3.0 + 0.1 * k};
                e = v[m];
                ++cases[m];
            }
            const double x = q[static_cast<std::size_t>(bc.index(c))];
            if (!std::isfinite(x)) { ++bad; }
            if (e == 0.0 && x == 0.0) { ++closed; }
            err = std::max(err, std::abs(x - e));
        }
        EXPECT_LT(aerr, kTol) << "face type " << d << ": the nest's areas differ from the parent's on flat ground";
        EXPECT_EQ(bad, 0) << "face type " << d << ": a non-finite restricted velocity";
        EXPECT_GT(closed, 0) << "face type " << d << ": no closed parent face restricted to 0";
        for (int m = 0; m < (d < 2 ? 4 : 2); ++m) {
            EXPECT_GT(cases[m], 0) << "face type " << d << ": case " << m << " not met";
        }
        EXPECT_LT(err, 10.0 * kTol) << "face type " << d << ": the restricted velocity differs from the open-area value";
    }
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

// The change of S between two sets: the second R differs from the first by a field that is largest (7.5) at one
// point of the inner box and larger still (100) outside it, where S is not set: the change is 7.5 per field, the
// size the largest |R_2| on the inner box. Beyond the plan (tracker row N-114, item 6).
TEST(NestDefect, ChangeOfTheSourceBetweenTwoSets)
{
    const amrex::Box dom(amrex::IntVect(0, 0, 0), amrex::IntVect(19, 15, 7));
    const gust::DefectRegion g = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, kOffset, dom);
    DefectCase p{parent_inputs(false)};
    p.make();
    gust::ResidualFields R1, R2;
    R1.define(p.ba, p.dm);
    R2.define(p.ba, p.dm);
    auto rv = [] (int i, int j, int k, int c) { return 1.0 + i - 2.0 * j + 3.0 * k + 0.25 * c; };
    const amrex::IntVect peak(7, 6, 3);   // inside the inner box (5..12, 4..9, all z) at every index type
    const amrex::IntVect out(1, 1, 1);    // outside the footprint
    int c = 0;
    for (auto& f : R1.fields()) {
        const int cc = c++;
        fill_all(*f.mf, [&] (int i, int j, int k) { return rv(i, j, k, cc); });
    }
    c = 0;
    for (auto& f : R2.fields()) {
        const int cc = c++;
        fill_all(*f.mf, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double bump = iv == peak ? 7.5 : iv == out ? 100.0 : 0.01 * ((i + j + k) % 5);
            return rv(i, j, k, cc) + bump;
        });
    }
    gust::NestDefect S(p.ba, p.dm);
    S.set(R1, {g});
    EXPECT_TRUE(S.change().first);
    S.set(R2, {g});
    EXPECT_FALSE(S.change().first);
    const amrex::MultiFab* r2[5] = {&R2.u, &R2.v, &R2.w, &R2.theta, &R2.tke};
    for (int n = 0; n < 5; ++n) {
        const amrex::Box in = amrex::convert(g.inner, r2[n]->ixType());
        ASSERT_TRUE(in.contains(peak));
        double big = 0;
        each_valid(*r2[n], [&] (int i, int j, int k, double v) {
            if (in.contains(amrex::IntVect(i, j, k))) { big = std::max(big, std::abs(v)); }
        });
        amrex::ParallelAllReduce::Max(big, amrex::ParallelDescriptor::Communicator());
        EXPECT_DOUBLE_EQ(static_cast<double>(S.change().change[static_cast<std::size_t>(n)]), 7.5) << "field " << n;
        EXPECT_DOUBLE_EQ(static_cast<double>(S.change().size[static_cast<std::size_t>(n)]), big) << "field " << n;
    }
    // The same R again: no change.
    S.set(R2, {g});
    for (int n = 0; n < 5; ++n) { EXPECT_EQ(static_cast<double>(S.change().change[static_cast<std::size_t>(n)]), 0.0); }
}

TEST(NestDefect, RelativeFeedbackChangeClosedForm)
{
    gust::FeedbackChange c;
    const std::array<amrex::Real, 5> nofloor{};
    EXPECT_TRUE(std::isinf(static_cast<double>(gust::relative_feedback_change(c, nofloor))));   // first
    c.first = false;
    c.change = {1, 2, 0, 3, 0.5};
    c.size = {10, 100, 0, 1, 50};
    // max(1/10, 2/100, 0/0 = 0, 3/1, 0.5/50) = 3
    const double eps = 4.0 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon());   // relative
    EXPECT_NEAR(static_cast<double>(gust::relative_feedback_change(c, nofloor)), 3.0, 3.0 * eps);
    // A floor above the size measures absolutely: theta 3 / 30.
    const std::array<amrex::Real, 5> floor{0, 0, 0, 30, 0};
    EXPECT_NEAR(static_cast<double>(gust::relative_feedback_change(c, floor)), 0.1, 0.1 * eps);
    c.change[2] = 1;   // a change over a size and floor of 0
    EXPECT_TRUE(std::isinf(static_cast<double>(gust::relative_feedback_change(c, floor))));
    c.change[2] = 0;
    {
        const gust::test::ScopedNoInvalidTrap no_trap;   // the NaN is planted on purpose
        c.change[4] = std::numeric_limits<amrex::Real>::quiet_NaN();
        EXPECT_TRUE(std::isinf(static_cast<double>(gust::relative_feedback_change(c, floor))));
    }
}

// The window test on the per-exchange changes (FeedbackWindow): W = 4, tol = 1e-3.
TEST(NestDefect, FeedbackWindowSumAndAitken)
{
    const amrex::Real tol = 1.0e-3;
    gust::FeedbackWindow w(4, tol);
    const amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    w.add(inf);   // the first exchange
    EXPECT_FALSE(w.settled());
    for (int n = 0; n < 3; ++n) { w.add(amrex::Real(1.0e-5)); }
    EXPECT_TRUE(std::isinf(static_cast<double>(w.sum())));   // 3 of 4
    EXPECT_FALSE(w.settled());
    // A steady drift of 1e-4 per exchange: sum 4e-4 <= tol, but rho = 1.
    gust::FeedbackWindow drift(4, tol);
    for (int n = 0; n < 4; ++n) { drift.add(amrex::Real(1.0e-4)); }
    EXPECT_NEAR(static_cast<double>(drift.sum()), 4.0e-4,
                4.0e-4 * 8.0 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()));
    EXPECT_FALSE(drift.settled());
    // Halving each exchange: d1 = 1.5e-4, d2 = 3.75e-5, rho = 0.25, remainder 1.25e-5.
    gust::FeedbackWindow fast(4, tol);
    for (const double d : {1.0e-4, 5.0e-5, 2.5e-5, 1.25e-5}) { fast.add(amrex::Real(d)); }
    EXPECT_TRUE(fast.settled());
    // Sum 9.5e-4 <= tol, but rho = 0.9: remainder 4.5e-4 * 0.9 / 0.1 = 4.05e-3 > tol.
    gust::FeedbackWindow slow(4, tol);
    for (const double d : {2.5e-4, 2.5e-4, 2.25e-4, 2.25e-4}) { slow.add(amrex::Real(d)); }
    EXPECT_FALSE(slow.settled());
    // The same with a 1.2e-3 sum (rho = 0.2) fails on the sum alone.
    gust::FeedbackWindow big(4, tol);
    for (const double d : {5.0e-4, 5.0e-4, 1.0e-4, 1.0e-4}) { big.add(amrex::Real(d)); }
    EXPECT_FALSE(big.settled());
    // The oldest leaves the window: after 1e-2 and four zeros the window is all zeros (round-off: settled).
    gust::FeedbackWindow old(4, tol);
    old.add(amrex::Real(1.0e-2));
    for (int n = 0; n < 4; ++n) { old.add(0); }
    EXPECT_EQ(static_cast<double>(old.sum()), 0.0);
    EXPECT_TRUE(old.settled());
    // An infinite change (a positive change over a zero size) empties the window.
    old.add(inf);
    EXPECT_FALSE(old.settled());
}
