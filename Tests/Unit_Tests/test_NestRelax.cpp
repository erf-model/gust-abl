// The patches' feedback by relaxation over their footprints (NestRelax.H;
// plan rev. 6 tracker, PR-71, nest.feedback = relax, beyond the plan):
//   UT-8.3.1: on a flat parent with uneven boxes and dx != dy, two
//     overlapping footprints: U is the mean over the footprint's cells of
//     the horizontal speed of q* (cell value from its faces); c = f U /
//     min(dx, dy) times w(s) = min(1, s / 2), s the distance in cells from
//     the nearest footprint face, computed here from the points' positions;
//     overlapping regions take the larger rate; 0 outside the footprints and
//     on their faces; R gains c (q - q*) for u, v, w, theta and k on a
//     non-uniform field, exactly 0 elsewhere; differences() is max |q - q*|
//     where c > 0. A footprint that spans a periodic direction has no faces
//     there.
//   UT-8.3.2: on the hill of HillTestHelpers.H as the application sets it up
//     (wall law, wall-k rows, Dirichlet ground and top w faces): an inactive
//     source leaves the residual and the implicit rows bit for bit; an active
//     one changes the residual by c (q - q*) on the footprint's rows that
//     keep their residual, nothing on the rows that replace theirs, and
//     nothing (by value) outside the footprint.
//   UT-8.3.3: implicit_rows gains c on the diagonal of those rows, no other
//     coefficient changes, and the change of the residual along a direction
//     is that diagonal times it, by finite differences.
//   UT-8.3.4: the rates and U do not depend on the rank map (bit for bit),
//     and agree with a one-box layout; with UT-8.3.1's closed form on 1, 2
//     and 4 ranks (the unit_np* tests).
//   nest.feedback_factor's check names the key.
//   The change of the target between two sets (the joint run's stop test,
//     JointRun.H; tracker row N-114, item 6, beyond the plan): max |q*_new -
//     q*_old| where c > 0 only (a larger change on a footprint face, where
//     c = 0, and outside is not counted), the sizes (one velocity scale, 0 for
//     theta, the largest k), against values computed here; first on the
//     first set only.
//   Boxes outside every footprint (tracker row N-114, item 14, beyond the
//     plan): on 4 x 4-cell boxes with two small footprints, the source and
//     the diagonal reach every point where the closed-form rate is > 0 (the
//     boxes B's margin alone reaches among them) and leave every other point
//     bit for bit; 12 of the 20 boxes hold no footprint point. The test
//     shows no box with a footprint point is skipped; the saving itself is
//     not visible in values.
#include "HillTestHelpers.H"
#include "ImplicitOperator.H"
#include "Inputs.H"
#include "MassFlux.H"
#include "Mesh.H"
#include "NestDefect.H"
#include "NestRelax.H"
#include "Residual.H"
#include "State.H"
#include "SurfaceLayer.H"
#include "TestHelpers.H"
#include "TransportTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>

#include <gtest/gtest.h>

#include <algorithm>
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
constexpr double kTol = kDouble ? 1.0e-12 : 2.0e-5;
constexpr double kMargin = 2;   // kDefectMargin, written out so a change of the constant shows here

// --- UT-8.3.1 and 8.3.4: a flat parent -----------------------------------------------------------

struct Parent
{
    gust::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    gust::State q;

    explicit Parent (bool periodic_y, int mgs_x = 7, int mgs_y = 6, bool reverse = false)
    {
        // 20 x 16 x 8 cells of 100 m x 80 m (dx != dy), uneven boxes.
        in.mesh.n_cell = amrex::IntVect(20, 16, 8);
        in.mesh.prob_lo = {-1000, 500, 0};
        in.mesh.prob_hi = {1000, 1780, 800};
        in.mesh.is_periodic = {0, periodic_y ? 1 : 0, 0};
        in.mesh.max_grid_size = {mgs_x, mgs_y};
        in.mesh.dz0 = 40;
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
    }
};

double pos (const amrex::Geometry& g, amrex::IndexType t, int i, int d)
{
    const double off = t.nodeCentered(d) ? 0.0 : 0.5;
    return static_cast<double>(g.ProbLo(d)) + (i + off) * static_cast<double>(g.CellSize(d));
}

// A smooth non-uniform field, different for each component c (and for the state, s = 0, and the target, s = 1).
double field (double x, double y, int k, int c, int s)
{
    return (4.0 + 0.5 * c - 0.7 * s) + 0.2 * k + 0.003 * x * (1 + c) - 0.002 * y * (1 + s)
         + 1.0e-6 * (1 + 0.1 * k) * x * y + 0.4 * std::sin(0.004 * x + 0.3 * c + s) * std::cos(0.003 * y);
}

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
        amrex::FArrayBox h(b, 1, amrex::The_Pinned_Arena());
        h.copy<amrex::RunOn::Device>(mf[mfi], b);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { fn(i, j, k, static_cast<double>(h(amrex::IntVect(i, j, k)))); });
    }
}

// u, v, w, theta, k of q filled with field(., s).
void fill_state (const Parent& p, gust::State& q, int s)
{
    amrex::MultiFab* f[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    for (int c = 0; c < 5; ++c) {
        const amrex::IndexType t = f[c]->ixType();
        fill_all(*f[c], [&] (int i, int j, int k) { return field(pos(p.geom, t, i, 0), pos(p.geom, t, j, 1), k, c, s); });
    }
}

// U of a footprint from the field (s = 1) directly: the mean over its cells of |(u_c, v_c)|.
double speed_of (const Parent& p, const amrex::Box& fp)
{
    const amrex::IndexType xf(amrex::IntVect(1, 0, 0)), yf(amrex::IntVect(0, 1, 0));
    double sum = 0;
    for (amrex::BoxIterator it(fp); it.ok(); ++it) {
        const int i = it()[0], j = it()[1], k = it()[2];
        const double yc = pos(p.geom, xf, j, 1), xc = pos(p.geom, yf, i, 0);
        const double uc = 0.5 * (field(pos(p.geom, xf, i, 0), yc, k, 0, 1) + field(pos(p.geom, xf, i + 1, 0), yc, k, 0, 1));
        const double vc = 0.5 * (field(xc, pos(p.geom, yf, j, 1), k, 1, 1) + field(xc, pos(p.geom, yf, j + 1, 1), k, 1, 1));
        sum += std::sqrt(uc * uc + vc * vc);
    }
    return sum / static_cast<double>(fp.numPts());
}

// The rate at point (i, j) of index type t from the regions, from the points' positions: w = min(1, s / 2),
// s = the distance from the nearest footprint face over dx or dy (none across a direction the region spans).
double expect_rate (const Parent& p, const std::vector<gust::DefectRegion>& regions, const std::vector<double>& c0,
                    amrex::IndexType t, int i, int j)
{
    double c = 0;
    for (std::size_t r = 0; r < regions.size(); ++r) {
        const amrex::Box& fp = regions[r].footprint;
        double s = std::numeric_limits<double>::infinity();
        const int ij[2] = {i, j};
        bool inside = true;
        for (int d = 0; d < 2; ++d) {
            const double h = static_cast<double>(p.geom.CellSize(d));
            const double lo = static_cast<double>(p.geom.ProbLo(d)) + fp.smallEnd(d) * h;
            const double hi = static_cast<double>(p.geom.ProbLo(d)) + (fp.bigEnd(d) + 1) * h;
            const double x = pos(p.geom, t, ij[d], d);
            if (x < lo - 1e-9 * h || x > hi + 1e-9 * h) { inside = false; }
            if (regions[r].inner.length(d) != fp.length(d)) { s = std::min(s, std::min(x - lo, hi - x) / h); }
        }
        if (!inside) { continue; }
        c = std::max(c, c0[r] * std::min(1.0, std::max(0.0, s / kMargin)));
    }
    return c;
}

struct Tally
{
    double err{0};      // max |value - expect| / scale
    long on{0};         // points with expect > 0
    long edge{0};       // points with 0 < w < 1 (the margin)
    long off_nonzero{0};  // points with expect == 0 and value != 0
};

void reduce (Tally& t)
{
    auto comm = amrex::ParallelDescriptor::Communicator();
    amrex::ParallelAllReduce::Max(t.err, comm);
    amrex::ParallelAllReduce::Sum(t.on, comm);
    amrex::ParallelAllReduce::Sum(t.edge, comm);
    amrex::ParallelAllReduce::Sum(t.off_nonzero, comm);
}

constexpr double kFactor = 3;

// UT-8.3.1's checks of a set of regions on parent p.
void check_closed_form (Parent& p, const std::vector<gust::DefectRegion>& regions)
{
    gust::State qstar;
    qstar.define(p.ba, p.dm);
    fill_state(p, qstar, 1);
    fill_state(p, p.q, 0);
    gust::NestRelax S(p.ba, p.dm, p.geom, amrex::Real(kFactor));
    EXPECT_FALSE(S.active());
    S.set(qstar, regions);
    EXPECT_TRUE(S.active());
    const double dxh = std::min(static_cast<double>(p.geom.CellSize(0)), static_cast<double>(p.geom.CellSize(1)));
    EXPECT_EQ(dxh, 80.0);
    std::vector<double> c0;
    ASSERT_EQ(S.nests().size(), regions.size());
    for (std::size_t r = 0; r < regions.size(); ++r) {
        const double U = speed_of(p, regions[r].footprint);
        EXPECT_EQ(S.nests()[r].name, regions[r].name);
        EXPECT_NEAR(static_cast<double>(S.nests()[r].speed), U, kTol * U) << regions[r].name;
        c0.push_back(kFactor * U / dxh);
        EXPECT_NEAR(static_cast<double>(S.nests()[r].rate), c0.back(), kTol * c0.back()) << regions[r].name;
    }
    // The rates.
    const Loc locs[4] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell};
    const char* lname[4] = {"xface", "yface", "zface", "cell"};
    for (int n = 0; n < 4; ++n) {
        const amrex::MultiFab& c = S.rate(locs[n]);
        const amrex::IndexType t = c.ixType();
        Tally tl;
        each_valid(c, [&] (int i, int j, int, double v) {
            const double e = expect_rate(p, regions, c0, t, i, j);
            if (e == 0) { if (v != 0) { ++tl.off_nonzero; } return; }
            ++tl.on;
            const double w = e / *std::max_element(c0.begin(), c0.end());
            if (w < 1 - 1e-9) { ++tl.edge; }
            tl.err = std::max(tl.err, std::abs(v - e) / e);
        });
        reduce(tl);
        EXPECT_LT(tl.err, kTol) << lname[n] << ": the rate is not f U / dx w(s)";
        EXPECT_EQ(tl.off_nonzero, 0) << lname[n] << ": a rate outside the footprints or on their faces";
        EXPECT_GT(tl.on, 0) << lname[n];
        EXPECT_GT(tl.edge, 0) << lname[n] << ": no point in the margin: the test is vacuous";
    }
    // The source.
    gust::ResidualFields R;
    R.define(p.ba, p.dm);
    for (auto& f : R.fields()) { f.mf->setVal(0.0); }
    S.add_to_residual(p.q, R);
    const amrex::MultiFab* rf[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    const std::vector<amrex::Real> diff = S.differences(p.q);
    ASSERT_EQ(diff.size(), 10u);
    for (int n = 0; n < 5; ++n) {
        const amrex::IndexType t = rf[n]->ixType();
        Tally tl;
        double dmax = 0;
        each_valid(*rf[n], [&] (int i, int j, int k, double v) {
            const double e = expect_rate(p, regions, c0, t, i, j);
            const double x = pos(p.geom, t, i, 0), y = pos(p.geom, t, j, 1);
            const double a = field(x, y, k, n, 0), b = field(x, y, k, n, 1);
            if (e == 0) { if (v != 0) { ++tl.off_nonzero; } return; }
            ++tl.on;
            const double expect = e * (a - b);
            tl.err = std::max(tl.err, std::abs(v - expect) / (e * (std::abs(a) + std::abs(b))));
            dmax = std::max(dmax, std::abs(a - b));
        });
        reduce(tl);
        amrex::ParallelAllReduce::Max(dmax, amrex::ParallelDescriptor::Communicator());
        EXPECT_LT(tl.err, kTol) << names[n] << ": R is not c (q - q*)";
        EXPECT_EQ(tl.off_nonzero, 0) << names[n] << ": a source outside the footprints";
        EXPECT_GT(tl.on, 0) << names[n];
        EXPECT_GT(dmax, 0.1) << names[n] << ": q - q* is too small for the test";
        EXPECT_NEAR(static_cast<double>(diff[static_cast<std::size_t>(n)]), dmax, kTol * dmax * 10) << names[n];
        EXPECT_GT(static_cast<double>(diff[static_cast<std::size_t>(5 + n)]), 0.0) << names[n];
        EXPECT_LE(static_cast<double>(diff[static_cast<std::size_t>(5 + n)]), dmax * (1 + kTol)) << names[n];
    }
}

// --- UT-8.3.2 and 8.3.3: the hill as the application sets it up ------------------------------------

using gust::hilltest::HillRun;

amrex::MultiFab& hfield (gust::State& q, int f)
{
    amrex::MultiFab* p[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    return *p[f];
}
const amrex::MultiFab& rfield (const gust::ResidualFields& R, int f)
{
    const amrex::MultiFab* p[6] = {&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke};
    return *p[f];
}

// q = q0 + s dq, dq smooth, non-uniform and nonzero everywhere; ghosts refilled.
void perturb (HillRun& h, const gust::State& q0, gust::State& q, double s)
{
    q.define(h.c.ba, h.c.dm);
    auto t = q.fields();
    const amrex::MultiFab* f0[7] = {&q0.u, &q0.v, &q0.w, &q0.p, &q0.theta, &q0.tke, &q0.z_phys_nd};
    for (std::size_t n = 0; n < t.size(); ++n) {
        amrex::MultiFab::Copy(*t[n].mf, *f0[n], 0, 0, 1, t[n].mf->nGrowVect());
    }
    const double L = HillRun::Lx, H = HillRun::H, pi = HillRun::kPi;
    const Loc locs[6] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell, Loc::cell};
    const double amp[6] = {0.3, 0.2, 0.05, 0.4, 0.5, 0.02};
    for (int f = 0; f < 6; ++f) {
        amrex::MultiFab d(hfield(q, f).boxArray(), hfield(q, f).DistributionMap(), 1, hfield(q, f).nGrow());
        const double a = amp[f] * s;
        gust::trtest::fill_loc(d, locs[f], h.c, [=] (double x, double, double z) {
            return a * (1.2 + std::sin(2 * pi * x / L + f) * std::cos(pi * z / H));
        });
        amrex::MultiFab::Add(hfield(q, f), d, 0, 0, 1, hfield(q, f).nGrow());
    }
    h.ls->bc.fill_ghosts(q, h.ls->surface(), h.s.fopt.buoyancy);
}

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

// The hill (16 x 2 x 12 cells, four boxes in x, periodic in x and y) and a footprint over cells 3 to 12 in x
// spanning y (no y faces).
HillRun make_hill () { return HillRun(0.003, 16, 12, 4); }
gust::DefectRegion hill_region ()
{
    return gust::DefectRegion{"nest.a", amrex::Box(amrex::IntVect(3, 0, 0), amrex::IntVect(12, 1, 11)),
                              amrex::Box(amrex::IntVect(5, 0, 0), amrex::IntVect(10, 1, 11))};
}

long footprint_columns () { return 10 * 2; }

} // namespace

TEST(NestRelax, SourceIsItsClosedFormOnUnevenBoxes)   // UT-8.3.1
{
    ASSERT_EQ(gust::kDefectMargin, 2);
    Parent p(false);
    const amrex::Box dom = p.geom.Domain();
    // A: parent cells 3-14 x 2-11; B: 12-18 x 8-14, overlapping A.
    const gust::DefectRegion a = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, amrex::IntVect(3, 2, 0), dom);
    const gust::DefectRegion b = gust::defect_region("nest.b", amrex::IntVect(14, 14, 8), 2, amrex::IntVect(12, 8, 0), dom);
    check_closed_form(p, {a, b});
}

TEST(NestRelax, PeriodicSpanHasNoFaces)   // UT-8.3.1
{
    Parent p(true);
    const amrex::Box dom = p.geom.Domain();
    const gust::DefectRegion a = gust::defect_region("nest.a", amrex::IntVect(24, 32, 8), 2, amrex::IntVect(3, 0, 0), dom,
                                                     {false, true});
    ASSERT_EQ(a.inner.length(1), a.footprint.length(1));
    check_closed_form(p, {a});
}

TEST(NestRelax, FootprintOnAPeriodicSeam)   // UT-8.3.1, beyond the plan
{
    // A nest from parent face 0 in a periodic y that does not span it: it has y faces, and the seam face (0 and its
    // periodic image 16) gets 0 on both copies.
    Parent p(true);
    const gust::DefectRegion a = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, amrex::IntVect(3, 0, 0),
                                                     p.geom.Domain(), {false, false});
    ASSERT_EQ(a.footprint.smallEnd(1), 0);
    ASSERT_LT(a.footprint.length(1), p.geom.Domain().length(1));
    check_closed_form(p, {a});
    gust::State qstar;
    qstar.define(p.ba, p.dm);
    fill_state(p, qstar, 1);
    gust::NestRelax S(p.ba, p.dm, p.geom, amrex::Real(kFactor));
    S.set(qstar, {a});
    double seam = 0;
    long seen = 0;
    each_valid(S.rate(Loc::yface), [&] (int, int j, int, double v) {
        if (j == 0 || j == 16) { seam = std::max(seam, std::abs(v)); ++seen; }
    });
    amrex::ParallelAllReduce::Max(seam, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(seen, amrex::ParallelDescriptor::Communicator());
    EXPECT_GT(seen, 0);
    EXPECT_EQ(seam, 0.0) << "a rate on the periodic seam face";
}

TEST(NestRelax, InactiveLeavesTheHillBitForBit)   // UT-8.3.2
{
    HillRun h = make_hill();
    gust::NestRelax S(h.c.ba, h.c.dm, h.c.geom, amrex::Real(kFactor));
    gust::State q;
    perturb(h, h.q, q, 1.0);
    gust::ResidualInputs off = h.ls->res_in;
    ASSERT_EQ(off.relax, nullptr);
    gust::ResidualInputs on = off;
    on.relax = &S;
    gust::ResidualFields Ron, Roff;
    gust::residual(q, on, h.c.mesh, h.c.geom, Ron);
    gust::residual(q, off, h.c.mesh, h.c.geom, Roff);
    const auto a = Ron.fields();
    const auto b = Roff.fields();
    for (std::size_t n = 0; n < a.size(); ++n) {
        EXPECT_EQ(gust::test::count_bit_differences(*a[n].mf, *b[n].mf, 0), 0) << a[n].name;
    }
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    for (const Equation eq : eqs) {
        amrex::MultiFab Aon, Aoff;
        rows(h, q, on, eq, Aon);
        rows(h, q, off, eq, Aoff);
        EXPECT_EQ(gust::test::count_bit_differences(Aon, Aoff, 0), 0) << static_cast<int>(eq);
    }
}

TEST(NestRelax, ResidualGainsTheSourceOnTheFootprintRows)   // UT-8.3.2
{
    HillRun h = make_hill();
    gust::NestRelax S(h.c.ba, h.c.dm, h.c.geom, amrex::Real(kFactor));
    S.set(h.q, {hill_region()});   // q* = the unperturbed state
    gust::State q;
    perturb(h, h.q, q, 1.0);
    gust::ResidualInputs off = h.ls->res_in;
    gust::ResidualInputs on = off;
    on.relax = &S;
    gust::ResidualFields Ron, Roff;
    gust::residual(q, on, h.c.mesh, h.c.geom, Ron);
    gust::residual(q, off, h.c.mesh, h.c.geom, Roff);
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    const Loc locs[5] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell};
    const int fi[5] = {0, 1, 2, 4, 5};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    // The rows that replace their residual: the ground and top w faces, the first cells' wall-k rows.
    const long expect_replaced[5] = {0, 0, 2 * footprint_columns(), 0, footprint_columns()};
    const double tol = kDouble ? 1e-11 : 1e-4;
    for (int n = 0; n < 5; ++n) {
        const amrex::MultiFab& c = S.rate(locs[n]);
        const amrex::MultiFab& t = S.target(eqs[n]);
        long outside = 0, replaced = 0, kept = 0;
        double err = 0;
        for (amrex::MFIter mfi(c, false); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox();
            auto host = [&] (const amrex::MultiFab& mf) {
                amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
                f.copy<amrex::RunOn::Device>(mf[mfi], b, 0, b, 0, 1);
                return f;
            };
            const auto hc = host(c), ht = host(t), ha = host(rfield(Ron, fi[n])), hb = host(rfield(Roff, fi[n]));
            const auto hq = host(hfield(q, fi[n]));
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                const double dR = static_cast<double>(ha(iv)) - static_cast<double>(hb(iv));
                const double cc = static_cast<double>(hc(iv));
                if (cc == 0) { if (dR != 0) { ++outside; } return; }
                if (dR == 0) { ++replaced; return; }
                ++kept;
                const double e = cc * (static_cast<double>(hq(iv)) - static_cast<double>(ht(iv)));
                const double scale = std::abs(static_cast<double>(ha(iv))) + std::abs(static_cast<double>(hb(iv))) + std::abs(e);
                err = std::max(err, std::abs(dR - e) / scale);
            });
        }
        auto comm = amrex::ParallelDescriptor::Communicator();
        amrex::ParallelAllReduce::Sum(outside, comm);
        amrex::ParallelAllReduce::Sum(replaced, comm);
        amrex::ParallelAllReduce::Sum(kept, comm);
        amrex::ParallelAllReduce::Max(err, comm);
        EXPECT_EQ(outside, 0) << names[n] << ": the residual changed outside the footprint";
        EXPECT_EQ(replaced, expect_replaced[n]) << names[n] << ": rows without the source";
        EXPECT_GT(kept, 0) << names[n];
        EXPECT_LT(err, tol) << names[n] << ": R does not gain c (q - q*)";
    }
    // Continuity has no source.
    EXPECT_EQ(gust::test::count_bit_differences(Ron.p, Roff.p, 0), 0);
}

TEST(NestRelax, DiagonalIsTheDerivativeOfTheSource)   // UT-8.3.3
{
    HillRun h = make_hill();
    gust::NestRelax S(h.c.ba, h.c.dm, h.c.geom, amrex::Real(kFactor));
    S.set(h.q, {hill_region()});
    gust::State q, q2;
    perturb(h, h.q, q, 1.0);
    const double eps = kDouble ? 1e-3 : 1e-1;   // linear source: the step only sets the round-off
    perturb(h, q, q2, eps);
    gust::ResidualInputs off = h.ls->res_in;
    gust::ResidualInputs on = off;
    on.relax = &S;
    gust::ResidualFields a1, b1, a2, b2;
    gust::residual(q, on, h.c.mesh, h.c.geom, a1);
    gust::residual(q, off, h.c.mesh, h.c.geom, b1);
    gust::residual(q2, on, h.c.mesh, h.c.geom, a2);
    gust::residual(q2, off, h.c.mesh, h.c.geom, b2);
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    const Loc locs[5] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell};
    const int fi[5] = {0, 1, 2, 4, 5};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    const double tol_fd = kDouble ? 1e-6 : 1e-2;
    for (int n = 0; n < 5; ++n) {
        amrex::MultiFab Aon, Aoff;
        rows(h, q, on, eqs[n], Aon);
        rows(h, q, off, eqs[n], Aoff);
        const amrex::MultiFab& c = S.rate(locs[n]);
        long others = 0, plus = 0, identity = 0, wrong = 0;
        double fd_err = 0;
        for (amrex::MFIter mfi(Aon, false); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox();
            auto host = [&] (const amrex::MultiFab& mf, int nc) {
                amrex::FArrayBox f(b, nc, amrex::The_Pinned_Arena());
                f.copy<amrex::RunOn::Device>(mf[mfi], b, 0, b, 0, nc);
                return f;
            };
            const auto ha = host(Aon, gust::n_implicit_coefs), hb = host(Aoff, gust::n_implicit_coefs);
            const auto hc = host(c, 1);
            const auto r1a = host(rfield(a1, fi[n]), 1), r1b = host(rfield(b1, fi[n]), 1);
            const auto r2a = host(rfield(a2, fi[n]), 1), r2b = host(rfield(b2, fi[n]), 1);
            const auto x1 = host(hfield(q, fi[n]), 1), x2 = host(hfield(q2, fi[n]), 1);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                for (int cc = 1; cc < gust::n_implicit_coefs; ++cc) {
                    if (ha(iv, cc) != hb(iv, cc)) { ++others; }
                }
                const double cv = static_cast<double>(hc(iv));
                const double dP = static_cast<double>(ha(iv, gust::coef_P)) - static_cast<double>(hb(iv, gust::coef_P));
                if (cv == 0) { if (dP != 0) { ++wrong; } return; }
                if (dP == 0) { ++identity; return; }   // a row that replaces its own (identity or wall-k row)
                if (std::abs(dP - cv) > (kDouble ? 1e-12 : 1e-4) * (std::abs(static_cast<double>(ha(iv, gust::coef_P))) + cv)) {
                    ++wrong;
                    return;
                }
                ++plus;
                const double ddR = (static_cast<double>(r2a(iv)) - static_cast<double>(r2b(iv)))
                                   - (static_cast<double>(r1a(iv)) - static_cast<double>(r1b(iv)));
                const double dx = static_cast<double>(x2(iv)) - static_cast<double>(x1(iv));
                if (dx == 0) { ++wrong; return; }
                fd_err = std::max(fd_err, std::abs(ddR - dP * dx) / (std::abs(dx) * cv));
            });
        }
        auto comm = amrex::ParallelDescriptor::Communicator();
        amrex::ParallelAllReduce::Sum(others, comm);
        amrex::ParallelAllReduce::Sum(plus, comm);
        amrex::ParallelAllReduce::Sum(identity, comm);
        amrex::ParallelAllReduce::Sum(wrong, comm);
        amrex::ParallelAllReduce::Max(fd_err, comm);
        const long expect_identity = (n == 2) ? 2 * footprint_columns() : (n == 4) ? footprint_columns() : 0;
        EXPECT_EQ(others, 0) << names[n] << ": off-diagonal coefficients changed";
        EXPECT_EQ(wrong, 0) << names[n] << ": a_P changed by neither 0 nor c";
        EXPECT_EQ(identity, expect_identity) << names[n] << ": rows that replace theirs";
        EXPECT_GT(plus, 0) << names[n];
        EXPECT_LT(fd_err, tol_fd) << names[n] << ": the diagonal is not the derivative of the source";
    }
}

TEST(NestRelax, SameSourceOnAnyRankMapAndLayout)   // UT-8.3.4
{
    std::vector<std::vector<double>> rates;
    std::vector<double> speeds;
    for (int layout = 0; layout < 3; ++layout) {
        // 0: uneven boxes; 1: the same boxes on reversed ranks; 2: one box.
        Parent p(false, layout == 2 ? 64 : 7, layout == 2 ? 64 : 6, layout == 1);
        const gust::DefectRegion a = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, amrex::IntVect(3, 2, 0),
                                                         p.geom.Domain());
        gust::State qstar;
        qstar.define(p.ba, p.dm);
        fill_state(p, qstar, 1);
        gust::NestRelax S(p.ba, p.dm, p.geom, amrex::Real(kFactor));
        S.set(qstar, {a});
        speeds.push_back(static_cast<double>(S.nests()[0].speed));
        std::vector<double> all;
        for (const Loc loc : {Loc::xface, Loc::yface, Loc::zface, Loc::cell}) {
            const amrex::MultiFab& c = S.rate(loc);
            const amrex::Box b = amrex::convert(p.geom.Domain(), c.ixType());
            std::vector<double> g(static_cast<std::size_t>(b.numPts()), 0.0), cnt(g.size(), 0.0);
            each_valid(c, [&] (int i, int j, int k, double v) {
                const auto m = static_cast<std::size_t>(b.index(amrex::IntVect(i, j, k)));
                g[m] += v;   // a face shared by two boxes is a valid point of both: averaged below
                cnt[m] += 1;
            });
            amrex::ParallelDescriptor::ReduceRealSum(g.data(), static_cast<int>(g.size()));
            amrex::ParallelDescriptor::ReduceRealSum(cnt.data(), static_cast<int>(cnt.size()));
            for (std::size_t m = 0; m < g.size(); ++m) { if (cnt[m] > 0) { g[m] /= cnt[m]; } }
            all.insert(all.end(), g.begin(), g.end());
        }
        rates.push_back(all);
    }
    EXPECT_EQ(speeds[0], speeds[1]) << "U depends on the rank map";
    EXPECT_NEAR(speeds[2], speeds[0], kTol * speeds[0]) << "U depends on the box layout";
    long differ = 0;
    double rel = 0;
    for (std::size_t m = 0; m < rates[0].size(); ++m) {
        if (rates[0][m] != rates[1][m]) { ++differ; }
        rel = std::max(rel, std::abs(rates[2][m] - rates[0][m]) / (std::abs(rates[0][m]) + 1e-30));
    }
    EXPECT_EQ(differ, 0) << "the rates depend on the rank map";
    EXPECT_LT(rel, kTol) << "the rates depend on the box layout";
}

TEST(NestRelax, FactorCheckNamesTheKey)
{
    EXPECT_TRUE(gust::check_feedback_factor(amrex::Real(1)).empty());
    EXPECT_TRUE(gust::check_feedback_factor(amrex::Real(16)).empty());
    const amrex::Real bad[] = {0, amrex::Real(-1), std::numeric_limits<amrex::Real>::infinity(),
                               std::numeric_limits<amrex::Real>::quiet_NaN()};
    for (const amrex::Real f : bad) {
        const std::string err = gust::check_feedback_factor(f);
        EXPECT_NE(err.find("nest.feedback_factor"), std::string::npos) << f;
    }
}

TEST(NestRelax, ChangeOfTheTargetBetweenTwoSets)
{
    Parent p(false);
    const gust::DefectRegion a = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, amrex::IntVect(3, 2, 0),
                                                     p.geom.Domain());
    gust::State q1, q2;
    q1.define(p.ba, p.dm);
    q2.define(p.ba, p.dm);
    fill_state(p, q1, 1);
    // q2 = q1 + a bump: 0.5 at one point deep in the footprint (c > 0 there at every index type), 50 on the
    // footprint's low x face (x-face 3: c = 0) and outside the footprint, 0.001 (i + j + k) mod 3 elsewhere.
    const amrex::IntVect deep(8, 6, 4), face(3, 6, 4), outside(17, 14, 2);
    amrex::MultiFab* f2[5] = {&q2.u, &q2.v, &q2.w, &q2.theta, &q2.tke};
    for (int c = 0; c < 5; ++c) {
        const amrex::IndexType t = f2[c]->ixType();
        fill_all(*f2[c], [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double bump = iv == deep ? 0.5 : (iv == outside || (c == 0 && iv == face)) ? 50.0 : 0.001 * ((i + j + k) % 3);
            return field(pos(p.geom, t, i, 0), pos(p.geom, t, j, 1), k, c, 1) + bump;
        });
    }
    gust::NestRelax S(p.ba, p.dm, p.geom, amrex::Real(kFactor));
    S.set(q1, {a});
    EXPECT_TRUE(S.change().first);
    S.set(q2, {a});
    EXPECT_FALSE(S.change().first);
    // The sizes from q2 where c > 0.
    const Loc locs[5] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell};
    double big[5] = {0, 0, 0, 0, 0};
    for (int c = 0; c < 5; ++c) {
        amrex::MultiFab both(f2[c]->boxArray(), f2[c]->DistributionMap(), 2, 0);
        amrex::MultiFab::Copy(both, S.rate(locs[c]), 0, 0, 1, 0);
        amrex::MultiFab::Copy(both, *f2[c], 0, 1, 1, 0);
        for (amrex::MFIter mfi(both, false); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox();
            amrex::FArrayBox h(b, 2, amrex::The_Pinned_Arena());
            h.copy<amrex::RunOn::Device>(both[mfi], b, 0, b, 0, 2);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                if (h(iv, 0) > 0) { big[c] = std::max(big[c], std::abs(static_cast<double>(h(iv, 1)))); }
            });
        }
        amrex::ParallelAllReduce::Max(big[c], amrex::ParallelDescriptor::Communicator());
    }
    const double U = std::max(big[0], std::max(big[1], big[2]));
    const double want[5] = {U, U, U, 0.0, big[4]};
    for (int c = 0; c < 5; ++c) {
        const auto n = static_cast<std::size_t>(c);
        EXPECT_NEAR(static_cast<double>(S.change().change[n]), 0.5, kTol) << "field " << c;
        EXPECT_DOUBLE_EQ(static_cast<double>(S.change().size[n]), want[c]) << "field " << c;
    }
    EXPECT_GT(U, 1.0);
}

TEST(NestRelax, BoxesOutsideEveryFootprintAreSkippedExactly)   // tracker row N-114, item 14, beyond the plan
{
    // 4 x 4-cell boxes (5 x 4 of them); A over cells 1-6 x 1-6, B over 13-17 x 9-13. add_to_residual and
    // add_to_implicit skip a box with no point of a footprint: every point where the closed-form rate is > 0 still
    // gains it (in every box that touches a footprint, among them the boxes B's margin alone reaches), and every
    // other point keeps its value bit for bit.
    Parent p(false, 4, 4);
    const amrex::Box dom = p.geom.Domain();
    const gust::DefectRegion a = gust::defect_region("nest.a", amrex::IntVect(12, 12, 8), 2, amrex::IntVect(1, 1, 0), dom);
    const gust::DefectRegion b = gust::defect_region("nest.b", amrex::IntVect(10, 10, 8), 2, amrex::IntVect(13, 9, 0), dom);
    const std::vector<gust::DefectRegion> regions{a, b};
    // The cell boxes: those that touch a footprint, those that touch B's footprint but not its inner box (the
    // margin), and the others (skipped).
    long touch = 0, margin_only = 0, skipped = 0;
    for (int m = 0; m < p.ba.size(); ++m) {
        const amrex::Box& bx = p.ba[m];
        const bool t = bx.intersects(a.footprint) || bx.intersects(b.footprint);
        touch += t ? 1 : 0;
        skipped += t ? 0 : 1;
        margin_only += (bx.intersects(b.footprint) && !bx.intersects(b.inner)) ? 1 : 0;
    }
    ASSERT_EQ(touch, 8);
    ASSERT_EQ(skipped, 12);
    ASSERT_EQ(margin_only, 3);

    gust::State qstar;
    qstar.define(p.ba, p.dm);
    fill_state(p, qstar, 1);
    fill_state(p, p.q, 0);
    gust::NestRelax S(p.ba, p.dm, p.geom, amrex::Real(kFactor));
    S.set(qstar, regions);
    const std::vector<double> c0{static_cast<double>(S.nests()[0].rate), static_cast<double>(S.nests()[1].rate)};
    constexpr double r0 = 1.5, a0 = 2.0;
    gust::ResidualFields R;
    R.define(p.ba, p.dm);
    for (auto& f : R.fields()) { f.mf->setVal(r0); }
    S.add_to_residual(p.q, R);
    const amrex::MultiFab* rf[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    for (int n = 0; n < 5; ++n) {
        const amrex::IndexType t = rf[n]->ixType();
        amrex::MultiFab A(rf[n]->boxArray(), rf[n]->DistributionMap(), gust::n_implicit_coefs, 0);
        A.setVal(a0);
        S.add_to_implicit(eqs[n], A);
        long on = 0, off_changed = 0, on_missed = 0, other_coef = 0;
        double err_r = 0, err_a = 0;
        for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.validbox();
            amrex::FArrayBox hr(bx, 1, amrex::The_Pinned_Arena()), ha(bx, gust::n_implicit_coefs, amrex::The_Pinned_Arena());
            hr.copy<amrex::RunOn::Device>((*rf[n])[mfi], bx);
            ha.copy<amrex::RunOn::Device>(A[mfi], bx, 0, bx, 0, gust::n_implicit_coefs);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                const double e = expect_rate(p, regions, c0, t, i, j);
                const double vr = static_cast<double>(hr(iv)), va = static_cast<double>(ha(iv, gust::coef_P));
                for (int c = 1; c < gust::n_implicit_coefs; ++c) { if (ha(iv, c) != amrex::Real(a0)) { ++other_coef; } }
                if (e == 0) {
                    if (hr(iv) != amrex::Real(r0) || ha(iv, gust::coef_P) != amrex::Real(a0)) { ++off_changed; }
                    return;
                }
                ++on;
                if (hr(iv) == amrex::Real(r0) || ha(iv, gust::coef_P) == amrex::Real(a0)) { ++on_missed; }
                const double x = pos(p.geom, t, i, 0), y = pos(p.geom, t, j, 1);
                const double qa = field(x, y, k, n, 0), qb = field(x, y, k, n, 1);
                err_r = std::max(err_r, std::abs(vr - r0 - e * (qa - qb)) / (r0 + e * (std::abs(qa) + std::abs(qb))));
                err_a = std::max(err_a, std::abs(va - a0 - e) / (a0 + e));
            });
        }
        auto comm = amrex::ParallelDescriptor::Communicator();
        amrex::ParallelAllReduce::Sum(on, comm);
        amrex::ParallelAllReduce::Sum(off_changed, comm);
        amrex::ParallelAllReduce::Sum(on_missed, comm);
        amrex::ParallelAllReduce::Sum(other_coef, comm);
        amrex::ParallelAllReduce::Max(err_r, comm);
        amrex::ParallelAllReduce::Max(err_a, comm);
        EXPECT_GT(on, 0) << names[n];
        EXPECT_EQ(on_missed, 0) << names[n] << ": a point with c > 0 lost its source (its box was skipped)";
        EXPECT_EQ(off_changed, 0) << names[n] << ": a point with c = 0 changed";
        EXPECT_EQ(other_coef, 0) << names[n] << ": an off-diagonal coefficient changed";
        EXPECT_LT(err_r, kTol) << names[n] << ": R is not R0 + c (q - q*)";
        EXPECT_LT(err_a, kTol) << names[n] << ": A_P is not A0 + c";
    }
}
