// The pressure-consistent relaxation (NestRelax.H, nest.feedback = relax_consistent):
//   on a flat parent with uneven boxes and dx != dy, two overlapping footprints: the added term is
//     w(s) R* with w = min(1, s / 2) from the points' positions (the larger weight where the footprints overlap),
//     0 on the footprints' faces and outside them whatever R* holds there; the rates are those of relax bit for bit.
//   on the hill of HillTestHelpers.H as the application sets it up, with q* a non-uniform state and R*
//     the residual at q* without the source: the residual with the source at q = q* is 0 on the rows where w = 1
//     that keep their residual (there |R*| is not small: the plain relaxation leaves R*), (1 - w) R* on the margin,
//     R* (bit for bit) where c = 0 and on the rows that replace theirs; continuity has no source.
//   the implicit rows are those of relax, bit for bit (the added term is a constant). A guard: the
//     implicit rows have no relax_consistent path today, so no mutant of tools/nest_relax_mutants.py reaches it.
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
#include <cstring>
#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

namespace {

using gust::Equation;
using gust::Loc;

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double kMargin = 2;   // kDefectMargin, written out so a change of the constant shows here
constexpr double kFactor = 3;

bool bits_differ (amrex::Real a, amrex::Real b) { return std::memcmp(&a, &b, sizeof(amrex::Real)) != 0; }

double pos (const amrex::Geometry& g, amrex::IndexType t, int i, int d)
{
    const double off = t.nodeCentered(d) ? 0.0 : 0.5;
    return static_cast<double>(g.ProbLo(d)) + (i + off) * static_cast<double>(g.CellSize(d));
}

double field (double x, double y, int k, int c, int s)
{
    return (4.0 + 0.5 * c - 0.7 * s) + 0.2 * k + 0.003 * x * (1 + c) - 0.002 * y * (1 + s)
         + 1.0e-6 * (1 + 0.1 * k) * x * y + 0.4 * std::sin(0.004 * x + 0.3 * c + s) * std::cos(0.003 * y);
}

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

// fn(i, j, k, a, b) on every valid point of a and b (same layout), through host copies.
void each_valid2 (const amrex::MultiFab& a, const amrex::MultiFab& b,
                  const std::function<void(int, int, int, double, double)>& fn)
{
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        amrex::FArrayBox ha(bx, 1, amrex::The_Pinned_Arena()), hb(bx, 1, amrex::The_Pinned_Arena());
        ha.copy<amrex::RunOn::Device>(a[mfi], bx);
        hb.copy<amrex::RunOn::Device>(b[mfi], bx);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            fn(i, j, k, static_cast<double>(ha(iv)), static_cast<double>(hb(iv)));
        });
    }
}

// --- A flat parent --------------------------------------------------------------------------------

struct Parent
{
    gust::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;

    Parent ()
    {
        // 20 x 16 x 8 cells of 100 m x 80 m (dx != dy), uneven boxes.
        in.mesh.n_cell = amrex::IntVect(20, 16, 8);
        in.mesh.prob_lo = {-1000, 500, 0};
        in.mesh.prob_hi = {1000, 1780, 800};
        in.mesh.is_periodic = {0, 0, 0};
        in.mesh.max_grid_size = {7, 6};
        in.mesh.dz0 = 40;
        geom = gust::make_geometry(in.mesh);
        ba = gust::make_box_array(in.mesh);
        dm = amrex::DistributionMapping(ba);
    }
};

// w at point (i, j) of index type t from the points' positions: the largest min(1, s / 2) over the footprints
// holding the point (s: the distance from the nearest face, over dx or dy); -1 outside every footprint.
double expect_weight (const Parent& p, const std::vector<gust::DefectRegion>& regions, amrex::IndexType t, int i, int j)
{
    double w = -1;
    for (const auto& g : regions) {
        const amrex::Box& fp = g.footprint;
        double s = std::numeric_limits<double>::infinity();
        const int ij[2] = {i, j};
        bool inside = true;
        for (int d = 0; d < 2; ++d) {
            const double h = static_cast<double>(p.geom.CellSize(d));
            const double lo = static_cast<double>(p.geom.ProbLo(d)) + fp.smallEnd(d) * h;
            const double hi = static_cast<double>(p.geom.ProbLo(d)) + (fp.bigEnd(d) + 1) * h;
            const double x = pos(p.geom, t, ij[d], d);
            if (x < lo - 1e-9 * h || x > hi + 1e-9 * h) { inside = false; }
            s = std::min(s, std::min(x - lo, hi - x) / h);
        }
        if (inside) { w = std::max(w, std::min(1.0, std::max(0.0, s / kMargin))); }
    }
    return w;
}

// --- The hill as the application sets it up ------------------------------------------------------------

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

// q = q0 + s dq, dq smooth, non-uniform and nonzero everywhere (p too); ghosts refilled.
void perturb (HillRun& h, const gust::State& q0, gust::State& q, double s)
{
    q.define(h.c.ba, h.c.dm);
    auto t = q.fields();
    const amrex::MultiFab* f0[7] = {&q0.u, &q0.v, &q0.w, &q0.p, &q0.theta, &q0.tke, &q0.z_phys_nd};
    for (std::size_t n = 0; n < static_cast<std::size_t>(t.size()); ++n) {
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

TEST(NestRelaxConsistent, AddedTermIsTheWeightedResidualOnUnevenBoxes)
{
    ASSERT_EQ(gust::kDefectMargin, 2);
    Parent p;
    const amrex::Box dom = p.geom.Domain();
    // A: parent cells 3-14 x 2-11; B: 12-18 x 8-14, overlapping A.
    const std::vector<gust::DefectRegion> regions{
        gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), 2, amrex::IntVect(3, 2, 0), dom),
        gust::defect_region("nest.b", amrex::IntVect(14, 14, 8), 2, amrex::IntVect(12, 8, 0), dom)};
    gust::State qstar;
    qstar.define(p.ba, p.dm);
    amrex::MultiFab* qf[5] = {&qstar.u, &qstar.v, &qstar.w, &qstar.theta, &qstar.tke};
    for (int c = 0; c < 5; ++c) {
        const amrex::IndexType t = qf[c]->ixType();
        fill_all(*qf[c], [&] (int i, int j, int k) { return field(pos(p.geom, t, i, 0), pos(p.geom, t, j, 1), k, c, 1); });
    }
    // R*: non-uniform on the footprints, 1e30 outside them (it must not reach the added term).
    gust::ResidualFields Rs;
    Rs.define(p.ba, p.dm);
    amrex::MultiFab* rf[5] = {&Rs.u, &Rs.v, &Rs.w, &Rs.theta, &Rs.tke};
    for (int c = 0; c < 5; ++c) {
        const amrex::IndexType t = rf[c]->ixType();
        fill_all(*rf[c], [&] (int i, int j, int k) {
            return expect_weight(p, regions, t, i, j) < 0 ? 1.0e30 : field(pos(p.geom, t, i, 0), pos(p.geom, t, j, 1), k, c, 2);
        });
    }
    gust::NestRelax plain(p.ba, p.dm, p.geom, amrex::Real(kFactor));
    gust::NestRelax S(p.ba, p.dm, p.geom, amrex::Real(kFactor), true);
    plain.set(qstar, regions);
    S.set(qstar, regions, Rs);
    ASSERT_TRUE(S.consistent());
    const Loc locs[4] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell};
    for (const Loc l : locs) {
        EXPECT_EQ(gust::test::count_bit_differences(S.rate(l), plain.rate(l), 0), 0) << "the rates are not relax's";
    }
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    const double tol = kDouble ? 1e-14 : 1e-6;
    for (int n = 0; n < 5; ++n) {
        const amrex::IndexType t = rf[n]->ixType();
        double err = 0, rmax = 0;   // |a - w R*| against the largest |R*| (a relative error per point divides by 0 where R* crosses 0)
        long edge = 0, full = 0, off_nonzero = 0;
        each_valid2(S.added(eqs[n]), *rf[n], [&] (int i, int j, int, double a, double r) {
            const double w = expect_weight(p, regions, t, i, j);
            if (w <= 0) { if (a != 0) { ++off_nonzero; } return; }
            if (w < 1) { ++edge; } else { ++full; }
            err = std::max(err, std::abs(a - w * r));
            rmax = std::max(rmax, std::abs(r));
        });
        auto comm = amrex::ParallelDescriptor::Communicator();
        amrex::ParallelAllReduce::Max(err, comm);
        amrex::ParallelAllReduce::Max(rmax, comm);
        amrex::ParallelAllReduce::Sum(edge, comm);
        amrex::ParallelAllReduce::Sum(full, comm);
        amrex::ParallelAllReduce::Sum(off_nonzero, comm);
        EXPECT_GT(rmax, 0.0) << names[n];
        EXPECT_LE(err, tol * rmax) << names[n] << ": the added term is not w R*";
        EXPECT_EQ(off_nonzero, 0) << names[n] << ": an added term on a footprint face or outside the footprints";
        EXPECT_GT(edge, 0) << names[n];
        EXPECT_GT(full, 0) << names[n];
    }
}

TEST(NestRelaxConsistent, TargetIsAFixedPointWhereTheWeightIsOne)
{
    HillRun h = make_hill();
    gust::State qstar;
    perturb(h, h.q, qstar, 1.0);   // q*: a non-uniform state, not a solution
    gust::ResidualInputs off = h.ls->res_in;
    ASSERT_EQ(off.relax, nullptr);
    gust::ResidualFields Rs;
    gust::residual(qstar, off, h.c.mesh, h.c.geom, Rs);   // R* = R_c(q*)
    gust::NestRelax S(h.c.ba, h.c.dm, h.c.geom, amrex::Real(kFactor), true);
    S.set(qstar, {hill_region()}, Rs);
    gust::NestRelax plain(h.c.ba, h.c.dm, h.c.geom, amrex::Real(kFactor));
    plain.set(qstar, {hill_region()});
    gust::ResidualInputs on = off;
    on.relax = &S;
    gust::ResidualInputs on_plain = off;
    on_plain.relax = &plain;
    gust::ResidualFields Ron, Rplain;
    gust::residual(qstar, on, h.c.mesh, h.c.geom, Ron);
    gust::residual(qstar, on_plain, h.c.mesh, h.c.geom, Rplain);
    ASSERT_EQ(S.nests().size(), 1U);
    const double c0 = static_cast<double>(S.nests()[0].rate);
    ASSERT_GT(c0, 0.0);
    const Loc locs[5] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell};
    const int fi[5] = {0, 1, 2, 4, 5};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    // The rows that replace their residual: the ground and top w faces, the first cells' wall-k rows.
    const long expect_replaced[5] = {0, 0, 2 * footprint_columns(), 0, footprint_columns()};
    const double tol = kDouble ? 1e-12 : 1e-5;
    for (int n = 0; n < 5; ++n) {
        const amrex::MultiFab& c = S.rate(locs[n]);
        long outside_changed = 0, replaced = 0, full = 0, edge = 0, plain_moved = 0;
        double full_err = 0, edge_err = 0, rs_full = 0, rs_edge = 0;
        for (amrex::MFIter mfi(c, false); mfi.isValid(); ++mfi) {
            const amrex::Box b = mfi.validbox();
            auto host = [&] (const amrex::MultiFab& mf) {
                amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
                f.copy<amrex::RunOn::Device>(mf[mfi], b, 0, b, 0, 1);
                return f;
            };
            const auto hc = host(c), ha = host(rfield(Ron, fi[n])), hs = host(rfield(Rs, fi[n]));
            const auto hp = host(rfield(Rplain, fi[n]));
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                const double cc = static_cast<double>(hc(iv));
                const double ron = static_cast<double>(ha(iv)), rs = static_cast<double>(hs(iv));
                if (cc == 0) {
                    if (bits_differ(ha(iv), hs(iv))) { ++outside_changed; }
                    return;
                }
                // A replaced row: the boundary row of q*, equal with and without the source (elsewhere the source
                // changes it: by -R* where w = 1, by -w R* on the margin).
                if (!bits_differ(ha(iv), hs(iv))) {
                    ++replaced;
                    return;
                }
                if (bits_differ(hp(iv), hs(iv))) { ++plain_moved; }   // relax adds c (q* - q*) = 0: it leaves R*
                const double w = cc / c0;
                if (cc == c0) {
                    ++full;
                    rs_full = std::max(rs_full, std::abs(rs));
                    full_err = std::max(full_err, std::abs(ron));
                } else {
                    ++edge;
                    edge_err = std::max(edge_err, std::abs(ron - (1 - w) * rs));
                    rs_edge = std::max(rs_edge, std::abs(rs));
                }
            });
        }
        auto comm = amrex::ParallelDescriptor::Communicator();
        amrex::ParallelAllReduce::Sum(outside_changed, comm);
        amrex::ParallelAllReduce::Sum(replaced, comm);
        amrex::ParallelAllReduce::Sum(full, comm);
        amrex::ParallelAllReduce::Sum(edge, comm);
        amrex::ParallelAllReduce::Sum(plain_moved, comm);
        amrex::ParallelAllReduce::Max(full_err, comm);
        amrex::ParallelAllReduce::Max(edge_err, comm);
        amrex::ParallelAllReduce::Max(rs_full, comm);
        amrex::ParallelAllReduce::Max(rs_edge, comm);
        EXPECT_EQ(outside_changed, 0) << names[n] << ": the residual changed where c = 0";
        EXPECT_EQ(replaced, expect_replaced[n]) << names[n] << ": rows without the source";
        EXPECT_EQ(plain_moved, 0) << names[n] << ": relax at q = q* is not R*";
        EXPECT_GT(full, 0) << names[n];
        EXPECT_GT(edge, 0) << names[n];
        // Not vacuous: R* is not small on those rows, so q* is not a fixed point of relax there.
        EXPECT_GT(rs_full, 1e-6) << names[n] << ": R* is too small for the test";
        EXPECT_LE(full_err, tol * rs_full) << names[n] << ": q* is not a fixed point where w = 1";
        EXPECT_LE(edge_err, tol * rs_edge) << names[n] << ": the margin's residual at q* is not (1 - w) R*";
    }
    // Continuity has no source.
    EXPECT_EQ(gust::test::count_bit_differences(Ron.p, Rs.p, 0), 0);
}

TEST(NestRelaxConsistent, ImplicitRowsAreThoseOfRelax)
{
    HillRun h = make_hill();
    gust::State qstar, q;
    perturb(h, h.q, qstar, 1.0);
    perturb(h, qstar, q, 0.5);
    gust::ResidualFields Rs;
    gust::residual(qstar, h.ls->res_in, h.c.mesh, h.c.geom, Rs);
    gust::NestRelax S(h.c.ba, h.c.dm, h.c.geom, amrex::Real(kFactor), true);
    S.set(qstar, {hill_region()}, Rs);
    gust::NestRelax plain(h.c.ba, h.c.dm, h.c.geom, amrex::Real(kFactor));
    plain.set(qstar, {hill_region()});
    gust::ResidualInputs on = h.ls->res_in, on_plain = h.ls->res_in;
    on.relax = &S;
    on_plain.relax = &plain;
    const Equation eqs[5] = {Equation::u, Equation::v, Equation::w, Equation::theta, Equation::tke};
    for (const Equation eq : eqs) {
        amrex::MultiFab A, Ap;
        rows(h, q, on, eq, A);
        rows(h, q, on_plain, eq, Ap);
        EXPECT_EQ(gust::test::count_bit_differences(A, Ap, 0), 0) << static_cast<int>(eq);
    }
}
