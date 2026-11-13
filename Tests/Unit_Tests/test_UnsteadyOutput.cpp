// Unit tests of the unsteady outputs (UnsteadyOutput.H, sub-phase 6.3, beyond
// the plan): the running averages against a two-pass weighted series, the
// probes' cells and values, the probe file on a restart, the checkpoint record
// and the averages' checkpoint on another layout, and the time of a resumed
// TimeDriver.

#include "Driver.H"
#include "LaminarTestHelpers.H"
#include "Simplec.H"
#include "TimeDriver.H"
#include "TimeTerm.H"
#include "TransportTestHelpers.H"
#include "UnsteadyOutput.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace {

using gust::lamtest::define_state;
using gust::lamtest::flat_case;
using gust::lamtest::Laminar;
using gust::ptest::PCase;
using gust::trtest::tol;

// A non-uniform state at time t: u = sin(0.7 i + 0.3 j + t) + 0.1 k on x-faces, v = cos(...) on y-faces,
// w = 0.01 k t on z-faces, theta = 300 + 0.5 sin(i + 2 t), k = 0.2 + 0.05 (j + k) t^2.
void set_state (gust::State& q, double t)
{
    const auto tt = static_cast<amrex::Real>(t);
    for (amrex::MFIter mfi(q.theta, false); mfi.isValid(); ++mfi) {
        const auto u = q.u.array(mfi);
        const auto v = q.v.array(mfi);
        const auto w = q.w.array(mfi);
        const auto th = q.theta.array(mfi);
        const auto k = q.tke.array(mfi);
        amrex::ParallelFor(amrex::surroundingNodes(mfi.validbox(), 0), [=] AMREX_GPU_DEVICE (int i, int j, int kk) {
            u(i, j, kk) = std::sin(amrex::Real(0.7) * i + amrex::Real(0.3) * j + tt) + amrex::Real(0.1) * kk;
        });
        amrex::ParallelFor(amrex::surroundingNodes(mfi.validbox(), 1), [=] AMREX_GPU_DEVICE (int i, int j, int kk) {
            v(i, j, kk) = std::cos(amrex::Real(0.7) * i + amrex::Real(0.3) * j + tt);
        });
        amrex::ParallelFor(amrex::surroundingNodes(mfi.validbox(), 2), [=] AMREX_GPU_DEVICE (int i, int j, int kk) {
            w(i, j, kk) = amrex::Real(0.01) * kk * tt;
        });
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int kk) {
            th(i, j, kk) = amrex::Real(300) + amrex::Real(0.5) * std::sin(amrex::Real(i) + 2 * tt);
            k(i, j, kk) = amrex::Real(0.2) + amrex::Real(0.05) * (j + kk) * tt * tt;
        });
    }
}

// The cell-centred fields (as the averages see them), into comps 0..4 of out.
void centre_fields (const gust::State& q, amrex::MultiFab& out)
{
    for (amrex::MFIter mfi(out, false); mfi.isValid(); ++mfi) {
        const auto o = out.array(mfi);
        const auto u = q.u.const_array(mfi);
        const auto v = q.v.const_array(mfi);
        const auto w = q.w.const_array(mfi);
        const auto th = q.theta.const_array(mfi);
        const auto k = q.tke.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int kk) {
            constexpr amrex::Real half = amrex::Real(0.5);
            o(i, j, kk, 0) = half * (u(i, j, kk) + u(i + 1, j, kk));
            o(i, j, kk, 1) = half * (v(i, j, kk) + v(i, j + 1, kk));
            o(i, j, kk, 2) = half * (w(i, j, kk) + w(i, j, kk + 1));
            o(i, j, kk, 3) = th(i, j, kk);
            o(i, j, kk, 4) = k(i, j, kk);
        });
    }
}

// acc += wt (x - shift)^2 comp by comp (shift = 0 for the first pass's acc += wt x with power 1).
void accumulate (amrex::MultiFab& acc, const amrex::MultiFab& x, const amrex::MultiFab* shift, amrex::Real wt)
{
    for (amrex::MFIter mfi(acc, false); mfi.isValid(); ++mfi) {
        const auto a = acc.array(mfi);
        const auto xx = x.const_array(mfi);
        if (shift == nullptr) {
            amrex::ParallelFor(mfi.validbox(), gust::kUnsteadyFields, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                a(i, j, k, n) += wt * xx(i, j, k, n);
            });
        } else {
            const auto s = shift->const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), gust::kUnsteadyFields, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) {
                const amrex::Real d = xx(i, j, k, n) - s(i, j, k, n);
                a(i, j, k, n) += wt * d * d;
            });
        }
    }
}

// max over cells and comps of |a - b| / (rel |b| + abs_floor).
double max_rel_diff (const amrex::MultiFab& a, const amrex::MultiFab& b, int acomp, int bcomp, double rel, double abs_floor)
{
    amrex::ReduceOps<amrex::ReduceOpMax> op;
    amrex::ReduceData<amrex::Real> data(op);
    using T = amrex::ReduceData<amrex::Real>::Type;
    const auto r = amrex::Real(rel);
    const auto f = amrex::Real(abs_floor);
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        const auto aa = a.const_array(mfi);
        const auto bb = b.const_array(mfi);
        op.eval(mfi.validbox(), gust::kUnsteadyFields, data, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n) -> T {
            return {amrex::Math::abs(aa(i, j, k, acomp + n) - bb(i, j, k, bcomp + n))
                    / (r * amrex::Math::abs(bb(i, j, k, bcomp + n)) + f)};
        });
    }
    amrex::Real m = amrex::get<0>(data.value(op));
    amrex::ParallelAllReduce::Max(m, amrex::ParallelDescriptor::Communicator());
    return double(m);
}

// The centre heights above the ground of column (i, j), from the node heights
// on the host (the mean of the 8 nodes minus the mean of the 4 bottom nodes).
std::vector<double> column_heights (const PCase& c, int i, int j)
{
    const int nz = c.geom.Domain().length(2);
    std::vector<double> h(nz, 0);
    for (amrex::MFIter mfi(c.z, false); mfi.isValid(); ++mfi) {
        const amrex::Box cells = amrex::enclosedCells(mfi.validbox());
        if (!cells.contains(amrex::IntVect(i, j, 0))) { continue; }
        amrex::FArrayBox f(c.z[mfi].box(), 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(c.z[mfi], c.z[mfi].box());
        amrex::Gpu::streamSynchronize();
        const auto z = f.const_array();
        const double g = 0.25 * (double(z(i, j, 0)) + z(i + 1, j, 0) + z(i, j + 1, 0) + z(i + 1, j + 1, 0));
        for (int k = 0; k < nz; ++k) {
            double s = 0;
            for (int dk = 0; dk <= 1; ++dk) {
                s += double(z(i, j, k + dk)) + z(i + 1, j, k + dk) + z(i, j + 1, k + dk) + z(i + 1, j + 1, k + dk);
            }
            h[k] = 0.125 * s - g;
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(h.data(), nz);
    return h;
}

bool io_rank () { return amrex::ParallelDescriptor::IOProcessor(); }

// true on every rank when ok is true on the I/O rank (EXPECTs on other ranks do not fail a test).
bool from_io (bool ok)
{
    int v = ok ? 1 : 0;
    amrex::ParallelDescriptor::Bcast(&v, 1, amrex::ParallelDescriptor::IOProcessorNumber());
    return v != 0;
}

std::string slurp (const std::string& file)
{
    std::ifstream is(file);
    std::stringstream ss;
    ss << is.rdbuf();
    return ss.str();
}

} // namespace

// UT-6.3.1: the running averages after irregular steps, the averaging start in
// the middle of the second step, equal the two-pass weighted mean and variance
// with the weights t^{n+1} - max(t^n, start), cell by cell, on several boxes.
TEST(UnsteadyOutput, AveragesAreTheWeightedSeries)
{
    const PCase c = flat_case(6, 4, 5, 1.0, 1.0, 1.0, true, true, 3);
    gust::State q;
    define_state(q, c);
    const std::vector<double> t = {0.0, 0.7, 1.9, 2.0, 3.5, 6.0, 6.25};
    const double start = 1.0;
    gust::RunningAverages avg(amrex::Real(start), c.ba, c.dm);
    amrex::MultiFab x(c.ba, c.dm, gust::kUnsteadyFields, 0);
    amrex::MultiFab sum(c.ba, c.dm, gust::kUnsteadyFields, 0);
    sum.setVal(0);
    double total = 0;
    std::vector<double> weights;
    for (std::size_t n = 1; n < t.size(); ++n) {
        set_state(q, t[n]);
        avg.add(q, t[n - 1], t[n]);
        const double w = t[n] - std::max(t[n - 1], start);
        weights.push_back(w);
        if (w > 0) {
            centre_fields(q, x);
            accumulate(sum, x, nullptr, amrex::Real(w));
            total += w;
        }
    }
    EXPECT_NEAR(weights[0], -0.3, 1e-15);
    EXPECT_NEAR(weights[1], 0.9, 1e-15);
    EXPECT_EQ(avg.weight(), total);
    amrex::MultiFab mean(c.ba, c.dm, gust::kUnsteadyFields, 0);
    amrex::MultiFab::Copy(mean, sum, 0, 0, gust::kUnsteadyFields, 0);
    mean.mult(amrex::Real(1.0 / total));
    amrex::MultiFab ss(c.ba, c.dm, gust::kUnsteadyFields, 0);
    ss.setVal(0);
    for (std::size_t n = 1; n < t.size(); ++n) {
        if (!(weights[n - 1] > 0)) { continue; }
        set_state(q, t[n]);
        centre_fields(q, x);
        accumulate(ss, x, &mean, amrex::Real(weights[n - 1]));
    }
    ss.mult(amrex::Real(1.0 / total));   // the variance

    amrex::MultiFab plot;
    std::vector<std::string> names;
    avg.plot_fields(plot, names);
    ASSERT_EQ(names.size(), 10U);
    EXPECT_EQ(names[0], "u_mean");
    EXPECT_EQ(names[4], "tke_mean");
    EXPECT_EQ(names[5], "u_var");
    EXPECT_EQ(names[9], "tke_var");
    EXPECT_LE(max_rel_diff(plot, mean, 0, 0, 1.0, 1.0), tol(1e-14, 1e-6));
    EXPECT_LE(max_rel_diff(plot, ss, gust::kUnsteadyFields, 0, 1.0, 1e-12), tol(1e-11, 1e-3));
}

// No sample before the averaging start: the mean and the variance stay 0 and W = 0.
TEST(UnsteadyOutput, NothingAveragedBeforeTheStart)
{
    const PCase c = flat_case(4, 4, 4, 1.0, 1.0, 1.0, true, true, 2);
    gust::State q;
    define_state(q, c);
    gust::RunningAverages avg(amrex::Real(5), c.ba, c.dm);
    set_state(q, 1.0);
    avg.add(q, 0.0, 2.5);
    avg.add(q, 2.5, 5.0);
    EXPECT_EQ(avg.weight(), 0.0);
    amrex::MultiFab plot;
    std::vector<std::string> names;
    avg.plot_fields(plot, names);
    for (int n = 0; n < plot.nComp(); ++n) { EXPECT_EQ(plot.norm0(n), amrex::Real(0)) << names[n]; }
}

// The averages' checkpoint read on another box layout, and the Unsteady record
// with values that have no short decimal form, come back bit for bit; a steady
// checkpoint has no record.
TEST(UnsteadyOutput, CheckpointRoundTrip)
{
    const PCase c = flat_case(6, 4, 5, 1.0, 1.0, 1.0, true, true, 3);
    gust::State q;
    define_state(q, c);
    gust::RunningAverages avg(amrex::Real(0), c.ba, c.dm);
    for (int n = 1; n <= 3; ++n) {
        set_state(q, 0.37 * n);
        avg.add(q, 0.37 * (n - 1), 0.37 * n);
    }
    const std::string dir = "ut_unsteady_output_chk";
    amrex::UtilCreateCleanDirectory(dir, true);
    amrex::UtilCreateCleanDirectory(dir + "/Level_0", true);
    avg.write_checkpoint(dir);
    gust::UnsteadyRecord r;
    r.time = 1.0 / 3.0;
    r.steps = 123456789012L;
    r.dt = 0.1;
    r.average_start = 2.0 / 7.0;
    r.weight = avg.weight();
    r.anchor_time = 1.0 / 7.0;
    r.anchor_steps = 98765432101L;
    gust::write_unsteady_record(dir, r);
    amrex::ParallelDescriptor::Barrier();

    amrex::BoxArray other(c.geom.Domain());
    other.maxSize(amrex::IntVect(2, 3, 5));
    const amrex::DistributionMapping odm(other);
    gust::RunningAverages back(amrex::Real(0), other, odm);
    gust::UnsteadyRecord rb;
    ASSERT_TRUE(gust::read_unsteady_record(dir, rb));
    back.read_checkpoint(dir, rb.weight);
    EXPECT_EQ(rb.time, r.time);
    EXPECT_EQ(rb.steps, r.steps);
    EXPECT_EQ(rb.dt, r.dt);
    EXPECT_EQ(rb.average_start, r.average_start);
    EXPECT_EQ(rb.weight, avg.weight());
    EXPECT_EQ(rb.anchor_time, r.anchor_time);
    EXPECT_EQ(rb.anchor_steps, r.anchor_steps);
    EXPECT_EQ(back.weight(), avg.weight());
    amrex::MultiFab m(c.ba, c.dm, gust::kUnsteadyFields, 0);
    amrex::MultiFab s(c.ba, c.dm, gust::kUnsteadyFields, 0);
    m.ParallelCopy(back.mean(), 0, 0, gust::kUnsteadyFields);
    s.ParallelCopy(back.m2(), 0, 0, gust::kUnsteadyFields);
    EXPECT_EQ(max_rel_diff(m, avg.mean(), 0, 0, 0.0, 1.0), 0.0);
    EXPECT_EQ(max_rel_diff(s, avg.m2(), 0, 0, 0.0, 1.0), 0.0);

    const std::string steady = "ut_unsteady_output_steady";
    amrex::UtilCreateCleanDirectory(steady, true);
    gust::UnsteadyRecord none;
    EXPECT_FALSE(gust::read_unsteady_record(steady, none));
    amrex::ParallelDescriptor::Barrier();
    if (io_rank()) {
        std::filesystem::remove_all(dir);
        std::filesystem::remove_all(steady);
    }
}

// The probe file's lines: comments and blank lines skipped, every bad line named.
TEST(UnsteadyOutput, ParseProbes)
{
    std::vector<std::string> err;
    const auto p = gust::parse_probes("# header\n\n 1 2 3   # a point\n-4.5 0 1e1\n", "f.txt", err);
    EXPECT_TRUE(err.empty());
    ASSERT_EQ(p.size(), 2U);
    EXPECT_EQ(p[1].x, amrex::Real(-4.5));
    EXPECT_EQ(p[1].h, amrex::Real(10));

    err.clear();
    gust::parse_probes("1 2\n1 2 x\n1 2 nan\n1 2 0\n1 2 3 4\n", "f.txt", err);
    ASSERT_EQ(err.size(), 5U);
    EXPECT_NE(err[0].find("unsteady.probes_file = f.txt, line 1: expected three numbers"), std::string::npos) << err[0];
    EXPECT_NE(err[1].find("line 2: 'x' is not a finite number"), std::string::npos) << err[1];
    EXPECT_NE(err[2].find("line 3: 'nan' is not a finite number"), std::string::npos) << err[2];
    EXPECT_NE(err[3].find("line 4: h = 0: the height above the ground must be > 0"), std::string::npos) << err[3];
    EXPECT_NE(err[4].find("line 5: expected three numbers"), std::string::npos) << err[4];

    err.clear();
    gust::parse_probes("# nothing\n", "f.txt", err);
    ASSERT_EQ(err.size(), 1U);
    EXPECT_NE(err[0].find("no probe"), std::string::npos);
}

// UT-6.3.2: each probe takes its column's cell whose centre height above the
// ground is nearest h, the lower one on a tie, and samples that cell's values.
// Flat ground at 0, 8 layers of 0.125 m: centres 0.0625 + 0.125 k. Then a
// bump, where h counts from the ground under the probe.
TEST(UnsteadyOutput, ProbesTakeTheNearestCell)
{
    {
        const PCase c = flat_case(4, 4, 8, 1.0, 1.0, 1.0, true, true, 2);
        std::vector<gust::Probe> list(4);
        list[0].x = amrex::Real(0.1); list[0].y = amrex::Real(0.6); list[0].h = amrex::Real(0.3);   // k = 2 (0.3125)
        list[1].x = amrex::Real(0.9); list[1].y = amrex::Real(0.3); list[1].h = amrex::Real(0.25);  // tie 0.1875 / 0.3125 -> k = 1
        list[2].x = amrex::Real(1.0); list[2].y = amrex::Real(1.0); list[2].h = amrex::Real(1.0);   // the high corner, the top cell
        list[3].x = amrex::Real(0.0); list[3].y = amrex::Real(0.0); list[3].h = amrex::Real(0.01);  // the low corner, k = 0
        const gust::Probes probes(list, c.z, c.geom, "f.txt");
        const auto& p = probes.probes();
        EXPECT_EQ(p[0].cell, amrex::IntVect(0, 2, 2));
        EXPECT_EQ(p[0].h_cell, amrex::Real(0.3125));
        EXPECT_EQ(p[1].cell, amrex::IntVect(3, 1, 1));
        EXPECT_EQ(p[1].h_cell, amrex::Real(0.1875));
        EXPECT_EQ(p[2].cell, amrex::IntVect(3, 3, 7));
        EXPECT_EQ(p[3].cell, amrex::IntVect(0, 0, 0));

        gust::State q;
        define_state(q, c);
        set_state(q, 0.5);
        const std::vector<amrex::Real> s = probes.sample(q);
        ASSERT_EQ(s.size(), 4U * gust::kUnsteadyFields);
        // Probe 1, cell (3, 1, 1): the formulas of set_state.
        const double t = 0.5;
        const double u = 0.5 * (std::sin(0.7 * 3 + 0.3 + t) + std::sin(0.7 * 4 + 0.3 + t)) + 0.1;
        const double v = 0.5 * (std::cos(0.7 * 3 + 0.3 + t) + std::cos(0.7 * 3 + 0.6 + t));
        const double w = 0.5 * 0.01 * (1 + 2) * t;
        const double th = 300 + 0.5 * std::sin(3 + 2 * t);
        const double k = 0.2 + 0.05 * (1 + 1) * t * t;
        const amrex::Real* s1 = s.data() + gust::kUnsteadyFields;
        EXPECT_NEAR(s1[0], u, tol(1e-15, 1e-6));
        EXPECT_NEAR(s1[1], v, tol(1e-15, 1e-6));
        EXPECT_NEAR(s1[2], w, tol(1e-15, 1e-7));
        EXPECT_NEAR(s1[3], th, tol(1e-13, 1e-4));
        EXPECT_NEAR(s1[4], k, tol(1e-15, 1e-6));
    }
    {
        // A bump of 0.3 m: the layers follow the ground (not uniformly), so the
        // expected cell comes from a host scan of the column's node heights.
        // Measured from z = 0 instead of the ground, h = 0.3 would take a cell
        // near the bottom of the column under the bump's top.
        const PCase c = flat_case(4, 4, 8, 1.0, 1.0, 1.0, true, true, 2,
                                  [] (double x, double y) { return 0.3 * std::exp(-20 * ((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5))); });
        std::vector<gust::Probe> list(1);
        list[0].x = amrex::Real(0.4); list[0].y = amrex::Real(0.6); list[0].h = amrex::Real(0.3);
        const gust::Probes probes(list, c.z, c.geom, "f.txt");
        const std::vector<double> hc = column_heights(c, 1, 2);
        int kbest = 0;
        for (int k = 1; k < int(hc.size()); ++k) {
            if (std::abs(hc[k] - 0.3) < std::abs(hc[kbest] - 0.3)) { kbest = k; }
        }
        EXPECT_GT(hc[0], 0.0);
        EXPECT_GE(kbest, 2);
        EXPECT_EQ(probes.probes()[0].cell, amrex::IntVect(1, 2, kbest));
        EXPECT_NEAR(probes.probes()[0].h_cell, hc[kbest], tol(1e-14, 1e-6));
    }
}

// The probe file: a new run writes the header and its rows; a restart at t = 1
// keeps the rows up to t = 1 and appends after them; a restart that finds no
// file starts a new one (open returns false: the caller writes the first row).
TEST(UnsteadyOutput, ProbeFileOnRestart)
{
    const PCase c = flat_case(4, 4, 8, 1.0, 1.0, 1.0, true, true, 2);
    std::vector<gust::Probe> list(1);
    list[0].x = amrex::Real(0.5); list[0].y = amrex::Real(0.5); list[0].h = amrex::Real(0.5);
    const gust::Probes probes(list, c.z, c.geom, "f.txt");
    gust::State q;
    define_state(q, c);
    const std::string out = "ut_unsteady_output_probes.csv";
    EXPECT_FALSE(probes.open(out, "f.txt", -1.0));
    for (int n = 0; n <= 3; ++n) {
        set_state(q, double(n));
        probes.write_row(out, double(n), q);
    }
    const std::string straight = io_rank() ? slurp(out) : std::string();
    EXPECT_TRUE(probes.open(out, "f.txt", 1.0));
    const std::string kept = io_rank() ? slurp(out) : std::string();
    for (int n = 2; n <= 3; ++n) {
        set_state(q, double(n));
        probes.write_row(out, double(n), q);
    }
    const std::string again = io_rank() ? slurp(out) : std::string();
    bool ok = true;
    if (io_rank()) {
        const std::string head = probes.header();
        ok = straight.rfind(head, 0) == 0 && kept.rfind(head, 0) == 0
             && std::count(straight.begin(), straight.end(), '\n') == std::count(head.begin(), head.end(), '\n') + 4
             && std::count(kept.begin(), kept.end(), '\n') == std::count(head.begin(), head.end(), '\n') + 2
             && again == straight;
        std::filesystem::remove(out);
    }
    EXPECT_TRUE(from_io(ok));
    amrex::ParallelDescriptor::Barrier();
    EXPECT_FALSE(probes.open(out, "f.txt", 1.0));
    bool fresh = true;
    if (io_rank()) {
        fresh = slurp(out) == probes.header();
        std::filesystem::remove(out);
    }
    EXPECT_TRUE(from_io(fresh));
}

// UT-6.3.3: a resumed TimeDriver numbers its steps on from the checkpoint's and,
// with the checkpoint's fixed step and time anchor, gives the times of the run
// without the restart, bit for bit. With another fixed step the time is
// t_ref + n dt. After a clipped last step the next step is one dt later, and a
// chain of restarts gives the times of a single restart.
TEST(TimeDriver, ResumeContinuesTheTime)
{
    auto make = [] (double stop, double dt) {
        gust::UnsteadyOptions uo;
        uo.stop_time = amrex::Real(stop);
        uo.dt = amrex::Real(dt);
        uo.max_inner = 1;
        uo.given = {"stop_time", "dt", "max_inner"};
        return uo;
    };
    // A channel at rest (no flow, one inner iteration a step): the driver's
    // time bookkeeping alone.
    const PCase c = flat_case(4, 4, 8, 1.0, 1.0, 1.0, true, true, 2);
    Laminar lam;
    lam.make(c, 0.01, gust::FaceType::no_slip_wall, gust::FaceType::no_slip_wall);
    lam.in.jvp_scales.velocity = 1;
    lam.in.jvp_scales.pressure = 1;
    // What a checkpoint records of the driver (UnsteadyRecord) and the steps run.
    struct Ran {
        std::vector<gust::StepRecord> records;
        double time{0};
        long steps{0};
        double dt{0};
        double anchor_time{0};
        long anchor_steps{0};
    };
    auto run = [&] (const gust::UnsteadyOptions& uo, const Ran* from) {
        gust::State q;
        define_state(q, c);
        q.theta.setVal(300.0);
        q.tke.setVal(1.0);
        gust::TimeTerm term(gust::TimeTermOptions{uo.dt}, c.ba, c.dm);
        lam.in.time = &term;
        gust::SteadyOptions so;
        so.stage = gust::StageKind::unsteady;
        so.zline_sweeps = 2;
        so.criterion.residual_tol = amrex::Real(1e-6);
        gust::SimplecOptions opt;
        auto stage = std::make_unique<gust::SimplecStage>(opt, so.zline_sweeps, lam.in, nullptr, lam.popt, c.mesh, c.geom,
                                                              c.z, 0, c.ba, c.dm);
        const gust::SimplecStage& simplec = *stage;
        std::vector<std::unique_ptr<gust::Stage>> stages;
        stages.push_back(std::move(stage));
        gust::Driver driver(so, std::move(stages), c.mesh, c.geom, &simplec.nu());
        gust::BoundaryConditions* bc = lam.bc.get();
        gust::TimeDriver td(uo, term, driver, c.mesh, c.geom, amrex::Real(1e-6), [bc] (gust::State& s) {
            gust::BuoyancyOptions buoy;
            buoy.gravity = 0;
            bc->fill_ghosts(s, nullptr, buoy);
        });
        if (from != nullptr) { td.resume(from->time, from->steps, from->dt, from->anchor_time, from->anchor_steps); }
        while (!td.finished()) { td.advance(q); }
        lam.in.time = nullptr;
        return Ran{td.records(), td.time_exact(), td.total_steps(), double(uo.dt), td.anchor_time(), td.anchor_steps()};
    };
    const double dt1 = double(amrex::Real(0.1));
    const double dt2 = double(amrex::Real(0.05));
    const Ran straight = run(make(1.0, 0.1), nullptr);
    ASSERT_EQ(straight.records.size(), 10U);
    // The checkpoint of step 3: t = 3 dt as the straight run summed it, anchor (0, 0).
    const Ran chk3{{}, 3 * dt1, 3, dt1, 0, 0};
    const Ran resumed = run(make(1.0, 0.1), &chk3);
    ASSERT_EQ(resumed.records.size(), 7U);
    for (std::size_t n = 0; n < 7; ++n) {
        EXPECT_EQ(resumed.records[n].step, int(n) + 4);
        EXPECT_EQ(resumed.records[n].time, straight.records[n + 3].time) << "step " << n + 4;
    }
    EXPECT_EQ(resumed.time, straight.time);
    // Another fixed step after the checkpoint: t_ref + n dt.
    const Ran other = run(make(0.5, 0.05), &chk3);
    ASSERT_GE(other.records.size(), 2U);
    EXPECT_EQ(other.records[0].time, amrex::Real(3 * dt1 + 1 * dt2));
    EXPECT_EQ(other.records[1].time, amrex::Real(3 * dt1 + 2 * dt2));
    EXPECT_EQ(other.records.back().time, amrex::Real(0.5));
    // A clipped last step (dt 0.3 to 1.0: 0.3, 0.6, 0.9, 1.0), then a restart to 2.0: the next step
    // ends at 1.0 + dt, not at 5 dt = 1.5.
    const Ran clipped = run(make(1.0, 0.3), nullptr);
    ASSERT_EQ(clipped.records.size(), 4U);
    EXPECT_EQ(clipped.time, 1.0);
    EXPECT_EQ(clipped.anchor_time, 1.0);
    EXPECT_EQ(clipped.anchor_steps, 4);
    const Ran after = run(make(2.0, 0.3), &clipped);
    ASSERT_EQ(after.records.size(), 4U);
    const double dt3 = double(amrex::Real(0.3));
    EXPECT_EQ(after.records[0].step, 5);
    EXPECT_EQ(after.records[0].time, amrex::Real(1.0 + 1 * dt3));
    EXPECT_EQ(after.records[2].time, amrex::Real(1.0 + 3 * dt3));
    EXPECT_EQ(after.records.back().time, amrex::Real(2.0));
    // Chained restarts: dt 0.1 to 0.3, then dt 0.05 to 1.0 in one restart, or to 0.6 and on to 1.0
    // in two; the second leg of the chain has the single restart's times.
    const Ran first = run(make(0.3, 0.1), nullptr);
    EXPECT_EQ(first.anchor_steps, 0);   // the last step is a full one (to round-off): no new anchor
    const Ran single = run(make(1.0, 0.05), &first);
    const Ran leg1 = run(make(0.6, 0.05), &first);
    EXPECT_EQ(leg1.anchor_steps, first.steps);
    const Ran leg2 = run(make(1.0, 0.05), &leg1);
    ASSERT_EQ(leg1.records.size() + leg2.records.size(), single.records.size());
    for (std::size_t n = 0; n < leg2.records.size(); ++n) {
        const std::size_t m = n + leg1.records.size();
        EXPECT_EQ(leg2.records[n].step, single.records[m].step);
        EXPECT_EQ(leg2.records[n].time, single.records[m].time) << "step " << single.records[m].step;
    }
    EXPECT_EQ(leg2.time, single.time);
}
