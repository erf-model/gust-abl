// The nest faces of a nested patch (NestFaces.H), sub-phase 9.2 (beyond the plan).
//   UT-9.2.1  the net flux through the nest faces is zero to round-off after
//             the correction, for a non-uniform parent over a hill; the
//             correction is the same delta on every face point (the closed
//             form net / sum of the weights); nothing is removed with an
//             outflow face; the sums do not depend on the rank map
//   UT-9.2.2  the blended ground is the parent's on the nest faces (bit for
//             bit), the patch's own from the blend width inwards (bit for
//             bit) and the closed form between
//   UT-9.2.3  the relaxation zone's rate and source equal the closed form at
//             every point of u, v, w, theta and k: zero where the parent's
//             flow leaves, zero outside the zone, the two faces' rates added
//             in a corner
//   UT-9.2.4  the implicit diagonal equals the finite-difference derivative
//             of the source
//   UT-9.2.5  with width 0 the residual and the implicit side are left bit
//             for bit
//   also: the ghost fills and the Dirichlet rows through BoundaryConditions
//             (closed form, corners included; the homogeneous fill; the abort
//             without the parent's values); each input fault is refused
//             naming its key; the restart record.
// The parent's fields are a + b x + c y + d x y on each level (NestParent's
// interpolation reproduces them to round-off), with normal velocities that
// change sign along every face.
#include "BoundaryConditions.H"
#include "ImplicitOperator.H"
#include "Inputs.H"
#include "Mesh.H"
#include "NestFaces.H"
#include "NestParent.H"
#include "Residual.H"
#include "State.H"
#include "StateIO.H"
#include "ZLevels.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Utility.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double kTol = kDouble ? 1.0e-11 : 2.0e-5;
constexpr double kBig = 1.0e30;

struct ScopedDir
{
    std::string name;
    ~ScopedDir () {
        amrex::ParallelDescriptor::Barrier();
        if (amrex::ParallelDescriptor::IOProcessor()) { amrex::FileSystem::RemoveAll(name); }
        amrex::ParallelDescriptor::Barrier();
    }
};

struct NestRun
{
    gust::Inputs in;
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::DistributionMapping dm;
    gust::Mesh mesh;
    gust::State q;

    void make (bool reverse = false, const gust::Mesh::GroundEdit& edit = nullptr)
    {
        geom = gust::make_geometry(in.mesh);
        ba = gust::make_box_array(in.mesh);
        if (dm.empty()) { dm = amrex::DistributionMapping(ba); }   // else: the map set by the test
        if (reverse) {
            amrex::Vector<int> pmap = dm.ProcessorMap();
            const int np = amrex::ParallelDescriptor::NProcs();
            for (auto& p : pmap) { p = np - 1 - p; }
            dm = amrex::DistributionMapping(pmap);
        }
        q.define(ba, dm);
        mesh.define(ba, dm);
        mesh.build(in, geom, q.z_phys_nd, edit);
    }
};

// The parent: 20 x 16 x 8 cells of 100 m, stretched layers, uneven boxes.
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
        in.terrain.smoothing = 1;
    }
    return in;
}

// The patch over parent cells 3 to 14 in x and 2 to 11 in y: 24 x 20 cells of
// 50 m, uneven boxes; its own hill, centred in the patch (not the parent's).
gust::Inputs patch_inputs (const gust::Inputs& parent)
{
    gust::Inputs in = parent;
    in.mesh.n_cell = amrex::IntVect(24, 20, 8);
    in.mesh.prob_lo = {-700, 700, 0};
    in.mesh.prob_hi = {500, 1700, 800};
    in.mesh.max_grid_size = {9, 7};
    return in;
}

// Field c on level k: bilinear in x and y. u (c = 0) and v (c = 1) change
// sign over the patch's faces at the lower levels.
double bilinear (double x, double y, int k, int c)
{
    return (1.0 + c) + 0.3 * k + (0.002 + 1.0e-4 * c) * x - (0.0015 - 1.0e-4 * k) * y
         + 1.0e-6 * (1.0 + 0.1 * k + 0.05 * c) * x * y;
}

double pos (const amrex::Geometry& g, const amrex::IndexType& t, int i, int d)
{
    return static_cast<double>(g.ProbLo(d)) + (i + (t.nodeCentered(d) ? 0.0 : 0.5)) * static_cast<double>(g.CellSize(d));
}

// fn(i, j, k, value) at the points of region(mfi) of every box, on the host.
void on_host (const amrex::MultiFab& mf, int comp, const std::function<amrex::Box(const amrex::MFIter&)>& region,
              const std::function<void(int, int, int, double)>& fn)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = region(mfi) & mfi.fabbox();
        if (b.isEmpty()) { continue; }
        amrex::FArrayBox h(b, 1, amrex::The_Pinned_Arena());
        h.copy<amrex::RunOn::Device>(mf[mfi], b, comp, b, 0, 1);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { fn(i, j, k, static_cast<double>(h(amrex::IntVect(i, j, k)))); });
    }
}

amrex::Box valid (const amrex::MFIter& mfi) { return mfi.validbox(); }

// mf = fn(i, j, k) at its valid points, ghost points `ghost`.
void set_valid (amrex::MultiFab& mf, const std::function<double(int, int, int)>& fn, double ghost = kBig)
{
    mf.setVal(static_cast<amrex::Real>(ghost));
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(i, j, k)); });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

struct Sum
{
    double err{0}, sum{0};
    long n{0}, bad{0};
    void reduce ()
    {
        const auto comm = amrex::ParallelDescriptor::Communicator();
        amrex::ParallelAllReduce::Max(err, comm);
        amrex::ParallelAllReduce::Sum(sum, comm);
        amrex::ParallelAllReduce::Sum(n, comm);
        amrex::ParallelAllReduce::Sum(bad, comm);
    }
    void add (double d) { err = std::max(err, std::isfinite(d) ? std::abs(d) : kBig); ++n; }
};

std::vector<std::pair<amrex::MultiFab*, int>> fields6 (gust::State& q)
{
    return {{&q.u, 0}, {&q.v, 1}, {&q.w, 2}, {&q.p, 3}, {&q.theta, 4}, {&q.tke, 5}};
}

// The parent's checkpoint with the bilinear fields.
std::string write_parent (NestRun& p, const std::string& prefix)
{
    for (auto [mf, c] : fields6(p.q)) {
        const amrex::IndexType t = mf->ixType();
        const int cc = c;
        set_valid(*mf, [&] (int i, int j, int k) { return bilinear(pos(p.geom, t, i, 0), pos(p.geom, t, j, 1), k, cc); });
    }
    return gust::write_checkpoint(prefix, p.q, p.in.mesh, &p.mesh.zlevels());
}

// Faces xlo, ylo, zlo, xhi, yhi, zhi (amrex::Orientation's order).
gust::BoundaryInputs bc_inputs (const std::string& parent, gust::FaceType xlo, gust::FaceType xhi, gust::FaceType ylo,
                                gust::FaceType yhi)
{
    gust::BoundaryInputs b;
    b.face = {xlo, ylo, gust::FaceType::slip_wall, xhi, yhi, gust::FaceType::slip_wall};
    b.nest_parent = parent;
    return b;
}

gust::BoundaryInputs all_nest (const std::string& parent)
{
    const auto n = gust::FaceType::nest;
    return bc_inputs(parent, n, n, n, n);
}

bool is_nest (const gust::BoundaryInputs& b, int d, bool lo)
{
    return b.face[amrex::Orientation(d, lo ? amrex::Orientation::low : amrex::Orientation::high)] == gust::FaceType::nest;
}

// The parent's normal velocity of direction d at face point (index m along
// the face, level k) of face (d, lo) of the patch, after the correction:
// the closed form, with the rows below the ground and above the top holding
// the nearest level (NestParent.H).
double face_velocity (const NestRun& f, int d, bool lo, int m, int k, double delta)
{
    const amrex::Box& dom = f.geom.Domain();
    const int e = 1 - d;
    const double xn = lo ? f.geom.ProbLo(d) : f.geom.ProbHi(d);
    const double xt = static_cast<double>(f.geom.ProbLo(e)) + (m - dom.smallEnd(e) + 0.5) * static_cast<double>(f.geom.CellSize(e));
    const int kc = std::max(dom.smallEnd(2), std::min(dom.bigEnd(2), k));
    const double V = (d == 0) ? bilinear(xn, xt, kc, 0) : bilinear(xt, xn, kc, 1);
    return V + (lo ? delta : -delta);   // V - s delta, s = -1 on a low face
}

// The zone's rate at point (i, j, k) of type t: the closed form of NestFaces.H.
double expected_rate (const NestRun& f, const gust::BoundaryInputs& b, const amrex::IndexType& t, int i, int j, int k, int N,
                      double factor, double delta)
{
    const amrex::Box& dom = f.geom.Domain();
    double rate = 0;
    for (int d = 0; d < 2; ++d) {
        for (const bool lo : {true, false}) {
            if (!is_nest(b, d, lo)) { continue; }
            const int face = lo ? dom.smallEnd(d) : dom.bigEnd(d) + 1;
            const bool nodal = t.nodeCentered(d);
            const double s = std::abs(((d == 0) ? i : j) + (nodal ? 0.0 : 0.5) - face);
            if (nodal ? (s < 0.5 || s > N - 0.5) : (s > N)) { continue; }
            const int e = 1 - d;
            const int oe = t.nodeCentered(e) ? 1 : 0, oz = t.nodeCentered(2) ? 1 : 0;
            const int m = (d == 0) ? j : i;
            const double U = 0.25 * (face_velocity(f, d, lo, m, k, delta) + face_velocity(f, d, lo, m - oe, k, delta)
                                     + face_velocity(f, d, lo, m, k - oz, delta) + face_velocity(f, d, lo, m - oe, k - oz, delta));
            const double Uin = std::max(0.0, lo ? U : -U);
            rate += (1.0 - s / N) * factor * Uin / (N * static_cast<double>(f.geom.CellSize(d)));
        }
    }
    return rate;
}

std::string abort_message (const std::function<void()>& fn)
{
    try { fn(); } catch (const std::runtime_error& e) { return e.what(); }
    return "";
}

bool has (const std::vector<std::string>& err, const std::string& what)
{
    return std::any_of(err.begin(), err.end(), [&] (const std::string& e) { return e.find(what) != std::string::npos; });
}

// sum s W V over the nest faces' valid points, and sum W max(0, -s V).
void host_flux (const NestRun& f, const gust::BoundaryInputs& b, const gust::State& par, double& net, double& inflow, long& n)
{
    Sum s, sin;
    const amrex::MultiFab* vn[2] = {&par.u, &par.v};
    for (int d = 0; d < 2; ++d) {
        amrex::MultiFab W(vn[d]->boxArray(), vn[d]->DistributionMap(), 1, 0);
        amrex::MultiFab::Copy(W, f.mesh.flux_area(d), 0, 0, 1, 0);
        W.mult(f.geom.InvCellSize(d));
        amrex::MultiFab::Multiply(W, *vn[d], 0, 0, 1, 0);   // W V
        for (const bool lo : {true, false}) {
            if (!is_nest(b, d, lo)) { continue; }
            const int face = lo ? f.geom.Domain().smallEnd(d) : f.geom.Domain().bigEnd(d) + 1;
            const double sg = lo ? -1.0 : 1.0;
            on_host(W, 0, [&] (const amrex::MFIter& mfi) { return amrex::Box(mfi.validbox()).setSmall(d, face).setBig(d, face); },
                    [&] (int, int, int, double wv) { s.sum += sg * wv; sin.sum += std::max(0.0, -sg * wv); ++s.n; });
        }
    }
    s.reduce();
    sin.reduce();
    net = s.sum;
    inflow = sin.sum;
    n = s.n;
}

// The patch over a hill with its parent read: what every test below starts from.
struct Case
{
    NestRun p, f;
    ScopedDir dir;
    std::unique_ptr<gust::NestParent> parent;
    Case (const std::string& prefix, bool hill = true, bool reverse = false)
    {
        p.in = parent_inputs(hill);
        p.make();
        dir.name = write_parent(p, prefix);
        f.in = patch_inputs(p.in);
        f.make(reverse);
        parent = std::make_unique<gust::NestParent>(dir.name, f.in.mesh, f.mesh.zlevels(), f.ba, f.dm);
    }
};

// UT-9.2.1
TEST (NestFaces, NetFluxIsZeroAfterTheCorrection)
{
    Case c("ut_nestfaces_flux");
    const gust::BoundaryInputs b = all_nest(c.dir.name);
    const gust::NestFaces nf(*c.parent, b, c.f.geom, c.f.mesh);
    gust::State raw;
    raw.define(c.f.ba, c.f.dm);
    c.parent->interpolate(raw);

    // The interpolated parent has a net flux (not a field that balances by itself).
    double net0 = 0, in0 = 0, net1 = 0, in1 = 0;
    long n = 0;
    host_flux(c.f, b, raw, net0, in0, n);
    EXPECT_EQ(n, 2L * 8 * (24 + 20));
    EXPECT_GT(std::abs(net0), 1.0e-3 * in0);
    EXPECT_GT(in0, 0.0);
    EXPECT_NEAR(static_cast<double>(nf.flux().net), net0, (kDouble ? 1.0e-12 : 1.0e-5) * in0);
    EXPECT_NEAR(static_cast<double>(nf.flux().inflow), in0, (kDouble ? 1.0e-12 : 1.0e-5) * in0);
    EXPECT_TRUE(nf.flux().applied);
    const double delta = nf.flux().delta;
    EXPECT_NEAR(delta, net0 / static_cast<double>(nf.flux().weight), (kDouble ? 1.0e-13 : 1.0e-6) * std::abs(delta));

    // After: zero to round-off of the flux entering.
    host_flux(c.f, b, nf.parent_fields(), net1, in1, n);
    EXPECT_LT(std::abs(net1), (kDouble ? 1.0e-13 : 1.0e-5) * in0);

    // The same delta off every face point's outward velocity; nothing elsewhere.
    const amrex::MultiFab* cor[2] = {&nf.parent_fields().u, &nf.parent_fields().v};
    const amrex::MultiFab* rw[2] = {&raw.u, &raw.v};
    Sum on_face, off_face;
    for (int d = 0; d < 2; ++d) {
        amrex::MultiFab diff(cor[d]->boxArray(), c.f.dm, 1, 0);
        amrex::MultiFab::Copy(diff, *cor[d], 0, 0, 1, 0);
        amrex::MultiFab::Subtract(diff, *rw[d], 0, 0, 1, 0);
        const int flo = c.f.geom.Domain().smallEnd(d), fhi = c.f.geom.Domain().bigEnd(d) + 1;
        on_host(diff, 0, valid, [&] (int i, int j, int, double v) {
            const int m = (d == 0) ? i : j;
            if (m == flo) { on_face.add(v - delta); }
            else if (m == fhi) { on_face.add(v + delta); }
            else { off_face.add(v); }
        });
    }
    on_face.reduce();
    off_face.reduce();
    EXPECT_GT(on_face.n, 0);
    EXPECT_LT(on_face.err, (kDouble ? 1.0e-13 : 1.0e-5) * std::max(1.0, std::abs(delta)));
    EXPECT_EQ(off_face.err, 0.0);
    // theta, k and w are the interpolation itself.
    amrex::MultiFab d(c.f.ba, c.f.dm, 1, 0);
    amrex::MultiFab::Copy(d, nf.parent_fields().theta, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, raw.theta, 0, 0, 1, 0);
    EXPECT_EQ(d.norminf(), amrex::Real(0));
}

TEST (NestFaces, NothingIsRemovedWithAnOutflowFace)
{
    Case c("ut_nestfaces_outflow");
    const auto n = gust::FaceType::nest;
    const gust::BoundaryInputs b = bc_inputs(c.dir.name, n, gust::FaceType::outflow, n, n);
    const gust::NestFaces nf(*c.parent, b, c.f.geom, c.f.mesh);
    EXPECT_FALSE(nf.flux().applied);
    EXPECT_EQ(nf.flux().delta, amrex::Real(0));
    EXPECT_FALSE(nf.is_nest(0, false));
    EXPECT_TRUE(nf.is_nest(0, true));
    gust::State raw;
    raw.define(c.f.ba, c.f.dm);
    c.parent->interpolate(raw);
    amrex::MultiFab d(nf.parent_fields().u.boxArray(), c.f.dm, 1, 0);
    amrex::MultiFab::Copy(d, nf.parent_fields().u, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(d, raw.u, 0, 0, 1, 0);
    EXPECT_EQ(d.norminf(), amrex::Real(0));
    EXPECT_NE(nf.summary().find("not removed"), std::string::npos);
}

TEST (NestFaces, TheNetFluxDoesNotDependOnTheRankMap)
{
    Case a("ut_nestfaces_map_a", true, false), r("ut_nestfaces_map_r", true, true);
    const gust::NestFaces na(*a.parent, all_nest(a.dir.name), a.f.geom, a.f.mesh);
    const gust::NestFaces nr(*r.parent, all_nest(r.dir.name), r.f.geom, r.f.mesh);
    EXPECT_EQ(na.flux().net, nr.flux().net);
    EXPECT_EQ(na.flux().weight, nr.flux().weight);
    EXPECT_EQ(na.flux().delta, nr.flux().delta);
    EXPECT_NE(na.flux().delta, amrex::Real(0));
}

TEST (NestFaces, NoNestFaceAborts)
{
    Case c("ut_nestfaces_none", false);
    const auto in = gust::FaceType::inflow, out = gust::FaceType::outflow;
    const std::string msg = abort_message([&] { gust::NestFaces nf(*c.parent, bc_inputs(c.dir.name, in, out, in, out), c.f.geom, c.f.mesh); });
    EXPECT_NE(msg.find("no lateral face is a nest face"), std::string::npos) << msg;
}

// The log's lines on the net flux: its share of the flux entering, a warning
// above kNestNetWarning of it, and other words when nothing enters.
TEST (NestFaces, SummaryWarnsAboveTheNetFluxThreshold)
{
    {
        Case c("ut_nestfaces_sum_hill");   // the bilinear fields over a hill: a large net
        const gust::NestFaces nf(*c.parent, all_nest(c.dir.name), c.f.geom, c.f.mesh);
        ASSERT_GT(nf.flux().inflow, 0.0);
        ASSERT_GT(std::abs(nf.flux().net), 10 * gust::kNestNetWarning * nf.flux().inflow);
        const std::string s = nf.summary();
        EXPECT_NE(s.find("of the flux entering; removed before the pressure solve"), std::string::npos) << s;
        EXPECT_NE(s.find("WARNING: nest faces: the net flux removed is above 0.001 of the flux entering"), std::string::npos) << s;
        EXPECT_NE(s.find("terrain.nest_blend_width"), std::string::npos) << s;
        // With an outflow face nothing is removed: no warning.
        const gust::NestFaces no(*c.parent, bc_inputs(c.dir.name, gust::FaceType::nest, gust::FaceType::outflow, gust::FaceType::nest,
                                                      gust::FaceType::nest), c.f.geom, c.f.mesh);
        EXPECT_EQ(no.summary().find("WARNING"), std::string::npos) << no.summary();
        EXPECT_NE(no.summary().find("not removed (an outflow face carries it)"), std::string::npos) << no.summary();
    }
    // A uniform wind (speed 5 m/s along x, or none) over flat ground: no net.
    for (const double speed : {5.0, 0.0}) {
        NestRun p, f;
        p.in = parent_inputs(false);
        p.make();
        const double val[6] = {speed, 0.0, 0.0, 0.0, 300.0, 0.1};
        for (auto [mf, n] : fields6(p.q)) {
            const double v = val[n];
            set_valid(*mf, [v] (int, int, int) { return v; });
        }
        const ScopedDir dir{gust::write_checkpoint("ut_nestfaces_sum_flat", p.q, p.in.mesh, &p.mesh.zlevels())};
        f.in = patch_inputs(p.in);
        f.make();
        const gust::NestParent parent(dir.name, f.in.mesh, f.mesh.zlevels(), f.ba, f.dm);
        const gust::NestFaces nf(parent, all_nest(dir.name), f.geom, f.mesh);
        const std::string s = nf.summary();
        EXPECT_TRUE(nf.flux().applied);
        EXPECT_EQ(s.find("WARNING"), std::string::npos) << s;
        if (speed > 0) {
            EXPECT_GT(nf.flux().inflow, 0.0);
            EXPECT_LT(std::abs(nf.flux().net), 1.0e-6 * nf.flux().inflow);
            EXPECT_NE(s.find("of the flux entering; removed before the pressure solve"), std::string::npos) << s;
        } else {
            EXPECT_EQ(nf.flux().inflow, 0.0);
            EXPECT_EQ(nf.flux().delta, amrex::Real(0));
            EXPECT_NE(s.find("no flux enters; mean outward velocity 0.000e+00 m/s; removed before the pressure solve"),
                      std::string::npos) << s;
        }
    }
}

// The ghost fills and rows through BoundaryConditions.
TEST (NestFaces, GhostsAndFacePointsHoldTheParentsValues)
{
    Case c("ut_nestfaces_fill");
    const gust::BoundaryInputs b = all_nest(c.dir.name);
    const gust::NestFaces nf(*c.parent, b, c.f.geom, c.f.mesh);
    const double delta = nf.flux().delta;
    gust::BoundaryConditions bc(b, c.f.geom, c.f.mesh, gust::CoriolisOptions{});
    const gust::BuoyancyOptions buoy;
    for (auto [mf, n] : fields6(c.f.q)) { mf->setVal(amrex::Real(7 + n)); }

    // Without the parent's values the fill is refused.
    const std::string msg = abort_message([&] { bc.fill_ghosts(c.f.q, nullptr, buoy); });
    EXPECT_NE(msg.find("set_nest_faces"), std::string::npos) << msg;

    bc.set_nest_faces(&nf);
    gust::ResidualFields R;
    R.define(c.f.ba, c.f.dm);
    for (auto& r : R.fields()) { r.mf->setVal(amrex::Real(-3)); }
    bc.apply_rows(c.f.q, R);   // before the fill: q = 7 (u), 8 (v) on the face points
    bc.fill_ghosts(c.f.q, nullptr, buoy);

    const amrex::Box& dom = c.f.geom.Domain();
    // Outside the domain in x or y, on the valid levels: the closed form at
    // the point itself (corners included); the normal velocity's face points
    // and its ghost points outside them carry the correction. Inside: untouched.
    Sum ghost, inside, pfield;
    for (auto [mf, n] : fields6(c.f.q)) {
        const amrex::IndexType t = mf->ixType();
        const amrex::Box pts = amrex::convert(dom, t);
        const int fld = n;
        on_host(*mf, 0, [&] (const amrex::MFIter& mfi) { return amrex::grow(mfi.validbox(), amrex::IntVect(2, 2, 0)); },
                [&] (int i, int j, int k, double v) {
            const bool xin = i >= pts.smallEnd(0) && i <= pts.bigEnd(0), yin = j >= pts.smallEnd(1) && j <= pts.bigEnd(1);
            const bool xface = fld == 0 && (i == pts.smallEnd(0) || i == pts.bigEnd(0));
            const bool yface = fld == 1 && (j == pts.smallEnd(1) || j == pts.bigEnd(1));
            if (fld == 3) {
                // p': zero gradient, so 7 + 3 everywhere.
                pfield.add(v - 10.0);
            } else if (fld == 2 && (k == pts.smallEnd(2) || k == pts.bigEnd(2))) {
                // w on the ground and on the top: the walls' rows (the z fill), not the nest's.
            } else if (xin && yin && !xface && !yface) {
                inside.add(v - (7.0 + fld));
            } else {
                double e = bilinear(pos(c.f.geom, t, i, 0), pos(c.f.geom, t, j, 1), k, fld);
                if (fld == 0 && i <= pts.smallEnd(0)) { e += delta; }
                if (fld == 0 && i >= pts.bigEnd(0)) { e -= delta; }
                if (fld == 1 && j <= pts.smallEnd(1)) { e += delta; }
                if (fld == 1 && j >= pts.bigEnd(1)) { e -= delta; }
                ghost.add(v - e);
            }
        });
    }
    ghost.reduce();
    inside.reduce();
    pfield.reduce();
    EXPECT_GT(ghost.n, 0);
    EXPECT_LT(ghost.err, kTol);
    EXPECT_EQ(inside.err, 0.0);
    EXPECT_EQ(pfield.err, 0.0);

    // The rows: R = q - V on the faces' normal-velocity points, untouched elsewhere.
    Sum rows, others;
    const amrex::MultiFab* rf[2] = {&R.u, &R.v};
    for (int d = 0; d < 2; ++d) {
        const int flo = dom.smallEnd(d), fhi = dom.bigEnd(d) + 1;
        on_host(*rf[d], 0, valid, [&] (int i, int j, int k, double v) {
            const int m = (d == 0) ? i : j, mt = (d == 0) ? j : i;
            if (m == flo || m == fhi) { rows.add(v - ((7.0 + d) - face_velocity(c.f, d, m == flo, mt, k, delta))); }
            else { others.add(v + 3.0); }
        });
    }
    rows.reduce();
    others.reduce();
    EXPECT_EQ(rows.n, 2L * 8 * (24 + 20));
    EXPECT_LT(rows.err, kTol);
    EXPECT_EQ(others.err, 0.0);

    // The identity rows of the implicit side are on the same points.
    amrex::MultiFab mask;
    bc.dirichlet_rows(gust::Equation::u, gust::Loc::xface, c.f.mesh, mask);
    EXPECT_EQ(static_cast<long>(std::lround(static_cast<double>(mask.sum(0)))), 2L * 8 * 20);

    // The homogeneous fill: zero in the ghosts and on the face points.
    for (auto [mf, n] : fields6(c.f.q)) { mf->setVal(amrex::Real(7 + n)); }
    bc.fill_ghosts(c.f.q, nullptr, buoy, true);
    Sum hom;
    for (auto [mf, n] : fields6(c.f.q)) {
        if (n == 3) { continue; }
        const amrex::Box pts = amrex::convert(dom, mf->ixType());
        const int fld = n;
        on_host(*mf, 0, [&] (const amrex::MFIter& mfi) { return amrex::grow(mfi.validbox(), amrex::IntVect(2, 2, 0)); },
                [&] (int i, int j, int, double v) {
            const bool xin = i > pts.smallEnd(0) - (fld == 0 ? 0 : 1) && i < pts.bigEnd(0) + (fld == 0 ? 0 : 1);
            const bool yin = j > pts.smallEnd(1) - (fld == 1 ? 0 : 1) && j < pts.bigEnd(1) + (fld == 1 ? 0 : 1);
            if (!(xin && yin)) { hom.add(v); }
        });
    }
    hom.reduce();
    EXPECT_GT(hom.n, 0);
    EXPECT_EQ(hom.err, 0.0);
}

// UT-9.2.2
TEST (NestFaces, BlendedGroundIsTheParentsOnTheFaceAndThePatchsInside)
{
    NestRun p;
    p.in = parent_inputs(true);
    p.make();
    const ScopedDir dir{write_parent(p, "ut_nestfaces_blend")};
    NestRun own;   // the patch's own ground
    own.in = patch_inputs(p.in);
    own.make();
    const double W = 200;
    for (const bool all : {true, false}) {
        const auto n = gust::FaceType::nest;
        // all: four nest faces. Else xlo and yhi only (the distance is to the nest faces alone).
        const gust::BoundaryInputs b = all ? all_nest(dir.name) : bc_inputs(dir.name, n, gust::FaceType::outflow, gust::FaceType::inflow, n);
        NestRun f;
        f.in = own.in;
        f.in.terrain.nest_blend_width = W;
        const amrex::Geometry g = gust::make_geometry(f.in.mesh);
        const amrex::BoxArray ba = gust::make_box_array(f.in.mesh);
        const amrex::DistributionMapping dm = own.dm;   // one map for the three grounds compared
        const gust::NestParent parent(dir.name, f.in.mesh, gust::make_zlevels(f.in.mesh, g), ba, dm, true);
        f.dm = dm;   // the map the parent was read for
        f.make(false, [&] (amrex::MultiFab& z) { gust::blend_nest_ground(parent, b, W, g, z); });
        EXPECT_GT(static_cast<double>(f.mesh.detJ.min(0)), 0.0);

        amrex::MultiFab hp(f.q.z_phys_nd.boxArray(), f.dm, 1, amrex::IntVect(1, 1, 0));
        parent.interpolate_ground(hp);
        // Per node of the ground: the blended height, the parent's and the patch's own.
        amrex::MultiFab three(f.q.z_phys_nd.boxArray(), f.dm, 3, 0);
        amrex::MultiFab::Copy(three, f.q.z_phys_nd, 0, 0, 1, 0);
        amrex::MultiFab::Copy(three, hp, 0, 1, 1, 0);
        amrex::MultiFab::Copy(three, own.q.z_phys_nd, 0, 2, 1, 0);
        const amrex::Box nodes = amrex::surroundingNodes(g.Domain());
        Sum face, inner, between, differ;
        for (amrex::MFIter mfi(three, false); mfi.isValid(); ++mfi) {
            amrex::Box bx = mfi.validbox();
            if (bx.smallEnd(2) > 0) { continue; }
            bx.setBig(2, 0);
            amrex::FArrayBox h(bx, 3, amrex::The_Pinned_Arena());
            h.copy<amrex::RunOn::Device>(three[mfi], bx, 0, bx, 0, 3);
            amrex::Gpu::streamSynchronize();
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                const double z = h(iv, 0), zp = h(iv, 1), zo = h(iv, 2);
                const double dx = g.CellSize(0), dy = g.CellSize(1);
                double dist = 1.0e30;
                if (is_nest(b, 0, true)) { dist = std::min(dist, (i - nodes.smallEnd(0)) * dx); }
                if (is_nest(b, 0, false)) { dist = std::min(dist, (nodes.bigEnd(0) - i) * dx); }
                if (is_nest(b, 1, true)) { dist = std::min(dist, (j - nodes.smallEnd(1)) * dy); }
                if (is_nest(b, 1, false)) { dist = std::min(dist, (nodes.bigEnd(1) - j) * dy); }
                if (dist == 0) {
                    if (z != zp) { ++face.bad; }
                    face.add(0);
                    differ.err = std::max(differ.err, std::abs(zp - zo));
                } else if (dist >= W) {
                    if (z != zo) { ++inner.bad; }
                    inner.add(0);
                } else {
                    const double t = dist / W, s = t * t * (3 - 2 * t);
                    between.add(z - ((1 - s) * zp + s * zo));
                }
            });
        }
        face.reduce();
        inner.reduce();
        between.reduce();
        differ.reduce();
        EXPECT_GE(face.n, all ? 2L * (25 + 21) - 4 : 25L + 21 - 1);   // the boxes share nodes
        EXPECT_EQ(face.bad, 0);
        EXPECT_GT(inner.n, 0);
        EXPECT_EQ(inner.bad, 0);
        EXPECT_GT(between.n, 0);
        EXPECT_LT(between.err, kDouble ? 1.0e-11 : 1.0e-3);
        EXPECT_GT(differ.err, 1.0) << "the parent's ground and the patch's differ on the face by more than 1 m";
    }
}

// UT-9.2.3 to UT-9.2.5 share this: the zone's rate, source and diagonal.
struct Zone
{
    Case c{"ut_nestfaces_zone"};
    gust::BoundaryInputs b;
    std::unique_ptr<gust::NestFaces> nf;
    Zone (int N, double factor, bool all = true)
    {
        const auto n = gust::FaceType::nest;
        b = all ? all_nest(c.dir.name) : bc_inputs(c.dir.name, n, gust::FaceType::outflow, n, n);
        b.nest_relax_width = N;
        b.nest_relax_factor = static_cast<amrex::Real>(factor);
        nf = std::make_unique<gust::NestFaces>(*c.parent, b, c.f.geom, c.f.mesh);
        // A state unlike the parent's.
        for (auto [mf, m] : fields6(c.f.q)) {
            const int fld = m;
            set_valid(*mf, [&] (int i, int j, int k) { return 3.0 + 0.5 * fld + 0.07 * i - 0.05 * j + 0.11 * k; });
        }
    }
};

TEST (NestFaces, ZoneRateAndSourceEqualTheClosedForm)
{
    for (const bool all : {true, false}) {
        const int N = 4;
        const double factor = 1.5;
        Zone z(N, factor, all);
        ASSERT_TRUE(z.nf->relax_on());
        const double delta = z.nf->flux().delta;
        gust::ResidualFields R;
        R.define(z.c.f.ba, z.c.f.dm);
        for (auto& r : R.fields()) { r.mf->setVal(amrex::Real(0)); }
        z.nf->add_relaxation(z.c.f.q, R);

        const gust::Loc locs[5] = {gust::Loc::xface, gust::Loc::yface, gust::Loc::zface, gust::Loc::cell, gust::Loc::cell};
        amrex::MultiFab* qf[5] = {&z.c.f.q.u, &z.c.f.q.v, &z.c.f.q.w, &z.c.f.q.theta, &z.c.f.q.tke};
        amrex::MultiFab* rf[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
        const int fld[5] = {0, 1, 2, 4, 5};
        Sum rate, src, pos_rate, zero_in_zone, at_face;
        for (int n = 0; n < 5; ++n) {
            const amrex::MultiFab& rt = z.nf->rate(locs[n]);
            const amrex::IndexType t = rt.ixType();
            ASSERT_EQ(rt.boxArray(), qf[n]->boxArray());
            // rate, q and R side by side.
            amrex::MultiFab three(rt.boxArray(), z.c.f.dm, 3, 0);
            amrex::MultiFab::Copy(three, rt, 0, 0, 1, 0);
            amrex::MultiFab::Copy(three, *qf[n], 0, 1, 1, 0);
            amrex::MultiFab::Copy(three, *rf[n], 0, 2, 1, 0);
            const amrex::Box& dom = z.c.f.geom.Domain();
            for (amrex::MFIter mfi(three, false); mfi.isValid(); ++mfi) {
                const amrex::Box bx = mfi.validbox();
                amrex::FArrayBox h(bx, 3, amrex::The_Pinned_Arena());
                h.copy<amrex::RunOn::Device>(three[mfi], bx, 0, bx, 0, 3);
                amrex::Gpu::streamSynchronize();
                amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                    const amrex::IntVect iv(i, j, k);
                    const double e = expected_rate(z.c.f, z.b, t, i, j, k, N, factor, delta);
                    rate.add(h(iv, 0) - e);
                    const double par = bilinear(pos(z.c.f.geom, t, i, 0), pos(z.c.f.geom, t, j, 1), k, fld[n]);
                    // The normal velocity's face points are the Dirichlet rows (apply_rows
                    // replaces what is added there): the source is not compared on them.
                    const bool on_x = n == 0 && (i == dom.smallEnd(0) || i == dom.bigEnd(0) + 1);
                    const bool on_y = n == 1 && (j == dom.smallEnd(1) || j == dom.bigEnd(1) + 1);
                    if (on_x || on_y) { ++at_face.n; return; }
                    src.add(h(iv, 2) - e * (h(iv, 1) - par));
                    if (e > 0) { ++pos_rate.n; }
                    // In the zone of xlo alone with the parent's flow leaving: no source.
                    if (i - dom.smallEnd(0) < N - 1 && j - dom.smallEnd(1) > N && dom.bigEnd(1) - j > N && e == 0) { ++zero_in_zone.n; }
                });
            }
        }
        for (Sum* s : {&rate, &src, &pos_rate, &zero_in_zone, &at_face}) { s->reduce(); }
        EXPECT_GT(pos_rate.n, 1000) << all;
        EXPECT_GT(zero_in_zone.n, 10) << "part of xlo is an outflow part" << all;
        EXPECT_LT(rate.err, kDouble ? 1.0e-14 : 1.0e-7) << all;
        EXPECT_LT(src.err, kDouble ? 1.0e-12 : 1.0e-5) << all;
        EXPECT_GT(at_face.n, 0);
    }
}

// UT-9.2.4
TEST (NestFaces, DiagonalIsTheSourcesDerivative)
{
    Zone z(3, 0.8);
    const double eps = kDouble ? 1.0e-3 : 0.25;
    gust::ResidualFields R0, R1;
    R0.define(z.c.f.ba, z.c.f.dm);
    R1.define(z.c.f.ba, z.c.f.dm);
    for (auto& r : R0.fields()) { r.mf->setVal(amrex::Real(0)); }
    for (auto& r : R1.fields()) { r.mf->setVal(amrex::Real(0)); }
    z.nf->add_relaxation(z.c.f.q, R0);
    amrex::MultiFab* qf[5] = {&z.c.f.q.u, &z.c.f.q.v, &z.c.f.q.w, &z.c.f.q.theta, &z.c.f.q.tke};
    for (auto* mf : qf) { mf->plus(static_cast<amrex::Real>(eps), 0, 1, 0); }
    z.nf->add_relaxation(z.c.f.q, R1);
    amrex::MultiFab* r0[5] = {&R0.u, &R0.v, &R0.w, &R0.theta, &R0.tke};
    amrex::MultiFab* r1[5] = {&R1.u, &R1.v, &R1.w, &R1.theta, &R1.tke};
    const gust::Equation eqs[5] = {gust::Equation::u, gust::Equation::v, gust::Equation::w, gust::Equation::theta, gust::Equation::tke};
    const gust::Loc locs[5] = {gust::Loc::xface, gust::Loc::yface, gust::Loc::zface, gust::Loc::cell, gust::Loc::cell};
    for (int n = 0; n < 5; ++n) {
        amrex::MultiFab A(r0[n]->boxArray(), z.c.f.dm, gust::n_implicit_coefs, 0);
        A.setVal(amrex::Real(0));
        z.nf->add_relaxation_to_implicit(eqs[n], locs[n], A);
        EXPECT_GT(A.norminf(gust::coef_P), amrex::Real(0)) << n;
        for (int cc = 1; cc < gust::n_implicit_coefs; ++cc) { EXPECT_EQ(A.norminf(cc), amrex::Real(0)) << n; }
        amrex::MultiFab fd(r0[n]->boxArray(), z.c.f.dm, 1, 0);
        amrex::MultiFab::Copy(fd, *r1[n], 0, 0, 1, 0);
        amrex::MultiFab::Subtract(fd, *r0[n], 0, 0, 1, 0);
        fd.mult(static_cast<amrex::Real>(1.0 / eps));
        amrex::MultiFab::Subtract(fd, A, gust::coef_P, 0, 1, 0);
        EXPECT_LT(static_cast<double>(fd.norminf()), (kDouble ? 1.0e-11 : 1.0e-4)) << n;
    }

    // The residual and the implicit side reach the zone through the boundary
    // conditions' damping hooks (no other damping here): the same bits.
    gust::BoundaryConditions bc(z.b, z.c.f.geom, z.c.f.mesh, gust::CoriolisOptions{});
    bc.set_nest_faces(z.nf.get());
    gust::ResidualFields Rb;
    Rb.define(z.c.f.ba, z.c.f.dm);
    for (auto& r : Rb.fields()) { r.mf->setVal(amrex::Real(0)); }
    bc.add_damping(z.c.f.q, Rb);
    amrex::MultiFab* rb[5] = {&Rb.u, &Rb.v, &Rb.w, &Rb.theta, &Rb.tke};
    for (int n = 0; n < 5; ++n) {
        EXPECT_GT(rb[n]->norminf(), amrex::Real(0)) << n;
        amrex::MultiFab::Subtract(*rb[n], *r1[n], 0, 0, 1, 0);
        EXPECT_EQ(rb[n]->norminf(), amrex::Real(0)) << n;
        amrex::MultiFab A(r0[n]->boxArray(), z.c.f.dm, gust::n_implicit_coefs, 0), Ab(r0[n]->boxArray(), z.c.f.dm, gust::n_implicit_coefs, 0);
        A.setVal(amrex::Real(0));
        Ab.setVal(amrex::Real(0));
        z.nf->add_relaxation_to_implicit(eqs[n], locs[n], A);
        bc.add_damping_to_implicit(eqs[n], locs[n], z.c.f.q.z_phys_nd, Ab);
        amrex::MultiFab::Subtract(Ab, A, 0, 0, gust::n_implicit_coefs, 0);
        EXPECT_EQ(Ab.norminf(gust::coef_P), amrex::Real(0)) << n;
    }
}

// UT-9.2.5
TEST (NestFaces, WidthZeroLeavesTheResidualAndTheDiagonal)
{
    Zone z(0, 1.0);
    EXPECT_FALSE(z.nf->relax_on());
    EXPECT_TRUE(z.nf->rate(gust::Loc::cell).empty());
    gust::ResidualFields R;
    R.define(z.c.f.ba, z.c.f.dm);
    int n = 0;
    for (auto& r : R.fields()) {
        const int m = n++;
        set_valid(*r.mf, [&] (int i, int j, int k) { return 0.1 * m + 0.013 * i - 0.007 * j + 0.003 * k; }, 0.0);
    }
    gust::ResidualFields R0;
    R0.define(z.c.f.ba, z.c.f.dm);
    {
        auto a = R.fields();
        auto b0 = R0.fields();
        for (int m = 0; m < static_cast<int>(a.size()); ++m) { amrex::MultiFab::Copy(*b0[m].mf, *a[m].mf, 0, 0, 1, 0); }
    }
    z.nf->add_relaxation(z.c.f.q, R);
    {
        auto a = R.fields();
        auto b0 = R0.fields();
        for (int m = 0; m < static_cast<int>(a.size()); ++m) {
            amrex::MultiFab::Subtract(*b0[m].mf, *a[m].mf, 0, 0, 1, 0);
            EXPECT_EQ(b0[m].mf->norminf(), amrex::Real(0)) << a[m].name;
        }
    }
    amrex::MultiFab A(z.c.f.ba, z.c.f.dm, gust::n_implicit_coefs, 0);
    A.setVal(amrex::Real(0.25));
    z.nf->add_relaxation_to_implicit(gust::Equation::theta, gust::Loc::cell, A);
    for (int cc = 0; cc < gust::n_implicit_coefs; ++cc) {
        EXPECT_EQ(A.min(cc), amrex::Real(0.25));
        EXPECT_EQ(A.max(cc), amrex::Real(0.25));
    }
}

// Each input fault is refused, naming its key.
TEST (NestFaces, InputFaultsAreRefusedNamingTheKey)
{
    gust::Inputs in = patch_inputs(parent_inputs(false));
    const amrex::Geometry g = gust::make_geometry(in.mesh);
    const gust::CoriolisOptions cor;
    const gust::BoundaryInputs good = all_nest("chk");
    EXPECT_TRUE(gust::check_boundary_inputs(good, g, cor).empty());
    EXPECT_TRUE(gust::has_nest_face(good));
    EXPECT_STREQ(gust::face_type_name(gust::FaceType::nest), "nest");
    EXPECT_EQ(gust::face_type_from_name("nest", "bc.xlo"), gust::FaceType::nest);

    auto errs = [&] (const std::function<void(gust::BoundaryInputs&)>& make) {
        gust::BoundaryInputs b = good;
        make(b);
        return gust::check_boundary_inputs(b, g, cor);
    };
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.face[5] = gust::FaceType::nest; }), "bc.zhi = nest: only the lateral faces"));
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.face[2] = gust::FaceType::nest; }), "bc.zlo = nest: only the lateral faces"));
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.face[0] = gust::FaceType::inflow; }), "a nest face with an inflow face needs an outflow face"));
    EXPECT_TRUE(errs([] (gust::BoundaryInputs& b) { b.face[0] = gust::FaceType::inflow; b.face[3] = gust::FaceType::outflow; b.profile_needed = true; }).empty());
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.nest_relax_width = -1; }), "bc.nest_relax_width = -1: must be >= 0"));
    // 24 x 20 cells: below half the cells across every nest face.
    EXPECT_TRUE(errs([] (gust::BoundaryInputs& b) { b.nest_relax_width = 9; }).empty());
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.nest_relax_width = 10; }), "bc.nest_relax_width = 10: must be below half the patch's cells in y (20 cells)"));
    EXPECT_FALSE(has(errs([] (gust::BoundaryInputs& b) { b.nest_relax_width = 10; }), "cells in x"));
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.nest_relax_width = 12; }), "cells in x (24 cells)"));
    // y outflow faces: only x counts.
    EXPECT_TRUE(errs([] (gust::BoundaryInputs& b) { b.face[1] = b.face[4] = gust::FaceType::outflow; b.nest_relax_width = 11; }).empty());
    // A lateral wall beside a nest face: the parent's flow crosses it.
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.face[1] = gust::FaceType::slip_wall; }),
                    "bc.ylo = slip_wall is not supported with a nest face"));
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.face[3] = gust::FaceType::no_slip_wall; }),
                    "bc.xhi = no_slip_wall is not supported with a nest face"));
    EXPECT_FALSE(has(errs([] (gust::BoundaryInputs& b) { b.face[1] = gust::FaceType::slip_wall; }), "bc.zhi"));
    // The damping layers relax towards the bc.inflow_* profile: with a nest
    // face and no inflow face the Rayleigh layer damps w and theta only, and
    // there is no sponge.
    auto rayleigh = [] (gust::BoundaryInputs& b, std::array<bool, 5> fields) {
        b.rayleigh_zdamp = 500;
        b.rayleigh_dampcoef = amrex::Real(0.1);
        b.rayleigh_fields = fields;
    };
    const char* const ray = "bc.rayleigh_fields: with a nest face and no inflow face the Rayleigh layer may damp w and theta only";
    EXPECT_TRUE(has(errs([&] (gust::BoundaryInputs& b) { rayleigh(b, {true, true, true, true, true}); }), ray));
    for (const int f : {0, 1, 4}) {
        std::array<bool, 5> one{false, false, true, false, false};
        one[static_cast<std::size_t>(f)] = true;
        EXPECT_TRUE(has(errs([&] (gust::BoundaryInputs& b) { rayleigh(b, one); }), ray)) << f;
    }
    EXPECT_TRUE(errs([&] (gust::BoundaryInputs& b) { rayleigh(b, {false, false, true, false, false}); }).empty());
    EXPECT_TRUE(errs([&] (gust::BoundaryInputs& b) { rayleigh(b, {false, false, true, true, false}); }).empty());
    // With an inflow face the profile is the run's own choice.
    EXPECT_FALSE(has(errs([&] (gust::BoundaryInputs& b) {
        b.face[0] = gust::FaceType::inflow; b.face[3] = gust::FaceType::outflow; rayleigh(b, {true, true, true, true, true}); }), ray));
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) {
        b.face[3] = gust::FaceType::outflow; b.sponge_width = 100; b.sponge_strength = amrex::Real(0.1); }),
                    "bc.sponge_width > 0 is not supported with a nest face and no inflow face"));
    for (const double bad : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
        EXPECT_TRUE(has(errs([bad] (gust::BoundaryInputs& b) { b.nest_relax_width = 2; b.nest_relax_factor = static_cast<amrex::Real>(bad); }),
                        "bc.nest_relax_factor = ")) << bad;
    }
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) { b.given.emplace_back("nest_relax_factor"); }),
                    "bc.nest_relax_factor has no effect without bc.nest_relax_width > 0"));
    EXPECT_TRUE(has(errs([] (gust::BoundaryInputs& b) {
        b.face[0] = gust::FaceType::inflow; b.face[3] = gust::FaceType::outflow; b.face[1] = b.face[4] = gust::FaceType::slip_wall;
        b.profile_needed = true; b.given.emplace_back("nest_relax_width"); }), "bc.nest_relax_width has no effect without a nest face"));

    // The checks of the whole run.
    auto run = [&] (const gust::BoundaryInputs& b, const gust::Inputs& i, bool guess, bool seq) { return gust::check_nest_run(b, i, guess, seq); };
    EXPECT_TRUE(run(good, in, true, false).empty());
    gust::BoundaryInputs no_parent = good;
    no_parent.nest_parent.clear();
    EXPECT_TRUE(has(run(no_parent, in, true, false), "bc.xlo = nest needs bc.nest_parent"));
    EXPECT_TRUE(has(run(no_parent, in, true, false), "init.first_guess = nest needs bc.nest_parent"));
    EXPECT_TRUE(has(run(good, in, false, false), "bc.xlo = nest needs init.first_guess = nest"));
    gust::Inputs restart = in;
    restart.io.restart = "chk00010";
    EXPECT_TRUE(run(good, restart, false, false).empty());
    EXPECT_TRUE(has(run(good, in, true, true), "bc.xlo = nest is not supported with steady.sequencing = 1"));
    gust::BoundaryInputs none = bc_inputs("chk", gust::FaceType::inflow, gust::FaceType::outflow, gust::FaceType::inflow, gust::FaceType::outflow);
    EXPECT_TRUE(run(none, in, true, false).empty());
    EXPECT_TRUE(has(run(none, in, false, false), "bc.nest_parent has no effect without init.first_guess = nest"));
    gust::BoundaryInputs yhi_only = none;
    yhi_only.face[4] = gust::FaceType::nest;
    yhi_only.nest_parent.clear();
    EXPECT_TRUE(has(run(yhi_only, in, true, false), "bc.yhi = nest needs bc.nest_parent"));

    // terrain.nest_blend_width: the patch is 1200 m x 1000 m.
    auto blend = [&] (double W, const gust::BoundaryInputs& b, bool given = true) {
        gust::Inputs i = in;
        i.terrain.nest_blend_width = static_cast<amrex::Real>(W);
        if (given) { i.terrain.given.emplace_back("nest_blend_width"); }
        return run(b, i, true, false);
    };
    EXPECT_TRUE(blend(499, good).empty());
    EXPECT_TRUE(has(blend(500, good), "terrain.nest_blend_width = 500: must be below half the patch's width in y (500 m)"));
    EXPECT_TRUE(has(blend(-1, good), "terrain.nest_blend_width = "));
    EXPECT_TRUE(has(blend(std::numeric_limits<double>::quiet_NaN(), good), "terrain.nest_blend_width = "));
    EXPECT_TRUE(has(blend(std::numeric_limits<double>::infinity(), good), "terrain.nest_blend_width = "));
    EXPECT_TRUE(has(blend(100, none), "terrain.nest_blend_width has no effect without a nest face"));
    EXPECT_TRUE(has(blend(0, none), "terrain.nest_blend_width has no effect without a nest face"));
    EXPECT_TRUE(blend(0, none, false).empty());
    // x faces only: half the width in x (600 m) is the bound.
    gust::BoundaryInputs xonly = good;
    xonly.face[1] = xonly.face[4] = gust::FaceType::outflow;
    EXPECT_TRUE(blend(550, xonly).empty());
    EXPECT_TRUE(has(blend(600, xonly), "in x (600 m)"));
}

// The restart record.
TEST (NestFaces, RestartRecordIsCheckedAgainstTheInputs)
{
    Case c("ut_nestfaces_rst", false);
    gust::BoundaryInputs b = all_nest(c.dir.name);
    const std::string h0 = gust::NestFaces::hash_inputs(*c.parent, b, 0);
    EXPECT_EQ(h0.size(), 64U);
    EXPECT_EQ(h0, gust::NestFaces::hash_inputs(*c.parent, b, 0));
    // Every setting of the record changes it; the parent's path does not.
    gust::BoundaryInputs b2 = b;
    b2.nest_relax_width = 3;
    EXPECT_NE(h0, gust::NestFaces::hash_inputs(*c.parent, b2, 0));
    gust::BoundaryInputs b3 = b2;
    b3.nest_relax_factor = 2;
    EXPECT_NE(gust::NestFaces::hash_inputs(*c.parent, b2, 0), gust::NestFaces::hash_inputs(*c.parent, b3, 0));
    EXPECT_NE(h0, gust::NestFaces::hash_inputs(*c.parent, b, 150));
    gust::BoundaryInputs b4 = b;
    b4.face[3] = gust::FaceType::outflow;
    EXPECT_NE(h0, gust::NestFaces::hash_inputs(*c.parent, b4, 0));
    // A face that is not nest counts too (an outflow face decides whether the net flux is removed).
    gust::BoundaryInputs b6 = b4;
    b6.face[3] = gust::FaceType::inflow;
    EXPECT_NE(gust::NestFaces::hash_inputs(*c.parent, b4, 0), gust::NestFaces::hash_inputs(*c.parent, b6, 0));
    gust::BoundaryInputs b5 = b;
    b5.nest_parent = "elsewhere";
    EXPECT_EQ(h0, gust::NestFaces::hash_inputs(*c.parent, b5, 0));
    // The parent: read again it gives the same record; another step, other
    // fields at the same step, or the patch one parent cell further, another.
    {
        const gust::NestParent again(c.dir.name, c.f.in.mesh, c.f.mesh.zlevels(), c.f.ba, c.f.dm);
        EXPECT_EQ(c.parent->content_digest().size(), 64U);
        EXPECT_EQ(again.content_digest(), c.parent->content_digest());
        EXPECT_EQ(h0, gust::NestFaces::hash_inputs(again, b, 0));

        c.p.q.step = 5;
        const ScopedDir later{gust::write_checkpoint("ut_nestfaces_rst_later", c.p.q, c.p.in.mesh, &c.p.mesh.zlevels())};
        const gust::NestParent p_later(later.name, c.f.in.mesh, c.f.mesh.zlevels(), c.f.ba, c.f.dm);
        EXPECT_EQ(p_later.parent_step(), 5);
        EXPECT_NE(h0, gust::NestFaces::hash_inputs(p_later, b, 0));

        c.p.q.step = c.parent->parent_step();
        c.p.q.u.plus(amrex::Real(0.5), 0, 1, 0);   // another flow, same step and cells
        const ScopedDir flow{gust::write_checkpoint("ut_nestfaces_rst_flow", c.p.q, c.p.in.mesh, &c.p.mesh.zlevels())};
        const gust::NestParent p_flow(flow.name, c.f.in.mesh, c.f.mesh.zlevels(), c.f.ba, c.f.dm);
        EXPECT_EQ(p_flow.parent_step(), c.parent->parent_step());
        EXPECT_EQ(p_flow.parent_cells(), c.parent->parent_cells());
        EXPECT_EQ(p_flow.offset(), c.parent->offset());
        EXPECT_NE(p_flow.content_digest(), c.parent->content_digest());
        EXPECT_NE(h0, gust::NestFaces::hash_inputs(p_flow, b, 0));

        gust::MeshInputs moved = c.f.in.mesh;   // one parent cell (100 m) further in x
        moved.prob_lo[0] += 100;
        moved.prob_hi[0] += 100;
        const gust::NestParent p_moved(c.dir.name, moved, c.f.mesh.zlevels(), c.f.ba, c.f.dm);
        EXPECT_EQ(p_moved.offset()[0], c.parent->offset()[0] + 1);
        EXPECT_EQ(p_moved.content_digest(), c.parent->content_digest());
        EXPECT_NE(h0, gust::NestFaces::hash_inputs(p_moved, b, 0));
    }

    // A checkpoint with the record, and one without (the parent's own).
    const ScopedDir with{"ut_nestfaces_rst_with"};
    amrex::UtilCreateCleanDirectory(with.name, true);
    gust::NestFaces::write_checkpoint(with.name, h0);
    gust::NestFaces::check_restart(with.name, h0);
    const std::string other = abort_message([&] { gust::NestFaces::check_restart(with.name, gust::NestFaces::hash_inputs(*c.parent, b2, 0)); });
    EXPECT_NE(other.find("io.restart = " + with.name), std::string::npos) << other;
    EXPECT_NE(other.find("bc.nest_relax_width"), std::string::npos) << other;
    const std::string dropped = abort_message([&] { gust::NestFaces::check_restart(with.name, ""); });
    EXPECT_NE(dropped.find("written with nest faces"), std::string::npos) << dropped;
    gust::NestFaces::check_restart(c.dir.name, "");
    const std::string added = abort_message([&] { gust::NestFaces::check_restart(c.dir.name, h0); });
    EXPECT_NE(added.find("written without nest faces"), std::string::npos) << added;
    // The record's text is pinned: a checkpoint written today is read by a later version.
    const ScopedDir fixed{"ut_nestfaces_rst_fixed"};
    amrex::UtilCreateCleanDirectory(fixed.name, true);
    auto put = [&] (const std::string& text) {
        if (amrex::ParallelDescriptor::IOProcessor()) { std::ofstream(fixed.name + "/NestFaces") << text; }
        amrex::ParallelDescriptor::Barrier();
    };
    put("gust_nest_faces 1\nsha256 " + h0 + "\n");
    gust::NestFaces::check_restart(fixed.name, h0);
    for (const std::string& bad : {"gust_nest_faces 2\nsha256 " + h0 + "\n", "gust_faces 1\nsha256 " + h0 + "\n",
                                   "gust_nest_faces 1\nmd5 " + h0 + "\n", std::string("gust_nest_faces 1\n")}) {
        put(bad);
        const std::string unread = abort_message([&] { gust::NestFaces::check_restart(fixed.name, h0); });
        EXPECT_NE(unread.find("its NestFaces file cannot be read"), std::string::npos) << bad << unread;
    }
    // No record is written without nest faces.
    const ScopedDir none{"ut_nestfaces_rst_none"};
    amrex::UtilCreateCleanDirectory(none.name, true);
    gust::NestFaces::write_checkpoint(none.name, "");
    gust::NestFaces::check_restart(none.name, "");
}

} // namespace
