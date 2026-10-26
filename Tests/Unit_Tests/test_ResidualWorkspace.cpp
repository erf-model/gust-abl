// The residual's persistent workspace and the wall state handed from the
// ghost fill: over a 3-D hill that cuts the inflow face, with an outflow
// face, the y sides periodic or open (a slip wall and an outflow, so two
// non-periodic faces meet at corners), the wall model on a z0 map and the
// k–L closure, on uneven boxes,
// (a) the wall state fill_ghosts returns is the one compute_surface_state
// gives on the filled state, bit for bit (ghost columns included), so the
// residual with it is the residual without it; with z_ref in the top
// interval (the wall model reads the top w face, which the fill writes) the
// two differ, and the residual computes its own;
// (b) a workspace reused across states and layouts, poisoned with NaN before
// every call, gives the bits of fresh fields;
// (c) jvp with the refill and a workspace gives the bits of the difference
// quotient built from public calls: the refilled state's wall state reaches
// the perturbed residual.
#include "BoundaryConditions.H"
#include "Closure.H"
#include "HybridConvection.H"
#include "Reductions.H"
#include "Residual.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace {

using gust::FaceType;
using gust::Loc;
using gust::ptest::PCase;

// 16 x 12 x 8 cells on 1000 m x 800 m x 500 m over a 3-D cosine hill of
// 120 m, shifted so that it cuts the inflow face (x = 0); x open, y periodic
// unless yopen; geometric layers; boxes of mgx x mgy cells.
PCase hill_case (int mgx, int mgy, bool yopen = false)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 12, 8);
    c.in.mesh.prob_hi = {1000, 800, 500};
    c.in.mesh.is_periodic = {0, yopen ? 0 : 1, 0};
    c.in.mesh.max_grid_size = {mgx, mgy};
    c.in.mesh.dz0 = amrex::Real(10);
    c.in.mesh.stretch_ratio = amrex::Real(1.3);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([] (double x, double y) {
        return 120 * 0.25 * (1 - std::cos(2 * M_PI * (x + 250) / 1000)) * (1 - std::cos(2 * M_PI * (y + 100) / 800));
    });
    return c;
}

// A non-uniform state: the manufactured flow, theta and k varying; the ghost
// cells start at values the fill must replace. scale multiplies the velocity
// (another state on the same mesh).
void hill_state (gust::State& q, const PCase& c, double scale)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
    const auto fl = gust::trtest::make_flow(c, 120);
    const amrex::Array<amrex::MultiFab*, 3> vel{&q.u, &q.v, &q.w};
    for (int d = 0; d < 3; ++d) {
        gust::trtest::fill_loc(*vel[d], gust::velocity_loc(d), c,
                               [&] (double x, double y, double z) { return scale * fl.vel(d, x, y, z); });
    }
    // w off zero on the top face, which the slip-wall fill resets.
    q.w.plus(amrex::Real(0.3), 0, 1, 0);
    gust::trtest::fill_loc(q.theta, Loc::cell, c, [&] (double x, double y, double z) { return 300 + 0.003 * z + fl.q(x, y, z); });
    gust::trtest::fill_loc(q.tke, Loc::cell, c, [&] (double x, double y, double z) { return 0.3 + 0.1 * fl.q(x, y, z); });
    gust::trtest::fill_loc(q.p, Loc::cell, c, [&] (double x, double y, double z) { return 2 * fl.q(x, y, z); });
    for (auto* mf : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { mf->setBndry(amrex::Real(37)); }
}

// A perturbation of every field.
void perturbation (gust::State& v, const PCase& c, double phase)
{
    v.define(c.ba, c.dm);
    amrex::MultiFab::Copy(v.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
    const amrex::Array<amrex::MultiFab*, 6> f{&v.u, &v.v, &v.w, &v.p, &v.theta, &v.tke};
    const Loc loc[6] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell, Loc::cell};
    for (int n = 0; n < 6; ++n) {
        const double a = (n == 5) ? 0.02 : 1.0;
        gust::trtest::fill_loc(*f[n], loc[n], c, [&] (double x, double y, double z) {
            return a * std::sin(x / 150 + n + phase) * std::cos(y / 130 - phase) * std::cos(z / 90);
        });
    }
}

// z0 varying from column to column (every column the fill and the fluxes
// read, ghost columns included). A free function: nvcc refuses a device
// lambda in a constructor.
void fill_z0_map (const PCase& c, amrex::MultiFab& z0)
{
    z0.define(gust::ground_slab(c.ba), c.dm, 1, amrex::IntVect(1, 1, 0));
    for (amrex::MFIter mfi(z0, false); mfi.isValid(); ++mfi) {
        const auto a = z0.array(mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k) = amrex::Real(0.05) * (1 + amrex::Real(0.5) * std::sin(amrex::Real(0.7) * i + amrex::Real(0.3) * j));
        });
    }
    z0.FillBoundary(c.geom.periodicity());
}

// The level's set-up: boundary conditions, the wall model on a z0 map, the
// closure with the realizability cap and the stratification gate, the
// upwind blend.
struct Level
{
    gust::BoundaryInputs bin;
    std::unique_ptr<gust::BoundaryConditions> bc;
    gust::SurfaceLayerOptions so;
    amrex::MultiFab z0, dist, sigma;
    std::unique_ptr<gust::Closure> closure;
    gust::ResidualInputs in;
    explicit Level (const PCase& c, bool yopen = false)
    {
        const FaceType ylo = yopen ? FaceType::slip_wall : FaceType::periodic;
        const FaceType yhi = yopen ? FaceType::outflow : FaceType::periodic;
        bin.face = {FaceType::inflow, ylo, FaceType::surface_layer, FaceType::outflow, yhi, FaceType::slip_wall};
        bin.inflow.direction = 15;
        bin.theta_ref.z = {0, 500};
        bin.theta_ref.theta = {300, amrex::Real(301.5)};
        bc = std::make_unique<gust::BoundaryConditions>(bin, c.geom, c.mesh, in.coriolis);
        fill_z0_map(c, z0);
        so.stability = gust::SurfaceLayerOptions::Stability::heat_flux;
        so.heat_flux = amrex::Real(-0.004);
        so.z0_map = &z0;
        gust::height_above_ground(c.z, dist);
        gust::ClosureOptions copt;
        copt.realizability = amrex::Real(0.309);
        copt.ri_gate = amrex::Real(0.25);
        closure = std::make_unique<gust::Closure>(copt);
        sigma.define(c.ba, c.dm, 1, 1);
        sigma.setVal(amrex::Real(0.4));
        in.closure = closure.get();
        in.wall_distance = &dist;
        in.theta_ref.z = bin.theta_ref.z;
        in.theta_ref.theta = bin.theta_ref.theta;
        in.momentum.limiter = gust::Limiter::van_albada;
        in.theta.limiter = gust::Limiter::van_albada;
        in.hybrid_sigma = &sigma;
        in.surface = &so;
        in.bc = bc.get();
    }
};

const char* const kSurfaceNames[] = {"ustar", "speed", "zref", "d0", "inv_L", "tau", "drag", "kwall", "area"};

amrex::Vector<const amrex::MultiFab*> surface_fields (const gust::SurfaceState& s)
{
    return {&s.ustar, &s.speed, &s.zref, &s.d0, &s.inv_L, &s.tau, &s.drag, &s.kwall, &s.area};
}

amrex::Vector<amrex::MultiFab*> surface_fields (gust::SurfaceState& s)
{
    return {&s.ustar, &s.speed, &s.zref, &s.d0, &s.inv_L, &s.tau, &s.drag, &s.kwall, &s.area};
}

amrex::Long surface_differences (const gust::SurfaceState& a, const gust::SurfaceState& b)
{
    const auto fa = surface_fields(a);
    const auto fb = surface_fields(b);
    amrex::Long n = 0;
    for (int f = 0; f < 9; ++f) { n += gust::count_bit_differences(*fa[f], *fb[f], amrex::IntVect(1, 1, 0)); }
    return n;
}

// The points where two residuals differ, per field.
void expect_same (const gust::ResidualFields& a, const gust::ResidualFields& b, const std::string& what)
{
    const auto fa = a.fields();
    const auto fb = b.fields();
    const char* names[] = {"u", "v", "w", "p", "theta", "tke"};
    for (int n = 0; n < 6; ++n) {
        EXPECT_EQ(gust::count_bit_differences(*fa[n], *fb[n], amrex::IntVect(0)), 0) << what << ": R_" << names[n];
    }
}

amrex::Long total_differences (const gust::ResidualFields& a, const gust::ResidualFields& b)
{
    amrex::Long n = 0;
    const auto fa = a.fields();
    const auto fb = b.fields();
    for (int f = 0; f < 6; ++f) { n += gust::count_bit_differences(*fa[f], *fb[f], amrex::IntVect(0)); }
    return n;
}

// Every field the workspace holds set to NaN (ghosts included): a field read
// before it is written then shows in the result.
void poison (amrex::MultiFab& mf)
{
    if (!mf.empty()) { mf.setVal(std::numeric_limits<amrex::Real>::quiet_NaN(), 0, mf.nComp(), mf.nGrowVect()); }
}

void poison (gust::ResidualWorkspace& w)
{
    for (auto& f : w.F) { poison(f); }
    for (auto& f : w.tmp) { poison(f); }
    for (auto& g : w.Gm) { for (auto& f : g) { poison(f); } }
    for (auto& f : w.Gt) { poison(f); }
    for (auto& f : w.Gk) { poison(f); }
    for (auto& g : w.Gc) { for (auto& f : g) { poison(f); } }
    for (auto* f : {&w.theta_ref, &w.d, &w.omega, &w.jvp_theta_ref}) { poison(*f); }
    for (auto* f : surface_fields(w.ss)) { poison(*f); }
    for (auto* f : {&w.cf.nu, &w.cf.gamma_theta, &w.cf.gamma_tke, &w.cf.nut, &w.cf.length, &w.cf.Rt, &w.cf.Ps, &w.cf.Pb,
                    &w.cf.eps, &w.cf.dSdk, &w.cf.cap_ratio, &w.cf.dSdk_exact, &w.cf.dnutdk}) { poison(*f); }
    for (auto* f : {&w.closure.N2, &w.closure.S2, &w.closure.N2g, &w.closure.cap}) { poison(*f); }
    if (!w.closure.dec.empty()) { w.closure.dec.setVal(0x5a5a5a5a); }
    for (auto* f : {&w.qp.u, &w.qp.v, &w.qp.w, &w.qp.p, &w.qp.theta, &w.qp.tke, &w.qp.z_phys_nd}) { poison(*f); }
    for (auto& f : w.base.fields()) { poison(*f.mf); }
}

// The lowest last-centre height of the mesh, from the range check of z_ref.
double lowest_last_centre (const PCase& c, const Level& l, const gust::State& q)
{
    gust::SurfaceLayerOptions so = l.so;
    so.zref = amrex::Real(1e4);
    gust::SurfaceState s;
    std::string msg;
    try {
        gust::compute_surface_state(q, so, l.in.buoyancy, c.mesh, c.geom, s);
    } catch (const std::runtime_error& e) {
        msg = e.what();
    }
    std::smatch m;
    if (!std::regex_search(msg, m, std::regex("\\(([-0-9.eE+]+) to ([-0-9.eE+]+) m\\)"))) {
        ADD_FAILURE() << "no range in: " << msg;
        return 0;
    }
    return std::stod(m[2].str());
}

} // namespace

// (a) The wall state of the fill is the residual's: equal bit for bit to
// compute_surface_state on the filled state, ghost columns included, on six
// uneven boxes, with y periodic and with y open; the residual with it equals
// the residual that computes its own. A wall state of another state changes
// the residual (the given state is the one used). With z_ref at the lowest
// last-centre height, that column's reference interval is the top one: the
// fill's state is not the filled state's, and the residual still equals the
// one that computes its own.
TEST(ResidualWorkspace, FillGivesTheResidualsWallState)
{
    for (const bool yopen : {false, true}) {
        const std::string tag = yopen ? "y open" : "y periodic";
        const PCase c = hill_case(8, 4, yopen);
        Level su(c, yopen);
        gust::State q, q2;
        hill_state(q, c, 1.0);
        hill_state(q2, c, 1.1);
        gust::SurfaceState from_fill, fresh, other;
        su.bc->fill_ghosts(q, &su.so, su.in.buoyancy, false, &from_fill);
        gust::compute_surface_state(q, su.so, su.in.buoyancy, c.mesh, c.geom, fresh);
        EXPECT_FALSE(from_fill.reads_top_face) << tag;
        const auto a = surface_fields(from_fill);
        const auto b = surface_fields(fresh);
        for (int n = 0; n < 9; ++n) {
            EXPECT_EQ(gust::count_bit_differences(*a[n], *b[n], amrex::IntVect(1, 1, 0)), 0) << tag << ": " << kSurfaceNames[n];
        }
        gust::ResidualFields R0, R1, R2;
        gust::residual(q, su.in, c.mesh, c.geom, R0);
        gust::ResidualWorkspace ws;
        gust::residual(q, su.in, c.mesh, c.geom, R1, nullptr, &from_fill, &ws);
        expect_same(R0, R1, tag + ", with the fill's wall state");
#ifndef GUST_DEBUG_CHECKS
        // (Debug builds refuse a wall state that is not q's.)
        su.bc->fill_ghosts(q2, &su.so, su.in.buoyancy, false, &other);
        gust::residual(q, su.in, c.mesh, c.geom, R2, nullptr, &other, &ws);
        EXPECT_GT(total_differences(R0, R2), 0) << tag << ": another state's wall state must show";
#endif
        // z_ref in the top interval of the lowest column.
        gust::State qt;
        hill_state(qt, c, 1.0);
        Level st(c, yopen);
        st.so.zref = static_cast<amrex::Real>(lowest_last_centre(c, st, qt));
        gust::SurfaceState top_fill, top_fresh;
        st.bc->fill_ghosts(qt, &st.so, st.in.buoyancy, false, &top_fill);
        gust::compute_surface_state(qt, st.so, st.in.buoyancy, c.mesh, c.geom, top_fresh);
        EXPECT_TRUE(top_fill.reads_top_face) << tag;
        EXPECT_TRUE(top_fresh.reads_top_face) << tag;
        EXPECT_GT(surface_differences(top_fill, top_fresh), 0) << tag << ": the top w face the fill wrote must show";
        gust::ResidualFields Rt0, Rt1;
        gust::residual(qt, st.in, c.mesh, c.geom, Rt0);
        gust::residual(qt, st.in, c.mesh, c.geom, Rt1, nullptr, &top_fill, &ws);
        expect_same(Rt0, Rt1, tag + ", z_ref in the top interval");
    }
}

// (b) One workspace through two states, another layout and back, poisoned
// with NaN before every call: every residual equals the one with fresh
// fields, bit for bit, and so do the closure's fields. The closure's gate,
// the realizability cap and the upwind blend are on, so every scratch field
// is written.
TEST(ResidualWorkspace, ReuseGivesFreshBits)
{
    const PCase ca = hill_case(8, 4);
    const PCase cb = hill_case(16, 6);
    Level sa(ca);
    Level sb(cb);
    gust::ResidualWorkspace ws;
    struct Run { const PCase* c; Level* s; double scale; const char* what; };
    for (const Run r : {Run{&ca, &sa, 1.0, "first"}, Run{&ca, &sa, 1.1, "second state"},
                        Run{&cb, &sb, 0.9, "other layout"}, Run{&ca, &sa, 1.05, "back"}}) {
        gust::State q;
        hill_state(q, *r.c, r.scale);
        r.s->bc->fill_ghosts(q, &r.s->so, r.s->in.buoyancy);
        gust::ResidualFields fresh, reused;
        gust::ClosureFields cf_fresh, cf_reused;
        gust::residual(q, r.s->in, r.c->mesh, r.c->geom, fresh, &cf_fresh);
        poison(ws);
        gust::residual(q, r.s->in, r.c->mesh, r.c->geom, reused, &cf_reused, nullptr, &ws);
        expect_same(fresh, reused, r.what);
        EXPECT_EQ(gust::count_bit_differences(cf_fresh.nu, cf_reused.nu, amrex::IntVect(1)), 0) << r.what;
        EXPECT_EQ(gust::count_bit_differences(cf_fresh.dSdk_exact, cf_reused.dSdk_exact, amrex::IntVect(0)), 0) << r.what;
        // Without the caller's closure fields the workspace holds them.
        gust::ResidualFields own;
        poison(ws);
        gust::residual(q, r.s->in, r.c->mesh, r.c->geom, own, nullptr, nullptr, &ws);
        expect_same(fresh, own, r.what);
    }
}

// (c) jvp with the refill (the Newton core's form) and a poisoned workspace,
// for two perturbations in turn, with the base residual computed inside and
// given: the bits of (R(fill(q + eps v)) - R(q)) / eps built from public
// calls with fresh fields (the perturbed residual computing its own wall
// state), at jvp's eps. In single precision jvp divides by its eps in
// double, which the returned Real does not carry: there the quotients agree
// to a few float ulps of each field's largest value instead.
TEST(ResidualWorkspace, JvpWithTheRefillGivesFreshBits)
{
    const PCase c = hill_case(8, 4);
    Level su(c);
    gust::State q;
    hill_state(q, c, 1.0);
    su.bc->fill_ghosts(q, &su.so, su.in.buoyancy);
    gust::ResidualWorkspace ws;
    gust::ResidualFields R0;
    gust::residual(q, su.in, c.mesh, c.geom, R0);
    for (const double phase : {0.0, 0.7}) {
        gust::State v;
        perturbation(v, c, phase);
        su.bc->fill_ghosts(v, &su.so, su.in.buoyancy, true);
        gust::ResidualFields J_ws, J_given;
        poison(ws);
        const amrex::Real eps = gust::jvp(q, v, su.in, c.mesh, c.geom, J_ws, nullptr, true, &ws);
        poison(ws);
        EXPECT_EQ(gust::jvp(q, v, su.in, c.mesh, c.geom, J_given, &R0, true, &ws), eps);
        // The reference: q + eps v (ghosts included), its full fill, its residual, minus R(q), over eps.
        gust::State qp;
        qp.define(c.ba, c.dm);
        amrex::MultiFab::Copy(qp.z_phys_nd, q.z_phys_nd, 0, 0, 1, q.z_phys_nd.nGrow());
        const amrex::Array<const amrex::MultiFab*, 6> qf{&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
        const amrex::Array<const amrex::MultiFab*, 6> vf{&v.u, &v.v, &v.w, &v.p, &v.theta, &v.tke};
        const amrex::Array<amrex::MultiFab*, 6> pf{&qp.u, &qp.v, &qp.w, &qp.p, &qp.theta, &qp.tke};
        for (int n = 0; n < 6; ++n) {
            amrex::MultiFab::LinComb(*pf[n], 1, *qf[n], 0, eps, *vf[n], 0, 0, 1, qf[n]->nGrow());
        }
        su.bc->fill_ghosts(qp, &su.so, su.in.buoyancy);
        gust::ResidualFields ref;
        gust::residual(qp, su.in, c.mesh, c.geom, ref);
        auto rf = ref.fields();
        const auto r0 = R0.fields();
        for (int n = 0; n < 6; ++n) {
            amrex::MultiFab::Subtract(*rf[n].mf, *r0[n].mf, 0, 0, 1, 0);
            rf[n].mf->mult(static_cast<amrex::Real>(1 / static_cast<double>(eps)), 0, 1, 0);
        }
        if constexpr (std::is_same_v<amrex::Real, double>) {
            expect_same(ref, J_ws, "base inside");
            expect_same(ref, J_given, "base given");
        } else {
            const auto fr = ref.fields();
            const auto f1 = J_ws.fields();
            const auto f2 = J_given.fields();
            for (int n = 0; n < 6; ++n) {
                const double scale = static_cast<double>(fr[n].mf->norm0(0, 0));
                for (const auto* J : {f1[n].mf, f2[n].mf}) {
                    amrex::MultiFab d(J->boxArray(), J->DistributionMap(), 1, 0);
                    amrex::MultiFab::LinComb(d, 1, *J, 0, -1, *fr[n].mf, 0, 0, 1, 0);
                    EXPECT_LE(static_cast<double>(d.norm0(0, 0)), 8 * std::numeric_limits<float>::epsilon() * scale) << fr[n].name;
                }
            }
        }
        EXPECT_GT(J_ws.u.norm0(0, 0), 0.0);   // not vacuous
    }
}
