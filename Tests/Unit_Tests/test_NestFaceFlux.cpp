// The conservative variant of the defect source (NestFaceFlux.H;
// nest.feedback_flux = nest), on a flat parent of 20 x 16 x 8 cells with
// uneven boxes and its nest at ratio 2 (24 x 20 cells, boxes that do not line
// up with the parent's):
//   - the nest's fluxes summed onto the parent faces: given fluxes G = a x^2 +
//     b y^2 + linear terms on every face location (the same function on both
//     domains), the correction on each ring point is the closed form
//       C = -s dxi^-1 (e_c - e_d + w (e_u - e_c)) / J,
//     with e the error of the summation rule for x^2 and y^2 (0 on a nest
//     face that lies on the parent face, a dx_f^2 / 4 for the mean of the
//     two nest faces either side of it or of r cell-centred faces, a dx_f^2 / 2
//     for the trapezoid rule over r + 1 nodes), s = +1 below chi_s, -1 above,
//     w the blend weight; 0 off the rings. Every equation (u, v, w, theta,
//     k), both directions, both sides; the nest's fluxes are set on their
//     valid points only, so a sum that reads across a nest box needs the
//     ghost exchange;
//   - from residual(): theta = s x, constant u = U, upwind convection: the
//     parent's upwind face value of q* lies dx_c / 2 upwind of the face, the
//     nest's dx_f / 2, so C_theta = -U s / 4 on the ring below chi_s along x
//     and +U s / 4 above (dx_c = 2 dx_f); 0 on the y rings (v = 0, theta
//     uniform in y) and for u and k (uniform). With van Leer convection and
//     the blend at weight w the high-order fluxes agree (exact for a linear
//     field) and C_theta = -+ w U s / 4.
#include "CVFluxes.H"
#include "Inputs.H"
#include "Mesh.H"
#include "NestDefect.H"
#include "NestFaceFlux.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <cmath>
#include <functional>
#include <type_traits>
#include <vector>

namespace {

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double kTol = kDouble ? 1.0e-10 : 2.0e-3;

struct Domain
{
    gust::Inputs in;
    amrex::Geometry geom{};
    amrex::BoxArray ba{};
    amrex::DistributionMapping dm{};
    gust::Mesh mesh{};
    gust::State q{};

    void make ()
    {
        geom = gust::make_geometry(in.mesh);
        ba = gust::make_box_array(in.mesh);
        dm = amrex::DistributionMapping(ba);
        q.define(ba, dm);
        mesh.define(ba, dm);
        mesh.build(in, geom, q.z_phys_nd);
    }
};

gust::Inputs parent_inputs ()
{
    gust::Inputs in;
    in.mesh.n_cell = amrex::IntVect(20, 16, 8);
    in.mesh.prob_lo = {-1000, 500, 0};
    in.mesh.prob_hi = {1000, 2100, 800};
    in.mesh.is_periodic = {0, 0, 0};
    in.mesh.max_grid_size = {7, 6};
    in.mesh.dz0 = 40;
    return in;
}

// Parent cells 3 to 14 in x and 2 to 11 in y at ratio 2.
gust::Inputs nest_inputs (const gust::Inputs& parent)
{
    gust::Inputs in = parent;
    in.mesh.n_cell = amrex::IntVect(24, 20, 8);
    in.mesh.prob_lo = {-700, 700, 0};
    in.mesh.prob_hi = {500, 1700, 800};
    in.mesh.max_grid_size = {8, 6};
    return in;
}

const amrex::IntVect kOffset(3, 2, 0);
constexpr int kR = 2;
constexpr int kLoc[5] = {1, 2, 4, 0, 0};

double pos (const amrex::Geometry& g, const amrex::MultiFab& mf, int i, int d)
{
    const double off = mf.ixType().nodeCentered(d) ? 0.0 : 0.5;
    return static_cast<double>(g.ProbLo(d)) + (i + off) * static_cast<double>(g.CellSize(d));
}

// fn on the points of b of every fab of mf (through a host copy); the other points keep their values.
void fill_on (amrex::MultiFab& mf, bool ghosts, const std::function<double(int, int, int)>& fn)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box fb = mfi.fabbox();
        const amrex::Box b = ghosts ? fb : mfi.validbox();
        amrex::FArrayBox f(fb, 1, amrex::The_Pinned_Arena());
        f.copy<amrex::RunOn::Device>(mf[mfi], fb);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) { f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(fn(i, j, k)); });
        mf[mfi].copy<amrex::RunOn::Device>(f, fb);
        amrex::Gpu::streamSynchronize();
    }
}

// The whole of mf (valid points) on every rank, indexed by its domain box.
std::vector<double> gather (const amrex::MultiFab& mf, const amrex::Box& dom)
{
    const amrex::Box b = amrex::convert(dom, mf.ixType());
    std::vector<double> out(static_cast<std::size_t>(b.numPts()), 0.0);
    std::vector<double> cnt(out.size(), 0.0);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        amrex::FArrayBox h(mf[mfi].box(), 1, amrex::The_Pinned_Arena());
        h.copy<amrex::RunOn::Device>(mf[mfi], mf[mfi].box());
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const auto n = static_cast<std::size_t>(b.index(amrex::IntVect(i, j, k)));
            out[n] += static_cast<double>(h(amrex::IntVect(i, j, k)));
            cnt[n] += 1.0;
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(out.data(), static_cast<int>(out.size()));
    amrex::ParallelDescriptor::ReduceRealSum(cnt.data(), static_cast<int>(cnt.size()));
    for (std::size_t n = 0; n < out.size(); ++n) { if (cnt[n] > 0.0) { out[n] /= cnt[n]; } }
    return out;
}

// The control-volume Jacobian at a point of type t (the mean of the cells it touches), from the gathered detJ.
double jac (const std::vector<double>& J, const amrex::Box& dom, int t, int i, int j, int k)
{
    double s = 0;
    int n = 0;
    for (int c = 0; c <= ((t >> 2) & 1); ++c) {
        for (int b = 0; b <= ((t >> 1) & 1); ++b) {
            for (int a = 0; a <= (t & 1); ++a) {
                const amrex::IntVect p(i - a, j - b, k - c);
                if (dom.contains(p)) { s += J[static_cast<std::size_t>(dom.index(p))]; ++n; }
            }
        }
    }
    return s / n;
}

// Which ring a parent point of type t is on: 0 none, else (d, side) packed as 1 + 2 d + side.
int ring_of (const gust::DefectRegion& g, int t, const amrex::IntVect& p)
{
    const amrex::Box base = amrex::convert(g.inner, amrex::IntVect(t & 1, (t >> 1) & 1, (t >> 2) & 1));
    for (int d = 0; d < 2; ++d) {
        for (int s = 0; s < 2; ++s) {
            amrex::Box b = base;
            b.setRange(d, s == 0 ? base.smallEnd(d) - 1 : base.bigEnd(d) + 1, 1);
            if (b.contains(p)) { return 1 + 2 * d + s; }
        }
    }
    return 0;
}

// A flux part: G = c0 + c1 x + c2 y + a x^2 + b y^2 + 0.1 k.
struct Quad
{
    double c0, c1, c2, a, b;
    double at (double x, double y, int k) const { return c0 + c1 * x + c2 * y + a * x * x + b * y * y + 0.1 * k; }
};

// The summation rule's error for the parent face of G-type tg along normal d (nest spacing h): see the header.
double rule_error (const Quad& q, int tg, int d, double hx, double hy)
{
    const double h[2] = {hx, hy};
    const double cf[2] = {q.a, q.b};
    double e = 0;
    for (int m = 0; m < 2; ++m) {
        const bool nod = ((tg >> m) & 1) != 0;
        if (m == d) { e += nod ? 0.0 : cf[m] * h[m] * h[m] / 4; }
        else { e += nod ? cf[m] * h[m] * h[m] / 2 : cf[m] * h[m] * h[m] / 4; }
    }
    return e;
}

// G[d] of a field at loc t filled from q (parent: all points; nest: valid points only, the ghosts 1e30).
void make_fluxes (const Domain& D, int t, const Quad& q, bool all_points,
                  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& G)
{
    gust::define_cv_fluxes(static_cast<gust::Loc>(t), D.ba, D.dm, G);
    for (int d = 0; d < 3; ++d) {
        G[d].setVal(amrex::Real(1.0e30));
        const amrex::MultiFab* m = &G[d];
        fill_on(G[d], all_points, [&] (int i, int j, int k) {
            return q.at(pos(D.geom, *m, i, 0), pos(D.geom, *m, j, 1), k);
        });
    }
}

// Passes the given parts of every equation to the probe.
struct GivenFluxes
{
    const Domain* D;
    Quad conv, up, diff;
    const amrex::MultiFab* sigma{nullptr};
    bool all_points{true};
    void run (gust::FluxProbe& p) const
    {
        for (int eq = 0; eq < 5; ++eq) {
            amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> G;
            make_fluxes(*D, kLoc[eq], conv, all_points, G);
            p.take(eq, gust::FluxPart::convective, G, nullptr);
            if (sigma != nullptr && eq < 4) {
                make_fluxes(*D, kLoc[eq], up, all_points, G);
                p.take(eq, gust::FluxPart::convective_upwind, G, sigma);
            }
            make_fluxes(*D, kLoc[eq], diff, all_points, G);
            p.take(eq, gust::FluxPart::diffusive, G, nullptr);
        }
    }
};

void summed_fluxes_case (double w)
{
    Domain p{parent_inputs()};
    p.make();
    Domain n{nest_inputs(p.in)};
    n.make();
    const amrex::Box dom = p.geom.Domain();
    const gust::DefectRegion g = gust::defect_region("nest.a", n.geom.Domain().size(), kR, kOffset, dom);
    amrex::MultiFab sig_c(p.ba, p.dm, 1, 1), sig_f(n.ba, n.dm, 1, 1);
    sig_c.setVal(amrex::Real(w));
    sig_f.setVal(amrex::Real(0.9));   // the nest's weight is not used
    const Quad conv{1.0, 0.01, -0.02, 1.0e-4, 3.0e-4};
    const Quad up{2.0, -0.01, 0.03, -2.0e-4, 0.5e-4};
    const Quad diff{-1.0, 0.02, 0.01, 0.7e-4, -1.1e-4};
    GivenFluxes gc{&p, conv, up, diff, w > 0 ? &sig_c : nullptr, true};
    GivenFluxes gf{&n, conv, up, diff, w > 0 ? &sig_f : nullptr, false};
    gust::NestFaceFlux ff(p.ba, p.dm, p.mesh, p.geom);
    ff.compute([&] (gust::FluxProbe& pr) { gc.run(pr); },
               {gust::FaceFluxNest{&g, &n.geom, kR, kOffset, [&] (gust::FluxProbe& pr) { gf.run(pr); }}});
    const std::vector<double> J = gather(p.mesh.detJ, dom);
    const double hx = n.geom.CellSize(0), hy = n.geom.CellSize(1);
    const amrex::MultiFab* C[5] = {&ff.correction().u, &ff.correction().v, &ff.correction().w,
                                   &ff.correction().theta, &ff.correction().tke};
    int on_ring = 0;
    double worst = 0, largest = 0;
    for (int eq = 0; eq < 5; ++eq) {
        const int t = kLoc[eq];
        const std::vector<double> c = gather(*C[eq], dom);
        const amrex::Box all = amrex::convert(dom, C[eq]->ixType());
        for (amrex::BoxIterator it(all); it.ok(); ++it) {
            const amrex::IntVect pt = it();
            if (t == 4 && (pt[2] == all.smallEnd(2) || pt[2] == all.bigEnd(2))) { continue; }   // w's boundary rows
            const int r = ring_of(g, t, pt);
            double expect = 0;
            if (r > 0) {
                const int d = (r - 1) / 2, s = (r - 1) % 2;
                const int tg = t ^ (1 << d);
                const double sgn = s == 0 ? 1.0 : -1.0;
                const double ww = (eq < 4) ? w : 0.0;
                const double e = rule_error(conv, tg, d, hx, hy) - rule_error(diff, tg, d, hx, hy)
                               + ww * (rule_error(up, tg, d, hx, hy) - rule_error(conv, tg, d, hx, hy));
                expect = -sgn * static_cast<double>(p.geom.InvCellSize(d)) * e / jac(J, dom, t, pt[0], pt[1], pt[2]);
                ++on_ring;
                largest = std::max(largest, std::abs(expect));
            }
            const double got = c[static_cast<std::size_t>(all.index(pt))];
            worst = std::max(worst, std::abs(got - expect));
            ASSERT_NEAR(got, expect, kTol) << "eq " << eq << " at " << pt << " ring " << r;
        }
    }
    EXPECT_GT(on_ring, 100);
    EXPECT_GT(largest, 1.0e-4) << "the closed form is zero: the test is vacuous";
}

TEST(NestFaceFlux, SummedNestFluxesOnEveryFaceLocation) { summed_fluxes_case(0.0); }

TEST(NestFaceFlux, SummedNestFluxesWithTheBlend) { summed_fluxes_case(0.3); }

// theta = s x and u = U on both domains; the probed residual() of each.
void linear_theta_case (bool blend)
{
    constexpr double kU = 7.0, kS = 0.004, kW = 0.5;
    Domain p{parent_inputs()};
    p.make();
    Domain n{nest_inputs(p.in)};
    n.make();
    for (Domain* D : {&p, &n}) {
        D->q.u.setVal(amrex::Real(kU));
        D->q.v.setVal(0.0);
        D->q.w.setVal(0.0);
        D->q.p.setVal(0.0);
        D->q.tke.setVal(0.3);
        const amrex::MultiFab* m = &D->q.theta;
        fill_on(D->q.theta, true, [&] (int i, int, int) { return 300.0 + kS * pos(D->geom, *m, i, 0); });
    }
    const amrex::Box dom = p.geom.Domain();
    const gust::DefectRegion g = gust::defect_region("nest.a", n.geom.Domain().size(), kR, kOffset, dom);
    amrex::MultiFab nu_c(p.ba, p.dm, 1, 1), nu_f(n.ba, n.dm, 1, 1), sig_c(p.ba, p.dm, 1, 1), sig_f(n.ba, n.dm, 1, 1);
    nu_c.setVal(2.0);
    nu_f.setVal(2.0);
    sig_c.setVal(amrex::Real(kW));
    sig_f.setVal(amrex::Real(kW));
    auto inputs = [&] (const amrex::MultiFab& nu, const amrex::MultiFab& sig) {
        gust::ResidualInputs in;
        in.nu = &nu;
        in.gamma_theta = &nu;
        in.gamma_tke = &nu;
        in.theta_ref.z = {0};
        in.theta_ref.theta = {300};
        if (blend) { in.hybrid_sigma = &sig; }
        else { in.theta.limiter = gust::Limiter::upwind; in.momentum.limiter = gust::Limiter::upwind; }
        return in;
    };
    gust::ResidualInputs in_c = inputs(nu_c, sig_c), in_f = inputs(nu_f, sig_f);
    gust::ResidualFields Rc, Rf;
    gust::NestFaceFlux ff(p.ba, p.dm, p.mesh, p.geom);
    ff.compute([&] (gust::FluxProbe& pr) { in_c.flux_probe = &pr; gust::residual(p.q, in_c, p.mesh, p.geom, Rc); },
               {gust::FaceFluxNest{&g, &n.geom, kR, kOffset, [&] (gust::FluxProbe& pr) {
                   in_f.flux_probe = &pr;
                   gust::residual(n.q, in_f, n.mesh, n.geom, Rf);
               }}});
    const double scale = blend ? kW : 1.0;
    const amrex::MultiFab* C[5] = {&ff.correction().u, &ff.correction().v, &ff.correction().w,
                                   &ff.correction().theta, &ff.correction().tke};
    int n_ring = 0;
    for (int eq = 0; eq < 5; ++eq) {
        const int t = kLoc[eq];
        const std::vector<double> c = gather(*C[eq], dom);
        const amrex::Box all = amrex::convert(dom, C[eq]->ixType());
        for (amrex::BoxIterator it(all); it.ok(); ++it) {
            const amrex::IntVect pt = it();
            if (t == 4 && (pt[2] == all.smallEnd(2) || pt[2] == all.bigEnd(2))) { continue; }
            const int r = ring_of(g, t, pt);
            double expect = 0;
            if (eq == 3 && (r == 1 || r == 2)) {
                expect = (r == 1 ? -1.0 : 1.0) * scale * kU * kS / 4;
                ++n_ring;
            }
            ASSERT_NEAR(c[static_cast<std::size_t>(all.index(pt))], expect, kTol * 100)
                << "eq " << eq << " at " << pt << " ring " << r;
        }
    }
    EXPECT_EQ(n_ring, 2 * 6 * 8);   // two x rings of 6 x 8 points
}

TEST(NestFaceFlux, LinearThetaUpwind) { linear_theta_case(false); }

TEST(NestFaceFlux, LinearThetaVanLeerWithTheBlend) { linear_theta_case(true); }

// NestDefect::add adds the correction to the source after set; before set it aborts (not tested: abort).
TEST(NestFaceFlux, DefectAddsTheCorrection)
{
    Domain p{parent_inputs()};
    p.make();
    const gust::DefectRegion g = gust::defect_region("nest.a", amrex::IntVect(24, 20, 8), kR, kOffset, p.geom.Domain());
    gust::ResidualFields R, C;
    R.define(p.ba, p.dm);
    C.define(p.ba, p.dm);
    for (auto& f : R.fields()) { f.mf->setVal(1.0); }
    for (auto& f : C.fields()) { f.mf->setVal(0.25); }
    gust::NestDefect S(p.ba, p.dm);
    S.set(R, {g});
    S.add(C);
    const std::vector<double> s = gather(S.source().theta, p.geom.Domain());
    const amrex::Box all = p.geom.Domain();
    for (amrex::BoxIterator it(all); it.ok(); ++it) {
        const double expect = g.inner.contains(it()) ? 1.25 : 0.25;
        ASSERT_NEAR(s[static_cast<std::size_t>(all.index(it()))], expect, kTol) << it();
    }
}

} // namespace
