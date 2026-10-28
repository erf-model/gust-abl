// The canopy's drag on the momentum equations (CanopyDrag.H), plan rev. 6
// beyond the plan, sub-phase 8.2 (UT-8.2.x):
//   the sink is f_c |U| u_c with f_c the face mean of the drag field and |U|
//   from the four-point means, against its closed form on linear velocity
//   fields, on two box layouts, ground faces included (UT-8.2.1);
//   the derivative diagonal is the finite-difference derivative of the sink,
//   picard is f_c |U| and none adds nothing, canopy_rows the same (UT-8.2.2);
//   the residual and implicit_rows carry the drag, and implicit_rows refuses
//   a canopy without its rows;
//   the momentum predictor's a_P carries the rows of the linearisation state,
//   also when it steps an increment from zero (the Newton preconditioner).
#include "CanopyDrag.H"
#include "CVFluxes.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"
#include "PoissonTestHelpers.H"
#include "Predictor.H"
#include "PseudoTime.H"
#include "Residual.H"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <string>
#include <type_traits>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-13 : 1e-5;

// The patch: centre (200, 200), diameter 300, height 17.3 m, Cd 0.25, LAI 4,
// uniform profile: f = Cd LAI / h in the cells wholly below the top.
constexpr double kH = 17.3, kCd = 0.25, kLai = 4.0, kR = 150.0;

// 16 x 16 x 24 cells over 400 m x 400 m, 1 m first layer stretched by 1.1,
// boxes of mgs cells, periodic in x and y, flat.
PCase make_case (int mgs)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 16, 24);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.mesh.dz0 = amrex::Real(1);
    c.in.mesh.stretch_ratio = amrex::Real(1.1);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

gust::Canopy make_canopy (const PCase& c, gust::CanopyDiagonal d)
{
    const fs::path dir = fs::temp_directory_path() / ("gust_canopy_drag_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / "patch.txt").string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        os << "1 200 200 " << kH << " " << 2 * kR << " " << kCd << " " << kLai << " 0.8\n";
    }
    amrex::ParallelDescriptor::Barrier();
    gust::CanopyOptions o;
    o.type = gust::CanopyType::patches;
    o.type_name = "patches";
    o.file = file;
    o.given = {"type", "file"};
    o.drag_diagonal = d;
    gust::Canopy cn;
    cn.build(o, c.geom, c.z, c.ba, c.dm);
    return cn;
}

// Linear functions of the index-space position (X, Y, Z), so that every
// four-point mean is the value at the face.
double gu (double X, double Y, double Z) { return 3 + 0.2 * X - 0.1 * Y + 0.05 * Z; }
double gv (double X, double Y, double Z) { return -1 + 0.1 * X + 0.15 * Y - 0.03 * Z; }
double gw (double X, double Y, double Z) { return 0.2 + 0.01 * X - 0.02 * Y + 0.04 * Z; }

// Component c's points (ghosts included) from g at its index-space position.
void fill_component (amrex::MultiFab& mf, int c, const std::function<double(double, double, double)>& g)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.fabbox();
        amrex::FArrayBox f(b, 1, amrex::The_Pinned_Arena());
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const double X = i + (c == 0 ? 0.0 : 0.5), Y = j + (c == 1 ? 0.0 : 0.5), Z = k + (c == 2 ? 0.0 : 0.5);
            f(amrex::IntVect(i, j, k)) = static_cast<amrex::Real>(g(X, Y, Z));
        });
        mf[mfi].copy<amrex::RunOn::Device>(f, b);
        amrex::Gpu::streamSynchronize();
    }
}

void fill_linear (gust::State& q)
{
    fill_component(q.u, 0, gu);
    fill_component(q.v, 1, gv);
    fill_component(q.w, 2, gw);
}

// The cell's drag in closed form: Cd LAI / h inside the patch and wholly
// below the top, 0 inside the patch and wholly above it or outside the
// patch; NaN for a part-filled cell.
double cell_f (int i, int j, int k, const amrex::Vector<amrex::Real>& zl)
{
    const double x = (i + 0.5) * 25 - 200, y = (j + 0.5) * 25 - 200;
    if (x * x + y * y > kR * kR) { return 0; }
    const double zb = zl[static_cast<std::size_t>(k)], zt = zl[static_cast<std::size_t>(k + 1)];
    if (zt <= kH) { return kCd * kLai / kH; }
    if (zb >= kH) { return 0; }
    return std::nan("");
}

// The largest |R_c - closed form| over component c's valid faces (faces
// next to a part-filled cell skipped) relative to the largest closed form,
// and the number of faces checked.
std::pair<double, amrex::Long> sink_error (const amrex::MultiFab& R, int c, const amrex::Vector<amrex::Real>& zl, int nz)
{
    double err = 0, big = 0;
    amrex::Long n = 0;
    for (amrex::MFIter mfi(R, false); mfi.isValid(); ++mfi) {
        const amrex::Box b = mfi.validbox();
        amrex::FArrayBox h(b, 1, amrex::The_Pinned_Arena());
        h.copy<amrex::RunOn::Device>(R[mfi], b);
        amrex::Gpu::streamSynchronize();
        amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
            const int di = (c == 0), dj = (c == 1), dk = (c == 2);
            const int ka = (c == 2 && k == 0) ? 0 : k - dk, kb = (c == 2 && k == nz) ? nz - 1 : k;
            const double fa = cell_f(i - di, j - dj, ka, zl), fb = cell_f(i, j, kb, zl);
            if (std::isnan(fa) || std::isnan(fb)) { return; }
            const double fc = 0.5 * (fa + fb);
            const double X = i + (c == 0 ? 0.0 : 0.5), Y = j + (c == 1 ? 0.0 : 0.5), Z = k + (c == 2 ? 0.0 : 0.5);
            const double U[3] = {gu(X, Y, Z), gv(X, Y, Z), gw(X, Y, Z)};
            const double want = fc * std::sqrt(U[0] * U[0] + U[1] * U[1] + U[2] * U[2]) * U[c];
            err = std::max(err, std::abs(static_cast<double>(h(amrex::IntVect(i, j, k))) - want));
            big = std::max(big, std::abs(want));
            ++n;
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(err);
    amrex::ParallelDescriptor::ReduceRealMax(big);
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return {big > 0 ? err / big : err, n};
}

// A host copy of mf's fab at mfi (valid box and ghosts).
amrex::FArrayBox host (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    const amrex::Box b = mf[mfi].box();
    amrex::FArrayBox h(b, mf.nComp(), amrex::The_Pinned_Arena());
    h.copy<amrex::RunOn::Device>(mf[mfi], b, 0, b, 0, mf.nComp());
    amrex::Gpu::streamSynchronize();
    return h;
}

// Smooth velocities, periodic over the 16 x 16 columns, sign-changing.
constexpr double kw = 2 * M_PI / 16;
void fill_periodic (gust::State& q)
{
    fill_component(q.u, 0, [] (double X, double Y, double Z) { return 2 * std::sin(kw * X + 0.3 * Z) + 0.3 * std::cos(kw * Y); });
    fill_component(q.v, 1, [] (double X, double Y, double Z) { return std::cos(kw * Y - kw * X) - 0.1 * Z; });
    fill_component(q.w, 2, [] (double X, double Y, double Z) { return 0.4 * std::sin(kw * (X + 2 * Y) + Z); });
}

void define_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
    q.p.setVal(0.0);
    q.theta.setVal(300.0);
    q.tke.setVal(0.1);
}

std::string abort_message (const std::function<void()>& f)
{
    try {
        f();
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return "";
}

} // namespace

// UT-8.2.1: R_c = f_c |U| u_c on every face, against the closed form, with
// 1 box and with 16 boxes (the faces on box edges take the ghost column),
// the ground and top z-faces taking the cell inside.
TEST(CanopyDrag, SinkHasItsClosedForm)
{
    for (const int mgs : {16, 4}) {
        PCase c = make_case(mgs);
        const gust::Canopy cn = make_canopy(c, gust::CanopyDiagonal::derivative);
        gust::CanopyFaceDrag fd;
        fd.build(cn, c.geom);
        gust::State q;
        define_state(q, c);
        fill_linear(q);
        amrex::MultiFab R[3];
        for (int d = 0; d < 3; ++d) {
            R[d].define(fd.f(d).boxArray(), c.dm, 1, 0);
            R[d].setVal(0.0);
        }
        gust::add_canopy_drag(fd, q.u, q.v, q.w, R[0], R[1], R[2]);
        const auto& zl = c.mesh.zlevels();
        for (int d = 0; d < 3; ++d) {
            const auto [err, n] = sink_error(R[d], d, zl, 24);
            EXPECT_LE(err, rtol) << "component " << d << ", max_grid_size " << mgs;
            EXPECT_GT(n, 16L * 16 * 20) << "component " << d;
        }
        // Some faces have drag, and the ground faces of w take the first cell's.
        EXPECT_GT(static_cast<double>(R[0].max(0)), 0.1);
        EXPECT_GT(static_cast<double>(R[2].max(0)), 0.0);
    }
}

// UT-8.2.2: the derivative diagonal equals the central difference of the
// sink in the face's own component (random velocities), picard is f_c |U|
// (the sink over u_c) and none adds nothing; every diagonal is >= 0.
TEST(CanopyDrag, DiagonalIsTheDerivativeOfTheSink)
{
    PCase c = make_case(8);
    gust::State q;
    define_state(q, c);
    fill_periodic(q);
    amrex::MultiFab* comps[3] = {&q.u, &q.v, &q.w};
    for (const auto d : {gust::CanopyDiagonal::derivative, gust::CanopyDiagonal::picard, gust::CanopyDiagonal::none}) {
        const gust::Canopy cn = make_canopy(c, d);
        gust::CanopyFaceDrag fd;
        fd.build(cn, c.geom);
        for (int cc = 0; cc < 3; ++cc) {
            const amrex::BoxArray& fba = fd.f(cc).boxArray();
            amrex::MultiFab A(fba, c.dm, gust::n_implicit_coefs, 0), S0(fba, c.dm, 1, 0), Sp(fba, c.dm, 1, 0),
                Sm(fba, c.dm, 1, 0), other[3];
            for (int e = 0; e < 3; ++e) {
                other[e].define(fd.f(e).boxArray(), c.dm, 1, 0);
                other[e].setVal(0.0);
            }
            A.setVal(0.0);
            gust::add_canopy_drag_to_implicit(fd, cc, q.u, q.v, q.w, A);
            // canopy_rows holds the same diagonal, for every option.
            gust::CanopyRows rows;
            gust::canopy_rows(fd, q.u, q.v, q.w, rows);
            amrex::MultiFab diff(fba, c.dm, 1, 0);
            amrex::MultiFab::Copy(diff, A, gust::coef_P, 0, 1, 0);
            amrex::MultiFab::Subtract(diff, rows.ap[static_cast<std::size_t>(cc)], 0, 0, 1, 0);
            EXPECT_EQ(static_cast<double>(diff.norm0()), 0.0) << "canopy_rows, component " << cc;
            const double h = 1e-6;
            auto sink = [&] (double shift, amrex::MultiFab& S) {
                comps[cc]->plus(static_cast<amrex::Real>(shift), 0, 1, comps[cc]->nGrow());
                S.setVal(0.0);
                amrex::MultiFab* R[3] = {&other[0], &other[1], &other[2]};
                R[cc] = &S;
                gust::add_canopy_drag(fd, q.u, q.v, q.w, *R[0], *R[1], *R[2]);
                comps[cc]->plus(static_cast<amrex::Real>(-shift), 0, 1, comps[cc]->nGrow());
            };
            sink(0, S0);
            sink(h, Sp);
            sink(-h, Sm);
            double err = 0, big = 0, neg = 0;
            for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
                const amrex::Box b = mfi.validbox();
                const amrex::FArrayBox ha = host(A, mfi), h0 = host(S0, mfi), hp = host(Sp, mfi), hm = host(Sm, mfi),
                                       hu = host(*comps[cc], mfi);
                const auto a = ha.const_array();
                const auto s0 = h0.const_array(), sp = hp.const_array(), sm = hm.const_array();
                const auto uc = hu.const_array();
                amrex::LoopOnCpu(b, [&] (int i, int j, int k) {
                    const double diag = a(i, j, k, gust::coef_P);
                    neg = std::min(neg, diag);
                    double want = 0;
                    if (d == gust::CanopyDiagonal::derivative) {
                        want = (sp(i, j, k) - sm(i, j, k)) / (2 * h);
                    } else if (d == gust::CanopyDiagonal::picard) {
                        if (std::abs(uc(i, j, k)) < 0.05) { return; }
                        want = s0(i, j, k) / uc(i, j, k);
                    }
                    err = std::max(err, std::abs(diag - want));
                    big = std::max(big, std::abs(want));
                });
            }
            amrex::ParallelDescriptor::ReduceRealMax(err);
            amrex::ParallelDescriptor::ReduceRealMax(big);
            amrex::ParallelDescriptor::ReduceRealMin(neg);
            const std::string what = gust::canopy_diagonal_name(d) + ", component " + std::to_string(cc);
            if (d == gust::CanopyDiagonal::none) {
                EXPECT_EQ(err, 0.0) << what;
            } else {
                EXPECT_GT(big, 0.05) << what;
                EXPECT_LE(err, (d == gust::CanopyDiagonal::derivative ? 1e-8 : rtol) * big) << what;
            }
            EXPECT_GE(neg, 0.0) << what;
        }
    }
}

// The residual's momentum rows gain the sink and implicit_rows' a_P the
// diagonal (no boundary conditions: periodic sides, every row); without the
// state implicit_rows refuses the canopy, and theta's rows are unchanged.
TEST(CanopyDrag, ResidualAndImplicitRowsCarryTheDrag)
{
    PCase c = make_case(8);
    gust::State q;
    define_state(q, c);
    fill_periodic(q);
    amrex::MultiFab nu(c.ba, c.dm, 1, 1);
    nu.setVal(0.5);
    gust::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.theta_ref.z = {0};
    in.theta_ref.theta = {300};
    const gust::Canopy cn = make_canopy(c, gust::CanopyDiagonal::derivative);
    gust::CanopyFaceDrag fd;
    fd.build(cn, c.geom);
    gust::ResidualFields R0, R1;
    gust::residual(q, in, c.mesh, c.geom, R0);
    gust::ResidualInputs inc = in;
    inc.canopy = &fd;
    gust::residual(q, inc, c.mesh, c.geom, R1);
    amrex::Array<amrex::MultiFab, 3> F;
    gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
    const amrex::Array<const amrex::MultiFab*, 3> Fp{&F[0], &F[1], &F[2]};
    gust::CanopyRows rows;
    gust::canopy_rows(fd, q.u, q.v, q.w, rows);
    amrex::MultiFab* r0[3] = {&R0.u, &R0.v, &R0.w};
    amrex::MultiFab* r1[3] = {&R1.u, &R1.v, &R1.w};
    const gust::Equation eqs[3] = {gust::Equation::u, gust::Equation::v, gust::Equation::w};
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab S(fd.f(d).boxArray(), c.dm, 1, 0), z[3];
        S.setVal(0.0);
        for (int e = 0; e < 3; ++e) {
            z[e].define(fd.f(e).boxArray(), c.dm, 1, 0);
            z[e].setVal(0.0);
        }
        amrex::MultiFab* R[3] = {&z[0], &z[1], &z[2]};
        R[d] = &S;
        gust::add_canopy_drag(fd, q.u, q.v, q.w, *R[0], *R[1], *R[2]);
        amrex::MultiFab::Subtract(*r1[d], *r0[d], 0, 0, 1, 0);
        amrex::MultiFab::Subtract(*r1[d], S, 0, 0, 1, 0);
        EXPECT_LE(static_cast<double>(r1[d]->norm0()), 1e-12 * static_cast<double>(S.norm0())) << "component " << d;
        EXPECT_GT(static_cast<double>(S.norm0()), 1e-3);

        amrex::MultiFab A0, A1, Dg(fd.f(d).boxArray(), c.dm, gust::n_implicit_coefs, 0);
        gust::implicit_rows(Fp, q.z_phys_nd, in, eqs[d], c.mesh, c.geom, A0);
        gust::implicit_rows(Fp, q.z_phys_nd, inc, eqs[d], c.mesh, c.geom, A1, nullptr, nullptr, &rows);
        Dg.setVal(0.0);
        gust::add_canopy_drag_to_implicit(fd, d, q.u, q.v, q.w, Dg);
        amrex::MultiFab::Subtract(A1, A0, 0, 0, gust::n_implicit_coefs, 0);
        amrex::MultiFab::Subtract(A1, Dg, 0, 0, gust::n_implicit_coefs, 0);
        EXPECT_LE(static_cast<double>(A1.norm0(gust::coef_P)), 1e-12 * static_cast<double>(Dg.norm0(gust::coef_P)));
        for (int n = 1; n < gust::n_implicit_coefs; ++n) { EXPECT_EQ(static_cast<double>(A1.norm0(n)), 0.0); }
        EXPECT_NE(abort_message([&] { gust::implicit_rows(Fp, q.z_phys_nd, inc, eqs[d], c.mesh, c.geom, A1); })
                      .find("with a canopy the momentum rows need the drag's diagonal"),
                  std::string::npos);
    }
    amrex::MultiFab T0, T1;
    gust::implicit_rows(Fp, q.z_phys_nd, in, gust::Equation::theta, c.mesh, c.geom, T0);
    gust::implicit_rows(Fp, q.z_phys_nd, inc, gust::Equation::theta, c.mesh, c.geom, T1);
    amrex::MultiFab::Subtract(T1, T0, 0, 0, gust::n_implicit_coefs, 0);
    EXPECT_EQ(static_cast<double>(T1.norm0()), 0.0);
    amrex::MultiFab::Subtract(R1.theta, R0.theta, 0, 0, 1, 0);
    EXPECT_EQ(static_cast<double>(R1.theta.norm0()), 0.0);
}

// The momentum predictor's a_P is that of the run without a canopy plus the
// rows computed at q before the step (canopy_rows), on all three components:
// the v and w rows do not see the u the step has just updated, and a step
// on an increment from zero (the Newton preconditioner's) has the same a_P
// with the rows of q. Without the rows the predictor refuses the canopy.
TEST(CanopyDrag, PredictorTakesTheRowsOfTheLinearisationState)
{
    PCase c = make_case(8);
    gust::State q;
    define_state(q, c);
    fill_periodic(q);
    amrex::MultiFab nu(c.ba, c.dm, 1, 1);
    nu.setVal(0.5);
    gust::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.theta_ref.z = {0};
    in.theta_ref.theta = {300};
    in.buoyancy.gravity = 0;
    const gust::Canopy cn = make_canopy(c, gust::CanopyDiagonal::derivative);
    gust::CanopyFaceDrag fd;
    fd.build(cn, c.geom);
    gust::ResidualInputs inc = in;
    inc.canopy = &fd;
    gust::CanopyRows rows;
    gust::canopy_rows(fd, q.u, q.v, q.w, rows);

    amrex::MultiFab dtau;
    gust::local_dtau(q, &nu, c.mesh, c.geom, 4, dtau);
    amrex::Array<amrex::MultiFab, 3> df, F;
    gust::face_dtau(dtau, c.mesh, c.geom, df);
    gust::transport_fluxes(q.u, q.v, q.w, q.z_phys_nd, c.mesh, c.geom, F);
    const amrex::Array<const amrex::MultiFab*, 3> Fp{&F[0], &F[1], &F[2]};
    auto step = [&] (const gust::ResidualInputs& ri, bool increment, const gust::CanopyRows* r, gust::PredictorFields& pf) {
        gust::State q1;
        define_state(q1, c);
        for (const auto& p : {std::make_pair(&q1.u, &q.u), std::make_pair(&q1.v, &q.v), std::make_pair(&q1.w, &q.w)}) {
            amrex::MultiFab::Copy(*p.first, *p.second, 0, 0, 1, p.second->nGrow());
        }
        gust::ResidualFields R;
        gust::residual(q1, ri, c.mesh, c.geom, R);
        if (increment) {
            for (auto* m : {&q1.u, &q1.v, &q1.w}) { m->setVal(0.0); }
        }
        gust::momentum_predictor(q1, R, Fp, ri, df, c.mesh, c.geom, 2, pf, nullptr, nullptr, 1,
                                 gust::LineSweepOptions(), 0, 0, r);
    };
    gust::PredictorFields p0, p1, p2;
    step(in, false, nullptr, p0);
    step(inc, false, &rows, p1);
    step(inc, true, &rows, p2);
    const amrex::MultiFab* a0[3] = {&p0.ap_u, &p0.ap_v, &p0.ap_w};
    const amrex::MultiFab* a1[3] = {&p1.ap_u, &p1.ap_v, &p1.ap_w};
    const amrex::MultiFab* a2[3] = {&p2.ap_u, &p2.ap_v, &p2.ap_w};
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab e(a1[d]->boxArray(), c.dm, 1, 0);
        amrex::MultiFab::Copy(e, *a1[d], 0, 0, 1, 0);
        amrex::MultiFab::Subtract(e, *a0[d], 0, 0, 1, 0);
        amrex::MultiFab::Subtract(e, rows.ap[static_cast<std::size_t>(d)], 0, 0, 1, 0);
        const double big = static_cast<double>(rows.ap[static_cast<std::size_t>(d)].norm0());
        EXPECT_GT(big, 0.0) << "component " << d;
        EXPECT_LE(static_cast<double>(e.norm0()), 1e-12 * static_cast<double>(a1[d]->norm0())) << "component " << d;
        amrex::MultiFab::Copy(e, *a2[d], 0, 0, 1, 0);
        amrex::MultiFab::Subtract(e, *a1[d], 0, 0, 1, 0);
        EXPECT_EQ(static_cast<double>(e.norm0()), 0.0) << "component " << d << " (increment)";
    }
    gust::PredictorFields p3;
    EXPECT_NE(abort_message([&] { step(inc, false, nullptr, p3); }).find("need the drag's diagonal"), std::string::npos);
}
