// The buildings' momentum terms and masked pressure correction (BuildingForcing.H, Poisson.H),
// beyond the plan, sub-phase 5.2 (UT-5.2.x):
//   an exact log law gives the exact u* at beta 0, 0.5 and 0.9 (UT-5.2.1); the blend, its
//   10 z0 floor and no reference cell (UT-5.2.2); the sink, its a_P and the closed wall face,
//   the fixed rows on blocked faces (UT-5.2.3, UT-5.2.5); the masked Poisson on a slab equals
//   the discrete cosine of the fluid channel, on both GMRES back ends, and an enclosed
//   block has phi = 0 and no flux through it (UT-5.2.6); GMRES iterations for 0, 1 and 9
//   buildings (UT-5.2.7); a pocket under a lifted building is counted as sealed.
#include "BuildingForcing.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"
#include "PoissonTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#ifdef _WIN32
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-12 : 2e-4;
constexpr double kKappa = 0.41;
const double kPi = std::acos(-1.0);

std::string write_text (const std::string& name, const std::function<void(std::ostream&)>& fill)
{
    static const long id = [] () {
        long pid = static_cast<long>(::getpid());
        amrex::ParallelDescriptor::Bcast(&pid, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        return pid;
    }();
    const fs::path dir = fs::temp_directory_path() / ("gust_building_forcing_" + std::to_string(id));
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(dir / name);
        os.precision(17);
        fill(os);
    }
    amrex::ParallelDescriptor::Barrier();
    return (dir / name).string();
}

// 8 x 8 x 24 cells of 10 m; a block on tiles 30..50 m in x and y, 25 m high.
PCase roof_case (int mgs)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(8, 8, 24);
    c.in.mesh.prob_hi = {80, 80, 240};
    c.in.mesh.is_periodic = {0, 0, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

gust::BuildingsOptions roof_options (const std::string& form = "fraction_stress")
{
    gust::BuildingsOptions o;
    o.type = gust::BuildingsType::immersed_forcing;
    o.type_name = "immersed_forcing";
    o.small_volfrac = 0;
    o.wall_form_name = form;
    o.wall_form = form == "drag" ? gust::BuildingWallForm::drag : gust::BuildingWallForm::fraction_stress;
    o.file = write_text("roof.txt", [] (std::ostream& os) {
        os << "5\n5\n";
        for (int r = 0; r < 2; ++r) { for (int t = 0; t < 5; ++t) { os << 20 * t << "\n"; } }
        for (int i = 0; i < 5; ++i) { for (int j = 0; j < 5; ++j) { os << ((i == 2 && j == 2) ? 25 : 0) << "\n"; } }
    });
    return o;
}

// The value at one point, on every rank, or -1e30 outside the field's points.
double value_at (const amrex::MultiFab& mf, int i, int j, int k, int n = 0)
{
    double v = -1e30;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        if (mfi.validbox().contains(amrex::IntVect(i, j, k))) { v = static_cast<double>(mf.const_array(mfi)(i, j, k, n)); }
    }
    amrex::ParallelDescriptor::ReduceRealMax(v);
    return v;
}

void set_uniform (amrex::MultiFab& mf, amrex::Real v) { mf.setVal(v); }

// A mesh's fractions from solid(i, j, k): faces touching a solid cell are blocked.
void close_cells (PCase& c, const std::function<bool(int, int, int)>& solid)
{
    amrex::MultiFab beta(c.mesh.detJ.boxArray(), c.dm, 1, 0);
    amrex::MultiFab fx(c.mesh.ax.boxArray(), c.dm, 1, 0), fy(c.mesh.ay.boxArray(), c.dm, 1, 0), fz(c.mesh.az.boxArray(), c.dm, 1, 0);
    amrex::MultiFab* f[3] = {&fx, &fy, &fz};
    for (amrex::MFIter mfi(beta, false); mfi.isValid(); ++mfi) {
        const auto b = beta.array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { b(i, j, k) = solid(i, j, k) ? 1 : 0; });
    }
    for (int d = 0; d < 3; ++d) {
        const int di = d == 0, dj = d == 1, dk = d == 2;
        for (amrex::MFIter mfi(*f[d], false); mfi.isValid(); ++mfi) {
            const auto a = f[d]->array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                a(i, j, k) = (solid(i, j, k) || solid(i - di, j - dj, k - dk)) ? 1 : 0;
            });
        }
    }
    c.mesh.set_open_faces(fx, fy, fz, beta, c.geom);
}

gust::PoissonBC periodic_xy_bc ()
{
    gust::PoissonBC bc;
    for (const int d : {0, 1}) {
        bc.face[amrex::Orientation(d, amrex::Orientation::low)] = gust::PoissonFace::periodic;
        bc.face[amrex::Orientation(d, amrex::Orientation::high)] = gust::PoissonFace::periodic;
    }
    return bc;
}

PCase slab_case (int nx, int nz, int mgs)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(nx, nx, nz);
    c.in.mesh.prob_hi = {amrex::Real(10 * nx), amrex::Real(10 * nx), amrex::Real(10 * nz)};
    c.in.mesh.is_periodic = {1, 1, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

// UT-5.2.1: the log law through d1 and d2 gives its own u*.
TEST (BuildingLaw, ExactLogLawGivesItsUstar)
{
    const double ustar = 0.37, z0 = 0.02;
    for (const double beta : {0.0, 0.5, 0.9}) {
        const double d1 = 0.5 * (1 - beta) * 10, d2 = (1 - beta) * 10 + 5;
        const double uo = ustar / kKappa * std::log(d1 / z0), ur = ustar / kKappa * std::log(d2 / z0);
        const double got = gust::form7_ustar2(amrex::Real(uo), amrex::Real(ur), amrex::Real(beta), amrex::Real(d1), amrex::Real(d2),
                                              amrex::Real(z0), amrex::Real(kKappa));
        EXPECT_NEAR(got, ustar * ustar, 10 * rtol) << "beta " << beta;
    }
}

// UT-5.2.2: the blend by hand; the floor; no reference cell.
TEST (BuildingLaw, BlendFloorAndNoReference)
{
    const double z0 = 0.1, uo = 3, ur = 5, d1 = 2.5, d2 = 10, beta = 0.5;
    const auto law = [&] (double b, double a, double c, double zz) {
        return static_cast<double>(gust::form7_ustar2(amrex::Real(uo), amrex::Real(ur), amrex::Real(b), amrex::Real(a), amrex::Real(c),
                                                      amrex::Real(zz), amrex::Real(kKappa)));
    };
    const double so = kKappa * uo / std::log(d1 / z0), sr = kKappa * ur / std::log(d2 / z0);
    EXPECT_NEAR(law(beta, d1, d2, z0), 0.5 * so * so + 0.5 * sr * sr, 10 * rtol);
    EXPECT_NEAR(law(0.9, d1, d2, z0), 0.1 * so * so + 0.9 * sr * sr, 10 * rtol);
    // Own height below 10 z0: read at 10 z0, finite as beta -> 1.
    const double sf = kKappa * uo / std::log(10.0);
    EXPECT_NEAR(law(0.999, 0.005, 0, z0), sf * sf, 10 * rtol);
    // d2 = 0: the own value alone, whatever beta.
    EXPECT_NEAR(law(beta, d1, 0, z0), so * so, 10 * rtol);
}

// UT-5.2.3, UT-5.2.5: on a partial roof cell (beta 0.5, k = 2).
TEST (BuildingForcing, SinkDiagonalAndClosedFace)
{
    for (const int mgs : {8, 4, 3}) {
        PCase c = roof_case(mgs);
        const gust::BuildingsOptions o = roof_options();
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1), w(c.mesh.az.boxArray(), c.dm, 1, 1);
        const double U = 3;
        set_uniform(u, amrex::Real(U));
        set_uniform(v, 0);
        set_uniform(w, 0);
        f.update(u, v);
        // d1 = 2.5, d2 = 10 (the next layer's centre), z0 = 0.1, a layer of 10 m, the whole z-face blocked.
        const double so = kKappa * U / std::log(25.0), sr = kKappa * U / std::log(100.0);
        const double r = (0.5 * so * so + 0.5 * sr * sr) / (10 * U);
        EXPECT_NEAR(value_at(f.rate(), 3, 3, 2), r, 100 * rtol * r);
        EXPECT_EQ(value_at(f.rate(), 3, 3, 3), 0.0) << "nothing above the partial cell";
        EXPECT_EQ(value_at(f.rate(), 3, 3, 1), 0.0) << "nothing in the solid";
        EXPECT_EQ(value_at(f.rate(), 2, 3, 2), 0.0) << "nothing beside the block";

        amrex::MultiFab Ru(c.mesh.ax.boxArray(), c.dm, 1, 0), Rv(c.mesh.ay.boxArray(), c.dm, 1, 0), Rw(c.mesh.az.boxArray(), c.dm, 1, 0);
        Ru.setVal(0);
        Rv.setVal(0);
        Rw.setVal(7);
        f.add_sink(u, v, Ru, Rv);
        EXPECT_NEAR(value_at(Ru, 4, 3, 2), r * U, 100 * rtol * r * U);
        EXPECT_NEAR(value_at(Ru, 3, 3, 2), 0.5 * r * U, 100 * rtol * r * U) << "half its control volume is over the roof";
        EXPECT_EQ(value_at(Rv, 3, 4, 2), 0.0) << "v is 0";

        // The fixed rows: R = q on a blocked face, others keep theirs.
        f.apply_rows(u, v, w, Ru, Rv, Rw);
        EXPECT_EQ(value_at(Ru, 4, 3, 1), U) << "u inside the block";
        EXPECT_EQ(value_at(Rw, 3, 3, 2), 0.0) << "the roof's z-face";
        EXPECT_EQ(value_at(Rw, 3, 3, 3), 7.0) << "the partial cell's top is open";
        EXPECT_EQ(value_at(Rw, 3, 3, 0), 7.0) << "w on the ground stays the boundary conditions' row";

        // a_P gains r; the closed face's coupling goes; the identity on a blocked face.
        amrex::MultiFab A(c.mesh.ax.boxArray(), c.dm, gust::n_implicit_coefs, 0);
        A.setVal(1);
        amrex::MultiFab open_z(c.mesh.az.boxArray(), c.dm, 1, 1);
        open_z.setVal(1);
        amrex::MultiFab::LinComb(open_z, 1, open_z, 0, -1, b.beta_face(2), 0, 0, 1, 0);
        open_z.FillBoundary(c.geom.periodicity());
        f.add_to_implicit(0, open_z, A);
        EXPECT_NEAR(value_at(A, 4, 3, 2, gust::coef_P), 1 + r, 100 * rtol);
        EXPECT_EQ(value_at(A, 4, 3, 2, gust::coef_ZM), 0.0);
        EXPECT_EQ(value_at(A, 4, 3, 2, gust::coef_ZP), 1.0);
        EXPECT_EQ(value_at(A, 4, 3, 1, gust::coef_P), 1.0);
        EXPECT_EQ(value_at(A, 4, 3, 1, gust::coef_ZP), 0.0);
        EXPECT_EQ(value_at(A, 3, 3, 2, gust::coef_ZM), 0.5) << "half the control volume's lower face is open";

        amrex::MultiFab mask(c.mesh.az.boxArray(), c.dm, 1, 0);
        mask.setVal(0);
        f.add_fixed_rows(2, mask);
        EXPECT_EQ(value_at(mask, 3, 3, 2), 1.0);
        EXPECT_EQ(value_at(mask, 3, 3, 3), 0.0);
        EXPECT_EQ(value_at(mask, 3, 3, 0), 0.0) << "nor is the ground's w a fixed row here";

        amrex::MultiFab Gu(amrex::convert(c.ba, amrex::IntVect(1, 0, 1)), c.dm, 1, 0), Gv(amrex::convert(c.ba, amrex::IntVect(0, 1, 1)), c.dm, 1, 0);
        Gu.setVal(2);
        Gv.setVal(2);
        f.close_viscous_fluxes(open_z, Gu, Gv);
        EXPECT_EQ(value_at(Gu, 4, 3, 2), 0.0);
        EXPECT_EQ(value_at(Gu, 3, 3, 2), 1.0);
        EXPECT_EQ(value_at(Gv, 3, 4, 3), 2.0);
    }
}

// UT-5.2.3: a roof or underside over half a cell's footprint, and partial faces in the mass fluxes.
TEST (BuildingForcing, PartCoveredWallsAndOpenFluxes)
{
    PCase c = roof_case(4);
    gust::BuildingsOptions o = roof_options();
    o.format = "xyz";
    const auto map = [] (const std::string& name, double h, double zs) {
        return write_text(name, [=] (std::ostream& os) {
            for (const double x : {0.0, 30.0, 40.0, 50.0, 80.0}) {
                for (const double y : {0.0, 20.0, 40.0, 60.0, 80.0}) {
                    const bool on = x == 40.0 && y == 40.0;   // the tile 35..45 m by 30..50 m
                    os << x << " " << y << " " << (on ? h : 0.0) << " " << (on ? zs : 0.0) << "\n";
                }
            }
        });
    };
    const double U = 3;
    const double su = kKappa * U / std::log(5 / 0.1);
    const double want = 0.5 * su * su / (10 * U);
    amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1);
    u.setVal(amrex::Real(U));
    v.setVal(0);
    {
        o.file = map("half_roof.txt", 20, 0);
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        f.update(u, v);
        EXPECT_NEAR(value_at(f.rate(), 3, 3, 2), want, 100 * rtol * want) << "half the footprint has a roof";
    }
    {
        o.file = map("half_under.txt", 40, 20);
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        f.update(u, v);
        EXPECT_NEAR(value_at(f.rate(), 3, 3, 1), want, 100 * rtol * want) << "half the footprint has an underside";
        amrex::MultiFab A(c.mesh.ax.boxArray(), c.dm, gust::n_implicit_coefs, 0);
        A.setVal(1);
        amrex::MultiFab open_z(c.mesh.az.boxArray(), c.dm, 1, 1);
        open_z.setVal(1);
        amrex::MultiFab::LinComb(open_z, 1, open_z, 0, -1, b.beta_face(2), 0, 0, 1, 0);
        open_z.FillBoundary(c.geom.periodicity());
        f.add_to_implicit(0, open_z, A);
        EXPECT_EQ(value_at(A, 4, 3, 1, gust::coef_ZP), 0.5) << "the row's upper face is half closed";
    }
    // Every face a quarter blocked: the transport fluxes carry three quarters.
    amrex::MultiFab bx(c.mesh.ax.boxArray(), c.dm, 1, 0), by(c.mesh.ay.boxArray(), c.dm, 1, 0), bz(c.mesh.az.boxArray(), c.dm, 1, 0);
    amrex::MultiFab beta(c.ba, c.dm, 1, 0);
    bx.setVal(0.25);
    by.setVal(0.25);
    bz.setVal(0.25);
    beta.setVal(0);
    c.mesh.set_open_faces(bx, by, bz, beta, c.geom);
    amrex::MultiFab uu(c.mesh.ax.boxArray(), c.dm, 1, 2), vv(c.mesh.ay.boxArray(), c.dm, 1, 2), ww(c.mesh.az.boxArray(), c.dm, 1, 1);
    uu.setVal(2);
    vv.setVal(3);
    ww.setVal(4);
    amrex::Array<amrex::MultiFab, 3> F;
    gust::transport_fluxes(uu, vv, ww, c.z, c.mesh, c.geom, F);
    EXPECT_NEAR(value_at(F[0], 4, 3, 2), 1.5, 10 * rtol);
    EXPECT_NEAR(value_at(F[1], 4, 3, 2), 2.25, 10 * rtol);
    EXPECT_NEAR(value_at(F[2], 4, 3, 2), 3.0, 10 * rtol);
}

TEST (BuildingForcing, DragFormHasFixedRowsOnly)
{
    PCase c = roof_case(8);
    const gust::BuildingsOptions o = roof_options("drag");
    gust::Buildings b;
    b.build(o, c.geom, c.z, c.ba, c.dm);
    gust::BuildingForcing f;
    f.build(b, o, c.geom, c.mesh.detJ);
    amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1);
    u.setVal(3);
    v.setVal(0);
    f.update(u, v);
    EXPECT_EQ(f.rate().norm0(), 0);
    amrex::MultiFab A(c.mesh.ax.boxArray(), c.dm, gust::n_implicit_coefs, 0);
    A.setVal(1);
    amrex::MultiFab open_z(c.mesh.az.boxArray(), c.dm, 1, 1);
    open_z.setVal(1);
    f.add_to_implicit(0, open_z, A);
    EXPECT_EQ(value_at(A, 4, 3, 2, gust::coef_P), 1.0);
    EXPECT_EQ(value_at(A, 4, 3, 1, gust::coef_ZP), 0.0) << "the fixed row stays";
}

// UT-5.2.6: the top half of the domain solid. In the fluid the cosine is an eigenvector of the
// discrete Laplacian with zero flux at the slab.
TEST (MaskedPoisson, SlabEqualsTheChannelsDiscreteCosine)
{
    for (const char* backend : {"gmres_fft", "gmres_mlmg"}) {
        for (const int mgs : {16, 7}) {
            const int nx = 16, nz = 32, nf = 16;
            const double dx = 10, dz = 10, L = nx * dx, h = nf * dz;
            PCase c = slab_case(nx, nz, mgs);
            close_cells(c, [=] (int, int, int k) { return k >= nf; });
            gust::PoissonOptions opt;
            opt.backend = backend;
            gust::PoissonSolver s(c.mesh, c.geom, c.z, opt);
            amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1);
            // In the solid rhs is left non-zero: its row must not take it.
            gust::ptest::fill_cells(rhs, c, [=] (double x, double, double z) {
                return z < h ? std::cos(2 * kPi * x / L) * std::cos(kPi * z / h) : 5.0;
            });
            phi.setVal(0);
            gust::FaceCoeffs fc;
            const gust::PoissonResult res = s.solve(rhs, phi, fc, periodic_xy_bc());
            const double lam = -(2 - 2 * std::cos(2 * kPi * dx / L)) / (dx * dx) - (2 - 2 * std::cos(kPi * dz / h)) / (dz * dz);
            EXPECT_TRUE(res.singular);
            amrex::MultiFab want(c.ba, c.dm, 1, 0);
            gust::ptest::fill_cells(want, c, [=] (double x, double, double z) {
                return z < h ? std::cos(2 * kPi * x / L) * std::cos(kPi * z / h) / lam : 0.0;
            });
            amrex::MultiFab::Subtract(want, phi, 0, 0, 1, 0);
            EXPECT_LT(want.norm0(0, 0), gust::ptest::tol(1e-6, 1e-1) * (1 / std::abs(lam))) << backend << " boxes " << mgs;
        }
    }
}

// UT-5.2.6: an enclosed block.
TEST (MaskedPoisson, EnclosedBlockHasNoFluxThroughIt)
{
    for (const char* backend : {"gmres_fft", "gmres_mlmg"}) {
        PCase c = slab_case(16, 16, 7);
        const auto solid = [] (int i, int j, int k) { return i >= 5 && i < 9 && j >= 5 && j < 9 && k < 6; };
        close_cells(c, solid);
        gust::PoissonOptions opt;
        opt.backend = backend;
        gust::PoissonSolver s(c.mesh, c.geom, c.z, opt);
        amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1);
        gust::ptest::fill_cells(rhs, c, [] (double x, double y, double z) { return std::cos(kPi * x / 80) * std::sin(kPi * y / 40) + 0.01 * z; });
        amrex::MultiFab keep(c.ba, c.dm, 1, 0);
        amrex::MultiFab::Copy(keep, rhs, 0, 0, 1, 0);
        phi.setVal(0);
        amrex::Array<amrex::MultiFab, 3> F;
        F[0].define(c.mesh.ax.boxArray(), c.dm, 1, 0);
        F[1].define(c.mesh.ay.boxArray(), c.dm, 1, 0);
        F[2].define(c.mesh.az.boxArray(), c.dm, 1, 0);
        gust::FaceCoeffs fc;
        const gust::PoissonResult res = s.solve(rhs, phi, fc, periodic_xy_bc(), &F);
        amrex::MultiFab div(c.ba, c.dm, 1, 0);
        gust::face_divergence(c.mesh, c.geom, {&F[0], &F[1], &F[2]}, div);
        // A phi = -div F in the fluid; rhs without its fluid mean; the block neither.
        amrex::MultiFab::Add(div, rhs, 0, 0, 1, 0);
        EXPECT_LT(div.norm0(0, 0), gust::ptest::tol(1e-7, 1e-1) * rhs.norm0(0, 0));
        for (amrex::MFIter mfi(phi, false); mfi.isValid(); ++mfi) {
            const auto p = phi.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (solid(i, j, k)) { EXPECT_EQ(p(i, j, k), 0); } });
        }
        (void)res;
        // The mean over the fluid is the one removed.
        double sum = 0, n = 0;
        for (amrex::MFIter mfi(keep, false); mfi.isValid(); ++mfi) {
            const auto a = keep.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { if (!solid(i, j, k)) { sum += a(i, j, k); n += 1; } });
        }
        amrex::ParallelDescriptor::ReduceRealSum(sum);
        amrex::ParallelDescriptor::ReduceRealSum(n);
        EXPECT_NEAR(res.removed_mean, sum / n, gust::ptest::tol(1e-11, 1e-3));
    }
}

// UT-5.2.7: the preconditioner without the mask; the iterations are reported and bounded.
TEST (MaskedPoisson, IterationsForZeroOneAndNineBuildings)
{
    for (const char* backend : {"gmres_fft", "gmres_mlmg"}) {
        for (const int count : {0, 1, 9}) {
            PCase c = slab_case(32, 16, 16);
            const auto solid = [=] (int i, int j, int k) {
                if (k >= 4 || count == 0) { return false; }
                const int n = count == 1 ? 1 : 3, s = 32 / n;
                const int ii = ((i % 32) + 32) % 32, jj = ((j % 32) + 32) % 32;
                return ii % s >= s / 2 - 2 && ii % s < s / 2 + 2 && jj % s >= s / 2 - 2 && jj % s < s / 2 + 2;
            };
            close_cells(c, solid);
            gust::PoissonOptions opt;
            opt.backend = backend;
            opt.reltol = gust::ptest::tol(1e-8, 1e-4);
            gust::PoissonSolver s(c.mesh, c.geom, c.z, opt);
            amrex::MultiFab rhs(c.ba, c.dm, 1, 0), phi(c.ba, c.dm, 1, 1);
            gust::ptest::fill_cells(rhs, c, [] (double x, double y, double z) { return std::cos(kPi * x / 40) * std::sin(kPi * y / 80) * std::cos(kPi * z / 160); });
            phi.setVal(0);
            gust::FaceCoeffs fc;
            const gust::PoissonResult res = s.solve(rhs, phi, fc, periodic_xy_bc());
            amrex::Print() << "UT-5.2.7: " << backend << ", " << count << " buildings: " << res.iterations << " GMRES iterations\n";
            EXPECT_LE(res.iterations, 60);
        }
    }
}

TEST (BuildingForcing, PocketUnderALiftedBuildingIsSealed)
{
    PCase c = roof_case(4);
    gust::BuildingsOptions o = roof_options();
    o.format = "xyz";
    const auto map = [&] (const std::string& name, bool closed) {
        return write_text(name, [=] (std::ostream& os) {
            for (int i = 0; i < 5; ++i) {
                for (int j = 0; j < 5; ++j) {
                    const bool inner = i >= 1 && i <= 3 && j >= 1 && j <= 3;
                    const bool centre = i == 2 && j == 2;
                    const bool gap = !closed && i == 2 && j == 1;
                    double h = inner ? 40 : 0, zs = 0;
                    if (centre || gap) { zs = 20; }
                    os << 20 * i << " " << 20 * j << " " << h << " " << zs << "\n";
                }
            }
        });
    };
    o.file = map("sealed.txt", true);
    {
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        // The pocket under the centre tile (20 x 20 m, 20 m): two cells across, two up.
        EXPECT_EQ(gust::count_sealed_fluid_cells(b, c.geom), 8);
    }
    o.file = map("open.txt", false);
    {
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        EXPECT_EQ(gust::count_sealed_fluid_cells(b, c.geom), 0) << "the gap opens it to the side";
    }
}

} // namespace
