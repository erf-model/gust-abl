// The buildings' momentum terms and masked pressure correction (BuildingForcing.H, Poisson.H),
// beyond the plan, sub-phase 5.2 (UT-5.2.x):
//   an exact log law gives the exact u* at beta 0, 0.5 and 0.9 (UT-5.2.1); the blend, its
//   10 z0 floor and no reference cell (UT-5.2.2); the sink, its a_P and the closed wall face,
//   the fixed rows on blocked faces (UT-5.2.3, UT-5.2.5); the masked Poisson on a slab equals
//   the discrete cosine of the fluid channel, on both GMRES back ends, and an enclosed
//   block has phi = 0 and no flux through it (UT-5.2.6); GMRES iterations for 0, 1 and 9
//   buildings (UT-5.2.7); a pocket under a lifted building is counted as sealed.
// Sub-phase 5.3 (UT-5.3.x), the sides of buildings:
//   an exact log law beside a side gives its stress on both components parallel to it, for
//   sides facing each way, and half of it where half the face is blocked (UT-5.3.1); an edge
//   cell (roof and side) and a corner cell (two sides) take the sum of their walls on each
//   component (UT-5.3.2); a partial side cell blends its reference cell; the viscous fluxes and
//   the implicit couplings through a side are closed, and roof_stress leaves the sides as they
//   were (UT-5.3.3); the same values on every box layout, with a side on a box edge and a
//   building across a periodic seam (UT-5.3.4).
// After #87: the speed of a cell with two walls from its faces that are not fixed rows
// (buildings.wall_speed), in the cell and in its reference cell, with the ground's w and partly blocked
// faces counted; the count of wall cells below the law's 10 z0 floor.
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
        // Binary: the same bytes, and so the same hash of the file, on every platform.
        std::ofstream os(dir / name, std::ios::binary);
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
    o.wall_form = form == "drag" ? gust::BuildingWallForm::drag
                : form == "roof_stress" ? gust::BuildingWallForm::roof_stress : gust::BuildingWallForm::fraction_stress;
    // As the reader's defaults: every face in the speed with the forms from before the open faces' mean.
    o.wall_speed = gust::default_wall_speed(o.wall_form);
    o.file = write_text("roof.txt", [] (std::ostream& os) {
        os << "5\n5\n";
        for (int r = 0; r < 2; ++r) { for (int t = 0; t < 5; ++t) { os << 20 * t << "\n"; } }
        for (int i = 0; i < 5; ++i) { for (int j = 0; j < 5; ++j) { os << ((i == 2 && j == 2) ? 25 : 0) << "\n"; } }
    });
    return o;
}

// A map in ERF's format on n x n points 20 m apart from 0: tile (ix, iy) covers 20 ix - 10 .. 20 ix + 10
// along x, so tile 2 is the cells 3 and 4 of roof_case. height(ix, iy) in m.
gust::BuildingsOptions tile_options (const std::string& name, const std::function<double(int, int)>& height,
                                     const std::string& form = "fraction_stress", int n = 5)
{
    gust::BuildingsOptions o = roof_options(form);
    o.file = write_text(name, [&] (std::ostream& os) {
        os << n << "\n" << n << "\n";
        for (int r = 0; r < 2; ++r) { for (int p = 0; p < n; ++p) { os << 20 * p << "\n"; } }
        for (int i = 0; i < n; ++i) { for (int j = 0; j < n; ++j) { os << height(i, j) << "\n"; } }
    });
    return o;
}

// f(i, j, k) at every point of mf, ghost points included.
void fill_by_index (amrex::MultiFab& mf, const std::function<double(int, int, int)>& f)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const auto a = mf.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { a(i, j, k) = static_cast<amrex::Real>(f(i, j, k)); });
    }
}

// u, v, w with one ghost layer and the residuals' fields, on a case's mesh.
struct Fields
{
    amrex::MultiFab u, v, w, Ru, Rv, Rw;
    explicit Fields (const PCase& c)
        : u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1), w(c.mesh.az.boxArray(), c.dm, 1, 1),
          Ru(c.mesh.ax.boxArray(), c.dm, 1, 0), Rv(c.mesh.ay.boxArray(), c.dm, 1, 0), Rw(c.mesh.az.boxArray(), c.dm, 1, 0)
    {
        for (amrex::MultiFab* mf : {&u, &v, &w, &Ru, &Rv, &Rw}) { mf->setVal(0); }
    }
};

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

using OpenArray = amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>;
using FluxArray = amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, AMREX_SPACEDIM>;

// The open fractions 1 - beta_face of a geometry, as Mesh::set_open_faces makes them: one ghost
// layer, periodic images, else 1.
OpenArray open_of (const PCase& c, const gust::Buildings& b, amrex::MultiFab& ox, amrex::MultiFab& oy, amrex::MultiFab& oz)
{
    amrex::MultiFab* o[3] = {&ox, &oy, &oz};
    for (int d = 0; d < 3; ++d) {
        o[d]->define(b.beta_face(d).boxArray(), c.dm, 1, 1);
        o[d]->setVal(1);
        amrex::MultiFab::LinComb(*o[d], 1, *o[d], 0, -1, b.beta_face(d), 0, 0, 1, 0);
        o[d]->FillBoundary(c.geom.periodicity());
    }
    return {&ox, &oy, &oz};
}

// G[c][d] on the d-faces of component c's control volumes, every value val.
FluxArray fluxes_of (const PCase& c, amrex::Real val)
{
    FluxArray G;
    for (int n = 0; n < 3; ++n) {
        gust::define_cv_fluxes(gust::velocity_loc(n), c.ba, c.dm, G[n]);
        for (auto& g : G[n]) { g.setVal(val); }
    }
    return G;
}

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
        f.update(u, v, w);
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
        f.add_sink(u, v, w, Ru, Rv, Rw);
        EXPECT_EQ(value_at(f.rate(), 3, 3, 2, 2), 0.0) << "a roof is no sink on w";
        EXPECT_EQ(value_at(f.rate(), 3, 3, 2, 1), value_at(f.rate(), 3, 3, 2, 0)) << "the same rate on u and v";
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
        amrex::MultiFab open_x, open_y, open_z;
        f.add_to_implicit(0, open_of(c, b, open_x, open_y, open_z), A);
        EXPECT_NEAR(value_at(A, 4, 3, 2, gust::coef_P) - 1, r, 1e3 * rtol * r);
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

        FluxArray G = fluxes_of(c, 2);
        f.close_viscous_fluxes({&open_x, &open_y, &open_z}, G);
        EXPECT_EQ(value_at(G[0][2], 4, 3, 2), 0.0);
        EXPECT_EQ(value_at(G[0][2], 3, 3, 2), 1.0);
        EXPECT_EQ(value_at(G[1][2], 3, 4, 3), 2.0);
        EXPECT_EQ(value_at(G[2][2], 3, 3, 2), 2.0) << "w's own z-flux is not a wall-parallel stress";
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
    amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1), w(c.mesh.az.boxArray(), c.dm, 1, 1);
    u.setVal(amrex::Real(U));
    v.setVal(0);
    w.setVal(0);
    {
        o.file = map("half_roof.txt", 20, 0);
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        f.update(u, v, w);
        EXPECT_NEAR(value_at(f.rate(), 3, 3, 2), want, 100 * rtol * want) << "half the footprint has a roof";
    }
    {
        o.file = map("half_under.txt", 40, 20);
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        f.update(u, v, w);
        EXPECT_NEAR(value_at(f.rate(), 3, 3, 1), want, 100 * rtol * want) << "half the footprint has an underside";
        amrex::MultiFab A(c.mesh.ax.boxArray(), c.dm, gust::n_implicit_coefs, 0);
        A.setVal(1);
        amrex::MultiFab open_x, open_y, open_z;
        f.add_to_implicit(0, open_of(c, b, open_x, open_y, open_z), A);
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
    amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1), w(c.mesh.az.boxArray(), c.dm, 1, 1);
    u.setVal(3);
    v.setVal(2);
    w.setVal(1);
    f.update(u, v, w);
    for (int n = 0; n < 3; ++n) { EXPECT_EQ(f.rate().norm0(n), 0); }
    amrex::MultiFab A(c.mesh.ax.boxArray(), c.dm, gust::n_implicit_coefs, 0);
    A.setVal(1);
    amrex::MultiFab open_x, open_y, open_z;
    f.add_to_implicit(0, open_of(c, b, open_x, open_y, open_z), A);
    EXPECT_EQ(value_at(A, 4, 2, 1, gust::coef_YP), 1.0) << "the side's coupling stays: viscous flux to the zero velocity";
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

// UT-5.3.1, UT-5.3.4: a block on the cells 3, 4 in x and y, 45 m high (four layers, and half of the
// fifth). The speed parallel to a side follows a log law in the distance from it, at 37 degrees to
// the horizontal: the sink on each parallel component is the stress's component over the cell's width.
TEST (BuildingSides, LogLawGivesItsStressOnBothParallelComponents)
{
    const double ustar = 0.37, z0 = 0.1, ch = 0.8, cv = 0.6, dn = 10;
    // The distance of cell n's centre from the block's nearer side along one direction.
    const auto dist = [] (int n) { return n <= 2 ? 30.0 - (10.0 * n + 5) : n >= 5 ? (10.0 * n + 5) - 50 : 1.0; };
    const auto U = [&] (int n) { return ustar / kKappa * std::log(dist(n) / z0); };
    const double stress = ustar * ustar / dn, tol = 100 * rtol * stress;
    double norm_first[3] = {0, 0, 0};
    for (const int mgs : {8, 5, 3}) {   // 5 and 3: the high and the low side lie on box edges
        PCase c = roof_case(mgs);
        const gust::BuildingsOptions o = tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; });
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        {
            // Sides normal to x: v and w.
            Fields q(c);
            fill_by_index(q.v, [&] (int i, int, int) { return ch * U(i); });
            fill_by_index(q.w, [&] (int i, int, int) { return cv * U(i); });
            f.update(q.u, q.v, q.w);
            for (const int i : {5, 2}) {
                const double r = ustar * ustar / (dn * U(i));
                EXPECT_NEAR(value_at(f.rate(), i, 3, 1, 1), r, 100 * rtol * r) << "the side at cell " << i << ", boxes " << mgs;
                EXPECT_NEAR(value_at(f.rate(), i, 3, 1, 2), r, 100 * rtol * r);
                EXPECT_EQ(value_at(f.rate(), i, 3, 1, 0), 0.0) << "u is normal to this side";
                EXPECT_NEAR(value_at(f.rate(), i, 3, 4, 1), 0.5 * r, 100 * rtol * r) << "half of the face is blocked in the fifth layer";
                EXPECT_EQ(value_at(f.rate(), i, 3, 5, 1), 0.0) << "above the block";
            }
            EXPECT_EQ(value_at(f.rate(), 6, 3, 1, 1), 0.0) << "the reference cell is no wall cell";
            f.add_sink(q.u, q.v, q.w, q.Ru, q.Rv, q.Rw);
            for (const int i : {5, 2}) {
                EXPECT_NEAR(value_at(q.Rv, i, 4, 1), ch * stress, tol) << "v's row between two wall cells, side at cell " << i;
                EXPECT_NEAR(value_at(q.Rw, i, 3, 1), cv * stress, tol) << "w's row between two wall cells";
                EXPECT_NEAR(value_at(q.Rv, i, 5, 1), 0.5 * ch * stress, tol) << "half the row's control volume is beside the side";
            }
            EXPECT_EQ(q.Ru.norm0(), 0.0) << "u is 0";
            for (int n = 0; n < 3; ++n) {
                const double norm = static_cast<double>(f.rate().norm1(n));
                if (mgs == 8) { norm_first[n] = norm; }
                EXPECT_NEAR(norm, norm_first[n], 100 * rtol * norm_first[n]) << "component " << n << ", boxes " << mgs;
            }
        }
        {
            // Sides normal to y: u and w.
            Fields q(c);
            fill_by_index(q.u, [&] (int, int j, int) { return ch * U(j); });
            fill_by_index(q.w, [&] (int, int j, int) { return cv * U(j); });
            f.update(q.u, q.v, q.w);
            f.add_sink(q.u, q.v, q.w, q.Ru, q.Rv, q.Rw);
            for (const int j : {5, 2}) {
                const double r = ustar * ustar / (dn * U(j));
                EXPECT_NEAR(value_at(f.rate(), 3, j, 1, 0), r, 100 * rtol * r) << "the side at cell " << j << ", boxes " << mgs;
                EXPECT_NEAR(value_at(f.rate(), 3, j, 1, 2), r, 100 * rtol * r);
                EXPECT_EQ(value_at(f.rate(), 3, j, 1, 1), 0.0) << "v is normal to this side";
                EXPECT_NEAR(value_at(q.Ru, 4, j, 1), ch * stress, tol);
                EXPECT_NEAR(value_at(q.Rw, 3, j, 1), cv * stress, tol);
            }
            EXPECT_EQ(q.Rv.norm0(), 0.0) << "v is 0";
        }
    }
}

// UT-5.3.2: an L of three 45 m tiles has an inner corner at the cell (2, 2); a 20 m tile beside the
// L makes a step, whose roof cell (5, 3, 2) lies against the taller tile's side.
TEST (BuildingSides, EdgeAndCornerCellsSumTheirWalls)
{
    const double a = 3, bb = 2, cc = 1, z0 = 0.1, dn = 10;
    const auto rate_of = [&] (double s) { const double us = kKappa * s / std::log(5 / z0); return us * us / (dn * s); };
    const double rx = rate_of(std::hypot(bb, cc)), ry = rate_of(std::hypot(a, cc)), rz = rate_of(std::hypot(a, bb));
    for (const int mgs : {8, 3}) {
        PCase c = roof_case(mgs);
        const gust::BuildingsOptions o = tile_options("ell.txt", [] (int i, int j) {
            if ((i == 2 && j == 2) || (i == 1 && j == 2) || (i == 2 && j == 1)) { return 45; }
            // Two 20 m steps: beside the side at x = 50 m, and beside the one at y = 50 m.
            return ((i == 3 && j == 2) || (i == 2 && j == 3)) ? 20 : 0;
        });
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        Fields q(c);
        q.u.setVal(amrex::Real(a));
        q.v.setVal(amrex::Real(bb));
        q.w.setVal(amrex::Real(cc));
        f.update(q.u, q.v, q.w);
        // The corner: sides at its high x and high y.
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 2, 2, 1, gust::kWallMask)), gust::kWallXhi | gust::kWallYhi);
        EXPECT_NEAR(value_at(f.rate(), 2, 2, 1, 0), ry, 100 * rtol * ry) << "u: the side normal to y";
        EXPECT_NEAR(value_at(f.rate(), 2, 2, 1, 1), rx, 100 * rtol * rx) << "v: the side normal to x";
        EXPECT_NEAR(value_at(f.rate(), 2, 2, 1, 2), rx + ry, 100 * rtol * (rx + ry)) << "w: both sides";
        // The edge: a roof below and a side at its low x.
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 5, 3, 2, gust::kWallMask)), gust::kWallRoof | gust::kWallXlo);
        EXPECT_NEAR(value_at(f.rate(), 5, 3, 2, 0), rz, 100 * rtol * rz) << "u: the roof";
        EXPECT_NEAR(value_at(f.rate(), 5, 3, 2, 1), rz + rx, 100 * rtol * (rz + rx)) << "v: the roof and the side";
        EXPECT_NEAR(value_at(f.rate(), 5, 3, 2, 2), rx, 100 * rtol * rx) << "w: the side";
        // Diagonally off the block's top corner: no wall.
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 5, 5, 4, gust::kWallMask)), 0);
        for (int n = 0; n < 3; ++n) { EXPECT_EQ(value_at(f.rate(), 5, 5, 4, n), 0.0); }
        f.add_sink(q.u, q.v, q.w, q.Ru, q.Rv, q.Rw);
        EXPECT_NEAR(value_at(q.Rw, 2, 2, 2), (rx + ry) * cc, 100 * rtol * (rx + ry) * cc) << "w's row between two corner cells";
        EXPECT_NEAR(value_at(q.Rv, 5, 4, 2), (rz + rx) * bb, 100 * rtol * (rz + rx) * bb) << "v's row between two edge cells";
        // With the blocked faces at rest, as in a run. The corner cell's u and v each have one face on the
        // other side, a fixed row: the speed takes the open face's value (buildings.wall_speed = open_faces),
        // so the rates are those of the uniform field above. all_faces takes the mean with the zero, half.
        const amrex::MultiFab* qs[3] = {&q.u, &q.v, &q.w};
        for (int n = 0; n < 3; ++n) {
            amrex::MultiFab mask(b.beta_face(n).boxArray(), c.dm, 1, 1);
            mask.setVal(0);
            amrex::MultiFab::Copy(mask, b.beta_face(n), 0, 0, 1, 0);
            mask.FillBoundary(c.geom.periodicity());
            amrex::MultiFab& x = const_cast<amrex::MultiFab&>(*qs[n]);
            for (amrex::MFIter mfi(x, false); mfi.isValid(); ++mfi) {
                const auto a4 = x.array(mfi);
                const auto m4 = mask.const_array(mfi);
                amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { if (m4(i, j, k) == amrex::Real(1)) { a4(i, j, k) = 0; } });
            }
        }
        f.update(q.u, q.v, q.w);
        EXPECT_NEAR(value_at(f.rate(), 2, 2, 1, 1), rx, 100 * rtol * rx) << "the corner's v is the open face's";
        EXPECT_NEAR(value_at(f.rate(), 2, 2, 1, 0), ry, 100 * rtol * ry) << "the corner's u is the open face's";
        EXPECT_NEAR(value_at(f.rate(), 2, 2, 1, 2), rx + ry, 100 * rtol * (rx + ry));
        // The edge cell: its low x face is on the side, so the roof's speed takes u from the high face alone;
        // its low z face is on the roof, so the side's speed takes w from the high face alone.
        const double rzo = rate_of(std::hypot(a, bb)), rxe = rate_of(std::hypot(bb, 0.0 + cc));
        EXPECT_NEAR(value_at(f.rate(), 5, 3, 2, 0), rzo, 100 * rtol * rzo) << "the edge's u is the open face's";
        EXPECT_NEAR(value_at(f.rate(), 5, 3, 2, 2), rxe, 100 * rtol * rxe) << "the edge's w is the open face's";
        gust::BuildingsOptions all = o;
        all.wall_speed = gust::BuildingWallSpeed::all_faces;
        gust::BuildingForcing g;
        g.build(b, all, c.geom, c.mesh.detJ);
        g.update(q.u, q.v, q.w);
        const double rxh = rate_of(std::hypot(0.5 * bb, cc)), ryh = rate_of(std::hypot(0.5 * a, cc));
        EXPECT_NEAR(value_at(g.rate(), 2, 2, 1, 1), rxh, 100 * rtol * rxh) << "all_faces: the corner's v is half the open face's";
        EXPECT_NEAR(value_at(g.rate(), 2, 2, 1, 0), ryh, 100 * rtol * ryh) << "all_faces: the corner's u is half the open face's";
        EXPECT_NEAR(value_at(g.rate(), 2, 2, 1, 2), rxh + ryh, 100 * rtol * (rxh + ryh));
        const double rzh = rate_of(std::hypot(0.5 * a, bb)), rxeh = rate_of(std::hypot(bb, 0.5 * cc));
        EXPECT_NEAR(value_at(g.rate(), 5, 3, 2, 0), rzh, 100 * rtol * rzh) << "all_faces: the edge's u is half";
        EXPECT_NEAR(value_at(g.rate(), 5, 3, 2, 2), rxeh, 100 * rtol * rxeh) << "all_faces: the edge's w is half";
        // The other step's edge cell: a roof below and a side at its low y. The roof's v and the side's w each
        // have one fixed face; u has none. v feels the roof alone, w the side alone, u both.
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 3, 5, 2, gust::kWallMask)), gust::kWallRoof | gust::kWallYlo);
        const double ryo = rate_of(std::hypot(a, cc));
        EXPECT_NEAR(value_at(f.rate(), 3, 5, 2, 1), rzo, 100 * rtol * rzo) << "the roof's v is the open face's";
        EXPECT_NEAR(value_at(f.rate(), 3, 5, 2, 2), ryo, 100 * rtol * ryo) << "the y side's w is the open face's";
        EXPECT_NEAR(value_at(f.rate(), 3, 5, 2, 0), rzo + ryo, 100 * rtol * (rzo + ryo));
        const double rzv = rate_of(std::hypot(a, 0.5 * bb)), ryh2 = rate_of(std::hypot(a, 0.5 * cc));
        EXPECT_NEAR(value_at(g.rate(), 3, 5, 2, 1), rzv, 100 * rtol * rzv) << "all_faces: the roof's v is half";
        EXPECT_NEAR(value_at(g.rate(), 3, 5, 2, 2), ryh2, 100 * rtol * ryh2) << "all_faces: the y side's w is half";
        // A plain side cell and a plain roof cell have no fixed face among theirs: the same in both forms.
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 3, 0, 1, gust::kWallMask)), gust::kWallYhi);
        EXPECT_EQ(value_at(g.rate(), 3, 0, 1, 0), value_at(f.rate(), 3, 0, 1, 0));
        EXPECT_GT(value_at(f.rate(), 3, 0, 1, 0), 0.0);
        // The first edge cell's v, which runs along the edge: the roof's rate and the side's.
        EXPECT_NEAR(value_at(f.rate(), 5, 3, 2, 1), rzo + rxe, 100 * rtol * (rzo + rxe));
        EXPECT_NEAR(value_at(g.rate(), 5, 3, 2, 1), rzh + rxeh, 100 * rtol * (rzh + rxeh));
    }
}

// Cells of 10 m along x and 5 m along y: each side's rate has its own normal's cell size and distance.
TEST (BuildingSides, RatesUseTheCellSizeAlongEachNormal)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(8, 16, 24);
    c.in.mesh.prob_hi = {80, 80, 240};
    c.in.mesh.is_periodic = {0, 0, 0};
    c.in.mesh.max_grid_size = {4, 4};
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    const gust::BuildingsOptions o = tile_options("block45.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; });
    gust::Buildings b;
    b.build(o, c.geom, c.z, c.ba, c.dm);
    gust::BuildingForcing f;
    f.build(b, o, c.geom, c.mesh.detJ);
    Fields q(c);
    const double a = 3, bb = 2, cc = 1, z0 = 0.1;
    q.u.setVal(amrex::Real(a));
    q.v.setVal(amrex::Real(bb));
    q.w.setVal(amrex::Real(cc));
    f.update(q.u, q.v, q.w);
    const auto rate_of = [&] (double s, double dn) { const double us = kKappa * s / std::log(0.5 * dn / z0); return us * us / (dn * s); };
    const double rx = rate_of(std::hypot(bb, cc), 10), ry = rate_of(std::hypot(a, cc), 5);
    // The block is the cells 3, 4 along x and 6..9 along y.
    ASSERT_EQ(static_cast<int>(value_at(b.wall(), 5, 7, 1, gust::kWallMask)), gust::kWallXlo);
    ASSERT_EQ(static_cast<int>(value_at(b.wall(), 3, 10, 1, gust::kWallMask)), gust::kWallYlo);
    ASSERT_EQ(static_cast<int>(value_at(b.wall(), 3, 5, 1, gust::kWallMask)), gust::kWallYhi);
    EXPECT_NEAR(value_at(f.rate(), 5, 7, 1, 1), rx, 100 * rtol * rx);
    EXPECT_NEAR(value_at(f.rate(), 3, 10, 1, 0), ry, 100 * rtol * ry);
    EXPECT_NEAR(value_at(f.rate(), 3, 5, 1, 2), ry, 100 * rtol * ry);
}

// Fixed faces beyond the cell's own: in its reference cell, on the ground, and partly blocked ones.
TEST (BuildingSides, OpenFacesInTheReferenceCellOnTheGroundAndPartlyBlocked)
{
    const auto zero_fixed = [] (const PCase& c, const gust::Buildings& b, Fields& q) {
        amrex::MultiFab* qs[3] = {&q.u, &q.v, &q.w};
        for (int n = 0; n < 3; ++n) {
            amrex::MultiFab mask(b.beta_face(n).boxArray(), c.dm, 1, 1);
            mask.setVal(0);
            amrex::MultiFab::Copy(mask, b.beta_face(n), 0, 0, 1, 0);
            mask.FillBoundary(c.geom.periodicity());
            for (amrex::MFIter mfi(*qs[n], false); mfi.isValid(); ++mfi) {
                const auto a4 = qs[n]->array(mfi);
                const auto m4 = mask.const_array(mfi);
                amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { if (m4(i, j, k) == amrex::Real(1)) { a4(i, j, k) = 0; } });
            }
        }
    };
    const auto build_two = [] (const PCase& c, const gust::BuildingsOptions& o, const gust::Buildings& b, gust::BuildingForcing& open,
                               gust::BuildingForcing& every) {
        open.build(b, o, c.geom, c.mesh.detJ);
        gust::BuildingsOptions all = o;
        all.wall_speed = gust::BuildingWallSpeed::all_faces;
        every.build(b, all, c.geom, c.mesh.detJ);
    };
    for (const int mgs : {4, 3}) {   // 3: a box edge between the wall cell and its reference cell
        // The reference cell: a wall 35..45 m along x leaves the cell 3 half solid (side at its high x, reference
        // the cell 2); a second building fills the cell 2 for y 10..30 m, so the reference cell (2, 3) has its low
        // y face on that building, a fixed row, and the wall cell (3, 3) has none.
        PCase c = roof_case(mgs);
        gust::BuildingsOptions o = roof_options();
        o.format = "xyz";
        o.file = write_text("reference_face.txt", [] (std::ostream& os) {
            for (const double x : {12.0, 28.0, 32.0, 38.0, 52.0, 80.0}) {   // tiles 20..30, 30..35, 35..45 m
                for (const double y : {0.0, 20.0, 40.0, 60.0, 80.0}) {      // tiles 10..30, 30..50 m
                    const bool wall = x == 38.0 && y == 40.0, second = x == 28.0 && y == 20.0;
                    os << x << " " << y << " " << ((wall || second) ? 45.0 : 0.0) << "\n";
                }
            }
        });
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 3, 3, 1, gust::kWallMask)), gust::kWallXhi);
        ASSERT_EQ(value_at(b.beta(), 3, 3, 1), 0.5);
        ASSERT_EQ(value_at(b.beta_face(1), 2, 3, 1), 1.0) << "the reference cell's low y face is on the second building";
        ASSERT_EQ(value_at(b.beta_face(1), 3, 3, 1), 0.5) << "the wall cell's own low y face is half blocked: no fixed row";
        gust::BuildingForcing f, g;
        build_two(c, o, b, f, g);
        Fields q(c);
        const double own = 3, ref = 5;
        fill_by_index(q.v, [&] (int i, int, int) { return i == 3 ? own : i == 2 ? ref : 7.0; });
        zero_fixed(c, b, q);
        f.update(q.u, q.v, q.w);
        g.update(q.u, q.v, q.w);
        const auto rate_of = [&] (double r) {
            const double so = kKappa * own / std::log(25.0), sr = kKappa * r / std::log(100.0);
            return (0.5 * so * so + 0.5 * sr * sr) / (10 * own);
        };
        EXPECT_NEAR(value_at(f.rate(), 3, 3, 1, 1), rate_of(ref), 100 * rtol * rate_of(ref)) << "the reference speed is its open face's, boxes " << mgs;
        EXPECT_NEAR(value_at(g.rate(), 3, 3, 1, 1), rate_of(0.5 * ref), 100 * rtol * rate_of(ref)) << "all_faces: half of it";
    }
    {
        // The ground: a 5 m annex (half the first layer) against a 45 m wall at x = 50 m. Its roof cell (5, 3, 0)
        // has the wall at its low x. The cell's ground face is covered, but w on the ground is the boundary
        // conditions' row, not a fixed one: it stays in the mean in both forms.
        PCase c = roof_case(4);
        const gust::BuildingsOptions o = tile_options("annex5.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 45 : (i == 3 && j == 2) ? 5 : 0; });
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 5, 3, 0, gust::kWallMask)), gust::kWallRoof | gust::kWallXlo);
        ASSERT_EQ(value_at(b.beta(), 5, 3, 0), 0.5);
        ASSERT_EQ(value_at(b.beta_face(2), 5, 3, 0), 1.0);
        gust::BuildingForcing f, g;
        build_two(c, o, b, f, g);
        Fields q(c);
        const double bb = 2, cc = 1;
        q.v.setVal(amrex::Real(bb));
        q.w.setVal(amrex::Real(cc));
        zero_fixed(c, b, q);   // the covered ground face too: w = 0 there, as on flat ground
        f.update(q.u, q.v, q.w);
        g.update(q.u, q.v, q.w);
        // The side at x = 50 m: (v, w) = (2, 1/2) in the cell and in its reference cell (6, 3, 0), both half solid.
        const double s = std::hypot(bb, 0.5 * cc), so = kKappa * s / std::log(25.0), sr = kKappa * s / std::log(100.0);
        // The wall faces the cell's upper half alone (the lower half is the annex): its area is the difference of
        // the blocked parts of the cell's two x faces, 1 - 0.5.
        ASSERT_EQ(value_at(b.beta_face(0), 5, 3, 0), 1.0);
        ASSERT_EQ(value_at(b.beta_face(0), 6, 3, 0), 0.5);
        const double r = 0.5 * (0.5 * so * so + 0.5 * sr * sr) / (10 * s);
        EXPECT_NEAR(value_at(f.rate(), 5, 3, 0, 2), r, 100 * rtol * r) << "the ground's w is counted";
        EXPECT_EQ(value_at(g.rate(), 5, 3, 0, 2), value_at(f.rate(), 5, 3, 0, 2)) << "no fixed face among these: the same in both forms";
    }
    {
        // A street one cell wide between two 45 m walls (the cells 2 and 4 along x): a fluid cell with a side on
        // each x face. Each takes its face's whole blocked part, so the rate is twice one side's.
        PCase c = roof_case(4);
        gust::BuildingsOptions o = roof_options();
        o.format = "xyz";
        o.file = write_text("street.txt", [] (std::ostream& os) {
            for (const double x : {5.0, 15.0, 25.0, 35.0, 45.0, 55.0}) {   // tiles of 10 m: the cells 1..4
                for (const double y : {0.0, 20.0, 40.0, 60.0, 80.0}) {
                    os << x << " " << y << " " << (((x == 25.0 || x == 45.0) && y == 40.0) ? 45.0 : 0.0) << "\n";
                }
            }
        });
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 3, 3, 1, gust::kWallMask)), gust::kWallXlo | gust::kWallXhi);
        ASSERT_EQ(value_at(b.beta(), 3, 3, 1), 0.0);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        Fields q(c);
        q.v.setVal(2);
        f.update(q.u, q.v, q.w);
        const double us = kKappa * 2 / std::log(50.0), r = 2 * us * us / (10 * 2.0);
        EXPECT_NEAR(value_at(f.rate(), 3, 3, 1, 1), r, 100 * rtol * r);
    }
    {
        // Partly blocked faces count as open: the roof cells of a 25 m block have their x and y faces half blocked.
        PCase c = roof_case(4);
        const gust::BuildingsOptions o = roof_options();
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        ASSERT_EQ(value_at(b.beta_face(0), 3, 3, 2), 0.5);
        ASSERT_EQ(value_at(b.beta_face(0), 4, 3, 2), 0.5);
        gust::BuildingForcing f, g;
        build_two(c, o, b, f, g);
        Fields q(c);
        q.u.setVal(3);
        f.update(q.u, q.v, q.w);
        g.update(q.u, q.v, q.w);
        EXPECT_GT(value_at(f.rate(), 3, 3, 2, 0), 0.0);
        EXPECT_EQ(value_at(f.rate(), 3, 3, 2, 0), value_at(g.rate(), 3, 3, 2, 0));
    }
}

// On a slope a roof's cells are partly solid by a fraction that changes from column to column, so each has a
// side bit although the roof is continuous. The side's rate takes the area facing that way, the difference of the
// two opposite faces' blocked parts, not the whole face's (which is the building's own cross-section).
TEST (BuildingSides, SlopedRoofsFalseSideTakesItsOwnSmallArea)
{
    for (const int mgs : {8, 3}) {
        PCase c;
        c.in.mesh.n_cell = amrex::IntVect(8, 8, 24);
        c.in.mesh.prob_hi = {80, 80, 240};
        c.in.mesh.is_periodic = {0, 0, 0};
        c.in.mesh.max_grid_size = {mgs, mgs};
        c.in.terrain.smoothing = 0;
        c.make([] (double x, double) { return 0.1 * x; });   // the ground rises 1 m per cell along x
        const gust::BuildingsOptions o = tile_options("block25.txt", [] (int i, int j) { return (i == 2 && j == 2) ? 25 : 0; });
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        // The roof's cells, the two columns of the block along x: partly solid, with a side bit.
        int k_roof = -1;
        for (int k = 0; k < 6; ++k) { if (value_at(b.beta(), 3, 3, k) > 0 && value_at(b.beta(), 3, 3, k) < 1) { k_roof = k; } }
        ASSERT_GE(k_roof, 1);
        const int mask = static_cast<int>(value_at(b.wall(), 3, 3, k_roof, gust::kWallMask));
        ASSERT_NE(mask & gust::kWallRoof, 0);
        ASSERT_NE(mask & (gust::kWallXlo | gust::kWallXhi), 0) << "the slope gives the roof cell a side bit";
        const bool lo = (mask & gust::kWallXlo) != 0;
        const double f_lo = value_at(b.beta_face(0), 3, 3, k_roof), f_hi = value_at(b.beta_face(0), 4, 3, k_roof);
        const double area = lo ? f_lo - f_hi : f_hi - f_lo, face = lo ? f_lo : f_hi;
        ASSERT_GT(area, 0.0);
        ASSERT_LT(area, 0.5 * face) << "the area facing that way is a small part of the face's blocked part";
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        Fields q(c);
        q.w.setVal(2);   // w alone: parallel to the sides, normal to the roof
        f.update(q.u, q.v, q.w);
        // Only the side acts on w. Its speed is 2 m/s in the cell and in its reference cell (no face is zeroed here).
        const double beta = value_at(b.beta(), 3, 3, k_roof);
        const double d1 = value_at(b.wall(), 3, 3, k_roof, gust::kWallD1x), d2 = value_at(b.wall(), 3, 3, k_roof, gust::kWallD2x);
        const double us2 = static_cast<double>(gust::form7_ustar2(amrex::Real(2), amrex::Real(2), amrex::Real(beta), amrex::Real(d1), amrex::Real(d2),
                                                                  amrex::Real(0.1), amrex::Real(kKappa)));
        const double want = area * us2 / (10 * 2.0);
        EXPECT_NEAR(value_at(f.rate(), 3, 3, k_roof, 2), want, 1e3 * rtol * want) << "boxes " << mgs;
        EXPECT_LT(value_at(f.rate(), 3, 3, k_roof, 2), 0.5 * face * us2 / (10 * 2.0)) << "not the whole face's";
    }
}

// A side in the middle of a cell: a tile 10 m wide leaves the cells 3 and 4 half solid along one
// direction, each with its centroid 2.5 m and its reference cell 10 m from its side. The four
// cells across have four speeds, so each side must read its own reference cell.
TEST (BuildingSides, PartialSideCellsBlendTheirOwnReference)
{
    const double val[8] = {9, 9, 5, 3, 4, 6, 9, 9}, wfac = 0.5;
    const auto rate_of = [&] (int own, int ref) {
        const double so = kKappa * std::hypot(1.0, wfac) * val[own] / std::log(25.0), sr = kKappa * std::hypot(1.0, wfac) * val[ref] / std::log(100.0);
        return (0.5 * so * so + 0.5 * sr * sr) / (10 * std::hypot(1.0, wfac) * val[own]);
    };
    for (const bool along_x : {true, false}) {
        PCase c = roof_case(4);
        gust::BuildingsOptions o = roof_options();
        o.format = "xyz";
        o.file = write_text(along_x ? "half_side_x.txt" : "half_side_y.txt", [=] (std::ostream& os) {
            for (const double p : {0.0, 30.0, 40.0, 50.0, 80.0}) {       // across the sides: the tile 35..45 m
                for (const double s : {0.0, 20.0, 40.0, 60.0, 80.0}) {   // along them: the tile 30..50 m
                    const double h = (p == 40.0 && s == 40.0) ? 45.0 : 0.0;
                    if (along_x) { os << p << " " << s << " " << h << "\n"; } else { os << s << " " << p << " " << h << "\n"; }
                }
            }
        });
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        Fields q(c);
        const auto clampi = [] (int n) { return n < 0 ? 0 : n > 7 ? 7 : n; };
        // The component along the side and w, both functions of the index across it.
        fill_by_index(along_x ? q.v : q.u, [&] (int i, int j, int) { return val[clampi(along_x ? i : j)]; });
        fill_by_index(q.w, [&] (int i, int j, int) { return wfac * val[clampi(along_x ? i : j)]; });
        f.update(q.u, q.v, q.w);
        const int lo_bit = along_x ? gust::kWallXlo : gust::kWallYlo, hi_bit = along_x ? gust::kWallXhi : gust::kWallYhi;
        const int tang = along_x ? 1 : 0, normal = along_x ? 0 : 1;
        const auto at = [&] (const amrex::MultiFab& mf, int n, int comp) { return along_x ? value_at(mf, n, 3, 1, comp) : value_at(mf, 3, n, 1, comp); };
        ASSERT_EQ(static_cast<int>(at(b.wall(), 3, gust::kWallMask)), hi_bit) << "the cell 3 has the side at its high end";
        ASSERT_EQ(static_cast<int>(at(b.wall(), 4, gust::kWallMask)), lo_bit);
        ASSERT_EQ(at(b.beta(), 3, 0), 0.5);
        const double r3 = rate_of(3, 2), r4 = rate_of(4, 5);
        EXPECT_NEAR(at(f.rate(), 3, tang), r3, 100 * rtol * r3) << (along_x ? "x" : "y") << ": the reference is the cell 2";
        EXPECT_NEAR(at(f.rate(), 3, 2), r3, 100 * rtol * r3);
        EXPECT_NEAR(at(f.rate(), 4, tang), r4, 100 * rtol * r4) << (along_x ? "x" : "y") << ": the reference is the cell 5";
        EXPECT_NEAR(at(f.rate(), 4, 2), r4, 100 * rtol * r4);
        EXPECT_EQ(at(f.rate(), 3, normal), 0.0);
    }
}

// UT-5.3.3: no viscous flux of a parallel component through a side, in the residual's fluxes and in the rows.
TEST (BuildingSides, ViscousFluxesAndCouplingsThroughASideAreClosed)
{
    for (const int mgs : {8, 5, 3}) {
        PCase c = roof_case(mgs);
        const auto height = [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; };
        const gust::BuildingsOptions o = tile_options("block45.txt", height);
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        amrex::MultiFab open_x, open_y, open_z;
        const OpenArray open = open_of(c, b, open_x, open_y, open_z);
        FluxArray G = fluxes_of(c, 2);
        f.close_viscous_fluxes(open, G);
        EXPECT_EQ(value_at(G[1][0], 5, 4, 1), 0.0) << "v through the side at x = 50 m, boxes " << mgs;
        EXPECT_EQ(value_at(G[1][0], 3, 4, 1), 0.0) << "v through the side at x = 30 m";
        EXPECT_EQ(value_at(G[1][0], 5, 5, 1), 1.0) << "half of this flux's face is on the side";
        EXPECT_EQ(value_at(G[1][0], 5, 6, 1), 2.0) << "past the block";
        EXPECT_EQ(value_at(G[2][0], 5, 3, 1), 0.0) << "w through the side at x = 50 m";
        EXPECT_EQ(value_at(G[2][0], 5, 3, 5), 1.5) << "the mean of a half blocked face and an open one: the row's convention where a side ends mid-cell";
        EXPECT_EQ(value_at(G[0][1], 4, 5, 1), 0.0) << "u through the side at y = 50 m";
        EXPECT_EQ(value_at(G[0][1], 4, 3, 1), 0.0) << "u through the side at y = 30 m";
        EXPECT_EQ(value_at(G[2][1], 3, 5, 1), 0.0) << "w through the side at y = 50 m";
        EXPECT_EQ(value_at(G[0][2], 4, 3, 4), 0.0) << "the roof's, as before";
        EXPECT_EQ(value_at(G[0][0], 5, 3, 1), 2.0) << "the normal stress of u stays";
        EXPECT_EQ(value_at(G[1][1], 3, 5, 1), 2.0);
        EXPECT_EQ(value_at(G[2][2], 3, 3, 5), 2.0);

        // The rows: a_P keeps its value without a rate, and the coupling through the side goes.
        amrex::MultiFab Av(c.mesh.ay.boxArray(), c.dm, gust::n_implicit_coefs, 0), Aw(c.mesh.az.boxArray(), c.dm, gust::n_implicit_coefs, 0);
        amrex::MultiFab Au(c.mesh.ax.boxArray(), c.dm, gust::n_implicit_coefs, 0);
        for (amrex::MultiFab* A : {&Au, &Av, &Aw}) { A->setVal(1); }
        f.add_to_implicit(0, open, Au);
        f.add_to_implicit(1, open, Av);
        f.add_to_implicit(2, open, Aw);
        EXPECT_EQ(value_at(Av, 5, 4, 1, gust::coef_XM), 0.0);
        EXPECT_EQ(value_at(Av, 5, 4, 1, gust::coef_XP), 1.0);
        EXPECT_EQ(value_at(Av, 5, 4, 1, gust::coef_P), 1.0);
        EXPECT_EQ(value_at(Av, 5, 4, 1, gust::coef_ZM), 1.0);
        EXPECT_EQ(value_at(Av, 2, 4, 1, gust::coef_XP), 0.0);
        EXPECT_EQ(value_at(Av, 2, 4, 1, gust::coef_XM), 1.0);
        EXPECT_EQ(value_at(Av, 5, 5, 1, gust::coef_XM), 0.5);
        EXPECT_EQ(value_at(Aw, 5, 3, 1, gust::coef_XM), 0.0);
        EXPECT_EQ(value_at(Aw, 2, 3, 1, gust::coef_XP), 0.0);
        EXPECT_EQ(value_at(Aw, 3, 5, 1, gust::coef_YM), 0.0);
        EXPECT_EQ(value_at(Aw, 3, 2, 1, gust::coef_YP), 0.0);
        EXPECT_EQ(value_at(Aw, 5, 3, 1, gust::coef_ZM), 1.0) << "w's own direction";
        EXPECT_EQ(value_at(Au, 4, 5, 1, gust::coef_YM), 0.0);
        EXPECT_EQ(value_at(Au, 4, 2, 1, gust::coef_YP), 0.0);
        EXPECT_EQ(value_at(Au, 4, 5, 1, gust::coef_XM), 1.0) << "u's own direction";
        // With a rate, a_P gains the mean of the row's two cells, on w as well.
        Fields q(c);
        q.v.setVal(2);
        q.w.setVal(1);
        f.update(q.u, q.v, q.w);
        Aw.setVal(1);
        f.add_to_implicit(2, open, Aw);
        const double s = std::hypot(2.0, 1.0), us = kKappa * s / std::log(50.0), r = us * us / (10 * s);
        EXPECT_NEAR(value_at(Aw, 5, 3, 1, gust::coef_P) - 1, r, 1e3 * rtol * r);
        EXPECT_NEAR(value_at(Aw, 5, 3, 5, gust::coef_P) - 1, 0.25 * r, 1e3 * rtol * r) << "half a face below, none above";
        EXPECT_GT(value_at(f.rate(), 5, 3, 0, 2), 0.0);
        EXPECT_EQ(value_at(Aw, 5, 3, 0, gust::coef_P), 1.0) << "w on the ground is the boundary conditions' row: no rate in it";

        // roof_stress: the sides as they were, with the fixed rows alone.
        const gust::BuildingsOptions old = tile_options("block45.txt", height, "roof_stress");
        gust::BuildingForcing g;
        g.build(b, old, c.geom, c.mesh.detJ);
        FluxArray H = fluxes_of(c, 2);
        g.close_viscous_fluxes(open, H);
        EXPECT_EQ(value_at(H[1][0], 5, 4, 1), 2.0);
        EXPECT_EQ(value_at(H[2][1], 3, 5, 1), 2.0);
        EXPECT_EQ(value_at(H[0][2], 4, 3, 4), 0.0) << "the roof's is closed in both forms";
        q.u.setVal(3);
        f.update(q.u, q.v, q.w);
        g.update(q.u, q.v, q.w);
        EXPECT_EQ(value_at(g.rate(), 3, 3, 4, 0), value_at(f.rate(), 3, 3, 4, 0)) << "the roof's rate is the same in both forms";
        EXPECT_EQ(g.rate().norm0(2), 0.0);
        EXPECT_EQ(value_at(g.rate(), 5, 3, 1, 1), 0.0);
        EXPECT_GT(value_at(g.rate(), 3, 3, 4, 0), 0.0) << "the roof's rate";
        EXPECT_EQ(value_at(g.rate(), 3, 3, 4, 0), value_at(g.rate(), 3, 3, 4, 1));
        q.Rw.setVal(7);
        g.add_sink(q.u, q.v, q.w, q.Ru, q.Rv, q.Rw);
        EXPECT_EQ(q.Rw.min(0), 7.0);
        EXPECT_EQ(q.Rw.max(0), 7.0);
        Av.setVal(1);
        g.add_to_implicit(1, open, Av);
        EXPECT_EQ(value_at(Av, 5, 4, 1, gust::coef_XM), 1.0);
        EXPECT_EQ(value_at(Av, 3, 4, 4, gust::coef_ZM), 0.0) << "the roof's coupling is closed in both forms";
        Aw.setVal(1);
        g.add_to_implicit(2, open, Aw);
        EXPECT_EQ(value_at(Aw, 5, 3, 1, gust::coef_P), 1.0);
        EXPECT_EQ(value_at(Aw, 5, 3, 1, gust::coef_XM), 1.0);
    }
}

// UT-5.3.4: a block across the periodic seam in x (the tile -10..10 m: cells 7 and 0).
TEST (BuildingSides, BlockAcrossAPeriodicSeam)
{
    const double a = 3, bb = 2, cc = 1, z0 = 0.1, dn = 10;
    const auto rate_of = [&] (double s) { const double us = kKappa * s / std::log(5 / z0); return us * us / (dn * s); };
    const double rx = rate_of(std::hypot(bb, cc)), ry = rate_of(std::hypot(a, cc));
    for (const int mgs : {8, 4, 3}) {
        PCase c;
        c.in.mesh.n_cell = amrex::IntVect(8, 8, 24);
        c.in.mesh.prob_hi = {80, 80, 240};
        c.in.mesh.is_periodic = {1, 1, 0};
        c.in.mesh.max_grid_size = {mgs, mgs};
        c.in.terrain.smoothing = 0;
        c.make([] (double, double) { return 0.0; });
        const gust::BuildingsOptions o = tile_options("seam.txt", [] (int i, int j) { return (i == 0 && j == 2) ? 45 : 0; }, "fraction_stress", 4);
        gust::Buildings b;
        b.build(o, c.geom, c.z, c.ba, c.dm);
        gust::BuildingForcing f;
        f.build(b, o, c.geom, c.mesh.detJ);
        Fields q(c);
        q.u.setVal(amrex::Real(a));
        q.v.setVal(amrex::Real(bb));
        q.w.setVal(amrex::Real(cc));
        f.update(q.u, q.v, q.w);
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 1, 3, 1, gust::kWallMask)), gust::kWallXlo) << "boxes " << mgs;
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 6, 3, 1, gust::kWallMask)), gust::kWallXhi);
        EXPECT_NEAR(value_at(f.rate(), 1, 3, 1, 1), rx, 100 * rtol * rx);
        EXPECT_NEAR(value_at(f.rate(), 6, 3, 1, 1), rx, 100 * rtol * rx);
        f.add_sink(q.u, q.v, q.w, q.Ru, q.Rv, q.Rw);
        // The x-face on the seam, between the cells 7 and 0 beside the block's side at y = 50 m.
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 0, 5, 1, gust::kWallMask)), gust::kWallYlo);
        ASSERT_EQ(static_cast<int>(value_at(b.wall(), 7, 5, 1, gust::kWallMask)), gust::kWallYlo);
        EXPECT_NEAR(value_at(q.Ru, 0, 5, 1), ry * a, 100 * rtol * ry * a) << "the rate's periodic image";
        EXPECT_NEAR(value_at(q.Rw, 0, 5, 1), ry * cc, 100 * rtol * ry * cc);
        amrex::MultiFab open_x, open_y, open_z;
        const OpenArray open = open_of(c, b, open_x, open_y, open_z);
        FluxArray G = fluxes_of(c, 2);
        f.close_viscous_fluxes(open, G);
        EXPECT_EQ(value_at(G[0][1], 0, 5, 1), 0.0) << "u through the side, on the seam";
        EXPECT_EQ(value_at(G[1][0], 1, 4, 1), 0.0);
        EXPECT_EQ(value_at(G[1][0], 7, 4, 1), 0.0);
    }
}

// The checkpoint's hash of the inputs. Before the sides had the law, fraction_stress was today's
// roof_stress: roof_stress keeps that hash (pinned from main at bba99be for this file), so such a
// checkpoint continues with it, and today's fraction_stress has another.
TEST (BuildingSides, RestartHashFollowsTheLawNotItsName)
{
    const auto height = [] (int i, int j) { return (i == 2 && j == 2) ? 45 : 0; };
    const std::string roof = gust::Buildings::hash_inputs(tile_options("block45.txt", height, "roof_stress"));
    const std::string all = gust::Buildings::hash_inputs(tile_options("block45.txt", height, "fraction_stress"));
    const std::string drag = gust::Buildings::hash_inputs(tile_options("block45.txt", height, "drag"));
    EXPECT_EQ(roof, "eec5cda82396300ca69b7beb0b92bed2523f63c361c54227734df2dc0288c15a");
    EXPECT_NE(all, roof);
    EXPECT_NE(drag, roof);
    EXPECT_NE(drag, all);
    // The speed from every face has the hash from before the open faces' mean; the open faces' has another.
    gust::BuildingsOptions every = tile_options("block45.txt", height, "fraction_stress");
    every.wall_tke = false;
    const std::string open = gust::Buildings::hash_inputs(every);
    every.wall_speed = gust::BuildingWallSpeed::all_faces;
    EXPECT_NE(gust::Buildings::hash_inputs(every), open);
    // roof_stress may take the open faces too: then its hash is not the earlier one.
    gust::BuildingsOptions roofs = tile_options("block45.txt", height, "roof_stress");
    roofs.wall_speed = gust::BuildingWallSpeed::open_faces;
    roofs.wall_tke = false;
    EXPECT_NE(gust::Buildings::hash_inputs(roofs), roof);
    EXPECT_TRUE(gust::check_buildings_options(roofs).empty());
    static_assert(gust::default_wall_speed(gust::BuildingWallForm::roof_stress) == gust::BuildingWallSpeed::all_faces);
    static_assert(gust::default_wall_speed(gust::BuildingWallForm::fraction_stress) == gust::BuildingWallSpeed::open_faces);
}

// The mean over a cell's two faces with their weights, and the count of wall cells below the 10 z0 floor.
TEST (BuildingLaw, OpenFaceMeanAndCellsBelowTheFloor)
{
    const auto mean = [] (double a, double b, double wa, double wb) {
        return static_cast<double>(gust::open_face_mean(amrex::Real(a), amrex::Real(b), amrex::Real(wa), amrex::Real(wb)));
    };
    EXPECT_EQ(mean(3, 5, 1, 1), 4.0);
    EXPECT_EQ(mean(3, 5, 1, 0), 3.0);
    EXPECT_EQ(mean(3, 5, 0, 1), 5.0);
    EXPECT_EQ(mean(3, 5, 0, 0), 0.0) << "both faces fixed: at rest";
    // Both weights 1 give the bits of the earlier half sum.
    const amrex::Real a = amrex::Real(0.1), b = amrex::Real(0.7);
    EXPECT_EQ(gust::open_face_mean(a, b, 1, 1), amrex::Real(0.5) * (a + b));

    // A block 25 m high on 10 m cells: its four roof cells are half solid, 2.5 m from the roof; the side cells
    // are 5 m from their sides. With z0 = 0.3 m the floor is 3 m: the roof cells are below it, the sides not.
    PCase c = roof_case(4);
    gust::BuildingsOptions o = roof_options();
    gust::Buildings bd;
    bd.build(o, c.geom, c.z, c.ba, c.dm);
    gust::BuildingForcing f;
    f.build(bd, o, c.geom, c.mesh.detJ);
    // Wall cells: 4 on the roof; beside each of the 4 sides 2 cells in each of 3 layers.
    EXPECT_EQ(f.wall_cells(), 4 + 4 * 2 * 3);
    EXPECT_EQ(f.cells_below_floor(), 0) << "z0 = 0.1 m: the floor is 1 m";
    o.z0 = 0.3;
    f.build(bd, o, c.geom, c.mesh.detJ);
    EXPECT_EQ(f.cells_below_floor(), 4);
    o.z0 = 0.6;
    f.build(bd, o, c.geom, c.mesh.detJ);
    EXPECT_EQ(f.cells_below_floor(), 4 + 4 * 2 * 3) << "z0 = 0.6 m: the floor is 6 m, above the sides' 5 m too";
    gust::BuildingsOptions roofs = roof_options("roof_stress");
    roofs.z0 = 0.6;
    f.build(bd, roofs, c.geom, c.mesh.detJ);
    EXPECT_EQ(f.wall_cells(), 4) << "roof_stress: the law is on the roof alone";
    EXPECT_EQ(f.cells_below_floor(), 4);
    f.build(bd, roof_options("drag"), c.geom, c.mesh.detJ);
    EXPECT_EQ(f.wall_cells(), 0) << "drag: no wall law";
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
