// The buildings' terms of theta (BuildingTheta.H), beyond the plan, sub-phase 5.5 (UT-5.5.x):
//   the walls' heat exchange equals the closed form of the form 7 law with a given surface temperature, on a
//   non-uniform theta: a roof's partly solid cell with its reference cell, a side with the flow along it, a side
//   half blocked, a side the flow is normal to (UT-5.5.1); with the reference profile as the walls' temperature a
//   state at that profile, stratified, exchanges nothing and its solid rows are zero (UT-5.5.2); the implicit
//   diagonal is the derivative of the exchange to the cell's own theta, by a finite difference (UT-5.5.3); the
//   diffusive flux of theta and its couplings are closed on blocked faces, so no heat is conducted out of a
//   building (UT-5.5.4). Also: solid cells hold theta at the walls' temperature; the heat the walls take, summed;
//   the solid cells out of the norm of theta; the options, their defaults and the restart hash.
#include "BoundaryConditions.H"
#include "BuildingTheta.H"
#include "PoissonTestHelpers.H"
#include "PseudoTime.H"
#include "Residual.H"
#include "Simplec.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

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
#include <string>
#include <type_traits>
#include <vector>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-12 : 2e-4;
// A temperature near 300 K, or a difference of two: round-off of 300 K (3e-5 K in single precision).
constexpr double ttol = std::is_same_v<amrex::Real, double> ? 1e-10 : 1e-3;
// A row's change on a base of 9: relative to the change in double; in single the round-off of 9 and of the 300 K
// in theta - theta_b (a few 1e-5 of the change's coefficient times 5 K), so 2 % of the change.
double rowtol (double change) { return (std::is_same_v<amrex::Real, double> ? 1e-9 : 0.02) * std::abs(change); }
constexpr double kKappa = 0.41;

std::string write_text (const std::string& name, const std::function<void(std::ostream&)>& fill)
{
    static const long id = [] () {
        long pid = static_cast<long>(::getpid());
        amrex::ParallelDescriptor::Bcast(&pid, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        return pid;
    }();
    const fs::path dir = fs::temp_directory_path() / ("gust_building_theta_" + std::to_string(id));
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(dir / name, std::ios::binary);
        os.precision(17);
        fill(os);
    }
    amrex::ParallelDescriptor::Barrier();
    return (dir / name).string();
}

// 8 x 8 x 24 cells of 10 m.
PCase mesh_case (int mgs, bool periodic = false)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(8, 8, 24);
    c.in.mesh.prob_hi = {80, 80, 240};
    c.in.mesh.is_periodic = {periodic ? 1 : 0, periodic ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.terrain.smoothing = 0;
    c.make([] (double, double) { return 0.0; });
    return c;
}

// One block on the tile 30 .. 50 m in x and y (the cells 3 and 4), `top` m high, in ERF's format, on n x n points
// 20 m apart (4: the tiles span the 80 m of a periodic domain).
gust::BuildingsOptions block_options (double top, gust::BuildingThetaMode mode, double surface = 0, int n = 5)
{
    gust::BuildingsOptions o;
    o.type = gust::BuildingsType::immersed_forcing;
    o.type_name = "immersed_forcing";
    o.small_volfrac = 0;
    o.z0 = 0.1;
    o.theta_mode = mode;
    o.theta_surface = surface;
    o.file = write_text("block" + std::to_string(static_cast<int>(top)) + "_" + std::to_string(n) + ".txt", [&] (std::ostream& os) {
        os << n << "\n" << n << "\n";
        for (int r = 0; r < 2; ++r) { for (int p = 0; p < n; ++p) { os << 20 * p << "\n"; } }
        for (int i = 0; i < n; ++i) { for (int j = 0; j < n; ++j) { os << ((i == 2 && j == 2) ? top : 0.0) << "\n"; } }
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

// theta differs from cell to cell, so a row reading the wrong cell fails.
double tval (int i, int j, int k) { return 300 + 0.1 * i + 0.2 * j + 0.3 * k; }

// A case's buildings with their wall law at u = U along x, and their theta terms.
struct Set
{
    PCase c;
    gust::BuildingsOptions o;
    gust::Buildings b;
    gust::BuildingForcing law;
    gust::BuildingTheta heat;
    gust::ReferenceProfile ref;
    amrex::MultiFab ox, oy, oz, theta;
    Set (int mgs, const gust::BuildingsOptions& opt, double U, const gust::ReferenceProfile& profile, bool periodic = false)
        : c(mesh_case(mgs, periodic)), o(opt), ref(profile)
    {
        b.build(o, c.geom, c.z, c.ba, c.dm);
        law.build(b, o, c.geom, c.mesh.detJ);
        heat.build(b, law, o, c.geom, c.z, c.mesh.detJ, ref);
        amrex::MultiFab* op[3] = {&ox, &oy, &oz};
        for (int d = 0; d < 3; ++d) {
            op[d]->define(b.beta_face(d).boxArray(), c.dm, 1, 1);
            op[d]->setVal(1);
            amrex::MultiFab::LinComb(*op[d], 1, *op[d], 0, -1, b.beta_face(d), 0, 0, 1, 0);
            op[d]->FillBoundary(c.geom.periodicity());
        }
        amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1), w(c.mesh.az.boxArray(), c.dm, 1, 1);
        u.setVal(amrex::Real(U));
        v.setVal(0);
        w.setVal(0);
        law.update(u, v, w);
        theta.define(c.ba, c.dm, 1, 1);
        for (amrex::MFIter mfi(theta, false); mfi.isValid(); ++mfi) {
            const auto a4 = theta.array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) { a4(i, j, k) = static_cast<amrex::Real>(tval(i, j, k)); });
        }
    }
    [[nodiscard]] amrex::Array<const amrex::MultiFab*, 3> open () const { return {&ox, &oy, &oz}; }
};

gust::ReferenceProfile uniform_profile ()
{
    gust::ReferenceProfile p;
    p.z = {0};
    p.theta = {300};
    return p;
}

} // namespace

// UT-5.5.1: the exchange of a roof, of a side along the flow, of a half blocked side and of a side across the flow.
TEST (BuildingTheta, WallExchangeIsTheForm7ClosedForm)
{
    const double U = 3, tb = 305;
    for (const int mgs : {8, 4, 3}) {
        // A roof at 25 m: the cell (3, 3, 2) is half solid, its fluid centroid 2.5 m from the roof, its reference
        // cell (3, 3, 3) 10 m.
        Set s(mgs, block_options(25, gust::BuildingThetaMode::fixed, tb), U, uniform_profile());
        const double so = kKappa * U / std::log(25.0), sr = kKappa * U / std::log(100.0);
        const double us = std::sqrt(0.5 * so * so + 0.5 * sr * sr);
        const double own = 0.5 * kKappa * us / std::log(25.0) / 10, ref = 0.5 * kKappa * us / std::log(100.0) / 10;
        EXPECT_NEAR(value_at(s.law.rate(), 3, 3, 2, gust::kBuildingHeatOwn), own, 100 * rtol * own) << "boxes " << mgs;
        EXPECT_NEAR(value_at(s.law.rate(), 3, 3, 2, gust::kBuildingHeatRefZ), ref, 100 * rtol * ref);
        EXPECT_EQ(value_at(s.law.rate(), 3, 3, 2, gust::kBuildingHeatRefX), 0.0);
        amrex::MultiFab R(s.c.ba, s.c.dm, 1, 0);
        R.setVal(9);
        s.heat.apply_rows(s.theta, R);
        const double roof = own * (tval(3, 3, 2) - tb) + ref * (tval(3, 3, 3) - tb);
        EXPECT_NEAR(value_at(R, 3, 3, 2) - 9, roof, rowtol(roof)) << "boxes " << mgs;
        // The side at y = 50 m, u along it: a fluid cell, 5 m from the side, its own term alone.
        const double ss = kKappa * U / std::log(50.0), side = kKappa * ss / std::log(50.0) / 10;
        EXPECT_NEAR(value_at(s.law.rate(), 3, 5, 1, gust::kBuildingHeatOwn), side, 100 * rtol * side);
        EXPECT_EQ(value_at(s.law.rate(), 3, 5, 1, gust::kBuildingHeatRefY), 0.0);
        EXPECT_NEAR(value_at(R, 3, 5, 1) - 9, side * (tval(3, 5, 1) - tb), rowtol(side * (tval(3, 5, 1) - tb)));
        // The same side in the third layer, where the block's top leaves half the face blocked: half the area.
        EXPECT_NEAR(value_at(R, 3, 5, 2) - 9, 0.5 * side * (tval(3, 5, 2) - tb), rowtol(0.5 * side * (tval(3, 5, 2) - tb)));
        // A side the flow is normal to (x = 50 m): no parallel speed, so u* = 0 and no exchange.
        EXPECT_EQ(value_at(R, 5, 3, 1), 9.0);
        // A solid cell holds theta_b; a cell with no wall keeps its row.
        EXPECT_NEAR(value_at(R, 3, 3, 1), tval(3, 3, 1) - tb, ttol);
        EXPECT_EQ(value_at(R, 0, 0, 0), 9.0);
        EXPECT_EQ(value_at(R, 3, 3, 4), 9.0);
        EXPECT_EQ(value_at(s.heat.surface_theta(), 3, 5, 1), tb);
    }
}

// UT-5.5.2: with the reference profile as the walls' temperature a state at that profile exchanges nothing, also
// where the profile changes with height (a roof's reference cell is one layer up), and the solid rows are zero.
TEST (BuildingTheta, AmbientWallsExchangeNothingAtTheReferenceProfile)
{
    gust::ReferenceProfile stratified;
    stratified.z = {0, 240};
    stratified.theta = {300, 312};
    Set s(3, block_options(25, gust::BuildingThetaMode::ambient), 3, stratified);
    amrex::MultiFab at_ref(s.c.ba, s.c.dm, 1, 1), R(s.c.ba, s.c.dm, 1, 0);
    gust::fill_theta_ref(stratified, s.c.z, at_ref);
    EXPECT_NEAR(value_at(s.heat.surface_theta(), 3, 3, 2), 300 + 12 * 25.0 / 240, ttol) << "the profile at the cell's height";
    ASSERT_GT(value_at(s.law.rate(), 3, 3, 2, gust::kBuildingHeatRefZ), 0.0);
    R.setVal(9);
    s.heat.apply_rows(at_ref, R);
    EXPECT_EQ(value_at(R, 3, 3, 2), 9.0) << "a roof's wall cell and its reference cell, each against its own height's";
    EXPECT_EQ(value_at(R, 3, 5, 1), 9.0);
    EXPECT_EQ(value_at(R, 3, 3, 1), 0.0) << "a solid cell at the profile";
    EXPECT_EQ(static_cast<double>(s.heat.wall_heat(at_ref)), 0.0);
    // Away from the profile the walls act: 1 K warmer air gives heat to them.
    at_ref.plus(amrex::Real(1), 0, 1, 1);
    R.setVal(9);
    s.heat.apply_rows(at_ref, R);
    EXPECT_GT(value_at(R, 3, 3, 2), 9.0);
    EXPECT_NEAR(value_at(R, 3, 3, 1), 1.0, ttol);
}

// UT-5.5.3: the implicit rows. a_P gains the derivative of the exchange to the cell's own theta (a finite
// difference of the row), solid cells are identity rows, and the reference cell's part stays in the residual.
TEST (BuildingTheta, ImplicitDiagonalIsTheDerivativeOfTheExchange)
{
    const double tb = 305;
    Set s(4, block_options(25, gust::BuildingThetaMode::fixed, tb), 3, uniform_profile());
    amrex::MultiFab A(s.c.ba, s.c.dm, gust::n_implicit_coefs, 0), R0(s.c.ba, s.c.dm, 1, 0), R1(s.c.ba, s.c.dm, 1, 0);
    A.setVal(2);
    s.heat.add_to_implicit(s.open(), A);
    R0.setVal(0);
    s.heat.apply_rows(s.theta, R0);
    amrex::MultiFab warmer(s.c.ba, s.c.dm, 1, 1);
    amrex::MultiFab::Copy(warmer, s.theta, 0, 0, 1, 1);
    const double eps = 0.5;
    for (const auto& cell : {amrex::IntVect(3, 3, 2), amrex::IntVect(3, 5, 1), amrex::IntVect(3, 5, 2)}) {
        for (amrex::MFIter mfi(warmer, false); mfi.isValid(); ++mfi) {
            if (!mfi.fabbox().contains(cell)) { continue; }
            const auto a4 = warmer.array(mfi);
            amrex::LoopOnCpu(amrex::Box(cell, cell), [&] (int i, int j, int k) { a4(i, j, k) = static_cast<amrex::Real>(tval(i, j, k) + eps); });
        }
    }
    R1.setVal(0);
    s.heat.apply_rows(warmer, R1);
    for (const auto& cell : {amrex::IntVect(3, 3, 2), amrex::IntVect(3, 5, 1), amrex::IntVect(3, 5, 2)}) {
        const int i = cell[0], j = cell[1], k = cell[2];
        const double dRdT = (value_at(R1, i, j, k) - value_at(R0, i, j, k)) / eps;
        EXPECT_GT(dRdT, 1e-4);
        // The row is linear in theta, so the difference is the derivative to round-off (single precision: of 300 K).
        const bool dbl = std::is_same_v<amrex::Real, double>;
        const double tol = (dbl ? 1e-9 : 0.05) * dRdT;
        EXPECT_NEAR(value_at(A, i, j, k, gust::coef_P) - 2, dRdT, tol) << cell;
    }
    // The roof's cell: its face on the roof is closed, so the coupling below is gone; a_P keeps what it had.
    EXPECT_EQ(value_at(A, 3, 3, 2, gust::coef_ZM), 0.0);
    EXPECT_EQ(value_at(A, 3, 3, 2, gust::coef_ZP), 2.0) << "the reference cell's part is not in the rows";
    EXPECT_EQ(value_at(A, 3, 3, 1, gust::coef_P), 1.0) << "a solid cell: the identity";
    for (const int n : {gust::coef_XM, gust::coef_XP, gust::coef_YM, gust::coef_YP, gust::coef_ZM, gust::coef_ZP}) {
        EXPECT_EQ(value_at(A, 3, 3, 1, n), 0.0);
    }
    EXPECT_EQ(value_at(A, 0, 0, 0, gust::coef_P), 2.0);
}

// UT-5.5.4: no heat is conducted out of a building. The diffusive flux through a blocked face is zero, half of it
// passes a half blocked face, and the implicit couplings follow.
TEST (BuildingTheta, DiffusiveFluxOfThetaIsClosedOnBlockedFaces)
{
    for (const int mgs : {8, 3}) {
        Set s(mgs, block_options(25, gust::BuildingThetaMode::fixed, 305), 3, uniform_profile());
        amrex::Array<amrex::MultiFab, 3> G;
        for (int d = 0; d < 3; ++d) {
            G[static_cast<std::size_t>(d)].define(s.b.beta_face(d).boxArray(), s.c.dm, 1, 0);
            G[static_cast<std::size_t>(d)].setVal(7);
        }
        s.heat.close_fluxes(s.open(), G);
        EXPECT_EQ(value_at(G[1], 3, 5, 1), 0.0) << "the side at y = 50 m";
        EXPECT_EQ(value_at(G[1], 3, 5, 2), 3.5) << "half blocked";
        EXPECT_EQ(value_at(G[1], 3, 5, 3), 7.0) << "above the block";
        EXPECT_EQ(value_at(G[0], 5, 3, 0), 0.0) << "the side at x = 50 m";
        EXPECT_EQ(value_at(G[2], 3, 3, 3), 7.0);
        EXPECT_EQ(value_at(G[2], 3, 3, 2), 0.0) << "under the roof's cell: inside the block";
        EXPECT_EQ(value_at(G[0], 1, 1, 1), 7.0);
        amrex::MultiFab A(s.c.ba, s.c.dm, gust::n_implicit_coefs, 0);
        A.setVal(2);
        s.heat.add_to_implicit(s.open(), A);
        EXPECT_EQ(value_at(A, 3, 5, 1, gust::coef_YM), 0.0);
        EXPECT_EQ(value_at(A, 3, 5, 2, gust::coef_YM), 1.0);
        EXPECT_EQ(value_at(A, 3, 5, 1, gust::coef_YP), 2.0);
    }
}

// The heat the walls take from the air, summed over the wall cells: the cells' exchange times their volume, the
// same on every layout. With the air 2 K above the walls everywhere it is 2 K times the sum of the coefficients.
TEST (BuildingTheta, WallHeatIsTheSumOfTheCellsExchange)
{
    double first = 0;
    for (const int mgs : {8, 4, 3}) {
        Set s(mgs, block_options(25, gust::BuildingThetaMode::fixed, 305), 3, uniform_profile());
        amrex::MultiFab warm(s.c.ba, s.c.dm, 1, 1);
        warm.setVal(307);
        double coef = 0;
        for (amrex::MFIter mfi(s.law.rate(), false); mfi.isValid(); ++mfi) {
            const auto r = s.law.rate().const_array(mfi);
            const auto beta = s.b.beta().const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (beta(i, j, k) >= 1) { return; }
                for (int n = gust::kBuildingHeatOwn; n <= gust::kBuildingHeatRefZ; ++n) { coef += static_cast<double>(r(i, j, k, n)); }
            });
        }
        amrex::ParallelDescriptor::ReduceRealSum(coef);
        ASSERT_GT(coef, 0.01);
        const double q = static_cast<double>(s.heat.wall_heat(warm));
        EXPECT_NEAR(q, 2 * coef * 1000, 1e3 * rtol * q) << "boxes " << mgs;
        if (mgs == 8) { first = q; }
        EXPECT_NEAR(q, first, 100 * rtol * q);
    }
}

// The scaled norm of R_theta counts the solid cells as zeros, as that of R_k does.
TEST (BuildingTheta, NormOfThetaLeavesTheSolidCellsOut)
{
    Set s(3, block_options(45, gust::BuildingThetaMode::fixed, 305), 3, uniform_profile());
    gust::ResidualFields R;
    R.define(s.c.ba, s.c.dm);
    for (amrex::MultiFab* f : {&R.u, &R.v, &R.w, &R.p, &R.tke}) { f->setVal(2); }
    R.theta.setVal(3);
    gust::State q;
    q.define(s.c.ba, s.c.dm);
    for (amrex::MultiFab* f : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { f->setVal(1.0); }
    const gust::JvpScales sc;
    const double st = static_cast<double>(sc.velocity * sc.theta);   // U theta_s / H with H = 1
    std::vector<amrex::Real> all, fluid, mon;
    amrex::MultiFab work;
    gust::stage_norms(R, q, nullptr, sc, 1, s.c.geom, all, mon);
    gust::stage_norms(R, q, nullptr, sc, 1, s.c.geom, fluid, mon, nullptr, nullptr, &s.heat, &work);
    const double solid = 16.0 / (8 * 8 * 24);   // 2 x 2 columns, 4 solid layers
    EXPECT_NEAR(static_cast<double>(all[4]), 3 / st, 100 * rtol * 3 / st);
    EXPECT_NEAR(static_cast<double>(fluid[4]), 3 * std::sqrt(1 - solid) / st, 100 * rtol * 3 / st);
    for (const int n : {0, 1, 2, 3, 5}) { EXPECT_EQ(fluid[static_cast<std::size_t>(n)], all[static_cast<std::size_t>(n)]) << n; }
    EXPECT_EQ(value_at(R.theta, 3, 3, 1), 3.0) << "R_theta itself is not changed";
}

// The same in the count of steady.confidence (SimplecStage::confidence) and in the norm the stage reports: after
// one iteration from air at rest at 300 K with walls at 305 K, the residual of theta is theta - theta_b = -5 K in
// the solid cells and zero elsewhere (no flow, so no exchange and no convection). Every cell meets the criterion
// and the norm is zero. (With a wind the first state is not divergence-free around the block, and the residual
// of theta in the cells beside it holds theta times that divergence, which hides the solid cells' part.)
TEST (BuildingTheta, ConfidenceLeavesTheSolidCellsOut)
{
    Set s(4, block_options(45, gust::BuildingThetaMode::fixed, 305, 4), 3, uniform_profile(), true);
    s.c.mesh.set_open_faces(s.b.beta_face(0), s.b.beta_face(1), s.b.beta_face(2), s.b.beta(), s.c.geom);
    gust::State q;
    q.define(s.c.ba, s.c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, s.c.z, 0, 0, 1, s.c.z.nGrow());
    q.u.setVal(0);
    q.v.setVal(0);
    q.w.setVal(0);
    q.p.setVal(0);
    q.theta.setVal(300);
    q.tke.setVal(1);
    amrex::MultiFab nu(s.c.ba, s.c.dm, 1, 1);
    nu.setVal(0.5);
    using gust::FaceType;
    gust::BoundaryInputs bin;
    bin.face = {FaceType::periodic, FaceType::periodic, FaceType::slip_wall, FaceType::periodic, FaceType::periodic, FaceType::slip_wall};
    gust::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.theta_ref = uniform_profile();
    in.buoyancy.gravity = 0;
    in.buildings = &s.law;
    in.building_theta = &s.heat;
    const gust::BoundaryConditions bc(bin, s.c.geom, s.c.mesh, in.coriolis);
    in.bc = &bc;
    gust::PoissonOptions popt;
    popt.reltol = amrex::Real(0.1);
    gust::SimplecOptions opt;
    gust::SimplecStage stage(opt, 2, in, nullptr, popt, s.c.mesh, s.c.geom, s.c.z, 0, s.c.ba, s.c.dm);
    stage.init(q);
    bc.fill_ghosts(q, nullptr, in.buoyancy);
    amrex::MultiFab dtau;
    gust::local_dtau(q, &stage.nu(), s.c.mesh, s.c.geom, amrex::Real(1), dtau);
    stage.iterate(q, dtau);
    EXPECT_NEAR(value_at(stage.residual_fields().theta, 3, 3, 1), -5.0, ttol) << "a solid cell's row";
    // |R_theta| <= 2.5 K/s as a scaled threshold; the other fields always pass.
    const gust::JvpScales sc;
    const double st = static_cast<double>(sc.velocity * sc.theta) / static_cast<double>(stage.norm_height());
    const std::vector<amrex::Real> threshold = {amrex::Real(1e20), amrex::Real(1e20), amrex::Real(1e20), amrex::Real(1e20),
                                                static_cast<amrex::Real>(2.5 / st), amrex::Real(1e20)};
    amrex::Real share = 0, height = -1;
    ASSERT_TRUE(stage.confidence(threshold, share, height));
    EXPECT_EQ(share, amrex::Real(1)) << "the 16 solid cells, with |R_theta| = 5, do not count against the criterion";
    EXPECT_EQ(height, amrex::Real(0));
    // The norm leaves them out too: they alone would give 5 sqrt(16 / 1536) / scale.
    EXPECT_LT(static_cast<double>(stage.residual_norms()[4]), 1e-6 * 5 * std::sqrt(16.0 / 1536) / st);
}

// Through residual(): the rows are applied, and the heat does not leave the block by diffusion.
TEST (BuildingTheta, ResidualCarriesTheRows)
{
    const double tb = 305;
    Set s(3, block_options(25, gust::BuildingThetaMode::fixed, tb), 3, uniform_profile());
    s.c.mesh.set_open_faces(s.b.beta_face(0), s.b.beta_face(1), s.b.beta_face(2), s.b.beta(), s.c.geom);
    gust::State q;
    q.define(s.c.ba, s.c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, s.c.z, 0, 0, 1, s.c.z.nGrow());
    q.u.setVal(3);
    q.v.setVal(0);
    q.w.setVal(0);
    q.p.setVal(0);
    q.tke.setVal(1);
    q.theta.setVal(300);
    amrex::MultiFab nu(s.c.ba, s.c.dm, 1, 1);
    nu.setVal(0.5);
    gust::ResidualInputs in;
    in.nu = &nu;
    in.gamma_theta = &nu;
    in.gamma_tke = &nu;
    in.theta_ref = uniform_profile();
    in.buildings = &s.law;
    gust::ResidualFields R0, R1;
    gust::residual(q, in, s.c.mesh, s.c.geom, R0);
    in.building_theta = &s.heat;
    gust::residual(q, in, s.c.mesh, s.c.geom, R1);
    // Uniform theta: no diffusion anywhere, so the rows are the whole difference.
    const double ss = kKappa * 3 / std::log(50.0), side = kKappa * ss / std::log(50.0) / 10;
    EXPECT_NEAR(value_at(R1.theta, 3, 5, 1) - value_at(R0.theta, 3, 5, 1), side * (300 - tb), rowtol(side * (300 - tb)));
    EXPECT_NEAR(value_at(R1.theta, 3, 3, 1), 300 - tb, ttol) << "a solid cell";
    EXPECT_EQ(value_at(R1.theta, 0, 0, 5), value_at(R0.theta, 0, 0, 5));
    EXPECT_EQ(value_at(R1.tke, 3, 5, 1), value_at(R0.tke, 3, 5, 1)) << "k is not touched";
}

// buildings.theta_mode and buildings.theta_surface: defaults, refusals and the restart hash.
TEST (BuildingTheta, OptionsDefaultsAndRefusals)
{
    gust::BuildingsOptions o = block_options(45, gust::BuildingThetaMode::none);
    EXPECT_TRUE(gust::check_buildings_options(o).empty());
    const std::string none = gust::Buildings::hash_inputs(o);
    o.theta_mode = gust::BuildingThetaMode::ambient;
    const std::string ambient = gust::Buildings::hash_inputs(o);
    EXPECT_NE(ambient, none) << "a checkpoint from before theta's terms continues with none alone";
    // fixed needs its temperature, finite and > 0, and the temperature is part of the hash.
    o.theta_mode = gust::BuildingThetaMode::fixed;
    o.theta_mode_name = "fixed";
    auto err = gust::check_buildings_options(o);
    ASSERT_EQ(err.size(), 1U);
    EXPECT_NE(err[0].find("buildings.theta_surface is required with buildings.theta_mode = fixed"), std::string::npos);
    o.given.emplace_back("theta_surface");
    o.theta_surface = 0;
    err = gust::check_buildings_options(o);
    ASSERT_EQ(err.size(), 1U);
    EXPECT_NE(err[0].find("buildings.theta_surface must be finite and > 0"), std::string::npos);
    o.theta_surface = 305;
    EXPECT_TRUE(gust::check_buildings_options(o).empty());
    const std::string fixed305 = gust::Buildings::hash_inputs(o);
    o.theta_surface = 306;
    EXPECT_NE(gust::Buildings::hash_inputs(o), fixed305);
    EXPECT_NE(fixed305, ambient);
    // The temperature without fixed is refused, as a key that would be ignored.
    o.theta_mode = gust::BuildingThetaMode::ambient;
    o.theta_mode_name = "ambient";
    err = gust::check_buildings_options(o);
    ASSERT_EQ(err.size(), 1U);
    EXPECT_NE(err[0].find("buildings.theta_surface has no effect without buildings.theta_mode = fixed"), std::string::npos);
    // The exchange reads the wall law's u* on every wall.
    o.given.clear();
    o.wall_form = gust::BuildingWallForm::roof_stress;
    o.wall_form_name = "roof_stress";
    o.wall_tke = false;
    err = gust::check_buildings_options(o);
    ASSERT_EQ(err.size(), 1U);
    EXPECT_NE(err[0].find("buildings.theta_mode = ambient needs buildings.wall_form = fraction_stress"), std::string::npos);
    o.theta_mode_name = "warm";
    err = gust::check_buildings_options(o);
    ASSERT_EQ(err.size(), 1U);
    EXPECT_NE(err[0].find("buildings.theta_mode must be none, ambient or fixed (got warm)"), std::string::npos);
    // The reader's defaults: ambient with the law on every wall, none with the other forms.
    {
        amrex::ParmParse pp("buildings");
        pp.add("type", std::string("immersed_forcing"));
        pp.add("file", o.file);
        const gust::BuildingsOptions read = gust::read_buildings_options();
        EXPECT_TRUE(read.theta_mode == gust::BuildingThetaMode::ambient);
        pp.add("wall_form", std::string("roof_stress"));
        EXPECT_TRUE(gust::read_buildings_options().theta_mode == gust::BuildingThetaMode::none);
        pp.remove("wall_form");
        pp.add("theta_mode", std::string("fixed"));
        pp.add("theta_surface", 303.5);
        const gust::BuildingsOptions fixed = gust::read_buildings_options();
        EXPECT_TRUE(fixed.theta_mode == gust::BuildingThetaMode::fixed);
        EXPECT_EQ(fixed.theta_surface, 303.5);
        EXPECT_TRUE(gust::check_buildings_options(fixed).empty());
        for (const char* k : {"type", "file", "theta_mode", "theta_surface"}) { pp.remove(k); }
    }
}
