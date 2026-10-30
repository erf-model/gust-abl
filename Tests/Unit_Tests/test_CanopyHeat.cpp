// The leaves' heat exchange (CanopyHeat.H), plan rev. 6 beyond the plan,
// sub-phase 8.4 (UT-8.4.x):
//   the leaf area density is the drag field's cell means with Cd = 1, a patch
//   with Cd = 0 still has leaves, and the coarse levels keep each layer's
//   integral of a (UT-8.4.0, beyond the plan);
//   the source lambda (theta_l - theta), lambda = 2 a / r_b, against its closed
//   form with ERF's constants on a sheared, stratified field, the wind floor
//   included, and the pinned value lambda = a / 5 at |U| = 5 m/s (UT-8.4.1);
//   the theta rows' diagonal is lambda, the derivative of the residual by
//   central differences, and the residual gains exactly lambda (theta - theta_l)
//   (UT-8.4.2);
//   no leaf heat leaves the residual and the theta rows bit for bit (UT-8.4.3);
//   the frozen-theta blend: theta_f exactly outside the canopy, the closed-form
//   weight inside, the midpoint where lambda k = eps (UT-8.4.4, beyond the
//   plan);
//   the option checks name the input.
#include "Canopy.H"
#include "CanopyHeat.H"
#include "Closure.H"
#include "ImplicitOperator.H"
#include "MassFlux.H"
#include "MeshHierarchy.H"
#include "PoissonTestHelpers.H"
#include "Reductions.H"
#include "Residual.H"
#include "TestHelpers.H"

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <type_traits>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;

// Two patches: one at (200, 200), diameter 300, height 17.3 m, Cd 0.25, LAI 4
// (uniform); one at (60, 340), diameter 80, height 9 m, Cd 0, LAI 3
// (Lalic–Mihailovic, laimax 0.7): no drag, but leaves.
constexpr double kCd = 0.25, kLai = 4.0, kH = 17.3, kThetaLeaf = 303.5;

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

std::string write_file (const std::string& name, const std::function<void(std::ofstream&)>& body)
{
    const fs::path dir = fs::temp_directory_path() / ("gust_canopy_heat_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / name).string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        body(os);
    }
    amrex::ParallelDescriptor::Barrier();
    return file;
}

gust::CanopyOptions patch_options (bool heat, gust::CanopyFrozenHeat frozen = gust::CanopyFrozenHeat::skip)
{
    const std::string file = write_file("patches.txt", [] (std::ofstream& os) {
        os << "1 200 200 " << kH << " 300 " << kCd << " " << kLai << " 0.8\n";
        os << "2 60 340 9 80 0 3 0.7\n";
    });
    gust::CanopyOptions o;
    o.type = gust::CanopyType::patches;
    o.type_name = "patches";
    o.file = file;
    o.given = {"type", "file"};
    if (heat) {
        o.heat.name = "fixed";
        o.heat.mode = gust::CanopyHeatMode::fixed;
        o.heat.leaf_theta = kThetaLeaf;
        o.heat.frozen = frozen;
        o.heat.frozen_name = (frozen == gust::CanopyFrozenHeat::blend) ? "blend" : "skip";
        o.given.insert(o.given.end(), {"heat", "leaf_theta"});
    }
    return o;
}

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

amrex::FArrayBox host (const amrex::MultiFab& mf, const amrex::MFIter& mfi)
{
    const amrex::Box b = mf[mfi].box();
    amrex::FArrayBox h(b, mf.nComp(), amrex::The_Pinned_Arena());
    h.copy<amrex::RunOn::Device>(mf[mfi], b, 0, b, 0, mf.nComp());
    amrex::Gpu::streamSynchronize();
    return h;
}

constexpr double kw = 2 * M_PI / 16;

// A periodic, sheared, stably stratified state; the wind is below the 0.01
// m/s floor of the leaf resistance in a band of x.
void fill_state (gust::State& q, const PCase& c)
{
    q.define(c.ba, c.dm);
    amrex::MultiFab::Copy(q.z_phys_nd, c.z, 0, 0, 1, c.z.nGrow());
    auto slow = [] (double X) { return (X >= 3 && X < 5) ? 1e-3 : 1.0; };
    fill_component(q.u, 0, [slow] (double X, double Y, double Z) { return slow(X) * (2 + 0.3 * Z + 0.5 * std::sin(kw * Y)); });
    fill_component(q.v, 1, [slow] (double X, double, double Z) { return slow(X) * (0.4 * std::cos(kw * X) + 0.05 * Z); });
    fill_component(q.w, 2, [slow] (double X, double Y, double) { return slow(X) * 0.05 * std::sin(kw * (X + Y)); });
    q.p.setVal(0.0);
    fill_component(q.theta, 3, [] (double X, double, double Z) { return 300 + 0.5 * Z + 0.2 * std::sin(kw * X); });
    fill_component(q.tke, 3, [] (double X, double Y, double) { return 0.5 * (1 + 0.3 * std::sin(kw * X) * std::cos(kw * Y)); });
}

// lambda by hand: ERF's r_b = 100 sqrt(0.05 / max(|U|, 0.01)), lambda = 2 a / r_b.
double lambda_by_hand (double a, double U)
{
    return 2 * a / (100 * std::sqrt(0.05 / std::max(U, 0.01)));
}

double cell_wind (const amrex::FArrayBox& u, const amrex::FArrayBox& v, const amrex::FArrayBox& w, int i, int j, int k)
{
    const double uc = 0.5 * (u(amrex::IntVect(i, j, k)) + u(amrex::IntVect(i + 1, j, k)));
    const double vc = 0.5 * (v(amrex::IntVect(i, j, k)) + v(amrex::IntVect(i, j + 1, k)));
    const double wc = 0.5 * (w(amrex::IntVect(i, j, k)) + w(amrex::IntVect(i, j, k + 1)));
    return std::sqrt(uc * uc + vc * vc + wc * wc);
}

struct HeatSetup
{
    PCase c;
    gust::Canopy canopy;
    gust::CanopyHeat heat;
    gust::State q;
    amrex::MultiFab dist;
    explicit HeatSetup (int mgs = 8, gust::CanopyFrozenHeat frozen = gust::CanopyFrozenHeat::skip) : c(make_case(mgs))
    {
        canopy.build(patch_options(true, frozen), c.geom, c.z, c.ba, c.dm);
        heat.build(canopy);
        fill_state(q, c);
        gust::height_above_ground(q.z_phys_nd, dist);
    }
};

std::vector<double> layer_integrals (const amrex::MultiFab& f, const amrex::MultiFab& detJ, const amrex::Geometry& geom)
{
    const amrex::Box& dom = geom.Domain();
    std::vector<double> out(static_cast<std::size_t>(dom.length(2)), 0.0);
    const double dxdy = static_cast<double>(geom.CellSize(0)) * static_cast<double>(geom.CellSize(1));
    const double dzeta = static_cast<double>(geom.CellSize(2));
    for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
        const auto ff = host(f, mfi), J = host(detJ, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            out[static_cast<std::size_t>(k - dom.smallEnd(2))] += ff(iv) * J(iv) * dxdy * dzeta;
        });
    }
    amrex::ParallelDescriptor::ReduceRealSum(out.data(), static_cast<int>(out.size()));
    return out;
}

} // namespace

// UT-8.4.0: a = f / Cd where Cd > 0 (the same cell means), a > 0 under the
// Cd = 0 patch where f = 0, a = 0 above the patches; without canopy.heat the
// field is not built.
TEST(CanopyHeat, LeafAreaDensityIsTheDragWithUnitCd)
{
    PCase c = make_case(8);
    gust::Canopy with, without;
    with.build(patch_options(true), c.geom, c.z, c.ba, c.dm);
    without.build(patch_options(false), c.geom, c.z, c.ba, c.dm);
    EXPECT_TRUE(without.lad().empty());
    ASSERT_FALSE(with.lad().empty());
    EXPECT_EQ(gust::count_bit_differences(with.drag(), without.drag(), amrex::IntVect(0)), 0);
    double worst = 0;
    amrex::Long leaves_no_drag = 0, crowns = 0, wrong_above = 0;
    for (amrex::MFIter mfi(with.lad(), false); mfi.isValid(); ++mfi) {
        const auto a = host(with.lad(), mfi), f = host(with.drag(), mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (f(iv) > 0) {
                worst = std::max(worst, std::abs(kCd * a(iv) - f(iv)) / f(iv));
                if (std::abs(a(iv) - kLai / kH) < (kDouble ? 1e-12 : 1e-5) * kLai / kH) { ++crowns; }
            } else if (a(iv) > 0) {
                ++leaves_no_drag;
            }
            if (k >= 20 && a(iv) != 0) { ++wrong_above; }
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    for (amrex::Long* n : {&leaves_no_drag, &crowns, &wrong_above}) { amrex::ParallelDescriptor::ReduceLongSum(*n); }
    EXPECT_LE(worst, kDouble ? 1e-14 : 1e-6);
    EXPECT_GT(leaves_no_drag, 0L);
    EXPECT_GT(crowns, 100L);
    EXPECT_EQ(wrong_above, 0L);
}

// UT-8.4.0, canopy.type = map with a Cd map: a = f / Cd column by column
// (the Cd map is replaced by 1, not the height or the LAI), with the
// Lalic–Mihailovic profile.
TEST(CanopyHeat, LeafAreaDensityOfAMap)
{
    PCase c = make_case(8);
    auto points = [] (const std::string& name, const std::function<double(int, int)>& v) {
        return write_file(name, [&] (std::ofstream& os) {
            for (int j = 0; j < 9; ++j) {
                for (int i = 0; i < 9; ++i) { os << 50 * i << " " << 50 * j << " " << v(i, j) << "\n"; }
            }
        });
    };
    gust::CanopyOptions o;
    o.type = gust::CanopyType::map;
    o.type_name = "map";
    o.height_file = points("h.xyz", [] (int i, int j) { return 10.0 + i + 0.5 * j; });
    o.lai_file = points("lai.xyz", [] (int i, int) { return 2.0 + 0.25 * i; });
    o.cd_file = points("cd.xyz", [] (int, int j) { return 0.1 + 0.02 * j; });
    o.profile = gust::CanopyProfile::lalic_mihailovic;
    o.profile_name = "lalic_mihailovic";
    o.laimax = 0.7;
    o.idw_neighbors = 4;
    o.given = {"type", "height_file", "lai_file", "cd_file", "profile", "laimax", "heat", "leaf_theta"};
    o.heat.name = "fixed";
    o.heat.mode = gust::CanopyHeatMode::fixed;
    o.heat.leaf_theta = kThetaLeaf;
    gust::Canopy cn;
    cn.build(o, c.geom, c.z, c.ba, c.dm);
    // The column's Cd by the same IDW: f / a in the canopy cells of one column is one value, and it varies over y.
    double spread = 0, cd_lo = 1e30, cd_hi = 0;
    amrex::Long cells = 0;
    for (amrex::MFIter mfi(cn.lad(), false); mfi.isValid(); ++mfi) {
        const auto a = host(cn.lad(), mfi), f = host(cn.drag(), mfi);
        const amrex::Box b = mfi.validbox();
        for (int j = b.smallEnd(1); j <= b.bigEnd(1); ++j) {
            for (int i = b.smallEnd(0); i <= b.bigEnd(0); ++i) {
                double first = -1;
                for (int k = b.smallEnd(2); k <= b.bigEnd(2); ++k) {
                    const amrex::IntVect iv(i, j, k);
                    // Cells with a share of the canopy too small for the precision's digits (a part-filled top
                    // cell) are left out.
                    if (!(a(iv) > (kDouble ? 0.0 : 1e-4))) { continue; }
                    const double r = f(iv) / a(iv);
                    if (first < 0) { first = r; }
                    spread = std::max(spread, std::abs(r - first) / first);
                    cd_lo = std::min(cd_lo, r);
                    cd_hi = std::max(cd_hi, r);
                    ++cells;
                }
            }
        }
    }
    amrex::ParallelDescriptor::ReduceRealMax(spread);
    amrex::ParallelDescriptor::ReduceRealMax(cd_hi);
    amrex::ParallelDescriptor::ReduceRealMin(cd_lo);
    amrex::ParallelDescriptor::ReduceLongSum(cells);
    EXPECT_GT(cells, 1000L);
    EXPECT_LE(spread, kDouble ? 1e-13 : 1e-5);
    EXPECT_GT(cd_lo, 0.1 - 1e-9);
    EXPECT_LT(cd_hi, 0.26 + 1e-9);
    EXPECT_GT(cd_hi - cd_lo, 0.05);
}

// UT-8.4.0: the coarse levels of the grid hierarchy keep each layer's
// integral of a (as of f), and take the fine level's heat options.
TEST(CanopyHeat, RestrictionKeepsLayerIntegralsOfTheLeafArea)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(32, 16, 12);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {0, 1, 0};
    c.in.mesh.max_grid_size = {8, 8};
    c.in.mesh.dz0 = amrex::Real(2);
    c.in.mesh.stretch_ratio = amrex::Real(1.3);
    c.in.mesh.stretch_ratio_given = true;
    c.in.mesh.hierarchy_levels = 3;
    c.in.terrain.type = "agnesi";
    c.in.terrain.hill_shape = "round";
    c.in.terrain.hill_height = 60;
    c.in.terrain.hill_half_width = 80;
    c.in.terrain.smoothing = 1;
    c.in.terrain.hierarchy_passes = 2;
    c.geom = gust::make_geometry(c.in.mesh);
    c.ba = gust::hierarchy_box_array(c.in.mesh);
    c.dm = amrex::DistributionMapping(c.ba);
    c.z.define(amrex::convert(c.ba, amrex::IntVect(1)), c.dm, 1, gust::ngrow_nodal);
    c.z.setVal(0.0);
    c.mesh.define(c.ba, c.dm);
    c.mesh.build(c.in, c.geom, c.z);
    gust::MeshHierarchy h(c.in, c.geom, c.z, c.mesh, 0);
    gust::Canopy fine;
    fine.build(patch_options(true, gust::CanopyFrozenHeat::blend), c.geom, c.z, c.ba, c.dm);
    const std::vector<double> want = layer_integrals(fine.lad(), c.mesh.detJ, c.geom);
    double total = 0;
    for (const double v : want) { total += v; }
    ASSERT_GT(total, 0.0);
    for (int l = 1; l < 3; ++l) {
        gust::Canopy coarse;
        coarse.restrict_from(fine, c.mesh.detJ, h.geom(l), h.mesh(l), h.boxArray(l), h.DistributionMap(l));
        EXPECT_TRUE(coarse.heat().on());
        EXPECT_EQ(coarse.heat().frozen, gust::CanopyFrozenHeat::blend);
        EXPECT_EQ(coarse.heat().leaf_theta, kThetaLeaf);
        const std::vector<double> got = layer_integrals(coarse.lad(), h.mesh(l).detJ, h.geom(l));
        ASSERT_EQ(got.size(), want.size());
        for (std::size_t k = 0; k < want.size(); ++k) {
            EXPECT_NEAR(got[k], want[k], (kDouble ? 1e-12 : 1e-5) * total) << "level " << l << " layer " << k;
        }
    }
}

// UT-8.4.1: the pinned value: at |U| = 5 m/s, r_b = 100 sqrt(0.05 / 5) = 10 s/m
// and lambda = a / 5; below 0.01 m/s, r_b = 100 sqrt(5) s/m.
TEST(CanopyHeat, ExchangeRatePinnedValues)
{
    const double t = kDouble ? 1e-15 : 1e-6;
    EXPECT_NEAR(static_cast<double>(gust::leaf_exchange_rate(amrex::Real(0.4), amrex::Real(5))), 0.08, t);
    EXPECT_NEAR(static_cast<double>(gust::leaf_exchange_rate(amrex::Real(0.4), amrex::Real(1e-4))), 0.8 / (100 * std::sqrt(5.0)), t);
    EXPECT_EQ(static_cast<double>(gust::leaf_exchange_rate(amrex::Real(0), amrex::Real(3))), 0.0);
}

// UT-8.4.1: R_theta from add_canopy_heat alone is lambda (theta - theta_l),
// lambda from the hand formula at each cell's own wind, on two layouts, the
// slow band (the wind floor) included.
TEST(CanopyHeat, SourceHasItsClosedForm)
{
    for (const int mgs : {8, 16}) {
        HeatSetup s(mgs);
        amrex::MultiFab R(s.c.ba, s.c.dm, 1, 0);
        R.setVal(0.0);
        gust::add_canopy_heat(s.heat, s.q, R);
        double worst = 0, scale = 0;
        amrex::Long floored = 0, heated = 0;
        for (amrex::MFIter mfi(R, false); mfi.isValid(); ++mfi) {
            const auto r = host(R, mfi), a = host(s.canopy.lad(), mfi), th = host(s.q.theta, mfi);
            const auto u = host(s.q.u, mfi), v = host(s.q.v, mfi), w = host(s.q.w, mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const amrex::IntVect iv(i, j, k);
                const double U = cell_wind(u, v, w, i, j, k);
                const double want = lambda_by_hand(a(iv), U) * (th(iv) - kThetaLeaf);
                worst = std::max(worst, std::abs(r(iv) - want));
                scale = std::max(scale, std::abs(want));
                if (a(iv) > 0) {
                    ++heated;
                    if (U < 0.01) { ++floored; }
                }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceRealMax(scale);
        amrex::ParallelDescriptor::ReduceLongSum(floored);
        amrex::ParallelDescriptor::ReduceLongSum(heated);
        EXPECT_LE(worst, (kDouble ? 1e-13 : 1e-5) * scale) << "max_grid_size " << mgs;
        EXPECT_GT(heated, 100L);
        EXPECT_GT(floored, 10L);
    }
}

// UT-8.4.2 and UT-8.4.3: the residual with the leaf heat minus without it is
// lambda (theta - theta_l); its derivative in theta by central differences is
// the theta rows' added diagonal lambda; without in.canopy_heat the theta
// rows and the residual are those of a canopy run without heat, bit for bit.
TEST(CanopyHeat, ResidualAndDiagonal)
{
    HeatSetup s(8);
    gust::Closure cl(gust::ClosureOptions{});
    gust::CanopyFaceDrag faces;
    faces.build(s.canopy, s.c.geom);
    gust::ResidualInputs in0;
    in0.closure = &cl;
    in0.wall_distance = &s.dist;
    in0.theta_ref.z = {0};
    in0.theta_ref.theta = {300};
    in0.canopy = &faces;
    gust::ResidualInputs in1 = in0;
    in1.canopy_heat = &s.heat;
    gust::ResidualFields R0, R1;
    gust::residual(s.q, in0, s.c.mesh, s.c.geom, R0);
    gust::residual(s.q, in1, s.c.mesh, s.c.geom, R1);
    for (const auto& p : {std::make_pair(&R0.u, &R1.u), std::make_pair(&R0.v, &R1.v), std::make_pair(&R0.w, &R1.w),
                    std::make_pair(&R0.p, &R1.p), std::make_pair(&R0.tke, &R1.tke)}) {
        EXPECT_EQ(gust::count_bit_differences(*p.first, *p.second, amrex::IntVect(0)), 0);
    }
    amrex::MultiFab want(s.c.ba, s.c.dm, 1, 0);
    want.setVal(0.0);
    gust::add_canopy_heat(s.heat, s.q, want);
    amrex::MultiFab diff(s.c.ba, s.c.dm, 1, 0);
    amrex::MultiFab::Copy(diff, R1.theta, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(diff, R0.theta, 0, 0, 1, 0);
    amrex::MultiFab::Subtract(diff, want, 0, 0, 1, 0);
    ASSERT_GT(static_cast<double>(want.norm0()), 0.0);
    EXPECT_LE(static_cast<double>(diff.norm0()), (kDouble ? 1e-12 : 1e-4) * static_cast<double>(want.norm0()));

    // The rows: implicit_rows' coef_P with the heat minus without is lambda.
    amrex::Array<amrex::MultiFab, AMREX_SPACEDIM> F;
    gust::transport_fluxes(s.q.u, s.q.v, s.q.w, s.q.z_phys_nd, s.c.mesh, s.c.geom, F);
    const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM> Fp{&F[0], &F[1], &F[2]};
    gust::ClosureFields cf;
    cl.evaluate(s.q, s.dist, gust::Roughness{}, gust::BuoyancyOptions{}, s.c.mesh, s.c.geom, cf);
    gust::CanopyRows rows;
    gust::canopy_rows(faces, s.q.u, s.q.v, s.q.w, rows);
    amrex::MultiFab A0, A1;
    gust::implicit_rows(Fp, s.q.z_phys_nd, in0, gust::Equation::theta, s.c.mesh, s.c.geom, A0, nullptr, &cf, &rows);
    EXPECT_TRUE(rows.theta.empty());
    gust::canopy_heat_rows(s.heat, s.q.u, s.q.v, s.q.w, rows);
    gust::implicit_rows(Fp, s.q.z_phys_nd, in1, gust::Equation::theta, s.c.mesh, s.c.geom, A1, nullptr, &cf, &rows);
    // Central differences of the leaf term in theta (it is linear in theta).
    const double h = 1e-3;
    gust::State qp, qm;
    for (gust::State* qq : {&qp, &qm}) {
        qq->define(s.c.ba, s.c.dm);
        for (const auto& p : {std::make_pair(&qq->u, &s.q.u), std::make_pair(&qq->v, &s.q.v), std::make_pair(&qq->w, &s.q.w),
                        std::make_pair(&qq->theta, &s.q.theta), std::make_pair(&qq->tke, &s.q.tke),
                        std::make_pair(&qq->p, &s.q.p), std::make_pair(&qq->z_phys_nd, &s.q.z_phys_nd)}) {
            amrex::MultiFab::Copy(*p.first, *p.second, 0, 0, 1, p.second->nGrow());
        }
    }
    qp.theta.plus(static_cast<amrex::Real>(h), 0, 1, qp.theta.nGrow());
    qm.theta.plus(static_cast<amrex::Real>(-h), 0, 1, qm.theta.nGrow());
    amrex::MultiFab Sp(s.c.ba, s.c.dm, 1, 0), Sm(s.c.ba, s.c.dm, 1, 0);
    Sp.setVal(0.0);
    Sm.setVal(0.0);
    gust::add_canopy_heat(s.heat, qp, Sp);
    gust::add_canopy_heat(s.heat, qm, Sm);
    double worst_rows = 0, worst_fd = 0, lmax = 0;
    for (amrex::MFIter mfi(A0, false); mfi.isValid(); ++mfi) {
        const auto a0 = host(A0, mfi), a1 = host(A1, mfi), lam = host(rows.theta, mfi), sp = host(Sp, mfi), sm = host(Sm, mfi);
        const auto a = host(s.canopy.lad(), mfi);
        const auto u = host(s.q.u, mfi), v = host(s.q.v, mfi), w = host(s.q.w, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            const double l = lambda_by_hand(a(iv), cell_wind(u, v, w, i, j, k));
            lmax = std::max(lmax, l);
            worst_rows = std::max({worst_rows, std::abs((a1(iv, gust::coef_P) - a0(iv, gust::coef_P)) - l), std::abs(lam(iv) - l)});
            for (int n = 0; n < gust::n_implicit_coefs; ++n) {
                if (n != gust::coef_P && a1(iv, n) != a0(iv, n)) { worst_rows = std::max(worst_rows, 1.0); }
            }
            worst_fd = std::max(worst_fd, std::abs((sp(iv) - sm(iv)) / (2 * h) - lam(iv)));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst_rows);
    amrex::ParallelDescriptor::ReduceRealMax(worst_fd);
    amrex::ParallelDescriptor::ReduceRealMax(lmax);
    ASSERT_GT(lmax, 0.0);
    EXPECT_LE(worst_rows, (kDouble ? 1e-13 : 1e-5) * lmax);
    EXPECT_LE(worst_fd, (kDouble ? 1e-9 : 1e-2) * lmax);
    // UT-8.4.3: without in.canopy_heat the rows ignore canopy->theta.
    amrex::MultiFab A2;
    gust::implicit_rows(Fp, s.q.z_phys_nd, in0, gust::Equation::theta, s.c.mesh, s.c.geom, A2, nullptr, &cf, &rows);
    EXPECT_EQ(gust::count_bit_differences(A0, A2, amrex::IntVect(0)), 0);
}

// UT-8.4.4: canopy.heat_frozen = blend: theta = theta_f exactly where a = 0,
// theta_f + w (theta_l - theta_f) with w = lambda k / (eps + lambda k) inside,
// 0 <= w < 1; w = 1/2 where lambda k = eps; a second blend from the same
// fields changes nothing (it sets theta from theta_f, not from theta).
TEST(CanopyHeat, FrozenBlend)
{
    EXPECT_EQ(static_cast<double>(gust::frozen_blend_weight(amrex::Real(0.2), amrex::Real(0.5), amrex::Real(0.1))), 0.5);
    EXPECT_EQ(static_cast<double>(gust::frozen_blend_weight(amrex::Real(0), amrex::Real(0.5), amrex::Real(0.1))), 0.0);
    EXPECT_EQ(static_cast<double>(gust::frozen_blend_weight(amrex::Real(0), amrex::Real(0), amrex::Real(0))), 0.0);
    HeatSetup s(8, gust::CanopyFrozenHeat::blend);
    EXPECT_EQ(s.heat.frozen(), gust::CanopyFrozenHeat::blend);
    amrex::MultiFab theta_f(s.c.ba, s.c.dm, 1, 0), eps(s.c.ba, s.c.dm, 1, 0);
    amrex::MultiFab::Copy(theta_f, s.q.theta, 0, 0, 1, 0);
    fill_component(eps, 3, [] (double X, double, double Z) { return 0.01 * (1 + 0.5 * std::cos(kw * X)) / (1 + 0.1 * Z); });
    // k below the floor (0.3) in part of the field: the weight takes the floored k.
    const double kfloor = 0.3;
    const amrex::Real moved = gust::blend_frozen_theta(s.heat, theta_f, eps, amrex::Real(kfloor), s.q);
    EXPECT_GT(static_cast<double>(moved), 0.0);
    double worst = 0, wmax = 0;
    amrex::Long outside_changed = 0, inside = 0;
    for (amrex::MFIter mfi(theta_f, false); mfi.isValid(); ++mfi) {
        const auto tf = host(theta_f, mfi), th = host(s.q.theta, mfi), e = host(eps, mfi), kk = host(s.q.tke, mfi);
        const auto a = host(s.canopy.lad(), mfi);
        const auto u = host(s.q.u, mfi), v = host(s.q.v, mfi), w = host(s.q.w, mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
            const amrex::IntVect iv(i, j, k);
            if (!(a(iv) > 0)) {
                if (th(iv) != tf(iv)) { ++outside_changed; }
                return;
            }
            ++inside;
            const double lk = lambda_by_hand(a(iv), cell_wind(u, v, w, i, j, k)) * std::max(static_cast<double>(kk(iv)), kfloor);
            const double wt = lk / (e(iv) + lk);
            wmax = std::max(wmax, wt);
            worst = std::max(worst, std::abs(th(iv) - (tf(iv) + wt * (kThetaLeaf - tf(iv)))));
        });
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    amrex::ParallelDescriptor::ReduceRealMax(wmax);
    amrex::ParallelDescriptor::ReduceLongSum(outside_changed);
    amrex::ParallelDescriptor::ReduceLongSum(inside);
    EXPECT_EQ(outside_changed, 0L);
    EXPECT_GT(inside, 100L);
    EXPECT_LE(worst, kDouble ? 1e-11 : 1e-3);
    EXPECT_GT(wmax, 0.05);
    EXPECT_LT(wmax, 1.0);
    EXPECT_EQ(static_cast<double>(gust::blend_frozen_theta(s.heat, theta_f, eps, amrex::Real(kfloor), s.q)), 0.0);
}

// The option checks name the input.
TEST(CanopyHeat, OptionChecks)
{
    auto errors = [] (const gust::CanopyOptions& o) {
        std::string all;
        for (const auto& e : gust::check_canopy_options(o)) { all += e + "\n"; }
        return all;
    };
    gust::CanopyOptions o = patch_options(true);
    EXPECT_EQ(errors(o), "");
    for (const double bad : {-1.0, 25.0, 149.9, 400.1, HUGE_VAL}) {
        o.heat.leaf_theta = bad;
        EXPECT_NE(errors(o).find("canopy.leaf_theta must be in [150, 400] K"), std::string::npos) << bad;
    }
    {
        // A NaN raises FE_INVALID in the check's comparisons; the CI builds trap it.
        const gust::test::ScopedNoInvalidTrap no_trap;
        o.heat.leaf_theta = std::nan("");
        EXPECT_NE(errors(o).find("canopy.leaf_theta must be in [150, 400] K"), std::string::npos);
    }
    o.heat.leaf_theta = 150;
    EXPECT_EQ(errors(o), "");
    o.heat.leaf_theta = 400;
    EXPECT_EQ(errors(o), "");
    o = patch_options(true);
    o.given = {"type", "file", "heat"};
    EXPECT_NE(errors(o).find("canopy.heat = fixed needs canopy.leaf_theta"), std::string::npos);
    o = patch_options(true);
    o.heat.frozen_name = "relax";
    EXPECT_NE(errors(o).find("canopy.heat_frozen must be skip or blend (got relax)"), std::string::npos);
    o = patch_options(false);
    o.heat.name = "on";
    EXPECT_NE(errors(o).find("canopy.heat must be off or fixed (got on)"), std::string::npos);
    o = patch_options(false);
    o.given.insert(o.given.end(), {"leaf_theta", "heat_frozen"});
    EXPECT_NE(errors(o).find("canopy.leaf_theta has no effect with canopy.heat = off"), std::string::npos);
    EXPECT_NE(errors(o).find("canopy.heat_frozen has no effect with canopy.heat = off"), std::string::npos);
}
