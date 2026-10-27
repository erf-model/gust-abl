// The canopy's drag field f = Cd a(z) (Canopy.H, Canopy_K.H), plan rev. 6
// beyond the plan, sub-phase 8.1 (UT-8.1.x):
//   the Lalic–Mihailovic profile integrates to one, its N pinned to quadrature
//   (UT-8.1.2); a thin cell gives ERF's point value of the profile;
//   every column's integral of f is Cd LAI, part-filled top cells included,
//   on non-uniform maps, both profiles (UT-8.1.1);
//   on a hill the canopy follows the local ground (UT-8.1.3);
//   patches give the same field on every box layout, and one across a
//   periodic seam reaches the columns on the other side (UT-8.1.4);
//   the coarse levels keep the integral of f over each layer, on a smoothed
//   hierarchy, gathered or not (UT-8.1.5);
//   a restart with another canopy is refused (UT-8.1.7);
//   the inputs and files are checked.
#include "Canopy.H"
#include "Canopy_K.H"
#include "MeshHierarchy.H"
#include "PoissonTestHelpers.H"

#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

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

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-12 : 3e-5;

// A flat or hilly mesh, n x n x 24 cells over 400 m x 400 m, 1 m first
// layer stretched by 1.1 (cells part-filled by a 17.3 m canopy), boxes of
// mgs cells, periodic in x (and y when told).
PCase make_case (int n, int mgs, bool periodic_y, const std::function<double(double, double)>& ground)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(n, n, 24);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {1, periodic_y ? 1 : 0, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.mesh.dz0 = amrex::Real(1);
    c.in.mesh.stretch_ratio = amrex::Real(1.1);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make(ground);
    return c;
}

double flat (double, double) { return 0.0; }

// Write a file on the I/O rank (precision 17); every rank waits.
std::string write_file (const std::string& name, const std::function<void(std::ofstream&)>& body)
{
    const fs::path dir = fs::temp_directory_path() / ("gust_canopy_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / name).string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        os.precision(17);
        body(os);
    }
    amrex::ParallelDescriptor::Barrier();
    return file;
}

// "x y value" at every column centre of an n x n mesh over 400 m (IDW hits
// each point exactly).
std::string write_map (const std::string& name, int n, const std::function<double(double, double)>& v)
{
    return write_file(name, [&] (std::ofstream& os) {
        const double d = 400.0 / n;
        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                const double x = (i + 0.5) * d, y = (j + 0.5) * d;
                os << x << " " << y << " " << v(x, y) << "\n";
            }
        }
    });
}

// The message of the abort f raises ("" if none).
std::string abort_message (const std::function<void()>& f)
{
    try {
        f();
    } catch (const std::runtime_error& e) {
        return e.what();
    }
    return "";
}

// The largest |column integral of f - cd lai(y)| / (cd lai(y)) over the
// columns, with lai(y) = lai0 + lai_y y at the column centre, and the
// number of columns whose integral is not finite.
struct ColumnError { double max_rel{0}; long bad{0}; };
ColumnError column_error (const amrex::MultiFab& f, const gust::Mesh& mesh, const amrex::Geometry& geom, double cd,
                          double lai0, double lai_y)
{
    const amrex::Box& dom = geom.Domain();
    const int k0 = dom.smallEnd(2), k1 = dom.bigEnd(2);
    const double dz = static_cast<double>(geom.CellSize(2)), dy = static_cast<double>(geom.CellSize(1));
    const double ylo = static_cast<double>(geom.ProbLo(1));
    amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum> ops;
    amrex::ReduceData<double, amrex::Long> data(ops);
    for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
        amrex::Box b = mfi.validbox();
        b.setBig(2, b.smallEnd(2));
        const auto a = f.const_array(mfi);
        const auto J = mesh.detJ.const_array(mfi);
        ops.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int) -> amrex::GpuTuple<double, amrex::Long> {
            double s = 0;
            for (int k = k0; k <= k1; ++k) { s += static_cast<double>(a(i, j, k)) * static_cast<double>(J(i, j, k)) * dz; }
            const double want = cd * (lai0 + lai_y * (ylo + (j + 0.5) * dy));
            if (!std::isfinite(s)) { return {0.0, amrex::Long(1)}; }
            return {std::abs(s - want) / want, amrex::Long(0)};
        });
    }
    const auto r = data.value();
    ColumnError e{amrex::get<0>(r), amrex::get<1>(r)};
    amrex::ParallelAllReduce::Max(e.max_rel, amrex::ParallelDescriptor::Communicator());
    amrex::ParallelAllReduce::Sum(e.bad, amrex::ParallelDescriptor::Communicator());
    return e;
}

// Cells with f > 0 whose lower face is at or above the canopy top (h(x) =
// h0 + h_x x), plus cells with f <= 0 whose lower face is below it.
long misplaced_cells (const amrex::MultiFab& f, const amrex::MultiFab& z, const amrex::Geometry& geom, double h0, double h_x)
{
    const double dx = static_cast<double>(geom.CellSize(0)), xlo = static_cast<double>(geom.ProbLo(0));
    const int k0 = geom.Domain().smallEnd(2);
    amrex::ReduceOps<amrex::ReduceOpSum> ops;
    amrex::ReduceData<amrex::Long> data(ops);
    for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
        const auto a = f.const_array(mfi);
        const auto zn = z.const_array(mfi);
        ops.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long> {
            auto face = [&] (int kk) {
                return 0.25 * (static_cast<double>(zn(i, j, kk)) + static_cast<double>(zn(i + 1, j, kk))
                               + static_cast<double>(zn(i, j + 1, kk)) + static_cast<double>(zn(i + 1, j + 1, kk)));
            };
            const double zb = face(k) - face(k0);
            const double h = h0 + h_x * (xlo + (i + 0.5) * dx);
            const bool inside = zb < h;
            return {(inside != (a(i, j, k) > 0)) ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    amrex::Long n = amrex::get<0>(data.value());
    amrex::ParallelAllReduce::Sum(n, amrex::ParallelDescriptor::Communicator());
    return n;
}

// The largest |f(k0) - want| / want over the columns (the first cell lies
// inside the canopy everywhere).
double first_cell_error (const amrex::MultiFab& f, const amrex::Geometry& geom, double want)
{
    const int k0 = geom.Domain().smallEnd(2);
    amrex::ReduceOps<amrex::ReduceOpMax> ops;
    amrex::ReduceData<double> data(ops);
    for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
        amrex::Box b = mfi.validbox();
        b.setBig(2, k0);
        const auto a = f.const_array(mfi);
        ops.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<double> {
            return {std::abs(static_cast<double>(a(i, j, k)) - want) / want};
        });
    }
    double e = amrex::get<0>(data.value());
    amrex::ParallelAllReduce::Max(e, amrex::ParallelDescriptor::Communicator());
    return e;
}

// sum of f detJ over each layer k (index k - k0), all ranks.
std::vector<double> layer_integrals (const amrex::MultiFab& f, const amrex::MultiFab& detJ, const amrex::Geometry& geom)
{
    const amrex::Box& dom = geom.Domain();
    const double area = static_cast<double>(geom.CellSize(0)) * static_cast<double>(geom.CellSize(1));
    std::vector<double> out;
    for (int kk = dom.smallEnd(2); kk <= dom.bigEnd(2); ++kk) {
        amrex::ReduceOps<amrex::ReduceOpSum> ops;
        amrex::ReduceData<double> data(ops);
        for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
            amrex::Box b = mfi.validbox();
            b.setSmall(2, kk);
            b.setBig(2, kk);
            const auto a = f.const_array(mfi);
            const auto J = detJ.const_array(mfi);
            ops.eval(b, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<double> {
                return {static_cast<double>(a(i, j, k)) * static_cast<double>(J(i, j, k)) * area};
            });
        }
        double s = amrex::get<0>(data.value());
        amrex::ParallelAllReduce::Sum(s, amrex::ParallelDescriptor::Communicator());
        out.push_back(s);
    }
    return out;
}

// One cell's value (valid or ghost), on every rank.
double cell_value (const amrex::MultiFab& f, int i, int j, int k)
{
    const amrex::Box b(amrex::IntVect(i, j, k), amrex::IntVect(i, j, k));
    // max over the boxes that hold it (ghosts included); -inf where none does.
    return static_cast<double>(f.max(b, 0, f.nGrow(), false));
}

// f on another layout, copied to a's layout (valid cells): max |a - b|.
double layout_difference (const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    amrex::MultiFab c(a.boxArray(), a.DistributionMap(), 1, 0);
    c.ParallelCopy(b, 0, 0, 1);
    amrex::MultiFab::Subtract(c, a, 0, 0, 1, 0);
    return static_cast<double>(c.norm0(0));
}

gust::CanopyOptions map_options (const std::string& h, const std::string& lai, double cd, bool lm, double laimax = 0.8)
{
    gust::CanopyOptions o;
    o.type = gust::CanopyType::map;
    o.type_name = "map";
    o.height_file = h;
    o.lai_file = lai;
    o.cd = cd;
    o.given = {"type", "height_file", "lai_file", "cd"};
    if (lm) {
        o.profile = gust::CanopyProfile::lalic_mihailovic;
        o.profile_name = "lalic_mihailovic";
        o.laimax = laimax;
        o.given.emplace_back("profile");
        o.given.emplace_back("laimax");
    }
    return o;
}

gust::CanopyOptions patch_options (const std::string& file)
{
    gust::CanopyOptions o;
    o.type = gust::CanopyType::patches;
    o.type_name = "patches";
    o.file = file;
    o.given = {"type", "file"};
    return o;
}

} // namespace

// UT-8.1.2: N in closed form against quadrature (scipy.integrate.quad, 1e-12),
// the shape integrates to N, the cumulative fraction is continuous at
// laimax, rises from 0 to 1 and its slope is P / N.
TEST(Canopy, LalicMihailovicProfileIntegratesToOne)
{
    // N from quad over [0, laimax] and [laimax, 1].
    const std::vector<std::pair<double, double>> pinned{
        {0.0, 0.688640915162}, {0.3, 0.747688277971}, {0.8, 0.313841713572}, {0.95, 0.078941679300}};
    for (const auto& [m, n_quad] : pinned) {
        const gust::LmProfile p = gust::make_lm_profile(m);
        EXPECT_NEAR(p.norm, n_quad, 1e-11) << "laimax " << m;
        // Simpson on [0, m] and [m, 1 - 1e-12] (the shape has a kink at m).
        auto simpson = [&p] (double a, double b, int n) {
            if (b <= a) { return 0.0; }
            const double h = (b - a) / n;
            double s = gust::lm_shape(a, p) + gust::lm_shape(b, p);
            for (int i = 1; i < n; ++i) { s += (i % 2 == 1 ? 4.0 : 2.0) * gust::lm_shape(a + i * h, p); }
            return s * h / 3;
        };
        const double integral = simpson(0, m, 20000) + simpson(m, 1 - 1e-12, 200000);
        EXPECT_NEAR(integral / p.norm, 1.0, 1e-8) << "laimax " << m;
        EXPECT_EQ(gust::lm_cdf(0.0, p), 0.0);
        EXPECT_EQ(gust::lm_cdf(1.0, p), 1.0);
        EXPECT_NEAR(gust::lm_cdf(m * (1 - 1e-15), p), gust::lm_cdf(m * (1 + 1e-15) + 1e-15, p), 1e-12) << "laimax " << m;
        double prev = 0;
        for (double eta = 0.01; eta < 1; eta += 0.01) {
            const double F = gust::lm_cdf(eta, p);
            EXPECT_GE(F, prev);
            const double d = 1e-6;
            const double slope = (gust::lm_cdf(eta + d, p) - gust::lm_cdf(eta - d, p)) / (2 * d);
            if (std::abs(eta - m) > 2 * d) { EXPECT_NEAR(slope, gust::lm_shape(eta, p) / p.norm, 1e-6) << "eta " << eta; }
            prev = F;
        }
        // The fraction between two heights is F(b) - F(a), and stays > 0 just
        // below the top, where F is within rounding of 1.
        for (const double a : {0.0, 0.1, m * 0.5, m, 0.5 * (m + 1), 0.9}) {
            for (const double b : {a + 0.01, a + 0.2, 1.0, 1.5}) {
                EXPECT_NEAR(gust::lm_fraction(a, b, p), gust::lm_cdf(b, p) - gust::lm_cdf(a, p), 1e-13) << a << " " << b;
            }
        }
        EXPECT_GT(gust::lm_fraction(1 - 1e-3, 1, p), 0.0);
        EXPECT_LT(gust::lm_fraction(1 - 1e-3, 1, p), 1e-6);
        EXPECT_EQ(gust::lm_fraction(0.5, 0.5, p), 0.0);
    }
}

// A thin cell gives ERF's point value cd (LAI / h) P(z / h) / N (ERF_ForestDrag.cpp);
// above the top 0; a cell across the top takes the part below it.
TEST(Canopy, ThinCellGivesErfPointValue)
{
    const double h = 22, lai = 5, cd = 0.2, d = 1e-5;
    const gust::LmProfile p = gust::make_lm_profile(0.7);
    for (const double z : {0.5, 5.0, 15.0, 15.4, 21.0}) {
        const double erf_value = cd * lai / h * gust::lm_shape(z / h, p) / p.norm;
        const double ours = gust::cell_drag(z - d, z + d, h, lai, cd, [&p] (double a, double b) { return gust::lm_fraction(a, b, p); });
        EXPECT_NEAR(ours / erf_value, 1.0, 1e-6) << "z " << z;
        const double uni = gust::cell_drag(z - d, z + d, h, lai, cd, [] (double a, double b) { return gust::uniform_fraction(a, b); });
        EXPECT_NEAR(uni / (cd * lai / h), 1.0, 1e-9) << "z " << z;   // eta of a 2e-5 m cell: ~1e-10 rounding
    }
    auto uni = [] (double a, double b) { return gust::uniform_fraction(a, b); };
    EXPECT_EQ(gust::cell_drag(22, 30, h, lai, cd, uni), 0.0);
    EXPECT_NEAR(gust::cell_drag(20, 30, h, lai, cd, uni), cd * lai / h * 2 / 10, 1e-15);
    EXPECT_EQ(gust::cell_drag(0, 1, 0, lai, cd, uni), 0.0);
    EXPECT_EQ(gust::cell_drag(0, 1, h, 0, cd, uni), 0.0);
}

// UT-8.1.1: every column's sum of f detJ dzeta is Cd LAI, with the canopy
// top inside a stretched cell, a height and an LAI that vary over the
// domain, both profiles, Cd constant and from a map; f is 0 exactly in the
// cells above the canopy and > 0 below.
TEST(Canopy, ColumnIntegralIsCdLai)
{
    const int n = 16;
    PCase c = make_case(n, 8, true, flat);
    const std::string hf = write_map("h.xyz", n, [] (double x, double) { return 17.3 + 0.01 * x; });
    const std::string lf = write_map("lai.xyz", n, [] (double, double y) { return 2.0 + 0.01 * y; });
    const std::string cf = write_map("cd.xyz", n, [] (double, double) { return 0.15; });
    for (const bool lm : {false, true}) {
        for (const bool cd_map : {false, true}) {
            gust::CanopyOptions o = map_options(hf, lf, 0.25, lm, 0.65);
            double cd = 0.25;
            if (cd_map) {
                o.cd_file = cf;
                cd = 0.15;
            }
            gust::Canopy can;
            can.build(o, c.geom, c.z, c.ba, c.dm);
            const ColumnError e = column_error(can.drag(), c.mesh, c.geom, cd, 2.0, 0.01);
            EXPECT_EQ(e.bad, 0);
            EXPECT_LT(e.max_rel, rtol) << "profile " << (lm ? "lalic_mihailovic" : "uniform") << ", cd map " << cd_map;
            EXPECT_EQ(misplaced_cells(can.drag(), c.z, c.geom, 17.3, 0.01), 0);
        }
    }
}

// UT-8.1.3: on a hill (50 m, BTF) a 20 m canopy sits on the local ground:
// every column integrates to Cd LAI and the first cell (inside the canopy
// everywhere) is cd LAI / h; heights from sea level would leave the hill bare.
TEST(Canopy, FollowsLocalGroundOnHill)
{
    const int n = 16;
    PCase c = make_case(n, 8, true, [] (double x, double y) {
        return 25 + 25 * std::cos(2 * M_PI * x / 400) * std::sin(2 * M_PI * y / 400);
    });
    const std::string hf = write_map("h20.xyz", n, [] (double, double) { return 20.0; });
    const std::string lf = write_map("lai3.xyz", n, [] (double, double) { return 3.0; });
    gust::Canopy can;
    can.build(map_options(hf, lf, 0.25, false), c.geom, c.z, c.ba, c.dm);
    const ColumnError e = column_error(can.drag(), c.mesh, c.geom, 0.25, 3.0, 0.0);
    EXPECT_EQ(e.bad, 0);
    EXPECT_LT(e.max_rel, rtol);
    EXPECT_LT(first_cell_error(can.drag(), c.geom, 0.25 * 3.0 / 20.0), rtol);
    EXPECT_EQ(misplaced_cells(can.drag(), c.z, c.geom, 20.0, 0.0), 0);
}

// UT-8.1.4: three patches (two overlapping across box edges, one across the
// periodic seam x = 0, one Lalic–Mihailovic) give the same field on one box,
// 16 boxes and uneven boxes; the seam patch reaches the last columns; the
// overlap adds; the periodic ghost column is the image, the non-periodic one
// the nearest column.
TEST(Canopy, PatchesSameOnEveryLayoutAcrossSeam)
{
    const std::string pf = write_file("patches.txt", [] (std::ofstream& os) {
        os << "# type x y height diameter cd lai laimax\n"
           << "1 200 200 18 90 0.2 4 0.8\n"
           << "2 215 190 25 60 0.3 6 0.7   # overlaps the first\n"
           << "\n"
           << "1 5 100 12 60 0.25 3 0\n";
    });
    const int n = 32;
    PCase one = make_case(n, 32, false, flat);
    PCase many = make_case(n, 8, false, flat);
    PCase uneven = make_case(n, 12, false, flat);
    gust::Canopy a, b, u;
    a.build(patch_options(pf), one.geom, one.z, one.ba, one.dm);
    b.build(patch_options(pf), many.geom, many.z, many.ba, many.dm);
    u.build(patch_options(pf), uneven.geom, uneven.z, uneven.ba, uneven.dm);
    EXPECT_EQ(layout_difference(a.drag(), b.drag()), 0.0);
    EXPECT_EQ(layout_difference(a.drag(), u.drag()), 0.0);

    // The seam patch (centre x = 5, radius 30): column i = 31 (x = 393.75)
    // is 11.25 m from it through the seam; row j = 8 (y = 106.25).
    const amrex::Vector<amrex::Real>& zl = one.mesh.zlevels();
    const double z0 = static_cast<double>(zl[0]), z1 = static_cast<double>(zl[1]);
    auto uni = [] (double a, double b) { return gust::uniform_fraction(a, b); };
    const double want_seam = gust::cell_drag(0, z1 - z0, 12, 3, 0.25, uni);
    EXPECT_NEAR(cell_value(a.drag(), 31, 8, 0), want_seam, rtol * want_seam);
    EXPECT_NEAR(cell_value(a.drag(), 0, 8, 0), want_seam, rtol * want_seam);
    // (200, 200) is in both overlapping patches: they add.
    const gust::LmProfile p = gust::make_lm_profile(0.7);
    const double want_both = gust::cell_drag(0, z1 - z0, 18, 4, 0.2, uni)
                             + gust::cell_drag(0, z1 - z0, 25, 6, 0.3, [&p] (double a, double b) { return gust::lm_fraction(a, b, p); });
    EXPECT_NEAR(cell_value(a.drag(), 15, 15, 0), want_both, rtol * want_both);
    // Layer 11 (18.53 to 21.38 m) is above the 18 m patch: the 25 m one alone.
    const double z11 = static_cast<double>(zl[11] - zl[0]), z12 = static_cast<double>(zl[12] - zl[0]);
    ASSERT_GT(z11, 18.0);
    const double want_top = gust::cell_drag(z11, z12, 25, 6, 0.3, [&p] (double a, double b) { return gust::lm_fraction(a, b, p); });
    EXPECT_NEAR(cell_value(a.drag(), 15, 15, 11), want_top, rtol * want_top);
    // Ghosts: periodic in x (i = -1 is i = 31), nearest column in y (j = -1 is j = 0).
    EXPECT_EQ(cell_value(b.drag(), -1, 8, 0), cell_value(b.drag(), 31, 8, 0));
    EXPECT_EQ(cell_value(b.drag(), 15, -1, 0), cell_value(b.drag(), 15, 0, 0));
    EXPECT_EQ(cell_value(b.drag(), 15, n, 0), cell_value(b.drag(), 15, n - 1, 0));
}

// UT-8.1.5: the coarse levels of a smoothed hierarchy (round hill, levels
// 1 and 2, gathered on one box or not) keep each layer's integral of f.
TEST(Canopy, RestrictionKeepsLayerIntegrals)
{
    for (const int gather : {0, 1}) {
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
        gust::MeshHierarchy h(c.in, c.geom, c.z, c.mesh, gather);
        const std::string hf = write_file("h_hier.xyz", [] (std::ofstream& os) {
            for (int j = 0; j < 9; ++j) {
                for (int i = 0; i < 9; ++i) { os << 50 * i << " " << 50 * j << " " << 10 + 2 * i + j << "\n"; }
            }
        });
        const std::string lf = write_file("lai_hier.xyz", [] (std::ofstream& os) {
            for (int j = 0; j < 9; ++j) {
                for (int i = 0; i < 9; ++i) { os << 50 * i << " " << 50 * j << " " << ((i + j) % 3 == 0 ? 0 : 2 + i) << "\n"; }
            }
        });
        gust::Canopy fine;
        fine.build(map_options(hf, lf, 0.2, true), c.geom, c.z, c.ba, c.dm);
        const std::vector<double> want = layer_integrals(fine.drag(), c.mesh.detJ, c.geom);
        double total = 0;
        for (const double v : want) { total += v; }
        ASSERT_GT(total, 0.0);
        for (int l = 1; l < 3; ++l) {
            gust::Canopy coarse;
            coarse.restrict_from(fine, c.mesh.detJ, h.geom(l), h.mesh(l), h.boxArray(l), h.DistributionMap(l));
            EXPECT_EQ(coarse.input_hash(), fine.input_hash());
            const std::vector<double> got = layer_integrals(coarse.drag(), h.mesh(l).detJ, h.geom(l));
            ASSERT_EQ(got.size(), want.size());
            for (std::size_t k = 0; k < want.size(); ++k) {
                EXPECT_NEAR(got[k], want[k], rtol * total) << "level " << l << " layer " << k << " gather " << gather;
            }
        }
    }
}

// UT-8.1.7: the checkpoint's canopy hash: the same inputs restart; a changed
// file, a canopy where the checkpoint has none, and none where it has one,
// are refused.
TEST(Canopy, RestartRefusesAnotherCanopy)
{
    PCase c = make_case(8, 8, true, flat);
    const std::string pa = write_file("restart_a.txt", [] (std::ofstream& os) { os << "1 200 200 18 90 0.2 4 0.8\n"; });
    const std::string pb = write_file("restart_b.txt", [] (std::ofstream& os) { os << "1 200 200 18 90 0.2 4.5 0.8\n"; });
    const fs::path base = fs::temp_directory_path() / ("gust_canopy_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string with = (base / "chk_with").string(), without = (base / "chk_without").string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(with);
        fs::create_directories(without);
        fs::remove(gust::canopy_checkpoint_file(without));
    }
    amrex::ParallelDescriptor::Barrier();
    gust::Canopy a, b, none;
    a.build(patch_options(pa), c.geom, c.z, c.ba, c.dm);
    b.build(patch_options(pb), c.geom, c.z, c.ba, c.dm);
    a.write_checkpoint(with);
    none.write_checkpoint(without);
    auto check = [] (const std::string& dir, const gust::Canopy& c) { gust::Canopy::check_restart(dir, c.input_hash()); };
    EXPECT_EQ(abort_message([&] { check(with, a); }), "");
    EXPECT_EQ(abort_message([&] { check(without, none); }), "");
    EXPECT_NE(abort_message([&] { check(with, b); }).find("the canopy inputs differ"), std::string::npos);
    EXPECT_NE(abort_message([&] { check(with, none); }).find("written with a canopy"), std::string::npos);
    EXPECT_NE(abort_message([&] { check(without, a); }).find("written without a canopy"), std::string::npos);
    EXPECT_NE(a.input_hash(), b.input_hash());
    // The hash follows the contents, not the name.
    const std::string pa2 = write_file("restart_a_copy.txt", [] (std::ofstream& os) { os << "1 200 200 18 90 0.2 4 0.8\n"; });
    EXPECT_EQ(gust::Canopy::hash_inputs(patch_options(pa2)), a.input_hash());
}

// The option checks name the input.
TEST(Canopy, OptionChecks)
{
    const std::string ok = write_file("ok.xyz", [] (std::ofstream& os) { os << "0 0 1\n"; });
    auto first = [] (const gust::CanopyOptions& o) {
        const auto e = gust::check_canopy_options(o);
        return e.empty() ? std::string() : e.front();
    };
    gust::CanopyOptions none;
    EXPECT_TRUE(gust::check_canopy_options(none).empty());
    EXPECT_FALSE(none.active());
    none.given = {"cd"};
    EXPECT_EQ(first(none), "canopy.cd has no effect with canopy.type = none");
    gust::CanopyOptions bad;
    bad.type_name = "forest";
    EXPECT_EQ(first(bad), "canopy.type must be none, patches or map (got forest)");

    gust::CanopyOptions p = patch_options("");
    EXPECT_EQ(first(p), "canopy.file is required with canopy.type = patches");
    p = patch_options("no_such_canopy_file.txt");
    EXPECT_EQ(first(p), "canopy.file = no_such_canopy_file.txt: cannot be read");
    p = patch_options(ok);
    p.given.emplace_back("profile");
    EXPECT_NE(first(p).find("canopy.profile has no effect with canopy.type = patches"), std::string::npos);

    gust::CanopyOptions m = map_options(ok, ok, 0.2, false);
    EXPECT_TRUE(gust::check_canopy_options(m).empty());
    m.given = {"type", "height_file", "lai_file"};
    EXPECT_NE(first(m).find("needs canopy.cd or canopy.cd_file"), std::string::npos);
    m.given = {"type", "height_file", "lai_file", "cd", "cd_file"};
    EXPECT_NE(first(m).find("cannot both be given"), std::string::npos);
    m = map_options(ok, ok, 0.0, false);
    EXPECT_EQ(first(m), "canopy.cd must be finite and > 0");
    m = map_options(ok, ok, std::nan(""), false);
    EXPECT_EQ(first(m), "canopy.cd must be finite and > 0");
    m = map_options(ok, ok, 0.2, true, 1.0);
    EXPECT_NE(first(m).find("canopy.laimax must satisfy 0 <= laimax < 1"), std::string::npos);
    m = map_options(ok, ok, 0.2, true, -0.1);
    EXPECT_NE(first(m).find("canopy.laimax must satisfy 0 <= laimax < 1"), std::string::npos);
    m = map_options(ok, ok, 0.2, false);
    m.given.emplace_back("laimax");
    EXPECT_EQ(first(m), "canopy.laimax has no effect with canopy.profile = uniform");
    m.given.pop_back();
    m.profile_name = "cone";
    EXPECT_EQ(first(m), "canopy.profile must be uniform or lalic_mihailovic (got cone)");
    m = map_options(ok, "", 0.2, false);
    EXPECT_EQ(first(m), "canopy.lai_file is required with canopy.type = map");
    m = map_options(ok, ok, 0.2, false);
    m.file = ok;
    m.given.emplace_back("file");
    EXPECT_NE(first(m).find("canopy.file has no effect with canopy.type = map"), std::string::npos);

    // ParmParse: a repeated key, a value that is not a number, two values.
    {
        amrex::ParmParse pp("canopy");
        pp.add("type", std::string("map"));
        pp.add("cd", std::string("0.2x"));
        pp.addarr("laimax", std::vector<amrex::Real>{amrex::Real(0.5), amrex::Real(0.6)});
        const gust::CanopyOptions r = gust::read_canopy_options(6, 2.0);
        EXPECT_EQ(r.type, gust::CanopyType::map);
        const auto e = gust::check_canopy_options(r);
        auto found = [&e] (const std::string& s) {
            for (const auto& x : e) { if (x.find(s) != std::string::npos) { return true; } }
            return false;
        };
        EXPECT_TRUE(found("canopy.laimax takes exactly one value"));
        EXPECT_TRUE(found("canopy.cd = 0.2x is not a number"));
        pp.add("type", std::string("map"));
        EXPECT_TRUE(gust::read_canopy_options(6, 2.0).parse_errors.front().find("canopy.type is given more than once")
                    != std::string::npos);
        for (const char* k : {"type", "cd", "laimax"}) { while (pp.remove(k) > 0) {} }
    }
}

// The patch file's lines are checked, naming the file and the line.
TEST(Canopy, PatchFileChecks)
{
    const std::string f = write_file("bad_patches.txt", [] (std::ofstream& os) {
        os << "1 0 0 10 20 0.2 3 0.8\n"     // line 1: fine
           << "1 0 0 10 20 0.2 3\n"         // 2: seven values
           << "3 0 0 10 20 0.2 3 0.8\n"     // 3: type 3
           << "1 0 0 0 20 0.2 3 0.8\n"      // 4: height 0
           << "1 0 0 10 -1 0.2 3 0.8\n"     // 5: diameter < 0
           << "1 0 0 10 20 -0.2 3 0.8\n"    // 6: cd < 0
           << "1 0 0 10 20 0.2 -3 0.8\n"    // 7: lai < 0
           << "2 0 0 10 20 0.2 3 1.0\n"     // 8: laimax 1 on a type 2 line
           << "1 0 0 10 20 0.2 3 nan\n"     // 9: not finite
           << "1 0 0 10 20 0.2 3 0.8x\n";   // 10: not a number
    });
    const std::string msg = abort_message([&] { (void)gust::read_canopy_patches(f, "canopy.file"); });
    for (const char* s : {"line 2: 7 values", "line 3: type must be 1", "line 4: height must be > 0",
                          "line 5: diameter must be > 0", "line 6: cd must be >= 0", "line 7: lai must be >= 0",
                          "line 8: laimax must satisfy", "line 9: 'nan' is not a finite number",
                          "line 10: '0.8x' is not a finite number"}) {
        EXPECT_NE(msg.find(s), std::string::npos) << s << "\n" << msg;
    }
    EXPECT_EQ(msg.find("line 1:"), std::string::npos) << msg;
    const std::string empty = write_file("empty_patches.txt", [] (std::ofstream& os) { os << "# nothing\n\n"; });
    EXPECT_NE(abort_message([&] { (void)gust::read_canopy_patches(empty, "canopy.file"); }).find("holds no patch"),
              std::string::npos);
    // A map with a negative value is refused when the field is built.
    PCase c = make_case(8, 8, true, flat);
    const std::string neg = write_file("neg.xyz", [] (std::ofstream& os) { os << "0 0 10\n100 0 -1\n0 100 10\n"; });
    const std::string ok = write_file("pos.xyz", [] (std::ofstream& os) { os << "0 0 1\n100 0 1\n0 100 1\n"; });
    gust::Canopy can;
    gust::CanopyOptions o = map_options(neg, ok, 0.2, false);
    o.idw_neighbors = 3;
    EXPECT_NE(abort_message([&] { can.build(o, c.geom, c.z, c.ba, c.dm); }).find("canopy.height_file = " + neg),
              std::string::npos);
    // A patch wider than the periodic domain (400 m) would reach columns
    // through two images: refused.
    const std::string wide = write_file("wide_patch.txt", [] (std::ofstream& os) { os << "1 200 200 10 450 0.2 3 0\n"; });
    EXPECT_NE(abort_message([&] { can.build(patch_options(wide), c.geom, c.z, c.ba, c.dm); }).find("wider than the periodic domain in x"),
              std::string::npos);
    o = map_options(ok, ok, 0.2, false);
    o.idw_neighbors = 4;
    EXPECT_NE(abort_message([&] { can.build(o, c.geom, c.z, c.ba, c.dm); }).find("terrain.idw_neighbors = 4 must be between 1"),
              std::string::npos);
}
