// The immersed ground as a continuous surface (Buildings::build_from_surface, ImmersedSurface_K.H):
//   UT-14.2.1  the solid fraction, the face fractions, the fluid's mean height above the ground and the
//              wall cells of a tilted plane at 10, 30 and 45 degrees against the closed form (the cell
//              clipped by the plane), with the slope along x and along y;
//   UT-14.2.2  the height above the ground on an Agnesi hill is z - h at every cell, and the parts of
//              the solver that need it read the helper (a search of the sources);
//   beyond the plan: a flat surface gives the fields of PR-87's flat ground.
#include "BuildingForcing.H"
#include "Buildings_K.H"
#include "Closure.H"
#include "ImmersedSurface_K.H"
#include "ImmersedTerrain.H"
#include "LaminarTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using gust::ptest::PCase;

constexpr bool kDouble = std::is_same_v<amrex::Real, double>;
constexpr double tol = kDouble ? 1e-11 : 2e-4;
constexpr double kZ0 = 0.1;
constexpr double kPi = 3.14159265358979323846;

// The node heights with the ground h(x, y) as their second component, at every node (ghosts too).
amrex::MultiFab with_ground (const PCase& c, const std::function<double(double, double)>& h)
{
    amrex::MultiFab z(c.z.boxArray(), c.dm, 2, c.z.nGrow());
    amrex::MultiFab::Copy(z, c.z, 0, 0, 1, c.z.nGrow());
    const auto dx = c.geom.CellSizeArray();
    const auto plo = c.geom.ProbLoArray();
    for (amrex::MFIter mfi(z, false); mfi.isValid(); ++mfi) {
        const auto a = z.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
            a(i, j, k, 1) = static_cast<amrex::Real>(h(static_cast<double>(plo[0]) + i * static_cast<double>(dx[0]),
                                                       static_cast<double>(plo[1]) + j * static_cast<double>(dx[1])));
        });
    }
    return z;
}

gust::TerrainInputs immersed ()
{
    gust::TerrainInputs t;
    t.representation = "immersed_forcing";
    t.flat_height = 25;
    return t;
}

// A file of points on the 10 m nodes of the 160 m x 80 m cases: the ramp slope * i + j (m).
std::string ramp_file (const std::string& name, double slope)
{
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / ("gust_immersed_surface_" + std::to_string(amrex::ParallelDescriptor::NProcs()));
    const std::string file = (dir / name).string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        for (int j = 0; j <= 8; ++j) {
            for (int i = 0; i <= 16; ++i) { os << 10.0 * i << " " << 10.0 * j << " " << slope * i + j << "\n"; }
        }
    }
    amrex::ParallelDescriptor::Barrier();
    return file;
}

// The surface's geometry and its wall law, as App.cpp builds them.
struct Surface
{
    gust::BuildingsOptions opt;
    gust::Buildings b;
    gust::BuildingForcing law;
    amrex::MultiFab z;
    Surface (const PCase& c, const std::function<double(double, double)>& h)
    {
        opt = gust::immersed_terrain_options(immersed(), kZ0);
        z = with_ground(c, h);
        b.build_from_surface(opt, "test", c.geom, z, c.ba, c.dm);
        law.build(b, opt, c.geom, c.mesh.detJ);
    }
};

// The largest value of f over the valid points of mf, on every rank (an EXPECT on one rank does
// not fail the test).
double max_over (const amrex::MultiFab& mf, const std::function<double(const amrex::Array4<amrex::Real const>&, int, int, int)>& f)
{
    double m = 0;
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const auto a = mf.const_array(mfi);
        amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { m = std::max(m, f(a, i, j, k)); });
    }
    amrex::ParallelDescriptor::ReduceRealMax(m);
    return m;
}

// The fluid of the rectangle [s0, s1] x [zb, zt] above the line z = h0 + m s: its area and centroid,
// from the rectangle clipped by the line (Sutherland-Hodgman, one edge) and the shoelace formulas.
struct Clip
{
    double area{0}, sc{0}, zc{0};
};

Clip fluid_above (double s0, double s1, double zb, double zt, double h0, double m)
{
    const std::array<std::array<double, 2>, 4> r{{{s0, zb}, {s1, zb}, {s1, zt}, {s0, zt}}};
    auto side = [&] (const std::array<double, 2>& p) { return p[1] - (h0 + m * p[0]); };
    std::vector<std::array<double, 2>> poly;
    for (std::size_t n = 0; n < 4; ++n) {
        const auto& p = r[n];
        const auto& q = r[(n + 1) % 4];
        const double fp = side(p), fq = side(q);
        if (fp >= 0) { poly.push_back(p); }
        if ((fp > 0 && fq < 0) || (fp < 0 && fq > 0)) {
            const double t = fp / (fp - fq);
            poly.push_back({p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])});
        }
    }
    Clip c;
    if (poly.size() < 3) { return c; }
    double a2 = 0, cs = 0, cz = 0;
    for (std::size_t n = 0; n < poly.size(); ++n) {
        const auto& p = poly[n];
        const auto& q = poly[(n + 1) % poly.size()];
        const double cr = p[0] * q[1] - q[0] * p[1];
        a2 += cr;
        cs += (p[0] + q[0]) * cr;
        cz += (p[1] + q[1]) * cr;
    }
    c.area = 0.5 * a2;
    if (c.area > 0) {
        c.sc = cs / (3 * a2);
        c.zc = cz / (3 * a2);
    }
    return c;
}

double clamp01 (double x) { return std::min(1.0, std::max(0.0, x)); }

std::string source (const std::string& name)
{
    std::ifstream f(std::string(GUST_SOURCE_DIR) + "/Source/" + name);
    std::stringstream s;
    s << f.rdbuf();
    return s.str();
}

// The closed form of a plane h = h0 + m s on cells of ds x dz, s along x (dir 0) or y (dir 1), m > 0.
struct Plane
{
    double h0, m, ds, dz, small;
    [[nodiscard]] double h (double s) const { return h0 + m * s; }
    [[nodiscard]] Clip fluid (int n, int k) const { return fluid_above(n * ds, (n + 1) * ds, k * dz, (k + 1) * dz, h0, m); }
    [[nodiscard]] double beta (int n, int k) const { return gust::clip_fraction(clamp01(1 - fluid(n, k).area / (ds * dz)), small); }
    // The fluid's mean height above the ground: z - h is linear, so its mean is its value at the centroid.
    [[nodiscard]] double height (int n, int k) const
    {
        const Clip c = fluid(n, k);
        return beta(n, k) < 1 ? c.zc - h(c.sc) : 0.5 * dz;
    }
    // A face between cells of fractions lo and hi whose own blocked part is raw (Buildings.H).
    [[nodiscard]] double face (double lo, double hi, double raw) const
    {
        if (lo == 1 || hi == 1) { return 1; }
        if (lo == 0 && hi == 0) { return 0; }
        return gust::clip_fraction(clamp01(raw), small);
    }
    // The face across the slope at s = n ds: the ground's height on it.
    [[nodiscard]] double along (int n, int k) const { return face(beta(n - 1, k), beta(n, k), (h(n * ds) - k * dz) / dz); }
    // The face beside the slope: the cell's own section.
    [[nodiscard]] double beside (int n, int k) const { return face(beta(n, k), beta(n, k), 1 - fluid(n, k).area / (ds * dz)); }
    // The low z face of cell (n, k): the part of the footprint whose ground lies above it.
    [[nodiscard]] double below (int n, int k) const
    {
        return face(k > 0 ? beta(n, k - 1) : 0.0, beta(n, k), (h((n + 1) * ds) - k * dz) / (m * ds));
    }
};

}   // namespace

// The closed-form moments of a triangle and a segment against a quadrature.
TEST (ImmersedSurface, MomentsMatchAQuadrature)
{
    const std::array<double, 3> v{3.0, 7.5, 5.0};
    for (const double t : {2.0, 4.0, 5.0, 6.0, 7.5, 9.0}) {
        // Midpoints of n^2 sub-triangles in barycentric coordinates.
        const int n = 600;
        double m0 = 0, m1 = 0, m2 = 0;
        long count = 0;
        for (int a = 0; a < n; ++a) {
            for (int b = 0; b < 2 * (n - a) - 1; ++b) {
                // Up and down sub-triangles of row a: centroids.
                const double l1 = (a + (b % 2 == 0 ? 1.0 : 2.0) / 3) / n;
                const double l2 = (b / 2 + (b % 2 == 0 ? 1.0 : 2.0) / 3) / n;
                const double h = v[0] * l1 + v[1] * l2 + v[2] * (1 - l1 - l2);
                const double d = std::max(0.0, t - h);
                m0 += d > 0 ? 1 : 0;
                m1 += d;
                m2 += d * d;
                ++count;
            }
        }
        const gust::SurfaceMoments s = gust::triangle_moments(v[0], v[1], v[2], t);
        EXPECT_NEAR(s.m0, m0 / count, 5e-3) << t;
        EXPECT_NEAR(s.m1, m1 / count, 1e-4) << t;
        EXPECT_NEAR(s.m2, m2 / count, 1e-4) << t;
        double seg = 0;
        const int ns = 200000;
        for (int a = 0; a < ns; ++a) { seg += std::max(0.0, t - (v[0] + (v[1] - v[0]) * (a + 0.5) / ns)); }
        EXPECT_NEAR(gust::segment_depth(v[0], v[1], t), seg / ns, 1e-8) << t;
    }
    // Two equal heights, low and high.
    for (const auto& q : {std::array<double, 3>{3.0, 3.0, 6.0}, std::array<double, 3>{3.0, 6.0, 6.0}}) {
        for (const double t : {3.0, 4.5, 6.0, 7.0}) {
            const int n = 400;
            double m0 = 0, m1 = 0, m2 = 0;
            long count = 0;
            for (int a = 0; a < n; ++a) {
                for (int b = 0; b < 2 * (n - a) - 1; ++b) {
                    const double l1 = (a + (b % 2 == 0 ? 1.0 : 2.0) / 3) / n;
                    const double l2 = (b / 2 + (b % 2 == 0 ? 1.0 : 2.0) / 3) / n;
                    const double d = std::max(0.0, t - (q[0] * l1 + q[1] * l2 + q[2] * (1 - l1 - l2)));
                    m0 += d > 0 ? 1 : 0;
                    m1 += d;
                    m2 += d * d;
                    ++count;
                }
            }
            const gust::SurfaceMoments s = gust::triangle_moments(q[0], q[1], q[2], t);
            EXPECT_NEAR(s.m0, m0 / count, 5e-3) << t;
            EXPECT_NEAR(s.m1, m1 / count, 1e-4) << t;
            EXPECT_NEAR(s.m2, m2 / count, 1e-4) << t;
        }
    }
    // A twisted footprint (corner heights 0, 1, 5, 2 around it): the surface is the four triangles
    // about the centre, whose height is the mean of the four corners (2), sampled on a 400 x 400 grid.
    {
        const double h00 = 0, h10 = 1, h11 = 5, h01 = 2, hc = 2;
        const int n = 400;
        for (const double t : {0.5, 1.5, 2.0, 3.0, 6.0}) {
            double m0 = 0, m1 = 0, m2 = 0;
            for (int a = 0; a < n; ++a) {
                for (int b = 0; b < n; ++b) {
                    const double x = (a + 0.5) / n - 0.5, y = (b + 0.5) / n - 0.5;   // from the centre
                    // The triangle the point lies in: the centre and the side nearest to it.
                    double h;
                    if (std::abs(y) >= std::abs(x)) {
                        const double r = 2 * std::abs(y), sx = 0.5 + x / r;   // along the side, 0 to 1
                        const double side = y < 0 ? h00 + (h10 - h00) * sx : h01 + (h11 - h01) * sx;
                        h = hc + r * (side - hc);
                    } else {
                        const double r = 2 * std::abs(x), sy = 0.5 + y / r;
                        const double side = x < 0 ? h00 + (h01 - h00) * sy : h10 + (h11 - h10) * sy;
                        h = hc + r * (side - hc);
                    }
                    const double d = std::max(0.0, t - h);
                    m0 += d > 0 ? 1 : 0;
                    m1 += d;
                    m2 += d * d;
                }
            }
            const gust::SurfaceMoments s = gust::footprint_moments(h00, h10, h11, h01, t);
            EXPECT_NEAR(s.m0, m0 / (n * n), 1e-2) << t;
            EXPECT_NEAR(s.m1, m1 / (n * n), 1e-3) << t;
            EXPECT_NEAR(s.m2, m2 / (n * n), 1e-3) << t;
        }
    }
    // A flat triangle or segment on the level: no fluid below, and no invalid operation.
    EXPECT_EQ(gust::triangle_moments(4, 4, 4, 4).m1, 0.0);
    EXPECT_EQ(gust::segment_depth(4, 4, 4), 0.0);
    EXPECT_EQ(gust::fluid_height({}, {}), 0.0);
}

// UT-14.2.1.
TEST (ImmersedSurface, TiltedPlaneMatchesTheClosedForm)
{
    for (const int mgs : {8, 4}) {
        // 8 x 8 x 24 cells of 10 m x 10 m x 10 m; no periodic seam cuts the plane.
        const PCase c = gust::lamtest::flat_case(8, 8, 24, 80, 80, 240, false, false, mgs);
        for (const double deg : {10.0, 30.0, 45.0}) {
            for (const int dir : {0, 1}) {
                SCOPED_TRACE("boxes of " + std::to_string(mgs) + ", " + std::to_string(deg) + " degrees along " + (dir == 0 ? "x" : "y"));
                const double m = std::tan(deg * kPi / 180);
                Surface g(c, [=] (double x, double y) { return 23.0 + m * (dir == 0 ? x : y); });
                const Plane p{23.0, m, 10.0, 10.0,
                              std::max(g.opt.small_volfrac, 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()))};
                auto along = [dir] (int i, int j) { return dir == 0 ? i : j; };

                EXPECT_LT(max_over(g.b.beta(), [&] (const auto& a, int i, int j, int k) {
                    return std::abs(static_cast<double>(a(i, j, k)) - p.beta(along(i, j), k));
                }), tol);
                EXPECT_LT(max_over(g.b.surface(), [&] (const auto& a, int i, int j, int k) {
                    const double e = std::abs(static_cast<double>(a(i, j, k, gust::kSurfaceHeight)) - p.height(along(i, j), k));
                    const double ex = std::abs(static_cast<double>(a(i, j, k, gust::kSurfaceSlopeX)) - (dir == 0 ? m : 0.0));
                    const double ey = std::abs(static_cast<double>(a(i, j, k, gust::kSurfaceSlopeY)) - (dir == 1 ? m : 0.0));
                    return std::max(e, std::max(ex, ey));
                }), 100 * tol);
                // The faces inside the domain (a face on the domain's side has one cell). A z face is
                // cut at its level less the round-off of a height (64 ulp of the top of the mesh, 240 m:
                // 1.8 mm in single precision), which moves the cut by that over the rise of a cell.
                const double zface_shift = 2 * 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()) * 240 / (m * 10);
                for (int d = 0; d < 3; ++d) {
                    EXPECT_LT(max_over(g.b.beta_face(d), [&] (const auto& a, int i, int j, int k) {
                        const int n = along(i, j);
                        if (k > 23 || i > 7 + (d == 0 ? 1 : 0) || j > 7 + (d == 1 ? 1 : 0)) { return 0.0; }
                        double want = 0;
                        if (d == 2) {
                            want = p.below(n, k);
                        } else if (d == dir) {
                            if (n == 0 || n == 8) { return 0.0; }
                            want = p.along(n, k);
                        } else {
                            if ((d == 0 ? i : j) == 0 || (d == 0 ? i : j) == 8) { return 0.0; }
                            want = p.beside(n, k);
                        }
                        return std::abs(static_cast<double>(a(i, j, k)) - want);
                    }), tol + (d == 2 ? zface_shift : 0.0)) << "faces " << d;
                }
                // The wall cells: the roof bit only (no side walls on a slope), the wall law's two
                // distances, and the roofs' areas of every column add to its footprint.
                EXPECT_LT(max_over(g.b.wall(), [&] (const auto& a, int i, int j, int k) {
                    const int n = along(i, j);
                    const double own = p.beta(n, k);
                    // The closed form, written out rather than through gust::ground_roof: a roof wherever the
                    // low face is more blocked than the high one.
                    const bool wall = own < 1 && p.below(n, k) > (k < 23 ? p.below(n, k + 1) : 0.0);
                    const double mask = std::abs(static_cast<double>(a(i, j, k, gust::kWallMask)) - (wall ? gust::kWallRoof : 0));
                    if (!wall) { return mask; }
                    const double d1 = std::abs(static_cast<double>(a(i, j, k, gust::kWallD1z)) - p.height(n, k));
                    const double d2 = std::abs(static_cast<double>(a(i, j, k, gust::kWallD2z)) - (p.beta(n, k + 1) < 1 ? p.height(n, k + 1) : 0.0));
                    return std::max(mask, std::max(d1, d2) / 100);
                }), tol);
                double worst = 0, cut = 0;
                for (amrex::MFIter mfi(g.b.beta(), false); mfi.isValid(); ++mfi) {
                    const auto be = g.b.beta().const_array(mfi);
                    const auto fz = g.b.beta_face(2).const_array(mfi);
                    const auto wl = g.b.wall().const_array(mfi);
                    const amrex::Box bx = mfi.validbox();
                    for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                        for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                            double sum = 0;
                            for (int k = 0; k < 24; ++k) {
                                if (wl(i, j, k, gust::kWallMask) > 0) {
                                    sum += static_cast<double>(gust::wall_area_fraction(fz(i, j, k), fz(i, j, k + 1), be(i, j, k)));
                                }
                                if (be(i, j, k) > 0 && be(i, j, k) < 1) { cut += 1; }
                            }
                            worst = std::max(worst, std::abs(sum - 1));
                        }
                    }
                }
                amrex::ParallelDescriptor::ReduceRealMax(worst);
                amrex::ParallelDescriptor::ReduceRealSum(cut);
                EXPECT_LT(worst, tol);
                EXPECT_GE(cut, 64.0);   // a cut cell in every column at least
            }
        }
    }
}

// Beyond the plan: a slope that rises through many layers across a column keeps its whole footprint
// as roof with any terrain.small_volfrac. 0.15 over 80 m columns of 1 m layers: 12 cut cells, about
// 1/12 of the footprint each. The buildings' test (wall_on_low_side) with a tolerance of 0.1 keeps a
// third of the footprint as roof (the test fails with it at every terrain.small_volfrac here).
TEST (ImmersedSurface, RoofsAddToTheFootprintOnASlopeThroughManyLayers)
{
    for (const int mgs : {8, 4}) {
        const PCase c = gust::lamtest::flat_case(8, 8, 128, 640, 640, 128, false, false, mgs);
        for (const double small : {0.005, 0.1, 0.3}) {
            SCOPED_TRACE("boxes of " + std::to_string(mgs) + ", terrain.small_volfrac " + std::to_string(small));
            gust::TerrainInputs t = immersed();
            t.small_volfrac = static_cast<amrex::Real>(small);
            const gust::BuildingsOptions opt = gust::immersed_terrain_options(t, kZ0);
            const amrex::MultiFab z = with_ground(c, [] (double x, double) { return 4.0 + 0.15 * x; });
            gust::Buildings b;
            b.build_from_surface(opt, "test", c.geom, z, c.ba, c.dm);
            double worst = 0, steep = 0;
            for (amrex::MFIter mfi(b.beta(), false); mfi.isValid(); ++mfi) {
                const auto be = b.beta().const_array(mfi);
                const auto fz = b.beta_face(2).const_array(mfi);
                const auto wl = b.wall().const_array(mfi);
                const amrex::Box bx = mfi.validbox();
                for (int j = bx.smallEnd(1); j <= bx.bigEnd(1); ++j) {
                    for (int i = bx.smallEnd(0); i <= bx.bigEnd(0); ++i) {
                        double sum = 0;
                        int cut = 0;
                        for (int k = 0; k < 128; ++k) {
                            if (wl(i, j, k, gust::kWallMask) > 0) {
                                sum += static_cast<double>(gust::wall_area_fraction(fz(i, j, k), fz(i, j, k + 1), be(i, j, k)));
                            }
                            if (be(i, j, k) > 0 && be(i, j, k) < 1) { ++cut; }
                        }
                        worst = std::max(worst, std::abs(sum - 1));
                        if (cut >= 4) { steep += 1; }
                    }
                }
            }
            amrex::ParallelDescriptor::ReduceRealMax(worst);
            amrex::ParallelDescriptor::ReduceRealSum(steep);
            EXPECT_LT(worst, tol);
            EXPECT_EQ(steep, 64.0) << "a column with fewer than 4 cut cells: the case does not spread the roof";
        }
    }
}

// Beyond the plan: a flat surface gives the cells, faces and wall cells of PR-87's flat ground
// (Buildings::build_from_map), with the ground inside a layer and on a cell face.
TEST (ImmersedSurface, FlatSurfaceGivesTheFlatGroundsFields)
{
    const PCase c = gust::lamtest::flat_case(4, 4, 24, 40, 40, 240, true, true, 2);
    for (const double height : {23.0, 20.0, 29.5}) {
        SCOPED_TRACE(height);
        Surface g(c, [=] (double, double) { return height; });
        gust::TerrainInputs t = immersed();
        t.flat_height = static_cast<amrex::Real>(height);
        gust::Buildings flat;
        flat.build_from_map(g.opt, gust::immersed_terrain_map(t, c.geom), "terrain.flat_height", c.geom, c.z, c.ba, c.dm);
        auto differ = [] (const amrex::MultiFab& a, const amrex::MultiFab& b, int comp, int ng) {
            amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, ng);
            amrex::MultiFab::Copy(d, a, comp, 0, 1, ng);
            amrex::MultiFab::Subtract(d, b, comp, 0, 1, ng);
            return static_cast<double>(d.norminf(0, ng));
        };
        EXPECT_LT(differ(g.b.beta(), flat.beta(), 0, 1), tol);
        for (int d = 0; d < 3; ++d) { EXPECT_LT(differ(g.b.beta_face(d), flat.beta_face(d), 0, 0), tol) << d; }
        EXPECT_EQ(differ(g.b.wall(), flat.wall(), gust::kWallMask, 1), 0.0);
        EXPECT_LT(differ(g.b.wall(), flat.wall(), gust::kWallD1z, 1), 100 * tol);
        EXPECT_LT(differ(g.b.wall(), flat.wall(), gust::kWallD2z, 1), 100 * tol);
        EXPECT_TRUE(g.b.from_surface());
        EXPECT_FALSE(flat.from_surface());
    }
}

// UT-14.2.2, the helper: on an Agnesi hill filled as the solver fills it, the height above the
// ground is the cell centre's height less the ground's mean over the four nodes of its column.
TEST (ImmersedSurface, HeightAboveTheGroundOnAnAgnesiHill)
{
    for (const int mgs : {16, 4}) {
        const PCase c = gust::lamtest::flat_case(16, 8, 24, 160, 80, 240, true, true, mgs);
        gust::TerrainInputs t = immersed();
        t.type = "agnesi";
        t.hill_shape = "round";
        t.hill_height = 60;
        t.hill_half_width = 40;
        amrex::MultiFab z = with_ground(c, [] (double, double) { return -1.0; });
        gust::fill_immersed_ground(z, c.geom, t);
        amrex::MultiFab dist;
        gust::height_above_ground(z, dist);
        auto hill = [] (int i, int j) {
            const double x = 10.0 * i - 80, y = 10.0 * j - 40;
            return 25 + 60 / (1 + (x * x + y * y) / (40.0 * 40.0));
        };
        const double worst = max_over(dist, [&] (const auto& a, int i, int j, int k) {
            const double ground = 0.25 * (hill(i, j) + hill(i + 1, j) + hill(i, j + 1) + hill(i + 1, j + 1));
            return std::abs(static_cast<double>(a(i, j, k)) - ((k + 0.5) * 10 - ground));
        });
        EXPECT_LT(worst, 100 * tol) << mgs;
        // The ground does not move the mesh: the first component is still the flat levels.
        EXPECT_EQ(static_cast<double>(z.min(0)), 0.0);
        EXPECT_NEAR(static_cast<double>(z.max(1)), 85.0, 100 * tol);
        // The extent check names the two limits.
        EXPECT_EQ(gust::check_immersed_ground_extent(z, 10, 230), "");
        EXPECT_NE(gust::check_immersed_ground_extent(z, 10, 80).find("terrain.flat_height"), std::string::npos);
        EXPECT_NE(gust::check_immersed_ground_extent(z, 40, 230).find("first layer"), std::string::npos);
    }
}

// UT-14.2.2, the consumers: each part that needs the height above the ground reads the helper, the
// helper reads the ground's component, and the parts not yet on it are refused with an immersed ground.
TEST (ImmersedSurface, EveryConsumerReadsTheHelper)
{
    for (const char* file : {"Closure.cpp", "Sequencing.cpp", "Simplec.cpp", "BoundaryConditions.cpp", "BoundaryDamping.cpp", "Initializer.cpp"}) {
        const std::string s = source(file);
        ASSERT_FALSE(s.empty()) << file;
        EXPECT_NE(s.find("height_above_ground"), std::string::npos) << file;
    }
    EXPECT_NE(source("Heights_K.H").find("z.nComp() - 1"), std::string::npos);
    EXPECT_NE(source("Closure.cpp").find("z_nd.nComp() - 1"), std::string::npos);
    // The wall distance of the closure and of the k rows comes from the surface's mean fluid height.
    EXPECT_NE(source("Sequencing.cpp").find("from_surface()"), std::string::npos);
    EXPECT_NE(source("BuildingTke.cpp").find("kSurfaceHeight"), std::string::npos);
    const std::string refusals = source("ImmersedTerrain.cpp");
    for (const char* key : {"canopy.heat = fixed", "steady.probes_file", "init.first_guess = column", "nests"}) {
        EXPECT_NE(refusals.find(key), std::string::npos) << key;
    }
}

// The restart record follows every input of the ground.
TEST (ImmersedSurface, HashFollowsTheGroundsInputs)
{
    if constexpr (kDouble) {
        // A flat ground at its default terrain.small_volfrac: the SHA-256 of
        //   "gust_immersed_terrain 3\nflat_height 0x1.9p+4\nsmall_volfrac 0x1.47ae147ae147bp-8\nz0 0x1.999999999999ap-4\n"
        //   "slope vertical\ntype flat\n",
        // the record of the checkpoints written before the wall-law form was fixed, so that they still restart.
        EXPECT_EQ(gust::immersed_terrain_hash(immersed(), kZ0), "29935c7fa3b12a59feb291a49854a5e532f3962290046cec9c4543b5f66c58b3");
    }
    gust::TerrainInputs t = immersed();
    t.type = "agnesi";
    t.hill_height = 60;
    t.hill_half_width = 40;
    const std::string base = gust::immersed_terrain_hash(t, kZ0);
    EXPECT_EQ(base, gust::immersed_terrain_hash(t, kZ0));
    auto changed = [&] (const std::function<void(gust::TerrainInputs&)>& edit) {
        gust::TerrainInputs e = t;
        edit(e);
        return gust::immersed_terrain_hash(e, kZ0) != base;
    };
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.type = "flat"; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.hill_height = 61; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.hill_half_width = 41; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.hill_shape = "round"; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.flat_height = 26; }));
    // The hybrid mesh: its own line, with the width in metres.
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.representation = "hybrid"; e.hybrid_half_width = 128; }));
    {
        gust::TerrainInputs a = t, b = t;
        a.representation = b.representation = "hybrid";
        a.hybrid_half_width = 128;
        b.hybrid_half_width = 256;
        EXPECT_NE(gust::immersed_terrain_hash(a, kZ0), gust::immersed_terrain_hash(b, kZ0));
    }
    EXPECT_NE(base, gust::immersed_terrain_hash(t, 0.2));

    // A file terrain: the file's contents and the interpolation's inputs.
    const std::string a = ramp_file("hash_a.xyz", 2.0), b = ramp_file("hash_b.xyz", 3.0);
    for (const std::string type : {"xyz", "erf_grid"}) {
        gust::TerrainInputs f = immersed();
        f.type = type;
        f.file = a;
        const std::string hf = gust::immersed_terrain_hash(f, kZ0);
        EXPECT_NE(hf, base) << type;
        gust::TerrainInputs e = f;
        e.file = b;
        EXPECT_NE(hf, gust::immersed_terrain_hash(e, kZ0)) << type << ": the file's contents";
        if (type == "xyz") {
            e = f;
            e.idw_neighbors += 1;
            EXPECT_NE(hf, gust::immersed_terrain_hash(e, kZ0)) << "terrain.idw_neighbors";
            e = f;
            e.idw_power = amrex::Real(3);
            EXPECT_NE(hf, gust::immersed_terrain_hash(e, kZ0)) << "terrain.idw_power";
        }
    }
    gust::TerrainInputs x = immersed(), g = immersed();
    x.type = "xyz";
    g.type = "erf_grid";
    x.file = g.file = a;
    EXPECT_NE(gust::immersed_terrain_hash(x, kZ0), gust::immersed_terrain_hash(g, kZ0));
}

// A ground from a file (terrain.type = xyz): the second component is terrain.flat_height plus the
// file's heights, which inverse-distance weighting returns exactly at the file's own points.
TEST (ImmersedSurface, GroundFromAnXyzFile)
{
    const PCase c = gust::lamtest::flat_case(16, 8, 24, 160, 80, 240, true, true, 8);
    gust::TerrainInputs t = immersed();
    t.type = "xyz";
    t.file = ramp_file("ground.xyz", 2.0);
    amrex::MultiFab z = with_ground(c, [] (double, double) { return -1.0; });
    gust::fill_immersed_ground(z, c.geom, t);
    // Away from the periodic seam (nodes 16 and 8 are the images of node 0).
    const double worst = max_over(z, [&] (const auto& a, int i, int j, int k) {
        if (i >= 16 || j >= 8) { return 0.0; }
        return std::abs(static_cast<double>(a(i, j, k, 1)) - (25 + 2.0 * i + j));
    });
    EXPECT_LT(worst, 100 * tol);
    EXPECT_EQ(static_cast<double>(z.min(0)), 0.0) << "the mesh's levels stay flat";
}

// The hybrid mesh (terrain.representation = hybrid, HybridTerrain.H): the same tilted plane on a mesh whose
// levels slope at s, z_k = s x + 10 k. The shear z' = z - s x keeps areas and vertical distances, and maps
// each cell to a 10 m x 10 m square and the ground to the plane h0 + (m - s) x: the solid fraction, the
// fluid's height above the ground and the face fractions are those of the flat mesh's closed form with the
// slope m - s. The slope fields stay m, the ground's own. A kernel that measures the ground from the
// level's mean (the flat mesh's) fails this.
TEST (ImmersedSurface, SlopedMeshMatchesTheShearedClosedForm)
{
    constexpr double s = 0.15;
    for (const int mgs : {8, 4}) {
        const PCase c = gust::lamtest::flat_case(8, 8, 24, 80, 80, 240, false, false, mgs);
        for (const double deg : {20.0, 30.0, 45.0}) {
            for (const int dir : {0, 1}) {
                SCOPED_TRACE("boxes of " + std::to_string(mgs) + ", " + std::to_string(deg) + " degrees along " + (dir == 0 ? "x" : "y"));
                const double m = std::tan(deg * kPi / 180);
                gust::TerrainInputs t = immersed();
                t.representation = "hybrid";
                t.hybrid_half_width = 10;
                const gust::BuildingsOptions opt = gust::immersed_terrain_options(t, kZ0);
                ASSERT_TRUE(opt.hybrid_mesh);
                amrex::MultiFab z = with_ground(c, [=] (double x, double y) { return 23.0 + m * (dir == 0 ? x : y); });
                const auto dx = c.geom.CellSizeArray();
                for (amrex::MFIter mfi(z, false); mfi.isValid(); ++mfi) {
                    const auto a = z.array(mfi);
                    amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
                        a(i, j, k, 0) = static_cast<amrex::Real>(s * (dir == 0 ? i : j) * static_cast<double>(dx[dir]) + 10.0 * k);
                    });
                }
                gust::Buildings b;
                b.build_from_surface(opt, "test", c.geom, z, c.ba, c.dm);
                const Plane p{23.0, m - s, 10.0, 10.0,
                              std::max(opt.small_volfrac, 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()))};
                auto along = [dir] (int i, int j) { return dir == 0 ? i : j; };

                EXPECT_LT(max_over(b.beta(), [&] (const auto& a, int i, int j, int k) {
                    return std::abs(static_cast<double>(a(i, j, k)) - p.beta(along(i, j), k));
                }), tol);
                EXPECT_LT(max_over(b.surface(), [&] (const auto& a, int i, int j, int k) {
                    const double e = std::abs(static_cast<double>(a(i, j, k, gust::kSurfaceHeight)) - p.height(along(i, j), k));
                    const double ex = std::abs(static_cast<double>(a(i, j, k, gust::kSurfaceSlopeX)) - (dir == 0 ? m : 0.0));
                    const double ey = std::abs(static_cast<double>(a(i, j, k, gust::kSurfaceSlopeY)) - (dir == 1 ? m : 0.0));
                    return std::max(e, std::max(ex, ey));
                }), 100 * tol);
                const double zface_shift = 2 * 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()) * 252 / ((m - s) * 10);
                for (int d = 0; d < 3; ++d) {
                    EXPECT_LT(max_over(b.beta_face(d), [&] (const auto& a, int i, int j, int k) {
                        const int n = along(i, j);
                        if (k > 23 || i > 7 + (d == 0 ? 1 : 0) || j > 7 + (d == 1 ? 1 : 0)) { return 0.0; }
                        double want = 0;
                        if (d == 2) {
                            want = p.below(n, k);
                        } else if (d == dir) {
                            if (n == 0 || n == 8) { return 0.0; }
                            want = p.along(n, k);
                        } else {
                            if ((d == 0 ? i : j) == 0 || (d == 0 ? i : j) == 8) { return 0.0; }
                            want = p.beside(n, k);
                        }
                        return std::abs(static_cast<double>(a(i, j, k)) - want);
                    }), tol + (d == 2 ? zface_shift : 0.0)) << "faces " << d;
                }
                double cut = 0;
                for (amrex::MFIter mfi(b.beta(), false); mfi.isValid(); ++mfi) {
                    const auto be = b.beta().const_array(mfi);
                    amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) { cut += (be(i, j, k) > 0 && be(i, j, k) < 1) ? 1 : 0; });
                }
                amrex::ParallelDescriptor::ReduceRealSum(cut);
                EXPECT_GE(cut, 64.0);   // a cut cell in every column at least
            }
        }
    }
}

// The hybrid kernels on a flat mesh give the flat mesh's geometry (the level is one height, so the ground
// seen from it is the ground itself, to round-off).
TEST (ImmersedSurface, HybridKernelsOnAFlatMeshMatchImmersed)
{
    const PCase c = gust::lamtest::flat_case(8, 8, 24, 80, 80, 240, false, false, 4);
    auto hill = [] (double x, double y) { return 23.0 + 40.0 / (1 + ((x - 40) * (x - 40) + (y - 35) * (y - 35)) / 400.0); };
    Surface flat(c, hill);
    gust::TerrainInputs t = immersed();
    t.representation = "hybrid";
    gust::BuildingsOptions opt = gust::immersed_terrain_options(t, kZ0);
    const amrex::MultiFab z = with_ground(c, hill);
    gust::Buildings b;
    b.build_from_surface(opt, "test", c.geom, z, c.ba, c.dm);
    auto diff = [&] (const amrex::MultiFab& x, const amrex::MultiFab& y) {
        amrex::MultiFab dd(x.boxArray(), x.DistributionMap(), x.nComp(), 0);
        amrex::MultiFab::Copy(dd, x, 0, 0, x.nComp(), 0);
        amrex::MultiFab::Subtract(dd, y, 0, 0, x.nComp(), 0);
        return static_cast<double>(dd.norminf(0, x.nComp(), amrex::IntVect(0)));
    };
    EXPECT_LT(diff(b.beta(), flat.b.beta()), tol);
    EXPECT_LT(diff(b.surface(), flat.b.surface()), 100 * tol);
    for (int d = 0; d < 3; ++d) { EXPECT_LT(diff(b.beta_face(d), flat.b.beta_face(d)), tol) << d; }
    EXPECT_LT(diff(b.wall(), flat.b.wall()), 100 * tol);
}
