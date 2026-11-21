// The immersed ground as a continuous surface (Buildings::build_from_surface, ImmersedSurface_K.H,
// ImmersedSlope.H; terrain.immersed_slope), tracker row PR-88:
//   UT-14.2.1  the solid fraction, the face fractions, the fluid's mean height above the ground and the
//              wall cells of a tilted plane at 10, 30 and 45 degrees against the closed form (the cell
//              clipped by the plane), with the slope along x and along y;
//   UT-14.2.2  the height above the ground on an Agnesi hill is z - h at every cell, and the parts of
//              the solver that need it read the helper (a search of the sources);
//   UT-14.2.3  terrain.immersed_slope = normal on a flat ground is bit-identical to vertical;
//   beyond the plan: a flat surface gives the fields of PR-87's flat ground; normal on a slope takes
//   the distance along the normal, the speed parallel to the ground (w included) and the slope's area.
#include "BuildingForcing.H"
#include "BuildingTke.H"
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
constexpr double kKappa = 0.41;
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

gust::TerrainInputs immersed (const std::string& slope = "vertical", const std::string& side = "off",
                              const std::string& wall_distance = "mean_height")
{
    gust::TerrainInputs t;
    t.representation = "immersed_forcing";
    t.flat_height = 25;
    t.immersed_slope = slope;
    t.immersed_side_walls = side;
    t.immersed_wall_distance = wall_distance;
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
    Surface (const PCase& c, const std::function<double(double, double)>& h, const std::string& slope = "vertical",
             const std::string& side = "off", const std::string& wall_distance = "mean_height")
    {
        opt = gust::immersed_terrain_options(immersed(slope, side, wall_distance), kZ0);
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
                    const bool wall = own < 1 && gust::wall_on_low_side(p.below(n, k), k < 23 ? p.below(n, k + 1) : 0.0, own, p.small);
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
    for (const char* key : {"canopy.type", "steady.probes_file", "init.first_guess = column", "bc.surface_z0_file", "nests"}) {
        EXPECT_NE(refusals.find(key), std::string::npos) << key;
    }
}

namespace {

// A wind that differs at every face, w included.
void fill_wind (const PCase& c, amrex::MultiFab& u, amrex::MultiFab& v, amrex::MultiFab& w)
{
    u.define(c.mesh.ax.boxArray(), c.dm, 1, 1);
    v.define(c.mesh.ay.boxArray(), c.dm, 1, 1);
    w.define(c.mesh.az.boxArray(), c.dm, 1, 1);
    amrex::MultiFab* f[3] = {&u, &v, &w};
    for (int d = 0; d < 3; ++d) {
        for (amrex::MFIter mfi(*f[d], false); mfi.isValid(); ++mfi) {
            const auto a = f[d]->array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
                const int ii = ((i % 4) + 4) % 4, jj = ((j % 4) + 4) % 4;   // periodic over the 4 cells
                a(i, j, k) = static_cast<amrex::Real>((d == 2 ? 0.3 : 2.0 + d) + 0.11 * ii - 0.07 * jj + 0.05 * k + 0.013 * ii * jj);
            });
        }
    }
}

}   // namespace

// UT-14.2.3.
TEST (ImmersedSurface, NormalOnFlatGroundIsBitIdenticalToVertical)
{
    for (const int mgs : {4, 2}) {
        const PCase c = gust::lamtest::flat_case(4, 4, 24, 40, 40, 240, true, true, mgs);
        for (const double height : {23.0, 20.0}) {
            Surface vert(c, [=] (double, double) { return height; }, "vertical");
            Surface norm(c, [=] (double, double) { return height; }, "normal");
            amrex::MultiFab u, v, w;
            fill_wind(c, u, v, w);
            vert.law.update(u, v, w);
            norm.law.update(u, v, w);
            const amrex::MultiFab& a = vert.law.rate();
            const amrex::MultiFab& b = norm.law.rate();
            ASSERT_EQ(a.nComp(), b.nComp());
            for (int n = 0; n < a.nComp(); ++n) {
                amrex::MultiFab d(a.boxArray(), a.DistributionMap(), 1, 0);
                amrex::MultiFab::Copy(d, a, n, 0, 1, 0);
                amrex::MultiFab::Subtract(d, b, n, 0, 1, 0);
                EXPECT_EQ(static_cast<double>(d.norminf()), 0.0) << "component " << n << " height " << height << " boxes " << mgs;
            }
            EXPECT_GT(static_cast<double>(a.norminf(gust::kBuildingUstar2)), 0.0);
        }
    }
}

// Beyond the plan: on a plane tilted along x and a wind of one speed S parallel to it (u = S cos a,
// w = S sin a), normal gives the law at the distances along the normal (d cos a) with the speed S, and
// the rate on the slope's area (the roof's area over cos a); vertical gives the law at the vertical
// distances with the horizontal speed S cos a. A wind along the normal gives no stress with normal.
TEST (ImmersedSurface, NormalTakesTheSlopesDistanceSpeedAndArea)
{
    const PCase c = gust::lamtest::flat_case(8, 8, 24, 80, 80, 240, false, false, 4);
    const double deg = 30, m = std::tan(deg * kPi / 180), cosa = std::cos(deg * kPi / 180), sina = std::sin(deg * kPi / 180);
    const double S = 4;
    // The slope along x, along y, and oblique (the unit vector (ex, ey) points up the slope).
    for (const auto& e : {std::array<double, 2>{1.0, 0.0}, std::array<double, 2>{0.0, 1.0}, std::array<double, 2>{0.6, 0.8}}) {
    const double ex = e[0], ey = e[1];
    auto h = [=] (double x, double y) { return 23.0 + m * (ex * x + ey * y); };
    Surface vert(c, h, "vertical"), norm(c, h, "normal");
    amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1), w(c.mesh.az.boxArray(), c.dm, 1, 1);
    u.setVal(static_cast<amrex::Real>(S * cosa * ex));
    v.setVal(static_cast<amrex::Real>(S * cosa * ey));
    w.setVal(static_cast<amrex::Real>(S * sina));
    vert.law.update(u, v, w);
    norm.law.update(u, v, w);
    const double dzeta = 10;
    double worst = 0, walls = 0;
    for (const Surface* g : {&vert, &norm}) {
        const bool normal = g == &norm;
        for (amrex::MFIter mfi(g->b.beta(), false); mfi.isValid(); ++mfi) {
            const auto be = g->b.beta().const_array(mfi);
            const auto fz = g->b.beta_face(2).const_array(mfi);
            const auto wl = g->b.wall().const_array(mfi);
            const auto r = g->law.rate().const_array(mfi);
            const auto J = c.mesh.detJ.const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                // Away from the domain's sides, where a face of the cell is not a wind face of the run.
                if (wl(i, j, k, gust::kWallMask) == 0 || i == 0 || i == 7 || j == 0 || j == 7) { return; }
                const double f = normal ? cosa : 1.0;
                const double speed = normal ? S : S * cosa;
                const double d1 = static_cast<double>(wl(i, j, k, gust::kWallD1z)) * f, d2 = static_cast<double>(wl(i, j, k, gust::kWallD2z)) * f;
                const double us2 = static_cast<double>(gust::form7_ustar2_exact(amrex::Real(speed), amrex::Real(speed), be(i, j, k), amrex::Real(d1),
                                                                                amrex::Real(d2), amrex::Real(kZ0), amrex::Real(kKappa)));
                const double area = static_cast<double>(gust::wall_area_fraction(fz(i, j, k), fz(i, j, k + 1), be(i, j, k))) / f;
                const double rate = area * us2 / (speed * static_cast<double>(J(i, j, k)) * dzeta);
                worst = std::max(worst, std::abs(static_cast<double>(r(i, j, k, gust::kBuildingUstar2)) - us2) / us2);
                worst = std::max(worst, std::abs(static_cast<double>(r(i, j, k, 0)) - rate) / rate);
                worst = std::max(worst, std::abs(static_cast<double>(r(i, j, k, 1)) - rate) / rate);
                // The sink on w: with normal only (the wind along the slope has a vertical part).
                worst = std::max(worst, std::abs(static_cast<double>(r(i, j, k, 2)) - (normal ? rate : 0.0)) / rate);
                walls += 1;
            });
        }
    }
    amrex::ParallelDescriptor::ReduceRealMax(worst);
    amrex::ParallelDescriptor::ReduceRealSum(walls);
    EXPECT_LT(worst, kDouble ? 1e-10 : 1e-3) << ex << " " << ey;
    EXPECT_GT(walls, 60.0) << ex << " " << ey;

    // A wind along the normal (-sin a (ex, ey), cos a): no speed parallel to the ground.
    u.setVal(static_cast<amrex::Real>(-S * sina * ex));
    v.setVal(static_cast<amrex::Real>(-S * sina * ey));
    w.setVal(static_cast<amrex::Real>(S * cosa));
    norm.law.update(u, v, w);
    const double us2 = max_over(norm.law.rate(), [&] (const auto& a, int i, int j, int k) {
        return (i == 0 || i == 7 || j == 0 || j == 7) ? 0.0 : static_cast<double>(a(i, j, k, gust::kBuildingUstar2));
    });
    EXPECT_LT(us2, kDouble ? 1e-20 : 1e-8) << ex << " " << ey;
    }
}

// The restart record follows every input of the ground.
// terrain.immersed_side_walls = steps. A ramp along x on 10 m cells and layers: at 0.5 m per metre the ground rises
// 5 m per column, less than a layer, so no cell gets a side wall (the wall fields equal those of off, bit for bit);
// on a flat ground the same. At 1.2 m per metre it rises 12 m per column: a cell that is not solid has a side wall
// exactly where the next column up (+x) is solid in its layer, with the step's distances, and keeps
// its roof bit as with off. The wall is on the high x face (Xhi). With off no cell has a side bit.
TEST (ImmersedSurface, StepSideWallsOnlyWhereTheGroundDropsMoreThanALayer)
{
    for (const int mgs : {8, 4}) {
        const PCase c = gust::lamtest::flat_case(16, 8, 24, 160, 80, 240, false, false, mgs);
        auto equal_walls = [] (const amrex::MultiFab& a, const amrex::MultiFab& b) {
            amrex::MultiFab d(a.boxArray(), a.DistributionMap(), a.nComp(), 0);
            amrex::MultiFab::Copy(d, a, 0, 0, a.nComp(), 0);
            amrex::MultiFab::Subtract(d, b, 0, 0, a.nComp(), 0);
            double m = 0;
            for (int n = 0; n < a.nComp(); ++n) { m = std::max(m, static_cast<double>(d.norminf(n))); }
            return m == 0;
        };
        for (const double slope : {0.0, 0.5}) {
            auto h = [=] (double x, double) { return 21.0 + slope * x; };
            Surface off(c, h), steps(c, h, "vertical", "steps");
            EXPECT_TRUE(equal_walls(off.b.wall(), steps.b.wall())) << "slope " << slope << " boxes " << mgs;
        }
        auto h = [] (double x, double) { return 21.0 + 1.2 * x; };
        Surface off(c, h), steps(c, h, "vertical", "steps");
        const int side = gust::kWallXlo | gust::kWallXhi | gust::kWallYlo | gust::kWallYhi;
        amrex::Long n_side = 0, n_bad = 0, n_off_side = 0, n_roof_changed = 0, n_dist_bad = 0, n_at_face = 0;
        const double dxd = 10;
        for (amrex::MFIter mfi(steps.b.wall(), false); mfi.isValid(); ++mfi) {
            const auto ws = steps.b.wall().const_array(mfi);
            const auto wo = off.b.wall().const_array(mfi);
            const auto br = steps.b.beta().const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const int ms = static_cast<int>(ws(i, j, k, gust::kWallMask)), mo = static_cast<int>(wo(i, j, k, gust::kWallMask));
                const double own = static_cast<double>(br(i, j, k));
                const bool want = own < 1 && br(i + 1, j, k) >= 1;
                if ((mo & side) != 0) { ++n_off_side; }
                if ((ms & gust::kWallRoof) != (mo & gust::kWallRoof)) { ++n_roof_changed; }
                if (((ms & gust::kWallXhi) != 0) != want || (ms & (gust::kWallXlo | gust::kWallYlo | gust::kWallYhi)) != 0) { ++n_bad; }
                if ((ms & gust::kWallXhi) != 0) {
                    ++n_side;
                    // The step's distances: half a column to the cell's fluid, 3/2 to the reference,
                    // none at the low x face (not periodic: the reference would be the ghost image).
                    const double d1 = 0.5 * dxd;
                    const double d2 = (i > 0 && br(i - 1, j, k) < 1) ? 1.5 * dxd : 0.0;
                    if (i == 0) { ++n_at_face; }
                    if (std::abs(static_cast<double>(ws(i, j, k, gust::kWallD1x)) - d1) > tol * dxd
                        || std::abs(static_cast<double>(ws(i, j, k, gust::kWallD2x)) - d2) > tol * dxd) { ++n_dist_bad; }
                }
            });
        }
        for (amrex::Long* v : {&n_side, &n_bad, &n_off_side, &n_roof_changed, &n_dist_bad, &n_at_face}) { amrex::ParallelDescriptor::ReduceLongSum(*v); }
        EXPECT_GT(n_side, 0) << "the steep ramp has steps: the case tests nothing otherwise";
        EXPECT_EQ(n_bad, 0) << "side bits where the next column is not solid, or missing, boxes " << mgs;
        EXPECT_EQ(n_off_side, 0) << "off has no side walls";
        EXPECT_EQ(n_roof_changed, 0) << "the roof bits are those of off";
        EXPECT_EQ(n_dist_bad, 0) << "the side's distances";
        EXPECT_GT(n_at_face, 0) << "no wall at the low x face: the domain-face case is not tested";
    }
}

// terrain.immersed_side_walls = steps_no_slip on the steep ramp of the test above. The step bits are those of steps.
// The marked faces are the x faces between a cell with the high x bit and the solid cell beyond it, and no domain
// face. The tangential fluxes across x faces (v and w) take the factor 2 on a marked face (the solid cell's zero half a
// column away) and the open part elsewhere; steps takes the open part (0 on a step face). The side cells take no wall
// law (their rate on v is 0, with steps it is not), and the implicit couplings take the same factors, a_P growing by
// the added part.
TEST (ImmersedSurface, NoSlipStepsDiffuseToTheSolidCell)
{
    for (const int mgs : {8, 4}) {
      for (const bool rising : {true, false}) {   // steps on the high x faces, then on the low ones
        const PCase c = gust::lamtest::flat_case(16, 8, 24, 160, 80, 240, false, false, mgs);
        auto h = [=] (double x, double) { return 21.0 + 1.2 * (rising ? x : 160 - x); };
        Surface ns(c, h, "vertical", "steps_no_slip"), st(c, h, "vertical", "steps");
        const int nx = c.geom.Domain().bigEnd(0) + 1;
        // The open parts (1 - beta on the faces) and the expected marks, from the wall bits, with one ghost layer.
        std::array<amrex::MultiFab, 3> open;
        amrex::MultiFab mark(ns.b.beta_face(0).boxArray(), c.dm, 1, 1);
        mark.setVal(0);
        for (int d = 0; d < 3; ++d) {
            open[d].define(ns.b.beta_face(d).boxArray(), c.dm, 1, 1);
            open[d].setVal(1);
            amrex::MultiFab::Subtract(open[d], ns.b.beta_face(d), 0, 0, 1, 0);
            open[d].FillBoundary(c.geom.periodicity());
        }
        amrex::Long n_mark = 0, n_mark_bad = 0, n_bits_differ = 0;
        for (amrex::MFIter mfi(mark, false); mfi.isValid(); ++mfi) {
            const auto m = mark.array(mfi);
            const auto f = ns.law.no_slip_faces(0).const_array(mfi);
            const auto wn = ns.b.wall().const_array(mfi);
            const auto ws = st.b.wall().const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                const bool want = (i > 0 && (static_cast<int>(wn(i - 1, j, k, gust::kWallMask)) & gust::kWallXhi) != 0)
                                  || (i < nx && (static_cast<int>(wn(i, j, k, gust::kWallMask)) & gust::kWallXlo) != 0);
                m(i, j, k) = want ? 1 : 0;
                if (want) { ++n_mark; }
                if ((static_cast<double>(f(i, j, k)) == 1.0) != want) { ++n_mark_bad; }
                if (i < nx && wn(i, j, k, gust::kWallMask) != ws(i, j, k, gust::kWallMask)) { ++n_bits_differ; }
            });
        }
        mark.FillBoundary(c.geom.periodicity());
        for (amrex::Long* v : {&n_mark, &n_mark_bad, &n_bits_differ}) { amrex::ParallelDescriptor::ReduceLongSum(*v); }
        EXPECT_GT(n_mark, 0) << "the steep ramp has steps";
        EXPECT_EQ(n_mark_bad, 0) << "the marked faces, boxes " << mgs;
        EXPECT_EQ(n_bits_differ, 0) << "the step bits are those of steps";
        const double y_marks = max_over(ns.law.no_slip_faces(1), [] (const auto& a, int i, int j, int k) { return std::abs(static_cast<double>(a(i, j, k))); });
        EXPECT_EQ(y_marks, 0.0) << "a ramp along x has no y steps";

        // The fluxes: all 1, then closed by each law.
        const amrex::Array<const amrex::MultiFab*, 3> op{&open[0], &open[1], &open[2]};
        auto closed = [&] (const Surface& s) {
            amrex::Array<amrex::Array<amrex::MultiFab, 3>, 3> G;
            for (int a = 0; a < 3; ++a) {
                for (int d = 0; d < 3; ++d) {
                    amrex::IntVect t(0);
                    t[a] = 1;
                    t[d] = 1;
                    if (a == d) { t = amrex::IntVect(0); }
                    G[a][d].define(amrex::convert(c.ba, t), c.dm, 1, 0);
                    G[a][d].setVal(1);
                }
            }
            s.law.close_viscous_fluxes(op, G);
            return G;
        };
        auto Gn = closed(ns);
        auto Gs = closed(st);
        amrex::Long n_two = 0, n_flux_bad = 0;
        for (const int a : {1, 2}) {
            const int di = (a == 1) ? 0 : 0, dj = (a == 1) ? 1 : 0, dk = (a == 2) ? 1 : 0;
            for (amrex::MFIter mfi(Gn[a][0], false); mfi.isValid(); ++mfi) {
                const auto gn = Gn[a][0].const_array(mfi);
                const auto gs = Gs[a][0].const_array(mfi);
                const auto o = open[0].const_array(mfi);
                const auto m = mark.const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    auto fac = [&] (int ii, int jj, int kk) { return m(ii, jj, kk) == 1 ? 2.0 : static_cast<double>(o(ii, jj, kk)); };
                    const double want_n = 0.5 * (fac(i - di, j - dj, k - dk) + fac(i, j, k));
                    const double want_s = 0.5 * (static_cast<double>(o(i - di, j - dj, k - dk)) + static_cast<double>(o(i, j, k)));
                    if (want_n == 2.0) { ++n_two; }
                    if (std::abs(static_cast<double>(gn(i, j, k)) - want_n) > tol || std::abs(static_cast<double>(gs(i, j, k)) - want_s) > tol) { ++n_flux_bad; }
                });
            }
        }
        for (amrex::Long* v : {&n_two, &n_flux_bad}) { amrex::ParallelDescriptor::ReduceLongSum(*v); }
        EXPECT_GT(n_two, 0) << "edges between two step faces: the factor 2 is not tested otherwise";
        EXPECT_EQ(n_flux_bad, 0) << "the tangential fluxes across x faces, boxes " << mgs;

        // The rates: the side cells take no law with no slip, so the rates are those of off (the roofs' alone);
        // with steps the side cells' rate on v exceeds off's.
        Surface off(c, h);
        amrex::MultiFab u(c.mesh.ax.boxArray(), c.dm, 1, 1), v(c.mesh.ay.boxArray(), c.dm, 1, 1), w(c.mesh.az.boxArray(), c.dm, 1, 1);
        u.setVal(0);
        v.setVal(2);
        w.setVal(1);
        ns.law.update(u, v, w);
        st.law.update(u, v, w);
        off.law.update(u, v, w);
        double ns_off = 0, st_gain = 1e300;
        amrex::Long n_side_cells = 0;
        for (amrex::MFIter mfi(off.law.rate(), false); mfi.isValid(); ++mfi) {
            const auto rn = ns.law.rate().const_array(mfi);
            const auto rs = st.law.rate().const_array(mfi);
            const auto ro = off.law.rate().const_array(mfi);
            const auto wl = st.b.wall().const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                for (int n = 0; n < 3; ++n) { ns_off = std::max(ns_off, std::abs(static_cast<double>(rn(i, j, k, n) - ro(i, j, k, n)))); }
                if ((static_cast<int>(wl(i, j, k, gust::kWallMask)) & (gust::kWallXlo | gust::kWallXhi)) != 0) {
                    ++n_side_cells;
                    st_gain = std::min(st_gain, static_cast<double>(rs(i, j, k, 1) - ro(i, j, k, 1)));
                }
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(ns_off);
        amrex::ParallelDescriptor::ReduceRealMin(st_gain);
        amrex::ParallelDescriptor::ReduceLongSum(n_side_cells);
        EXPECT_GT(n_side_cells, 0);
        EXPECT_EQ(ns_off, 0.0) << "no wall law on a no-slip step: the rates are off's";
        EXPECT_GT(st_gain, 0.0) << "steps: the law acts on every side cell";

        // The implicit couplings of v: XM and XP take the factors, a_P grows by the part above 1.
        auto implicit = [&] (const Surface& s) {
            amrex::MultiFab A(c.mesh.ay.boxArray(), c.dm, gust::n_implicit_coefs, 0);
            A.setVal(1);
            A.setVal(6, gust::coef_P, 1);
            s.law.add_to_implicit(1, op, A);
            return A;
        };
        const amrex::MultiFab An = implicit(ns), As = implicit(st);
        amrex::Long n_row_bad = 0, n_rows = 0;
        for (amrex::MFIter mfi(An, false); mfi.isValid(); ++mfi) {
            const auto an = An.const_array(mfi);
            const auto as = As.const_array(mfi);
            const auto o = open[0].const_array(mfi);
            const auto m = mark.const_array(mfi);
            const auto rn = ns.law.rate().const_array(mfi, 1);
            const auto rs = st.law.rate().const_array(mfi, 1);
            const auto bf = ns.b.beta_face(1).const_array(mfi);
            amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                if (static_cast<double>(bf(i, j, k)) == 1.0) { return; }   // a blocked face: the identity in both
                auto fac = [&] (int ii, int jj, int kk) { return m(ii, jj, kk) == 1 ? 2.0 : static_cast<double>(o(ii, jj, kk)); };
                const double flo = 0.5 * (fac(i, j - 1, k) + fac(i, j, k));
                const double fhi = 0.5 * (fac(i + 1, j - 1, k) + fac(i + 1, j, k));
                const double dr = 0.5 * (static_cast<double>(rn(i, j - 1, k) + rn(i, j, k)) - static_cast<double>(rs(i, j - 1, k) + rs(i, j, k)));
                const double dP = dr + std::max(flo - 1, 0.0) + std::max(fhi - 1, 0.0);
                if (flo > 1 || fhi > 1) { ++n_rows; }
                if (std::abs(static_cast<double>(an(i, j, k, gust::coef_XM)) - flo) > tol || std::abs(static_cast<double>(an(i, j, k, gust::coef_XP)) - fhi) > tol
                    || std::abs(static_cast<double>(an(i, j, k, gust::coef_P) - as(i, j, k, gust::coef_P)) - dP) > 10 * tol) { ++n_row_bad; }
            });
        }
        for (amrex::Long* v : {&n_rows, &n_row_bad}) { amrex::ParallelDescriptor::ReduceLongSum(*v); }
        EXPECT_GT(n_rows, 0) << "rows next to a step face";
        EXPECT_EQ(n_row_bad, 0) << "the implicit couplings of v, boxes " << mgs << (rising ? " rising" : " falling");
      }
    }
}

// terrain.immersed_wall_distance = column_centre. Ramps along x on 10 m cells and layers. At 0.5 m per metre the
// ground spans 5 m across a footprint, less than the layer: the distances equal mean_height's, bit for bit, and so
// on a flat ground. At 1.2 m per metre it spans 12 m: in every cell that is not solid the distance is
// min(d_mean, max(z_c - h_c, dz / 2)), with z_c the cell centre's height and h_c = 21 + 1.2 x_c the ground's at the
// column centre (closed form on a straight ramp); somewhere the column centre gives the smaller value. With and
// without step side walls, on one box and several.
TEST (ImmersedSurface, ColumnCentreDistanceOnlyWhereTheFootprintSpansMoreThanALayer)
{
    const amrex::Real floor = amrex::Real(1e-6), cmu0 = amrex::Real(0.5477);
    auto dist = [&] (const PCase& c, const Surface& s) {
        gust::BuildingTke t;
        t.build(s.b, s.law, s.opt, c.geom, s.z, floor, cmu0);
        amrex::MultiFab d(c.ba, c.dm, 1, 0);
        amrex::MultiFab::Copy(d, t.distance(), 0, 0, 1, 0);
        return d;
    };
    for (const int mgs : {8, 4}) {
        const PCase c = gust::lamtest::flat_case(16, 8, 24, 160, 80, 240, false, false, mgs);
        for (const double slope : {0.0, 0.5}) {
            auto h = [=] (double x, double) { return 21.0 + slope * x; };
            const Surface mh(c, h), cc(c, h, "vertical", "off", "column_centre");
            amrex::MultiFab a = dist(c, mh), b = dist(c, cc);
            amrex::MultiFab::Subtract(a, b, 0, 0, 1, 0);
            EXPECT_EQ(static_cast<double>(a.norminf(0)), 0.0) << "slope " << slope << " boxes " << mgs;
        }
        for (const std::string side : {"off", "steps"}) {
            auto h = [] (double x, double) { return 21.0 + 1.2 * x; };
            const Surface mh(c, h, "vertical", side), cc(c, h, "vertical", side, "column_centre");
            const amrex::MultiFab dm = dist(c, mh), dc = dist(c, cc);
            amrex::Long n_fluid = 0, n_bad = 0, n_smaller = 0;
            for (amrex::MFIter mfi(dc, false); mfi.isValid(); ++mfi) {
                const auto am = dm.const_array(mfi);
                const auto ac = dc.const_array(mfi);
                const auto bt = cc.b.beta().const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    if (static_cast<double>(bt(i, j, k)) >= 1) {
                        if (ac(i, j, k) != am(i, j, k)) { ++n_bad; }   // solid cells keep theirs
                        return;
                    }
                    ++n_fluid;
                    const double zc = 10.0 * (k + 0.5), hc = 21.0 + 1.2 * 10.0 * (i + 0.5);
                    const double want = std::min(static_cast<double>(am(i, j, k)), std::max(zc - hc, 5.0));
                    if (want < static_cast<double>(am(i, j, k)) - 1e-3) { ++n_smaller; }
                    if (std::abs(static_cast<double>(ac(i, j, k)) - want) > tol * 100) { ++n_bad; }
                });
            }
            for (amrex::Long* v : {&n_fluid, &n_bad, &n_smaller}) { amrex::ParallelDescriptor::ReduceLongSum(*v); }
            EXPECT_GT(n_fluid, 0);
            EXPECT_GT(n_smaller, 0) << "the column centre is nowhere nearer: the case tests nothing";
            EXPECT_EQ(n_bad, 0) << "the column-centre distance, side walls " << side << " boxes " << mgs;
        }
    }
}

TEST (ImmersedSurface, HashFollowsTheGroundsInputs)
{
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
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.immersed_slope = "normal"; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.immersed_side_walls = "steps"; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.immersed_side_walls = "steps_no_slip"; }));
    {
        gust::TerrainInputs a = t, b = t;
        a.immersed_side_walls = "steps";
        b.immersed_side_walls = "steps_no_slip";
        EXPECT_NE(gust::immersed_terrain_hash(a, kZ0), gust::immersed_terrain_hash(b, kZ0));
    }
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.immersed_wall_distance = "column_centre"; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.type = "flat"; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.hill_height = 61; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.hill_half_width = 41; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.hill_shape = "round"; }));
    EXPECT_TRUE(changed([] (gust::TerrainInputs& e) { e.flat_height = 26; }));
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
