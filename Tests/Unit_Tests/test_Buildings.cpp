// The building geometry (Buildings.H, Buildings_K.H), beyond the plan,
// sub-phase 5.1 (UT-5.1.x):
//   a box on flat ground: cell and face fractions are the exact overlap, part-filled
//   cells included (UT-5.1.1); on a slope, each column's solid volume is the
//   map's mean, from the local ground (UT-5.1.2); masks, d1 and d2 of roof, side,
//   edge and corner cells equal hand values (UT-5.1.3); the same bits on every
//   box layout, with a building across a periodic seam (UT-5.1.4); bare columns
//   hold nothing, the ground included (UT-5.1.5); ERF's own cube file: a core and
//   the rim's volume as ERF's interpolated heights (UT-5.1.6, closed form);
//   the two file formats agree; inputs and files are checked.
#include "Buildings.H"
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
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;
using HFn = std::function<double(double, double)>;

constexpr double rtol = std::is_same_v<amrex::Real, double> ? 1e-11 : 2e-4;

// n x n x 24 cells over L x L x 240 m: 10 m layers, boxes of mgs cells.
PCase make_case (int n, double L, int mgs, bool periodic_x, const HFn& ground)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(n, n, 24);
    c.in.mesh.prob_hi = {amrex::Real(L), amrex::Real(L), 240};
    c.in.mesh.is_periodic = {periodic_x ? 1 : 0, 0, 0};
    c.in.mesh.max_grid_size = {mgs, mgs};
    c.in.terrain.smoothing = 0;
    c.make(ground);
    return c;
}

double flat (double, double) { return 0.0; }

std::string file_path (const std::string& name)
{
    // One folder per run: CTest may run two unit jobs of the same size at once.
    static const long id = [] () {
        long pid = static_cast<long>(::getpid());
        amrex::ParallelDescriptor::Bcast(&pid, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        return pid;
    }();
    const fs::path dir = fs::temp_directory_path() / ("gust_buildings_" + std::to_string(id));
    if (amrex::ParallelDescriptor::IOProcessor()) { fs::create_directories(dir); }
    amrex::ParallelDescriptor::Barrier();
    return (dir / name).string();
}

std::string write_text (const std::string& name, const std::string& text)
{
    const std::string file = file_path(name);
    if (amrex::ParallelDescriptor::IOProcessor()) { std::ofstream(file) << text; }
    amrex::ParallelDescriptor::Barrier();
    return file;
}

// nt x nt points x0 + d t in ERF's format.
std::string write_erf (const std::string& name, int nt, double x0, double d, const HFn& h)
{
    const std::string file = file_path(name);
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(file);
        os.precision(17);
        os << nt << "\n" << nt << "\n";
        for (int rep = 0; rep < 2; ++rep) { for (int t = 0; t < nt; ++t) { os << x0 + d * t << "\n"; } }
        for (int i = 0; i < nt; ++i) { for (int j = 0; j < nt; ++j) { os << h(x0 + d * i, x0 + d * j) << "\n"; } }
    }
    amrex::ParallelDescriptor::Barrier();
    return file;
}

gust::BuildingsOptions options (const std::string& file, double small = 0.005)
{
    gust::BuildingsOptions o;
    o.type = gust::BuildingsType::immersed_forcing;
    o.type_name = "immersed_forcing";
    o.file = file;
    o.small_volfrac = small;
    return o;
}

gust::Buildings make_buildings (const PCase& c, const gust::BuildingsOptions& o)
{
    gust::Buildings b;
    b.build(o, c.geom, c.z, c.ba, c.dm);
    return b;
}

// A field's valid data on every rank, indexable by (i, j, k, n).
struct Gathered
{
    amrex::Box box;
    int ncomp{0};
    std::vector<double> v;
    [[nodiscard]] double operator() (int i, int j, int k, int n = 0) const
    {
        const auto l = box.length();
        return v[static_cast<std::size_t>((i - box.smallEnd(0)) + l[0] * ((j - box.smallEnd(1)) + static_cast<long>(l[1]) * ((k - box.smallEnd(2)) + static_cast<long>(l[2]) * n)))];
    }
};

Gathered gather (const amrex::MultiFab& f, const amrex::Geometry& geom)
{
    Gathered g;
    g.box = amrex::convert(geom.Domain(), f.ixType());
    g.ncomp = f.nComp();
    amrex::BoxArray ba(g.box);
    amrex::DistributionMapping dm(amrex::Vector<int>{amrex::ParallelDescriptor::IOProcessorNumber()});
    amrex::MultiFab one(ba, dm, g.ncomp, 0, amrex::MFInfo().SetArena(amrex::The_Pinned_Arena()));
    one.setVal(0);
    one.ParallelCopy(f, 0, 0, g.ncomp);
    amrex::Gpu::streamSynchronize();
    g.v.resize(static_cast<std::size_t>(g.box.numPts()) * static_cast<std::size_t>(g.ncomp));
    for (amrex::MFIter mfi(one, false); mfi.isValid(); ++mfi) {
        const amrex::Real* p = one[mfi].dataPtr();
        for (std::size_t n = 0; n < g.v.size(); ++n) { g.v[n] = static_cast<double>(p[n]); }
    }
    amrex::ParallelDescriptor::Bcast(g.v.data(), g.v.size(), amrex::ParallelDescriptor::IOProcessorNumber());
    return g;
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

double overlap (double a, double b, double c, double d) { return std::max(0.0, std::min(b, d) - std::max(a, c)); }

// The tiles of the tests: 5 m, covering 0 to 400 m.
HFn box_map (double xa, double xb, double ya, double yb, double h)
{
    return [=] (double x, double y) { return (x > xa && x < xb && y > ya && y < yb) ? h : 0.0; };
}

} // namespace

TEST (Buildings, FlatBoxFractionsAreTheExactOverlap)
{
    const PCase c = make_case(16, 400, 8, false, flat);
    const std::string file = write_erf("box.txt", 80, 2.5, 5, box_map(110, 190, 140, 235, 47));
    const gust::Buildings b = make_buildings(c, options(file));
    const Gathered beta = gather(b.beta(), c.geom);
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            for (int k = 0; k < 24; ++k) {
                const double want = overlap(i * 25., (i + 1) * 25., 110, 190) / 25 * overlap(j * 25., (j + 1) * 25., 140, 235) / 25
                                    * overlap(k * 10., (k + 1) * 10., 0, 47) / 10;
                EXPECT_NEAR(beta(i, j, k), want, rtol) << i << " " << j << " " << k;
            }
        }
    }
    // Hand values: x = 125 m is inside the footprint, 100 m outside; z = 40 m lies on a solid cell, z = 50 m above the roof.
    const Gathered fx = gather(b.beta_face(0), c.geom), fy = gather(b.beta_face(1), c.geom), fz = gather(b.beta_face(2), c.geom);
    EXPECT_NEAR(fx(5, 6, 4), 0.7, rtol);
    EXPECT_EQ(fx(5, 6, 2), 1.0);
    EXPECT_EQ(fx(4, 6, 2), 0.0);
    EXPECT_NEAR(fy(5, 6, 4), 0.7, rtol);
    EXPECT_NEAR(fy(4, 6, 4), 0.6 * 0.7, rtol);
    EXPECT_EQ(fz(5, 6, 4), 1.0);
    EXPECT_EQ(fz(5, 6, 5), 0.0);
    EXPECT_EQ(fz(5, 6, 0), 1.0);
    EXPECT_NEAR(fz(4, 6, 2), 0.6, rtol);
}

TEST (Buildings, SlopedColumnsHoldTheMapsVolumeFromTheLocalGround)
{
    const PCase c = make_case(16, 400, 6, false, [] (double x, double y) { return 0.05 * x + 0.02 * y; });
    const std::string file = write_erf("slope.txt", 80, 2.5, 5, box_map(110, 190, 140, 235, 47));
    const gust::Buildings b = make_buildings(c, options(file, 0));
    const Gathered beta = gather(b.beta(), c.geom), J = gather(c.mesh.detJ, c.geom);
    const double dzeta = static_cast<double>(c.geom.CellSize(2));
    for (int i = 0; i < 16; ++i) {
        for (int j = 0; j < 16; ++j) {
            double sum = 0;
            for (int k = 0; k < 24; ++k) { sum += beta(i, j, k) * J(i, j, k) * dzeta; }
            const double want = overlap(i * 25., (i + 1) * 25., 110, 190) / 25 * overlap(j * 25., (j + 1) * 25., 140, 235) / 25 * 47;
            EXPECT_NEAR(sum, want, 100 * rtol) << i << " " << j;
        }
    }
}

TEST (Buildings, WallMasksAndDistancesEqualHandValues)
{
    using namespace gust;
    const PCase c = make_case(16, 400, 8, false, flat);
    // A 100 m cube footprint on cell edges (i, j = 4 to 7), solid to k = 3, half of k = 4.
    Buildings b = make_buildings(c, options(write_erf("cube.txt", 80, 2.5, 5, box_map(100, 200, 100, 200, 45))));
    Gathered w = gather(b.wall(), c.geom);
    auto mask = [&] (int i, int j, int k) { return static_cast<int>(w(i, j, k, kWallMask)); };
    EXPECT_EQ(mask(3, 5, 2), kWallXhi);
    EXPECT_NEAR(w(3, 5, 2, kWallD1x), 12.5, rtol);
    EXPECT_NEAR(w(3, 5, 2, kWallD2x), 37.5, rtol);
    EXPECT_EQ(w(3, 5, 2, kWallD1y), 0.0);
    EXPECT_EQ(mask(8, 5, 2), kWallXlo);
    EXPECT_EQ(mask(5, 3, 2), kWallYhi);
    EXPECT_EQ(mask(5, 8, 2), kWallYlo);
    EXPECT_EQ(mask(5, 5, 4), kWallRoof);
    EXPECT_NEAR(w(5, 5, 4, kWallD1z), 2.5, rtol);
    EXPECT_NEAR(w(5, 5, 4, kWallD2z), 10.0, rtol);
    EXPECT_EQ(mask(3, 5, 4), kWallXhi);   // fluid beside the wall from 40 to 45 m
    EXPECT_EQ(mask(4, 5, 4), kWallRoof);  // the wall lies on its low face, so the fluid cell has it
    EXPECT_EQ(w(4, 5, 4, kWallD1x), 0.0);
    EXPECT_EQ(mask(5, 5, 2), 0);   // solid
    EXPECT_EQ(mask(5, 5, 5), 0);   // above a partial cell
    EXPECT_EQ(mask(4, 5, 0), 0);   // solid at the ground: no underside
    // A roof on a layer face: the fluid cell above the solid one is the wall cell.
    b = make_buildings(c, options(write_erf("cube40.txt", 80, 2.5, 5, box_map(100, 200, 100, 200, 40))));
    w = gather(b.wall(), c.geom);
    EXPECT_EQ(mask(5, 5, 4), kWallRoof);
    EXPECT_NEAR(w(5, 5, 4, kWallD1z), 5.0, rtol);
    EXPECT_NEAR(w(5, 5, 4, kWallD2z), 15.0, rtol);
    // A half-cell rim: x from 112.5 m (i = 4 half full).
    b = make_buildings(c, options(write_erf("rim.txt", 160, 1.25, 2.5, box_map(112.5, 200, 100, 200, 45))));
    w = gather(b.wall(), c.geom);
    EXPECT_EQ(mask(4, 5, 2), kWallXhi);
    EXPECT_NEAR(w(4, 5, 2, kWallD1x), 6.25, rtol);
    EXPECT_EQ(mask(3, 5, 2), 0);
    EXPECT_EQ(mask(4, 5, 0), kWallXhi);   // the rim at the ground: no underside
    EXPECT_EQ(mask(4, 5, 4), kWallRoof | kWallXhi);   // edge
    EXPECT_NEAR(w(4, 5, 4, kWallD1x), 9.375, rtol);
    // Rims in x and y: a corner.
    b = make_buildings(c, options(write_erf("corner.txt", 160, 1.25, 2.5, box_map(112.5, 200, 112.5, 200, 45))));
    w = gather(b.wall(), c.geom);
    EXPECT_EQ(mask(4, 4, 2), kWallXhi | kWallYhi);
    EXPECT_EQ(mask(4, 4, 4), kWallRoof | kWallXhi | kWallYhi);
    // A step in the roof: the wall from 45 to 50 m on the partial cell beside the solid one.
    b = make_buildings(c, options(write_erf("step.txt", 80, 2.5, 5, [] (double x, double y) {
        return (x > 100 && x < 200 && y > 100 && y < 200) ? (x < 150 ? 45.0 : 55.0) : 0.0;
    })));
    w = gather(b.wall(), c.geom);
    EXPECT_EQ(mask(5, 5, 4), kWallRoof | kWallXhi);
    EXPECT_NEAR(w(5, 5, 4, kWallD2x), 25.0, rtol);
    // A one-cell street: both sides, no usable reference cell.
    b = make_buildings(c, options(write_erf("street.txt", 80, 2.5, 5, [] (double x, double y) {
        return ((x > 100 && x < 125) || (x > 150 && x < 175)) && y > 100 && y < 200 ? 60.0 : 0.0;
    })));
    w = gather(b.wall(), c.geom);
    EXPECT_EQ(mask(5, 5, 2), kWallXlo | kWallXhi);
    EXPECT_NEAR(w(5, 5, 2, kWallD1x), 12.5, rtol);
    EXPECT_EQ(w(5, 5, 2, kWallD2x), 0.0);
    // A building thinner than a cell: no net wall, and the log says so.
    b = make_buildings(c, options(write_erf("fin.txt", 160, 1.25, 2.5, box_map(110, 115, 100, 200, 45))));
    EXPECT_NE(b.summary(c.mesh.detJ, c.geom).find("thinner than a cell"), std::string::npos);
    // A building on a non-periodic domain side is no wall there, and the log says so.
    b = make_buildings(c, options(write_erf("side.txt", 80, 2.5, 5, box_map(0, 100, 100, 200, 45))));
    w = gather(b.wall(), c.geom);
    EXPECT_EQ(mask(0, 5, 2), 0);
    EXPECT_EQ(mask(4, 5, 2), kWallXlo);
    EXPECT_NE(b.summary(c.mesh.detJ, c.geom).find("non-periodic side face"), std::string::npos);
}

TEST (Buildings, SameBitsOnEveryLayoutAcrossAPeriodicSeam)
{
    using namespace gust;
    // Blocks at both ends of x (the same building through the seam); one tile straddles x = 0.
    const HFn h = [] (double x, double y) { return ((x < 30 || x > 360) && y > 100 && y < 200) ? 45.0 : 0.0; };
    const std::string file = write_erf("seam.txt", 80, 0, 5, h);
    const PCase c1 = make_case(16, 400, 16, true, flat), c2 = make_case(16, 400, 5, true, flat);
    const Buildings b1 = make_buildings(c1, options(file)), b2 = make_buildings(c2, options(file));
    const Gathered a1 = gather(b1.beta(), c1.geom), a2 = gather(b2.beta(), c2.geom);
    EXPECT_EQ(a1.v, a2.v);
    EXPECT_EQ(gather(b1.wall(), c1.geom).v, gather(b2.wall(), c2.geom).v);
    for (int d = 0; d < 3; ++d) { EXPECT_EQ(gather(b1.beta_face(d), c1.geom).v, gather(b2.beta_face(d), c2.geom).v) << d; }
    // Tiles reach from x = 362.5 m through 400 m to 27.5 m.
    EXPECT_EQ(a1(15, 5, 2), 1.0);
    EXPECT_EQ(a1(0, 5, 2), 1.0);
    EXPECT_NEAR(a1(1, 5, 2), 2.5 / 25, rtol);
    EXPECT_NEAR(a1(14, 5, 2), 12.5 / 25, rtol);
    const Gathered w = gather(b1.wall(), c1.geom), fx = gather(b1.beta_face(0), c1.geom);
    EXPECT_EQ(static_cast<int>(w(0, 5, 2, kWallMask)), 0);   // no wall at the seam
    EXPECT_EQ(static_cast<int>(w(15, 5, 2, kWallMask)), 0);
    EXPECT_EQ(fx(0, 5, 2), 1.0);
    EXPECT_EQ(fx(16, 5, 2), 1.0);
    // Both ends of x listed: the seam twice.
    const std::string twice = write_erf("twice.txt", 81, 0, 5, h);
    EXPECT_NE(abort_message([&] () { make_buildings(c1, options(twice)); }).find("wider than the periodic domain"), std::string::npos);
    const std::string far = write_erf("far.txt", 80, 700, 5, [] (double, double) { return 10.0; });
    EXPECT_NE(abort_message([&] () { make_buildings(c1, options(far)); }).find("more than one period"), std::string::npos);
}

TEST (Buildings, BareColumnsHoldNothingTheGroundIncluded)
{
    const PCase c = make_case(16, 400, 8, false, [] (double x, double) { return 0.05 * x; });
    const gust::Buildings b = make_buildings(c, options(write_erf("bare.txt", 80, 2.5, 5, box_map(100, 200, 100, 200, 45))));
    const Gathered beta = gather(b.beta(), c.geom), fz = gather(b.beta_face(2), c.geom), w = gather(b.wall(), c.geom);
    for (int k = 0; k < 24; ++k) {
        EXPECT_EQ(beta(12, 12, k), 0.0) << k;
        EXPECT_EQ(fz(12, 12, k), 0.0) << k;
        EXPECT_EQ(w(12, 12, k, gust::kWallMask), 0.0) << k;
    }
    gust::BuildingsOptions none;
    EXPECT_FALSE(none.active());
    EXPECT_TRUE(gust::check_buildings_options(none).empty());
    none.given.emplace_back("small_volfrac");
    EXPECT_EQ(gust::check_buildings_options(none).size(), 1u);
}

// ERF's PartialCells cube (make_buildings.py): a 40 m cube on 10 m nodes, 140 to
// 180 m. ERF's reader takes heights at the cell corners: a 4-cell core, rims at half height,
// corners at a quarter. The columns here hold the same volumes.
TEST (Buildings, ErfCubeFileGivesTheColumnVolumesOfItsInterpolatedHeights)
{
    const PCase c = make_case(32, 320, 16, true, flat);
    const std::string file = write_erf("erf_cube.txt", 32, 0, 10, [] (double x, double y) {
        return (x >= 140 && x <= 180 && y >= 140 && y <= 180) ? 40.0 : 0.0;
    });
    const gust::Buildings b = make_buildings(c, options(file, 0));
    const Gathered beta = gather(b.beta(), c.geom);
    auto column = [&] (int i, int j) {
        double s = 0;
        for (int k = 0; k < 24; ++k) { s += 10 * beta(i, j, k); }
        return s;
    };
    EXPECT_NEAR(column(14, 14), 40, 100 * rtol);
    EXPECT_NEAR(column(17, 17), 40, 100 * rtol);
    EXPECT_NEAR(column(13, 15), 20, 100 * rtol);
    EXPECT_NEAR(column(18, 15), 20, 100 * rtol);
    EXPECT_NEAR(column(13, 13), 10, 100 * rtol);
    EXPECT_NEAR(column(12, 15), 0, 100 * rtol);
}

TEST (Buildings, TheTwoFileFormatsGiveTheSameField)
{
    const PCase c = make_case(16, 400, 8, false, flat);
    const HFn h = box_map(110, 190, 140, 235, 47);
    const std::string erf = write_erf("same.txt", 80, 2.5, 5, h);
    std::string text = "# x y height, in no order\n";
    {
        std::ostringstream os;
        for (int i = 79; i >= 0; --i) {
            for (int j = 0; j < 80; ++j) { os << 2.5 + 5 * i << " " << 2.5 + 5 * j << " " << h(2.5 + 5 * i, 2.5 + 5 * j) << "\n"; }
        }
        text += os.str();
    }
    gust::BuildingsOptions o = options(write_text("same.xyz", text));
    o.format = "xyz";
    EXPECT_EQ(gather(make_buildings(c, options(erf)).beta(), c.geom).v, gather(make_buildings(c, o).beta(), c.geom).v);
    auto refused = [&] (const std::string& body) {
        return abort_message([&] () { gust::read_building_map(write_text("bad.xyz", body), "xyz"); });
    };
    EXPECT_NE(refused("0 0 1\n1 0 1\n0 1 1\n").find("full grid of 4"), std::string::npos);
    EXPECT_NE(refused("0 0 1\n1 0 1\n0 1 1\n0 1 2\n").find("appears twice"), std::string::npos);
    EXPECT_NE(refused("0 0 1\n1 0 1\n0 1 1\n1 1 -2\n").find("must be finite and >= 0"), std::string::npos);
    EXPECT_NE(refused("0 0 1\n0 1 1\n").find("at least 2 x"), std::string::npos);
}

// x y height zstart: a lifted 100 m block, 10 m layers.
gust::Buildings lifted (const PCase& c, const std::string& name, double h, double zs)
{
    std::ostringstream os;
    for (int i = 0; i < 80; ++i) {
        for (int j = 0; j < 80; ++j) {
            const double x = 2.5 + 5 * i, y = 2.5 + 5 * j;
            const bool in = x > 100 && x < 200 && y > 100 && y < 200;
            os << x << ", " << y << ", " << (in ? h : 0.0) << ", " << (in ? zs : 0.0) << "\n";
        }
    }
    gust::BuildingsOptions o = options(write_text(name, os.str()));
    o.format = "xyz";
    return make_buildings(c, o);
}

TEST (Buildings, VoidDeckIsOpenWithAWallUnderTheBuilding)
{
    using namespace gust;
    const PCase c = make_case(16, 400, 8, false, flat);
    // Underside inside k = 2, roof inside k = 6.
    Buildings b = lifted(c, "void25.xyz", 65, 25);
    Gathered beta = gather(b.beta(), c.geom), w = gather(b.wall(), c.geom);
    Gathered fx = gather(b.beta_face(0), c.geom), fz = gather(b.beta_face(2), c.geom);
    const double want[] = {0, 0, 0.5, 1, 1, 1, 0.5, 0};
    for (int k = 0; k < 8; ++k) { EXPECT_NEAR(beta(5, 5, k), want[k], rtol) << k; }
    auto mask = [&] (int i, int j, int k) { return static_cast<int>(w(i, j, k, kWallMask)); };
    EXPECT_EQ(mask(5, 5, 2), kWallUnder);
    EXPECT_NEAR(w(5, 5, 2, kWallD1z), 2.5, rtol);
    EXPECT_NEAR(w(5, 5, 2, kWallD2z), 10.0, rtol);
    EXPECT_EQ(mask(5, 5, 1), 0);
    EXPECT_EQ(mask(5, 5, 6), kWallRoof);
    EXPECT_EQ(mask(3, 5, 0), 0);   // beside the void: no wall
    EXPECT_EQ(fz(5, 5, 0), 0.0);   // the ground and the deck stay open
    EXPECT_EQ(fz(5, 5, 1), 0.0);
    EXPECT_EQ(fz(5, 5, 2), 0.0);   // below the underside, inside the void
    EXPECT_EQ(fz(5, 5, 3), 1.0);
    EXPECT_EQ(fx(4, 5, 0), 0.0);
    EXPECT_NEAR(fx(4, 5, 2), 0.5, rtol);
    EXPECT_EQ(fx(4, 5, 3), 1.0);
    // The underside on a layer face: the fluid cell below the solid one is the wall cell.
    b = lifted(c, "void30.xyz", 65, 30);
    w = gather(b.wall(), c.geom);
    fz = gather(b.beta_face(2), c.geom);
    EXPECT_EQ(mask(5, 5, 2), kWallUnder);
    EXPECT_NEAR(w(5, 5, 2, kWallD1z), 5.0, rtol);
    EXPECT_NEAR(w(5, 5, 2, kWallD2z), 15.0, rtol);
    EXPECT_EQ(fz(5, 5, 3), 1.0);
    // A void inside the ground layer has an underside; a low building in that layer has a roof alone.
    b = lifted(c, "void5.xyz", 65, 5);
    w = gather(b.wall(), c.geom);
    EXPECT_EQ(mask(5, 5, 0), kWallUnder);
    EXPECT_NEAR(w(5, 5, 0, kWallD1z), 2.5, rtol);
    EXPECT_EQ(w(5, 5, 0, kWallD2z), 0.0);   // below it is the ground
    b = lifted(c, "low.xyz", 5, 0);
    w = gather(b.wall(), c.geom);
    EXPECT_EQ(mask(5, 5, 0), kWallRoof);
    // Refusals.
    auto refused = [&] (const std::string& body) {
        return abort_message([&] () { read_building_map(write_text("bad4.xyz", body), "xyz"); });
    };
    EXPECT_NE(refused("0 0 5 5\n").find("zstart must be"), std::string::npos);
    EXPECT_NE(refused("0 0 5 -1\n").find("zstart must be"), std::string::npos);
    EXPECT_NE(refused("0 0 5 1 2\n").find("holds 5 values"), std::string::npos);
    EXPECT_NE(refused("x y h\n").find("not a number"), std::string::npos);
}

// On a stretched mesh d2 reaches the centre of the next layer, not 3/2 of this one.
TEST (Buildings, ReferenceDistanceFollowsTheStretchedLayers)
{
    PCase c;
    c.in.mesh.n_cell = amrex::IntVect(16, 16, 24);
    c.in.mesh.prob_hi = {400, 400, 1000};
    c.in.mesh.is_periodic = {0, 0, 0};
    c.in.mesh.max_grid_size = {8, 8};
    c.in.mesh.dz0 = amrex::Real(2);
    c.in.mesh.stretch_ratio = amrex::Real(1.2);
    c.in.mesh.stretch_ratio_given = true;
    c.in.terrain.smoothing = 0;
    c.make(flat);
    const Gathered J = gather(c.mesh.detJ, c.geom);
    const double dzeta = static_cast<double>(c.geom.CellSize(2));
    double top = 0;
    for (int k = 0; k < 5; ++k) { top += J(5, 5, k) * dzeta; }
    const double dz5 = J(5, 5, 5) * dzeta, dz6 = J(5, 5, 6) * dzeta;
    const gust::Buildings b = make_buildings(c, options(write_erf("stretch.txt", 80, 2.5, 5, box_map(100, 200, 100, 200, top + 0.5 * dz5))));
    const Gathered w = gather(b.wall(), c.geom);
    EXPECT_EQ(static_cast<int>(w(5, 5, 5, gust::kWallMask)), gust::kWallRoof);
    EXPECT_NEAR(w(5, 5, 5, gust::kWallD1z), 0.25 * dz5, 100 * rtol);
    EXPECT_NEAR(w(5, 5, 5, gust::kWallD2z), 0.5 * dz5 + 0.5 * dz6, 100 * rtol);
}

// small_volfrac = 0 with a width that does not sum exactly: an aligned block is still solid.
TEST (Buildings, RoundOffLeavesSolidCellsSolid)
{
    const PCase c = make_case(48, 1000, 16, false, flat);
    const std::string file = write_erf("round.txt", 240, 1000.0 / 480, 1000.0 / 240, [] (double x, double y) {
        const double d = 1000.0 / 48;
        return (x > 10 * d && x < 20 * d && y > 10 * d && y < 20 * d) ? 45.0 : 0.0;
    });
    const Gathered beta = gather(make_buildings(c, options(file, 0)).beta(), c.geom);
    for (int i = 10; i < 20; ++i) {
        for (int j = 10; j < 20; ++j) { EXPECT_EQ(beta(i, j, 2), 1.0) << i << " " << j; }
    }
}

TEST (Buildings, RestartsAreRefusedByTheirHash)
{
    using gust::Buildings;
    const std::string file = write_erf("hash.txt", 80, 2.5, 5, box_map(100, 200, 100, 200, 45));
    const std::string h = Buildings::hash_inputs(options(file));
    EXPECT_NE(h, Buildings::hash_inputs(options(file, 0.01)));
    gust::BuildingsOptions o = options(file);
    o.format = "xyz";
    EXPECT_NE(h, Buildings::hash_inputs(o));
    const std::string dir = fs::path(file_path("chk")).string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::remove_all(dir);
        fs::create_directories(dir);
    }
    amrex::ParallelDescriptor::Barrier();
    EXPECT_NE(abort_message([&] () { Buildings::check_restart(dir, h); }).find("written without buildings"), std::string::npos);
    EXPECT_EQ(abort_message([&] () { Buildings::check_restart(dir, ""); }), "");
}

TEST (Buildings, FacesFollowTheClipAndTheCells)
{
    using namespace gust;
    const PCase c = make_case(16, 400, 8, false, flat);
    // small_volfrac 0.1. Nearly filling cell i = 4 in x (0.95): it clips to solid, and its open low
    // face is then blocked. The thin end at i = 8 (0.05) clips to fluid.
    const HFn clip = box_map(101.25, 201.25, 100, 200, 45);
    Buildings b = make_buildings(c, options(write_erf("clip.txt", 320, 0.625, 1.25, clip), 0.1));
    Gathered beta = gather(b.beta(), c.geom), fx = gather(b.beta_face(0), c.geom);
    EXPECT_EQ(beta(4, 5, 2), 1.0);
    EXPECT_EQ(beta(8, 5, 2), 0.0);
    EXPECT_EQ(fx(4, 5, 2), 1.0);
    // A wall 5 m thick, 10 m long across y = 125 m: both cells below 0.1 (0.2 x 0.2), so fluid; the
    // face between two fluid cells is open, though the wall crosses 0.2 of it.
    b = make_buildings(c, options(write_erf("thin.txt", 320, 0.625, 1.25, [] (double x, double y) {
        return ((x > 110 && x < 115 && y > 120 && y < 130) || (x > 300 && x < 350 && y > 300 && y < 350)) ? 45.0 : 0.0;
    }), 0.1));
    const Gathered fy = gather(b.beta_face(1), c.geom);
    EXPECT_EQ(fy(4, 5, 2), 0.0);
    // A roof on a layer face over partial cells (x fraction 0.6) blocks that part of the face.
    b = make_buildings(c, options(write_erf("roof.txt", 80, 2.5, 5, box_map(110, 190, 100, 200, 40))));
    EXPECT_NEAR(gather(b.beta_face(2), c.geom)(4, 5, 4), 0.6, rtol);
    // A wall on a tile edge blocks its face: only the tile on the high side holds the building.
    b = make_buildings(c, options(write_erf("edge.txt", 80, 2.5, 5, box_map(100, 200, 100, 200, 45))));
    EXPECT_NEAR(gather(b.beta_face(0), c.geom)(4, 5, 4), 0.5, rtol);
}

TEST (Buildings, FilesAndInputsAreChecked)
{
    auto refused = [&] (const std::string& body) {
        return abort_message([&] () { gust::read_building_map(write_text("bad.txt", body)); });
    };
    EXPECT_NE(refused("two\n2\n").find("integers nx and ny"), std::string::npos);
    EXPECT_NE(refused("1 2 0 0 1 1 1\n").find("between 2"), std::string::npos);
    EXPECT_NE(refused("2 2 0 1 0 1 1 1 1\n").find("need 10"), std::string::npos);
    EXPECT_NE(refused("2 2 1 0 0 1 1 1 1 1\n").find("x values must be"), std::string::npos);
    EXPECT_NE(refused("2 2 0 1 1 1 1 1 1 1\n").find("y values must be"), std::string::npos);
    EXPECT_NE(refused("2 2 0 1 0 1 1 1 1 -1\n").find("must be finite and >= 0"), std::string::npos);
    EXPECT_NE(refused("2 2 0 1 0 1 1 x 1 1\n").find("not a number"), std::string::npos);
    gust::BuildingsOptions o = options("no_such_file.txt", 0.5);
    o.format = "csv";
    EXPECT_EQ(gust::check_buildings_options(o).size(), 3u);
    o.type_name = "foo";
    o.type = gust::BuildingsType::none;
    EXPECT_NE(gust::check_buildings_options(o)[0].find("must be none or immersed_forcing"), std::string::npos);
    const PCase c = make_case(16, 400, 8, false, flat);
    const std::string tall = write_erf("tall.txt", 80, 2.5, 5, box_map(100, 200, 100, 200, 235));
    EXPECT_NE(abort_message([&] () { make_buildings(c, options(tall)); }).find("reaches the top layer"), std::string::npos);
    const std::string outside = write_erf("outside.txt", 80, 1002.5, 5, [] (double, double) { return 10.0; });
    EXPECT_NE(abort_message([&] () { make_buildings(c, options(outside)); }).find("no building lies in the domain"), std::string::npos);
}
