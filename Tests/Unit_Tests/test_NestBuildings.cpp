// Buildings inside a nested patch (NestBuildings.H), beyond the plan,
// sub-phase 9.4 (tracker row PR-76; sub-phase 7.4 in the tracker's
// numbering):
//   UT-9.4.1 a building inside a patch has the solid fractions, the face
//            fractions (the fixed rows) and the wall cells of a whole-domain
//            run on the same cells, on sloping ground, in several boxes;
//   UT-9.4.2 a building tile within the gap of a nest face is named, with the
//            face; the gap is 2 cells + max(bc.nest_relax_width cells,
//            terrain.nest_blend_width); a tile on the gap's edge is clear;
//            faces that are not nest, tiles beside a non-periodic patch and tiles
//            with no height are not counted; a periodic direction counts every
//            tile along the face.
#include "NestBuildings.H"
#include "PoissonTestHelpers.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#ifdef _WIN32
#include <process.h>
#define getpid _getpid
#else
#include <unistd.h>
#endif

#include <filesystem>
#include <fstream>
#include <functional>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

using gust::ptest::PCase;

// One tile per building: edges xe, ye, a height in tile (it, jt) only.
gust::BuildingHeightMap one_tile (std::vector<double> xe, std::vector<double> ye, int it, int jt, double h)
{
    gust::BuildingHeightMap m;
    m.nx = static_cast<int>(xe.size()) - 1;
    m.ny = static_cast<int>(ye.size()) - 1;
    m.xe = std::move(xe);
    m.ye = std::move(ye);
    m.h.assign(static_cast<std::size_t>(m.nx * m.ny), 0.0);
    m.zs.assign(m.h.size(), 0.0);
    m.h[static_cast<std::size_t>(it * m.ny + jt)] = h;
    return m;
}

// A patch of 40 x 30 cells of 10 m over [200, 600] x [100, 400] m.
amrex::Geometry patch_geom (bool periodic_y = false)
{
    const amrex::Box dom(amrex::IntVect(0, 0, 0), amrex::IntVect(39, 29, 7));
    const amrex::RealBox rb({AMREX_D_DECL(amrex::Real(200), amrex::Real(100), amrex::Real(0))},
                            {AMREX_D_DECL(amrex::Real(600), amrex::Real(400), amrex::Real(80))});
    return amrex::Geometry(dom, rb, amrex::CoordSys::cartesian, {AMREX_D_DECL(0, periodic_y ? 1 : 0, 0)});
}

gust::BoundaryInputs nest_all (bool periodic_y = false)
{
    gust::BoundaryInputs b;
    b.face[amrex::Orientation(0, amrex::Orientation::low)] = gust::FaceType::nest;
    b.face[amrex::Orientation(0, amrex::Orientation::high)] = gust::FaceType::nest;
    const gust::FaceType fy = periodic_y ? gust::FaceType::periodic : gust::FaceType::nest;
    b.face[amrex::Orientation(1, amrex::Orientation::low)] = fy;
    b.face[amrex::Orientation(1, amrex::Orientation::high)] = fy;
    return b;
}

bool names (const std::vector<std::string>& err, const std::string& face)
{
    for (const auto& e : err) {
        if (e.find("nest face " + face) != std::string::npos) { return true; }
    }
    return false;
}

} // namespace

TEST (NestBuildings, TileOnTheGapEdgeIsClearAndOneCellInsideIsNamed)
{
    const amrex::Geometry g = patch_geom();
    const gust::BoundaryInputs b = nest_all();
    // The gap is 2 cells (20 m). A tile from x = 220 m: clear.
    EXPECT_TRUE(gust::check_nest_buildings(one_tile({180, 220, 240, 260}, {200, 230, 260}, 1, 0, 15), "m", b, g, 0).empty());
    // From x = 210 m: within the gap of bc.xlo, and only of it.
    const auto e = gust::check_nest_buildings(one_tile({180, 210, 240, 260}, {200, 230, 260}, 1, 0, 15), "m.txt", b, g, 0);
    ASSERT_EQ(e.size(), 1u);
    EXPECT_TRUE(names(e, "bc.xlo"));
    EXPECT_NE(e[0].find("buildings.file = m.txt: the building tile x in [210, 240] m, y in [200, 230] m (top 15 m)"),
              std::string::npos) << e[0];
    // A tile across the face, and one outside the patch beyond it (not in the patch: not counted).
    EXPECT_TRUE(names(gust::check_nest_buildings(one_tile({190, 205, 240}, {200, 230}, 0, 0, 15), "m", b, g, 0), "bc.xlo"));
    EXPECT_TRUE(gust::check_nest_buildings(one_tile({150, 200, 240}, {200, 230}, 0, 0, 15), "m", b, g, 0).empty());
}

TEST (NestBuildings, EachFaceAndCorner)
{
    const amrex::Geometry g = patch_geom();
    const gust::BoundaryInputs b = nest_all();
    // High x: 600 - 20 = 580 m is the gap's edge.
    EXPECT_TRUE(gust::check_nest_buildings(one_tile({550, 580, 700}, {200, 230}, 0, 0, 5), "m", b, g, 0).empty());
    EXPECT_TRUE(names(gust::check_nest_buildings(one_tile({550, 581, 700}, {200, 230}, 0, 0, 5), "m", b, g, 0), "bc.xhi"));
    // Low and high y.
    EXPECT_TRUE(names(gust::check_nest_buildings(one_tile({300, 330}, {90, 119, 150}, 0, 1, 5), "m", b, g, 0), "bc.ylo"));
    EXPECT_TRUE(gust::check_nest_buildings(one_tile({300, 330}, {90, 120, 150}, 0, 1, 5), "m", b, g, 0).empty());
    EXPECT_TRUE(names(gust::check_nest_buildings(one_tile({300, 330}, {381, 390}, 0, 0, 5), "m", b, g, 0), "bc.yhi"));
    // A corner tile is named for both faces.
    const auto c = gust::check_nest_buildings(one_tile({205, 215}, {105, 115}, 0, 0, 5), "m", b, g, 0);
    EXPECT_EQ(c.size(), 2u);
    EXPECT_TRUE(names(c, "bc.xlo"));
    EXPECT_TRUE(names(c, "bc.ylo"));
    // No height: no building.
    EXPECT_TRUE(gust::check_nest_buildings(one_tile({205, 215}, {105, 115}, 0, 0, 0), "m", b, g, 0).empty());
}

TEST (NestBuildings, ZoneAndBlendWidenTheGap)
{
    const amrex::Geometry g = patch_geom();
    gust::BoundaryInputs b = nest_all();
    const auto tile = one_tile({235, 260}, {200, 230}, 0, 0, 12);   // 35 m from bc.xlo
    EXPECT_TRUE(gust::check_nest_buildings(tile, "m", b, g, 0).empty());
    // The gap is 2 cells beyond the larger of the zone and the blend.
    b.nest_relax_width = 1;   // 20 + 10 m
    EXPECT_TRUE(gust::check_nest_buildings(tile, "m", b, g, 0).empty());
    b.nest_relax_width = 2;   // 20 + 20 m
    const auto ez = gust::check_nest_buildings(tile, "m", b, g, 0);
    ASSERT_EQ(ez.size(), 1u);
    EXPECT_NE(ez[0].find("within 40 m"), std::string::npos) << ez[0];
    EXPECT_NE(ez[0].find("bc.nest_relax_width = 2"), std::string::npos) << ez[0];
    b.nest_relax_width = 0;
    EXPECT_TRUE(gust::check_nest_buildings(tile, "m", b, g, amrex::Real(15)).empty());
    const auto eb = gust::check_nest_buildings(tile, "m", b, g, amrex::Real(16));
    ASSERT_EQ(eb.size(), 1u);
    EXPECT_NE(eb[0].find("terrain.nest_blend_width = 16"), std::string::npos) << eb[0];
    // The zone and the blend do not add: the larger counts.
    b.nest_relax_width = 1;
    EXPECT_TRUE(gust::check_nest_buildings(tile, "m", b, g, amrex::Real(15)).empty());
}

TEST (NestBuildings, OnlyNestFacesAndTilesAlongThePatch)
{
    const amrex::Geometry g = patch_geom();
    gust::BoundaryInputs b = nest_all();
    b.face[amrex::Orientation(0, amrex::Orientation::low)] = gust::FaceType::inflow;
    EXPECT_TRUE(gust::check_nest_buildings(one_tile({205, 215}, {200, 230}, 0, 0, 5), "m", b, g, 0).empty());
    b = nest_all();
    // Near bc.xlo in x but beside the patch in y (y < 100 m): not in the patch.
    EXPECT_TRUE(gust::check_nest_buildings(one_tile({205, 215}, {40, 100}, 0, 0, 5), "m", b, g, 0).empty());
    // y periodic: an image lies in the patch, so it counts.
    const amrex::Geometry gp = patch_geom(true);
    EXPECT_TRUE(names(gust::check_nest_buildings(one_tile({205, 215}, {40, 90}, 0, 0, 5), "m", nest_all(true), gp, 0), "bc.xlo"));
}

TEST (NestBuildings, MessagesAreCappedWithACount)
{
    const amrex::Geometry g = patch_geom();
    gust::BuildingHeightMap m = one_tile({205, 215}, {100, 110, 120, 130, 140, 150, 160}, 0, 0, 5);
    for (double& h : m.h) { h = 5; }
    const auto e = gust::check_nest_buildings(m, "m", nest_all(), g, 0, 4);
    ASSERT_EQ(e.size(), 5u);   // 6 tiles near bc.xlo, the first two also near bc.ylo: 8, 4 named
    EXPECT_NE(e.back().find("4 more pairs of a building tile and a nest face"), std::string::npos) << e.back();
}

// UT-9.4.1: the patch builds its buildings as the whole domain does on the same cells.
namespace {

struct Plain
{
    amrex::Box box;
    int ncomp{0};
    std::vector<double> v;
    [[nodiscard]] double at (int i, int j, int k, int n) const
    {
        const auto l = box.length();
        return v[static_cast<std::size_t>((i - box.smallEnd(0)) + l[0] * ((j - box.smallEnd(1)) + static_cast<long>(l[1]) * ((k - box.smallEnd(2)) + static_cast<long>(l[2]) * n)))];
    }
};

Plain gather_plain (const amrex::MultiFab& f, const amrex::Geometry& geom)
{
    Plain g;
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

std::string map_file ()
{
    static const long id = [] () {
        long pid = static_cast<long>(::getpid());
        amrex::ParallelDescriptor::Bcast(&pid, 1, amrex::ParallelDescriptor::IOProcessorNumber());
        return pid;
    }();
    const fs::path dir = fs::temp_directory_path() / ("gust_nest_buildings_" + std::to_string(id));
    const std::string file = (dir / "map.txt").string();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        fs::create_directories(dir);
        std::ofstream os(file);
        os.precision(17);
        // 5 m tiles over [100, 300] x [100, 300] m: an L-shaped block (two tops) and a lifted deck, walls between cell faces.
        os << "# x y height zstart\n";
        for (int i = 0; i < 40; ++i) {
            for (int j = 0; j < 40; ++j) {
                const double x = 102.5 + 5 * i, y = 102.5 + 5 * j;
                double h = 0, zs = 0;
                if (x > 160 && x < 215 && y > 150 && y < 245) { h = 27; }
                if (x > 215 && x < 240 && y > 150 && y < 185) { h = 43; }
                if (x > 170 && x < 200 && y > 255 && y < 280) { h = 18; zs = 6; }
                os << x << " " << y << " " << h << " " << zs << "\n";
            }
        }
    }
    amrex::ParallelDescriptor::Barrier();
    return file;
}

} // namespace

TEST (NestBuildings, PatchHasTheWholeDomainsFractionsAndWalls)
{
    const auto ground = [] (double x, double y) { return 0.04 * x + 0.015 * y + 5e-5 * (x - 200) * (y - 180); };
    // The whole domain: 32 x 32 cells of 12.5 m over [0, 400]^2. The patch: 16 x 17 of its cells over
    // [125, 325] x [112.5, 325] m, in other boxes; the buildings lie 2.8 or more cells inside it.
    PCase whole;
    whole.in.mesh.n_cell = amrex::IntVect(32, 32, 24);
    whole.in.mesh.prob_hi = {400, 400, 240};
    whole.in.mesh.max_grid_size = {16, 16};
    whole.in.terrain.smoothing = 0;
    whole.make(ground);
    PCase patch;
    patch.in.mesh.n_cell = amrex::IntVect(16, 17, 24);
    patch.in.mesh.prob_lo = {125, 112.5, 0};
    patch.in.mesh.prob_hi = {325, 325, 240};
    patch.in.mesh.max_grid_size = {5, 6};
    patch.in.terrain.smoothing = 0;
    patch.make(ground);

    gust::BuildingsOptions o;
    o.type = gust::BuildingsType::immersed_forcing;
    o.type_name = "immersed_forcing";
    o.file = map_file();
    o.format = "xyz";
    // The check passes for the patch with all four faces nest.
    EXPECT_TRUE(gust::check_nest_buildings(gust::read_building_map(o.file, o.format), o.file, nest_all(), patch.geom, 0).empty());

    gust::Buildings bw, bp;
    bw.build(o, whole.geom, whole.z, whole.ba, whole.dm);
    bp.build(o, patch.geom, patch.z, patch.ba, patch.dm);
    const int oi = 10, oj = 9;   // the whole domain's cell under patch cell (0, 0)
    struct Field { const amrex::MultiFab* w; const amrex::MultiFab* p; const char* name; };
    const std::vector<Field> fields = {{&bw.beta(), &bp.beta(), "beta"},
                                       {&bw.beta_face(0), &bp.beta_face(0), "beta_x"},
                                       {&bw.beta_face(1), &bp.beta_face(1), "beta_y"},
                                       {&bw.beta_face(2), &bp.beta_face(2), "beta_z"},
                                       {&bw.wall(), &bp.wall(), "wall"}};
    long solid = 0, walls = 0;
    for (const Field& f : fields) {
        const Plain w = gather_plain(*f.w, whole.geom), p = gather_plain(*f.p, patch.geom);
        for (int n = 0; n < p.ncomp; ++n) {
            for (int k = p.box.smallEnd(2); k <= p.box.bigEnd(2); ++k) {
                for (int j = p.box.smallEnd(1); j <= p.box.bigEnd(1); ++j) {
                    for (int i = p.box.smallEnd(0); i <= p.box.bigEnd(0); ++i) {
                        const double a = w.at(i + oi, j + oj, k, n), b = p.at(i, j, k, n);
                        EXPECT_EQ(a, b) << f.name << " comp " << n << " at patch cell " << i << " " << j << " " << k;
                        if (std::string(f.name) == "beta" && b > 0) { ++solid; }
                        if (std::string(f.name) == "wall" && n == gust::kWallMask && b != 0) { ++walls; }
                    }
                }
            }
        }
    }
    // Not vacuous: the patch holds solid, partly solid and wall cells.
    EXPECT_GT(solid, 100);
    EXPECT_GT(walls, 50);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::error_code ec;
        fs::remove_all(fs::path(o.file).parent_path(), ec);
    }
}
