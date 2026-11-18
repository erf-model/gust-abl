// The nest inputs (Nests.H): nest.names and each nest's corners, the
// placement checks, and the settings a nest derives from the outer deck.
// The joint run and the abort messages are exercised by the deck tests in
// Tests/Canonical_RANS/Nests_One_Process.
#include "Forcing.H"
#include "Mesh.H"
#include "NestParent.H"
#include "Nests.H"
#include "State.H"
#include "SurfaceLayer.H"
#include "TerrainSurface.H"

#include <AMReX_MultiFab.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>

namespace {

bool any_mentions (const std::vector<std::string>& msgs, const std::string& s)
{
    return std::any_of(msgs.begin(), msgs.end(),
                       [&] (const std::string& m) { return m.find(s) != std::string::npos; });
}

// Adds keys the way an inputs file does (test_Inputs.cpp).
void add_inputs_text (const std::string& text)
{
    const std::string fname = "ut_nests_tmp";
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::ofstream os(fname);
        os << text;
    }
    amrex::ParallelDescriptor::Barrier();
    amrex::ParmParse::addfile(fname);
    amrex::ParallelDescriptor::Barrier();
    if (amrex::ParallelDescriptor::IOProcessor()) { std::remove(fname.c_str()); }
}

// Removes the keys a test added.
struct ScopedKeys
{
    std::vector<std::pair<std::string, std::string>> keys; // (prefix, name)
    ~ScopedKeys () {
        for (const auto& [prefix, name] : keys) {
            amrex::ParmParse pp(prefix);
            while (pp.remove(name) > 0) {}
        }
    }
};

// 16 x 12 x 8 cells of 100 m in x and y, from (0, 0).
gust::MeshInputs outer_mesh ()
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(16, 12, 8);
    m.prob_hi = {1600, 1200, 800};
    m.max_grid_size = {13, 9};
    m.dz0 = 5;
    return m;
}

gust::NestSpec nest (const std::string& name, amrex::Real xlo, amrex::Real ylo, amrex::Real xhi, amrex::Real yhi)
{
    gust::NestSpec s;
    s.name = name;
    s.lo = {xlo, ylo};
    s.hi = {xhi, yhi};
    return s;
}

std::vector<std::string> check (const std::vector<gust::NestSpec>& nests, const gust::MeshInputs& m)
{
    gust::NestOptions opt;
    opt.nests = nests;
    return gust::check_nests(opt, m);
}

// The largest |ground - hm / (1 + (r / L)^2)| over the valid ground nodes, r
// the distance (double) from (xc, yc).
double agnesi_ground_error (const gust::MeshInputs& mi, const gust::TerrainInputs& t, double xc, double yc)
{
    const amrex::Geometry geom = gust::make_geometry(mi);
    const amrex::BoxArray ba = gust::make_box_array(mi);
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab z(amrex::convert(ba, amrex::IntVect(1)), dm, 1, gust::ngrow_nodal);
    z.setVal(0);
    gust::fill_terrain_surface(z, geom, t);
    amrex::MultiFab e(z.boxArray(), dm, 1, 0);
    e.setVal(0);
    const auto dx = geom.CellSizeArray();
    const auto plo = geom.ProbLoArray();
    const double hm = t.hill_height, L = t.hill_half_width;
    for (amrex::MFIter mfi(e); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox();
        bx.setBig(2, 0);
        const auto ea = e.array(mfi);
        const auto za = z.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double x = static_cast<double>(plo[0]) + i * static_cast<double>(dx[0]) - xc;
            const double y = static_cast<double>(plo[1]) + j * static_cast<double>(dx[1]) - yc;
            const double h = hm / (1 + (x * x + y * y) / (L * L));
            ea(i, j, k) = static_cast<amrex::Real>(std::abs(static_cast<double>(za(i, j, k)) - h));
        });
    }
    return static_cast<double>(e.norm0(0, 0));
}

} // namespace

TEST (Nests, NoNamesIsNoNest)
{
    const auto opt = gust::read_nest_options();
    EXPECT_FALSE(opt.active());
    EXPECT_TRUE(opt.errors.empty());
}

TEST (Nests, ReadsNamesAndCorners)
{
    ScopedKeys guard{{{"nest", "names"}, {"nest.a", "lo"}, {"nest.a", "hi"}, {"nest.b2", "lo"}, {"nest.b2", "hi"}}};
    add_inputs_text("nest.names = a b2\n"
                    "nest.a.lo = 400 300\nnest.a.hi = 1200 900\n"
                    "nest.b2.lo = 200 100\nnest.b2.hi = 400 300\n");
    const auto opt = gust::read_nest_options();
    ASSERT_TRUE(opt.errors.empty()) << opt.errors.front();
    ASSERT_EQ(opt.nests.size(), 2u);
    EXPECT_EQ(opt.nests[0].name, "a");
    EXPECT_EQ(opt.nests[0].lo[0], 400);
    EXPECT_EQ(opt.nests[0].lo[1], 300);
    EXPECT_EQ(opt.nests[0].hi[0], 1200);
    EXPECT_EQ(opt.nests[0].hi[1], 900);
    EXPECT_EQ(opt.nests[1].name, "b2");
    EXPECT_EQ(opt.nests[1].hi[1], 300);
}

TEST (Nests, ReadErrorsNameTheKey)
{
    ScopedKeys guard{{{"nest", "names"}, {"nest.a", "lo"}, {"nest.a", "hi"}, {"nest.c", "lo"}, {"nest.c", "hi"}}};
    add_inputs_text("nest.names = a 9x names a c\n"
                    "nest.a.lo = 400 300\n"
                    "nest.c.lo = 400 300 0\nnest.c.hi = 1200 900\nnest.c.hi = 1200 900\n");
    const auto opt = gust::read_nest_options();
    EXPECT_TRUE(any_mentions(opt.errors, "nest.a.hi is missing"));
    EXPECT_TRUE(any_mentions(opt.errors, "'9x' is not a nest name"));
    EXPECT_TRUE(any_mentions(opt.errors, "'names' is not a nest name"));
    EXPECT_TRUE(any_mentions(opt.errors, "'a' is named more than once"));
    EXPECT_TRUE(any_mentions(opt.errors, "nest.c.lo takes two values"));
    EXPECT_TRUE(any_mentions(opt.errors, "nest.c.hi is given more than once"));
    EXPECT_TRUE(opt.nests.empty());
}

TEST (Nests, PlacementOnTheOuterCellFaces)
{
    const auto m = outer_mesh();
    EXPECT_TRUE(check({nest("a", 400, 300, 1200, 900)}, m).empty());
    const auto off = check({nest("a", 450, 300, 1200, 900)}, m);
    ASSERT_EQ(off.size(), 1u);
    EXPECT_TRUE(any_mentions(off, "nest.a.lo and nest.a.hi x must lie on the outer cell faces"));
    // Within the 1e-4-cell tolerance of check_nest_placement.
    EXPECT_TRUE(check({nest("a", 400.005, 300, 1200, 900)}, m).empty());
    EXPECT_TRUE(any_mentions(check({nest("a", 400, 300, 1200, 300)}, m), "nest.a.hi y must be above nest.a.lo y"));
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    EXPECT_TRUE(any_mentions(check({nest("a", 400, nan, 1200, 900)}, m), "y must be finite"));
}

TEST (Nests, MarginOfTwoOuterCells)
{
    auto m = outer_mesh();
    // 2 cells from each face: allowed; 1 cell: refused.
    EXPECT_TRUE(check({nest("a", 200, 200, 1400, 1000)}, m).empty());
    EXPECT_TRUE(any_mentions(check({nest("a", 100, 200, 1400, 1000)}, m), "at least 2 outer cells inside the outer x faces"));
    EXPECT_TRUE(any_mentions(check({nest("a", 200, 200, 1400, 1100)}, m), "at least 2 outer cells inside the outer y faces"));
    // Across a periodic direction the nest may reach the outer faces.
    m.is_periodic[0] = 1;
    EXPECT_TRUE(check({nest("a", 0, 200, 1600, 1000)}, m).empty());
    EXPECT_FALSE(check({nest("a", 0, 100, 1600, 1000)}, m).empty());
    // Periodic in both: the whole outer domain leaves the nest no nest face.
    m.is_periodic[1] = 1;
    EXPECT_TRUE(check({nest("a", 0, 0, 1600, 1100)}, m).empty());
    const auto all = check({nest("a", 0, 0, 1600, 1200)}, m);
    ASSERT_EQ(all.size(), 1u);
    EXPECT_TRUE(any_mentions(all, "nest.a covers the whole outer domain in x and y"));
}

TEST (Nests, NestsMayTouchNotOverlap)
{
    const auto m = outer_mesh();
    EXPECT_TRUE(check({nest("a", 200, 200, 800, 600), nest("b", 800, 200, 1200, 600)}, m).empty());
    // Touching at one corner only.
    EXPECT_TRUE(check({nest("a", 200, 200, 800, 600), nest("b", 800, 600, 1200, 1000)}, m).empty());
    const auto e = check({nest("a", 200, 200, 800, 600), nest("b", 700, 500, 1200, 1000)}, m);
    ASSERT_EQ(e.size(), 1u);
    EXPECT_TRUE(any_mentions(e, "nest.b overlaps nest.a"));
}

TEST (Nests, MeshIsTwiceTheOuterCellsInXY)
{
    auto m = outer_mesh();
    m.is_periodic[1] = 1;
    m.hierarchy_levels = 3;
    const auto n = gust::nest_mesh(nest("a", 400, 300, 1200, 900), m);
    EXPECT_EQ(n.n_cell, amrex::IntVect(16, 12, 8));
    EXPECT_EQ(n.prob_lo[0], 400);
    EXPECT_EQ(n.prob_lo[1], 300);
    EXPECT_EQ(n.prob_lo[2], 0);
    EXPECT_EQ(n.prob_hi[0], 1200);
    EXPECT_EQ(n.prob_hi[1], 900);
    EXPECT_EQ(n.prob_hi[2], 800);
    EXPECT_EQ(n.is_periodic[0], 0);
    EXPECT_EQ(n.is_periodic[1], 0);
    EXPECT_EQ(n.hierarchy_levels, 1);
    EXPECT_EQ(n.max_grid_size[0], 13);
    EXPECT_EQ(n.max_grid_size[1], 9);
    EXPECT_EQ(n.dz0, 5);
    // An outer domain that does not start at 0.
    m.prob_lo = {-800, 1000, 0};
    m.prob_hi = {800, 2200, 800};
    const auto k = gust::nest_mesh(nest("k", -600, 1300, 0, 1500), m);
    EXPECT_EQ(k.n_cell[0], 12);
    EXPECT_EQ(k.n_cell[1], 4);
    // A corner within the tolerance is placed on the outer cell face.
    const auto s = gust::nest_mesh(nest("s", -599.995, 1300, 0.004, 1500), m);
    EXPECT_EQ(s.prob_lo[0], -600);
    EXPECT_EQ(s.prob_hi[0], 0);
    // The whole width of a periodic outer direction: periodic there.
    const auto w = gust::nest_mesh(nest("w", -600, 1000, 0, 2200), m);
    EXPECT_EQ(w.n_cell[1], 24);
    EXPECT_EQ(w.is_periodic[0], 0);
    EXPECT_EQ(w.is_periodic[1], 1);
    EXPECT_EQ(w.prob_lo[1], 1000);
    EXPECT_EQ(w.prob_hi[1], 2200);
    // On the high face of a periodic direction, not from its low face: that face exactly (in double,
    // -800.1 + 12 * (1600.4 / 12) is 1.1e-13 above 800.3).
    m.prob_lo = {-800, -800.1, 0};
    m.prob_hi = {800, 800.3, 800};
    const auto h = gust::nest_mesh(nest("h", -600, 0.1, 0, 800.3), m);
    EXPECT_EQ(h.n_cell[1], 12);
    EXPECT_EQ(h.is_periodic[1], 0);
    EXPECT_EQ(h.prob_hi[1], m.prob_hi[1]);
}

TEST (Nests, InputsFromTheOuterDeck)
{
    gust::Inputs in;
    in.mesh = outer_mesh();
    in.io.plot_file = "plt";
    in.io.chk_file = "chk";
    in.io.turb_plot_file.clear();
    in.io.restart = "chk00010";
    in.terrain.given = {"terrain.type"};
    gust::apply_nest_inputs(nest("hill", 400, 300, 1200, 900), in);
    EXPECT_EQ(in.mesh.n_cell, amrex::IntVect(16, 12, 8));
    EXPECT_EQ(in.io.plot_file, "plt_hill");
    EXPECT_EQ(in.io.chk_file, "chk_hill");
    EXPECT_EQ(in.io.turb_plot_file, "");
    EXPECT_EQ(in.io.restart, "");
    EXPECT_TRUE(in.terrain.given.empty());
}

TEST (Nests, BoundaryInputsFromTheOuterDeck)
{
    gust::BoundaryInputs bin;
    bin.face[amrex::Orientation(0, amrex::Orientation::low)] = gust::FaceType::inflow;
    bin.face[amrex::Orientation(0, amrex::Orientation::high)] = gust::FaceType::outflow;
    const auto ground = bin.face[amrex::Orientation(2, amrex::Orientation::low)];
    const auto top = bin.face[amrex::Orientation(2, amrex::Orientation::high)];
    bin.sponge_width = 300;
    bin.sponge_strength = 0.01;
    bin.nest_relax_width = 3;
    bin.nest_relax_factor = 0.5;
    bin.given = {"bc.xlo"};
    gust::apply_nest_boundary_inputs("plt_chk_final", amrex::Array<int, 3>{0, 0, 0}, false, bin);
    for (int d = 0; d < 2; ++d) {
        EXPECT_EQ(bin.face[amrex::Orientation(d, amrex::Orientation::low)], gust::FaceType::nest);
        EXPECT_EQ(bin.face[amrex::Orientation(d, amrex::Orientation::high)], gust::FaceType::nest);
    }
    EXPECT_EQ(bin.face[amrex::Orientation(2, amrex::Orientation::low)], ground);
    EXPECT_EQ(bin.face[amrex::Orientation(2, amrex::Orientation::high)], top);
    EXPECT_EQ(bin.nest_parent, "plt_chk_final");
    EXPECT_EQ(bin.sponge_width, 0);
    EXPECT_EQ(bin.sponge_strength, 0);
    const std::array<bool, 5> fields{false, false, true, true, false};
    EXPECT_EQ(bin.rayleigh_fields, fields);
    EXPECT_FALSE(bin.rayleigh_from_outer);
    EXPECT_EQ(bin.nest_relax_width, 0);
    EXPECT_EQ(bin.nest_relax_factor, 1);
    EXPECT_TRUE(bin.given.empty());
    // A nest periodic in y keeps periodic faces there.
    gust::BoundaryInputs per;
    gust::apply_nest_boundary_inputs("c", amrex::Array<int, 3>{0, 1, 0}, false, per);
    EXPECT_EQ(per.face[amrex::Orientation(0, amrex::Orientation::low)], gust::FaceType::nest);
    EXPECT_EQ(per.face[amrex::Orientation(1, amrex::Orientation::low)], gust::FaceType::periodic);
    EXPECT_EQ(per.face[amrex::Orientation(1, amrex::Orientation::high)], gust::FaceType::periodic);
    // nest.damping = outer: the outer deck's Rayleigh fields as they are, and the nest-face check lets them pass.
    gust::BoundaryInputs od;
    od.rayleigh_fields = {true, true, true, false, true};
    gust::apply_nest_boundary_inputs("c", amrex::Array<int, 3>{0, 0, 0}, true, od);
    const std::array<bool, 5> kept{true, true, true, false, true};
    EXPECT_EQ(od.rayleigh_fields, kept);
    EXPECT_TRUE(od.rayleigh_from_outer);
    EXPECT_EQ(od.sponge_width, 0);
}

// A nest off the centre of the outer domain keeps the outer domain's hill
// (TerrainInputs::hill_centre, which App.cpp sets from the parent checkpoint).
TEST (Nests, AgnesiHillStaysAtTheOuterCentre)
{
    gust::MeshInputs outer = outer_mesh();   // centre (800, 600)
    gust::TerrainInputs t;
    t.type = "agnesi";
    t.hill_shape = "round";
    t.hill_height = 60;
    t.hill_half_width = 250;
    // The nest over outer cells 10 to 14 in x and 3 to 7 in y, centre (1200, 500).
    gust::MeshInputs m = gust::nest_mesh(nest("e", 1000, 300, 1400, 700), outer);
    const double tol = 1e-4;   // m: single precision of 60 m heights
    t.hill_centre = {800, 600};
    EXPECT_LT(agnesi_ground_error(m, t, 800, 600), tol);
    // Without it the hill is the nest's own: 60 m at (1200, 500).
    t.hill_centre = {std::numeric_limits<amrex::Real>::quiet_NaN(), std::numeric_limits<amrex::Real>::quiet_NaN()};
    EXPECT_LT(agnesi_ground_error(m, t, 1200, 500), tol);
    EXPECT_GT(agnesi_ground_error(m, t, 800, 600), 30);
}

TEST (Nests, OutputNames)
{
    EXPECT_EQ(gust::nest_output_name("plt", "n1"), "plt_n1");
    EXPECT_EQ(gust::nest_output_name("", "n1"), "");
}

// The outer run's forcing reaches each nest: the geostrophic wind, the two
// horizontal pressure-gradient components and the surface heat flux; the
// vertical component and the nest's other settings stay as the deck gives.
TEST (Nests, OuterForcingReachesNest)
{
    gust::ForcingOptions outer_f;
    gust::SurfaceLayerOptions outer_s;
    outer_f.coriolis.ug = amrex::Real(7.5);
    outer_f.coriolis.vg = amrex::Real(-2.25);
    outer_f.pressure_gradient = {amrex::Real(1.5e-4), amrex::Real(-3.0e-4), amrex::Real(9.0)};
    outer_s.heat_flux = amrex::Real(0.125);
    gust::OuterOutcome out;
    EXPECT_FALSE(out.forcing_set);

    gust::ForcingOptions nest_f;
    gust::SurfaceLayerOptions nest_s;
    nest_f.coriolis.ug = amrex::Real(1);
    nest_f.coriolis.vg = amrex::Real(2);
    nest_f.pressure_gradient = {amrex::Real(3), amrex::Real(4), amrex::Real(5)};
    nest_s.heat_flux = amrex::Real(6);
    // Nothing kept: the nest keeps its deck's values.
    gust::take_outer_forcing(out, nest_f, nest_s);
    EXPECT_EQ(nest_f.coriolis.ug, amrex::Real(1));
    EXPECT_EQ(nest_f.coriolis.vg, amrex::Real(2));
    EXPECT_EQ(nest_f.pressure_gradient[0], amrex::Real(3));
    EXPECT_EQ(nest_f.pressure_gradient[1], amrex::Real(4));
    EXPECT_EQ(nest_s.heat_flux, amrex::Real(6));

    gust::keep_outer_forcing(outer_f, outer_s, out);
    EXPECT_TRUE(out.forcing_set);
    gust::take_outer_forcing(out, nest_f, nest_s);
    EXPECT_EQ(nest_f.coriolis.ug, amrex::Real(7.5));
    EXPECT_EQ(nest_f.coriolis.vg, amrex::Real(-2.25));
    EXPECT_EQ(nest_f.pressure_gradient[0], amrex::Real(1.5e-4));
    EXPECT_EQ(nest_f.pressure_gradient[1], amrex::Real(-3.0e-4));
    EXPECT_EQ(nest_f.pressure_gradient[2], amrex::Real(5));
    EXPECT_EQ(nest_s.heat_flux, amrex::Real(0.125));
}
