// The hybrid mesh's envelope (HybridTerrain.H, terrain.representation = hybrid):
//   - on a slope a x + b y the envelope lies a p dx + b q dy below the ground away from the faces (the
//     minimum of the square is its low corner, and the mean of a plane over the square is its centre);
//   - in a valley the envelope is the mean of the minimum, not the minimum;
//   - a one-node spike is removed, and the envelope is at most the ground everywhere;
//   - across a periodic face the square wraps: a shifted ground gives the shifted envelope, and a square
//     wider than the line gives the line's minimum;
//   - the half-width in nodes is floor(width / dx);
//   - make_hybrid_envelope gives the host envelope at every node, ghosts included, whatever the boxes and
//     the rank count (the suite runs on 1, 2 and 4 ranks).
#include "HybridTerrain.H"
#include "LaminarTestHelpers.H"
#include "MassFlux.H"
#include "SurfaceRoughness.H"

#include <AMReX_ParallelDescriptor.H>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using gust::ptest::PCase;

std::vector<double> grid (int nx, int ny, double (*f)(int, int))
{
    std::vector<double> h(static_cast<std::size_t>(nx) * ny);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) { h[static_cast<std::size_t>(j) * nx + i] = f(i, j); }
    }
    return h;
}

}   // namespace

TEST (HybridTerrain, EnvelopeOfAPlaneIsTheShiftedPlane)
{
    const int nx = 40, ny = 30, px = 3, py = 2;
    const auto h = grid(nx, ny, [] (int i, int j) { return 0.7 * i - 0.4 * j + 5.0; });
    const auto e = gust::lower_envelope(h, nx, ny, px, py, false, false);
    double worst = 0, above = 0;
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            const std::size_t n = static_cast<std::size_t>(j) * nx + i;
            above = std::max(above, e[n] - h[n]);
            // The squares of every node of the square stay inside the domain.
            if (i >= 2 * px && i < nx - 2 * px && j >= 2 * py && j < ny - 2 * py) {
                worst = std::max(worst, std::abs(e[n] - (h[n] - 0.7 * px - 0.4 * py)));
            }
        }
    }
    EXPECT_LT(worst, 1e-12);
    EXPECT_LE(above, 0.0);
    // No half-width: the ground itself.
    EXPECT_EQ(gust::lower_envelope(h, nx, ny, 0, 0, false, false), h);
}

// A V-shaped valley a |i - c|: the minimum is a max(0, |i - c| - p), and the mean of that over the window of
// node c + p is a p (p + 1) / (2 (2 p + 1)) (the minimum alone would give 0 there).
TEST (HybridTerrain, EnvelopeOfAValleyIsTheMeanOfTheMinimum)
{
    const int nx = 41, c = 20, p = 4;
    std::vector<double> h(nx);
    for (int i = 0; i < nx; ++i) { h[static_cast<std::size_t>(i)] = 1.5 * std::abs(i - c); }
    const auto e = gust::lower_envelope(h, nx, 1, p, 0, false, false);
    EXPECT_NEAR(e[c], 0.0, 1e-14);
    EXPECT_NEAR(e[c + p], 1.5 * p * (p + 1) / (2.0 * (2 * p + 1)), 1e-14);
    EXPECT_NEAR(e[c - p], 1.5 * p * (p + 1) / (2.0 * (2 * p + 1)), 1e-14);
    // Away from the valley the slope's shift: a p.
    EXPECT_NEAR(e[c + 3 * p], h[c + 3 * p] - 1.5 * p, 1e-13);
}

TEST (HybridTerrain, EnvelopeRemovesASpikeAndStaysBelowTheGround)
{
    const int nx = 21, ny = 1;
    std::vector<double> h(nx, 3.0);
    h[10] = 50;
    const auto e = gust::lower_envelope(h, nx, ny, 1, 0, false, false);
    for (int i = 0; i < nx; ++i) { EXPECT_DOUBLE_EQ(e[static_cast<std::size_t>(i)], 3.0) << i; }
    // A rough ground: never above it.
    const auto r = grid(31, 17, [] (int i, int j) { return std::sin(1.3 * i) * std::cos(0.7 * j) * 20 + 0.1 * i * j; });
    for (const bool per : {false, true}) {
        const auto er = gust::lower_envelope(r, 31, 17, 4, 3, per, per);
        for (std::size_t n = 0; n < r.size(); ++n) { EXPECT_LE(er[n], r[n]) << n; }
    }
}

TEST (HybridTerrain, PeriodicEnvelopeWraps)
{
    const int nx = 32, ny = 4, shift = 11;
    auto f = [] (int i) { return 10 * std::sin(2 * 3.14159265358979323846 * i / 32.0) + 3 * std::cos(6 * 3.14159265358979323846 * i / 32.0); };
    std::vector<double> a(static_cast<std::size_t>(nx) * ny), b(a.size());
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            a[static_cast<std::size_t>(j) * nx + i] = f(i) + j;
            b[static_cast<std::size_t>(j) * nx + i] = f((i + shift) % nx) + j;
        }
    }
    const auto ea = gust::lower_envelope(a, nx, ny, 5, 1, true, true);
    const auto eb = gust::lower_envelope(b, nx, ny, 5, 1, true, true);
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            EXPECT_NEAR(eb[static_cast<std::size_t>(j) * nx + i], ea[static_cast<std::size_t>(j) * nx + (i + shift) % nx], 1e-12) << i << " " << j;
        }
    }
    // A square wider than the line: every node's minimum is the line's, and so is the mean.
    const std::vector<double> line{4.0, 1.0, 7.0, 2.0};
    const auto el = gust::lower_envelope(line, 4, 1, 5, 0, true, false);
    for (const double v : el) { EXPECT_DOUBLE_EQ(v, 1.0); }
}

TEST (HybridTerrain, HalfWidthInNodes)
{
    const PCase c = gust::lamtest::flat_case(64, 8, 4, 1024, 320, 400, true, true, 32);   // dx 16, dy 40
    EXPECT_EQ(gust::hybrid_envelope_nodes(128, c.geom), amrex::IntVect(8, 3, 0));
    EXPECT_EQ(gust::hybrid_envelope_nodes(127.9, c.geom), amrex::IntVect(7, 3, 0));
    EXPECT_EQ(gust::hybrid_envelope_nodes(0, c.geom), amrex::IntVect(0, 0, 0));
}

// Without the keys, terrain.flat_height is the first layer's thickness and terrain.hybrid_half_width is
// 2 cells (the larger lateral cell size; a direction of one cell does not count). Given keys and other
// representations are left alone.
TEST (HybridTerrain, DefaultsAreOneFirstLayerAndTwoCells)
{
    auto make = [] () {
        gust::Inputs in;
        in.mesh.n_cell = amrex::IntVect(64, 8, 20);
        in.mesh.prob_lo = {-512, 0, 0};
        in.mesh.prob_hi = {512, 320, 1000};   // dx 16, dy 40
        in.mesh.is_periodic = {1, 1, 0};
        in.mesh.dz0 = 2;
        in.mesh.stretch_ratio = 1.1;
        in.mesh.stretch_ratio_given = true;
        in.terrain.representation = "hybrid";
        in.terrain.type = "agnesi";
        return in;
    };
    {
        auto in = make();
        EXPECT_EQ(gust::apply_hybrid_defaults(in).size(), 2u);
        EXPECT_NEAR(static_cast<double>(in.terrain.flat_height), 2.0, 1e-12);        // z_1 = dz0
        EXPECT_NEAR(static_cast<double>(in.terrain.hybrid_half_width), 80.0, 1e-12); // 2 dy
    }
    {
        auto in = make();
        in.mesh.dz0 = 0;   // uniform layers: 1000 / 20
        in.mesh.stretch_ratio = 0;
        in.mesh.stretch_ratio_given = false;
        in.mesh.n_cell = amrex::IntVect(64, 1, 20);   // one cell in y: only dx counts
        gust::apply_hybrid_defaults(in);
        EXPECT_NEAR(static_cast<double>(in.terrain.flat_height), 50.0, 1e-12);
        EXPECT_NEAR(static_cast<double>(in.terrain.hybrid_half_width), 32.0, 1e-12);
    }
    {
        auto in = make();
        in.terrain.flat_height = 0;
        in.terrain.hybrid_half_width = 16;   // one node in x, none in y
        in.terrain.given = {"flat_height", "hybrid_half_width"};
        EXPECT_TRUE(gust::apply_hybrid_defaults(in).empty());
        EXPECT_EQ(in.terrain.flat_height, 0);
        EXPECT_EQ(in.terrain.hybrid_half_width, 16);
    }
    {
        // Warnings: a given width below one cell in both directions, and flat terrain.
        auto in = make();
        in.terrain.hybrid_half_width = 15.9;
        in.terrain.given = {"flat_height", "hybrid_half_width"};
        auto log = gust::apply_hybrid_defaults(in);
        ASSERT_EQ(log.size(), 1u);
        EXPECT_NE(log[0].find("below one cell"), std::string::npos) << log[0];
        in.terrain.type = "flat";
        in.terrain.hybrid_half_width = 0;   // 0 is the documented terrain-following choice: no width warning
        log = gust::apply_hybrid_defaults(in);
        ASSERT_EQ(log.size(), 1u);
        EXPECT_NE(log[0].find("flat terrain"), std::string::npos) << log[0];
    }
    {
        auto in = make();
        in.terrain.representation = "immersed_forcing";
        in.terrain.flat_height = 7;
        EXPECT_TRUE(gust::apply_hybrid_defaults(in).empty());
        EXPECT_EQ(in.terrain.flat_height, 7);
        EXPECT_EQ(in.terrain.hybrid_half_width, 0);
    }
}

// The envelope on the mesh against the host envelope of the same nodes, for periodic and cut faces and
// boxes of 16 and 8 cells.
TEST (HybridTerrain, MeshEnvelopeDoesNotDependOnTheBoxes)
{
    for (const bool per : {true, false}) {
        for (const int mgs : {16, 8}) {
            SCOPED_TRACE(std::string(per ? "periodic" : "cut") + ", boxes of " + std::to_string(mgs));
            const PCase c = gust::lamtest::flat_case(32, 16, 4, 320, 160, 400, per, per, mgs);
            const auto dx = c.geom.CellSizeArray();
            auto hill = [] (double x, double y) { return 60.0 / (1 + ((x - 120) * (x - 120) + (y - 70) * (y - 70)) / 900.0) + 0.02 * x; };
            amrex::MultiFab z(c.z.boxArray(), c.dm, 1, c.z.nGrow());
            z.setVal(0);
            for (amrex::MFIter mfi(z, false); mfi.isValid(); ++mfi) {
                const auto a = z.array(mfi);
                const amrex::Box dom = amrex::surroundingNodes(c.geom.Domain());
                amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
                    if (k != 0) { return; }
                    // Ghost nodes: periodic image, else the nearest domain node (fill_terrain_surface's).
                    const int ii = per ? ((i % 32) + 32) % 32 : std::max(0, std::min(dom.bigEnd(0), i));
                    const int jj = per ? ((j % 16) + 16) % 16 : std::max(0, std::min(dom.bigEnd(1), j));
                    a(i, j, k) = static_cast<amrex::Real>(hill(ii * static_cast<double>(dx[0]), jj * static_cast<double>(dx[1])));
                });
            }
            gust::make_hybrid_envelope(z, c.geom, 30.0);   // 3 x 3 nodes
            const int nx = per ? 32 : 33, ny = per ? 16 : 17;
            std::vector<double> h(static_cast<std::size_t>(nx) * ny);
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    h[static_cast<std::size_t>(j) * nx + i] = static_cast<double>(static_cast<amrex::Real>(hill(i * static_cast<double>(dx[0]), j * static_cast<double>(dx[1]))));
                }
            }
            const auto e = gust::lower_envelope(h, nx, ny, 3, 3, per, per);
            double worst = 0;
            for (amrex::MFIter mfi(z, false); mfi.isValid(); ++mfi) {
                const auto a = z.const_array(mfi);
                amrex::Box bx = mfi.fabbox();
                bx.setSmall(2, 0);
                bx.setBig(2, 0);
                amrex::LoopOnCpu(bx, [&] (int i, int j, int) {
                    const int ii = per ? ((i % nx) + nx) % nx : std::max(0, std::min(nx - 1, i));
                    const int jj = per ? ((j % ny) + ny) % ny : std::max(0, std::min(ny - 1, j));
                    worst = std::max(worst, std::abs(static_cast<double>(a(i, j, 0)) - e[static_cast<std::size_t>(jj) * nx + ii]));
                });
            }
            amrex::ParallelDescriptor::ReduceRealMax(worst);
            EXPECT_LT(worst, (std::is_same_v<amrex::Real, double> ? 1e-12 : 1e-4));
        }
    }
}

// The surface layer's immersed columns (hybrid_immersed_columns): 1 where the first cell's solid fraction is
// above 0, ghost columns included, on the ground slab of the cells' boxes; the count is that of the valid
// columns over all ranks.
TEST (HybridTerrain, ImmersedColumnsWhereTheFirstCellIsCut)
{
    const PCase c = gust::lamtest::flat_case(32, 16, 4, 320, 160, 400, true, true, 8);
    amrex::MultiFab beta(c.ba, c.dm, 1, 1);
    for (amrex::MFIter mfi(beta, false); mfi.isValid(); ++mfi) {
        const auto b = beta.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
            // k = 0: a fraction above 0 in columns with i + j odd (one in five of them 1e-12), else 0;
            // above: 1 everywhere (must not be read).
            b(i, j, k) = k > 0 ? 1 : (((i + j) % 2 != 0) ? (((i + 32) % 5 == 0) ? 1e-12 : 0.4) : 0);
        });
    }
    amrex::MultiFab mask;
    const amrex::Long n = gust::hybrid_immersed_columns(beta, gust::ground_slab(c.ba), c.dm, 1, mask);
    EXPECT_EQ(n, 32 * 16 / 2);
    EXPECT_EQ(mask.nGrowVect(), amrex::IntVect(1, 1, 0));
    amrex::Long bad = 0;
    for (amrex::MFIter mfi(mask, false); mfi.isValid(); ++mfi) {
        const auto m = mask.const_array(mfi);
        EXPECT_EQ(mfi.validbox().smallEnd(2), 0);
        EXPECT_EQ(mfi.validbox().bigEnd(2), 0);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
            bad += (m(i, j, k) != (((i + j) % 2 != 0) ? 1 : 0)) ? 1 : 0;
        });
    }
    amrex::ParallelDescriptor::ReduceLongSum(bad);
    EXPECT_EQ(bad, 0);
}

// The hybrid mesh's pressure fill (hybrid_solid_pressure) on a mesh that follows a ramp of slope 0.3: the
// solid cells (beta = 1) lie below a cut cell (beta 0.4) at a height that changes from column to column,
// and hold 0, as the Poisson solution's solid cells do. The fluid holds p = a x, which has no vertical
// change. The fill gives every solid cell the fluid's value above it, so the filled field is a x at every
// cell of the fab (ghost columns and the layer below the ground included), and the x component of the
// pressure gradient is a on every x face of the fluid to the metric's error on the smoothed levels (about
// 0.5 % here). The field without the fill gives a slope term of order 0.3 a x / dz at the faces of the cut
// cells, more than a itself. A second field, a x + k, checks that the value is the cut cell's (the first
// cell above that is not solid), not the top cell's or one above the cut cell.
TEST (HybridTerrain, SolidPressureMakesTheSlopedGradientExact)
{
    constexpr double a = 2.5e-3;
    for (const int mgs : {16, 8}) {
        SCOPED_TRACE("boxes of " + std::to_string(mgs));
        const PCase c = gust::lamtest::flat_case(32, 8, 16, 320, 80, 160, false, true, mgs,
                                                 [] (double x, double) { return 0.3 * x; });
        const auto dx = c.geom.CellSizeArray();
        auto ground = [] (int i, int j) { return 1 + ((i + 64) + 2 * (j + 64)) % 3; };   // the cut cell's k
        amrex::MultiFab beta(c.ba, c.dm, 1, 1), p(c.ba, c.dm, 1, 1), ps, pk(c.ba, c.dm, 1, 1), pks;
        for (amrex::MFIter mfi(p, false); mfi.isValid(); ++mfi) {
            const auto b = beta.array(mfi);
            const auto q = p.array(mfi);
            const auto qk = pk.array(mfi);
            amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
                const int g = ground(i, j);
                b(i, j, k) = k < g ? 1 : (k == g ? amrex::Real(0.4) : 0);
                q(i, j, k) = k < g ? 0 : static_cast<amrex::Real>(a * (i + 0.5) * static_cast<double>(dx[0]));
                qk(i, j, k) = k < g ? 0 : q(i, j, k) + static_cast<amrex::Real>(k);
            });
        }
        gust::hybrid_solid_pressure(beta, p, ps);
        gust::hybrid_solid_pressure(beta, pk, pks);
        ASSERT_EQ(ps.nGrowVect(), p.nGrowVect());
        double worst = 0;
        amrex::Long filled = 0;
        for (amrex::MFIter mfi(ps, false); mfi.isValid(); ++mfi) {
            const auto f = ps.const_array(mfi);
            const auto fk = pks.const_array(mfi);
            const amrex::Box vb = mfi.validbox();
            amrex::Box bx = amrex::grow(vb, amrex::IntVect(1, 1, 0));
            bx.growLo(2, 1);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                worst = std::max(worst, std::abs(static_cast<double>(f(i, j, k)) - static_cast<double>(static_cast<amrex::Real>(a * (i + 0.5) * static_cast<double>(dx[0])))));
                const int g = ground(i, j);
                filled += (k < g) ? 1 : 0;
                const double want = static_cast<double>(f(i, j, k)) + std::max(k, g);
                worst = std::max(worst, std::abs(static_cast<double>(fk(i, j, k)) - want));
            });
        }
        amrex::ParallelDescriptor::ReduceRealMax(worst);
        amrex::ParallelDescriptor::ReduceLongSum(filled);
        EXPECT_LT(worst, 8 * std::numeric_limits<amrex::Real>::epsilon() * 20);
        EXPECT_GT(filled, 0);

        auto gradient_error = [&] (const amrex::MultiFab& pf) {
            amrex::Array<amrex::MultiFab, 3> G;
            for (int d = 0; d < 3; ++d) {
                G[d].define(amrex::convert(c.ba, amrex::IntVect::TheDimensionVector(d)), c.dm, 1, 0);
            }
            gust::pressure_gradient(pf, c.z, c.mesh, c.geom, {&G[0], &G[1], &G[2]});
            double e = 0;
            for (amrex::MFIter mfi(G[0], false); mfi.isValid(); ++mfi) {
                const auto g = G[0].const_array(mfi);
                amrex::LoopOnCpu(mfi.validbox(), [&] (int i, int j, int k) {
                    // Interior x faces between two cells of the fluid or its cut cells.
                    if (i == 0 || i == 32 || k < std::max(ground(i - 1, j), ground(i, j))) { return; }
                    e = std::max(e, std::abs(static_cast<double>(g(i, j, k)) - a));
                });
            }
            amrex::ParallelDescriptor::ReduceRealMax(e);
            return e;
        };
        EXPECT_LT(gradient_error(ps), 0.01 * a);
        EXPECT_GT(gradient_error(p), a) << "the field without the fill must show the slope term";
    }
}
