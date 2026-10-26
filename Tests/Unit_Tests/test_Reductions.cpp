// Deterministic sum and max, identical to the last bit on 1, 2 and
// 4 ranks. CTest runs this file's tests at each rank count; every run must
// reproduce the same pinned bits.
#include "Reductions.H"
#include "Residual.H"
#include "Simplec.H"
#include "State.H"
#include "TestHelpers.H"

#include <gtest/gtest.h>

#include <cstdint>
#include <cmath>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace {

using gust::test::fill_pattern;
using gust::test::pattern;

// 24 x 20 x 16 cells cut 8 x 8 in x/y: 3 x 3 boxes (the y boxes are 8, 8, 4).
gust::MeshInputs reduction_mesh ()
{
    gust::MeshInputs m;
    m.n_cell = amrex::IntVect(24, 20, 16);
    m.prob_hi = {1, 1, 1};
    m.max_grid_size = {8, 8};
    return m;
}

std::uint64_t bits (double x)
{
    std::uint64_t b = 0;
    std::memcpy(&b, &x, sizeof(b));
    return b;
}

// -0x1.f9a008f38d09cp+19, pinned from the 1-rank run (macOS arm64, AppleClang);
// the 2- and 4-rank runs, and Linux GCC/Clang, must reproduce it. The addends
// are dyadic and the only operations are additions, so no FMA contraction or
// libm difference can change it.
constexpr std::uint64_t kSumBits = 0xc12f9a008f38d09cULL;

TEST (Reductions, SumIsIdenticalOnEveryRankCount)
{
    if constexpr (!std::is_same_v<amrex::Real, double>) {
        GTEST_SKIP() << "pinned value is for double precision";
    } else {
        const auto ba = gust::make_box_array(reduction_mesh());
        ASSERT_EQ(ba.size(), 9);
        const amrex::DistributionMapping dm(ba);
        amrex::MultiFab mf(ba, dm, 1, 1);
        fill_pattern(mf, 3);

        const double s = gust::deterministic_sum(mf, 0);
        EXPECT_EQ(bits(s), kSumBits) << std::hexfloat << s;

        // The same value with a different DistributionMapping of the same boxes.
        amrex::Vector<int> pmap(ba.size());
        for (int b = 0; b < ba.size(); ++b) {
            pmap[b] = (ba.size() - 1 - b) % amrex::ParallelDescriptor::NProcs();
        }
        const amrex::DistributionMapping dm_rev(pmap);
        amrex::MultiFab mf_rev(ba, dm_rev, 1, 1);
        fill_pattern(mf_rev, 3);
        EXPECT_EQ(bits(gust::deterministic_sum(mf_rev, 0)), bits(s));

        // Accuracy: within 16 eps sum|x| of the exact sum, pinned from Python's
        // math.fsum of the same addends (-0x1.f9a008f38d094p+19).
        double abs_sum = 0;
        for (int k = 0; k < 16; ++k) {
            for (int j = 0; j < 20; ++j) {
                for (int i = 0; i < 24; ++i) { abs_sum += std::abs(static_cast<double>(pattern(i, j, k, 3))); }
            }
        }
        const double exact = -0x1.f9a008f38d094p+19;
        EXPECT_LE(std::abs(s - exact), 16 * std::numeric_limits<double>::epsilon() * abs_sum);
    }
}

TEST (Reductions, SumIgnoresGhostCells)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab a(ba, dm, 1, 1), b(ba, dm, 1, 1);
    fill_pattern(a, 3);
    fill_pattern(b, 3);
    b.setBndry(1.0e30);
    EXPECT_EQ(gust::deterministic_sum(a, 0), gust::deterministic_sum(b, 0));
}

TEST (Reductions, MaxEqualsTheExactMaximum)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab mf(ba, dm, 1, 1);
    fill_pattern(mf, 5);
    mf.setBndry(1.0e30);   // ghost cells must not count

    amrex::Real ref = std::numeric_limits<amrex::Real>::lowest();
    for (int k = 0; k < 16; ++k) {
        for (int j = 0; j < 20; ++j) {
            for (int i = 0; i < 24; ++i) { ref = std::max(ref, pattern(i, j, k, 5)); }
        }
    }
    EXPECT_EQ(gust::global_max(mf, 0), ref);
}

TEST (Reductions, MaxReportsNaNAndPositiveZero)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab mf(ba, dm, 1, 1);

    mf.setVal(-0.0);
    const amrex::Real z = gust::global_max(mf, 0);
    EXPECT_EQ(z, 0);
    EXPECT_FALSE(std::signbit(z));

    mf.setVal(-1.0);
    gust::test::set_one_nan(mf, amrex::IntVect(17, 3, 9));
    EXPECT_TRUE(std::isnan(gust::global_max(mf, 0)));
}

namespace {

//! Write v at the point iv of whichever box holds it in its valid region
//! (ghost = false), or in its first ghost layer outside the domain (ghost = true).
void set_point (amrex::MultiFab& mf, const amrex::IntVect& iv, amrex::Real v, bool ghost)
{
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box vb = mfi.validbox();
        if (ghost ? !(amrex::grow(vb, 1).contains(iv) && !vb.contains(iv)) : !vb.contains(iv)) { continue; }
        const auto a = mf.array(mfi);
        amrex::ParallelFor(amrex::Box(iv, iv), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { a(i, j, k) = v; });
        if (ghost) { break; }   // one box's ghost only
    }
}

} // namespace

// The validity check that replaced a min and a max per field (audit A7, I3):
// every NaN, infinity and negative value is counted, zero only when zero is
// not allowed (-0 counts as zero), ghost layers only when asked, and every
// rank gets the full count. The bad points sit in different boxes, so on 2
// and 4 ranks they are on different ranks. No ordered compare touches the
// NaN: on Linux the invalid trap is on, and the count must not raise it.
TEST (Reductions, CountNotPositiveFindsEveryBadPoint)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab a(ba, dm, 1, 1), b(ba, dm, 1, 1);
    a.setVal(1.0);
    b.setVal(amrex::Real(1.0e-30));
    const auto count = [&] (const std::vector<const amrex::MultiFab*>& f, bool allow_zero, int ng) {
        const std::vector<amrex::Long> n = gust::count_not_positive(f, allow_zero, ng);
        EXPECT_EQ(n.size(), f.size());
        return n;
    };
    EXPECT_EQ(count({&a, &b}, false, 1), (std::vector<amrex::Long>{0, 0}));
    const amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    set_point(a, amrex::IntVect(0, 0, 0), nan, false);
    set_point(a, amrex::IntVect(23, 19, 15), inf, false);
    set_point(a, amrex::IntVect(12, 10, 3), -inf, false);
    set_point(a, amrex::IntVect(5, 17, 8), -std::numeric_limits<amrex::Real>::denorm_min(), false);
    set_point(b, amrex::IntVect(20, 2, 11), 0.0, false);
    set_point(b, amrex::IntVect(9, 9, 9), -0.0, false);
    set_point(b, amrex::IntVect(-1, 4, 4), -1.0, true);   // a ghost cell outside the domain
    // > 0: the NaN, both infinities, the negative, both zeros; ghosts only with ngrow = 1
    EXPECT_EQ(count({&a, &b}, false, 0), (std::vector<amrex::Long>{4, 2}));
    EXPECT_EQ(count({&a, &b}, false, 1), (std::vector<amrex::Long>{4, 3}));
    // >= 0: the zeros pass
    EXPECT_EQ(count({&a, &b}, true, 0), (std::vector<amrex::Long>{4, 0}));
    EXPECT_EQ(count({&a, &b}, true, 1), (std::vector<amrex::Long>{4, 1}));
    EXPECT_EQ(count({&b}, true, 0), (std::vector<amrex::Long>{0}));
    EXPECT_EQ(count({&b}, false, 0), (std::vector<amrex::Long>{2}));
    // require_positive refuses on every rank, naming the operator and the failing fields only
    try {
        gust::require_positive("op_x", {"field_a", "field_b"}, {&a, &b}, true);
        ADD_FAILURE() << "require_positive accepted a NaN";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("gust: op_x: field_a must be finite and >= 0 everywhere (4 points are not)"),
                  std::string::npos) << e.what();
    }
    try {
        gust::require_positive("op_x", {"field_a", "field_b"}, {&a, &b}, false);
        ADD_FAILURE() << "require_positive accepted a NaN";
    } catch (const std::runtime_error& e) {
        EXPECT_NE(std::string(e.what()).find("gust: op_x: field_a, field_b must be finite and > 0 everywhere (6 points are not)"),
                  std::string::npos) << e.what();
    }
    EXPECT_NO_THROW(gust::require_positive("op_x", {"field_b"}, {&b}, true));
}

// count_bit_differences compares bit patterns: +0 and -0 differ, a NaN
// equals the same NaN, a NaN differs from a number; every component counts,
// and ghost cells only up to ngrow. a against itself: 0.
TEST (Reductions, CountBitDifferencesComparesBitPatterns)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab a(ba, dm, 2, 1), b(ba, dm, 2, 1);
    a.setVal(1.5);
    b.setVal(1.5);
    const amrex::Real nan = std::numeric_limits<amrex::Real>::quiet_NaN();
    set_point(a, amrex::IntVect(3, 4, 5), nan, false);
    set_point(b, amrex::IntVect(3, 4, 5), nan, false);   // same NaN: equal
    set_point(a, amrex::IntVect(1, 1, 1), 0.0, false);
    set_point(b, amrex::IntVect(1, 1, 1), -0.0, false);  // +0 vs -0: differ
    set_point(a, amrex::IntVect(23, 19, 15), nan, false); // NaN vs 1.5: differ
    set_point(b, amrex::IntVect(-1, 0, 0), 2.0, true);   // a ghost cell
    EXPECT_EQ(gust::count_bit_differences(a, a, amrex::IntVect(1)), 0);
    EXPECT_EQ(gust::count_bit_differences(a, b, amrex::IntVect(0)), 2);
    EXPECT_EQ(gust::count_bit_differences(a, b, amrex::IntVect(1)), 3);
    // The second component: a point there counts once more.
    b.setVal(amrex::Real(1.75), 1, 1, 0);
    EXPECT_EQ(gust::count_bit_differences(a, b, amrex::IntVect(0)), 2 + ba.numPts());
}

TEST (Reductions, SumRejectsFaceData)
{
    const auto ba = gust::make_box_array(reduction_mesh());
    const amrex::DistributionMapping dm(ba);
    amrex::MultiFab u(amrex::convert(ba, amrex::IntVect(1, 0, 0)), dm, 1, 0);
    u.setVal(1.0);
    EXPECT_THROW(gust::deterministic_sum(u, 0), std::runtime_error);
}

} // namespace

// ---------------------------------------------------------------------------
// Audit A4 (I4, I5, U2): sums over face data count each face once and do not
// depend on the box layout or the ranks. An x-face field of a wrapped pattern
// (face n equals face 0 when x is periodic, as the solver keeps them), x
// periodic and not, on 3 x 3 uneven boxes, one box and 6-cell boxes, with the
// boxes on their ranks and in reverse.
namespace {

gust::MeshInputs face_mesh (int perx, int max_grid)
{
    gust::MeshInputs m = reduction_mesh();
    m.is_periodic = {perx, 0, 0};
    m.max_grid_size = {max_grid, max_grid};
    return m;
}

//! The pattern at x face i, wrapped when x is periodic (nx cells).
double face_value (int i, int j, int k, int nx, bool perx)
{
    return static_cast<double>(pattern(perx ? i % nx : i, j, k, 2));
}

void fill_faces (amrex::MultiFab& u, int nx, bool perx)
{
    for (amrex::MFIter mfi(u, false); mfi.isValid(); ++mfi) {
        const auto a = u.array(mfi);
        amrex::LoopOnCpu(mfi.fabbox(), [&] (int i, int j, int k) {
            a(i, j, k) = static_cast<amrex::Real>(face_value(i, j, k, nx, perx));
        });
    }
}

struct FaceCase
{
    amrex::Geometry geom;
    amrex::BoxArray ba;
    amrex::MultiFab u, u_rev;
    std::unique_ptr<amrex::iMultiFab> own, own_rev;
    FaceCase (int perx, int max_grid)
    {
        const auto m = face_mesh(perx, max_grid);
        geom = gust::make_geometry(m);
        ba = amrex::convert(gust::make_box_array(m), amrex::IntVect(1, 0, 0));
        const amrex::DistributionMapping dm(ba);
        amrex::Vector<int> pmap(ba.size());
        for (int b = 0; b < ba.size(); ++b) { pmap[b] = (ba.size() - 1 - b) % amrex::ParallelDescriptor::NProcs(); }
        const amrex::DistributionMapping dm_rev(pmap);
        u.define(ba, dm, 1, 0);
        u_rev.define(ba, dm_rev, 1, 0);
        fill_faces(u, 24, perx != 0);
        fill_faces(u_rev, 24, perx != 0);
        own = u.OwnerMask(geom.periodicity());
        own_rev = u_rev.OwnerMask(geom.periodicity());
    }
};

} // namespace

TEST (Reductions, DotCountsEachFaceOnceOnAnyLayout)
{
    for (const int perx : {0, 1}) {
        // The exact sum over the distinct faces: i = 0..24, or 0..23 periodic.
        double exact = 0;
        for (int k = 0; k < 16; ++k) {
            for (int j = 0; j < 20; ++j) {
                for (int i = 0; i <= (perx ? 23 : 24); ++i) { exact += face_value(i, j, k, 24, perx != 0) * face_value(i, j, k, 24, perx != 0); }
            }
        }
        for (const int mg : {8, 64, 6}) {
            FaceCase fc(perx, mg);
            const double s = gust::deterministic_dot(fc.u, 0, fc.u, 0, fc.own.get());
            EXPECT_EQ(bits(gust::deterministic_dot(fc.u_rev, 0, fc.u_rev, 0, fc.own_rev.get())), bits(s)) << perx << " " << mg;
            EXPECT_NEAR(s, exact, 1e-13 * exact) << perx << " " << mg;
        }
    }
}

TEST (Reductions, LevelSumsCountEachFaceOnceOnAnyLayout)
{
    for (const int perx : {0, 1}) {
        for (const int mg : {8, 64, 6}) {
            FaceCase fc(perx, mg);
            const auto s = gust::deterministic_level_sums(fc.u, 0, *fc.own, fc.geom.Domain());
            const auto s_rev = gust::deterministic_level_sums(fc.u_rev, 0, *fc.own_rev, fc.geom.Domain());
            ASSERT_EQ(s.size(), 16U);
            for (int k = 0; k < 16; ++k) {
                double exact = 0, mag = 0;
                for (int j = 0; j < 20; ++j) {
                    for (int i = 0; i <= (perx ? 23 : 24); ++i) {
                        exact += face_value(i, j, k, 24, perx != 0);
                        mag += std::abs(face_value(i, j, k, 24, perx != 0));
                    }
                }
                EXPECT_EQ(bits(s_rev[k]), bits(s[k])) << perx << " " << mg << " " << k;
                EXPECT_NEAR(s[k], exact, 1e-13 * mag) << perx << " " << mg << " " << k;
            }
        }
    }
}

// The scaled residual RMS of every field (stage_norms) over its distinct
// points: the same on the three layouts to 1e-13, and for u the exact RMS.
// Before, Dot counted a box-seam face twice and both copies of a periodic
// face, so the RMS changed with mesh.max_grid_size.
TEST (Reductions, ResidualRmsDoesNotDependOnTheLayout)
{
    for (const int perx : {0, 1}) {
        double exact = 0;
        for (int k = 0; k < 16; ++k) {
            for (int j = 0; j < 20; ++j) {
                for (int i = 0; i <= (perx ? 23 : 24); ++i) { exact += face_value(i, j, k, 24, perx != 0) * face_value(i, j, k, 24, perx != 0); }
            }
        }
        const double npts = (perx ? 24.0 : 25.0) * 20 * 16;
        const gust::JvpScales sc;
        const double su = static_cast<double>(sc.velocity * sc.velocity);   // U^2 / H with H = 1
        std::vector<amrex::Real> ref;
        for (const int mg : {8, 64, 6}) {
            const auto m = face_mesh(perx, mg);
            const auto geom = gust::make_geometry(m);
            const auto ba = gust::make_box_array(m);
            const amrex::DistributionMapping dm(ba);
            gust::ResidualFields R;
            R.define(ba, dm);
            fill_faces(R.u, 24, perx != 0);
            fill_pattern(R.v, 4);
            fill_pattern(R.w, 5);
            fill_pattern(R.p, 6);
            fill_pattern(R.theta, 7);
            fill_pattern(R.tke, 8);
            gust::State q;
            q.define(ba, dm);
            for (amrex::MultiFab* f : {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke}) { f->setVal(1.0); }
            std::vector<amrex::Real> res, mon;
            gust::stage_norms(R, q, nullptr, sc, 1, geom, res, mon);
            const double rt = std::is_same_v<amrex::Real, double> ? 1e-13 : 1e-5;   // the result is a Real
            EXPECT_NEAR(static_cast<double>(res[0]), std::sqrt(exact / npts) / su, rt * std::sqrt(exact / npts) / su) << perx << " " << mg;
            if (ref.empty()) { ref = res; continue; }
            for (int n = 0; n < 6; ++n) {
                EXPECT_NEAR(static_cast<double>(res[n]), static_cast<double>(ref[n]), rt * std::abs(static_cast<double>(ref[n])))
                    << perx << " " << mg << " " << n;
            }
        }
    }
}
