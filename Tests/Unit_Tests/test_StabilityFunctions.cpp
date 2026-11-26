// The stable side of MOST's stability functions (StabilityFunctions.H, bc.surface_stable_functions).
// Businger-Dyer (the default) is ERF's calc_phi/calc_psi exactly; QNSE (Sukoriansky et al. 2005, ERF's
// enable_qnse_stable_functions) is pinned to values computed independently in Python, its psi to a numerical
// integral of (1 - phi)/x, and its unstable side is Businger-Dyer's. The Jimenez pair (surface_temperature only) is
// pinned to Python values and to ERF's calc_psi_m2/calc_psi_h2 on both sides.

#include "ColumnGuess.H"
#include "ERF_MOSTUtils.H"
#include "Initializer.H"
#include "StabilityFunctions.H"
#include "SurfaceLayer.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

namespace {

using gust::StableCoefficients;
using gust::SurfaceLayerOptions;
using gust::trtest::tol;

const StableCoefficients kQnse{amrex::Real(2.5), amrex::Real(0.2), amrex::Real(2.5), amrex::Real(0.2)};

bool mentions (const std::vector<std::string>& err, const std::string& what)
{
    for (const auto& e : err) { if (e.find(what) != std::string::npos) { return true; } }
    return false;
}

} // namespace

// QNSE with ERF's coefficients: phi = (1 + 2.5 zeta)/(1 + 0.2 zeta) and psi = -11.5 ln(1 + 0.2 zeta) at
// zeta = 0.1, 1, 10, 100 (Python); phi tends to 12.5; distinct coefficients per function. psi is the integral of (1 - phi_m)/x (Simpson, 2000 panels).
TEST(StabilityFunctions, QnseStableSide)
{
    struct Row { double zeta, phi, psi; };
    const Row rows[] = {{0.1, 1.2254901960784315, -0.2277302139060669}, {1, 2.916666666666667, -2.0966979031304778},
                        {10, 8.666666666666666, -12.634041319683263}, {100, 11.952380952380953, -35.01200803381936}};
    for (const Row& r : rows) {
        const auto z = static_cast<amrex::Real>(r.zeta);
        amrex::Real pm, ph;
        gust::most_phi(z, kQnse, pm, ph);
        EXPECT_NEAR(pm, r.phi, tol(1e-14, 1e-6) * r.phi) << r.zeta;
        EXPECT_NEAR(ph, r.phi, tol(1e-14, 1e-6) * r.phi) << r.zeta;
        EXPECT_NEAR(gust::most_psi_m<true>(z, kQnse), r.psi, tol(1e-13, 1e-5) * std::abs(r.psi)) << r.zeta;
        EXPECT_NEAR(gust::most_psi_h<true>(z, kQnse), r.psi, tol(1e-13, 1e-5) * std::abs(r.psi)) << r.zeta;
    }
    // Distinct coefficients (a_m, b_m, a_h, b_h) = (3, 0.1, 4, 0.5) at zeta = 2: each function takes its own pair.
    {
        const StableCoefficients d{3, amrex::Real(0.1), 4, amrex::Real(0.5)};
        amrex::Real pm, ph;
        gust::most_phi(2, d, pm, ph);
        EXPECT_NEAR(pm, 5.833333333333334, tol(1e-14, 1e-6) * 6);
        EXPECT_NEAR(ph, 4.5, tol(1e-14, 1e-6) * 5);
        EXPECT_NEAR(gust::most_psi_m<true>(2, d), -5.287325147024682, tol(1e-13, 1e-5) * 6);
        EXPECT_NEAR(gust::most_psi_h<true>(2, d), -4.852030263919617, tol(1e-13, 1e-5) * 5);
    }
    for (const double zeta : {0.5, 3.0}) {
        const int n = 2000;
        const double hstep = zeta / n;
        double s = 0;
        for (int i = 0; i <= n; ++i) {
            const double x = i * hstep;
            amrex::Real pm, ph;
            gust::most_phi(static_cast<amrex::Real>(x), kQnse, pm, ph);
            const double fx = (i == 0) ? -(2.5 - 0.2) : (1 - static_cast<double>(pm)) / x;   // the limit at 0
            s += fx * ((i == 0 || i == n) ? 1 : (i % 2 ? 4 : 2));
        }
        s *= hstep / 3;
        EXPECT_NEAR(gust::most_psi_m<true>(static_cast<amrex::Real>(zeta), kQnse), s, tol(1e-10, 1e-4)) << zeta;
    }
}

// The unstable side is Businger-Dyer's whatever the stable side; the default coefficients give ERF's
// Businger-Dyer phi and psi on both sides (bit for bit for psi: the same function).
TEST(StabilityFunctions, BusingerDyerAndTheUnstableSide)
{
    const similarity_funs sf;
    for (double zeta = -20; zeta <= 20; zeta += 0.25) {
        const auto z = static_cast<amrex::Real>(zeta);
        amrex::Real pq, hq, pb, hb;
        gust::most_phi(z, kQnse, pq, hq);
        gust::most_phi(z, StableCoefficients{}, pb, hb);
        EXPECT_NEAR(pb, sf.calc_phi_m(z), tol(1e-14, 1e-6) * (1 + std::abs(double(pb)))) << zeta;
        EXPECT_NEAR(hb, sf.calc_phi_h(z), tol(1e-14, 1e-6) * (1 + std::abs(double(hb)))) << zeta;
        EXPECT_EQ(gust::most_psi_m<false>(z, kQnse), sf.calc_psi_m(z)) << zeta;
        EXPECT_EQ(gust::most_psi_h<false>(z, kQnse), sf.calc_psi_h(z)) << zeta;
        if (zeta <= 0) {
            EXPECT_EQ(pq, pb) << zeta;
            EXPECT_EQ(hq, hb) << zeta;
            EXPECT_EQ(gust::most_psi_m<true>(z, kQnse), sf.calc_psi_m(z)) << zeta;
            EXPECT_EQ(gust::most_psi_h<true>(z, kQnse), sf.calc_psi_h(z)) << zeta;
        } else {
            EXPECT_LT(pq, pb) << zeta;   // QNSE mixes more than Businger-Dyer in stable air
            EXPECT_GT(gust::most_psi_m<true>(z, kQnse), sf.calc_psi_m(z)) << zeta;
        }
    }
}

// The options: qnse needs b > 0 and a >= b; the coefficients without qnse, and qnse with nothing stable to act
// on, are refused naming the key.
TEST(StabilityFunctions, InputChecks)
{
    SurfaceLayerOptions o;
    o.stable_functions = SurfaceLayerOptions::StableFunctions::qnse;
    EXPECT_TRUE(gust::check_surface_layer_options(o).empty());
    const StableCoefficients c = o.stable_coefficients();
    EXPECT_EQ(c.am, amrex::Real(2.5));
    EXPECT_EQ(c.bh, amrex::Real(0.2));
    o.qnse_bm = 0;
    EXPECT_TRUE(mentions(gust::check_surface_layer_options(o), "qnse_bm"));
    o.qnse_bm = amrex::Real(0.2);
    o.qnse_ah = amrex::Real(0.1);
    EXPECT_TRUE(mentions(gust::check_surface_layer_options(o), "qnse_ah"));
    SurfaceLayerOptions bd;
    EXPECT_EQ(bd.stable_coefficients().am, amrex::Real(5));
    EXPECT_EQ(bd.stable_coefficients().bm, amrex::Real(0));
    bd.qnse_am = 3;
    EXPECT_TRUE(mentions(gust::check_surface_layer_options(bd), "qnse_am is set but stable_functions is businger_dyer"));

    gust::InitOptions io;
    SurfaceLayerOptions q;
    q.stable_functions = SurfaceLayerOptions::StableFunctions::qnse;
    EXPECT_TRUE(mentions(gust::check_stable_functions(io, q, true), "bc.surface_stable_functions = qnse has no effect"));
    io.first_guess = gust::InitOptions::FirstGuess::column;
    io.column_mol = 200;
    EXPECT_TRUE(gust::check_stable_functions(io, q, true).empty());
    gust::InitOptions none;
    q.stability = SurfaceLayerOptions::Stability::heat_flux;
    EXPECT_TRUE(gust::check_stable_functions(none, q, true).empty());
}

// The column's wall gradients take the stable side: u* phi_m/(kappa z) and theta* phi_h/(kappa z) with QNSE's phi.
TEST(StabilityFunctions, ColumnWallGradients)
{
    const double kappa = 0.41, us = 0.3, ts = 0.05, z = 2.0, il = 0.5;   // zeta = 1
    double dU = 0, dT = 0;
    gust::column_wall_gradients(us, ts, il, z, 0.0, kappa, dU, dT, kQnse);
    EXPECT_NEAR(dU, us * 2.916666666666667 / (kappa * z), tol(1e-12, 1e-6));
    EXPECT_NEAR(dT, ts * 2.916666666666667 / (kappa * z), tol(1e-12, 1e-6));
}

// bc.surface_stable_functions = jimenez, ERF's Jimenez et al. (2012) pair. Stable side Cheng and Brutsaert (2005): psi_m = -6.1 ln(zeta + (1 +
// zeta^2.5)^(1/2.5)), psi_h = -5.3 ln(zeta + (1 + zeta^1.1)^(1/1.1)) and phi = 1 - zeta psi', pinned at zeta = 0.1,
// 1, 10, 100 to values computed in Python (mpmath, 30 digits, phi by its numerical derivative), and equal to ERF's
// calc_psi_m2/calc_psi_h2 (extern_erf/ERF_MOSTUtils.H). The unstable side (the Businger-Dyer and free-convection
// blend) likewise; zeta = 0 gives psi = 0 and phi = 1.
TEST(StabilityFunctions, JimenezPair)
{
    struct Row { double zeta, psi_m, psi_h, phi_m, phi_h; };
    const Row rows[] = {
        {0.1, -0.5883959379291539, -0.84098277137607233, 1.5713921003566282, 1.8089697833695153},
        {1, -5.1322658401041618, -5.6023522548976726, 5.3649340796961089, 4.5708225288528175},
        {10, -18.277819976407348, -16.064719899140247, 7.090379385815271, 6.098220471076017},
        {100, -32.319748135894224, -28.096256286095339, 7.0999695002439979, 6.2833369639835132}};
    const similarity_funs sf;
    for (const Row& r : rows) {
        const auto z = static_cast<amrex::Real>(r.zeta);
        amrex::Real pm, ph;
        gust::most_phi_jimenez(z, pm, ph);
        EXPECT_NEAR(gust::most_psi_m_jimenez(z), r.psi_m, tol(1e-13, 1e-5) * std::abs(r.psi_m)) << r.zeta;
        EXPECT_NEAR(gust::most_psi_h_jimenez(z), r.psi_h, tol(1e-13, 1e-5) * std::abs(r.psi_h)) << r.zeta;
        EXPECT_NEAR(pm, r.phi_m, tol(1e-12, 1e-5) * r.phi_m) << r.zeta;
        EXPECT_NEAR(ph, r.phi_h, tol(1e-12, 1e-5) * r.phi_h) << r.zeta;
        EXPECT_NEAR(gust::most_psi_m_jimenez(z), sf.calc_psi_m2(z), tol(1e-13, 1e-5) * std::abs(r.psi_m)) << r.zeta;
        EXPECT_NEAR(gust::most_psi_h_jimenez(z), sf.calc_psi_h2(z), tol(1e-13, 1e-5) * std::abs(r.psi_h)) << r.zeta;
    }
    // Unstable side: the blend, pinned at zeta = -0.1, -1, -5, -50 (mpmath, ERF's formulas), and ERF's
    // calc_psi_m2/calc_psi_h2; phi by the numerical derivative of that psi.
    const Row urows[] = {
        {-0.1, 0.28336310869192893, 0.53513129873064347, 0.78806858654634916, 0.61839727881579906},
        {-1, 1.1226456679594171, 1.8882809525335583, 0.46464826857294017, 0.26706949856546559},
        {-5, 2.1662698828511077, 3.1259808982560125, 0.26457943850303298, 0.18465386239021944},
        {-50, 4.0376028880047439, 5.135011778146106, 0.12573219345306984, 0.083939696564407946}};
    for (const Row& r : urows) {
        const auto z = static_cast<amrex::Real>(r.zeta);
        amrex::Real pm, ph;
        gust::most_phi_jimenez(z, pm, ph);
        EXPECT_NEAR(gust::most_psi_m_jimenez(z), r.psi_m, tol(1e-13, 1e-5) * std::abs(r.psi_m)) << r.zeta;
        EXPECT_NEAR(gust::most_psi_h_jimenez(z), r.psi_h, tol(1e-13, 1e-5) * std::abs(r.psi_h)) << r.zeta;
        EXPECT_NEAR(pm, r.phi_m, tol(1e-12, 1e-5) * r.phi_m) << r.zeta;
        EXPECT_NEAR(ph, r.phi_h, tol(1e-12, 1e-5) * r.phi_h) << r.zeta;
        EXPECT_NEAR(gust::most_psi_m_jimenez(z), sf.calc_psi_m2(z), tol(1e-13, 1e-5) * std::abs(r.psi_m)) << r.zeta;
        EXPECT_NEAR(gust::most_psi_h_jimenez(z), sf.calc_psi_h2(z), tol(1e-13, 1e-5) * std::abs(r.psi_h)) << r.zeta;
    }
    amrex::Real pm, ph;
    gust::most_phi_jimenez(0, pm, ph);
    EXPECT_EQ(pm, 1);
    EXPECT_EQ(ph, 1);
    EXPECT_EQ(gust::most_psi_m_jimenez(0), 0);
    EXPECT_EQ(gust::most_psi_h_jimenez(0), 0);
}

// jimenez is refused with the other stabilities (its kernels are the surface-temperature wall model's) and
// accepted with surface_temperature.
TEST(StabilityFunctions, JimenezOnlyWithSurfaceTemperature)
{
    SurfaceLayerOptions o;
    o.stable_functions = SurfaceLayerOptions::StableFunctions::jimenez;
    EXPECT_TRUE(mentions(gust::check_surface_layer_options(o), "jimenez is available only"));
    o.stability = SurfaceLayerOptions::Stability::heat_flux;
    EXPECT_TRUE(mentions(gust::check_surface_layer_options(o), "jimenez is available only"));
    o.stability = SurfaceLayerOptions::Stability::surface_temperature;
    o.surface_temperature = 290;
    EXPECT_TRUE(gust::check_surface_layer_options(o).empty());
}
