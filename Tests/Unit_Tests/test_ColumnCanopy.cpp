// The column first guess with a forest canopy (ColumnGuess.H, ColumnCanopy), plan rev. 6 beyond the plan,
// sub-phase 8.5b (UT-8.5.1's 1-D part): the force-driven neutral column of 64 x 1 m layers (z0 = 0.1 m, a pressure
// gradient F = 2.5e-3 m/s^2, no rotation, a Rayleigh layer on u, v and k above 48 m) under a dense uniform canopy,
// f = 0.5 1/m (Cd LAI = 20) over the lowest 40 layers. Deep inside such a canopy the shear and the diffusion vanish, so
// the steady balances are closed forms:
//   momentum  F = f |W| |W|                               ->  |W| = sqrt(F / f) = 0.0707 m/s;
//   k         beta_p |W|^3 = beta_d |W| k + Cmu0^3 k^1.5 / alpha   (P_c = D_c + eps with L = alpha / f, f cancels),
// and with canopy.tke = off the k equation has no canopy sink or cap there, so k is what diffuses down from the
// canopy top against the dissipation Cmu0^3 k^1.5 / (kappa (z + z0)): far above the closed form with the terms on
// (8.6e-4 m^2/s^2, the root above, against 0.4 to 0.6 at 10 to 30 m, measured), and its diffusion then lifts the speed above
// sqrt(F/f), by 34 % at 10 m to about 110 % at 20 m (the layers checked). A canopy of f = 0 everywhere gives
// the no-canopy table bit for bit; a profile that is not one value per layer is refused.
#include "Closure.H"
#include "ColumnGuess.H"
#include "Initializer.H"
#include "TestHelpers.H"
#include "TransportTestHelpers.H"

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

namespace {

using gust::BoundaryInputs;
using gust::FaceType;
using gust::InitOptions;
using gust::trtest::tol;

constexpr double F = 2.5e-3, FDRAG = 0.5;
constexpr int NZ = 64, NCANOPY = 40;

struct Column
{
    gust::Inputs in;
    InitOptions io;
    gust::SteadyOptions so;
    BoundaryInputs bin;
    gust::SurfaceLayerOptions surf;
    gust::ClosureOptions copt;
    gust::ForcingOptions fopt;
    gust::ColumnCanopy canopy;
    Column ()
    {
        in.mesh.n_cell = amrex::IntVect(4, 4, NZ);
        in.mesh.prob_lo = {0, 0, 0};
        in.mesh.prob_hi = {256, 256, NZ};
        in.mesh.is_periodic = {1, 1, 0};
        in.mesh.max_grid_size = {4, 4};
        in.terrain.type = "flat";
        io.first_guess = InitOptions::FirstGuess::column;
        io.column_tol = amrex::Real(1e-10);
        so.stage = gust::StageKind::simplec;
        bin.face = {FaceType::periodic, FaceType::periodic, FaceType::surface_layer, FaceType::periodic,
                    FaceType::periodic, FaceType::slip_wall};
        bin.inflow.kind = gust::InflowProfileOptions::Kind::log_law;
        bin.inflow.ustar = amrex::Real(0.4);
        bin.inflow.z0 = amrex::Real(0.1);
        bin.theta_ref.z = {0};
        bin.theta_ref.theta = {300};
        bin.rayleigh_zdamp = 48;
        bin.rayleigh_dampcoef = amrex::Real(0.05);
        bin.rayleigh_fields = {true, true, true, false, true};
        surf.z0 = amrex::Real(0.1);
        fopt.pressure_gradient = {amrex::Real(F), 0, 0};
        fopt.buoyancy.gravity = amrex::Real(9.81);
        fopt.buoyancy.theta0 = 300;
    }
    void forest (double f, bool tke_on = true)
    {
        canopy.f.assign(NZ, 0.0);
        for (int k = 0; k < NCANOPY; ++k) { canopy.f[static_cast<std::size_t>(k)] = f; }
        canopy.tke.on = tke_on;
        canopy.tke.name = tke_on ? "on" : "off";
    }
    gust::ColumnGuessResult run () const
    {
        return gust::run_column_guess(in, io, so, bin, surf, true, copt, fopt, std::string(), canopy);
    }
};

// The k of the deep-canopy balance beta_p U^3 = beta_d U k + Cmu0^3 k^1.5 / alpha, by bisection on [0, beta_p U^2 / beta_d]
// (the right side grows with k; at the upper end it exceeds the left).
double deep_canopy_k (double U, double beta_p, double beta_d, double alpha, double cmu0)
{
    double lo = 0, hi = beta_p * U * U / beta_d;
    for (int it = 0; it < 200; ++it) {
        const double k = 0.5 * (lo + hi);
        const double rhs = beta_d * U * k + std::pow(cmu0, 3) * std::pow(k, 1.5) / alpha;
        if (rhs < beta_p * U * U * U) { lo = k; } else { hi = k; }
    }
    return 0.5 * (lo + hi);
}

} // namespace

TEST(ColumnCanopy, DeepCanopyBalancesAreClosedForms)
{
    Column s;
    s.forest(FDRAG);
    const auto r = s.run();
    ASSERT_TRUE(r.steady);
    const gust::ClosureConstants cc = gust::closure_constants(s.copt);
    const double U = std::sqrt(F / FDRAG);
    const double kdeep = deep_canopy_k(U, s.canopy.tke.beta_p, s.canopy.tke.beta_d, s.canopy.tke.length_alpha,
                                       static_cast<double>(cc.Cmu0));
    ASSERT_GT(kdeep, static_cast<double>(cc.tke_floor));
    // The table's row i is the cell i - 1 (row 0 is the ground): the layers at 10 to 30 m.
    int checked = 0;
    for (std::size_t i = 1; i < r.z.size(); ++i) {
        const double z = static_cast<double>(r.z[i]);
        if (z < 10 || z > 30) { continue; }
        const double speed = std::abs(std::complex<double>(static_cast<double>(r.u[i]), static_cast<double>(r.v[i])));
        EXPECT_NEAR(speed, U, tol(1e-3, 1e-2) * U) << "z " << z;
        EXPECT_EQ(r.v[i], amrex::Real(0)) << "z " << z;   // no rotation
        EXPECT_NEAR(static_cast<double>(r.tke[i]), kdeep, tol(1e-3, 1e-2) * kdeep) << "z " << z;
        ++checked;
    }
    EXPECT_EQ(checked, 20);
    // Above the canopy the wind is the open column's order (the forcing accelerates it): well above U.
    EXPECT_GT(static_cast<double>(r.u.back()), 10 * U);
}

TEST(ColumnCanopy, KTermsOffKeepsTheDragOnly)
{
    Column off;
    off.forest(FDRAG, false);
    const auto r = off.run();
    ASSERT_TRUE(r.steady);
    Column on;
    on.forest(FDRAG, true);
    const auto ron = on.run();
    ASSERT_TRUE(ron.steady);
    const double U = std::sqrt(F / FDRAG);
    int checked = 0;
    for (std::size_t i = 1; i < r.z.size(); ++i) {
        const double z = static_cast<double>(r.z[i]);
        if (z < 10 || z > 20) { continue; }
        // The larger k makes a larger nu_t, so the diffusion from above adds to the forcing: the speed is above
        // sqrt(F/f), by 34 % at 10 m to about 110 % at 20 m (measured; 240 % at 30 m), never below it.
        EXPECT_GT(static_cast<double>(r.u[i]), U) << "z " << z;
        EXPECT_LT(static_cast<double>(r.u[i]), 2.5 * U) << "z " << z;
        EXPECT_GT(static_cast<double>(r.tke[i]), 0.3) << "z " << z;
        EXPECT_GT(static_cast<double>(r.tke[i]), 100 * static_cast<double>(ron.tke[i])) << "z " << z;
        ++checked;
    }
    EXPECT_EQ(checked, 10);
}

TEST(ColumnCanopy, NoCanopyAndZeroCanopyAgreeBitForBit)
{
    Column bare;
    const auto a = bare.run();
    Column zero;
    zero.forest(0.0);
    const auto b = zero.run();
    ASSERT_TRUE(a.steady && b.steady);
    ASSERT_EQ(a.u.size(), b.u.size());
    for (std::size_t i = 0; i < a.u.size(); ++i) {
        EXPECT_EQ(a.u[i], b.u[i]) << i;
        EXPECT_EQ(a.v[i], b.v[i]) << i;
        EXPECT_EQ(a.tke[i], b.tke[i]) << i;
    }
    EXPECT_EQ(a.iterations, b.iterations);
    // The canopy slows the column: the dense forest's first cell against the bare column's.
    Column dense;
    dense.forest(FDRAG);
    const auto d = dense.run();
    EXPECT_LT(static_cast<double>(d.u[1]), 0.2 * static_cast<double>(a.u[1]));
}

TEST(ColumnCanopy, ProfileMustBeOnePerLayer)
{
    Column s;
    s.forest(FDRAG);
    s.canopy.f.pop_back();
    EXPECT_THROW(s.run(), std::runtime_error);
}
