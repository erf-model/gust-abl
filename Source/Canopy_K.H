#ifndef GUST_CANOPY_K_H_
#define GUST_CANOPY_K_H_

// Point relations of the canopy's leaf area density (Canopy.H).
//
// A canopy of height h and leaf area index LAI has the leaf area density
//   a(z) = (LAI / h) P(z / h) / N,   N = int_0^1 P(eta) d eta,   0 <= z < h,
// and a = 0 above h, so that int_0^h a dz = LAI. Two shapes P (ERF's
// erf.forest_tree_type, ERF_ForestDrag.cpp):
//   uniform:            P = 1, N = 1;
//   Lalic–Mihailovic:   r = (1 - eta_m) / (1 - eta),
//                       P = r^6 exp(6 (1 - r))        eta <  eta_m,
//                       P = r^(1/2) exp((1 - r) / 2)  eta >= eta_m,
//   eta_m = laimax, the height of the largest density as a fraction of h
//   (Lalic and Mihailovic 2004, doi:10.1175/1520-0450(2004)043<0641:AERDLD>2.0.CO;2).
//
// The drag field of a cell is the mean of Cd a(z) over the cell's height
// above the ground, z_b to z_t:
//   f = Cd LAI (F(z_t / h) - F(z_b / h)) / (z_t - z_b),
// with F the profile's cumulative fraction, F(eta) = int_0^eta P / N (F = 1
// above the top). F is in closed form. The difference is formed from the
// branch integrals themselves (lm_fraction), not as F - F, so a cell just
// below the top, where F is within rounding of 1, keeps its small share.
// With t = r(eta) and c = 1 - eta_m,
// d eta = c dt / t^2, so
//   lower branch: int_0^eta P = c [E(c) - E(t)],
//                 E(t) = exp(6 (1 - t)) (t^4/6 + 4 t^3/6^2 + 12 t^2/6^3 + 24 t/6^4 + 24/6^5),
//   upper branch: int_eta^1 P = c [2 t^(-1/2) exp((1 - t)/2) - sqrt(2 pi e) erfc(sqrt(t / 2))],
// and N is the sum of the two at eta = eta_m (t = 1). ERF samples a(z) at the
// cell centre with N from a 100-point sum instead (ERF's N is within 0.8 % of
// the exact value for eta_m <= 0.97); the centre value's column integral
// departs from Cd LAI by 10-50 % once the layers are h/4 or thicker.
// Everything is evaluated in double precision in every build.

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>

#include <cmath>

namespace gust {

//! The constants of one Lalic–Mihailovic profile (make_lm_profile).
struct LmProfile
{
    double eta_m{0.8};   //!< laimax
    double c{0.2};       //!< 1 - eta_m
    double e_c{0};       //!< E(c), the lower branch's value at the ground
    double norm{1};      //!< N
};

//! E(t) of the lower branch (see above).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double lm_lower_e (double t) noexcept
{
    const double poly = (((t / 6.0 + 4.0 / 36.0) * t + 12.0 / 216.0) * t + 24.0 / 1296.0) * t + 24.0 / 7776.0;
    return std::exp(6.0 * (1.0 - t)) * poly;
}

//! int_eta^1 P of the upper branch at t = r(eta) >= 1 (see above); 0 at t = inf.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double lm_upper_tail (double t, double c) noexcept
{
    constexpr double sqrt_2pie = 4.1327313541224930;   // sqrt(2 pi e)
    return c * (2.0 / std::sqrt(t) * std::exp(0.5 * (1.0 - t)) - sqrt_2pie * std::erfc(std::sqrt(0.5 * t)));
}

//! The profile constants for eta_m in [0, 1) (checked by the caller).
inline LmProfile make_lm_profile (double eta_m)
{
    LmProfile p;
    p.eta_m = eta_m;
    p.c = 1.0 - eta_m;
    p.e_c = lm_lower_e(p.c);
    p.norm = p.c * (p.e_c - lm_lower_e(1.0)) + lm_upper_tail(1.0, p.c);
    return p;
}

//! The shape P(eta) itself, 0 <= eta < 1 (for the tests).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double lm_shape (double eta, const LmProfile& p) noexcept
{
    const double r = p.c / (1.0 - eta);
    return eta < p.eta_m ? std::pow(r, 6) * std::exp(6.0 * (1.0 - r)) : std::sqrt(r) * std::exp(0.5 * (1.0 - r));
}

//! F(eta) of the uniform profile: eta clipped to [0, 1].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double uniform_cdf (double eta) noexcept
{
    return eta <= 0.0 ? 0.0 : (eta >= 1.0 ? 1.0 : eta);
}

//! F(b) - F(a) of the uniform profile (a <= b).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double uniform_fraction (double a, double b) noexcept
{
    return uniform_cdf(b) - uniform_cdf(a);
}

//! F(eta) of the Lalic–Mihailovic profile: 0 at and below the ground, 1 at
//! and above the top.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double lm_cdf (double eta, const LmProfile& p) noexcept
{
    if (eta <= 0.0) { return 0.0; }
    if (eta >= 1.0) { return 1.0; }
    const double t = p.c / (1.0 - eta);
    if (eta <= p.eta_m) { return p.c * (p.e_c - lm_lower_e(t)) / p.norm; }
    return 1.0 - lm_upper_tail(t, p.c) / p.norm;
}

//! int_eta^1 P for eta_m <= eta <= 1 (0 at the top).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double lm_tail_at (double eta, const LmProfile& p) noexcept
{
    return eta >= 1.0 ? 0.0 : lm_upper_tail(p.c / (1.0 - eta), p.c);
}

//! F(b) - F(a) of the Lalic–Mihailovic profile (a <= b), from the branch
//! integrals: the lower branch's E at both ends, the upper branch's tails at
//! both ends, or each branch's part when [a, b] spans eta_m.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double lm_fraction (double a, double b, const LmProfile& p) noexcept
{
    a = a <= 0.0 ? 0.0 : (a >= 1.0 ? 1.0 : a);
    b = b <= 0.0 ? 0.0 : (b >= 1.0 ? 1.0 : b);
    if (!(a < b)) { return 0.0; }
    if (b <= p.eta_m) { return p.c * (lm_lower_e(p.c / (1.0 - a)) - lm_lower_e(p.c / (1.0 - b))) / p.norm; }
    if (a >= p.eta_m) { return (lm_tail_at(a, p) - lm_tail_at(b, p)) / p.norm; }
    return (p.c * (lm_lower_e(p.c / (1.0 - a)) - lm_lower_e(1.0)) + lm_upper_tail(1.0, p.c) - lm_tail_at(b, p)) / p.norm;
}

//! The cell mean of Cd a(z) between z_b and z_t above the ground (z_t > z_b)
//! for a canopy of height h, leaf area index lai and drag coefficient cd;
//! 0 where h <= 0 or lai <= 0. fraction(a, b) is F(b) - F(a): uniform_fraction
//! or a profile's lm_fraction.
template <typename Fraction>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double cell_drag (double z_b, double z_t, double h, double lai, double cd, const Fraction& fraction) noexcept
{
    if (!(h > 0.0) || !(lai > 0.0) || !(z_b < h)) { return 0.0; }
    return cd * lai * fraction(z_b / h, z_t / h) / (z_t - z_b);
}

} // namespace gust

#endif
