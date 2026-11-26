#ifndef GUST_SURFACE_WSTAR_K_H_
#define GUST_SURFACE_WSTAR_K_H_

// The convective velocity scale of the surface layer (bc.surface_wstar = beljaars, SurfaceLayer.H): Beljaars's
// (1995) gustiness, as ERF's surface flux and surface_temp modes (ERF_MOSTStress.H with ERF_Wstar.H's calc_wstar,
// ERF SHA 151d6e0), re-implemented here with no copied code:
//
//   w* = (g / theta_0  z_i  q)^(1/3) for q > 0, else 0,
//   U_eff = sqrt(U^2 + (beta w*)^2),  beta = kBeljaarsBeta = 1.2,
//
// with q the kinematic heat flux (K m/s, upward > 0; theta_v's in a moist run), z_i the boundary layer's depth
// (m) and U the wind speed at the reference point. U_eff replaces U in the column's Monin-Obukhov solve and
// in its drag: u* = kappa U_eff / cm, tau = u*^2 / U_eff u_t. In a calm, heated column the eddies of the
// convective boundary layer keep a flux that the mean wind alone would not carry.
//
// ERF's g / theta_v is the cell's; GUST's goth = g / theta_0 is the reference value every buoyancy term uses.
// With a given flux (heat_flux) U_eff is direct. With a given ground temperature (surface_temperature) q
// depends on U_eff: ERF lags q by one step; GUST solves the column's fixed point
//   U_eff = F(U_eff) = sqrt(U^2 + (beta w*(q(U_eff)))^2)
// to round-off (surface_temperature_most_wstar below), so u* and q stay functions of the state.
// Where q grows at most linearly in U_eff, dF/dU_eff <= beta^2 w*^2 / (3 F U_eff) < 1/3 and the iteration
// contracts. The floor of MOST's denominator (den_min) puts a kink in q(U_eff), and past it q can grow faster:
// with z0 = 2, d = 5, kB^-1 = 2, z_i = 10, U = 0.3 and dtheta = -8 the iteration cycles between two values. So when
// it has not converged, bisection on F(x) - x finds the root: F(U) >= U, and F(x) < x for x large (w* grows as
// q^(1/3)).

#include "SurfaceTemperature_K.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace gust {

//! Beljaars's (1995) coefficient of w* in U_eff (ERF's Bjr_beta).
constexpr amrex::Real kBeljaarsBeta = amrex::Real(1.2);
//! The range of bc.surface_pblh and of bc.surface_pblh_file's values (m).
constexpr double kPblhMin = 10;
constexpr double kPblhMax = 10000;

//! w* (m/s) from goth = g / theta_0 (m/s^2/K), z_i (m) and q (K m/s): zero for q <= 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real convective_velocity (amrex::Real goth, amrex::Real zi, amrex::Real q) noexcept
{
    return (q > 0) ? std::cbrt(goth * zi * q) : amrex::Real(0);
}

//! U_eff = sqrt(U^2 + (kBeljaarsBeta w*)^2) (m/s).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real beljaars_speed (amrex::Real U, amrex::Real goth, amrex::Real zi, amrex::Real q) noexcept
{
    const amrex::Real bw = kBeljaarsBeta * convective_velocity(goth, zi, q);
    return std::sqrt(U * U + bw * bw);
}

//! surface_temperature_most with U_eff: the fixed point U_eff = beljaars_speed(U, goth, zi, q(U_eff)) from
//! U_eff = U (bisection if that does not converge in 100 steps), then the solve at U_eff (its u*, q, c and zeta).
//! zi is z_i (m); the other arguments are surface_temperature_most's. ueff returns U_eff; converged is false if
//! either the zeta solve or the search for U_eff did not converge.
template <bool JIMENEZ = false>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SurfaceTemperatureMost surface_temperature_most_wstar (amrex::Real kappa, amrex::Real U, amrex::Real d,
                                                       amrex::Real lz, amrex::Real den_min, amrex::Real goth,
                                                       amrex::Real dtheta, amrex::Real rib_min,
                                                       amrex::Real rib_max, amrex::Real kb, amrex::Real zi,
                                                       amrex::Real& ueff) noexcept
{
    const amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real u_eff = U;
    SurfaceTemperatureMost r = surface_temperature_most<JIMENEZ>(kappa, u_eff, d, lz, den_min, goth, dtheta, rib_min,
                                                                 rib_max, kb);
    bool fixed = (r.q <= 0);   // no w*: U_eff = U
    // The zeta solve stops at a step of 64 eps max(|zeta|, 1) (SurfaceTemperature_K.H), so F carries a noise of a
    // few tens of eps: a tighter test can cycle between two values for ever. Stop at 256 eps, or once the step at
    // round-off level (below 4096 eps) no longer shrinks; with dF/dU_eff < 1/3 the error is then of that order.
    amrex::Real step_old = std::numeric_limits<amrex::Real>::max();
    for (int it = 0; it < 100 && !fixed; ++it) {
        const amrex::Real next = beljaars_speed(U, goth, zi, r.q);
        const amrex::Real step = std::abs(next - u_eff);
        fixed = (step <= 256 * eps * next) || (step <= 4096 * eps * next && step >= step_old);
        step_old = step;
        u_eff = next;
        r = surface_temperature_most<JIMENEZ>(kappa, u_eff, d, lz, den_min, goth, dtheta, rib_min, rib_max, kb);
    }
    if (!fixed) {
        // Bisection on F(x) - x over [lo, hi] with F(lo) >= lo and F(hi) < hi (doubling hi from U_eff).
        amrex::Real lo = U;
        amrex::Real hi = amrex::max(u_eff, U);
        bool bracket = false;
        for (int it = 0; it < 64 && !bracket; ++it) {
            const SurfaceTemperatureMost rh = surface_temperature_most<JIMENEZ>(kappa, hi, d, lz, den_min, goth, dtheta,
                                                                                rib_min, rib_max, kb);
            bracket = (beljaars_speed(U, goth, zi, rh.q) < hi);
            if (!bracket) { lo = hi; hi *= 2; }
        }
        for (int it = 0; it < 200 && bracket && hi - lo > 256 * eps * hi; ++it) {
            const amrex::Real mid = amrex::Real(0.5) * (lo + hi);
            const SurfaceTemperatureMost rm = surface_temperature_most<JIMENEZ>(kappa, mid, d, lz, den_min, goth, dtheta,
                                                                                rib_min, rib_max, kb);
            if (beljaars_speed(U, goth, zi, rm.q) >= mid) { lo = mid; } else { hi = mid; }
        }
        fixed = bracket && (hi - lo <= 256 * eps * hi);
        u_eff = amrex::Real(0.5) * (lo + hi);
        r = surface_temperature_most<JIMENEZ>(kappa, u_eff, d, lz, den_min, goth, dtheta, rib_min, rib_max, kb);
    }
    r.converged = r.converged && fixed;
    ueff = u_eff;
    return r;
}

} // namespace gust

#endif
