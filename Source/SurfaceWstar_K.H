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
// The map is a contraction: q grows at most linearly in U_eff, so dF/dU_eff <= beta^2 w*^2 / (3 F U_eff) < 1/3.

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
//! U_eff = U, then the solve at U_eff (its u*, q, c and zeta). zi is z_i (m); the other arguments are
//! surface_temperature_most's. ueff returns U_eff; converged is false if either iteration did not converge.
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
    for (int it = 0; it < 100 && !fixed; ++it) {
        const amrex::Real next = beljaars_speed(U, goth, zi, r.q);
        fixed = std::abs(next - u_eff) <= 16 * eps * next;
        u_eff = next;
        r = surface_temperature_most<JIMENEZ>(kappa, u_eff, d, lz, den_min, goth, dtheta, rib_min, rib_max, kb);
    }
    r.converged = r.converged && fixed;
    ueff = u_eff;
    return r;
}

} // namespace gust

#endif
