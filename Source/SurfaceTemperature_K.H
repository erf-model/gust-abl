#ifndef GUST_SURFACE_TEMPERATURE_K_H_
#define GUST_SURFACE_TEMPERATURE_K_H_

// The surface layer of one column with a given ground temperature theta_s
// (bc.surface_stability = surface_temperature, SurfaceLayer.H): Monin-Obukhov
// similarity from the bulk Richardson number of the reference point, after
// ERF's surface_temp mode (ERF_MOSTStress.H, the surface_temp flux at ERF SHA
// 92410d5), re-implemented here with no copied code. ERF solves on L with
// Jimenez's psi; this solves on zeta with Businger-Dyer's:
//
//   Ri_b = g d dtheta / (theta_0 U^2),  dtheta = theta_r - theta_s,  bounded to [kRibMin, kRibMax],
//   G(zeta) = zeta - Ri_b cm(zeta)^2 / ch(zeta) = 0,
//   cm = max(lz - psi_m(zeta), m),  ch = max(lz + kB - psi_h(zeta), m),  lz = ln((d + z0)/z0),
//   u* = kappa U / cm,  theta* = kappa dtheta / ch,  q = -u* theta*,  1/L = zeta / d.
//
// kB = ln(z0 / z0t) >= 0 is the thermal roughness term (bc.surface_kb_inverse): the heat's roughness length z0t
// = z0 exp(-kB), so ln((d + z0)/z0t) = lz + kB. kB = 0 (the default) is z0t = z0 and the form before the input,
// bit for bit (lz + 0 is exact).
//
// zeta has the sign of Ri_b; the root is bracketed in [-100, 0] (unstable) or
// [0, 100] (stable), the wall model's clamp on d/L, and found by Newton's
// method with a bisection step whenever Newton leaves the bracket. Without a
// sign change in the bracket (|Ri_b| beyond the clamp's reach) zeta is the
// bracket's end. The tolerance is round-off, so u* and q are smooth
// functions of the state (as jvp needs). A Ri_b of zero gives zeta = 0 and the
// log law exactly.
//
// The bound on Ri_b: ERF's -4 below. Above, 0.19: Businger-Dyer's stable
// side, cm = ch = lz + 5 zeta, gives Ri_b = zeta ch / cm^2 = zeta / (lz + 5 zeta) < 0.2 for
// every zeta (the column guess's kRibMaxBusingerDyer, ColumnGuess.cpp), so a
// calm, cooled column keeps a small flux instead of none. The root then lies at a
// zeta of order 10, far past Businger-Dyer's fitted range (zeta up to about 1).
//
// With bc.surface_stable_functions = jimenez (JIMENEZ = true) both sides are ERF's surface_temp
// pair (calc_psi_m2/calc_psi_h2, StabilityFunctions.H), as WRF's revised MM5 surface layer: the
// stable side is Cheng and Brutsaert's (2005), psi grows as a log, Ri_b has no critical value and
// the bound above is ERF's 4 (kSurfaceTemperatureRibMaxJimenez); the unstable side is Jimenez et
// al.'s (2012) blend of Businger-Dyer and a free-convection form. The bracket is unchanged.

#include "StabilityFunctions.H"

#include <AMReX_Algorithm.H>
#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>
#include <limits>

namespace gust {

//! The bounds of the bulk Richardson number (see above).
constexpr amrex::Real kSurfaceTemperatureRibMin = amrex::Real(-4);
constexpr amrex::Real kSurfaceTemperatureRibMax = amrex::Real(0.19);
//! The upper bound with the Cheng-Brutsaert stable side: ERF's (ERF_MOSTStress.H, surface_temp).
constexpr amrex::Real kSurfaceTemperatureRibMaxJimenez = amrex::Real(4);

//! The result of the column's solve.
struct SurfaceTemperatureMost
{
    amrex::Real ustar;   //!< u* (m/s)
    amrex::Real q;       //!< kinematic heat flux, upward > 0 (K m/s)
    amrex::Real c;       //!< kappa u* / ch (m/s): q = c (theta_s - theta_r)
    amrex::Real zeta;    //!< d / L
    bool converged;
};

//! The column's surface layer: kappa, U (m/s, already floored), d (m), lz = ln((d + z0)/z0), den_min the floor m,
//! goth = g / theta_0 (m/s^2/K), dtheta = theta_r - theta_s (K), rib_min and rib_max the bounds of Ri_b,
//! kb = ln(z0 / z0t) >= 0. JIMENEZ: the Cheng-Brutsaert stable side (false: Businger-Dyer).
template <bool JIMENEZ = false>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SurfaceTemperatureMost surface_temperature_most (amrex::Real kappa, amrex::Real U, amrex::Real d, amrex::Real lz,
                                                 amrex::Real den_min, amrex::Real goth, amrex::Real dtheta,
                                                 amrex::Real rib_min, amrex::Real rib_max,
                                                 amrex::Real kb = 0) noexcept
{
    const StableCoefficients sc{};   // Businger-Dyer (JIMENEZ = false)
    auto psi_m = [&] (amrex::Real z) noexcept {
        if constexpr (JIMENEZ) { return most_psi_m_jimenez(z); } else { return most_psi_m<false>(z, sc); }
    };
    auto psi_h = [&] (amrex::Real z) noexcept {
        if constexpr (JIMENEZ) { return most_psi_h_jimenez(z); } else { return most_psi_h<false>(z, sc); }
    };
    const amrex::Real rib = amrex::max(rib_min, amrex::min(rib_max, goth * d * dtheta / (U * U)));
    auto cm_of = [&] (amrex::Real z) noexcept { return amrex::max(lz - psi_m(z), den_min); };
    const amrex::Real lzh = lz + kb;
    auto ch_of = [&] (amrex::Real z) noexcept { return amrex::max(lzh - psi_h(z), den_min); };
    auto G = [&] (amrex::Real z) noexcept {
        const amrex::Real cm = cm_of(z);
        return z - rib * cm * cm / ch_of(z);
    };
    const amrex::Real eps = std::numeric_limits<amrex::Real>::epsilon();
    amrex::Real zeta = 0;
    bool converged = true;
    if (rib != 0) {
        // G(0) = -Ri_b lz^2/lzh has the sign of -Ri_b; the root lies on Ri_b's side of 0.
        amrex::Real lo = (rib < 0) ? amrex::Real(-100) : amrex::Real(0);
        amrex::Real hi = (rib < 0) ? amrex::Real(0) : amrex::Real(100);
        const amrex::Real glo = G(lo), ghi = G(hi);
        if (rib > 0 && !(ghi > 0)) {
            zeta = hi;
        } else if (rib < 0 && !(glo < 0)) {
            zeta = lo;
        } else {
            converged = false;
            // From the first-order guess zeta = Ri_b lz (cm = ch = lz), inside the bracket.
            zeta = amrex::max(lo, amrex::min(hi, rib * amrex::max(lz, den_min)));
            const amrex::Real h = std::sqrt(eps);
            for (int it = 0; it < 200 && !converged; ++it) {
                const amrex::Real g = G(zeta);
                if (g < 0) { lo = zeta; } else { hi = zeta; }
                const amrex::Real dz = h * amrex::max(std::abs(zeta), amrex::Real(1));
                const amrex::Real gp = (G(zeta + dz) - g) / dz;
                amrex::Real next = (gp > 0) ? zeta - g / gp : amrex::Real(0.5) * (lo + hi);
                if (!(next > lo && next < hi)) { next = amrex::Real(0.5) * (lo + hi); }
                converged = std::abs(next - zeta) <= 64 * eps * amrex::max(std::abs(next), amrex::Real(1))
                            || !(hi - lo > 64 * eps * amrex::max(std::abs(zeta), amrex::Real(1)));
                zeta = next;
            }
        }
    }
    const amrex::Real cm = cm_of(zeta), ch = ch_of(zeta);
    const amrex::Real ustar = kappa * U / cm;
    const amrex::Real c = kappa * ustar / ch;
    return {ustar, -c * dtheta, c, zeta, converged};
}

} // namespace gust

#endif
