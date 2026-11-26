#ifndef GUST_STABILITYFUNCTIONS_H_
#define GUST_STABILITYFUNCTIONS_H_

// MOST's stability functions of zeta = z/L, selected by bc.surface_stable_functions.
//
// Unstable side (zeta < 0): Businger-Dyer, ERF's similarity_funs (extern_erf/ERF_MOSTUtils.H), except with jimenez.
// Stable side (zeta >= 0), one form with two coefficients per function:
//
//   phi(zeta) = (1 + a zeta) / (1 + b zeta),        psi(zeta) = -((a - b) / b) ln(1 + b zeta)   (b > 0),
//
// - businger_dyer (default): a = 5, b = 0, so phi = 1 + 5 zeta and psi = -5 zeta (the b -> 0 limit), ERF's
//   calc_phi_m/calc_psi_m (Dyer 1974);
// - qnse: the quasi-normal scale elimination functions (Sukoriansky, Galperin and Perov 2005), as ERF's
//   enable_qnse_stable_functions (ERF_ComputeDiffusivityYSUNew.cpp): a_m, b_m, a_h, b_h from
//   bc.surface_qnse_am, _bm, _ah, _bh, ERF's defaults 2.5, 0.2, 2.5, 0.2. phi tends to a/b = 12.5 for large
//   zeta (the turbulence is not cut off), where Businger-Dyer grows without bound.
// - jimenez (bc.surface_stability = surface_temperature only): ERF's calc_psi_m2/calc_psi_h2 (Jimenez et al. 2012,
//   WRF's revised MM5 surface layer). Stable side Cheng and Brutsaert (2005), fitted to very stable data:
//
//     psi(zeta) = -a ln(zeta + (1 + zeta^b)^(1/b)),
//     phi(zeta) = 1 - zeta psi'(zeta) = 1 + a (zeta + zeta^b (1 + zeta^b)^(1/b - 1)) / (zeta + (1 + zeta^b)^(1/b)),
//
//   a = 6.1, b = 2.5 (momentum) and a = 5.3, b = 1.1 (heat). psi grows as a log, so Ri_b has no critical value
//   (Businger-Dyer's is 0.2).
//
//   Its unstable side is the rest of the pair: psi = (psi_k + zeta^2 psi_c) / (1 + zeta^2), the Businger-Dyer
//   (Kansas) psi_k blended with a free-convection psi_c (Grachev, Fairall and Bradley 2000), y = (1 - c zeta)^(1/3)
//   with c = 10 (momentum) and 34 (heat). Near zeta = 0 it is Businger-Dyer; for large -zeta, the convective
//   form. most_psi_m_jimenez and the others below; this whole pair is jimenez.
//
// psi is the integral of (1 - phi(x))/x from 0 to zeta. The forms are branch-safe as most_phi's
// (SurfaceLayer.H): each arm of a select only sees an argument of its own sign, |zeta| forms rather than
// min/max, so neither the log nor ERF's square roots meet a negative argument in an arm that is not taken.

#include "ERF_MOSTUtils.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <cmath>

namespace gust {

//! The stable-side coefficients: phi_m = (1 + am zeta)/(1 + bm zeta), phi_h = (1 + ah zeta)/(1 + bh zeta).
//! The defaults are Businger-Dyer's (b = 0).
struct StableCoefficients
{
    amrex::Real am{5}, bm{0}, ah{5}, bh{0};
};

//! phi_m and phi_h of zeta with the given stable side (Businger-Dyer unstable side).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void most_phi (amrex::Real zeta, const StableCoefficients& c, amrex::Real& phi_m, amrex::Real& phi_h) noexcept
{
    const amrex::Real one = 1, gamma = 16;   // Dyer 1974, ERF's value
    const amrex::Real az = std::abs(zeta);
    const amrex::Real zs = amrex::Real(0.5) * (az + zeta);   // max(zeta, 0)
    const amrex::Real zu = amrex::Real(0.5) * (zeta - az);   // min(zeta, 0)
    const amrex::Real r = std::sqrt(one - gamma * zu);       // >= 1
    phi_m = (one + c.am * zs) / (one + c.bm * zs) + (one / std::sqrt(r) - one);
    phi_h = (one + c.ah * zs) / (one + c.bh * zs) + (one / r - one);
}

//! The stable psi of one function, -((a - b)/b) ln(1 + b zeta_s), zeta_s >= 0 and b > 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real qnse_psi (amrex::Real zeta_s, amrex::Real a, amrex::Real b) noexcept
{
    return -((a - b) / b) * std::log(amrex::Real(1) + b * zeta_s);
}

//! psi_m of zeta: QNSE's stable side (QNSE = true) or ERF's Businger-Dyer calc_psi_m (QNSE = false, the same values).
template <bool QNSE>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real most_psi_m (amrex::Real zeta, const StableCoefficients& c) noexcept
{
    const similarity_funs sf;
    if constexpr (QNSE) {
        const amrex::Real az = std::abs(zeta);
        const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
        return (zeta > 0) ? qnse_psi(zs, c.am, c.bm) : sf.calc_psi_m(zu);
    } else {
        // ERF's calc_psi_m, its arms fed by |zeta| forms: -5 zeta stable, the unstable form at min(zeta, 0).
        amrex::ignore_unused(c);
        const amrex::Real az = std::abs(zeta);
        const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
        return (zeta > 0) ? -amrex::Real(5) * zs : sf.calc_psi_m(zu);   // ERF's beta_m = 5 (private there)
    }
}

//! psi_h of zeta, as most_psi_m.
template <bool QNSE>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real most_psi_h (amrex::Real zeta, const StableCoefficients& c) noexcept
{
    const similarity_funs sf;
    if constexpr (QNSE) {
        const amrex::Real az = std::abs(zeta);
        const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
        return (zeta > 0) ? qnse_psi(zs, c.ah, c.bh) : sf.calc_psi_h(zu);
    } else {
        // ERF's calc_psi_h, its arms fed by |zeta| forms: -5 zeta stable, the unstable form at min(zeta, 0).
        amrex::ignore_unused(c);
        const amrex::Real az = std::abs(zeta);
        const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
        return (zeta > 0) ? -amrex::Real(5) * zs : sf.calc_psi_h(zu);   // ERF's beta_h = 5
    }
}

//! Cheng and Brutsaert's (2005) stable psi of one function, zeta_s >= 0 (see above).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cb_psi (amrex::Real zeta_s, amrex::Real a, amrex::Real b) noexcept
{
    return -a * std::log(zeta_s + std::pow(amrex::Real(1) + std::pow(zeta_s, b), amrex::Real(1) / b));
}

//! Cheng and Brutsaert's stable phi of one function, 1 - zeta psi'(zeta), zeta_s >= 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real cb_phi (amrex::Real zeta_s, amrex::Real a, amrex::Real b) noexcept
{
    const amrex::Real one = 1;
    const amrex::Real zb = std::pow(zeta_s, b);
    const amrex::Real r = std::pow(one + zb, one / b);   // (1 + zeta^b)^(1/b)
    return one + a * (zeta_s + zb * r / (one + zb)) / (zeta_s + r);
}

//! The coefficients of Cheng and Brutsaert (2005): psi_m (6.1, 2.5), psi_h (5.3, 1.1).
struct ChengBrutsaert
{
    static constexpr amrex::Real am = amrex::Real(6.1), bm = amrex::Real(2.5);
    static constexpr amrex::Real ah = amrex::Real(5.3), bh = amrex::Real(1.1);
};

//! The convective psi of the Jimenez unstable blend (Grachev, Fairall and Bradley 2000; ERF's psi_c_m, psi_c_h), y =
//! (1 - c zeta)^(1/3) >= 1, and its phi = 1 / y.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real jimenez_psi_c (amrex::Real y) noexcept
{
    const amrex::Real one = 1, two = 2, three = 3, s3 = std::sqrt(three);
    return amrex::Real(1.5) * std::log((y * y + y + one) / three) - s3 * std::atan((two * y + one) / s3)
           + amrex::Real(3.14159265358979323846) / s3;
}

//! The unstable phi of the blend psi = (psi_k + zeta^2 psi_c) / (1 + zeta^2), zeta <= 0: phi = 1 - zeta psi' =
//! (phi_k + zeta^2 phi_c) / (1 + zeta^2) - 2 zeta^2 (psi_c - psi_k) / (1 + zeta^2)^2 (from zeta psi_k' = 1 - phi_k,
//! zeta psi_c' = 1 - phi_c).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real jimenez_phi_unstable (amrex::Real zu, amrex::Real phi_k, amrex::Real psi_k, amrex::Real y) noexcept
{
    const amrex::Real one = 1, z2 = zu * zu, w = one + z2;
    return (phi_k + z2 / y) / w - amrex::Real(2) * z2 * (jimenez_psi_c(y) - psi_k) / (w * w);
}

//! psi_m of zeta, the Jimenez et al. (2012) pair (ERF's calc_psi_m2): Cheng-Brutsaert stable side, the blend of
//! Businger-Dyer and the convective form unstable side.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real most_psi_m_jimenez (amrex::Real zeta) noexcept
{
    const similarity_funs sf;
    const amrex::Real az = std::abs(zeta);
    const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
    return (zeta > 0) ? cb_psi(zs, ChengBrutsaert::am, ChengBrutsaert::bm) : sf.calc_psi_m2(zu);
}

//! psi_h of zeta, as most_psi_m_jimenez (ERF's calc_psi_h2).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real most_psi_h_jimenez (amrex::Real zeta) noexcept
{
    const similarity_funs sf;
    const amrex::Real az = std::abs(zeta);
    const amrex::Real zs = amrex::Real(0.5) * (az + zeta), zu = amrex::Real(0.5) * (zeta - az);
    return (zeta > 0) ? cb_psi(zs, ChengBrutsaert::ah, ChengBrutsaert::bh) : sf.calc_psi_h2(zu);
}

//! phi_m and phi_h of zeta, the companions of most_psi_m_jimenez and most_psi_h_jimenez (1 - zeta psi').
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void most_phi_jimenez (amrex::Real zeta, amrex::Real& phi_m, amrex::Real& phi_h) noexcept
{
    const similarity_funs sf;
    const amrex::Real one = 1, third = one / amrex::Real(3);
    const amrex::Real az = std::abs(zeta);
    const amrex::Real zs = amrex::Real(0.5) * (az + zeta);
    const amrex::Real zu = amrex::Real(0.5) * (zeta - az);
    // Unstable arms at zu <= 0 (at zu = 0 both give exactly 1: psi_k = psi_c = 0 is not needed, z2 = 0).
    const amrex::Real r = std::sqrt(one - 16 * zu);                      // >= 1
    const amrex::Real um = jimenez_phi_unstable(zu, one / std::sqrt(r), sf.calc_psi_m(zu),
                                                std::pow(one - 10 * zu, third));
    const amrex::Real uh = jimenez_phi_unstable(zu, one / r, sf.calc_psi_h(zu), std::pow(one - 34 * zu, third));
    // The stable arm gives exactly 1 at zs = 0, the unstable arm exactly 1 at zu = 0.
    phi_m = cb_phi(zs, ChengBrutsaert::am, ChengBrutsaert::bm) + (um - one);
    phi_h = cb_phi(zs, ChengBrutsaert::ah, ChengBrutsaert::bh) + (uh - one);
}

} // namespace gust

#endif
