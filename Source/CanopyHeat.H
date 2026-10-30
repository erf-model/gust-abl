#ifndef GUST_CANOPY_HEAT_H_
#define GUST_CANOPY_HEAT_H_

// The leaves' sensible heat exchange with the air (canopy.heat = fixed,
// docs/canopy.rst), the steady form of ERF's fixed-leaf-temperature canopy
// biophysics (ERF_CanopyBiophysics.H, AddCanopyBiophysicsHeatSources, ERF
// 1b00938, erf-model/ERF#3897, from the ICLASS canopy model):
//
//   S_theta = lambda (theta_l - theta),   lambda = 2 a / r_b,
//   r_b = C_b sqrt(l_w / max(|U|, U_min)),
//
// with a the leaf area density (1/m, Canopy::lad), theta_l = canopy.leaf_theta
// the leaves' potential temperature, |U| the cell-centred wind (the faces'
// means, as ERF and CanopyTke.H), r_b (s/m) the leaf boundary-layer
// resistance, and ERF's constants C_b = 100 s^(1/2)/m, l_w = 0.05 m (the leaf
// width), U_min = 0.01 m/s. The 2 counts both sides of a leaf. ERF adds
// 2 rho c_p (T_l - T) a / r_b / (c_p Pi) to rho theta, T = Pi theta: the
// density, c_p and the Exner function Pi cancel. ERF's transpiration (with
// moisture) is not here: GUST carries no vapour.
//
// The residual R_theta gains -S_theta = lambda (theta - theta_l). Implicit
// side: lambda does not depend on theta, so the theta rows' a_P gains
// lambda, the exact derivative; its dependence on the wind is lagged, as the
// drag's (CanopyRows, CanopyDrag.H). lambda >= 0, so the diagonal stays
// positive, and the source pulls theta towards theta_l: with no other
// source, theta stays between theta_l and the boundary values.
//
// simplec.theta = frozen (theta is data, no theta equation):
//   skip   the leaf heat is left out (a warning at start-up);
//   blend  in the canopy, theta is the steady point balance of the leaf
//          exchange with turbulent mixing towards the frozen profile theta_f,
//          the mixing a relaxation at the eddy turnover rate eps / k:
//            0 = lambda (theta_l - theta) + (eps / k) (theta_f - theta),
//            theta = theta_f + w (theta_l - theta_f),
//            w = Da / (1 + Da) = lambda k / (eps + lambda k),   Da = lambda k / eps,
//          from the closure's k (at its floor where below it, as eps takes it)
//          and eps at the iteration's state, set once per SIMPLEC iteration;
//          eps without the canopy's short-circuit sink D_c (that is the wake
//          turbulence's cascade, not the air's renewal). w = 0 where a = 0 (theta_f exactly), w -> 1 in
//          dense foliage at low turbulence, and 0 <= w < 1: no new input, and
//          theta moves smoothly from theta_f to theta_l with the foliage. The
//          rate eps / k is a modelling choice: in a canopy it lies between the
//          flushing rate u*/h (about 4 times slower) and nu_t / L^2 (about 3
//          times faster), so w is uncertain by about that factor in Da.
//          Needs the k–L closure.

#include "Canopy.H"
#include "CanopyDrag.H"
#include "State.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>

#include <cmath>

namespace gust {

//! ERF's constants of the leaf boundary-layer resistance (above).
inline constexpr double kLeafWidth = 0.05;           //!< l_w (m)
inline constexpr double kLeafMinWind = 0.01;         //!< U_min (m/s)
inline constexpr double kLeafBoundaryCoeff = 100.0;  //!< C_b (s^(1/2) / m)

//! lambda = 2 a / r_b (1/s) at leaf area density a (1/m) and wind speed |U| (m/s).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real leaf_exchange_rate (amrex::Real a, amrex::Real wind) noexcept
{
    const amrex::Real lw = amrex::Real(kLeafWidth), umin = amrex::Real(kLeafMinWind);
    const amrex::Real rb = amrex::Real(kLeafBoundaryCoeff) * std::sqrt(lw / (wind > umin ? wind : umin));
    return amrex::Real(2) * a / rb;
}

//! The blend weight w = lambda k / (eps + lambda k) (above); 0 where both vanish.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real frozen_blend_weight (amrex::Real lambda, amrex::Real k, amrex::Real eps) noexcept
{
    const amrex::Real lk = lambda * k;
    const amrex::Real den = eps + lk;
    return den > 0 ? lk / den : amrex::Real(0);
}

//! The level's canopy heat: the leaf area density and the options, from the
//! level's Canopy.
class CanopyHeat
{
public:
    //! From canopy.lad() and canopy.heat(). Aborts unless canopy.heat().on().
    void build (const Canopy& canopy);
    [[nodiscard]] const amrex::MultiFab& lad () const { return *m_lad; }
    [[nodiscard]] amrex::Real leaf_theta () const { return m_leaf_theta; }
    [[nodiscard]] CanopyFrozenHeat frozen () const { return m_frozen; }

private:
    const amrex::MultiFab* m_lad{nullptr};
    amrex::Real m_leaf_theta{0};
    CanopyFrozenHeat m_frozen{CanopyFrozenHeat::skip};
};

//! R_theta += lambda (theta - theta_l) on the valid cells (R on the cells,
//! no ghosts, on the leaf area density's layout); u, v, w of each cell's own
//! faces (no ghosts read). Aborts on another layout.
void add_canopy_heat (const CanopyHeat& heat, const State& q, amrex::MultiFab& R_theta);

//! rows.theta = lambda at (u, v, w) on the cells (defined here if empty): the
//! theta rows' diagonal at the solver's linearisation state.
void canopy_heat_rows (const CanopyHeat& heat, const amrex::MultiFab& u, const amrex::MultiFab& v,
                       const amrex::MultiFab& w, CanopyRows& rows);

//! canopy.heat_frozen = blend: theta = theta_f + w (theta_l - theta_f) on the
//! valid cells, from the frozen field theta_f (on q's cells, no ghosts read),
//! q's velocity and k (taken at tke_floor where below it, as the closure takes
//! it for eps), and eps (the closure's, on q's cells). Returns the largest
//! |change| of theta over the level (collective).
amrex::Real blend_frozen_theta (const CanopyHeat& heat, const amrex::MultiFab& theta_f, const amrex::MultiFab& eps,
                                amrex::Real tke_floor, State& q);

} // namespace gust

#endif
