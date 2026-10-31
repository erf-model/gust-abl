#ifndef GUST_BUILDING_FORCING_H_
#define GUST_BUILDING_FORCING_H_

// How buildings act on momentum (Buildings.H has the geometry; docs/buildings.rst):
//
// - a blocked face (fraction 1) is a fixed row, u = 0: R = u and the identity in the
//   implicit rows, no 1/dtau, d = 0 in SIMPLEC, the steady limit of ERF's
//   if_implicit_projection, with no drag coefficient;
// - partial faces carry the open part of the mass flux (Mesh::set_open_faces);
// - buildings.wall_form = fraction_stress: a roof or underside takes the form 7 wall
//   law (ERF erf.if_wall_form = fraction_stress, Kynema-SGF wall_form 7), neutral:
//     u*^2 = (1 - beta) u*_own^2 + beta u*_ref^2,  u*_x = kappa U_x / ln(d_x / z0),
//   with the own height floored at 10 z0 (Kynema-SGF 996cdc9e) and u*^2 = u*_own^2 where
//   there is no reference cell (d2 = 0). It is a sink on the horizontal velocity,
//     S = -r u,  r = f u*^2 / (Delta_z |U_own|),
//   f the blocked part of the wall's z-face, Delta_z the layer's height. The
//   viscous flux through the blocked part is closed, in the residual and in the
//   line solve. r stays out of the local pseudo-time step and goes into a_P
//   (its derivative with r frozen; the fixed point does not depend on it).
// - buildings.wall_form = drag: the fixed rows alone (Kynema-SGF's laminar pathway).
//
// Sides of buildings have the fixed rows here and their wall law in PR-40; k and theta
// in PR-41 and PR-42. Horizontal speeds: on a slope the roof's tilt is not taken out.

#include "Buildings.H"
#include "ImplicitOperator.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <cmath>

namespace gust {

//! The wall law of one cell: u*^2 of the own and reference speeds.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real form7_ustar2 (amrex::Real u_own, amrex::Real u_ref, amrex::Real beta, amrex::Real d1, amrex::Real d2,
                          amrex::Real z0, amrex::Real kappa) noexcept
{
    const amrex::Real floor = amrex::Real(10) * z0;
    const amrex::Real own = kappa * u_own / std::log(amrex::max(d1, floor) / z0);
    if (d2 <= amrex::Real(0)) { return own * own; }
    const amrex::Real ref = kappa * u_ref / std::log(amrex::max(d2, floor) / z0);
    return (amrex::Real(1) - beta) * own * own + beta * ref * ref;
}

//! Speeds below this do not divide the rate (m/s): u*^2 / U stays finite at rest.
inline constexpr amrex::Real kBuildingMinSpeed = amrex::Real(1.0e-3);

//! The momentum terms of one level's buildings.
class BuildingForcing
{
public:
    void build (const Buildings& b, const BuildingsOptions& opt, const amrex::Geometry& geom, const amrex::MultiFab& detJ);

    [[nodiscard]] bool wall_law () const { return m_wall_form == BuildingWallForm::fraction_stress; }
    //! The rates r at the cells (one ghost layer) from the last update.
    [[nodiscard]] const amrex::MultiFab& rate () const { return m_rate; }

    //! The rates at the state (u, v on their faces; valid points are read).
    void update (const amrex::MultiFab& u, const amrex::MultiFab& v);
    //! R_u, R_v += r u, r v with the rates of update.
    void add_sink (const amrex::MultiFab& u, const amrex::MultiFab& v, amrex::MultiFab& Ru, amrex::MultiFab& Rv) const;
    //! Gu, Gv: the viscous z-fluxes of u and v, times the open part of each flux's face.
    void close_viscous_fluxes (const amrex::MultiFab& open_z, amrex::MultiFab& Gu, amrex::MultiFab& Gv) const;
    //! R_c = q_c on the blocked faces (last: they replace every other row).
    void apply_rows (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w, amrex::MultiFab& Ru,
                     amrex::MultiFab& Rv, amrex::MultiFab& Rw) const;
    //! Component c's rows: a_P gains r, the wall faces' couplings go, and the blocked
    //! faces become identity rows.
    void add_to_implicit (int c, const amrex::MultiFab& open_z, amrex::MultiFab& A) const;
    //! mask = 1 on the blocked faces of component c as well.
    void add_fixed_rows (int c, amrex::MultiFab& mask) const;

private:
    const Buildings* m_b{nullptr};
    amrex::Geometry m_geom;
    const amrex::MultiFab* m_detJ{nullptr};
    amrex::MultiFab m_rate;
    amrex::Real m_z0{0.1};
    BuildingWallForm m_wall_form{BuildingWallForm::fraction_stress};
};

//! How many fluid cells have no open path to the largest-index region. Aborts are
//! the caller's. A sealed pocket has its own null vector in the pressure correction.
amrex::Long count_sealed_fluid_cells (const Buildings& b, const amrex::Geometry& geom);

} // namespace gust

#endif
