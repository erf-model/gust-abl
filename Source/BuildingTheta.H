#ifndef GUST_BUILDING_THETA_H_
#define GUST_BUILDING_THETA_H_

// How buildings act on the potential temperature (buildings.theta_mode = ambient | fixed;
// Buildings.H has the geometry, BuildingForcing.H the wall law; docs/buildings.rst). The
// form 7 wall treatment with a given surface temperature theta_b (ERF ib_wall::wall_state,
// surface_temp), neutral:
//
// - a solid cell (beta = 1) holds theta at theta_b: R_theta = theta - theta_b, an identity
//   row;
// - each wall with the law takes heat from the air of its cell at the kinematic flux
//     q = u* theta*,  theta* = (1 - beta) kappa (theta - theta_b) / ln(d1 / z0)
//                            +      beta  kappa (theta_ref - theta_b) / ln(d2 / z0),
//   u* the momentum law's of that wall, d1 and d2 its distances with their floor of 10 z0,
//   theta_ref the next cell out along the wall's normal (BuildingForcing.H:
//   form7_heat_velocities). Over the cell,
//     R_theta += sum over the walls of f q / Delta_n,
//   f the wall's part of its face's area and Delta_n the cell's size along its normal. The
//   coefficients come from BuildingForcing::rate (kBuildingHeatOwn and the three after it),
//   at the state of its last update. The own cell's coefficient goes into a_P; the reference
//   cell's part is in the residual alone. No stability functions: the exchange does not feel
//   the temperature difference it makes (psi_h = 0);
// - the diffusive flux of theta through the blocked part of a face is closed, in the residual
//   and in the implicit rows: no heat is conducted out of a building, the wall's flux above
//   is the exchange.
//
// theta_b: with buildings.theta_mode = ambient the reference profile at the cell's height
// (bc.theta_profile), each cell's theta against its own height's (the reference cell's too),
// so a state at that profile exchanges no heat; with fixed,
// buildings.theta_surface, the same on every wall. buildings.theta_mode = none leaves theta
// as it was before these terms. Not here: a given heat flux, stability functions on the
// walls, a temperature that differs from wall to wall, the heating of a building in time.
//
// Ghost cells: the rows read theta one cell out of a partly solid wall cell (its reference
// cell), so theta needs its one ghost layer filled; theta_b has one layer, as the reference
// profile's field.

#include "BuildingForcing.H"
#include "Forcing.H"
#include "ImplicitOperator.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace gust {

//! The theta terms of one level's buildings.
class BuildingTheta
{
public:
    //! law: its rates give the walls' exchange coefficients (updated by the residual before
    //! the rows here are applied). z_nd: the mesh's node heights; detJ: the cells' (for
    //! wall_heat). ref: the reference profile of buildings.theta_mode = ambient.
    void build (const Buildings& b, const BuildingForcing& law, const BuildingsOptions& opt, const amrex::Geometry& geom,
                const amrex::MultiFab& z_nd, const amrex::MultiFab& detJ, const ReferenceProfile& ref);

    //! theta_b at the cells (one ghost layer).
    [[nodiscard]] const amrex::MultiFab& surface_theta () const { return m_tb; }

    //! Gt[d], the diffusive flux of theta through the d-faces: times the open part of each face.
    void close_fluxes (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open,
                       amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& Gt) const;
    //! The rows of R_theta in solid cells and the walls' exchange in wall cells (last, after
    //! every other term of R_theta). theta with its ghost layer filled.
    void apply_rows (const amrex::MultiFab& theta, amrex::MultiFab& Rt) const;
    //! The same in the implicit coefficients of theta: the own coefficient in a_P, identity
    //! rows in solid cells, and the couplings through the closed parts of the faces.
    void add_to_implicit (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open, amrex::MultiFab& A) const;
    //! 0 in the solid cells of a field on the cells: their rows of R_theta are left out of the
    //! residual norm (Simplec.H: stage_norms), as the rows of k are.
    void zero_solid_cells (amrex::MultiFab& field) const;
    //! The heat the walls take from the air, sum of f q times the wall's area (K m^3/s), at
    //! the rates of the law's last update. Collective: the same value on every rank.
    [[nodiscard]] amrex::Real wall_heat (const amrex::MultiFab& theta) const;

private:
    const Buildings* m_b{nullptr};
    const BuildingForcing* m_law{nullptr};
    const amrex::MultiFab* m_detJ{nullptr};
    amrex::Geometry m_geom;
    amrex::MultiFab m_tb;
};

} // namespace gust

#endif
