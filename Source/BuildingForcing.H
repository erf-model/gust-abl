#ifndef GUST_BUILDING_FORCING_H_
#define GUST_BUILDING_FORCING_H_

// How buildings act on momentum (Buildings.H has the geometry; docs/buildings.rst):
//
// - a blocked face (fraction 1) is a fixed row, u = 0: R = u and the identity in the
//   implicit rows, no 1/dtau, d = 0 in SIMPLEC, the steady limit of ERF's
//   if_implicit_projection, with no drag coefficient; w on the ground stays the
//   boundary conditions' row, as they set it from the slope;
// - partial faces carry the open part of the mass flux (Mesh::set_open_faces);
// - buildings.wall_form = fraction_stress: every wall of a building takes the form 7
//   wall law (ERF erf.if_wall_form = fraction_stress, Kynema-SGF wall_form 7), neutral:
//     u*^2 = (1 - beta) u*_own^2 + beta u*_ref^2,  u*_x = kappa U_x / ln(d_x / z0),
//   with the own distance floored at 10 z0 (Kynema-SGF 996cdc9e) and u*^2 = u*_own^2 where
//   there is no reference cell (d2 = 0). U is the speed
//   parallel to the wall, at the cell centre: (u, v) at a roof or underside, (v, w) at a
//   side normal to x, (u, w) at one normal to y. Each component is the mean of the cell's
//   two faces that are not fixed rows (buildings.wall_speed = open_faces, the default with
//   this form), in the cell and in its reference cell: in an edge or corner cell one face
//   can lie on the other wall, and the mean with its zero would halve the speed. all_faces
//   takes both faces, as before. A wall is a sink on the two velocity
//   components parallel to it,
//     S = -r u_t,  r = f u*^2 / (Delta_n |U_own|),
//   f the wall's part of its face's area (wall_area_fraction: the blocked part of the face in
//   a fluid cell, the difference of the two opposite faces' in a partly solid one), Delta_n
//   the cell's size along the normal
//   (the layer's height, dx or dy). A cell with several walls (an edge: roof and side;
//   a corner: two sides) takes the sum of its walls' rates on each component, so the
//   component along an edge feels both walls. The viscous flux of a parallel component
//   through the blocked part of a wall's face is closed, in the residual and in the
//   implicit rows. r stays out of the local pseudo-time step and goes into a_P
//   (its derivative with r frozen; the fixed point does not depend on it).
// - buildings.wall_form = roof_stress: the law on roofs and undersides alone, the sides
//   with their fixed rows and viscous fluxes to the zero velocity (the form before the
//   sides had the law).
// - buildings.wall_form = drag: the fixed rows alone (Kynema-SGF's laminar pathway).
//
// k at buildings: BuildingTke.H. theta: BuildingTheta.H. On a slope a roof's speed is horizontal: its
// tilt is not taken out. A side is vertical, so (v, w) and (u, w) are parallel to it.

#include "Buildings.H"
#include "ImplicitOperator.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
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

//! The component of BuildingForcing::rate() after the three rates: the largest u*^2 of the
//! cell's walls (m^2/s^2), which the closure's wall value of k reads (BuildingTke.H).
inline constexpr int kBuildingUstar2 = 3;
//! The heat exchange of the cell's walls with the law (1/s; BuildingTheta.H): the sum over
//! the walls of the coefficient on the cell's own theta - theta_b, and, in a partly solid
//! cell, the coefficient on the reference cell's theta - theta_b for its wall along x, y, z
//! (such a cell has one wall along a direction at most; the reference cell is the next one
//! out from that wall).
inline constexpr int kBuildingHeatOwn = 4, kBuildingHeatRefX = 5, kBuildingHeatRefY = 6, kBuildingHeatRefZ = 7;
inline constexpr int kNumBuildingRates = 8;

//! The heat transfer velocities (m/s) of one wall under the form 7 law with a given surface
//! temperature, neutral (ERF ib_wall::wall_state, surface_temp, psi_h = 0):
//!   q = u* theta*,  theta* = (1 - beta) kappa (theta - theta_b) / ln(d1 / z0)
//!                          +      beta  kappa (theta_ref - theta_b) / ln(d2 / z0),
//! q the kinematic heat flux from the air to the wall (K m/s), with the momentum law's u*
//! (us2 = u*^2), its floor of 10 z0 on both distances, and the own term alone where there is
//! no reference cell (d2 = 0). own multiplies theta - theta_b, ref multiplies
//! theta_ref - theta_b. Both are 0 at rest.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void form7_heat_velocities (amrex::Real us2, amrex::Real beta, amrex::Real d1, amrex::Real d2, amrex::Real z0, amrex::Real kappa,
                            amrex::Real& own, amrex::Real& ref) noexcept
{
    const amrex::Real floor = amrex::Real(10) * z0;
    const amrex::Real c = kappa * std::sqrt(us2);
    own = c / std::log(amrex::max(d1, floor) / z0);
    ref = 0;
    if (d2 <= amrex::Real(0)) { return; }
    // beta = 0 (a fluid cell, which can have a wall on both sides of a direction) gives ref = 0
    // exactly: BuildingTheta reads no reference cell there.
    own *= amrex::Real(1) - beta;
    ref = beta * c / std::log(amrex::max(d2, floor) / z0);
}

//! The mean of a velocity component over a cell's two faces, with the faces' weights (1, or 0 on
//! a fixed row with buildings.wall_speed = open_faces): an edge or corner cell then takes the open
//! face's value, not half of it. Both weights 1 give (a + b) / 2; both 0 give 0.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real open_face_mean (amrex::Real a, amrex::Real b, amrex::Real wa, amrex::Real wb) noexcept
{
    return (wa * a + wb * b) / amrex::max(wa + wb, amrex::Real(1));
}

//! The part of a cell's face area that a wall takes: in a fluid cell (beta = 0) the blocked part
//! f_wall of the face the wall lies on; in a partly solid cell the difference from the opposite
//! face's blocked part, which is the area facing that way (Buildings_K.H: wall_on_low_side). On
//! a slope the solid fraction of a roof's cells changes a little from column to column, so they
//! carry a side bit with a small such area; the whole face's fraction would count the building's
//! own cross-section as a wall.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wall_area_fraction (amrex::Real f_wall, amrex::Real f_opposite, amrex::Real beta) noexcept
{
    return beta > amrex::Real(0) ? amrex::max(f_wall - f_opposite, amrex::Real(0)) : f_wall;
}

//! Speeds below this do not divide the rate (m/s): u*^2 / U stays finite at rest.
inline constexpr amrex::Real kBuildingMinSpeed = amrex::Real(1.0e-3);

//! The momentum terms of one level's buildings.
class BuildingForcing
{
public:
    void build (const Buildings& b, const BuildingsOptions& opt, const amrex::Geometry& geom, const amrex::MultiFab& detJ);

    [[nodiscard]] bool wall_law () const { return m_wall_form != BuildingWallForm::drag; }
    //! The law on the sides as well as on roofs and undersides.
    [[nodiscard]] bool side_law () const { return m_wall_form == BuildingWallForm::fraction_stress; }
    //! The rates at the cells (one ghost layer) from the last update: component c < 3 is the
    //! rate on velocity component c, the sum over the cell's walls parallel to it; component
    //! kBuildingUstar2 is the largest u*^2 of the cell's walls; kBuildingHeatOwn and the three
    //! after it are the walls' heat exchange coefficients.
    [[nodiscard]] const amrex::MultiFab& rate () const { return m_rate; }

    //! The wall cells with a wall under the law, and those among them whose fluid centroid is
    //! within 10 z0 of such a wall, where the law reads the speed as at 10 z0 (all ranks).
    [[nodiscard]] amrex::Long wall_cells () const { return m_wall_cells; }
    [[nodiscard]] amrex::Long cells_below_floor () const { return m_below_floor; }

    //! The rates at the state (u, v, w on their faces). Roofs and undersides read valid
    //! points; the sides read one ghost cell of each component, at their reference cells.
    void update (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w);
    //! R_c += r_c q_c with the rates of update, a row taking the mean of its two cells.
    void add_sink (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w, amrex::MultiFab& Ru,
                   amrex::MultiFab& Rv, amrex::MultiFab& Rw) const;
    //! G[c][d], the viscous flux of component c through its control volumes' d-faces: times
    //! the open part of that face, for each wall normal d with the law and c parallel to it.
    //! open[d]: Mesh::open_fraction(d), one ghost layer.
    void close_viscous_fluxes (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open,
                               amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, AMREX_SPACEDIM>& G) const;
    //! R_c = q_c on the blocked faces (last: they replace every other row).
    void apply_rows (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w, amrex::MultiFab& Ru,
                     amrex::MultiFab& Rv, amrex::MultiFab& Rw) const;
    //! Component c's rows: a_P gains r_c, the couplings through the closed parts of the
    //! wall faces go, and the blocked faces become identity rows.
    void add_to_implicit (int c, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open, amrex::MultiFab& A) const;
    //! mask = 1 on the blocked faces of component c as well.
    void add_fixed_rows (int c, amrex::MultiFab& mask) const;

private:
    const Buildings* m_b{nullptr};
    amrex::Geometry m_geom;
    const amrex::MultiFab* m_detJ{nullptr};
    amrex::MultiFab m_rate;
    std::array<amrex::MultiFab, 3> m_weight;   //!< the faces' weights in a cell's speed (one ghost layer)
    amrex::Long m_wall_cells{0}, m_below_floor{0};
    amrex::Real m_z0{0.1};
    BuildingWallForm m_wall_form{BuildingWallForm::fraction_stress};
};

//! How many fluid cells have no open path to the largest-index region. Aborts are
//! the caller's. A sealed pocket has its own null vector in the pressure correction.
amrex::Long count_sealed_fluid_cells (const Buildings& b, const amrex::Geometry& geom);

} // namespace gust

#endif
