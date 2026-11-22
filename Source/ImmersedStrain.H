#ifndef GUST_IMMERSED_STRAIN_H_
#define GUST_IMMERSED_STRAIN_H_

// The strain invariant S2 next to an immersed ground (terrain.representation = immersed_forcing;
// docs/immersed_terrain.rst). In the first full cell above a cut cell the closure's central
// difference puts the cut cell's speed at its cell centre, though it stands for its fluid centroid,
// d1 = (1 - beta) dz / 2 above the wall. Here the vertical gradients of u and v in that cell are the
// log-law gradient through the centroid and the next centre, times the error of the flat ground's
// central difference at the same cell sizes, so that the cell makes the flat ground's error:
//
//   du/dz = R (u_{k+1} - u_c) / [ (d2 + z0) ln((d2 + h + z0) / (d1 + z0)) ],
//   R = (f3 + z0) ln((f4 + z0) / (f2 + z0)) / (2 dz_k),
//
// with c the cut cell, k = c + 1 the cell above it, d2 its distance from the wall (kWallD2z),
// h = z_{k+1} - z_k, and f2 = dz_c / 2, f3 = dz_c + dz_k / 2, f4 = f3 + h the centres of the same
// cells over a flat ground at the cut cell's bottom face. With beta = 0 (d1 = f2, d2 = f3) it is the
// closure's central difference (u_{k+1} - u_c) / (2 dz_k) (Stencil_K.H; f4 - f2 = 2 dz_k on uniform
// layers), so a flat ground at a cell face is unchanged, on stretched layers too.
// In the cut cell itself (apply_cut_cell_strain) the vertical gradient is the wall law's at its fluid
// centroid, u* / (kappa (d1 + z0)), as the body-fitted ground's first cell has it: the closure's
// difference reads the solid cell below. Its k is the wall value (BuildingTke); S2 there acts through
// turb.realizability. S2 changes
// by the difference of the (u, z) and (v, z) terms of 2 S_ij S_ij, the other gradients as the
// closure has them. This follows the rule of ERF's erf.if_flat_wall_dissipation and kynema-sgf
// #2067 (an immersed wall repeats the flat ground's discretisation), applied to the shear.

#include "BuildingForcing.H"
#include "Closure.H"

namespace gust {

//! S2 in the cut cells of the immersed ground b: the (u, z) and (v, z) terms with the wall law's
//! gradient. rough: the ground's roughness, a uniform z0 (m, > 0) or a map on the ground slab of b's
//! cells with one lateral ghost column (SurfaceRoughness.H; each column its own z0, one kernel per
//! kind). On S2's valid cells and lateral ghost ring.
void apply_cut_cell_strain (const Buildings& b, const Roughness& rough, const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2);

//! Aborts (naming op) unless rough is a uniform z0 > 0 or a map on the ground slab of b's cells.
void check_ground_roughness (const char* op, const Buildings& b, const Roughness& rough);

class ImmersedStrain final : public StrainCorrection
{
public:
    //! law: the immersed ground's wall law (law.immersed_ground()), whose geometry must outlive
    //! this object. rough: the ground's roughness, a uniform z0 (m) > 0 or a map (as
    //! apply_cut_cell_strain's; it must outlive this object).
    ImmersedStrain (const BuildingForcing& law, const Roughness& rough);
    //! S2 in each cut cell and in the first full cell above it, on S2's valid cells and lateral ghost ring.
    void apply (const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2) const override;

private:
    const Buildings* m_b;
    Roughness m_rough;
};

}   // namespace gust

#endif
