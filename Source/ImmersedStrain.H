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
// layers), so a flat ground at a cell face is unchanged, on stretched layers too. The cut cell keeps the closure's S2: its k is the wall value (BuildingTke). S2 changes
// by the difference of the (u, z) and (v, z) terms of 2 S_ij S_ij, the other gradients as the
// closure has them. This follows the rule of ERF's erf.if_flat_wall_dissipation and kynema-sgf
// #2067 (an immersed wall repeats the flat ground's discretisation), applied to the shear.

#include "BuildingForcing.H"
#include "Closure.H"

namespace gust {

class ImmersedStrain final : public StrainCorrection
{
public:
    //! law: the immersed ground's wall law (law.immersed_ground()), whose geometry must outlive
    //! this object. z0: the ground's roughness length (m), > 0.
    ImmersedStrain (const BuildingForcing& law, amrex::Real z0);
    //! S2 in the first full cell above each cut cell, on S2's valid cells and lateral ghost ring.
    void apply (const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2) const override;

private:
    const Buildings* m_b;
    amrex::Real m_z0;
};

}   // namespace gust

#endif
