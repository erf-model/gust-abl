#ifndef GUST_LOG_SHEAR_H_
#define GUST_LOG_SHEAR_H_

// turb.shear_stencil = log_law: the vertical gradients of u and v in the strain invariant S2, exact
// for a logarithmic profile over the ground. The closure's central difference of a log profile
// over-estimates du/dz near the ground: by 21 % in the second cell of uniform layers (S2, and with
// it the shear production of k, by 46 % there; z0 -> 0), 6 % in the third, 3 % in the fourth.
// Here, in every cell k above the wall cell and below the top cell,
//
//   du/dz = (u_{k+1} - u_{k-1}) / [ (d_k + z0) ln((d_{k+1} + z0) / (d_{k-1} + z0)) ],
//
// d the closure's wall distance (the height above the ground; the fluid centroid's, d1, in the cut
// cell of an immersed ground, BuildingTke.cpp) and z0 the ground's roughness length (the column's,
// with a roughness map). With d_{k+1} - d_{k-1} = 2 dz << d_k it
// tends to the central difference (the factor is 1 + dz^2 / (3 d^2) + ...). S2 changes by the
// difference of the (u, z) and (v, z) terms of 2 S_ij S_ij, the other gradients as the closure has
// them. The wall cell keeps the closure's S2: its k is the wall value. With an immersed ground this
// takes the place of ImmersedStrain (its factor R, the flat ground's central-difference error, is 1
// here); solid cells are left as they are, and the cut cell takes the wall law's gradient
// (apply_cut_cell_strain, ImmersedStrain.H). A fluid cell with a side wall only (sloping ground) keeps
// the closure's S2.

#include "Buildings.H"
#include "Closure.H"
#include "SurfaceRoughness.H"

namespace gust {

class LogShear final : public StrainCorrection
{
public:
    //! dist: the height above the ground on the level's cells, one lateral ghost cell (read here
    //! only). rough: the ground's roughness, a uniform z0 > 0 (m) or a map on dist's layout with one
    //! lateral ghost column. ground: an immersed ground's geometry (it must outlive this object), or
    //! null with a body-fitted ground. Aborts on a uniform z0 that is not positive.
    LogShear (const amrex::MultiFab& dist, const Roughness& rough, const Buildings* ground);
    //! S2 on its valid cells and lateral ghost ring.
    void apply (const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2) const override;

private:
    amrex::MultiFab m_len;   //!< d + z0
    const Buildings* m_b;
    amrex::Real m_z0 = 0;    //!< the immersed ground's uniform z0 (its cut cells, ImmersedStrain.H)
};

}   // namespace gust

#endif
