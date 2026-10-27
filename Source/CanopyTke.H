#ifndef GUST_CANOPY_TKE_H_
#define GUST_CANOPY_TKE_H_

// The forest canopy's sources of k (docs/canopy.rst), the steady form of
// ERF's canopy terms (erf-model/ERF#4340, extern_erf/ERF_CanopyTurbulence.H):
//
//   S_k += P_c - D_c,   P_c = beta_p f |U|^3,   D_c = beta_d f |U| k,
//
// with f = Cd a the drag field (Canopy.H) and |U| the cell-centred wind
// (the faces' means, as ERF). The drag removes mean kinetic energy f |U|^3
// per unit mass; a fraction beta_p becomes turbulence in the wakes of the
// leaves, at scales below the mesh, and dissipates fast (the short-circuit
// sink). Where the canopy pair alone dominates, k balances at
// beta_p |U|^2 / beta_d, the balance P/c of ERF's exponential stage update;
// with the length cap the dissipation C_mu0^3 k^1.5 f / alpha joins it
// (docs/canopy.rst: k = kappa_c |U|^2).
//
// Implicit side: D_c is linear in k, so its Patankar diagonal D_c / k and its
// exact derivative are both beta_d f |U| (0 on the exact diagonal at the k
// floor, where the kernel holds k fixed). P_c does not depend on k and stays
// an explicit source; its dependence on the wind is lagged, as the shear
// production's.
//
// The length cap L <= alpha / f is part of the closure's point kernel
// (Closure_K.H).

#include "Closure.H"
#include "State.H"
#include "SurfaceRoughness.H"

#include <AMReX_MultiFab.H>

namespace gust {

//! Fill cf.Pc and cf.Dc (defined here if empty) on the valid cells from the
//! state q (u, v, w on each cell's own faces: no ghosts read),
//! the drag field f (cell-centred, on q's cells) and the coefficients; add
//! beta_d f |U| to cf.dSdk, and subtract it from cf.dSdk_exact where k is
//! above tke_floor. k below the floor is taken at the floor, as the closure
//! takes it.
void add_canopy_tke (const State& q, const amrex::MultiFab& f, amrex::Real beta_p, amrex::Real beta_d,
                     amrex::Real tke_floor, ClosureFields& cf);

//! A ground under the trees as rough as a forest: the canopy's drag and k
//! terms already stand for the trees, so the ground's z0 should be that of
//! the soil and litter (0.01 to 0.1 m), not the forest's (about h/10).
inline constexpr amrex::Real kForestSizedZ0 = amrex::Real(0.5);

//! The columns (valid cells of every rank) whose first cell has canopy
//! (f > 0) and whose ground roughness (rough: uniform or map) is z0_max or more.
amrex::Long canopy_columns_on_rough_ground (const amrex::MultiFab& f, const Roughness& rough, amrex::Real z0_max);

} // namespace gust

#endif
