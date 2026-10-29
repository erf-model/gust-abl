#ifndef GUST_PSEUDO_TIME_H_
#define GUST_PSEUDO_TIME_H_

// The local pseudo-time step: in every cell
//
//   dtau = CFL / ( |u_c|/dx + |v_c|/dy + |Omega_c|/(J dzeta) + 2 nu (1/dx^2 + 1/dy^2) ),
//
// with u_c, v_c the means of the cell's x- and y-face velocities, Omega_c the
// mean of its two z-face contravariant fluxes (MassFlux.H: the velocity
// through the terrain-following layers, in height units per second, so
// |Omega| / (J dzeta) is the rate at which it crosses the cell's physical
// thickness J dzeta) and nu the momentum diffusivity (cell-centred; null: 0).
// The vertical diffusive limit is left out: vertical diffusion is implicit in
// every stage (z-line solves, Predictor.H). The rate is floored at
// kMinRate so a cell at rest without diffusion gets a large finite step.
//
// The step floor (steady.dtau_floor = phi, 0.1 by default): the rate is also
// kept at or above
//
//   phi U_max / min(dx, dy),   U_max = max over the level's faces of |u|, |v|,
//
// the convective rate of the level's fastest cell. A cell at rest (a
// recirculation pocket, the air under a dense forest canopy) then steps at
// most as the fastest cell does, 1/phi times, instead of CFL / kMinRate.
// Why (the observation; the mechanism is not derived): with the plain
// local step such a cell's momentum rows have 1/dtau -> 0 and a vanishing
// drag derivative, so the SIMPLEC d = a_P - sum a_nb (Simplec.H) and with
// it the bound on the pressure correction's velocity update in the cell
// vanish; with a canopy the iteration cycles once the CFL cap passes
// about 3, and the floor removes the cycle (docs/canopy.rst: the forested
// escarpment and the same forest on flat ground, both steady at the
// default cap with the floor; a 2-D SIMPLEC model with the same terms did
// not reproduce the cycle, so what closes the loop in the 3-D solver is
// not identified). Not tried: a rate from the cell's own face speeds
// (the cell-centred |u_c| of a reversal cell is 0 while its faces move)
// or a neighbourhood rate, which would bound such cells without a
// collective. phi from 0.03 to 1 converges them; 0.1 is the default, the largest
// value that leaves the empty cases' paths close to what they were (the 2-D
// hill and the Newton switch on it unchanged, the steep ridge with the w cap
// 15 % more iterations); 1 is the fastest on the forest and slows the Newton
// hill and the ridge (docs/solution_algorithm.rst).
// phi = 0 is the plain local step (the form before the floor). U_max is one
// all-reduce per call, the same on every rank (a maximum).

#include "Mesh.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace gust {

//! 1/s; dtau <= CFL / kMinRate.
inline constexpr amrex::Real kMinRate = amrex::Real(1.0e-12);

//! dtau: cell-centred on the mesh's layout, no ghosts needed (defined here if
//! empty). q's u, v, w must have at least one ghost cell filled (Omega reads
//! the neighbouring faces of the slopes). phi: the step floor above (0:
//! none; collective when > 0). Aborts, naming the argument, on a wrong
//! layout, cfl <= 0 or a non-finite cfl, or a phi that is negative or not finite.
void local_dtau (const State& q, const amrex::MultiFab* nu, const Mesh& mesh, const amrex::Geometry& geom,
                 amrex::Real cfl, amrex::MultiFab& dtau, amrex::Real phi = 0);

//! The SER law of the CFL (steady.cfl0, cfl_max, ser_exponent p,
//! ser_max_growth g) after an iteration with relative residual r (> 0, or
//! 0 for a vanished residual), from the current CFL cfl:
//!
//!   min(cfl_max, g cfl, cfl0 r^-p).
//!
//! Below r_min = (cfl0 / cfl_max)^(1/p) the law already gives the cap, so the
//! power is never taken there (it would overflow for a vanished residual and
//! p > 1); with p = 0 the law is cfl0. The driver (Driver.H) and the
//! pseudo-transient continuation (Continuation.H) share it.
amrex::Real ser_cfl (amrex::Real cfl, amrex::Real r, amrex::Real cfl0, amrex::Real cfl_max, amrex::Real p,
                     amrex::Real growth);

} // namespace gust

#endif
