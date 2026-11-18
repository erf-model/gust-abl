#ifndef GUST_NESTS_H_
#define GUST_NESTS_H_

// Nests run in the same process as the outer domain (plan rev. 6 tracker,
// PR-69). The inputs give each nest's name and horizontal extent only:
//
//   nest.names     = nest1 nest2      names of the nests (identifiers)
//   nest.nest1.lo  = x y              the nest's low corner (m)
//   nest.nest1.hi  = x y              its high corner (m)
//   nest.feedback  = none | defect | relax | relax_consistent | inject
//                                     one way (default), or the joint run
//                                     with the defect source, the relaxation
//                                     (plain or consistent) or the injection
//   nest.feedback_factor = 1          the relaxation's factor (relax and
//                                     relax_consistent), or the
//                                     injection's fraction per exchange
//                                     (inject; default 0.5, at most 1)
//   nest.feedback_flux = outer | nest the flux the outer cells just outside
//                                     the defect source take through its
//                                     lateral faces (defect only; NestFaceFlux.H)
//   nest.damping   = w_theta | outer  the nests' Rayleigh layer: on w and
//                                     theta only (default), or the outer
//                                     deck's on every field it names
//
// Everything else is the outer domain's: the deck's physics, closure,
// solver and terrain keys, and these derived settings:
//   * mesh: refinement 2 2 1 (kNestRatio), the outer z-levels, one grid
//     level (mesh.hierarchy_levels = 1); periodic only in a periodic outer
//     direction the nest spans whole;
//   * faces: bc.xlo, bc.xhi, bc.ylo and bc.yhi = nest (periodic where the
//     nest is), ground and top as the outer domain's; bc.nest_parent = the
//     outer run's final checkpoint; no sponge, and a Rayleigh layer on w and
//     theta only, or with nest.damping = outer the outer deck's
//     (apply_nest_boundary_inputs);
//   * terrain: the outer deck's; a terrain.type = agnesi hill stays at the
//     outer domain's centre (TerrainInputs::hill_centre, set from the parent);
//   * forcing: the outer run's (forcing.ug, forcing.vg, the horizontal
//     forcing.pressure_gradient and bc.surface_heat_flux after a column
//     guess's adjustments, App.cpp OuterOutcome);
//   * first guess: init.first_guess = nest; no grid sequencing; no probes;
//   * outputs: io.plot_file, io.chk_file and io.turb_plot_file with
//     "_<name>" appended.
// The outer domain runs to steady and writes its checkpoint (it stops, and no
// nest runs, if it does not reach steady); then each nest runs from it, one
// after the other (and the run stops with an error after the last if one did
// not reach steady). STF (terrain.smoothing = 1) over terrain is refused: it
// smooths over each mesh's own ground. One way: the outer domain does not see the nests, and two
// nests that touch each take the outer values on their shared face but
// correct their own net flux, so their normal velocities there differ.
//
// nest.feedback = defect (JointRun.H): after the one-way runs above (whose
// outputs are written as usual), the outer domain and every nest iterate
// together in memory: each exchange the outer domain gets the nests' defect
// source (NestDefect.H) and each nest's faces the current outer domain
// (NestFaces::refresh), until every domain is steady in the same exchange;
// then each writes its outputs again, at its new step.
//
// nest.feedback = relax (NestRelax.H): the same joint run with the relaxation
// of the outer domain towards the restricted nests over their footprints in
// place of the defect source; nest.feedback_factor (default 1, finite and
// > 0) scales its rate.
//
// nest.feedback = relax_consistent (NestRelax.H): the relaxation with the
// residual of its target, weighted by the relaxation's weight, taken from the
// source as well, so the target satisfies the outer rows where the weight is
// 1 (the parent's pressure there mostly follows the nests'); same factor.
//
// nest.feedback_flux = nest (NestFaceFlux.H; defect only): the outer cells
// just outside the defect source's region take the nests' own fluxes through
// its lateral faces (summed over the nest faces) in place of the outer
// scheme's flux of the restricted nest state (the default, outer).
//
// nest.feedback = inject (NestInject.H): the same joint run with the outer
// state replaced by the restricted nests on the defect source's region at
// each exchange (a fraction nest.feedback_factor of the way, default 0.5,
// finite, > 0 and <= 1), held as data through the outer iteration.

#include "BoundaryConditions.H"
#include "Inputs.H"

#include <AMReX_REAL.H>

#include <array>
#include <string>
#include <vector>

namespace gust {

struct ForcingOptions;
struct SurfaceLayerOptions;

struct NestSpec
{
    std::string name;
    //! The corners in x and y (m), as given.
    std::array<amrex::Real, 2> lo{0, 0}, hi{0, 0};
};

//! nest.feedback: none (one way, the default), defect (the joint run with the defect source, NestDefect.H,
//! JointRun.H), relax (the joint run with the relaxation, NestRelax.H), relax_consistent (the relaxation with
//! the weighted residual of its target, NestRelax.H, so that the target solves the relaxed equations where the
//! weight is 1) or inject (the joint run with the injection, NestInject.H).
enum class NestFeedback { none, defect, relax, relax_consistent, inject };

//! nest.feedback_flux: outer (the default: the outer scheme's fluxes) or nest (the nests' own summed face fluxes
//! through the defect source's lateral faces, NestFaceFlux.H; nest.feedback = defect only).
enum class NestFeedbackFlux { outer, nest };

//! nest.damping: w_theta (the default: the nests' Rayleigh layer damps w and theta only) or outer (the outer deck's
//! layer on every field bc.rayleigh_fields names, towards the bc.inflow_* profile the outer domain relaxes to).
enum class NestDamping { w_theta, outer };

//! The value of nest.feedback as the deck gives it.
const char* nest_feedback_name (NestFeedback f);

struct NestOptions
{
    std::vector<NestSpec> nests;
    NestFeedback feedback{NestFeedback::none};
    //! nest.feedback_factor: relax's and relax_consistent's factor (finite, > 0; default 1) or inject's fraction (finite, > 0 and <= 1;
    //! default kInjectFactorDefault); not given with other feedbacks.
    amrex::Real feedback_factor{1};
    //! nest.feedback_flux (nest.feedback = defect only).
    NestFeedbackFlux feedback_flux{NestFeedbackFlux::outer};
    //! nest.damping.
    NestDamping damping{NestDamping::w_theta};
    //! Value-count and name problems found while reading.
    std::vector<std::string> errors;
    [[nodiscard]] bool active () const { return !nests.empty(); }
};

//! Read nest.names, nest.feedback, nest.feedback_factor, nest.feedback_flux, nest.damping and each nest's nest.<name>.lo and
//! nest.<name>.hi. No placement checks (check_nests).
NestOptions read_nest_options ();

//! The placement of every nest in the outer mesh: corners on the outer cell
//! faces, at least kNestMargin outer cells inside each non-periodic outer
//! face (inside the domain across a periodic one), not the whole outer
//! domain in x and y, no two nests overlapping. One message per problem,
//! each naming the key.
std::vector<std::string> check_nests (const NestOptions& opt, const MeshInputs& outer);

//! The nest's mesh inputs: its corners in x and y (placed on the outer cell
//! faces exactly), twice the outer cells there, the outer z (cells, ends,
//! dz0, stretch ratio, z-levels file), periodic only across the whole width
//! of a periodic outer direction, one grid level, the outer max_grid_size.
//! Call after check_nests reports no problem.
MeshInputs nest_mesh (const NestSpec& nest, const MeshInputs& outer);

//! The output name of a nest: base + "_" + name (empty base: empty).
std::string nest_output_name (const std::string& base, const std::string& name);

//! The nest's inputs from the outer deck's: nest_mesh; the io names with
//! "_<name>"; no io.restart; the terrain keys' "given" record cleared (the
//! outer domain's checks cover the deck).
void apply_nest_inputs (const NestSpec& nest, Inputs& in);

//! The nest's boundary inputs from the outer deck's: the four lateral faces
//! nest with parent the checkpoint parent_chk (periodic in a direction where
//! the nest's mesh, periodic, is); no sponge (it relaxes
//! towards the inflow profile, and a nest has no inflow face); a Rayleigh
//! layer on w and theta only (u, v and k would relax towards the inflow
//! profile, not the parent), or with outer_damping (nest.damping = outer)
//! the outer deck's layer on every field it names (rayleigh_from_outer: the
//! parent is damped the same way, towards the same bc.inflow_* profile, over
//! the footprint); the relaxation zone off; the "given" record cleared (the
//! outer domain's checks cover the deck).
void apply_nest_boundary_inputs (const std::string& parent_chk, const amrex::Array<int, 3>& periodic, bool outer_damping,
                                 BoundaryInputs& bin);

//! What the outer domain's run hands its nests: whether it reached steady, and
//! the forcing and surface heat flux it ran with (a column first guess with
//! a target wind or init.column_stability = mol replaces the deck's values;
//! a nest takes the run's, so its interior is forced as its faces' parent).
struct OuterOutcome
{
    bool steady{false};
    bool forcing_set{false};
    amrex::Real ug{0}, vg{0};
    std::array<amrex::Real, 2> pressure_gradient{0, 0};
    amrex::Real surface_heat_flux{0};
};

//! The outer run's forcing.ug, forcing.vg, horizontal
//! forcing.pressure_gradient and bc.surface_heat_flux, kept for its nests.
void keep_outer_forcing (const ForcingOptions& fopt, const SurfaceLayerOptions& surf, OuterOutcome& out);

//! A nest's forcing and surface heat flux: the outer run's (no change when
//! the outer run kept none).
void take_outer_forcing (const OuterOutcome& out, ForcingOptions& fopt, SurfaceLayerOptions& surf);

} // namespace gust

#endif
