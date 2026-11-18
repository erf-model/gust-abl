#ifndef GUST_NEST_INJECT_H_
#define GUST_NEST_INJECT_H_

// The patches' feedback on the parent by injection (nest.feedback = inject;
// plan rev. 6 tracker, row N-114 item 2, left from PR-70). The joint run
// (JointRun.H) is the same as with nest.feedback = defect or relax, with this
// in place of their source or relaxation. Each exchange:
//
//   up: q* = the parent state with I q_f on the nests' footprints
//       (restrict_patch, NestDefect.H), then
//
//     q_c := q_c + omega (q* - q_c)   on chi_s,
//
// for u, v, w, p, theta and k, with
//   q_c    the parent's state;
//   q*     the restricted nests (I q_f) on their footprints, the parent's
//          state elsewhere;
//   omega  nest.feedback_factor (0 < omega <= 1; default
//          kInjectFactorDefault): the fraction of the way to the nest taken
//          per exchange;
//   chi_s  the region of the defect source (DefectRegion::inner): the
//          parent's points at least kDefectMargin parent cells inside the
//          footprint's lateral faces, all z (cells; faces from L to H + 1
//          along their normal, so chi_s's boundary faces are injected too;
//          w on the ground and top faces is left to the boundary condition,
//          which on the terrain makes it follow u).
//          The ring between chi_s and the footprint's faces keeps the
//          parent's own state and equations, as with defect and relax: the
//          nest's face values come from the parent there.
//
// The parent's iteration holds chi_s as data, a Dirichlet region:
// - residual (Residual.cpp, the end of residual()): the rows of u, v, w, p
//   (continuity), theta and k are zero on chi_s (zero_rows), after every
//   other term and row, so the wall-k rows and the ground and top rows
//   under a nest too. The residual norms and the confidence stop see only
//   the points outside chi_s; the RMS is still over all the parent's points
//   (zeros on chi_s).
// - SIMPLEC (SimplecStage::iterate): the velocity on chi_s is reset to the
//   target after the momentum predictor, and every field at the end of the
//   iteration (hold). The pressure correction holds chi_s's cells and faces
//   out of its Poisson solve (FaceCoeffs::held_fluid, held_open; as
//   buildings' solid cells): their rows are p' = 0, no flux crosses their
//   faces (the operator's areas and the returned fluxes are both masked), so
//   the ring cells beside chi_s see the injected boundary fluxes as given
//   and chi_s's faces are not corrected. Without this the correction moves
//   chi_s's faces, the next injection resets them, and the parent's
//   iteration does not settle. The line sweeps still solve chi_s's points
//   with the rest (their rows are zero but coupled), so a ring row sees
//   its held neighbours move within a sweep until hold resets them; at the
//   fixed point those moves vanish. With simplec.theta = frozen, theta on
//   chi_s is still injected (the nest's theta).
// - The Newton stage is refused (Sequencing.cpp; App.cpp refuses Newton with
//   any joint feedback).
//
// Why omega < 1: with omega = 1 (plain injection) the parent and the nests
// are two Dirichlet-coupled regions overlapping by kDefectMargin parent
// cells, each taking one pseudo-time step per exchange. On the
// rt_nests_two_way hill omega = 0.25, 0.5 and 0.6 reach the same steady
// state (omega moves only the path: at the fixed point q_c = q* on chi_s
// for any omega); at 0.75 the difference grows after about 8 exchanges into
// an oscillation, and at 1 the oscillation is larger
// (docs/boundary_conditions.rst, two-way-inject).
//
// Pressure: p is injected with restrict_patch's level shift (the restricted
// p has the parent's mean over the footprint). chi_s's boundary faces are
// held, so no row the parent solves reads p on chi_s: the injected p is the
// parent's output there.
//
// After a steady joint run, run_joint injects once in full (omega = 1), so
// the parent's output under a nest is I q_f of the final nest. A run stopped
// before steady (steady.max_iterations) writes the iterate, so its restart
// continues the straight run bit for bit.
//
// Limit: with turb.nut_relax < 1 the closure's nu_t on chi_s blends with
// the parent's own history (Closure::accept), as with defect and relax; it
// enters the ring's diffusive fluxes on chi_s's boundary faces.

#include "NestDefect.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <array>
#include <string>
#include <vector>

namespace gust {

//! nest.feedback_factor's default with nest.feedback = inject (NestInject.H, Why omega < 1).
constexpr amrex::Real kInjectFactorDefault = amrex::Real(0.5);

//! nest.feedback_factor's check with nest.feedback = inject: finite, > 0 and <= 1. Empty: no problem.
std::string check_inject_factor (amrex::Real factor);

class NestInject
{
public:
    //! On the parent's layout and geometry (for the owner mask of differences()), with omega = factor
    //! (aborts unless check_inject_factor passes); inactive until set.
    NestInject (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::Geometry& geom,
                amrex::Real factor = kInjectFactorDefault);

    //! The regions whose inner boxes are chi_s (DefectRegion; a direction where inner spans the footprint has no
    //! faces). Sets the injection active.
    void set (const std::vector<DefectRegion>& regions);

    //! The target := q + omega (qstar - q) on chi_s (qstar itself with full; kept until the next call), then
    //! q := the target there for u, v, w, p, theta and k (valid points; the caller fills q's ghosts after).
    //! Nothing while inactive.
    void inject (const State& qstar, State& q, bool full = false);

    //! q := the target on chi_s for the first nfields of u, v, w, p, theta, k (3: the velocity; 6: all), valid
    //! points; the caller fills the ghosts. SimplecStage::iterate holds chi_s with it after the momentum
    //! predictor (the velocity) and at its end (all). Nothing while inactive;
    //! aborts if active before the first inject.
    void hold (State& q, int nfields) const;

    //! R := 0 on chi_s for u, v, w, p (continuity), theta and k. Nothing while inactive.
    void zero_rows (ResidualFields& R) const;

    //! f := 0 on chi_s's points of f's type (valid points; f on the parent's layout); zero_rows uses it.
    //! Nothing while inactive.
    void zero_points (amrex::MultiFab& f) const;

    //! max |q - qstar| and its RMS over chi_s (each point once), for u, v, w, theta and k: 10 values, the
    //! maxima first. Zeros while inactive. Every rank gets the same values.
    [[nodiscard]] std::vector<amrex::Real> differences (const State& q, const State& qstar) const;

    //! The pressure correction's masks (FaceCoeffs::held_fluid, held_open; set builds them): a cell mask, 0 on
    //! chi_s's cells and 1 elsewhere, and per direction a face mask, 0 on chi_s's faces (its boundary faces
    //! included) and 1 elsewhere. Valid points, no ghosts; on the parent's layout.
    [[nodiscard]] const amrex::MultiFab& fluid_mask () const { return m_fluid; }
    [[nodiscard]] const amrex::MultiFab& open_mask (int d) const { return m_open[static_cast<std::size_t>(d)]; }

    //! The change of the target at the last inject (FeedbackChange, NestDefect.H; the joint run's stop window): max
    //! |new target - old target| over chi_s's data points for u, v, w, theta and k, the old target being q there
    //! (the stage holds it); sizes as NestRelax's (the largest |u|, |v|, |w| of the new target for the velocity,
    //! 0 for theta, the largest k for k). Never first: q on chi_s is the old target, or the outer state at the
    //! first exchange, so a joint restart measures the change as the run without it does.
    [[nodiscard]] const FeedbackChange& change () const { return m_change; }
    //! The target of u, v, w, theta and k (a joint checkpoint stores them: JointRecord.H).
    [[nodiscard]] std::array<const amrex::MultiFab*, 5> feedback_fields () const
    {
        return {&m_target[0], &m_target[1], &m_target[2], &m_target[4], &m_target[5]};
    }
    //! A joint restart: the target of the last exchange before the checkpoint (valid points). The state's chi_s
    //! holds it as well, and inject reads that, so this only keeps feedback_fields the run's own.
    void restore_previous (const std::array<const amrex::MultiFab*, 5>& f);

    [[nodiscard]] bool active () const { return m_active; }
    [[nodiscard]] amrex::Real factor () const { return m_factor; }
    [[nodiscard]] const std::vector<DefectRegion>& regions () const { return m_regions; }

private:
    amrex::BoxArray m_ba;
    amrex::DistributionMapping m_dm;
    amrex::Geometry m_geom;
    std::vector<DefectRegion> m_regions;
    amrex::Real m_factor;
    std::array<amrex::MultiFab, 6> m_target;   //!< u, v, w, p, theta, tke of the last inject (on chi_s)
    FeedbackChange m_change;
    amrex::MultiFab m_fluid;                   //!< 0 on chi_s's cells, 1 elsewhere
    std::array<amrex::MultiFab, 3> m_open;     //!< 0 on chi_s's x, y, z faces, 1 elsewhere
    bool m_active{false};
    bool m_held{false};
};

} // namespace gust

#endif
