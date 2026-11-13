#ifndef GUST_TIME_DRIVER_H_
#define GUST_TIME_DRIVER_H_

// The time driver of unsteady RANS by dual time stepping (steady.stage =
// unsteady). An outer loop over physical steps t^n -> t^{n+1} = t^n + dt;
// each step is one backward-Euler (BDF1) step of the steady operator,
//
//   R(q^{n+1}) + (q^{n+1} - q^n) / dt = 0,
//
// solved by the pseudo-time iteration of the steady driver (Driver.H, the
// SIMPLEC stage with its local dtau and SER ramp), with the physical time
// term of TimeTerm.H in the residual and on the implicit diagonals.
//
// One step (TimeDriver::advance):
//   1. the ghost cells of q are filled (the caller's GhostFill) and dt is
//      chosen: unsteady.dt, or with unsteady.cfl > 0 the global step
//      convective_dt(q) = cfl / max over the cells of
//      (|u_c|/dx + |v_c|/dy + |Omega_c|/(J dzeta)), the convective rate of
//      local_dtau (PseudoTime.H) without the diffusive term or the floor;
//      the last step ends the run at unsteady.stop_time (a step within 0.1 %
//      of the time left takes all of it, so no sliver of a step is left);
//      a dt that is not finite, or a target CFL over a field at rest (the
//      largest rate at kMinRate), aborts naming unsteady.cfl;
//   2. q^n <- q (TimeTerm::store) and the term's dt is set;
//   3. inner iterations, Driver::step, until every scaled residual R_q
//      (the stage's norms, time term included) is at or below
//
//        max(10^-d R_q^1, steady.residual_tol),   d = unsteady.inner_drop_orders,
//
//      or unsteady.max_inner iterations have run. R_q^1 is the first inner
//      iteration's norm, the residual at q = q^n, where the time term is 0:
//      the steady residual of the last step's answer. The floor at
//      steady.residual_tol stops a step whose start is already steady
//      (and a field such as k, held at its floor in some cells, whose
//      residual does not fall in proportion) instead of asking for
//      round-off. The stage's norm of iteration m is the residual before
//      that iteration's update, so a converged step has had one more update
//      than its last norm shows.
//   A step that reaches max_inner is reported and the run goes on; a
//   residual that becomes NaN or infinite (the driver's Diverged) ends it.
//
// The driver's own verdict (steady.residual_drop_orders, the monitors, the
// confidence) is not a stop here, and with SteadyOptions::stage = unsteady
// the driver ramps the CFL whatever the verdict (Driver::step). Its SER ramp
// of the pseudo-time CFL and its histories run on across the steps: the
// reference hold at steady.cfl0 lasts the first steady.residual_ref_iterations
// inner iterations of the run, which may span several steps, and a restart
// continues the checkpoint's history (no new hold, its CFL and reference).
//
// The time: t^n is accumulated in double (n dt for a fixed step), so that
// the last-step decision does not drift in a single-precision build.
//
// MPI: the stage's norms and convective_dt are all-reduced, so every rank
// takes the same dt and the same inner count.
//
// Inputs (unsteady.*, read with steady.stage = unsteady): see
// UnsteadyOptions. Not in this module yet: the checkpoint of q^n and the
// time, averages and probes (the unsteady outputs).

#include "Driver.H"
#include "Mesh.H"
#include "State.H"
#include "TimeTerm.H"

#include <AMReX_Geometry.H>

#include <functional>
#include <string>
#include <vector>

namespace gust {

struct UnsteadyOptions
{
    //! unsteady.stop_time [s] (finite, > 0; required): the physical time the run ends at.
    amrex::Real stop_time{0};
    //! unsteady.dt [s] (finite, > 0): the physical step; exactly one of unsteady.dt and unsteady.cfl.
    amrex::Real dt{0};
    //! unsteady.cfl (0: use unsteady.dt; else finite and > 0): the target convective CFL of a global dt.
    amrex::Real cfl{0};
    //! unsteady.inner_drop_orders (finite, > 0 and <= 15; 3): the inner stop's drop below the residual at q^n.
    //! 3 keeps the algebraic error of a slowly converging inner solve below a tenth of the BDF1 error; 1 does not.
    amrex::Real inner_drop_orders{3};
    //! unsteady.max_inner (>= 1; 50): the inner iterations of a step at most.
    int max_inner{50};
    std::vector<std::string> given;   //!< the unsteady.* keys given
};

//! Range checks, each message naming unsteady.<key>.
std::vector<std::string> check_unsteady_options (const UnsteadyOptions& opt);
//! Read unsteady.*; aborts, naming the key, on a bad value, a repeated key or a key
//! that takes one value given several. Unknown unsteady.* keys are left to
//! the caller's unknown-key check (they stay unread).
UnsteadyOptions read_unsteady_options ();
//! The keys read_unsteady_options reads.
const std::vector<std::string>& unsteady_keys ();

//! cfl / (the largest convective rate of local_dtau over the valid cells), the
//! rate floored at kMinRate. q's u and v need one filled ghost cell. Collective.
amrex::Real convective_dt (const State& q, const Mesh& mesh, const amrex::Geometry& geom, amrex::Real cfl);

//! The check of a physical step dt from TimeDriver::advance: an empty string when dt is finite and > 0
//! and, with cfl > 0, below cfl / (2 kMinRate) (a field with a convective rate, not one at rest);
//! else the message, naming unsteady.cfl or unsteady.dt.
std::string check_step_dt (amrex::Real dt, amrex::Real cfl);
//! true when a step of dt, with left of the run's time to go, is the last: left <= dt (1 + 1e-3). The last
//! step then takes all of left (at most 0.1 % longer than dt).
bool is_last_step (double left, double dt);

//! One physical step as the time driver ran it.
struct StepRecord
{
    int step{0};                       //!< 1-based
    amrex::Real time{0};               //!< t^{n+1}, at the end of the step [s]
    amrex::Real dt{0};                 //!< [s]
    int inner{0};                      //!< inner iterations
    //! max_q R_q / max(R_q^1, steady.residual_tol) at the last inner iteration.
    amrex::Real inner_residual{0};
    bool converged{false};             //!< the inner stop met (false: unsteady.max_inner reached)
    bool diverged{false};              //!< a residual became NaN or infinite
};

class TimeDriver
{
public:
    //! Fills the ghost cells of q (the boundary conditions' fill_ghosts).
    using GhostFill = std::function<void (State&)>;

    //! term and driver are referenced: the driver's stages must read term through their
    //! ResidualInputs::time. residual_tol: steady.residual_tol, the inner stop's floor.
    //! Aborts on check_unsteady_options failures and on a residual_tol that is not finite and > 0.
    TimeDriver (const UnsteadyOptions& opt, TimeTerm& term, Driver& driver, const Mesh& mesh,
                const amrex::Geometry& geom, amrex::Real residual_tol, GhostFill fill);

    //! One physical step (the header's steps 1 to 3). Aborts when the run has already reached
    //! unsteady.stop_time or diverged.
    const StepRecord& advance (State& q);
    //! true once the time has reached unsteady.stop_time or a step diverged.
    [[nodiscard]] bool finished () const { return m_diverged || m_done; }
    [[nodiscard]] bool diverged () const { return m_diverged; }
    [[nodiscard]] amrex::Real time () const { return static_cast<amrex::Real>(m_time); }
    [[nodiscard]] const std::vector<StepRecord>& records () const { return m_records; }
    //! The steps that reached unsteady.max_inner without meeting the inner stop.
    [[nodiscard]] int unconverged_steps () const;
    //! The inner iterations over the run.
    [[nodiscard]] long inner_iterations () const;

private:
    UnsteadyOptions m_opt;
    TimeTerm& m_term;
    Driver& m_driver;
    const Mesh& m_mesh;
    const amrex::Geometry& m_geom;
    amrex::Real m_tol;
    GhostFill m_fill;
    double m_time{0};   //!< t^n, in double (n dt for a fixed step)
    bool m_done{false};
    bool m_diverged{false};
    std::vector<StepRecord> m_records;
};

//! One line for a step: "step 3: t = ... s, dt = ... s, 12 inner iterations, inner residual ...".
std::string format_step (const StepRecord& s);

} // namespace gust

#endif
