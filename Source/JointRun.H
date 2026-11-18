#ifndef GUST_JOINT_RUN_H_
#define GUST_JOINT_RUN_H_

// The joint two-level run (nest.feedback = defect; plan rev. 6 tracker,
// PR-70; nest.feedback = relax, PR-71, takes NestRelax::set in place of
// NestDefect::set in step 1, NestRelax.H; nest.feedback = inject, row N-114 item 2,
// replaces the outer state on the source's region in step 1 instead, NestInject.H). The outer domain and its nests have each run one way to steady
// (Nests.H) and stay in memory (Run.cpp holds every run_domain open through
// DomainRole::joint). Then, each exchange:
//
//   1. up:    q* = the outer state with every nest restricted onto its
//             footprint (restrict_patch); S = chi_s R_c(q*)
//             (NestDefect::set);
//   2. the outer domain takes one iteration with R_c - S;
//   3. down:  each nest's parent images and faces from the outer state now
//             (NestParent::refresh, NestFaces::refresh);
//   4. each nest takes one iteration.
//
// The exchange runs every iteration (fixed). The run is steady when every
// domain's driver is steady after the same exchange (each domain's steady.*
// criteria, with S and the faces as they are then) and the feedback the outer
// domain was stepped with has settled. Each exchange's relative change from
// the last,
//
//   d_n = max over u, v, w, theta, k of  max|F_new - F_old| / max(size, floor)
//
// (relative_feedback_change, NestDefect.H; F = S, or the relaxation's target
// q* where its rate is > 0; size and floor per JointMember::feedback_floor),
// is tested over the outer domain's steady.monitor_window exchanges as a
// monitor is over that many iterations (FeedbackWindow, NestDefect.H): the
// window's sum and Aitken's extrapolated remainder at or below its
// steady.monitor_rel_tol. S is a constant
// of a step, so a domain's own verdict does not see how S moves between
// exchanges; this test does. The run also stops when one domain diverges,
// or after steady.max_iterations exchanges over the run. A restart from the
// outer domain's last checkpoint goes on with the exchanges (JointRecord.H;
// first_exchange the exchanges its record counts) with the window and the
// last feedback the record and checkpoint hold, so it stops after the same
// exchange as the run without a restart. A line per exchange (the first
// ten, then every tenth) gives the largest |S| per equation (relax: the
// differences under the nests), the feedback's change and the largest change
// of a nest's face velocities.

#include "Domain.H"

#include <string>
#include <vector>

namespace gust {

//! What the joint run did.
struct JointResult
{
    int exchanges{0};             //!< over the run, a restart's earlier exchanges included
    bool steady{false};
    bool diverged{false};
    double wall_seconds{0};       //!< the whole loop
    double exchange_seconds{0};   //!< restriction, source and refresh
    //! Not steady: every domain was steady at the last exchange, the feedback had not settled over the window.
    bool waiting_on_feedback{false};
    double feedback_change{0};    //!< the last exchange's relative change of the feedback
    double feedback_window{0};    //!< the sum of the changes over the window (infinite until it is full)
    std::vector<amrex::Real> window;   //!< the window's changes at the end, oldest first (for the joint record)
};

//! Run the exchanges: members[0] is the outer domain (with its defect source),
//! the others its nests (with their parent images and faces); first_exchange: the exchanges an earlier run took
//! (a joint restart; 0: none), so this run takes at most max_exchanges - first_exchange; prior_window: a joint
//! restart's window from its record (JointRecord.H; empty: none), the outer domain's feedback already restored.
JointResult run_joint (const std::vector<JointMember*>& members, int max_exchanges, int first_exchange = 0,
                       const std::vector<amrex::Real>& prior_window = {});

//! The widest domain name the joint run's log lines print (characters).
constexpr int kJointLogNameWidth = 24;

//! name as the joint run's log lines print it: as it is up to kJointLogNameWidth characters; a longer one keeps its
//! first 12 and last 9 characters around "..." (kJointLogNameWidth in all), so a suffix of up to 9 characters stays readable.
std::string joint_log_name (const std::string& name);

} // namespace gust

#endif
