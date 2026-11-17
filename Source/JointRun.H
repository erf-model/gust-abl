#ifndef GUST_JOINT_RUN_H_
#define GUST_JOINT_RUN_H_

// The joint two-level run (nest.feedback = defect; plan rev. 6 tracker,
// PR-70; nest.feedback = relax, PR-71, takes NestRelax::set in place of
// NestDefect::set in step 1, NestRelax.H). The outer domain and its nests have each run one way to steady
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
// The exchange runs every iteration (fixed). The run stops when every domain's
// driver is steady after the same exchange (each domain's steady.* criteria,
// with S and the faces as they are then), when one diverges, or after
// steady.max_iterations exchanges over the run. A restart from the outer domain's
// last checkpoint goes on with the exchanges (JointRecord.H; first_exchange the
// exchanges its record counts). A line per exchange (the first ten, then
// every tenth) gives the largest |S| per equation and the largest change of a
// nest's face velocities.

#include "Domain.H"

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
};

//! Run the exchanges: members[0] is the outer domain (with its defect source),
//! the others its nests (with their parent images and faces); first_exchange: the exchanges an earlier run took
//! (a joint restart; 0: none), so this run takes at most max_exchanges - first_exchange.
JointResult run_joint (const std::vector<JointMember*>& members, int max_exchanges, int first_exchange = 0);

} // namespace gust

#endif
