#include "App.H"

#include "Domain.H"
#include "JointRun.H"
#include "Nests.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <functional>
#include <string>
#include <vector>

namespace gust {

void run ()
{
    // The nests (Nests.H): read first, so that the outer domain's unknown-key check sees their keys read.
    const NestOptions nests = read_nest_options();
    DomainRole outer;
    if (nests.active() || !nests.errors.empty()) {
        outer.nests = &nests;
        // Every nest's derived inputs are checked before the outer domain runs.
        outer.after_checks = [&nests] {
            for (const auto& s : nests.nests) {
                DomainRole r;
                r.nest = &s;
                r.parent_chk = "(the outer domain's checkpoint)";
                r.check_only = true;
                run_domain(r);
            }
        };
    }
    OuterOutcome outcome;
    if (nests.active()) { outer.outcome = &outcome; }
    const bool joint = nests.active() && nests.feedback != NestFeedback::none;
    // nest.feedback = defect or relax (JointRun.H): the outer domain, once steady, runs each nest inside its own run_domain
    // (DomainRole::joint), and the last nest runs the joint loop, so every domain stays in memory.
    std::vector<JointMember*> members;
    std::vector<OuterOutcome> nest_outcomes(nests.nests.size());
    std::vector<std::string> nest_chks(nests.nests.size());
    std::function<void(std::size_t)> chain;
    chain = [&] (std::size_t i) {
        if (i == nests.nests.size()) {
            run_joint(members, members[0]->max_iterations);
            return;
        }
        const NestSpec& s = nests.nests[i];
        amrex::Print() << "gust: nest." << s.name << ": " << s.lo[0] << " " << s.lo[1] << " to " << s.hi[0] << " " << s.hi[1]
                       << " m, from the outer domain's checkpoint " << members[0]->chk << " (joint run)\n";
        DomainRole r;
        r.nest = &s;
        r.parent_chk = members[0]->chk;
        r.from_outer = &outcome;
        r.outcome = &nest_outcomes[i];
        r.joint = [&chain, &members, i] (JointMember& m) {
            members.push_back(&m);
            chain(i + 1);
        };
        nest_chks[i] = run_domain(r);
        if (members.size() < i + 2) {
            // finish() stops a diverged nest before it joins, so this is an internal check.
            amrex::Abort("gust: nest." + s.name + ": internal: the nest did not join the joint run");
        }
    };
    if (joint) {
        outer.joint = [&] (JointMember& m) {
            if (m.verdict.status != SteadyStatus::Steady) { return; }   // reported after run_domain returns
            members.push_back(&m);
            // The outer run's files are complete on disk before a nest reads them.
            amrex::ParallelDescriptor::Barrier();
            chain(0);
        };
    }
    const std::string chk = run_domain(outer);
    if (nests.active() && !outcome.steady && members.empty()) {
        amrex::Abort("gust: nest.names: the outer domain did not reach steady (steady.max_iterations, steady.residual_tol, "
                     "steady.residual_drop_orders), "
                     "so the nests are not run; its last checkpoint is " + chk);
    }
    std::string not_steady;
    if (joint) {
        for (std::size_t i = 0; i < nests.nests.size(); ++i) {
            if (!nest_outcomes[i].steady) { not_steady += "\n  nest." + nests.nests[i].name + " (last checkpoint " + nest_chks[i] + ")"; }
        }
        if (!outcome.steady) { not_steady += "\n  the outer domain (last checkpoint " + chk + ")"; }
        if (!not_steady.empty()) {
            amrex::Abort(std::string("gust: nest.feedback = ") + nest_feedback_name(nests.feedback)
                         + ": the joint run did not reach steady (steady.max_iterations, "
                         "steady.residual_tol, steady.residual_drop_orders); a joint run does not restart, so rerun "
                         "it from the start with other settings. Not steady:" + not_steady);
        }
        return;
    }
    // The outer run's files are complete on disk before a nest reads them.
    amrex::ParallelDescriptor::Barrier();
    for (const auto& s : nests.nests) {
        if (chk.empty()) { amrex::Abort("gust: nest." + s.name + ": the outer domain wrote no checkpoint"); }
        amrex::Print() << "gust: nest." << s.name << ": " << s.lo[0] << " " << s.lo[1] << " to " << s.hi[0] << " " << s.hi[1]
                       << " m, from the outer domain's checkpoint " << chk << "\n";
        DomainRole r;
        r.nest = &s;
        r.parent_chk = chk;
        r.from_outer = &outcome;
        OuterOutcome nest_outcome;
        r.outcome = &nest_outcome;
        const std::string nest_chk = run_domain(r);
        if (!nest_outcome.steady) { not_steady += "\n  nest." + s.name + " (last checkpoint " + nest_chk + ")"; }
    }
    // Every nest runs; then the run stops with an error if any did not reach steady.
    if (!not_steady.empty()) {
        amrex::Abort("gust: nest.names: these nests did not reach steady (steady.max_iterations, steady.residual_tol, "
                     "steady.residual_drop_orders):" + not_steady);
    }
}

} // namespace gust
