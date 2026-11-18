#include "App.H"

#include "Domain.H"
#include "Driver.H"
#include "JointRecord.H"
#include "JointRun.H"
#include "NestDefect.H"
#include "NestInject.H"
#include "NestRelax.H"
#include "Nests.H"

#include <AMReX.H>
#include <AMReX_FileSystem.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <array>
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
    // nest.feedback = defect, relax or inject (JointRun.H): the outer domain, once steady, runs each nest inside its own run_domain
    // (DomainRole::joint), and the last nest runs the joint loop, so every domain stays in memory.
    std::vector<JointMember*> members;
    std::vector<OuterOutcome> nest_outcomes(nests.nests.size());
    std::vector<std::string> nest_chks(nests.nests.size());
    // A joint restart (JointRecord.H): io.restart holds the record of an earlier joint run. Without one, a restart
    // continues the outer domain's one-way run and the nests start from its checkpoint, as in a first run.
    JointRecord prior;
    bool resume = false;
    std::vector<DomainRole::JointResume> resumes(nests.nests.size() + 1);
    if (nests.active() && !joint) {
        // A joint checkpoint restarted one way: the outer domain goes on alone, its history already holding the joint
        // run's iterations.
        std::string restart;
        amrex::ParmParse("io").query("restart", restart);
        if (!restart.empty() && read_joint_record(restart, prior)) {
            amrex::Print() << "gust: io.restart = " << restart << ": a joint checkpoint (" << prior.exchanges
                           << " exchanges with nest.feedback = " << prior.feedback << "); with nest.feedback = none "
                           "the outer domain goes on one way and the nests start again from its checkpoint\n";
            prior = JointRecord{};
        }
    }
    if (joint) {
        std::string restart;
        amrex::ParmParse("io").query("restart", restart);
        if (!restart.empty() && read_joint_record(restart, prior)) {
            resume = true;
            std::string names, recorded;
            for (const auto& s : nests.nests) { names += " " + s.name; }
            for (const auto& n : prior.nests) { recorded += " " + n.name; }
            if (names != recorded) {
                amrex::Abort("gust: io.restart = " + restart + ": its joint record has the nests" + recorded
                             + ", nest.names has" + names + "; a joint restart keeps the nests and their order");
            }
            int max_iterations = SteadyOptions{}.max_iterations;
            amrex::ParmParse("steady").query("max_iterations", max_iterations);
            if (max_iterations <= prior.exchanges) {
                amrex::Abort("gust: io.restart = " + restart + ": the joint run took " + std::to_string(prior.exchanges)
                             + " exchanges; steady.max_iterations = " + std::to_string(max_iterations)
                             + " leaves none (the joint run counts its exchanges over the run)");
            }
            // Every recorded nest checkpoint is there before anything runs (the record names them relative to the
            // directory the joint run started in).
            std::string missing;
            if (amrex::ParallelDescriptor::IOProcessor()) {
                for (const auto& n : prior.nests) {
                    if (!amrex::FileSystem::Exists(n.chk + "/Header")) { missing += " " + n.chk; }
                }
            }
            amrex::Vector<char> mbuf(missing.begin(), missing.end());
            amrex::BroadcastArray(mbuf, amrex::ParallelDescriptor::MyProc(), amrex::ParallelDescriptor::IOProcessorNumber(),
                                  amrex::ParallelDescriptor::Communicator());
            missing.assign(mbuf.begin(), mbuf.end());
            if (!missing.empty()) {
                amrex::Abort("gust: io.restart = " + restart + ": its joint record names nest checkpoints that are not"
                             " here:" + missing + " (the record names them relative to the directory the joint run"
                             " started in; restart from that directory)");
            }
            const std::string fbk = nest_feedback_name(nests.feedback);
            if (fbk != prior.feedback) {
                amrex::Print() << "gust: io.restart = " << restart << ": the joint run took " << prior.exchanges
                               << " exchanges with nest.feedback = " << prior.feedback << "; it goes on with " << fbk << "\n";
            }
            resumes[0] = {restart, prior.outer_cfl_floor};
            for (std::size_t i = 0; i < prior.nests.size(); ++i) {
                resumes[i + 1] = {prior.nests[i].chk, prior.nests[i].cfl_floor};
            }
            outer.resume = &resumes[0];
        }
    }
    JointResult joint_result;
    // The CFL floors of the joint run's drivers, kept for its record. They are read here, right after the joint run:
    // a nest's member lives inside its run_domain, which has returned by the time chain(0) does.
    std::vector<amrex::Real> floors;
    // The outer domain's feedback at the last exchange (S or the target q*), kept for the checkpoint as the floors are.
    std::array<amrex::MultiFab, 5> last_feedback;
    std::function<void(std::size_t)> chain;
    chain = [&] (std::size_t i) {
        if (i == nests.nests.size()) {
            JointMember& P = *members[0];
            const auto fields = [&P] {
                return P.inject != nullptr ? P.inject->feedback_fields()
                       : P.relax != nullptr ? P.relax->feedback_fields() : P.defect->feedback_fields();
            };
            std::vector<amrex::Real> prior_window;
            if (resume) {
                // The stop test's window and the last feedback (JointRecord.H), when the record has them for this mode.
                const std::string fbk = nest_feedback_name(nests.feedback);
                if (prior.has_feedback && prior.feedback == fbk) {
                    std::array<amrex::MultiFab, 5> f;
                    read_joint_feedback(resumes[0].chk, fields(), f);
                    const std::array<const amrex::MultiFab*, 5> pf{&f[0], &f[1], &f[2], &f[3], &f[4]};
                    if (P.inject != nullptr) {
                        P.inject->restore_previous(pf);
                    } else if (P.relax != nullptr) {
                        P.relax->restore_previous(pf);
                    } else {
                        P.defect->restore_previous(pf);
                    }
                    prior_window = prior.window;
                } else {
                    amrex::Print() << "gust: io.restart = " << resumes[0].chk << ": "
                                   << (prior.has_feedback ? "its joint record holds the feedback of nest.feedback = "
                                                            + prior.feedback
                                                          : std::string("its joint record (GUST_JOINT 1) holds no feedback"))
                                   << ", so the joint stop test's window (steady.monitor_window exchanges) fills again: "
                                      "the run may stop later than it would have without the restart\n";
                }
            }
            joint_result = run_joint(members, P.max_iterations, resume ? prior.exchanges : 0, prior_window);
            for (const JointMember* j : members) { floors.push_back(j->cfl_floor); }
            const auto f = fields();
            for (std::size_t n = 0; n < 5; ++n) {
                last_feedback[n].define(f[n]->boxArray(), f[n]->DistributionMap(), 1, 0);
                amrex::MultiFab::Copy(last_feedback[n], *f[n], 0, 0, 1, 0);
            }
            return;
        }
        const NestSpec& s = nests.nests[i];
        amrex::Print() << "gust: nest." << s.name << ": " << s.lo[0] << " " << s.lo[1] << " to " << s.hi[0] << " " << s.hi[1]
                       << " m, " << (resume ? "from its joint checkpoint " + resumes[i + 1].chk + " and" : std::string("from"))
                       << " the outer domain's checkpoint " << members[0]->chk << " (joint run)\n";
        DomainRole r;
        r.nest = &s;
        r.parent_chk = members[0]->chk;
        r.from_outer = &outcome;
        r.outcome = &nest_outcomes[i];
        if (resume) { r.resume = &resumes[i + 1]; }
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
            // A joint restart joins whatever its verdict: the checkpoint is a joint one.
            if (!resume && m.verdict.status != SteadyStatus::Steady) { return; }   // reported after run_domain returns
            members.push_back(&m);
            // The outer run's files are complete on disk before a nest reads them.
            amrex::ParallelDescriptor::Barrier();
            chain(0);
        };
    }
    const std::string chk = run_domain(outer);
    if (joint && !members.empty() && joint_result.exchanges > 0) {
        // The record of the joint run in the outer domain's last checkpoint, beside the nests' (written by then).
        JointRecord rec;
        rec.feedback = nest_feedback_name(nests.feedback);
        rec.exchanges = joint_result.exchanges;
        rec.outer_cfl_floor = floors.at(0);
        rec.window = joint_result.window;
        for (std::size_t i = 0; i < nests.nests.size(); ++i) {
            rec.nests.push_back({nests.nests[i].name, floors.at(i + 1), nest_chks[i]});
        }
        write_joint_record(chk, rec);
        write_joint_feedback(chk, {&last_feedback[0], &last_feedback[1], &last_feedback[2], &last_feedback[3],
                                   &last_feedback[4]});
    }
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
        // The blend is refused with io.restart (App.cpp): the restart then runs without it.
        const std::string blend = amrex::ParmParse("convection").contains("hybrid_sigma_max")
            ? " and without convection.hybrid_sigma_max (the blend is not in the checkpoint)" : "";
        if (not_steady.empty() && joint_result.waiting_on_feedback) {
            // Every domain's own verdict is steady, but the feedback between them had not settled (JointRun.H).
            char what[300];
            std::snprintf(what, sizeof(what), "last change %.3e, sum over the window %.3e",
                          joint_result.feedback_change, joint_result.feedback_window);
            amrex::Abort(std::string("gust: nest.feedback = ") + nest_feedback_name(nests.feedback)
                         + ": the joint run did not reach steady: every domain is steady, but the "
                         + (nests.feedback == NestFeedback::relax    ? "relaxation's target"
                            : nests.feedback == NestFeedback::inject ? "injection's target"
                                                                     : "defect source")
                         + " had not settled over steady.monitor_window exchanges to steady.monitor_rel_tol (" + what
                         + ") within steady.max_iterations = " + std::to_string(joint_result.exchanges) + " exchanges; "
                         "to go on, restart with io.restart = " + chk + " (its joint record names the nests' "
                         "checkpoints) and a larger steady.max_iterations" + blend);
        }
        if (!not_steady.empty()) {
            amrex::Abort(std::string("gust: nest.feedback = ") + nest_feedback_name(nests.feedback)
                         + ": the joint run did not reach steady after " + std::to_string(joint_result.exchanges)
                         + " exchanges (steady.max_iterations, steady.residual_tol, steady.residual_drop_orders); "
                         "to go on, restart with io.restart = " + chk + " (its joint record names the nests' checkpoints), "
                         "a larger steady.max_iterations" + blend + ". Not steady:" + not_steady);
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
