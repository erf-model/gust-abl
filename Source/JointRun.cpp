#include "JointRun.H"

#include "Mesh.H"
#include "NestDefect.H"
#include "NestFaces.H"
#include "NestParent.H"
#include "NestRelax.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>
#include <AMReX_Utility.H>

#include <algorithm>
#include <cstdio>
#include <string>

namespace gust {

namespace {

//! The solution fields and the step; the node heights do not change in the joint run and are copied once.
void copy_state (const State& src, State& dst)
{
    const amrex::MultiFab* s[6] = {&src.u, &src.v, &src.w, &src.p, &src.theta, &src.tke};
    amrex::MultiFab* d[6] = {&dst.u, &dst.v, &dst.w, &dst.p, &dst.theta, &dst.tke};
    for (int n = 0; n < 6; ++n) { amrex::MultiFab::Copy(*d[n], *s[n], 0, 0, 1, s[n]->nGrow()); }
    dst.step = src.step;
}

const char* status_name (SteadyStatus s)
{
    return s == SteadyStatus::Steady ? "steady" : s == SteadyStatus::Diverged ? "diverged" : "not steady";
}

} // namespace

std::string joint_log_name (const std::string& name)
{
    constexpr std::size_t width = kJointLogNameWidth, head = 12, tail = 9;
    static_assert(head + 3 + tail == width, "joint_log_name keeps head + \"...\" + tail characters");
    if (name.size() <= width) { return name; }
    return name.substr(0, head) + "..." + name.substr(name.size() - tail);
}

JointResult run_joint (const std::vector<JointMember*>& members, int max_exchanges, int first_exchange,
                       const std::vector<amrex::Real>& prior_window)
{
    if (first_exchange < 0 || first_exchange >= max_exchanges) {
        amrex::Abort("gust: run_joint: " + std::to_string(first_exchange) + " exchanges taken of at most "
                     + std::to_string(max_exchanges) + " (steady.max_iterations)");
    }
    if (members.size() < 2 || (members[0]->defect == nullptr) == (members[0]->relax == nullptr)) {
        amrex::Abort("gust: run_joint: needs the outer domain with its defect source or its relaxation (one of them) and "
                     "at least one nest");
    }
    JointMember& P = *members[0];
    const bool relax = P.relax != nullptr;
    std::vector<DefectRegion> regions;
    for (std::size_t n = 1; n < members.size(); ++n) {
        const JointMember& N = *members[n];
        if (N.parent == nullptr || N.faces == nullptr) {
            amrex::Abort("gust: run_joint: " + N.name + " has no parent images or nest faces");
        }
        regions.push_back(defect_region(N.name, N.geom->Domain().size(), N.parent->ratio(), N.parent->offset(),
                                        P.geom->Domain(), {N.geom->isPeriodic(0), N.geom->isPeriodic(1)}));
    }
    if (relax) {
        amrex::Print() << "gust: joint run (nest.feedback = relax, nest.feedback_factor = " << P.relax->factor()
                       << "): the outer domain and " << members.size() - 1 << " nest(s) exchange every iteration, up to "
                       << max_exchanges << " exchanges" << (first_exchange > 0 ? " over the run" : "")
                       << "; the relaxation's weight rises from 0 on each nest's faces to 1 "
                       << kDefectMargin << " outer cells inside them\n";
    } else {
        amrex::Print() << "gust: joint run (nest.feedback = defect): the outer domain and " << members.size() - 1
                       << " nest(s) exchange every iteration, up to " << max_exchanges << " exchanges"
                       << (first_exchange > 0 ? " over the run" : "") << "; the source acts "
                       << kDefectMargin << " outer cells inside each nest's faces\n";
    }
    State qr;
    qr.define(P.state->p.boxArray(), P.state->p.DistributionMap());
    amrex::MultiFab::Copy(qr.z_phys_nd, P.state->z_phys_nd, 0, 0, 1, P.state->z_phys_nd.nGrow());
    // The nests' ratios and offsets are fixed after the one-way runs.
    std::vector<NestImage> images;
    for (std::size_t n = 1; n < members.size(); ++n) {
        const JointMember& N = *members[n];
        images.push_back({N.state, N.mesh, N.geom, N.parent->ratio(), N.parent->offset()});
    }
    ResidualFields R;
    JointResult out;
    FeedbackWindow window(P.feedback_window, P.feedback_tol);
    window.seed(prior_window);
    const double t_all = amrex::second();
    if (first_exchange > 0) {
        amrex::Print() << "gust: joint run: restarted after exchange " << first_exchange << "\n";
    }
    for (int it = first_exchange; it < max_exchanges; ++it) {
        // 1. Up: the nests restricted onto the outer state's copy, the source from its residual.
        double t0 = amrex::second();
        // q* = the outer state with I q_f on every footprint (two outer ghost fills for any number of nests).
        copy_state(*P.state, qr);
        for (std::size_t n = 1; n < members.size(); ++n) { members[n]->fill_ghosts(*members[n]->state); }
        restrict_patches(images, *P.mesh, *P.geom, P.fill_ghosts, qr);
        if (relax) {
            P.relax->set(qr, regions);
            if (it == first_exchange) {
                for (const auto& n : P.relax->nests()) {
                    char line[200];
                    std::snprintf(line, sizeof(line), "gust: joint: %s: mean horizontal speed %.4e m/s, rate f U / dx "
                                  "%.4e 1/s\n", joint_log_name(n.name).c_str(), double(n.speed), double(n.rate));
                    amrex::Print() << line;
                }
            }
        } else {
            P.residual_without_defect(qr, R);
            P.defect->set(R, regions);
        }
        // The feedback's change from the last exchange (infinite at the first).
        const amrex::Real fchange = relative_feedback_change(relax ? P.relax->change() : P.defect->change(), P.feedback_floor);
        window.add(fchange);
        double tx = amrex::second() - t0;
        // 2. The outer domain's iteration.
        t0 = amrex::second();
        P.verdict = P.step();
        ++P.iterations;
        P.wall_seconds += amrex::second() - t0;
        // 3. Down: every nest's faces from the outer state now.
        t0 = amrex::second();
        amrex::Real change = 0;
        for (std::size_t n = 1; n < members.size(); ++n) {
            JointMember& N = *members[n];
            N.parent->refresh(*P.state);
            change = std::max(change, N.faces->refresh(*N.parent, *N.mesh));
            // The ghosts from the new faces: the driver's local step reads them before the iteration fills them, and
            // a joint restart's first step fills them from these faces (Stage::init), so the run and its restart agree.
            N.fill_ghosts(*N.state);
        }
        tx += amrex::second() - t0;
        // 4. The nests' iterations.
        bool all_steady = P.verdict.status == SteadyStatus::Steady;
        bool diverged = P.verdict.status == SteadyStatus::Diverged;
        for (std::size_t n = 1; n < members.size(); ++n) {
            JointMember& N = *members[n];
            t0 = amrex::second();
            N.verdict = N.step();
            ++N.iterations;
            N.wall_seconds += amrex::second() - t0;
            all_steady = all_steady && N.verdict.status == SteadyStatus::Steady;
            diverged = diverged || N.verdict.status == SteadyStatus::Diverged;
        }
        // The exchange's time is the slowest rank's.
        amrex::ParallelDescriptor::ReduceRealMax(tx);
        out.exchange_seconds += tx;
        out.exchanges = it + 1;
        out.feedback_change = static_cast<double>(fchange);
        out.feedback_window = static_cast<double>(window.sum());
        // Steady: every domain's verdict, and the feedback they were stepped with has settled over the window.
        const bool settled = window.settled();
        const bool last = (all_steady && settled) || diverged || it + 1 == max_exchanges;
        if (it < first_exchange + 10 || it % 10 == 9 || last) {
            char line[512];
            if (relax) {
                // The outer state after its step against the target of this exchange, where the rate is > 0.
                const auto d = P.relax->differences(*P.state);
                std::snprintf(line, sizeof(line), "gust: joint %5d: outer - nest under the nests: max u %.3e v %.3e "
                              "w %.3e theta %.3e k %.3e, RMS u %.3e v %.3e w %.3e theta %.3e k %.3e; face change %.3e m/s; "
                              "target change %.3e (window %.3e); outer %s", it + 1, double(d[0]), double(d[1]), double(d[2]), double(d[3]),
                              double(d[4]), double(d[5]), double(d[6]), double(d[7]), double(d[8]), double(d[9]),
                              double(change), double(fchange), double(window.sum()), status_name(P.verdict.status));
            } else {
                const auto s = P.defect->max_norms();
                std::snprintf(line, sizeof(line), "gust: joint %5d: max |S| u %.3e v %.3e w %.3e theta %.3e k %.3e; S change "
                              "%.3e (window %.3e); face change %.3e m/s; outer %s", it + 1, double(s[0]), double(s[1]),
                              double(s[2]), double(s[3]), double(s[4]), double(fchange), double(window.sum()),
                              double(change), status_name(P.verdict.status));
            }
            std::string text = line;
            for (std::size_t n = 1; n < members.size(); ++n) {
                text += ", " + joint_log_name(members[n]->name) + " " + status_name(members[n]->verdict.status);
            }
            amrex::Print() << text << "\n";
        }
        if (diverged) { out.diverged = true; break; }
        if (all_steady && settled) { out.steady = true; break; }
        out.waiting_on_feedback = all_steady;
    }
    out.window = window.values();
    out.wall_seconds = amrex::second() - t_all;
    amrex::ParallelDescriptor::ReduceRealMax(out.wall_seconds);
    std::string per;
    for (const JointMember* m : members) {
        char b[160];
        std::snprintf(b, sizeof(b), "; %s %d iterations %.3f s", joint_log_name(m->name).c_str(), m->iterations,
                      m->wall_seconds);
        per += b;
    }
    char line[200];
    std::snprintf(line, sizeof(line), "gust: joint run: %s after %d exchanges, %.3f s (exchanges %.3f s)",
                  out.diverged ? "diverged" : out.steady ? "steady" : "not steady", out.exchanges, out.wall_seconds,
                  out.exchange_seconds);
    amrex::Print() << line << per << "\n";
    if (!out.steady && !out.diverged && out.waiting_on_feedback) {
        char note[400];
        std::snprintf(note, sizeof(note), "gust: joint run: every domain was steady at the last exchange, but the %s had "
                      "not settled: last change %.3e, sum over steady.monitor_window = %d exchanges %.3e, "
                      "steady.monitor_rel_tol = %.3e\n", relax ? "target" : "source S", out.feedback_change,
                      P.feedback_window, out.feedback_window, double(P.feedback_tol));
        amrex::Print() << note;
    }
    return out;
}

} // namespace gust
