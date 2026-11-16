#include "App.H"

#include "Domain.H"
#include "Nests.H"

#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Print.H>

#include <string>

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
    const std::string chk = run_domain(outer);
    if (nests.active() && !outcome.steady) {
        amrex::Abort("gust: nest.names: the outer domain did not reach steady (steady.max_iterations, steady.residual_tol, "
                     "steady.residual_drop_orders), "
                     "so the nests are not run; its last checkpoint is " + chk);
    }
    // The outer run's files are complete on disk before a nest reads them.
    amrex::ParallelDescriptor::Barrier();
    std::string not_steady;
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
