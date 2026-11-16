#ifndef GUST_DOMAIN_H_
#define GUST_DOMAIN_H_

// One domain's run (App.cpp run_domain) and the role it runs in: the outer
// domain, or a nest of it (Nests.H). run() (Run.cpp) runs the outer domain,
// then each nest.

#include "Nests.H"

#include <functional>
#include <string>

namespace gust {

//! Which domain run_domain runs.
struct DomainRole
{
    //! Filled with the run's steady verdict (and, for the outer domain, the forcing its nests take).
    OuterOutcome* outcome{nullptr};
    //! A nest: the outer run's outcome (null in a check-only pass).
    const OuterOutcome* from_outer{nullptr};
    //! The nest (Nests.H); null: the outer domain.
    const NestSpec* nest{nullptr};
    //! A nest's parent checkpoint (the outer run's last); in a check-only
    //! pass a stand-in (nothing is read).
    std::string parent_chk;
    //! The outer domain with nest.names: the nests to check with its inputs.
    const NestOptions* nests{nullptr};
    //! Called once the inputs pass their checks, before anything is built.
    std::function<void()> after_checks;
    //! Return once the inputs pass their checks.
    bool check_only{false};
};

//! Read and check the inputs, build the mesh and state (fresh or from
//! io.restart), run the stage, and write the plotfile and checkpoint of one
//! domain; returns the name of the last checkpoint written (empty: none).
std::string run_domain (const DomainRole& role);

} // namespace gust

#endif
