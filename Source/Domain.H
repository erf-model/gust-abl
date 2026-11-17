#ifndef GUST_DOMAIN_H_
#define GUST_DOMAIN_H_

// One domain's run (App.cpp run_domain) and the role it runs in: the outer
// domain, or a nest of it (Nests.H). run() (Run.cpp) runs the outer domain,
// then each nest; with nest.feedback = defect or relax, then the outer domain and the
// nests together (JointRun.H).

#include "Nests.H"
#include "SteadyCriterion.H"

#include <AMReX_Geometry.H>

#include <functional>
#include <string>

namespace gust {

class Mesh;
class NestDefect;
class NestRelax;
class NestFaces;
class NestParent;
struct ResidualFields;
struct State;

//! One domain of the joint run (nest.feedback = defect or relax, JointRun.H), handed
//! out by run_domain once its one-way run is done; its pointers and functions
//! are valid while the hook that receives it runs.
struct JointMember
{
    std::string name;                 //!< "outer" or "nest.<name>"
    State* state{nullptr};
    const Mesh* mesh{nullptr};
    const amrex::Geometry* geom{nullptr};
    std::string chk;                  //!< the checkpoint the one-way run wrote
    SteadyVerdict verdict;            //!< the last step's (the one-way run's at first)
    int iterations{0};                //!< joint iterations taken
    int max_iterations{0};            //!< steady.max_iterations as the run read it
    double wall_seconds{0};           //!< their wall time
    //! This domain's ghost-cell fill (BoundaryConditions::fill_ghosts) of a state on its layout.
    std::function<void(State&)> fill_ghosts;
    //! One iteration of this domain's driver (logged; the step count advances).
    std::function<SteadyVerdict()> step;
    // The outer domain: its residual without the source (with the stage's upwind blend), and the source.
    std::function<void(const State&, ResidualFields&)> residual_without_defect;
    NestDefect* defect{nullptr};
    //! The outer domain with nest.feedback = relax: the relaxation (NestRelax.H) in place of the defect source.
    NestRelax* relax{nullptr};
    // A nest: its parent's images and its faces.
    NestParent* parent{nullptr};
    NestFaces* faces{nullptr};
};

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
    //! nest.feedback = defect or relax: called with this domain once its one-way run is done and its outputs
    //! are written (not if it diverged); after it returns a domain that took joint iterations
    //! writes its outputs again. Empty: one way.
    std::function<void(JointMember&)> joint;
};

//! Read and check the inputs, build the mesh and state (fresh or from
//! io.restart), run the stage, and write the plotfile and checkpoint of one
//! domain; returns the name of the last checkpoint written (empty: none).
std::string run_domain (const DomainRole& role);

} // namespace gust

#endif
