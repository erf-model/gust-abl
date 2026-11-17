#ifndef GUST_NEST_DEFECT_H_
#define GUST_NEST_DEFECT_H_

// The patches' feedback on the parent in the joint two-level run
// (nest.feedback = defect; plan rev. 6 tracker, PR-70; the relaxation of
// nest.feedback = relax, PR-71, is NestRelax.H and reuses restrict_patch,
// DefectRegion and kDefectMargin from here). The parent and its
// nests iterate together in memory (JointRun.H); each exchange:
//
//   down: each nest's faces take the current parent (NestFaces::refresh);
//   up:   the parent's residual gets the defect source of every nest,
//
//     R_c(q_c) - S = 0,    S = chi_s R_c(q*),    q* = q_c with I q_f on the nests' footprints,
//
// with
//   R_c      the parent's residual (Residual.H), with the stage's upwind blend;
//   q_c      the parent's state;
//   q_f      a nest's state; I q_f its restriction onto the parent cells under
//            the nest (restrict_patch: cells the detJ-weighted mean, u and v
//            the face flux over the coarse area, w from the restricted
//            vertical flux Omega);
//   chi_s    1 on the parent's points at least kDefectMargin parent cells
//            inside the footprint's lateral faces (all z), 0 elsewhere. The
//            ring keeps the nest's face values out of S's own cells; the
//            upwind stencil at chi_s's edge still reads the face next to it.
// On chi_s the parent's equations become R_c(q_c) = R_c(q*): there the parent
// carries the nests' solution, and outside the footprint the parent's own
// equations see it through their stencils. q* takes I q_f on the ring too, so
// R_c(q*) on chi_s does not mix the nest's state with the parent's (with q*
// = I q_f on chi_s only, the Gorge test case settles into a cycle near a
// residual of 1e-4 and does not reach steady). At the joint steady state q_c
// on the ring is the parent's own, close to but not exactly I q_f. S is a constant
// of the parent's iteration (no implicit term; it changes only at an
// exchange). Continuity gets no source: I q_f is divergence-free.
//
// The rows a boundary condition or the wall law replaces (the wall-k rows, the
// Dirichlet rows) replace S too: the source is added right after the time
// term (Residual.cpp), before them.
//
// The ground: a nest keeps the parent's z-levels; its ground can have more
// detail. The restriction goes by zeta index (as NestParent's interpolation),
// with the fluxes conservative across the different ground.
//
// Limit: with turb.nut_relax < 1 the closure's nu_t blends with the parent's
// own history (Closure::accept), so R_c(q*) uses that history, not one of
// q*'s own; at a steady state the history is nu_t of q_c, so only the path
// changes.

#include "Mesh.H"
#include "Residual.H"
#include "State.H"

#include <array>
#include <AMReX_Box.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <functional>
#include <string>
#include <vector>

namespace gust {

//! Parent cells between the footprint's lateral faces and the defect source's region.
inline constexpr int kDefectMargin = 2;

//! I q_f: the nest state fine (on the layout of mesh_f; u and v with one
//! ghost cell filled, as BoundaryConditions::fill_ghosts leaves them) restricted
//! at r r 1 onto the parent cells from offset (parent cells, z: 0) in the
//! parent state coarse (mesh_c's layout). Only the footprint's valid points
//! change; outside them coarse keeps its values. w on the footprint's
//! z-faces is w from the restricted vertical flux Omega; coarse_ghosts then
//! resets the ground and top faces to their boundary values (it refills
//! them, and the ghosts, after u and v and again at the end).
void restrict_patch (const State& fine, const Mesh& mesh_f, const amrex::Geometry& geom_f, int r,
                     const amrex::IntVect& offset, const Mesh& mesh_c, const amrex::Geometry& geom_c,
                     const std::function<void(State&)>& coarse_ghosts, State& coarse);

//! The footprint of a nest on the parent (cells) and the source's region
//! inside it (kDefectMargin cells in from each lateral face, all z).
struct DefectRegion
{
    std::string name;
    amrex::Box footprint;   //!< parent cells under the nest
    amrex::Box inner;       //!< footprint grown by -kDefectMargin in x and y (not in a periodic direction the nest spans)
};

//! The region of a nest with fine_cells cells at ratio r r 1 from offset in a
//! parent domain dom_c; name is the domain's name in messages ("nest.<name>"). A nest that spans a periodic parent
//! direction has no faces there (periodic[d]): the inner region keeps the whole span. Aborts if the inner region is
//! empty.
DefectRegion defect_region (const std::string& name, const amrex::IntVect& fine_cells, int r, const amrex::IntVect& offset,
                            const amrex::Box& dom_c, std::array<bool, 2> periodic = {false, false});

class NestDefect
{
public:
    //! The source on the parent's layout, zero and inactive until set.
    NestDefect (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    //! S = R on the points of the regions' inner boxes (cells; faces from L to
    //! H + 1 along their normal), 0 elsewhere; R is the parent's residual at
    //! q* without the source. Sets the source active.
    void set (const ResidualFields& R, const std::vector<DefectRegion>& regions);

    //! R -= S for u, v, w, theta and k (nothing while inactive).
    void add_to_residual (ResidualFields& R) const;

    [[nodiscard]] bool active () const { return m_active; }
    //! The largest |S| of u, v, w, theta and k over the parent (after set).
    [[nodiscard]] std::vector<amrex::Real> max_norms () const;
    [[nodiscard]] const ResidualFields& source () const { return m_S; }

private:
    ResidualFields m_S;
    bool m_active{false};
};

} // namespace gust

#endif
