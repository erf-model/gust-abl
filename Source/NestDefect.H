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
#include <deque>
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

//! One nest for restrict_patches.
struct NestImage
{
    const State* fine = nullptr;
    const Mesh* mesh = nullptr;
    const amrex::Geometry* geom = nullptr;
    int ratio = 2;
    amrex::IntVect offset;
};

//! restrict_patch for several nests with two parent ghost fills in all: the
//! cell fields, u and v of every nest; coarse_ghosts; w of every nest;
//! coarse_ghosts. w reads the parent's u and v (and their ghosts), so it comes
//! after every nest's u and v: on a face two footprints share, w on both
//! sides sees the face's final flux, and the parent's Omega on every
//! footprint is the restricted nest Omega. That face keeps the restriction of
//! the later nest in the list; the cells next to it on the earlier nest's
//! side are then not divergence-free by the difference of the two
//! restrictions. restrict_patch is this with one nest.
void restrict_patches (const std::vector<NestImage>& nests, const Mesh& mesh_c, const amrex::Geometry& geom_c,
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

//! How much a feedback changed at its last set (the defect source S, NestDefect::set, or the relaxation's
//! target q*, NestRelax::set), per field u, v, w, theta, k: change = max |new - old| over the points the feedback
//! acts on (a NaN counts as infinite), size = the new field's size there (each class says which). first: no
//! earlier set to compare with. Every rank holds the same values.
struct FeedbackChange
{
    bool first{true};
    std::array<amrex::Real, 5> change{};
    std::array<amrex::Real, 5> size{};
};

//! The joint run's measure of a feedback's change between two exchanges (JointRun.H): the largest over the fields
//! of change / max(size, floor) (0 / 0 = 0, a positive change over 0 is infinite); infinite when c.first. floor:
//! the size below which a change is measured absolutely (as SteadyCriterion's monitor scales).
amrex::Real relative_feedback_change (const FeedbackChange& c, const std::array<amrex::Real, 5>& floor);

//! The joint run's test that a feedback has settled (JointRun.H), on the per-exchange changes d_n of
//! relative_feedback_change, as SteadyCriterion tests a monitor over steady.monitor_window iterations (one exchange
//! is one iteration): settled when the window holds W changes, their sum (which bounds the feedback's relative
//! change over the window) is at most tol, and Aitken's estimate of the change still to come, d2 rho / (1 - rho)
//! with d1, d2 the sums over the window's two halves and rho = d2 / d1, is at most tol (skipped when both sums are
//! below 1e-3 tol; refused when rho >= 1).
class FeedbackWindow
{
public:
    FeedbackWindow (int window, amrex::Real tol) : m_window(window), m_tol(tol) {}
    //! Adds an exchange's change (infinite, NaN or negative empties the window: nothing before it counts).
    void add (amrex::Real d);
    //! The sum of the changes in the window (infinite until it holds W).
    amrex::Real sum () const;
    bool settled () const;
    //! The changes in the window, oldest first (at most W; a joint checkpoint stores them).
    [[nodiscard]] std::vector<amrex::Real> values () const { return {m_d.begin(), m_d.end()}; }
    //! A joint restart: the window as values() gave it before the checkpoint (each added in order).
    void seed (const std::vector<amrex::Real>& d) { for (const amrex::Real x : d) { add(x); } }
private:
    int m_window;
    amrex::Real m_tol;
    std::deque<amrex::Real> m_d;
};

class NestDefect
{
public:
    //! The source on the parent's layout, zero and inactive until set.
    NestDefect (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    //! S = R on the points of the regions' inner boxes (cells; faces from L to
    //! H + 1 along their normal), 0 elsewhere; R is the parent's residual at
    //! q* without the source. Sets the source active. Records the change of
    //! S (change(): max |S_new - S_old| over the inner boxes' points, size the
    //! largest |S_new|; first on the first set).
    void set (const ResidualFields& R, const std::vector<DefectRegion>& regions);

    //! R -= S for u, v, w, theta and k (nothing while inactive).
    void add_to_residual (ResidualFields& R) const;

    [[nodiscard]] bool active () const { return m_active; }
    //! The largest |S| of u, v, w, theta and k over the parent (after set).
    [[nodiscard]] std::vector<amrex::Real> max_norms () const;
    [[nodiscard]] const ResidualFields& source () const { return m_S; }
    //! The change of S at the last set.
    [[nodiscard]] const FeedbackChange& change () const { return m_change; }
    //! S of u, v, w, theta and k (a joint checkpoint stores them: JointRecord.H).
    [[nodiscard]] std::array<const amrex::MultiFab*, 5> feedback_fields () const
    {
        return {&m_S.u, &m_S.v, &m_S.w, &m_S.theta, &m_S.tke};
    }
    //! A joint restart: S of the last exchange before the checkpoint (valid points), so that the first set
    //! measures its change as the run without a restart does. The source stays inactive until set.
    void restore_previous (const std::array<const amrex::MultiFab*, 5>& f);

private:
    ResidualFields m_S;
    FeedbackChange m_change;
    bool m_active{false};
    bool m_have_previous{false};
};

} // namespace gust

#endif
