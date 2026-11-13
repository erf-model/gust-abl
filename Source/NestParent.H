#ifndef GUST_NEST_PARENT_H_
#define GUST_NEST_PARENT_H_

// The parent run of a nested patch: a smaller, finer domain inside a
// finished (frozen) parent run, one way (the patch does not change the
// parent). The parent is read from its checkpoint (bc.nest_parent), not its
// plotfile: the plotfile holds cell-centred velocities only, the checkpoint
// the face velocities, theta, k and p'.
//
// Placement (every rule aborts naming bc.nest_parent):
//  * refinement 2 2 1: the patch's cells are half the parent's in x and y,
//    and its faces lie on parent cell faces, so each parent cell holds
//    2 x 2 patch cells, aligned (to 1e-4 of a parent cell, which covers
//    round-off in single precision with coordinates of order 1e6 m);
//  * the same z-levels: n_cell z, prob_lo z, prob_hi z and the reference
//    z-levels (ZLevels.H; in the checkpoint Header since format version 2)
//    equal bit for bit. A ratio of 8 is a chain of three nests, each run
//    the parent of the next; no level refines in z, so nothing is
//    interpolated between z-levels;
//  * inside the parent: in a non-periodic direction of the parent at least
//    kNestMargin parent cells from its domain faces, so every patch point,
//    its two ghost layers included, lies between parent cell centres and
//    is interpolated, never extrapolated (k and theta stay within the
//    parent's range); in a periodic direction anywhere inside the domain
//    (the parent's periodic images give the neighbours), and a periodic
//    patch must span the parent's whole width there.
//
// Interpolation: every field at its own staggered point, per zeta-level:
// horizontal bilinear interpolation of the parent values around it
// (Transfer_K.H, the operator of the grid hierarchy's prolong_fields with
// an offset). A field a + b x + c y + d x y on each zeta-level is reproduced
// to round-off. The patch's node heights are its own (its ground can have
// more detail than the parent's); the values are taken at the same zeta-level
// index, not at the same height (theta is then off by about
// lapse rate x ground difference near the ground; docs/initial_conditions.rst).
// The parent's ground flux is not kept: the zero-gradient ghost row below the
// ground is the patch's, so the caller sets w on the ground faces with
// zero_ground_flux (Initializer.H) before the correction.
//
// Layouts: the parent is read on its own BoxArray, then copied
// (ParallelCopy, with the parent's periodicity) to the parent cells under
// each patch box, on the patch's DistributionMapping, with kNestImageGhost
// ghost cells; the parent's fields are released after the copy. The result
// does not depend on the parent's or the patch's box layout or rank count.

#include "Inputs.H"
#include "State.H"
#include "StateIO.H"
#include "Transfer_K.H"

#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace gust {

//! The refinement of a nest in x and y (z: 1).
inline constexpr int kNestRatio = 2;
//! Parent cells between the patch and a non-periodic parent face.
inline constexpr int kNestMargin = 2;
//! Ghost cells of the parent images: the parent cells around the patch's
//! points and its ngrow_state ghost layers.
inline constexpr int kNestImageGhost = 2;

//! The problems with placing the patch (mesh inputs patch, reference z-levels
//! patch_zlevels) in the parent whose checkpoint Header is h (empty: none);
//! each message starts with "bc.nest_parent". offset, when given and there is
//! no problem, receives the parent cells between the parent's low corner and
//! the patch's in x and y (z: 0).
std::vector<std::string> check_nest_placement (const CheckpointHeader& h, const MeshInputs& patch,
                                               const amrex::Vector<amrex::Real>& patch_zlevels,
                                               amrex::IntVect* offset = nullptr);

class NestParent
{
public:
    //! Read the parent checkpoint dir for a patch with mesh inputs patch,
    //! reference z-levels patch_zlevels and layout (patch_ba, patch_dm).
    //! Aborts, naming bc.nest_parent, if dir has no Header, the Header is
    //! not a GUST checkpoint of format version 2 or later, or the placement
    //! breaks a rule above.
    NestParent (const std::string& dir, const MeshInputs& patch, const amrex::Vector<amrex::Real>& patch_zlevels,
                const amrex::BoxArray& patch_ba, const amrex::DistributionMapping& patch_dm);

    //! Set u, v, w, p', theta and k of patch (defined on the layout given to
    //! the constructor) from the parent: at the valid points and the
    //! ngrow_state ghost layers in x and y (interpolated as the valid
    //! points: neighbouring boxes agree bit for bit, and a periodic seam
    //! gets its images), then the ghost rows below the ground and above the
    //! top from the nearest valid row (zero gradient), so u and v are ready
    //! for the divergence. BoundaryConditions::fill_ghosts later sets the
    //! domain faces' ghost cells. z_phys_nd is not touched.
    void interpolate (State& patch) const;

    [[nodiscard]] const amrex::IntVect& offset () const { return m_offset; }
    [[nodiscard]] int parent_step () const { return m_step; }
    [[nodiscard]] const amrex::IntVect& parent_cells () const { return m_parent_cells; }
    //! One line for the log: the parent, its step and cells, the patch's place in it.
    [[nodiscard]] std::string summary () const;

private:
    std::string m_dir;
    int m_step{0};
    amrex::IntVect m_offset{0, 0, 0};
    amrex::IntVect m_parent_cells{0, 0, 0};
    amrex::IntVect m_patch_cells{0, 0, 0};
    amrex::BoxArray m_patch_ba;
    // The parent's values under the patch's boxes (parent indices).
    amrex::MultiFab m_u, m_v, m_w, m_p, m_theta, m_tke;
    transfer::Axis m_cx, m_cy, m_nx, m_ny;
};

} // namespace gust

#endif
