#ifndef GUST_STATE_H_
#define GUST_STATE_H_

#include "Inputs.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace gust {

//! Ghost cells on the solution fields (the MUSCL reconstruction of
//! Convection.H reads two cells upwind of a face) and on the nodal heights.
inline constexpr int ngrow_state = 2;
inline constexpr int ngrow_nodal = 1;

//! The solution on ERF's Arakawa C-grid: normal velocities on faces, p', theta
//! and k at cell centres, physical heights at nodes.
struct State
{
    amrex::MultiFab u;          //!< x-face normal velocity
    amrex::MultiFab v;          //!< y-face normal velocity
    amrex::MultiFab w;          //!< z-face normal velocity
    amrex::MultiFab p;          //!< perturbation pressure p' (cell centre)
    amrex::MultiFab theta;      //!< potential temperature (cell centre)
    amrex::MultiFab tke;        //!< turbulent kinetic energy k (cell centre)
    //! Physical height of each node (component 0). With an immersed ground (terrain.representation =
    //! immersed_forcing, ImmersedTerrain.H) a second component holds the ground's height in the node's
    //! column, the same at every k: Heights_K.H reads the last component as the ground, so every
    //! consumer of a height above the ground follows the immersed surface. The copies of the heights
    //! between states and the plotfile's z_phys_nd carry component 0 alone, and so do the checkpoint's: a
    //! restart fills the ground from the inputs, which the checkpoint's record pins (ImmersedTerrain.H).
    amrex::MultiFab z_phys_nd;
    int step{0};

    //! Allocate every field on the cell-centred BoxArray ba. Values are set
    //! to zero, ghost cells included. ground_component: z_phys_nd also holds the immersed ground's
    //! height (two components).
    void define (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, bool ground_component = false);

    //! The fields in checkpoint order, with their file names.
    struct Field { const char* name; amrex::MultiFab* mf; };
    amrex::Vector<Field> fields ();
};

//! Heights of a flat, uniformly spaced grid: z = prob_lo[2] + k dz, ghost
//! nodes included (extrapolated linearly), computed as fma(k, dz, z0). The
//! top node can differ from prob_hi[2] in the last bit; code that needs the
//! domain top uses Geometry::ProbHi(2). A terrain mesh (Mesh.H) replaces this.
void init_flat_z_phys_nd (amrex::MultiFab& z_phys_nd, const amrex::Geometry& geom);

//! Geometry and cell-centred BoxArray for the mesh inputs. Boxes are split in
//! x and y only, so every box spans the full z column.
amrex::Geometry make_geometry (const MeshInputs& mesh);
amrex::BoxArray make_box_array (const MeshInputs& mesh);

} // namespace gust

#endif
