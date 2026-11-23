#ifndef GUST_HYBRID_TERRAIN_H_
#define GUST_HYBRID_TERRAIN_H_

// The hybrid mesh (terrain.representation = hybrid; docs/immersed_terrain.rst, Hybrid mesh). The mesh
// follows a smoothed lower envelope h_s <= h of the terrain's height h, and the ground, at
// terrain.flat_height plus h, is immersed in it as on the flat mesh of immersed_forcing
// (ImmersedTerrain.H). The gap between the mesh's lower face and the ground is
// terrain.flat_height + h - h_s >= terrain.flat_height >= 0. A flat envelope gives immersed_forcing's
// mesh; terrain.hybrid_half_width = 0 gives h_s = h, a body-fitted mesh with the ground
// terrain.flat_height above its lower face. With terrain.flat_height = 0 the lower face is the ground
// wherever h_s = h: those columns take the surface layer at their lower face, and the others, whose
// lower face lies inside the terrain, the immersed ground's wall law (hybrid_immersed_columns).
//
// The envelope: the running minimum of h over a square of half-width p nodes in x and q in y, then the
// running mean of that minimum over the same square. Each node lies in the square of every node of its
// own square, so the mean of the minima is at most h there: h_s <= h everywhere. At a face that is not
// periodic the square is cut at the face (both passes); across a periodic face it wraps.

#include "Inputs.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Vector.H>

#include <string>
#include <vector>

namespace gust {

//! The envelope of a node array h (nx * ny values, x fastest): the running minimum over the square of
//! half-width px by py nodes, then the running mean of that over the same square, each at most h at
//! its node (min(mean, h): the mean of values at most h stays at most h after round-off). perx/pery:
//! the direction is periodic, node n standing for node 0 (h holds n distinct nodes), and the square
//! wraps; else it is cut at the array's ends. px = py = 0 returns h.
std::vector<double> lower_envelope (const std::vector<double>& h, int nx, int ny, int px, int py, bool perx, bool pery);

//! The defaults of terrain.representation = hybrid, set once the inputs have passed validate(): without
//! terrain.flat_height, the ground lies one first-layer thickness (z_1 of make_zlevels) above the
//! envelope; without terrain.hybrid_half_width, the half-width is 2 cells (twice the larger cell size of
//! the lateral directions with more than one cell; 0 when neither has). Another representation, or a
//! key given in the inputs, is left as it is. Returns one log line per default it set, and a warning line
//! when a given width is below one cell (no envelope) or the terrain is flat (the immersed_forcing mesh).
std::vector<std::string> apply_hybrid_defaults (Inputs& in);

//! The envelope's half-widths in nodes for terrain.hybrid_half_width = width (m): floor(width / dx + tol),
//! tol = max(1e-9, 64 eps(Real)), in each lateral direction.
amrex::IntVect hybrid_envelope_nodes (amrex::Real width, const amrex::Geometry& geom);

//! Replaces the terrain's height in z_nd's component 0 at the domain's lowest nodes (valid and ghost
//! columns; a ghost column takes its periodic image's or the nearest domain column's value, as
//! fill_terrain_surface's) with its envelope for terrain.hybrid_half_width = width. The surface is
//! gathered onto every rank (one value per domain node), so the envelope does not depend on the boxes
//! or the rank count. A collective call.
void make_hybrid_envelope (amrex::MultiFab& z_nd, const amrex::Geometry& geom, amrex::Real width);

//! The refusal of the hybrid mesh that needs its node heights: at every node the ground (z_nd's last
//! component) must lie at or below the bottom of the column's top layer. Empty when it does. A
//! collective call.
std::string check_hybrid_ground_extent (const amrex::MultiFab& z_nd, const amrex::Geometry& geom);

//! The gap h + flat_height - h_s between the mesh's lower face and the ground (z_nd's last component
//! minus component 0 at the lowest nodes): its smallest and largest values over the valid nodes, and the
//! height above the lower face to which zlevels keeps the first layer's thickness (within 1e-6 of it,
//! from the bottom): the layers through the gap are uniform where the largest gap lies below it.
struct HybridGap
{
    double min{0}, max{0};
    double uniform_top{0};
};
HybridGap hybrid_gap (const amrex::MultiFab& z_nd, const amrex::Vector<amrex::Real>& zlevels);

//! The columns whose lower face lies inside the terrain: 1 where the first cell's solid fraction
//! (beta, one ghost layer whose lateral ghosts are filled) is above 0, else 0, at the domain's lowest
//! layer of slab_ba (the surface state's layout, ground_slab) and its ghost columns (ng lateral ghosts,
//! at most beta's). The surface layer acts at the lower face where it is 0
//! (SurfaceLayerOptions::immersed_columns). Returns the number of valid columns of 1 (a collective call).
amrex::Long hybrid_immersed_columns (const amrex::MultiFab& beta, const amrex::BoxArray& slab_ba,
                                     const amrex::DistributionMapping& dm, int ng, amrex::MultiFab& mask);

//! p with each solid cell (beta = 1) given the pressure of the first cell above it in its column that is
//! not solid, into ps (defined on p's layout when it is not): the valid cells, the lateral ghost columns
//! (one) and the ghost layer below the domain. The pressure gradient's slope terms (pressure_gradient:
//! -z_x dp/dz on the x faces) then take no vertical difference across the terrain's solid boundary,
//! where the solid cells hold 0. The solid cells lie below the fluid in every column, and the boxes span
//! the column (aborts otherwise). beta and p need one ghost cell, filled laterally; aborts otherwise. Only one
//! lateral ghost ring is filled: a further ring keeps p as it is.
void hybrid_solid_pressure (const amrex::MultiFab& beta, const amrex::MultiFab& p, amrex::MultiFab& ps);

}   // namespace gust

#endif
