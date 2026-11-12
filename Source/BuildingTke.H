#ifndef GUST_BUILDING_TKE_H_
#define GUST_BUILDING_TKE_H_

// How buildings act on the turbulence closure (buildings.wall_tke = on; Buildings.H has the
// geometry, BuildingForcing.H the wall law; docs/buildings.rst):
//
// - a solid cell (beta = 1) holds k at the closure's floor: R_k = k - k_floor, an identity
//   row;
// - a wall cell's k row relaxes towards the wall value of the law's u*, as the ground's
//   first cell does (SurfaceLayer.H):
//     W = a (k - k_w),  k_w = max(u*^2 / Cmu0^2, k_floor),  a = max(u*, 0.01 m/s) / d_w,
//   u*^2 the largest of the cell's walls (BuildingForcing::rate, component kBuildingUstar2),
//   d_w the cell's distance below, and f the largest blocked part of the faces its walls
//   lie on. With bc.surface_wall_k = dirichlet (the default) or no surface layer the row
//   takes the place of the cell's equation by that part,
//     R_k = (1 - f) R_k + f W,
//   so a fully blocked face replaces it and a slightly blocked one leaves it nearly as it
//   was; in the ground's first cell, which has the ground's relaxation row, f W is added
//   to that row, so k lies between the two wall values. With bc.surface_wall_k = source,
//   f W is added to the cell's equation everywhere, as the ground's is;
// - the diffusive flux of k through the blocked part of a face is closed, in the residual
//   and in the implicit rows (k does not diffuse into a building);
// - the closure's geometric length l_g = kappa (d + z0) takes the nearer of the ground and
//   the nearest wall of a building, each with its own roughness:
//     d + z0 = min(d_ground + z0_ground, d_b + z0_b),  z0_b = buildings.z0
//   (limit_wall_distance), d_b the distance to the nearest wall of a building:
//     wall cell   d_b = max(min over its walls of d1, z0_b), d1 = (1 - beta) Delta_n / 2
//                 the distance of its fluid centroid (form 7);
//     elsewhere   the distance from the cell centre to the nearest wall, each wall a
//                 rectangle on its plane, placed from the cell's solid fraction and its
//                 faces' blocked fractions (BuildingTke.cpp: wall_depth, wall_rectangle; a
//                 roof over a slab beta Delta_z above the cell's floor, the side of a solid
//                 neighbour on the face, a building's edge inside the cell where it lies).
//                 It is found by passing each cell's nearest
//                 wall to its 26 neighbours until none changes (a vector distance
//                 transform): exact for a convex block, and close to it otherwise. A wall
//                 as far as the tallest column is not passed on: the ground is nearer.
//   Solid cells keep the ground's distance (the closure is not used there).
// - the solid cells' rows count as zeros in the residual norm of k (zero_solid_cells).
//
// Not here: the flat-wall dissipation length of ERF and Kynema-SGF (their wall cells keep
// the k equation, GUST's have the relaxation row, where the dissipation does not enter);
// stability functions on building walls (neutral: k_w has no buoyancy term); theta (BuildingTheta.H).
//
// Ghost cells: d_b has one layer, periodic images, else the nearest domain column's along x
// and y (as the ground's distance) and a large value below and above the domain. The rows
// read valid cells only.

#include "BuildingForcing.H"
#include "ImplicitOperator.H"
#include "SurfaceRoughness.H"

#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

#include <cmath>

namespace gust {

//! The distance from the point p to an axis-aligned rectangle of centre c and half sizes h
//! (one of them 0), with the nearest periodic image along x and y (lx, ly: the periods, 0:
//! not periodic).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real rectangle_distance (amrex::Real px, amrex::Real py, amrex::Real pz, amrex::Real cx, amrex::Real cy, amrex::Real cz,
                                amrex::Real hx, amrex::Real hy, amrex::Real hz, amrex::Real lx, amrex::Real ly) noexcept
{
    amrex::Real dx = px - cx, dy = py - cy;
    if (lx > 0) { dx -= lx * std::floor(dx / lx + amrex::Real(0.5)); }
    if (ly > 0) { dy -= ly * std::floor(dy / ly + amrex::Real(0.5)); }
    const amrex::Real ex = amrex::max(std::abs(dx) - hx, amrex::Real(0));
    const amrex::Real ey = amrex::max(std::abs(dy) - hy, amrex::Real(0));
    const amrex::Real ez = amrex::max(std::abs(pz - cz) - hz, amrex::Real(0));
    return std::sqrt(ex * ex + ey * ey + ez * ez);
}

//! Friction velocities below this do not set the relaxation rate (m/s), as the ground's.
inline constexpr amrex::Real kBuildingKUstarMin = amrex::Real(0.01);
//! "No wall found": larger than any distance in a domain.
inline constexpr amrex::Real kBuildingFar = amrex::Real(1.0e30);

//! How a wall cell's row meets the cell's own: blended in everywhere (replace), blended above
//! the ground's first cells and added to the ground's relaxation row in them (ground_row), or
//! added everywhere (source). BuildingTke.H's header has the rows.
enum class BuildingKRow { replace, ground_row, source };

//! The closure's terms of one level's buildings.
class BuildingTke
{
public:
    //! law: its rates give u*^2 (updated by the residual before the rows here are applied).
    //! z_nd: the mesh's node heights (one ghost node). tke_floor, Cmu0: the closure's
    //! (ClosureConstants::tke_floor, ClosureOptions::Cmu0).
    void build (const Buildings& b, const BuildingForcing& law, const BuildingsOptions& opt, const amrex::Geometry& geom,
                const amrex::MultiFab& z_nd, amrex::Real tke_floor, amrex::Real Cmu0);

    //! d_b at the cells (one ghost layer); kBuildingFar in solid cells and where the nearest wall
    //! is as far as the tallest column or farther (cutoff()), where the ground is nearer anyway.
    [[nodiscard]] const amrex::MultiFab& distance () const { return m_dist; }
    //! The passes the distance took, the last one, which changed nothing, included.
    [[nodiscard]] int distance_passes () const { return m_passes; }
    //! The tallest column of the mesh (m): walls at this distance or more are not looked for.
    [[nodiscard]] amrex::Real cutoff () const { return m_cutoff; }
    //! The closure's length reads d + z0 with the ground's z0 (Closure.H). Where a building's
    //! wall is the nearer one by that measure the sum must be d_b + z0_b, so
    //!   dist = min(dist, d_b + z0_b - z0)
    //! on the cells and their ghost layer, z0 the ground's (ground; its map needs a lateral ghost
    //! column). The value is the closure's alone: it is below zero where the ground's z0 exceeds
    //! d_b + z0_b. With z0_b = z0 it is min(dist, d_b), bit for bit. Returns the number of valid
    //! cells whose distance fell (all ranks).
    amrex::Long limit_wall_distance (amrex::MultiFab& dist, const Roughness& ground) const;
    //! 0 in the solid cells of a field on the cells: their rows of R_k hold k at its floor and
    //! are left out of the residual norm (Simplec.H: stage_norms).
    void zero_solid_cells (amrex::MultiFab& field) const;

    //! Gk[d], the diffusive flux of k through the d-faces: times the open part of each face.
    void close_fluxes (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open,
                       amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& Gk) const;
    //! The rows of R_k in solid and wall cells (last, after every other term of R_k).
    void apply_rows (const amrex::MultiFab& tke, amrex::MultiFab& Rk, BuildingKRow form) const;
    //! The same rows of the implicit coefficients of k, and the couplings through the
    //! closed parts of the faces.
    void add_to_implicit (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open, amrex::MultiFab& A,
                          BuildingKRow form) const;

private:
    const Buildings* m_b{nullptr};
    const BuildingForcing* m_law{nullptr};
    amrex::Geometry m_geom;
    amrex::MultiFab m_dist;
    amrex::Real m_floor{0}, m_cmu0{amrex::Real(0.5562)}, m_cutoff{0}, m_z0b{0};
    int m_passes{0};
};

} // namespace gust

#endif
