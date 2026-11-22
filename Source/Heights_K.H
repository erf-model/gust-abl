#ifndef GUST_HEIGHTS_K_H_
#define GUST_HEIGHTS_K_H_

// Heights of the points of the terrain-following mesh from its nodes
// (shared by the boundary fills and the initializer): the
// node heights with the lateral indices clamped to the stored ghost node and
// the vertical ones extrapolated, a point's height above its column's ground
// and its vertical height, for every location type of Stencil_K.H.
//
// The ground of a column is the last component of the node heights at the
// ground layer: the mesh's lowest nodes on a body-fitted mesh (one component),
// the immersed ground's height where the heights carry it as a second
// component (State.H; terrain.representation = immersed_forcing or hybrid). A point
// below an immersed ground has a negative height.

#include "Stencil_K.H"

#include <AMReX_Array4.H>
#include <AMReX_Dim3.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

namespace gust {

//! Node heights with the lateral indices clamped to the stored nodes and the
//! vertical ones extrapolated linearly beyond them (the second ghost layer).
struct NodeZ
{
    amrex::Array4<amrex::Real const> z;
    amrex::Dim3 lo, hi;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int k) const noexcept
    {
        i = amrex::max(lo.x, amrex::min(hi.x, i));
        j = amrex::max(lo.y, amrex::min(hi.y, j));
        if (k < lo.z) { return z(i, j, lo.z) + (lo.z - k) * (z(i, j, lo.z) - z(i, j, lo.z + 1)); }
        if (k > hi.z) { return z(i, j, hi.z) + (k - hi.z) * (z(i, j, hi.z) - z(i, j, hi.z - 1)); }
        return z(i, j, k);
    }
    //! The ground's height in the column of node (i, j), lateral indices clamped; klo: the ground layer.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real ground (int i, int j, int klo) const noexcept
    {
        i = amrex::max(lo.x, amrex::min(hi.x, i));
        j = amrex::max(lo.y, amrex::min(hi.y, j));
        return z(i, j, klo, z.nComp() - 1);
    }
};

//! Height of the point (i, j, k) of type T: its node-mean height minus the
//! node-mean ground under its footprint (NodeZ::ground at the ground layer klo).
//! A point centred in a layer that an immersed ground cuts takes the height of the centroid of the
//! layer's part above the ground, half the open height (the inflow profile, the first guess and the
//! damping heights then read the fluid's height, as the wall law does).
template <int T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real height_above_ground (const NodeZ& z, int i, int j, int k, int klo) noexcept
{
    using namespace stencil;
    constexpr int nx = nodal(T, 0) ? 1 : 2, ny = nodal(T, 1) ? 1 : 2, nz = nodal(T, 2) ? 1 : 2;
    amrex::Real zp = 0, zg = 0, zb = 0, zt = 0;
    for (int a = 0; a < nx; ++a) {
        for (int b = 0; b < ny; ++b) {
            zg += z.ground(i + a, j + b, klo);
            for (int c = 0; c < nz; ++c) { zp += z(i + a, j + b, k + c); }
            zb += z(i + a, j + b, k);
            zt += z(i + a, j + b, k + nz - 1);
        }
    }
    zg /= (nx * ny);
    // A point whose layer is cut by an immersed ground: the centroid of the layer's part above the ground.
    // (On a body-fitted mesh no layer is cut: the ground is the lowest nodes.)
    if (nz == 2 && zg * (nx * ny) > zb && zg * (nx * ny) < zt) { return amrex::Real(0.5) * (zt / (nx * ny) - zg); }
    return zp / (nx * ny * nz) - zg;
}

//! Vertical height of the point (its node mean), for the gradient fills.
template <int T>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real point_z (const NodeZ& z, int i, int j, int k) noexcept
{
    using namespace stencil;
    constexpr int nx = nodal(T, 0) ? 1 : 2, ny = nodal(T, 1) ? 1 : 2, nz = nodal(T, 2) ? 1 : 2;
    amrex::Real zp = 0;
    for (int a = 0; a < nx; ++a) {
        for (int b = 0; b < ny; ++b) {
            for (int c = 0; c < nz; ++c) { zp += z(i + a, j + b, k + c); }
        }
    }
    return zp / (nx * ny * nz);
}

} // namespace gust

#endif
