#ifndef GUST_BOUNDARY_CONDITIONS_K_H_
#define GUST_BOUNDARY_CONDITIONS_K_H_

// Internal to the boundary-condition files (BoundaryConditions.cpp and
// BoundaryDamping.cpp) and to NestFaces.cpp: the field numbering of the fills,
// the inflow profile's value per field and the ghost slabs of a face.

#include "InflowProfile.H"

#include <AMReX_Box.H>
#include <AMReX_Extension.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

#include <algorithm>

namespace gust::bc_detail {

//! The fields the fills act on, in the order of State::fields.
enum Field : int { f_u = 0, f_v, f_w, f_p, f_theta, f_tke, n_fields };

//! The inflow profile's value for field f at height h (0 when homogeneous).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real profile_value (const InflowProfile::View& p, int f, amrex::Real h) noexcept
{
    amrex::Real u, v, t, k, nut;
    p.eval(h, u, v, t, k, nut);
    switch (f) {
    case f_u: return u;
    case f_v: return v;
    case f_theta: return t;
    case f_tke: return k;
    default: return 0;
    }
}

//! The ghost slab of layer L (1 or 2) outside face (d, lo) of a box, grown
//! over the other directions as fill_bc_ghosts: fully over the directions
//! already filled (e < d), and over later ones only where values exist
//! (inside the domain, or across periodic faces).
inline amrex::Box ghost_slab (const amrex::Box& vb, const amrex::Box& dom, const amrex::Geometry& geom, int d, bool lo,
                              int layer, int ng)
{
    const bool nodal_d = vb.type(d) == amrex::IndexType::NODE;
    const int edge = lo ? dom.smallEnd(d) : (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
    amrex::Box gb = vb;
    for (int e = 0; e < AMREX_SPACEDIM; ++e) {
        if (e == d) { continue; }
        gb.grow(e, ng);
        if (e > d && !geom.isPeriodic(e)) {
            const bool nodal_e = vb.type(e) == amrex::IndexType::NODE;
            gb.setSmall(e, std::max(gb.smallEnd(e), dom.smallEnd(e)));
            gb.setBig(e, std::min(gb.bigEnd(e), nodal_e ? dom.bigEnd(e) + 1 : dom.bigEnd(e)));
        }
    }
    const int g = lo ? edge - layer : edge + layer;
    gb.setSmall(d, g);
    gb.setBig(d, g);
    return gb;
}

//! Whether box vb touches the domain face (d, lo).
inline bool touches (const amrex::Box& vb, const amrex::Box& dom, int d, bool lo)
{
    const bool nodal_d = vb.type(d) == amrex::IndexType::NODE;
    return lo ? vb.smallEnd(d) == dom.smallEnd(d) : vb.bigEnd(d) == (nodal_d ? dom.bigEnd(d) + 1 : dom.bigEnd(d));
}

} // namespace gust::bc_detail

#endif
