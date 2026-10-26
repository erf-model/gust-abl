#ifndef GUST_BOUNDARY_CONDITIONS_K_H_
#define GUST_BOUNDARY_CONDITIONS_K_H_

// Internal to the boundary-condition files (BoundaryConditions.cpp and
// BoundaryDamping.cpp): the field numbering of the fills and the inflow
// profile's value per field.

#include "InflowProfile.H"

#include <AMReX_Extension.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>

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

} // namespace gust::bc_detail

#endif
