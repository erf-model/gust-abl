#ifndef GUST_ERF_SHIM_INDEX_DEFINES_H_
#define GUST_ERF_SHIM_INDEX_DEFINES_H_

// Shim for ERF's ERF_IndexDefines.H, which the files in extern_erf/ include.
// They use only the boundary-condition enum below and the two state
// component indices that ERF_CanopyTurbulence.H's AddCanopyTKESources names
// (not called by GUST); both are ERF's (Source/ERF_IndexDefines.H at the
// commit in external/VERSIONS), the enumerators in ERF's order. ERF's
// licence: extern_erf/LICENSE.md.

enum struct ERF_BC {
    symmetry,
    inflow,
    outflow,
    inflow_outflow,
    ho_outflow,
    open,
    no_slip_wall,
    slip_wall,
    periodic,
    surface_layer,
    undefined
};

#ifndef Rho_comp
#define Rho_comp       0
#endif
#ifndef RhoKE_comp
#define RhoKE_comp     (Rho_comp+2)
#endif

#endif
