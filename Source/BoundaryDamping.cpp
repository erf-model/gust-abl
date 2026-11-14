// The damping layers of the boundary conditions: the Rayleigh layer under the
// top and the lateral sponge, in the residual and on the implicit diagonal
// (split from BoundaryConditions.cpp, audit A11).
#include "BoundaryConditions.H"
#include "BoundaryConditions_K.H"
#include "NestFaces.H"

#include "CVFluxes.H"
#include "Heights_K.H"
#include "ImplicitOperator.H"
#include "Stencil_K.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <cmath>

namespace gust {

using namespace bc_detail;

namespace {

//! The damping coefficient c(x, y, z) of both layers (1/s), device-usable.
struct DampingCoef
{
    bool rayleigh, sponge;
    amrex::Real zd, rc, ztop, sw, inv_w, ss;
    amrex::GpuArray<amrex::Real, 3> plo, phi;
    amrex::GpuArray<int, 4> lateral;   //!< xlo, xhi, ylo, yhi hold an inflow or outflow
    //! The Rayleigh part alone (1/s).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real rayleigh_part (amrex::Real z) const noexcept
    {
        const amrex::Real pi = amrex::Real(3.14159265358979323846);
        if (rayleigh && z > zd) {
            const amrex::Real s = std::sin(amrex::Real(0.5) * pi * (z - zd) / (ztop - zd));
            return rc * s * s;
        }
        return 0;
    }
    //! The sponge part alone (1/s).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real sponge_part (amrex::Real x, amrex::Real y) const noexcept
    {
        amrex::Real c = 0;
        if (sponge) {
            const amrex::Real dist[4] = {x - plo[0], phi[0] - x, y - plo[1], phi[1] - y};
            for (int n = 0; n < 4; ++n) {
                if (lateral[n] && dist[n] < sw) {
                    const amrex::Real r = 1 - dist[n] * inv_w;   // inv_w: 1/sw, finite when off
                    c += ss * r * r;
                }
            }
        }
        return c;
    }
    //! Both parts; with_rayleigh false leaves the Rayleigh part out (a
    //! field bc.rayleigh_fields does not name).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (amrex::Real x, amrex::Real y, amrex::Real z, bool with_rayleigh = true) const noexcept
    {
        return (with_rayleigh ? rayleigh_part(z) : amrex::Real(0)) + sponge_part(x, y);
    }
};

//! ztop: the mesh's top level, where the Rayleigh ramp ends.
DampingCoef damping_coef (const BoundaryInputs& in, const amrex::Geometry& geom, amrex::Real ztop)
{
    DampingCoef d{};
    d.zd = in.rayleigh_zdamp;
    d.rc = in.rayleigh_dampcoef;
    d.ztop = ztop;
    d.sw = in.sponge_width;
    d.ss = in.sponge_strength;
    d.rayleigh = d.zd > 0;
    d.sponge = d.sw > 0;
    d.inv_w = d.sponge ? 1 / d.sw : amrex::Real(0);   // no x/0 for the compiler to evaluate early
    d.plo = geom.ProbLoArray();
    d.phi = geom.ProbHiArray();
    for (int dd = 0; dd < 2; ++dd) {
        for (const bool lo : {true, false}) {
            const FaceType t = in.face[amrex::Orientation(dd, lo ? amrex::Orientation::low : amrex::Orientation::high)];
            d.lateral[2 * dd + (lo ? 0 : 1)] = (t == FaceType::inflow || t == FaceType::outflow) ? 1 : 0;
        }
    }
    return d;
}

//! Height above the ground and absolute height of the point (i, j, k) of type t.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void point_heights (int t, const NodeZ& nz, int i, int j, int k, int klo, amrex::Real& h, amrex::Real& zc) noexcept
{
    if (t == stencil::XF) { h = height_above_ground<stencil::XF>(nz, i, j, k, klo); zc = point_z<stencil::XF>(nz, i, j, k); }
    else if (t == stencil::YF) { h = height_above_ground<stencil::YF>(nz, i, j, k, klo); zc = point_z<stencil::YF>(nz, i, j, k); }
    else if (t == stencil::ZF) { h = height_above_ground<stencil::ZF>(nz, i, j, k, klo); zc = point_z<stencil::ZF>(nz, i, j, k); }
    else { h = height_above_ground<stencil::CC>(nz, i, j, k, klo); zc = point_z<stencil::CC>(nz, i, j, k); }
}

} // namespace

void BoundaryConditions::add_damping (const State& q, ResidualFields& R) const
{
    if (m_nest != nullptr) { m_nest->add_relaxation(q, R); }
    if (!damping_on()) { return; }
    const int klo = m_geom.Domain().smallEnd(2);
    const InflowProfile::View prof = m_profile.view();
    const DampingCoef coef = damping_coef(m_in, m_geom, m_mesh.zlevels().back());
    const auto plo = m_geom.ProbLoArray();
    const auto dx = m_geom.CellSizeArray();
    const amrex::MultiFab* qf[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    amrex::MultiFab* rf[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    const int ff[5] = {f_u, f_v, f_w, f_theta, f_tke};
    for (int n = 0; n < 5; ++n) {
        const int f = ff[n];
        const bool with_rayleigh = m_in.rayleigh_fields[static_cast<std::size_t>(n)];
        const bool above_ground = m_in.rayleigh_above_ground;
        for (amrex::MFIter mfi(*rf[n], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto r = rf[n]->array(mfi);
            const auto a = qf[n]->const_array(mfi);
            const amrex::Box zb = q.z_phys_nd[mfi].box();
            const NodeZ nz{q.z_phys_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
            const int tt = (n < 3) ? (1 << n) : 0;
            const bool nx = (tt & 1) != 0, ny = (tt & 2) != 0;
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                amrex::Real h, zc;
                point_heights(tt, nz, i, j, k, klo, h, zc);
                const amrex::Real x = plo[0] + (i + (nx ? 0 : amrex::Real(0.5))) * dx[0];
                const amrex::Real y = plo[1] + (j + (ny ? 0 : amrex::Real(0.5))) * dx[1];
                // The Rayleigh layer relaxes towards the free atmosphere:
                // theta towards the sounding at the point's absolute height
                // (the terrain-following profile height would impose a warm
                // anomaly over a crest), u and v towards the profile, w to 0.
                // The sponge relaxes towards the inflow profile (the height
                // above the ground, as the inflow itself).
                const amrex::Real zr = above_ground ? h : zc;   // the layer's height coordinate
                const amrex::Real cr = with_rayleigh ? coef.rayleigh_part(zr) : amrex::Real(0);
                const amrex::Real cs = coef.sponge_part(x, y);
                if (cr > 0) {
                    const amrex::Real ref = (f == f_w) ? amrex::Real(0) : profile_value(prof, f, (f == f_theta) ? zr : h);
                    r(i, j, k) += cr * (a(i, j, k) - ref);
                }
                if (cs > 0) {
                    const amrex::Real ref = (f == f_w) ? amrex::Real(0) : profile_value(prof, f, h);
                    r(i, j, k) += cs * (a(i, j, k) - ref);
                }
            });
        }
    }
}

BoundaryConditions::CanopyZoneCounts BoundaryConditions::canopy_in_boundary_zones (const amrex::MultiFab& f,
                                                                                   const amrex::MultiFab& z_nd) const
{
    const amrex::Box& dom = m_geom.Domain();
    const int klo = dom.smallEnd(2);
    // The first cell column inside each inflow face.
    const bool in_xlo = m_in.face[amrex::Orientation(0, amrex::Orientation::low)] == FaceType::inflow;
    const bool in_xhi = m_in.face[amrex::Orientation(0, amrex::Orientation::high)] == FaceType::inflow;
    const bool in_ylo = m_in.face[amrex::Orientation(1, amrex::Orientation::low)] == FaceType::inflow;
    const bool in_yhi = m_in.face[amrex::Orientation(1, amrex::Orientation::high)] == FaceType::inflow;
    const int x0 = dom.smallEnd(0), x1 = dom.bigEnd(0), y0 = dom.smallEnd(1), y1 = dom.bigEnd(1);
    const DampingCoef coef = damping_coef(m_in, m_geom, m_mesh.zlevels().back());
    const bool above_ground = m_in.rayleigh_above_ground;
    // The Rayleigh layer counts when it relaxes a field other than w (u, v, theta, k) towards the profile.
    const auto& rf = m_in.rayleigh_fields;
    const bool rayleigh_to_profile = rf[0] || rf[1] || rf[3] || rf[4];
    const auto plo = m_geom.ProbLoArray();
    const auto dx = m_geom.CellSizeArray();
    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> rop;
    amrex::ReduceData<amrex::Long, amrex::Long> data(rop);
    for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const auto fa = f.const_array(mfi);
        const amrex::Box zb = z_nd[mfi].box();
        const NodeZ nz{z_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
        rop.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long, amrex::Long>
        {
            if (!(fa(i, j, k) > 0)) { return {0, 0}; }
            const bool inflow = (in_xlo && i == x0) || (in_xhi && i == x1) || (in_ylo && j == y0) || (in_yhi && j == y1);
            amrex::Real h, zc;
            point_heights(stencil::CC, nz, i, j, k, klo, h, zc);
            const amrex::Real x = plo[0] + (i + amrex::Real(0.5)) * dx[0];
            const amrex::Real y = plo[1] + (j + amrex::Real(0.5)) * dx[1];
            const bool damped = (rayleigh_to_profile && coef.rayleigh_part(above_ground ? h : zc) > 0) || coef.sponge_part(x, y) > 0;
            return {inflow ? amrex::Long(1) : amrex::Long(0), damped ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    auto v = data.value(rop);
    amrex::Long n[2] = {amrex::get<0>(v), amrex::get<1>(v)};
    amrex::ParallelDescriptor::ReduceLongSum(n, 2);
    return {n[0], n[1]};
}

void BoundaryConditions::add_damping_to_implicit (Equation eq, Loc loc, const amrex::MultiFab& z_nd, amrex::MultiFab& A) const
{
    if (m_nest != nullptr) { m_nest->add_relaxation_to_implicit(eq, loc, A); }
    if (!damping_on()) { return; }
    const int klo = m_geom.Domain().smallEnd(2);
    const DampingCoef coef = damping_coef(m_in, m_geom, m_mesh.zlevels().back());
    const int n = (eq == Equation::u) ? 0 : (eq == Equation::v) ? 1 : (eq == Equation::w) ? 2 : (eq == Equation::theta) ? 3 : 4;
    const bool with_rayleigh = m_in.rayleigh_fields[static_cast<std::size_t>(n)];
    const bool above_ground = m_in.rayleigh_above_ground;
    const auto plo = m_geom.ProbLoArray();
    const auto dx = m_geom.CellSizeArray();
    const int tt = static_cast<int>(loc);
    const bool nx = (tt & 1) != 0, ny = (tt & 2) != 0;
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto a = A.array(mfi);
        const amrex::Box zb = z_nd[mfi].box();
        const NodeZ nz{z_nd.const_array(mfi), amrex::lbound(zb), amrex::ubound(zb)};
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real h, zc;
            point_heights(tt, nz, i, j, k, klo, h, zc);
            const amrex::Real x = plo[0] + (i + (nx ? 0 : amrex::Real(0.5))) * dx[0];
            const amrex::Real y = plo[1] + (j + (ny ? 0 : amrex::Real(0.5))) * dx[1];
            a(i, j, k, coef_P) += (with_rayleigh ? coef.rayleigh_part(above_ground ? h : zc) : amrex::Real(0)) + coef.sponge_part(x, y);
        });
    }
}

} // namespace gust
