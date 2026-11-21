#include "LogShear.H"

#include "ImmersedStrain.H"

#include "Stencil_K.H"

#include <cmath>

namespace gust {

namespace {

//! The (u, z) and (v, z) terms of S2 in one cell with the log-law gradients (LogShear.H).
//! d holds d + z0.
template <typename U, typename V, typename W, typename Z, typename D, typename S>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void log_shear_cell (int i, int j, int k, const U& u, const V& v, const W& w, const Z& z, const D& d, const S& s2,
                     const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxinv) noexcept
{
    using namespace stencil;
    const amrex::Real dm = d(i, j, k - 1), dc = d(i, j, k), dp = d(i, j, k + 1);
    // Not where the distance does not grow upward (under an overhang, next to another wall).
    if (!(dm > amrex::Real(0) && dc > dm && dp > dc)) { return; }
    const amrex::Real lg = dc * std::log(dp / dm);
    const amrex::Real gu = amrex::Real(0.5) * (u(i, j, k + 1) + u(i + 1, j, k + 1) - u(i, j, k - 1) - u(i + 1, j, k - 1)) / lg;
    const amrex::Real gv = amrex::Real(0.5) * (v(i, j, k + 1) + v(i, j + 1, k + 1) - v(i, j, k - 1) - v(i, j + 1, k - 1)) / lg;
    // The closure's (u, z) and (v, z) terms of S2 (Closure.cpp: stratification_and_strain).
    const Metric m = metric_at<CC>(z, i, j, k, dxinv);
    const auto Gu = grad_at<XF, CC>(u, i, j, k, m, dxinv);
    const auto Gv = grad_at<YF, CC>(v, i, j, k, m, dxinv);
    const auto Gw = grad_at<ZF, CC>(w, i, j, k, m, dxinv);
    const amrex::Real e_old = Gu[2] + Gw[0], f_old = Gv[2] + Gw[1];
    const amrex::Real e_new = gu + Gw[0], f_new = gv + Gw[1];
    // Not below zero: the closure's S2 holds e_old^2 + f_old^2 to round-off only.
    s2(i, j, k) = amrex::max(amrex::Real(0), s2(i, j, k) + (e_new * e_new - e_old * e_old) + (f_new * f_new - f_old * f_old));
}

//! len = dist + z0(i, j) on every cell of len, ghost cells included.
template <typename Z0>
void add_roughness (amrex::MultiFab& len, const amrex::MultiFab& dist, const Roughness& rough)
{
    for (amrex::MFIter mfi(len, false); mfi.isValid(); ++mfi) {
        const auto o = len.array(mfi);
        const auto d = dist.const_array(mfi);
        const Z0 zf = Z0::of(rough, mfi);
        amrex::ParallelFor(mfi.fabbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            o(i, j, k) = d(i, j, k) + zf(i, j);
        });
    }
}

}   // namespace

LogShear::LogShear (const amrex::MultiFab& dist, const Roughness& rough, const Buildings* ground) : m_b(ground)
{
    if (dist.nGrowVect()[0] < 1 || dist.nGrowVect()[1] < 1) { amrex::Abort("gust: LogShear: the wall distance needs one lateral ghost cell"); }
    // One lateral ghost cell (apply's cells), none in z: the stencil stays inside the domain.
    m_len.define(dist.boxArray(), dist.DistributionMap(), 1, amrex::IntVect(1, 1, 0));
    if (m_b != nullptr) {
        if (rough.map != nullptr) { amrex::Abort("gust: turb.shear_stencil = log_law: an immersed ground needs a uniform bc.surface_z0"); }
        m_z0 = rough.z0;
    }
    if (rough.map == nullptr) {
        if (!(rough.z0 > 0)) { amrex::Abort("gust: turb.shear_stencil = log_law needs a roughness length z0 > 0 (bc.surface_z0)"); }
        add_roughness<UniformZ0>(m_len, dist, rough);
    } else {
        check_roughness_map("turb.shear_stencil = log_law", *rough.map, dist);
        add_roughness<MapZ0>(m_len, dist, rough);
    }
}

void LogShear::apply (const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2) const
{
    const auto dxinv = geom.InvCellSizeArray();
    const int klo = geom.Domain().smallEnd(2), khi = geom.Domain().bigEnd(2);
    for (amrex::MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        amrex::Box bx = amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0));
        bx.setSmall(2, amrex::max(bx.smallEnd(2), klo + 1));
        bx.setBig(2, amrex::min(bx.bigEnd(2), khi - 1));
        if (!bx.ok()) { continue; }
        const auto s2 = S2.array(mfi);
        const auto u = q.u.const_array(mfi);
        const auto v = q.v.const_array(mfi);
        const auto w = q.w.const_array(mfi);
        const auto z = q.z_phys_nd.const_array(mfi);
        const auto d = m_len.const_array(mfi);
        if (m_b == nullptr) {
            // A body-fitted ground: every cell above the first.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                log_shear_cell(i, j, k, u, v, w, z, d, s2, dxinv);
            });
        } else {
            // An immersed ground: the fluid cells above the wall cell (the cut cell, or the first
            // fluid cell over a ground on a layer's face: the cell the roof bit marks). A cell with only
            // side bits (terrain.immersed_side_walls = steps) is not a cut cell and takes the log law.
            const auto beta = m_b->beta().const_array(mfi);
            const auto wl = m_b->wall().const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (beta(i, j, k) >= amrex::Real(1)) { return; }
                if ((static_cast<int>(wl(i, j, k, kWallMask)) & (kWallRoof | kWallUnder)) != 0) { return; }
                log_shear_cell(i, j, k, u, v, w, z, d, s2, dxinv);
            });
        }
    }
    // The cut cells of an immersed ground: the wall law's gradient (ImmersedStrain.H).
    if (m_b != nullptr) { apply_cut_cell_strain(*m_b, m_z0, q, geom, S2); }
}

}   // namespace gust
