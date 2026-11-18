#include "ImmersedStrain.H"

#include "Stencil_K.H"

#include <cmath>

namespace gust {

ImmersedStrain::ImmersedStrain (const BuildingForcing& law, amrex::Real z0) : m_b(law.geometry()), m_z0(z0)
{
    if (m_b == nullptr || !law.immersed_ground()) {
        amrex::Abort("gust: ImmersedStrain needs the wall law of an immersed ground (terrain.representation = immersed_forcing)");
    }
    if (!(z0 > 0)) { amrex::Abort("gust: ImmersedStrain needs z0 > 0"); }
}

void ImmersedStrain::apply (const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2) const
{
    using namespace stencil;
    const auto dxinv = geom.InvCellSizeArray();
    const amrex::Real z0 = m_z0;
    const int klo = geom.Domain().smallEnd(2), khi = geom.Domain().bigEnd(2);
    for (amrex::MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0));
        const auto s2 = S2.array(mfi);
        const auto u = q.u.const_array(mfi);
        const auto v = q.v.const_array(mfi);
        const auto w = q.w.const_array(mfi);
        const auto z = q.z_phys_nd.const_array(mfi);
        const auto beta = m_b->beta().const_array(mfi);
        const auto wl = m_b->wall().const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int c = k - 1;   // the cut cell below
            if (c < klo || k >= khi) { return; }
            auto roof = [&] (int kk) {
                return beta(i, j, kk) < amrex::Real(1) && (static_cast<int>(wl(i, j, kk, kWallMask)) & kWallRoof) != 0;
            };
            if (roof(k) || !roof(c)) { return; }
            auto zc = [&] (int kk) {
                return amrex::Real(0.125) * (z(i, j, kk) + z(i + 1, j, kk) + z(i, j + 1, kk) + z(i + 1, j + 1, kk) + z(i, j, kk + 1)
                                             + z(i + 1, j, kk + 1) + z(i, j + 1, kk + 1) + z(i + 1, j + 1, kk + 1));
            };
            auto height = [&] (int kk) {
                return amrex::Real(0.25) * (z(i, j, kk + 1) + z(i + 1, j, kk + 1) + z(i, j + 1, kk + 1) + z(i + 1, j + 1, kk + 1)
                                            - z(i, j, kk) - z(i + 1, j, kk) - z(i, j + 1, kk) - z(i + 1, j + 1, kk));
            };
            auto uc = [&] (int kk) { return amrex::Real(0.5) * (u(i, j, kk) + u(i + 1, j, kk)); };
            auto vc = [&] (int kk) { return amrex::Real(0.5) * (v(i, j, kk) + v(i, j + 1, kk)); };
            const amrex::Real h = zc(k + 1) - zc(k);
            const amrex::Real f2 = amrex::Real(0.5) * height(c), f3 = height(c) + amrex::Real(0.5) * height(k), f4 = f3 + h;
            const amrex::Real ratio = (f3 + z0) * std::log((f4 + z0) / (f2 + z0)) / (amrex::Real(2) * height(k));
            const amrex::Real d1 = wl(i, j, c, kWallD1z), d2 = wl(i, j, c, kWallD2z);
            const amrex::Real lg = (d2 + z0) * std::log((d2 + h + z0) / (d1 + z0));
            const amrex::Real gu = ratio * (uc(k + 1) - uc(c)) / lg;
            const amrex::Real gv = ratio * (vc(k + 1) - vc(c)) / lg;
            // The closure's (u, z) and (v, z) terms of S2 (Closure.cpp: stratification_and_strain).
            const Metric m = metric_at<CC>(z, i, j, k, dxinv);
            const auto Gu = grad_at<XF, CC>(u, i, j, k, m, dxinv);
            const auto Gv = grad_at<YF, CC>(v, i, j, k, m, dxinv);
            const auto Gw = grad_at<ZF, CC>(w, i, j, k, m, dxinv);
            const amrex::Real e_old = Gu[2] + Gw[0], f_old = Gv[2] + Gw[1];
            const amrex::Real e_new = gu + Gw[0], f_new = gv + Gw[1];
            // Not below zero: the closure's S2 holds e_old^2 + f_old^2 to round-off only.
            s2(i, j, k) = amrex::max(amrex::Real(0), s2(i, j, k) + (e_new * e_new - e_old * e_old) + (f_new * f_new - f_old * f_old));
        });
    }
}

}   // namespace gust
