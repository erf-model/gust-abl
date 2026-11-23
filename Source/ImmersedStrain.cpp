#include "ImmersedStrain.H"

#include "Stencil_K.H"

#include <cmath>
#include <string>

namespace gust {

namespace {

//! apply_cut_cell_strain with the accessor Z0.
template <typename Z0>
void cut_cell_strain (const Buildings& b, const Roughness& rough, const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2)
{
    using namespace stencil;
    const auto dxinv = geom.InvCellSizeArray();
    const int klo = geom.Domain().smallEnd(2), khi = geom.Domain().bigEnd(2);
    for (amrex::MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0));
        const auto s2 = S2.array(mfi);
        const auto u = q.u.const_array(mfi);
        const auto v = q.v.const_array(mfi);
        const auto w = q.w.const_array(mfi);
        const auto z = q.z_phys_nd.const_array(mfi);
        const auto beta = b.beta().const_array(mfi);
        const auto wl = b.wall().const_array(mfi);
        const Z0 zf = Z0::of(rough, mfi);
        // The cut cell itself: the wall law's gradient at its fluid centroid, u* / (kappa (d1 + z0)) along its
        // wind, as the body-fitted ground's first cell has it (BoundaryConditions.cpp: ground_gradients). The
        // closure's difference reads the solid cell below (u = 0). kappa cancels (u* = kappa U / ln).
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k < klo || k > khi) { return; }
            if (!(beta(i, j, k) < amrex::Real(1) && (static_cast<int>(wl(i, j, k, kWallMask)) & kWallRoof) != 0)) { return; }
            auto uc = [&] (int kk) { return amrex::Real(0.5) * (u(i, j, kk) + u(i + 1, j, kk)); };
            auto vc = [&] (int kk) { return amrex::Real(0.5) * (v(i, j, kk) + v(i, j + 1, kk)); };
            const amrex::Real d1 = wl(i, j, k, kWallD1z), d2 = wl(i, j, k, kWallD2z), z0 = zf(i, j);
            const amrex::Real uo = uc(k), vo = vc(k);
            const amrex::Real own = std::sqrt(uo * uo + vo * vo);
            const amrex::Real ref = k < khi ? std::sqrt(uc(k + 1) * uc(k + 1) + vc(k + 1) * vc(k + 1)) : amrex::Real(0);
            const amrex::Real us_k = std::sqrt(form7_ustar2_exact(own, ref, beta(i, j, k), d1, d2, z0, amrex::Real(1)));
            // No division by a zero speed: a cut cell at rest has no wind direction, and its gradient is 0
            // (us_k is not 0 there when the cell above has wind).
            const amrex::Real g = own > amrex::Real(0) ? us_k / ((d1 + z0) * own) : amrex::Real(0);
            const amrex::Real gu = g * uo, gv = g * vo;
            const Metric m = metric_at<CC>(z, i, j, k, dxinv);
            const auto Gu = grad_at<XF, CC>(u, i, j, k, m, dxinv);
            const auto Gv = grad_at<YF, CC>(v, i, j, k, m, dxinv);
            const auto Gw = grad_at<ZF, CC>(w, i, j, k, m, dxinv);
            const amrex::Real e_old = Gu[2] + Gw[0], f_old = Gv[2] + Gw[1];
            const amrex::Real e_new = gu + Gw[0], f_new = gv + Gw[1];
            s2(i, j, k) = amrex::max(amrex::Real(0), s2(i, j, k) + (e_new * e_new - e_old * e_old) + (f_new * f_new - f_old * f_old));
        });
    }
}

//! S2 in the first full cell above each cut cell (ImmersedStrain.H), with the accessor Z0.
template <typename Z0>
void above_cut_cell_strain (const Buildings& b, const Roughness& rough, const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2)
{
    using namespace stencil;
    const auto dxinv = geom.InvCellSizeArray();
    const int klo = geom.Domain().smallEnd(2), khi = geom.Domain().bigEnd(2);
    for (amrex::MFIter mfi(S2, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0));
        const auto s2 = S2.array(mfi);
        const auto u = q.u.const_array(mfi);
        const auto v = q.v.const_array(mfi);
        const auto w = q.w.const_array(mfi);
        const auto z = q.z_phys_nd.const_array(mfi);
        const auto beta = b.beta().const_array(mfi);
        const auto wl = b.wall().const_array(mfi);
        const Z0 zf = Z0::of(rough, mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int c = k - 1;   // the cut cell below
            if (c < klo || k >= khi) { return; }
            auto roof = [&] (int kk) {
                return beta(i, j, kk) < amrex::Real(1) && (static_cast<int>(wl(i, j, kk, kWallMask)) & kWallRoof) != 0;
            };
            if (roof(k) || !roof(c)) { return; }
            const amrex::Real z0 = zf(i, j);
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

}   // namespace

void check_ground_roughness (const char* op, const Buildings& b, const Roughness& rough)
{
    if (rough.map != nullptr) {
        check_roughness_map(op, *rough.map, b.beta());
    } else if (!(rough.z0 > 0)) {
        amrex::Abort(std::string("gust: ") + op + " needs z0 > 0");
    }
}

void apply_cut_cell_strain (const Buildings& b, const Roughness& rough, const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2)
{
    if (rough.map != nullptr) {
        cut_cell_strain<MapZ0>(b, rough, q, geom, S2);
    } else {
        cut_cell_strain<UniformZ0>(b, rough, q, geom, S2);
    }
}

ImmersedStrain::ImmersedStrain (const BuildingForcing& law, const Roughness& rough) : m_b(law.geometry()), m_rough(rough)
{
    if (m_b == nullptr || !law.immersed_ground()) {
        amrex::Abort("gust: ImmersedStrain needs the wall law of an immersed ground (terrain.representation = immersed_forcing or hybrid)");
    }
    check_ground_roughness("ImmersedStrain", *m_b, rough);
}

void ImmersedStrain::apply (const State& q, const amrex::Geometry& geom, amrex::MultiFab& S2) const
{
    if (m_rough.map != nullptr) {
        above_cut_cell_strain<MapZ0>(*m_b, m_rough, q, geom, S2);
    } else {
        above_cut_cell_strain<UniformZ0>(*m_b, m_rough, q, geom, S2);
    }
    apply_cut_cell_strain(*m_b, m_rough, q, geom, S2);
}

}   // namespace gust
