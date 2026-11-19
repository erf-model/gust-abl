#include "ImmersedSlope.H"

#include "BuildingForcing.H"
#include "Buildings_K.H"

#include <cmath>

namespace gust {

void immersed_slope_rates (const Buildings& b, const std::array<amrex::MultiFab, 3>& weight, const amrex::MultiFab& detJ,
                           const amrex::Geometry& geom, const amrex::MultiFab& u, const amrex::MultiFab& v,
                           const amrex::MultiFab& w, amrex::Real z0, amrex::Real kappa, amrex::Real umin, amrex::MultiFab& rate_mf)
{
    if (!b.from_surface()) { amrex::Abort("gust: immersed_slope_rates needs a ground built from a surface"); }
    const amrex::Real dzeta = geom.CellSize(2);
    const int khi = geom.Domain().bigEnd(2);
    for (amrex::MFIter mfi(rate_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto r = rate_mf.array(mfi);
        const auto uu = u.const_array(mfi);
        const auto vv = v.const_array(mfi);
        const auto ww = w.const_array(mfi);
        const auto beta = b.beta().const_array(mfi);
        const auto wx = weight[0].const_array(mfi);
        const auto wy = weight[1].const_array(mfi);
        const auto wz = weight[2].const_array(mfi);
        const auto wl = b.wall().const_array(mfi);
        const auto sf = b.surface().const_array(mfi);
        const auto fz = b.beta_face(2).const_array(mfi);
        const auto J = detJ.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int mask = static_cast<int>(wl(i, j, k, kWallMask));
            if (beta(i, j, k) >= amrex::Real(1) || (mask & kWallRoof) == 0) { return; }
            const amrex::Real hx = sf(i, j, k, kSurfaceSlopeX), hy = sf(i, j, k, kSurfaceSlopeY);
            const amrex::Real stretch = std::sqrt(amrex::Real(1) + hx * hx + hy * hy);
            const amrex::Real cosa = amrex::Real(1) / stretch;
            const amrex::Real nx = -hx * cosa, ny = -hy * cosa, nz = cosa;
            // The speed parallel to the ground. Without slope the w terms are exact zeros.
            auto speed = [&] (int kk) {
                const amrex::Real uc = open_face_mean(uu(i, j, kk), uu(i + 1, j, kk), wx(i, j, kk), wx(i + 1, j, kk));
                const amrex::Real vc = open_face_mean(vv(i, j, kk), vv(i, j + 1, kk), wy(i, j, kk), wy(i, j + 1, kk));
                const amrex::Real wc = open_face_mean(ww(i, j, kk), ww(i, j, kk + 1), wz(i, j, kk), wz(i, j, kk + 1));
                const amrex::Real un = uc * nx + vc * ny + wc * nz;
                const amrex::Real tx = uc - un * nx, ty = vc - un * ny, tz = wc - un * nz;
                const amrex::Real horizontal = tx * tx + ty * ty;
                return std::sqrt(horizontal + tz * tz);
            };
            const amrex::Real own = speed(k);
            const amrex::Real inv = amrex::Real(1) / (amrex::max(own, umin) * J(i, j, k) * dzeta);
            const amrex::Real inv_dn = amrex::Real(1) / (J(i, j, k) * dzeta);
            const amrex::Real d1 = wl(i, j, k, kWallD1z) * cosa, d2 = wl(i, j, k, kWallD2z) * cosa;
            const amrex::Real ref = k < khi ? speed(k + 1) : amrex::Real(0);
            const amrex::Real us2 = form7_ustar2_exact(own, ref, beta(i, j, k), d1, d2, z0, kappa);
            const amrex::Real area = wall_area_fraction(fz(i, j, k), fz(i, j, k + 1), beta(i, j, k)) * stretch;
            const amrex::Real rate = area * us2 * inv;
            amrex::Real ho = 0, hr = 0;
            form7_heat_velocities_exact(us2, beta(i, j, k), d1, d2, z0, kappa, ho, hr);
            r(i, j, k, 0) = rate;
            r(i, j, k, 1) = rate;
            r(i, j, k, 2) = (hx != amrex::Real(0) || hy != amrex::Real(0)) ? rate : amrex::Real(0);
            r(i, j, k, kBuildingUstar2) = us2;
            r(i, j, k, kBuildingHeatOwn) = area * ho * inv_dn;
            r(i, j, k, kBuildingHeatRefZ) = area * hr * inv_dn;
        });
    }
}

}   // namespace gust
