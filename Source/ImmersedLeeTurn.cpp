#include "ImmersedLeeTurn.H"

#include "BuildingForcing.H"

namespace gust {

void ImmersedLeeTurn::build (const Buildings& b, const amrex::Geometry& geom, const amrex::MultiFab& detJ, amrex::Real c0,
                             amrex::Real eps)
{
    m_b = &b;
    m_detJ = &detJ;
    m_geom = geom;
    m_c0 = c0;
    m_eps = eps;
    m_active = c0 > amrex::Real(0) && b.from_surface();
    m_rates.clear();
    if (m_active) {
        if (!(eps > amrex::Real(0))) { amrex::Abort("gust: ImmersedLeeTurn needs terrain.lee_turn_ramp > 0"); }
        m_rates.define(b.beta().boxArray(), b.beta().DistributionMap(), kNumLeeTurnRates, 1);
        m_rates.setVal(0.0);
    }
}

void ImmersedLeeTurn::update (const amrex::MultiFab& u, const amrex::MultiFab& v, const std::array<amrex::MultiFab, 3>& weight)
{
    if (!m_active) { return; }
    // Ghost cells stay 0: the w rows read the cells below and above their face, which lie in the box's column
    // or, above the top face, in the ghost layer, where there is no cut cell.
    m_rates.setVal(0.0);
    const amrex::Real c0 = m_c0, eps = m_eps, dzeta = m_geom.CellSize(2);
    const amrex::Real dx = m_geom.CellSize(0), dy = m_geom.CellSize(1);
    const int khi = m_geom.Domain().bigEnd(2);
    for (amrex::MFIter mfi(m_rates, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto L = m_rates.array(mfi);
        const auto uu = u.const_array(mfi);
        const auto vv = v.const_array(mfi);
        const auto beta = m_b->beta().const_array(mfi);
        const auto wx = weight[0].const_array(mfi);
        const auto wy = weight[1].const_array(mfi);
        const auto s = m_b->surface().const_array(mfi);
        const auto J = m_detJ->const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (beta(i, j, k) <= amrex::Real(0) || beta(i, j, k) >= amrex::Real(1)) { return; }
            // The reference cell: the first fully fluid cell above in the column (the top layer holds fluid:
            // check_immersed_ground_extent).
            int kr = k + 1;
            while (kr < khi && beta(i, j, kr) > amrex::Real(0)) { ++kr; }
            const amrex::Real ur = open_face_mean(uu(i, j, kr), uu(i + 1, j, kr), wx(i, j, kr), wx(i + 1, j, kr));
            const amrex::Real vr = open_face_mean(vv(i, j, kr), vv(i, j + 1, kr), wy(i, j, kr), wy(i, j + 1, kr));
            amrex::Real lam = 0, lam_wt = 0;
            lee_turn_rate(s(i, j, k, kSurfaceSlopeX), s(i, j, k, kSurfaceSlopeY), ur, vr, J(i, j, k) * dzeta, dx, dy, c0, eps,
                          kBuildingMinSpeed, lam, lam_wt);
            L(i, j, k, kLeeTurnRate) = lam;
            L(i, j, k, kLeeTurnTarget) = lam_wt;
        });
    }
}

void ImmersedLeeTurn::add_sink (const amrex::MultiFab& w, amrex::MultiFab& Rw) const
{
    if (!m_active) { return; }
    const int kg = m_geom.Domain().smallEnd(2);
    for (amrex::MFIter mfi(Rw, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto L = m_rates.const_array(mfi);
        const auto x = w.const_array(mfi);
        const auto a = Rw.array(mfi);
        const auto fz = m_b->beta_face(2).const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k == kg || fz(i, j, k) >= amrex::Real(1)) { return; }
            const amrex::Real lam = amrex::Real(0.5) * (L(i, j, k - 1, kLeeTurnRate) + L(i, j, k, kLeeTurnRate));
            const amrex::Real tgt = amrex::Real(0.5) * (L(i, j, k - 1, kLeeTurnTarget) + L(i, j, k, kLeeTurnTarget));
            a(i, j, k) += lam * x(i, j, k) - tgt;
        });
    }
}

void ImmersedLeeTurn::add_to_implicit (amrex::MultiFab& A) const
{
    if (!m_active) { return; }
    const int kg = m_geom.Domain().smallEnd(2);
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto a = A.array(mfi);
        const auto L = m_rates.const_array(mfi);
        const auto fz = m_b->beta_face(2).const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (k == kg || fz(i, j, k) >= amrex::Real(1)) { return; }
            a(i, j, k, coef_P) += amrex::Real(0.5) * (L(i, j, k - 1, kLeeTurnRate) + L(i, j, k, kLeeTurnRate));
        });
    }
}

}   // namespace gust
