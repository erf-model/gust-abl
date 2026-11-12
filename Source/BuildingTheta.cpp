#include "BuildingTheta.H"

#include "Reductions.H"

#include <AMReX_ParallelDescriptor.H>

namespace gust {

namespace {

//! What a cell's row needs: its kind, the walls' coefficients and theta_b.
struct ThetaRow
{
    amrex::Array4<const amrex::Real> beta, wl, rate, tb;
    //! 0: the cell has no row here; 1: solid; 2: a wall cell.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int kind (int i, int j, int k) const noexcept
    {
        if (beta(i, j, k) >= amrex::Real(1)) { return 1; }
        return static_cast<int>(wl(i, j, k, kWallMask)) != 0 ? 2 : 0;
    }
    //! The walls' exchange of the cell, sum of f q / Delta_n (K/s). A partly solid cell has one
    //! wall along a direction at most, and its reference cell is the next one out from it; in a
    //! fluid cell the reference coefficients are 0.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real exchange (const amrex::Array4<const amrex::Real>& t, int i, int j, int k) const noexcept
    {
        const int mask = static_cast<int>(wl(i, j, k, kWallMask));
        // Each temperature against theta_b of its own cell: with the reference profile as theta_b a
        // state at that profile exchanges nothing, also where the profile changes with height.
        const int di = (mask & kWallXlo) != 0 ? 1 : -1, dj = (mask & kWallYlo) != 0 ? 1 : -1, dk = (mask & kWallRoof) != 0 ? 1 : -1;
        amrex::Real s = rate(i, j, k, kBuildingHeatOwn) * (t(i, j, k) - tb(i, j, k));
        // A coefficient of 0 leaves its neighbour unread.
        if (rate(i, j, k, kBuildingHeatRefX) != amrex::Real(0)) { s += rate(i, j, k, kBuildingHeatRefX) * (t(i + di, j, k) - tb(i + di, j, k)); }
        if (rate(i, j, k, kBuildingHeatRefY) != amrex::Real(0)) { s += rate(i, j, k, kBuildingHeatRefY) * (t(i, j + dj, k) - tb(i, j + dj, k)); }
        if (rate(i, j, k, kBuildingHeatRefZ) != amrex::Real(0)) { s += rate(i, j, k, kBuildingHeatRefZ) * (t(i, j, k + dk) - tb(i, j, k + dk)); }
        return s;
    }
};

} // namespace

void BuildingTheta::build (const Buildings& b, const BuildingForcing& law, const BuildingsOptions& opt, const amrex::Geometry& geom,
                           const amrex::MultiFab& z_nd, const amrex::MultiFab& detJ, const ReferenceProfile& ref)
{
    if (!b.active()) { amrex::Abort("gust: BuildingTheta::build needs the buildings' geometry"); }
    if (!law.side_law()) { amrex::Abort("gust: BuildingTheta::build needs the wall law on every wall (buildings.wall_form = fraction_stress)"); }
    if (opt.theta_mode == BuildingThetaMode::none) { amrex::Abort("gust: BuildingTheta::build with buildings.theta_mode = none"); }
    m_b = &b;
    m_law = &law;
    m_detJ = &detJ;
    m_geom = geom;
    m_tb.define(b.beta().boxArray(), b.beta().DistributionMap(), 1, 1);
    if (opt.theta_mode == BuildingThetaMode::fixed) {
        m_tb.setVal(static_cast<amrex::Real>(opt.theta_surface));
    } else {
        fill_theta_ref(ref, z_nd, m_tb);
    }
}

void BuildingTheta::close_fluxes (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open,
                                  amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& Gt) const
{
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab& G = Gt[static_cast<std::size_t>(d)];
        for (amrex::MFIter mfi(G, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto g = G.array(mfi);
            const auto o = open[static_cast<std::size_t>(d)]->const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { g(i, j, k) *= o(i, j, k); });
        }
    }
}

void BuildingTheta::apply_rows (const amrex::MultiFab& theta, amrex::MultiFab& Rt) const
{
    if (theta.nGrow() < 1) { amrex::Abort("gust: BuildingTheta::apply_rows: theta needs a ghost layer (the walls' reference cells)"); }
    for (amrex::MFIter mfi(Rt, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto r = Rt.array(mfi);
        const auto t = theta.const_array(mfi);
        const ThetaRow c{m_b->beta().const_array(mfi), m_b->wall().const_array(mfi), m_law->rate().const_array(mfi), m_tb.const_array(mfi)};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int kind = c.kind(i, j, k);
            if (kind == 1) { r(i, j, k) = t(i, j, k) - c.tb(i, j, k); }
            if (kind == 2) { r(i, j, k) += c.exchange(t, i, j, k); }
        });
    }
}

void BuildingTheta::add_to_implicit (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open, amrex::MultiFab& A) const
{
    // The couplings through the closed parts of the faces, as close_fluxes: the whole
    // coupling is scaled and a_P keeps the closed part, so the row stays diagonally dominant.
    for (int d = 0; d < 3; ++d) {
        const int ei = (d == 0) ? 1 : 0, ej = (d == 1) ? 1 : 0, ek = (d == 2) ? 1 : 0;
        const int lo = coef_XM + 2 * d, hi = lo + 1;
        for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = A.array(mfi);
            const auto o = open[static_cast<std::size_t>(d)]->const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                a(i, j, k, lo) *= o(i, j, k);
                a(i, j, k, hi) *= o(i + ei, j + ej, k + ek);
            });
        }
    }
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto a = A.array(mfi);
        const ThetaRow c{m_b->beta().const_array(mfi), m_b->wall().const_array(mfi), m_law->rate().const_array(mfi), m_tb.const_array(mfi)};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int kind = c.kind(i, j, k);
            if (kind == 1) {
                for (int n = 0; n < n_implicit_coefs; ++n) { a(i, j, k, n) = 0; }
                a(i, j, k, coef_P) = 1;
            } else if (kind == 2) {
                a(i, j, k, coef_P) += c.rate(i, j, k, kBuildingHeatOwn);
            }
        });
    }
}

void BuildingTheta::zero_solid_cells (amrex::MultiFab& field) const
{
    if (field.boxArray() != m_tb.boxArray() || field.DistributionMap() != m_tb.DistributionMap()) {
        amrex::Abort("gust: BuildingTheta::zero_solid_cells: the field is not on the buildings' cells");
    }
    for (amrex::MFIter mfi(field, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto f = field.array(mfi);
        const auto beta = m_b->beta().const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (beta(i, j, k) >= amrex::Real(1)) { f(i, j, k) = 0; }
        });
    }
}

amrex::Real BuildingTheta::wall_heat (const amrex::MultiFab& theta) const
{
    if (theta.nGrow() < 1) { amrex::Abort("gust: BuildingTheta::wall_heat: theta needs a ghost layer (the walls' reference cells)"); }
    // The cells' exchange times their volume, summed in a fixed order (deterministic_sum).
    amrex::MultiFab q(m_tb.boxArray(), m_tb.DistributionMap(), 1, 0);
    const auto cs = m_geom.CellSizeArray();
    const amrex::Real cell = cs[0] * cs[1] * cs[2];
    for (amrex::MFIter mfi(q, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto h = q.array(mfi);
        const auto t = theta.const_array(mfi);
        const auto J = m_detJ->const_array(mfi);
        const ThetaRow c{m_b->beta().const_array(mfi), m_b->wall().const_array(mfi), m_law->rate().const_array(mfi), m_tb.const_array(mfi)};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            h(i, j, k) = c.kind(i, j, k) == 2 ? c.exchange(t, i, j, k) * J(i, j, k) * cell : amrex::Real(0);
        });
    }
    return deterministic_sum(q, 0);
}

} // namespace gust
