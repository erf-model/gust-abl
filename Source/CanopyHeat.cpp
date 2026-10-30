#include "CanopyHeat.H"

#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <cmath>

namespace gust {

namespace {

// what (named in the message) must be on the leaf area density's cells; with face set, on its faces.
void check_layout (const char* op, const char* what, const amrex::MultiFab& lad, const amrex::MultiFab& f, bool face = false)
{
    const amrex::BoxArray cells = face ? amrex::convert(f.boxArray(), amrex::IntVect::TheCellVector()) : f.boxArray();
    if (f.empty() || cells != lad.boxArray() || f.DistributionMap() != lad.DistributionMap()) {
        amrex::Abort(std::string("gust: ") + op + ": " + what + " is not on the leaf area density's layout");
    }
}

// The cell-centred wind of a cell from its own faces.
struct CellWind
{
    amrex::Array4<amrex::Real const> u, v, w;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real operator() (int i, int j, int k) const noexcept
    {
        const amrex::Real half = amrex::Real(0.5);
        const amrex::Real uc = half * (u(i, j, k) + u(i + 1, j, k));
        const amrex::Real vc = half * (v(i, j, k) + v(i, j + 1, k));
        const amrex::Real wc = half * (w(i, j, k) + w(i, j, k + 1));
        return std::sqrt(uc * uc + vc * vc + wc * wc);
    }
};

} // namespace

void CanopyHeat::build (const Canopy& canopy)
{
    if (!canopy.active() || !canopy.heat().on()) { amrex::Abort("gust: CanopyHeat::build needs a canopy with canopy.heat on"); }
    m_lad = &canopy.lad();
    m_leaf_theta = static_cast<amrex::Real>(canopy.heat().leaf_theta);
    m_frozen = canopy.heat().frozen;
}

void add_canopy_heat (const CanopyHeat& heat, const State& q, amrex::MultiFab& R_theta)
{
    check_layout("add_canopy_heat", "theta", heat.lad(), q.theta);
    check_layout("add_canopy_heat", "R_theta", heat.lad(), R_theta);
    check_layout("add_canopy_heat", "u", heat.lad(), q.u, true);
    check_layout("add_canopy_heat", "v", heat.lad(), q.v, true);
    check_layout("add_canopy_heat", "w", heat.lad(), q.w, true);
    const amrex::Real tl = heat.leaf_theta();
    for (amrex::MFIter mfi(R_theta, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto r = R_theta.array(mfi);
        const auto a = heat.lad().const_array(mfi);
        const auto th = q.theta.const_array(mfi);
        const CellWind wind{q.u.const_array(mfi), q.v.const_array(mfi), q.w.const_array(mfi)};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r(i, j, k) += leaf_exchange_rate(a(i, j, k), wind(i, j, k)) * (th(i, j, k) - tl);
        });
    }
}

void canopy_heat_rows (const CanopyHeat& heat, const amrex::MultiFab& u, const amrex::MultiFab& v,
                       const amrex::MultiFab& w, CanopyRows& rows)
{
    const amrex::MultiFab& lad = heat.lad();
    check_layout("canopy_heat_rows", "u", lad, u, true);
    check_layout("canopy_heat_rows", "v", lad, v, true);
    check_layout("canopy_heat_rows", "w", lad, w, true);
    if (rows.theta.empty() || rows.theta.boxArray() != lad.boxArray() || rows.theta.DistributionMap() != lad.DistributionMap()) {
        rows.theta.define(lad.boxArray(), lad.DistributionMap(), 1, 0);
    }
    for (amrex::MFIter mfi(rows.theta, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto out = rows.theta.array(mfi);
        const auto a = lad.const_array(mfi);
        const CellWind wind{u.const_array(mfi), v.const_array(mfi), w.const_array(mfi)};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            out(i, j, k) = leaf_exchange_rate(a(i, j, k), wind(i, j, k));
        });
    }
}

amrex::Real blend_frozen_theta (const CanopyHeat& heat, const amrex::MultiFab& theta_f, const amrex::MultiFab& eps,
                                amrex::Real tke_floor, State& q)
{
    check_layout("blend_frozen_theta", "theta", heat.lad(), q.theta);
    check_layout("blend_frozen_theta", "k", heat.lad(), q.tke);
    check_layout("blend_frozen_theta", "the frozen theta", heat.lad(), theta_f);
    check_layout("blend_frozen_theta", "eps (the closure's)", heat.lad(), eps);
    check_layout("blend_frozen_theta", "u", heat.lad(), q.u, true);
    check_layout("blend_frozen_theta", "v", heat.lad(), q.v, true);
    check_layout("blend_frozen_theta", "w", heat.lad(), q.w, true);
    const amrex::Real tl = heat.leaf_theta();
    amrex::ReduceOps<amrex::ReduceOpMax> rop;
    amrex::ReduceData<amrex::Real> data(rop);
    for (amrex::MFIter mfi(q.theta, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto th = q.theta.array(mfi);
        const auto tf = theta_f.const_array(mfi);
        const auto e = eps.const_array(mfi);
        const auto kk = q.tke.const_array(mfi);
        const auto a = heat.lad().const_array(mfi);
        const CellWind wind{q.u.const_array(mfi), q.v.const_array(mfi), q.w.const_array(mfi)};
        rop.eval(mfi.tilebox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real>
        {
            const amrex::Real lam = leaf_exchange_rate(a(i, j, k), wind(i, j, k));
            const amrex::Real kf = kk(i, j, k) > tke_floor ? kk(i, j, k) : tke_floor;
            const amrex::Real wt = frozen_blend_weight(lam, kf, e(i, j, k));
            const amrex::Real next = tf(i, j, k) + wt * (tl - tf(i, j, k));
            const amrex::Real change = std::abs(next - th(i, j, k));
            th(i, j, k) = next;
            return {change};
        });
    }
    amrex::Real m = amrex::get<0>(data.value(rop));
    amrex::ParallelDescriptor::ReduceRealMax(m);
    return m;
}

} // namespace gust
