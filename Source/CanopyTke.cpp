#include "CanopyTke.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <cmath>

namespace gust {

void add_canopy_tke (const State& q, const amrex::MultiFab& f, amrex::Real beta_p, amrex::Real beta_d,
                     amrex::Real tke_floor, ClosureFields& cf)
{
    if (f.boxArray() != q.p.boxArray() || f.DistributionMap() != q.p.DistributionMap()) {
        amrex::Abort("gust: add_canopy_tke: the drag field is not on the state's cells");
    }
    if (cf.Pc.empty()) {
        cf.Pc.define(q.p.boxArray(), q.p.DistributionMap(), 1, 0);
        cf.Dc.define(q.p.boxArray(), q.p.DistributionMap(), 1, 0);
    }
    for (amrex::MFIter mfi(cf.Pc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto pc = cf.Pc.array(mfi), dc = cf.Dc.array(mfi), ds = cf.dSdk.array(mfi), de = cf.dSdk_exact.array(mfi);
        const auto fa = f.const_array(mfi);
        const auto u = q.u.const_array(mfi), v = q.v.const_array(mfi), w = q.w.const_array(mfi);
        const auto kk = q.tke.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real half = amrex::Real(0.5);
            const amrex::Real uc = half * (u(i, j, k) + u(i + 1, j, k));
            const amrex::Real vc = half * (v(i, j, k) + v(i, j + 1, k));
            const amrex::Real wc = half * (w(i, j, k) + w(i, j, k + 1));
            const amrex::Real wind = std::sqrt(uc * uc + vc * vc + wc * wc);
            const amrex::Real fd = fa(i, j, k);
            const bool above = kk(i, j, k) > tke_floor;
            const amrex::Real kf = above ? kk(i, j, k) : tke_floor;
            const amrex::Real c = beta_d * fd * wind;   // the sink's rate (1/s)
            pc(i, j, k) = beta_p * fd * wind * wind * wind;
            dc(i, j, k) = c * kf;
            ds(i, j, k) += c;
            de(i, j, k) -= above ? c : amrex::Real(0);
        });
    }
}

namespace {

template <class Z0>
amrex::Long count_rough (const amrex::MultiFab& f, const Roughness& rough, amrex::Real z0_max)
{
    amrex::ReduceOps<amrex::ReduceOpSum> rop;
    amrex::ReduceData<amrex::Long> data(rop);
    for (amrex::MFIter mfi(f, false); mfi.isValid(); ++mfi) {
        amrex::Box bx = mfi.validbox();
        bx.setBig(2, bx.smallEnd(2));   // the first cell of each column (the boxes span the column)
        const auto fa = f.const_array(mfi);
        const Z0 z0 = Z0::of(rough, mfi);
        rop.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
        {
            return {(fa(i, j, k) > 0 && z0(i, j) >= z0_max) ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    amrex::Long n = amrex::get<0>(data.value(rop));
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

} // namespace

amrex::Long canopy_columns_on_rough_ground (const amrex::MultiFab& f, const Roughness& rough, amrex::Real z0_max)
{
    return (rough.map != nullptr) ? count_rough<MapZ0>(f, rough, z0_max) : count_rough<UniformZ0>(f, rough, z0_max);
}

} // namespace gust
