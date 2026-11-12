#include "ClosureWake.H"

#include <cmath>
#include <limits>

namespace gust {

void nut_cap_field (const State& q, const amrex::MultiFab& S2, const ClosureOptions& opt, const ClosureConstants& c,
                    amrex::MultiFab& cap, const amrex::MultiFab* weight)
{
    if (cap.empty()) { cap.define(q.p.boxArray(), q.p.DistributionMap(), 1, 1); }
    if (!(opt.realizability > 0)) {
        cap.setVal(std::numeric_limits<amrex::Real>::max());
        return;
    }
    for (amrex::MFIter mfi(cap, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0));
        const auto cp = cap.array(mfi);
        const auto s2 = S2.const_array(mfi);
        const auto tk = q.tke.const_array(mfi);
        const amrex::Real cr = opt.realizability, floor_k = c.tke_floor;
        const amrex::Real big = std::numeric_limits<amrex::Real>::max(), tiny = std::numeric_limits<amrex::Real>::min();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real s = std::sqrt(s2(i, j, k));
            cp(i, j, k) = (s > tiny) ? cr * amrex::max(tk(i, j, k), floor_k) / s : big;
        });
    }
    if (weight == nullptr) { return; }
    if (weight->boxArray() != cap.boxArray() || weight->DistributionMap() != cap.DistributionMap() || weight->nGrow() < 1) {
        amrex::Abort("gust: nut_cap_field: the cap's weight is not on the closure's cells with a ghost layer");
    }
    // The weight: its own kernel (no option branch in the cap's). cap / w, the largest Real at w = 0 (and where
    // the cap is already that).
    for (amrex::MFIter mfi(cap, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = amrex::grow(mfi.validbox(), amrex::IntVect(1, 1, 0));
        const auto cp = cap.array(mfi);
        const auto w = weight->const_array(mfi);
        const amrex::Real big = std::numeric_limits<amrex::Real>::max();
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real v = w(i, j, k);
            if (v >= amrex::Real(1)) { return; }
            cp(i, j, k) = (v > amrex::Real(0) && cp(i, j, k) < big * v) ? cp(i, j, k) / v : big;
        });
    }
}

} // namespace gust
