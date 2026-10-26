#include "Reductions.H"

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_GpuContainers.H>
#include <AMReX_ParallelDescriptor.H>

#include <bit>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace gust {

amrex::Real deterministic_sum (const amrex::MultiFab& mf, int comp)
{
    // Face and nodal boxes share their boundary points, which would be counted twice.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf.is_cell_centered(),
                                     "deterministic_sum needs cell-centred data");
    const int nboxes = static_cast<int>(mf.boxArray().size());
    amrex::Vector<amrex::Real> box_sum(nboxes, amrex::Real(0));

    // No tiling: one partial sum per box, whatever the tile size.
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const int klo = bx.smallEnd(2);
        const int khi = bx.bigEnd(2);
        amrex::Box cols = bx;
        cols.setBig(2, klo);

        amrex::FArrayBox colsum(cols, 1, amrex::The_Async_Arena());
        const auto cs = colsum.array();
        const auto a = mf.const_array(mfi, comp);
        amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int k0) noexcept
        {
            amrex::Real s = 0;
            for (int k = klo; k <= khi; ++k) { s += a(i, j, k); }
            cs(i, j, k0) = s;
        });

        const auto npts = static_cast<std::size_t>(cols.numPts());
        amrex::Vector<amrex::Real> h(npts);
        amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, colsum.dataPtr(),
                              colsum.dataPtr() + npts, h.begin());
        amrex::Gpu::streamSynchronize();

        amrex::Real s = 0;
        for (const auto x : h) { s += x; }
        box_sum[mfi.index()] = s;
    }

    amrex::ParallelDescriptor::ReduceRealSum(box_sum.dataPtr(), nboxes);

    amrex::Real total = 0;
    for (const auto x : box_sum) { total += x; }
    return total;
}

double box_dot (const amrex::MultiFab& a, int acomp, const amrex::MultiFab& b, int bcomp, const amrex::iMultiFab* owner,
                const amrex::MFIter& mfi)
{
    // mfi may iterate over another index type (the cells); the field's own box of the same index.
    const amrex::Box bx = amrex::convert(mfi.validbox(), a.ixType());
    const int klo = bx.smallEnd(2);
    const int khi = bx.bigEnd(2);
    amrex::Box cols = bx;
    cols.setBig(2, klo);
    const auto npts = static_cast<std::size_t>(cols.numPts());
    amrex::Gpu::DeviceVector<double> colsum(npts);
    double* cs = colsum.data();
    const auto lo = amrex::lbound(cols);
    const auto len = amrex::length(cols);
    const auto aa = a.const_array(mfi, acomp);
    const auto ba = b.const_array(mfi, bcomp);
    if (owner != nullptr) {
        const auto m = owner->const_array(mfi);
        amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int k0) noexcept
        {
            double s = 0;
            for (int k = klo; k <= khi; ++k) {
                if (m(i, j, k) != 0) { s += static_cast<double>(aa(i, j, k)) * static_cast<double>(ba(i, j, k)); }
            }
            cs[(i - lo.x) + (j - lo.y) * len.x + (k0 - lo.z) * len.x * len.y] = s;
        });
    } else {
        amrex::ParallelFor(cols, [=] AMREX_GPU_DEVICE (int i, int j, int k0) noexcept
        {
            double s = 0;
            for (int k = klo; k <= khi; ++k) { s += static_cast<double>(aa(i, j, k)) * static_cast<double>(ba(i, j, k)); }
            cs[(i - lo.x) + (j - lo.y) * len.x + (k0 - lo.z) * len.x * len.y] = s;
        });
    }
    amrex::Vector<double> h(npts);
    amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, colsum.begin(), colsum.end(), h.begin());
    amrex::Gpu::streamSynchronize();
    double s = 0;
    for (const double x : h) { s += x; }
    return s;
}

double deterministic_dot (const amrex::MultiFab& a, int acomp, const amrex::MultiFab& b, int bcomp,
                          const amrex::iMultiFab* owner)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(a.boxArray() == b.boxArray() && a.DistributionMap() == b.DistributionMap(),
                                     "deterministic_dot needs two fields on one layout");
    const int nboxes = static_cast<int>(a.boxArray().size());
    amrex::Vector<double> box_sum(nboxes, 0.0);
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        box_sum[mfi.index()] = box_dot(a, acomp, b, bcomp, owner, mfi);
    }
    amrex::ParallelDescriptor::ReduceRealSum(box_sum.dataPtr(), nboxes);
    double total = 0;
    for (const double x : box_sum) { total += x; }
    return total;
}

namespace {

//! The local count of count_not_positive for one field (no communication).
//! AllowZero a template argument: no option branch in the kernel.
template <bool AllowZero>
amrex::Long local_not_positive (const amrex::MultiFab& mf, int ngrow)
{
    // The IEEE 754 layout of amrex::Real: the sign bit and the exponent mask.
    using Bits = std::conditional_t<sizeof(amrex::Real) == 8, std::uint64_t, std::uint32_t>;
    static_assert(sizeof(Bits) == sizeof(amrex::Real));
    constexpr Bits kSign = Bits(1) << (8 * sizeof(Bits) - 1);
    constexpr Bits kExp = (sizeof(Bits) == 8) ? Bits(0x7ff0000000000000ULL) : Bits(0x7f800000UL);
    amrex::ReduceOps<amrex::ReduceOpSum> op;
    amrex::ReduceData<amrex::Long> data(op);
    for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(ngrow);
        const auto a = mf.const_array(mfi);
        op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
        {
            // From the bit pattern, integer operations only: no floating-point
            // compare ever sees a NaN, so nothing raises FE_INVALID (trapped in
            // CI). GCC's std::isinf is an ordered compare, evaluated even when
            // isnan is already true.
            const auto bits = std::bit_cast<Bits>(a(i, j, k));
            const Bits mag = bits & ~kSign;
            const bool zero = mag == 0;                            // +0 or -0
            bool bad = mag >= kExp || ((bits & kSign) != 0 && !zero);   // NaN, +-inf, or negative
            if constexpr (!AllowZero) { bad = bad || zero; }
            return {bad ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    return amrex::get<0>(data.value(op));
}

} // namespace

std::vector<amrex::Long> count_not_positive (const std::vector<const amrex::MultiFab*>& fields, bool allow_zero,
                                             int ngrow)
{
    std::vector<amrex::Long> n;
    for (const amrex::MultiFab* mf : fields) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf->nGrowVect().allGE(amrex::IntVect(ngrow)),
                                         "count_not_positive: the field has fewer ghost layers than asked");
        n.push_back(allow_zero ? local_not_positive<true>(*mf, ngrow) : local_not_positive<false>(*mf, ngrow));
    }
    if (!n.empty()) { amrex::ParallelDescriptor::ReduceLongSum(n.data(), static_cast<int>(n.size())); }
    return n;
}

void require_positive (const char* op, const std::vector<std::string>& names,
                       const std::vector<const amrex::MultiFab*>& fields, bool allow_zero, int ngrow)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(names.size() == fields.size(), "require_positive: one name per field");
    const std::vector<amrex::Long> n = count_not_positive(fields, allow_zero, ngrow);
    std::string bad;
    amrex::Long total = 0;
    for (std::size_t i = 0; i < n.size(); ++i) {
        if (n[i] == 0) { continue; }
        bad += (bad.empty() ? "" : ", ") + names[i];
        total += n[i];
    }
    if (total != 0) {
        amrex::Abort(std::string("gust: ") + op + ": " + bad + " must be finite and " + (allow_zero ? ">= 0" : "> 0")
                     + " everywhere (" + std::to_string(total) + " points are not)");
    }
}

std::vector<double> deterministic_level_sums (const amrex::MultiFab& mf, int comp, const amrex::iMultiFab& owner,
                                              const amrex::Box& domain)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(mf.ixType().cellCentered(2), "deterministic_level_sums needs cell data in z");
    const int klo_d = domain.smallEnd(2);
    const int nz = domain.length(2);
    const int nboxes = static_cast<int>(mf.boxArray().size());
    amrex::Vector<double> slots(static_cast<std::size_t>(nboxes) * static_cast<std::size_t>(nz), 0.0);
    for (amrex::MFIter mfi(mf, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox();
        const int ilo = bx.smallEnd(0);
        const int ihi = bx.bigEnd(0);
        amrex::Box rows = bx;
        rows.setBig(0, ilo);
        const auto npts = static_cast<std::size_t>(rows.numPts());
        amrex::Gpu::DeviceVector<double> rowsum(npts);
        double* rs = rowsum.data();
        const auto lo = amrex::lbound(rows);
        const auto len = amrex::length(rows);
        const auto a = mf.const_array(mfi, comp);
        const auto m = owner.const_array(mfi);
        amrex::ParallelFor(rows, [=] AMREX_GPU_DEVICE (int i0, int j, int k) noexcept
        {
            double s = 0;
            for (int i = ilo; i <= ihi; ++i) {
                if (m(i, j, k) != 0) { s += static_cast<double>(a(i, j, k)); }
            }
            rs[(i0 - lo.x) + (j - lo.y) * len.x + (k - lo.z) * len.x * len.y] = s;
        });
        amrex::Vector<double> h(npts);
        amrex::Gpu::copyAsync(amrex::Gpu::deviceToHost, rowsum.begin(), rowsum.end(), h.begin());
        amrex::Gpu::streamSynchronize();
        // rows has length 1 in x: index (j - lo.y) + (k - lo.z) len.y.
        for (int k = 0; k < len.z; ++k) {
            double s = 0;
            for (int j = 0; j < len.y; ++j) { s += h[static_cast<std::size_t>(j + k * len.y)]; }
            slots[static_cast<std::size_t>(mfi.index()) * nz + (lo.z + k - klo_d)] = s;
        }
    }
    amrex::ParallelDescriptor::ReduceRealSum(slots.dataPtr(), static_cast<int>(slots.size()));
    std::vector<double> sums(static_cast<std::size_t>(nz), 0.0);
    for (int b = 0; b < nboxes; ++b) {
        for (int k = 0; k < nz; ++k) { sums[k] += slots[static_cast<std::size_t>(b) * nz + k]; }
    }
    return sums;
}

namespace {
struct OwnerMaskEntry
{
    amrex::BoxArray ba;               // held, so its identity is not reused
    amrex::DistributionMapping dm;
    amrex::IntVect period;
    std::unique_ptr<amrex::iMultiFab> mask;
};
std::vector<OwnerMaskEntry>& owner_masks ()
{
    static std::vector<OwnerMaskEntry> cache;
    return cache;
}
} // namespace

const amrex::iMultiFab& owner_mask (const amrex::MultiFab& mf, const amrex::Periodicity& period)
{
    auto& cache = owner_masks();
    for (const auto& e : cache) {
        if (e.ba == mf.boxArray() && e.dm == mf.DistributionMap() && e.period == period.intVect()) { return *e.mask; }
    }
    // Four index types per level: eight entries hold the level in use and the
    // one before (grid sequencing); the oldest goes first.
    if (cache.size() >= 8) { cache.erase(cache.begin()); }
    // The masks live in AMReX's arenas: freed before amrex::Finalize (once per
    // filling of an empty cache; a second clear of an empty cache is harmless).
    if (cache.empty()) { amrex::ExecOnFinalize([] { owner_masks().clear(); }); }
    cache.push_back({mf.boxArray(), mf.DistributionMap(), period.intVect(), mf.OwnerMask(period)});
    return *cache.back().mask;
}

amrex::Real global_max (const amrex::MultiFab& mf, int comp)
{
    // max(x, NaN) keeps x, so a NaN would be kept or dropped depending on the
    // decomposition; test for it first (contains_nan is collective).
    if (mf.contains_nan(comp, 1, 0)) {
        return std::numeric_limits<amrex::Real>::quiet_NaN();
    }
    const amrex::Real m = mf.max(comp, 0);
    // max(+0, -0) returns its first argument; report zero as +0.
    return (m == amrex::Real(0)) ? amrex::Real(0) : m;
}

} // namespace gust
