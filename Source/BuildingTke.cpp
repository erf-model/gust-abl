#include "BuildingTke.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <cmath>

namespace gust {

namespace {

// The components of the working field of the distance: the distance and the nearest
// wall's rectangle (centre and half sizes).
enum { kD = 0, kCx, kCy, kCz, kHx, kHy, kHz, kNumWork };

//! A cell's geometry from the node heights: centre and the mean heights of its floor and top.
struct CellGeom
{
    amrex::Array4<const amrex::Real> z;
    amrex::Real x0, y0, dx, dy;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real zlo (int i, int j, int k) const noexcept
    {
        return amrex::Real(0.25) * (z(i, j, k) + z(i + 1, j, k) + z(i, j + 1, k) + z(i + 1, j + 1, k));
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real xc (int i) const noexcept { return x0 + (i + amrex::Real(0.5)) * dx; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real yc (int j) const noexcept { return y0 + (j + amrex::Real(0.5)) * dy; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real zc (int i, int j, int k) const noexcept { return amrex::Real(0.5) * (zlo(i, j, k) + zlo(i, j, k + 1)); }
};

//! The rectangle of the wall `bit` of the cell (i, j, k), solid fraction beta: its plane lies
//! beta of the cell's size in from the face the wall is on. A side covers the part f of its
//! face's height, from the cell's floor (from_floor: the top of a building) or down from its
//! top (the lower end of a lifted one).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wall_rectangle (const CellGeom& g, int i, int j, int k, int bit, amrex::Real beta, amrex::Real f, bool from_floor,
                     amrex::Real c[3], amrex::Real h[3]) noexcept
{
    const amrex::Real half = amrex::Real(0.5);
    const amrex::Real zl = g.zlo(i, j, k), zh = g.zlo(i, j, k + 1), dz = zh - zl;
    c[0] = g.xc(i);
    c[1] = g.yc(j);
    c[2] = from_floor ? zl + half * f * dz : zh - half * f * dz;
    h[0] = half * g.dx;
    h[1] = half * g.dy;
    h[2] = half * f * dz;
    if (bit == kWallRoof) { c[2] = zl + beta * dz; h[2] = 0; }
    else if (bit == kWallUnder) { c[2] = zh - beta * dz; h[2] = 0; }
    else if (bit == kWallXlo) { c[0] += (beta - half) * g.dx; h[0] = 0; }
    else if (bit == kWallXhi) { c[0] -= (beta - half) * g.dx; h[0] = 0; }
    else if (bit == kWallYlo) { c[1] += (beta - half) * g.dy; h[1] = 0; }
    else { c[1] -= (beta - half) * g.dy; h[1] = 0; }
}

} // namespace

void BuildingTke::build (const Buildings& b, const BuildingForcing& law, const BuildingsOptions& opt, const amrex::Geometry& geom,
                         const amrex::MultiFab& z_nd, amrex::Real tke_floor, amrex::Real Cmu0)
{
    if (!b.active()) { amrex::Abort("gust: BuildingTke::build needs the buildings' geometry"); }
    if (!law.side_law()) { amrex::Abort("gust: BuildingTke::build needs the wall law on every wall (buildings.wall_form = fraction_stress)"); }
    m_b = &b;
    m_law = &law;
    m_geom = geom;
    m_floor = tke_floor;
    m_cmu0 = Cmu0;
    const amrex::BoxArray& ba = b.beta().boxArray();
    const amrex::DistributionMapping& dm = b.beta().DistributionMap();
    const amrex::Box& dom = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto cs = geom.CellSizeArray();
    const amrex::Real lx = geom.isPeriodic(0) ? geom.ProbLength(0) : amrex::Real(0);
    const amrex::Real ly = geom.isPeriodic(1) ? geom.ProbLength(1) : amrex::Real(0);
    const amrex::Real z0b = static_cast<amrex::Real>(opt.z0), far = kBuildingFar;
    // Only min(d_ground, d_b) is used, and the ground is never farther than the tallest
    // column: a wall at that distance or more is not passed on, which bounds the passes by
    // that height in cells.
    amrex::Real cutoff = 0;
    {
        amrex::ReduceOps<amrex::ReduceOpMax> op;
        amrex::ReduceData<amrex::Real> data(op);
        const int k0 = dom.smallEnd(2), k1 = dom.bigEnd(2) + 1;
        for (amrex::MFIter mfi(b.beta(), false); mfi.isValid(); ++mfi) {
            const auto z = z_nd.const_array(mfi);
            amrex::Box col = amrex::surroundingNodes(mfi.validbox());
            col.setBig(2, col.smallEnd(2));
            op.eval(col, data, [=] AMREX_GPU_DEVICE (int i, int j, int) -> amrex::GpuTuple<amrex::Real> { return {z(i, j, k1) - z(i, j, k0)}; });
        }
        cutoff = amrex::get<0>(data.value(op));
        amrex::ParallelDescriptor::ReduceRealMax(cutoff);
    }
    m_cutoff = cutoff;
    // In a periodic direction every index has an image; else the domain's own.
    const int ilo = geom.isPeriodic(0) ? dom.smallEnd(0) - 1 : dom.smallEnd(0), ihi = geom.isPeriodic(0) ? dom.bigEnd(0) + 1 : dom.bigEnd(0);
    const int jlo = geom.isPeriodic(1) ? dom.smallEnd(1) - 1 : dom.smallEnd(1), jhi = geom.isPeriodic(1) ? dom.bigEnd(1) + 1 : dom.bigEnd(1);
    const int klo = dom.smallEnd(2), khi = dom.bigEnd(2);

    amrex::MultiFab work[2];
    for (auto& w : work) {
        w.define(ba, dm, kNumWork, 1);
        w.setVal(0);
        w.setVal(far, kD, 1, 1);
    }
    // The sides' blocked fractions with a ghost layer (Buildings::beta_face has none).
    amrex::MultiFab fside[2];
    for (int d = 0; d < 2; ++d) {
        fside[d].define(b.beta_face(d).boxArray(), dm, 1, 1);
        fside[d].setVal(0);
        amrex::MultiFab::Copy(fside[d], b.beta_face(d), 0, 0, 1, 0);
        fside[d].FillBoundary(geom.periodicity());
    }
    // The seeds: each cell's nearest wall among those of its 27 cells. A wall cell keeps
    // the distance of its own fluid centroid.
    for (amrex::MFIter mfi(work[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto w = work[0].array(mfi);
        const auto beta = b.beta().const_array(mfi);
        const auto wl = b.wall().const_array(mfi);
        const auto fx = fside[0].const_array(mfi);
        const auto fy = fside[1].const_array(mfi);
        const CellGeom g{z_nd.const_array(mfi), plo[0], plo[1], cs[0], cs[1]};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (beta(i, j, k) >= amrex::Real(1)) { return; }
            const amrex::Real px = g.xc(i), py = g.yc(j), pz = g.zc(i, j, k);
            amrex::Real best = far, bc[3] = {0, 0, 0}, bh[3] = {0, 0, 0};
            for (int kk = amrex::max(k - 1, klo); kk <= amrex::min(k + 1, khi); ++kk) {
                for (int jj = amrex::max(j - 1, jlo); jj <= amrex::min(j + 1, jhi); ++jj) {
                    for (int ii = amrex::max(i - 1, ilo); ii <= amrex::min(i + 1, ihi); ++ii) {
                        const int mask = static_cast<int>(wl(ii, jj, kk, kWallMask));
                        if (mask == 0 || beta(ii, jj, kk) >= amrex::Real(1)) { continue; }
                        for (int bit = kWallRoof; bit <= kWallUnder; bit *= 2) {
                            if ((mask & bit) == 0) { continue; }
                            // A side's blocked part: from the floor when the face below is blocked at
                            // least as much (or is the ground), else from the top.
                            amrex::Real f = 1, below = 1;
                            if (bit == kWallXlo || bit == kWallXhi) {
                                const int fi = (bit == kWallXlo) ? ii : ii + 1;
                                f = fx(fi, jj, kk);
                                below = kk > klo ? fx(fi, jj, kk - 1) : amrex::Real(1);
                            } else if (bit == kWallYlo || bit == kWallYhi) {
                                const int fj = (bit == kWallYlo) ? jj : jj + 1;
                                f = fy(ii, fj, kk);
                                below = kk > klo ? fy(ii, fj, kk - 1) : amrex::Real(1);
                            }
                            amrex::Real c[3], h[3];
                            wall_rectangle(g, ii, jj, kk, bit, beta(ii, jj, kk), f, below >= f, c, h);
                            const amrex::Real d = rectangle_distance(px, py, pz, c[0], c[1], c[2], h[0], h[1], h[2], lx, ly);
                            if (d < best && d < cutoff) {
                                best = d;
                                for (int n = 0; n < 3; ++n) { bc[n] = c[n]; bh[n] = h[n]; }
                            }
                        }
                    }
                }
            }
            if (static_cast<int>(wl(i, j, k, kWallMask)) != 0) {
                amrex::Real d1 = far;
                for (int n = kWallD1x; n <= kWallD1z; n += 2) {
                    if (wl(i, j, k, n) > amrex::Real(0)) { d1 = amrex::min(d1, wl(i, j, k, n)); }
                }
                best = amrex::max(d1, z0b);
            }
            w(i, j, k, kD) = best;
            for (int n = 0; n < 3; ++n) { w(i, j, k, kCx + n) = bc[n]; w(i, j, k, kHx + n) = bh[n]; }
        });
    }
    // Pass each cell's nearest wall to its neighbours until no distance falls.
    m_passes = 0;
    int cur = 0;
    for (;;) {
        amrex::MultiFab& src = work[cur];
        amrex::MultiFab& dst = work[1 - cur];
        src.FillBoundary(geom.periodicity());
        amrex::ReduceOps<amrex::ReduceOpSum> op;
        amrex::ReduceData<amrex::Long> data(op);
        for (amrex::MFIter mfi(dst, false); mfi.isValid(); ++mfi) {
            const auto s = src.const_array(mfi);
            const auto t = dst.array(mfi);
            const auto beta = b.beta().const_array(mfi);
            const auto wl = b.wall().const_array(mfi);
            const CellGeom g{z_nd.const_array(mfi), plo[0], plo[1], cs[0], cs[1]};
            op.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
            {
                for (int n = 0; n < kNumWork; ++n) { t(i, j, k, n) = s(i, j, k, n); }
                // Solid cells have no distance; wall cells keep their own.
                if (beta(i, j, k) >= amrex::Real(1) || static_cast<int>(wl(i, j, k, kWallMask)) != 0) { return {amrex::Long(0)}; }
                const amrex::Real px = g.xc(i), py = g.yc(j), pz = g.zc(i, j, k);
                amrex::Real best = s(i, j, k, kD);
                int bi = i, bj = j, bk = k;
                for (int kk = amrex::max(k - 1, klo); kk <= amrex::min(k + 1, khi); ++kk) {
                    for (int jj = j - 1; jj <= j + 1; ++jj) {
                        for (int ii = i - 1; ii <= i + 1; ++ii) {
                            if (s(ii, jj, kk, kD) >= far) { continue; }
                            const amrex::Real d = rectangle_distance(px, py, pz, s(ii, jj, kk, kCx), s(ii, jj, kk, kCy), s(ii, jj, kk, kCz),
                                                                     s(ii, jj, kk, kHx), s(ii, jj, kk, kHy), s(ii, jj, kk, kHz), lx, ly);
                            // A strict fall, beyond rounding, so that ties do not pass back and forth.
                            if (d < best * (amrex::Real(1) - amrex::Real(1.0e-6)) && d < cutoff) { best = d; bi = ii; bj = jj; bk = kk; }
                        }
                    }
                }
                if (bi == i && bj == j && bk == k) { return {amrex::Long(0)}; }
                t(i, j, k, kD) = best;
                for (int n = kCx; n < kNumWork; ++n) { t(i, j, k, n) = s(bi, bj, bk, n); }
                return {amrex::Long(1)};
            });
        }
        amrex::Long changed = amrex::get<0>(data.value(op));
        amrex::ParallelDescriptor::ReduceLongSum(changed);
        cur = 1 - cur;
        ++m_passes;
        if (changed == 0) { break; }
    }
    m_dist.define(ba, dm, 1, 1);
    m_dist.setVal(far);
    amrex::MultiFab::Copy(m_dist, work[cur], kD, 0, 1, 0);
    m_dist.FillBoundary(geom.periodicity());
    // Outside a non-periodic side: the nearest domain column's, as the ground's distance there.
    for (amrex::MFIter mfi(m_dist, false); mfi.isValid(); ++mfi) {
        const auto d = m_dist.array(mfi);
        const amrex::Box vb = mfi.validbox();
        amrex::Box gb = amrex::grow(vb, 0, 1);
        gb.grow(1, 1);
        const int i0 = dom.smallEnd(0), i1 = dom.bigEnd(0), j0 = dom.smallEnd(1), j1 = dom.bigEnd(1);
        const bool px = geom.isPeriodic(0), py = geom.isPeriodic(1);
        amrex::ParallelFor(gb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int ic = px ? i : amrex::min(amrex::max(i, i0), i1), jc = py ? j : amrex::min(amrex::max(j, j0), j1);
            if (ic != i || jc != j) { d(i, j, k) = d(ic, jc, k); }
        });
    }
}

amrex::Long BuildingTke::limit_wall_distance (amrex::MultiFab& dist) const
{
    if (dist.boxArray() != m_dist.boxArray() || dist.DistributionMap() != m_dist.DistributionMap() || dist.nGrow() < 1) {
        amrex::Abort("gust: BuildingTke::limit_wall_distance: the wall distance is not on the buildings' cells with a ghost layer");
    }
    amrex::ReduceOps<amrex::ReduceOpSum> op;
    amrex::ReduceData<amrex::Long> data(op);
    for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
        const auto d = dist.array(mfi);
        const auto db = m_dist.const_array(mfi);
        const amrex::Box vb = mfi.validbox();
        op.eval(amrex::grow(vb, 1), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
        {
            if (!(db(i, j, k) < d(i, j, k))) { return {amrex::Long(0)}; }
            d(i, j, k) = db(i, j, k);
            return {vb.contains(i, j, k) ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    amrex::Long n = amrex::get<0>(data.value(op));
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

void BuildingTke::close_fluxes (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open,
                                amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>& Gk) const
{
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab& G = Gk[static_cast<std::size_t>(d)];
        for (amrex::MFIter mfi(G, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto g = G.array(mfi);
            const auto o = open[static_cast<std::size_t>(d)]->const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { g(i, j, k) *= o(i, j, k); });
        }
    }
}

namespace {

//! The cells of bx in the ground's first layer (slab) and above it (rest); either may be empty.
void split_at_ground (const amrex::Box& bx, int klo, amrex::Box& slab, amrex::Box& rest)
{
    slab = bx;
    rest = bx;
    slab.setBig(2, amrex::min(bx.bigEnd(2), klo));
    rest.setSmall(2, amrex::max(bx.smallEnd(2), klo + 1));
}

//! What a cell's row needs: the faces' blocked fractions, the law's u*^2 and the distance.
struct RowData
{
    amrex::Array4<const amrex::Real> beta, wl, fx, fy, fz, us2, dw;
    amrex::Real floor, inv_cmu2, us_min;
    //! 0: the cell has no row here; 1: solid; 2: a wall cell.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int kind (int i, int j, int k) const noexcept
    {
        if (beta(i, j, k) >= amrex::Real(1)) { return 1; }
        return static_cast<int>(wl(i, j, k, kWallMask)) != 0 ? 2 : 0;
    }
    //! The largest blocked part of the faces the cell's walls lie on.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real blocked (int i, int j, int k) const noexcept
    {
        const int mask = static_cast<int>(wl(i, j, k, kWallMask));
        amrex::Real f = 0;
        if ((mask & kWallRoof) != 0) { f = amrex::max(f, fz(i, j, k)); }
        if ((mask & kWallUnder) != 0) { f = amrex::max(f, fz(i, j, k + 1)); }
        if ((mask & kWallXlo) != 0) { f = amrex::max(f, fx(i, j, k)); }
        if ((mask & kWallXhi) != 0) { f = amrex::max(f, fx(i + 1, j, k)); }
        if ((mask & kWallYlo) != 0) { f = amrex::max(f, fy(i, j, k)); }
        if ((mask & kWallYhi) != 0) { f = amrex::max(f, fy(i, j + 1, k)); }
        return f;
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real rate (int i, int j, int k) const noexcept { return amrex::max(std::sqrt(us2(i, j, k)), us_min) / dw(i, j, k); }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real k_wall (int i, int j, int k) const noexcept { return amrex::max(us2(i, j, k) * inv_cmu2, floor); }
};

} // namespace

void BuildingTke::apply_rows (const amrex::MultiFab& tke, amrex::MultiFab& Rk, BuildingKRow form) const
{
    const int klo = m_geom.Domain().smallEnd(2);
    for (amrex::MFIter mfi(Rk, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto r = Rk.array(mfi);
        const auto t = tke.const_array(mfi);
        const RowData c{m_b->beta().const_array(mfi), m_b->wall().const_array(mfi), m_b->beta_face(0).const_array(mfi),
                        m_b->beta_face(1).const_array(mfi), m_b->beta_face(2).const_array(mfi),
                        m_law->rate().const_array(mfi, kBuildingUstar2), m_dist.const_array(mfi), m_floor,
                        amrex::Real(1) / (m_cmu0 * m_cmu0), kBuildingKUstarMin};
        // blend: the wall's row in the place of the cell's equation, by the blocked part f.
        // add: the wall's row, times f, added to the row the cell has.
        amrex::Box blend = mfi.tilebox(), add;
        if (form == BuildingKRow::ground_row) { split_at_ground(mfi.tilebox(), klo, add, blend); }
        if (form == BuildingKRow::source) { add = mfi.tilebox(); blend = amrex::Box(); }
        if (blend.ok()) {
            amrex::ParallelFor(blend, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int kind = c.kind(i, j, k);
                if (kind == 1) { r(i, j, k) = t(i, j, k) - c.floor; }
                if (kind != 2) { return; }
                const amrex::Real f = c.blocked(i, j, k);
                r(i, j, k) = (amrex::Real(1) - f) * r(i, j, k) + f * c.rate(i, j, k) * (t(i, j, k) - c.k_wall(i, j, k));
            });
        }
        if (add.ok()) {
            amrex::ParallelFor(add, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int kind = c.kind(i, j, k);
                if (kind == 1) { r(i, j, k) = t(i, j, k) - c.floor; }
                if (kind != 2) { return; }
                r(i, j, k) += c.blocked(i, j, k) * c.rate(i, j, k) * (t(i, j, k) - c.k_wall(i, j, k));
            });
        }
    }
}

void BuildingTke::add_to_implicit (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open, amrex::MultiFab& A,
                                   BuildingKRow form) const
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
    const int klo = m_geom.Domain().smallEnd(2);
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto a = A.array(mfi);
        const RowData c{m_b->beta().const_array(mfi), m_b->wall().const_array(mfi), m_b->beta_face(0).const_array(mfi),
                        m_b->beta_face(1).const_array(mfi), m_b->beta_face(2).const_array(mfi),
                        m_law->rate().const_array(mfi, kBuildingUstar2), m_dist.const_array(mfi), m_floor,
                        amrex::Real(1) / (m_cmu0 * m_cmu0), kBuildingKUstarMin};
        amrex::Box blend = mfi.tilebox(), add;
        if (form == BuildingKRow::ground_row) { split_at_ground(mfi.tilebox(), klo, add, blend); }
        if (form == BuildingKRow::source) { add = mfi.tilebox(); blend = amrex::Box(); }
        if (blend.ok()) {
            amrex::ParallelFor(blend, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int kind = c.kind(i, j, k);
                if (kind == 0) { return; }
                const amrex::Real f = (kind == 1) ? amrex::Real(1) : c.blocked(i, j, k);
                for (int n = 0; n < n_implicit_coefs; ++n) { a(i, j, k, n) *= amrex::Real(1) - f; }
                a(i, j, k, coef_P) += (kind == 1) ? amrex::Real(1) : f * c.rate(i, j, k);
            });
        }
        if (add.ok()) {
            amrex::ParallelFor(add, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int kind = c.kind(i, j, k);
                if (kind == 1) {
                    for (int n = 0; n < n_implicit_coefs; ++n) { a(i, j, k, n) = 0; }
                    a(i, j, k, coef_P) = 1;
                } else if (kind == 2) {
                    a(i, j, k, coef_P) += c.blocked(i, j, k) * c.rate(i, j, k);
                }
            });
        }
    }
}

} // namespace gust
