#include "BuildingTke.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <cmath>
#include <limits>
#include <string>

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

//! The blocked fractions of a cell's six faces, f[d][0] the low and f[d][1] the high one.
struct FaceFractions
{
    amrex::Real f[3][2];
};

//! How far a wall's plane lies in from the face it is on, as a part of the cell's size along
//! the normal. The solid of a cell is the part that runs through to the opposite face
//! (f_opposite of the cross-section) and the wall's own part (f_wall - f_opposite), which
//! reaches this far:
//!   beta = f_opposite + (f_wall - f_opposite) depth.
//! A slab under a roof gives beta, the side of a neighbouring solid cell 0 (the wall is on the
//! face), whatever a roof in the same cell adds to beta. 0 in a fluid cell.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
amrex::Real wall_depth (amrex::Real beta, amrex::Real f_wall, amrex::Real f_opposite) noexcept
{
    const amrex::Real own = f_wall - f_opposite;
    if (!(beta > amrex::Real(0)) || !(own > amrex::Real(0))) { return 0; }
    return amrex::min(amrex::max((beta - f_opposite) / own, amrex::Real(0)), amrex::Real(1));
}

//! The side (0: low, 1: high) of the cell's one wall along d (0: x, 1: y, 2: z), or -1 with
//! none or both.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int wall_side (int mask, int d) noexcept
{
    const int lo = (d == 0) ? kWallXlo : (d == 1) ? kWallYlo : kWallRoof;
    const int hi = (d == 0) ? kWallXhi : (d == 1) ? kWallYhi : kWallUnder;
    const bool l = (mask & lo) != 0, h = (mask & hi) != 0;
    return (l == h) ? -1 : (l ? 0 : 1);
}

//! The part [a, b] of a cell's size, measured from its low (side 0) or high face, as a centre
//! and half size: lo and hi are the cell's faces.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void set_range (amrex::Real lo, amrex::Real hi, int side, amrex::Real a, amrex::Real b, amrex::Real& c, amrex::Real& h) noexcept
{
    const amrex::Real half = amrex::Real(0.5), size = hi - lo;
    c = (side == 0) ? lo + half * (a + b) * size : hi - half * (a + b) * size;
    h = half * (b - a) * size;
}

//! The rectangle of the wall `bit` of the cell (i, j, k): solid fraction beta, wall bits mask,
//! face fractions ff.
//!
//! Its plane lies wall_depth in from the face the wall is on. Over the plane it covers the
//! blocked part f_wall of that face. The wall's own part is f_wall - f_opposite of it; the
//! rest lies under the solid that runs through the cell (below a roof, behind a tower's
//! side) and is no cell's nearest wall, so it is not cut away:
//! - a side, along z, the blocked part of the heights, from the cell's floor (from_floor:
//!   the top of a building) or down from its top (the lower end of a lifted one). Where the
//!   cell has a wall along the other horizontal direction and no roof, the faces are blocked
//!   along that direction instead (a building's edge inside the cell), and the range is
//!   taken along it, from that wall;
//! - a roof or underside, the footprint: all of it where its own part is the whole face or
//!   the cell is fluid; with one side wall in the cell, the blocked part from that wall;
//! - the corner of a block inside the cell (walls along all three directions, their depths'
//!   product the wall's own part): the block's face.
//! Other cells with walls along three directions (a tower's corner over a podium) keep the
//! whole footprint and the heights' range.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void wall_rectangle (const CellGeom& g, int i, int j, int k, int bit, int mask, amrex::Real beta, const FaceFractions& ff,
                     bool from_floor, amrex::Real c[3], amrex::Real h[3]) noexcept
{
    // tol: the rounding of the fractions (made in double, held as Real) and of their products.
    const amrex::Real half = amrex::Real(0.5), tol = 64 * std::numeric_limits<amrex::Real>::epsilon();
    const amrex::Real zl = g.zlo(i, j, k), zh = g.zlo(i, j, k + 1);
    const amrex::Real lo[3] = {g.xc(i) - half * g.dx, g.yc(j) - half * g.dy, zl};
    const amrex::Real hi[3] = {g.xc(i) + half * g.dx, g.yc(j) + half * g.dy, zh};
    const amrex::Real size[3] = {g.dx, g.dy, zh - zl};
    c[0] = g.xc(i);
    c[1] = g.yc(j);
    c[2] = half * (zl + zh);
    for (int d = 0; d < 3; ++d) { h[d] = half * size[d]; }
    // The wall's normal n and side s; its faces' fractions.
    const int n = (bit == kWallRoof || bit == kWallUnder) ? 2 : (bit == kWallXlo || bit == kWallXhi) ? 0 : 1;
    const int s = (bit == kWallRoof || bit == kWallXlo || bit == kWallYlo) ? 0 : 1;
    const amrex::Real f_wall = ff.f[n][s], f_opp = ff.f[n][1 - s];
    const bool partly = beta > amrex::Real(0);
    // The plane.
    const amrex::Real depth = wall_depth(beta, f_wall, f_opp);
    if (n == 2) { c[2] = (s == 0) ? zl + depth * size[2] : zh - depth * size[2]; }
    else { c[n] += (s == 0) ? (depth - half) * size[n] : -(depth - half) * size[n]; }
    h[n] = 0;
    // The blocked part of the wall's face, from 0 to b along the direction the faces are blocked along.
    const amrex::Real a = 0, b = f_wall;
    const int t1 = (n == 0) ? 1 : 0, t2 = (n == 2) ? 1 : 2;   // the tangents: t2 is z for a side
    const int s1 = partly ? wall_side(mask, t1) : -1, s2 = partly ? wall_side(mask, t2) : -1;
    if (s1 >= 0 && s2 >= 0) {
        const amrex::Real d1 = wall_depth(beta, ff.f[t1][s1], ff.f[t1][1 - s1]);
        const amrex::Real d2 = wall_depth(beta, ff.f[t2][s2], ff.f[t2][1 - s2]);
        if (f_opp <= tol && std::abs(d1 * d2 - f_wall) <= tol) {   // the corner of a block
            set_range(lo[t1], hi[t1], s1, amrex::Real(0), d1, c[t1], h[t1]);
            set_range(lo[t2], hi[t2], s2, amrex::Real(0), d2, c[t2], h[t2]);
            return;
        }
    }
    if (n == 2) {
        if (!partly || f_wall - f_opp >= amrex::Real(1) - tol) { return; }
        if (s1 >= 0 && s2 < 0) { set_range(lo[0], hi[0], s1, a, b, c[0], h[0]); }
        else if (s2 >= 0 && s1 < 0) { set_range(lo[1], hi[1], s2, a, b, c[1], h[1]); }
        return;
    }
    // A side: along the other horizontal direction where the cell has a wall there and neither roof nor underside.
    if (s1 >= 0 && (mask & (kWallRoof | kWallUnder)) == 0) { set_range(lo[t1], hi[t1], s1, a, b, c[t1], h[t1]); }
    else { set_range(lo[2], hi[2], from_floor ? 0 : 1, a, b, c[2], h[2]); }
}

//! A wall on the face of a fluid cell is the outer face of the solid in the cell across it.
//! Where that cell is partly solid and its solid is one block (a building's edge or top
//! inside it), the rectangle is the block's cross-section on the shared face: along each
//! tangent with one wall in that cell, from that wall's side to its depth. mid, hs, lo, hi: the
//! centre, half size and faces of that cell along x and y, and of the fluid cell's own column along z
//! (the two share the face); n: the wall's normal; plane: the shared face's coordinate; beta,
//! mask, ff: the partly solid cell's; f: the shared face's blocked fraction.
//! false, and c, h untouched, where the depths' product is not f: the solid is no single
//! block there, and the caller draws the fluid cell's own rectangle.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool solid_face_rectangle (const amrex::Real mid[3], const amrex::Real hs[3], const amrex::Real lo[3], const amrex::Real hi[3], int n, amrex::Real plane, amrex::Real beta, int mask,
                           const FaceFractions& ff, amrex::Real f, amrex::Real c[3], amrex::Real h[3]) noexcept
{
    const amrex::Real tol = 64 * std::numeric_limits<amrex::Real>::epsilon();
    const int t[2] = {(n == 0) ? 1 : 0, (n == 2) ? 1 : 2};
    int side[2];
    amrex::Real depth[2], prod = 1;
    for (int m = 0; m < 2; ++m) {
        side[m] = wall_side(mask, t[m]);
        depth[m] = side[m] >= 0 ? wall_depth(beta, ff.f[t[m]][side[m]], ff.f[t[m]][1 - side[m]]) : amrex::Real(1);
        prod *= depth[m];
    }
    if (std::abs(prod - f) > tol) { return false; }
    c[n] = plane;
    h[n] = 0;
    for (int m = 0; m < 2; ++m) {
        if (side[m] >= 0) { set_range(lo[t[m]], hi[t[m]], side[m], amrex::Real(0), depth[m], c[t[m]], h[t[m]]); }
        else { c[t[m]] = mid[t[m]]; h[t[m]] = hs[t[m]]; }
    }
    return true;
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
    m_z0b = static_cast<amrex::Real>(opt.z0);
    const amrex::BoxArray& ba = b.beta().boxArray();
    const amrex::DistributionMapping& dm = b.beta().DistributionMap();
    const amrex::Box& dom = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto cs = geom.CellSizeArray();
    const amrex::Real lx = geom.isPeriodic(0) ? geom.ProbLength(0) : amrex::Real(0);
    const amrex::Real ly = geom.isPeriodic(1) ? geom.ProbLength(1) : amrex::Real(0);
    const amrex::Real z0b = m_z0b, far = kBuildingFar;
    // The closure takes the nearer of the ground and a wall (limit_wall_distance), and the
    // ground is never farther than the tallest column: a wall at that distance or more is
    // not passed on, which bounds the passes by that height in cells. (With a ground rougher
    // than the walls a wall just beyond it could still be the nearer by d + z0, by less than
    // the roughnesses' difference; it is not looked for.)
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
    // The faces' blocked fractions, the cells' solid fraction and their wall bits with two ghost
    // layers (Buildings has none and one): periodic images. The seeds below read the domain's
    // cells and those images alone: a wall of a cell one layer out, and the cell across its face.
    amrex::MultiFab fside[3];
    for (int d = 0; d < 3; ++d) {
        fside[d].define(b.beta_face(d).boxArray(), dm, 1, 2);
        fside[d].setVal(0);
        amrex::MultiFab::Copy(fside[d], b.beta_face(d), 0, 0, 1, 0);
        fside[d].FillBoundary(geom.periodicity());
    }
    amrex::MultiFab cells(ba, dm, 2, 2);   // 0: beta, 1: the wall bits
    cells.setVal(0);
    amrex::MultiFab::Copy(cells, b.beta(), 0, 0, 1, 0);
    amrex::MultiFab::Copy(cells, b.wall(), kWallMask, 1, 1, 0);
    cells.FillBoundary(geom.periodicity());
    // The cells with a partly solid neighbour to look at: two image layers in a periodic direction.
    const int ilo2 = geom.isPeriodic(0) ? ilo - 1 : ilo, ihi2 = geom.isPeriodic(0) ? ihi + 1 : ihi;
    const int jlo2 = geom.isPeriodic(1) ? jlo - 1 : jlo, jhi2 = geom.isPeriodic(1) ? jhi + 1 : jhi;
    const amrex::Real ftol = 64 * std::numeric_limits<amrex::Real>::epsilon();
    // The seeds: each cell's nearest wall among those of its 27 cells. A wall cell keeps
    // the distance of its own fluid centroid.
    for (amrex::MFIter mfi(work[0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto w = work[0].array(mfi);
        const auto beta = b.beta().const_array(mfi);
        const auto wl = b.wall().const_array(mfi);
        const auto cl = cells.const_array(mfi);
        const auto fx = fside[0].const_array(mfi);
        const auto fy = fside[1].const_array(mfi);
        const auto fz = fside[2].const_array(mfi);
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
                            // A side's part of its face: from the floor when the face below is blocked at
                            // least as much (or is the ground), else from the top. kk - 1 is two layers
                            // below the cell at most: in the box, as boxes span the column (Buildings::build).
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
                            const FaceFractions ff{{{fx(ii, jj, kk), fx(ii + 1, jj, kk)}, {fy(ii, jj, kk), fy(ii, jj + 1, kk)},
                                                    {fz(ii, jj, kk), fz(ii, jj, kk + 1)}}};
                            amrex::Real c[3], h[3];
                            // A fluid cell's wall on a partly blocked face: the cross-section of the solid in
                            // the cell across it, where that is one block.
                            bool drawn = false;
                            const int n = (bit == kWallRoof || bit == kWallUnder) ? 2 : (bit == kWallXlo || bit == kWallXhi) ? 0 : 1;
                            const int s = (bit == kWallRoof || bit == kWallXlo || bit == kWallYlo) ? 0 : 1;
                            const amrex::Real fw = ff.f[n][s];
                            if (!(beta(ii, jj, kk) > amrex::Real(0)) && fw < amrex::Real(1) - ftol) {
                                const int step = (s == 0) ? -1 : 1;
                                const int ni = ii + (n == 0 ? step : 0), nj = jj + (n == 1 ? step : 0), nk = kk + (n == 2 ? step : 0);
                                if (ni >= ilo2 && ni <= ihi2 && nj >= jlo2 && nj <= jhi2 && nk >= klo && nk <= khi
                                    && cl(ni, nj, nk, 0) > amrex::Real(0) && cl(ni, nj, nk, 0) < amrex::Real(1)) {
                                    const FaceFractions nf{{{fx(ni, nj, nk), fx(ni + 1, nj, nk)}, {fy(ni, nj, nk), fy(ni, nj + 1, nk)},
                                                            {fz(ni, nj, nk), fz(ni, nj, nk + 1)}}};
                                    const amrex::Real half = amrex::Real(0.5);
                                    const amrex::Real zl = g.zlo(ii, jj, kk), zh = g.zlo(ii, jj, kk + 1);
                                    const amrex::Real nmid[3] = {g.xc(ni), g.yc(nj), half * (zl + zh)};
                                    const amrex::Real nhs[3] = {half * g.dx, half * g.dy, half * (zh - zl)};
                                    const amrex::Real nlo[3] = {g.xc(ni) - half * g.dx, g.yc(nj) - half * g.dy, zl};
                                    const amrex::Real nhi[3] = {g.xc(ni) + half * g.dx, g.yc(nj) + half * g.dy, zh};
                                    // The shared face: the fluid cell's own, on the wall's side.
                                    const amrex::Real plane = (n == 2) ? (s == 0 ? zl : zh)
                                                              : (n == 0) ? g.xc(ii) + step * half * g.dx : g.yc(jj) + step * half * g.dy;
                                    drawn = solid_face_rectangle(nmid, nhs, nlo, nhi, n, plane, cl(ni, nj, nk, 0), static_cast<int>(cl(ni, nj, nk, 1)),
                                                                 nf, fw, c, h);
                                }
                            }
                            if (!drawn) { wall_rectangle(g, ii, jj, kk, bit, mask, beta(ii, jj, kk), ff, below >= f, c, h); }
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

namespace {

//! dist = min(dist, d_b + z0_b - z0) on the cells and their ghost layer, z0 the ground's by the
//! accessor Z0 (UniformZ0 or MapZ0: one kernel per kind). The number of valid cells that fell.
template <typename Z0>
amrex::Long limit_distance (amrex::MultiFab& dist, const amrex::MultiFab& building, const Roughness& rough, amrex::Real z0b)
{
    amrex::ReduceOps<amrex::ReduceOpSum> op;
    amrex::ReduceData<amrex::Long> data(op);
    for (amrex::MFIter mfi(dist, false); mfi.isValid(); ++mfi) {
        const auto d = dist.array(mfi);
        const auto db = building.const_array(mfi);
        const Z0 z0 = Z0::of(rough, mfi);
        const amrex::Box vb = mfi.validbox();
        op.eval(amrex::grow(vb, 1), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
        {
            // No wall found: the ground's distance stays. With equal roughnesses the shift is 0
            // and the cell takes d_b itself.
            if (!(db(i, j, k) < kBuildingFar)) { return {amrex::Long(0)}; }
            const amrex::Real de = db(i, j, k) + (z0b - z0(i, j));
            if (!(de < d(i, j, k))) { return {amrex::Long(0)}; }
            d(i, j, k) = de;
            return {vb.contains(i, j, k) ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    amrex::Long n = amrex::get<0>(data.value(op));
    amrex::ParallelDescriptor::ReduceLongSum(n);
    return n;
}

} // namespace

amrex::Long BuildingTke::limit_wall_distance (amrex::MultiFab& dist, const Roughness& ground) const
{
    const char* op = "BuildingTke::limit_wall_distance";
    if (dist.boxArray() != m_dist.boxArray() || dist.DistributionMap() != m_dist.DistributionMap() || dist.nGrow() < 1) {
        amrex::Abort(std::string("gust: ") + op + ": the wall distance is not on the buildings' cells with a ghost layer");
    }
    if (ground.map == nullptr) {
        if (std::isnan(ground.z0) || !(ground.z0 >= 0 && std::isfinite(ground.z0))) {
            amrex::Abort(std::string("gust: ") + op + ": the ground's z0 must be finite and >= 0");
        }
        return limit_distance<UniformZ0>(dist, m_dist, ground, m_z0b);
    }
    check_roughness_map(op, *ground.map, dist);   // on the ground slab, with a ghost column in x and y
    return limit_distance<MapZ0>(dist, m_dist, ground, m_z0b);
}

void BuildingTke::zero_solid_cells (amrex::MultiFab& field) const
{
    if (field.boxArray() != m_dist.boxArray() || field.DistributionMap() != m_dist.DistributionMap()) {
        amrex::Abort("gust: BuildingTke::zero_solid_cells: the field is not on the buildings' cells");
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
    //! The largest part of a face's area that one of the cell's walls takes (wall_area_fraction).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real blocked (int i, int j, int k) const noexcept
    {
        const int mask = static_cast<int>(wl(i, j, k, kWallMask));
        amrex::Real f = 0;
        if ((mask & kWallRoof) != 0) { f = amrex::max(f, wall_area_fraction(fz(i, j, k), fz(i, j, k + 1), beta(i, j, k))); }
        if ((mask & kWallUnder) != 0) { f = amrex::max(f, wall_area_fraction(fz(i, j, k + 1), fz(i, j, k), beta(i, j, k))); }
        if ((mask & kWallXlo) != 0) { f = amrex::max(f, wall_area_fraction(fx(i, j, k), fx(i + 1, j, k), beta(i, j, k))); }
        if ((mask & kWallXhi) != 0) { f = amrex::max(f, wall_area_fraction(fx(i + 1, j, k), fx(i, j, k), beta(i, j, k))); }
        if ((mask & kWallYlo) != 0) { f = amrex::max(f, wall_area_fraction(fy(i, j, k), fy(i, j + 1, k), beta(i, j, k))); }
        if ((mask & kWallYhi) != 0) { f = amrex::max(f, wall_area_fraction(fy(i, j + 1, k), fy(i, j, k), beta(i, j, k))); }
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
