#ifndef GUST_BUILDINGS_K_H_
#define GUST_BUILDINGS_K_H_

// The kernels of the building geometry (Buildings.H). The height map is a
// grid of rectangular tiles of constant height, in the column heights above
// the ground, so the mean over a cell's footprint is an exact sum over the
// tiles it overlaps.

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>

namespace gust {

//! The tile i of direction d runs from xe[i] to xe[i + 1]; the height of tile
//! (i, j) is h[i * ny + j], ERF's order, solid from zs[i * ny + j] up to it.
//! lx, ly: the period (0: not periodic).
struct HeightMapView
{
    const double* xe;
    const double* ye;
    const double* h;
    const double* zs;
    int nx, ny;
    double lx, ly;
};

//! The largest i in [0, n - 1] with e[i] <= x (x inside [e[0], e[n]]).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int tile_index (const double* e, int n, double x) noexcept
{
    int lo = 0, hi = n - 1;
    while (lo < hi) {
        const int mid = (lo + hi + 1) / 2;
        if (e[mid] <= x) { lo = mid; } else { hi = mid - 1; }
    }
    return lo;
}

//! The tile holding x or a periodic image of x; -1 when none does.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int locate_tile (const double* e, int n, double period, double x) noexcept
{
    const int m = period > 0 ? 1 : 0;
    for (int s = -m; s <= m; ++s) {
        const double xs = x + s * period;
        if (xs >= e[0] && xs <= e[n]) { return tile_index(e, n, xs); }
    }
    return -1;
}

//! The mean of f(height, start) over the rectangle [xa, xb] x [ya, yb], periodic
//! images included. Outside the tiles the height is 0, where f is taken as 0.
template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double area_mean (const HeightMapView& m, double xa, double xb, double ya, double yb, const F& f) noexcept
{
    double sum = 0;
    const int mx = m.lx > 0 ? 1 : 0, my = m.ly > 0 ? 1 : 0;
    for (int sx = -mx; sx <= mx; ++sx) {
        const double ax = xa + sx * m.lx, bx = xb + sx * m.lx;
        if (bx <= m.xe[0] || ax >= m.xe[m.nx]) { continue; }
        for (int i = tile_index(m.xe, m.nx, amrex::max(ax, m.xe[0])); i < m.nx && m.xe[i] < bx; ++i) {
            const double wx = amrex::min(bx, m.xe[i + 1]) - amrex::max(ax, m.xe[i]);
            for (int sy = -my; sy <= my; ++sy) {
                const double ay = ya + sy * m.ly, by = yb + sy * m.ly;
                if (by <= m.ye[0] || ay >= m.ye[m.ny]) { continue; }
                for (int j = tile_index(m.ye, m.ny, amrex::max(ay, m.ye[0])); j < m.ny && m.ye[j] < by; ++j) {
                    const double wy = amrex::min(by, m.ye[j + 1]) - amrex::max(ay, m.ye[j]);
                    sum += wx * wy * f(m.h[i * m.ny + j], m.zs[i * m.ny + j]);
                }
            }
        }
    }
    return sum / ((xb - xa) * (yb - ya));
}

//! The mean of f(height, start) along the part [ta, tb] of a vertical face at
//! xn. A face on the edge between two tiles takes the larger f: the side of a
//! building blocks the face it lies on (eps: the rounding of xn). ne, nn, ln:
//! tile edges, count and period normal to the face; tan_e, nt, lt: along it;
//! sn, st: the strides of the two directions in h.
template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double line_mean (const double* ne, int nn, double ln, const double* tan_e, int nt, double lt, const double* h, const double* zs, int sn, int st,
                  double xn, double ta, double tb, double eps, const F& f) noexcept
{
    const int na = locate_tile(ne, nn, ln, xn - eps);
    const int nb = locate_tile(ne, nn, ln, xn + eps);
    if (na < 0 && nb < 0) { return 0; }
    double sum = 0;
    const int mt = lt > 0 ? 1 : 0;
    for (int s = -mt; s <= mt; ++s) {
        const double a = ta + s * lt, b = tb + s * lt;
        if (b <= tan_e[0] || a >= tan_e[nt]) { continue; }
        for (int j = tile_index(tan_e, nt, amrex::max(a, tan_e[0])); j < nt && tan_e[j] < b; ++j) {
            const double w = amrex::min(b, tan_e[j + 1]) - amrex::max(a, tan_e[j]);
            const double fa = na >= 0 ? f(h[na * sn + j * st], zs[na * sn + j * st]) : 0.0;
            const double fb = nb >= 0 ? f(h[nb * sn + j * st], zs[nb * sn + j * st]) : 0.0;
            sum += w * amrex::max(fa, fb);
        }
    }
    return sum / (tb - ta);
}

//! The part of the heights [zb, zt] above the ground inside a tile's solid [zs, h].
struct VerticalOverlap
{
    double zb, zt;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() (double h, double zs) const noexcept
    {
        return amrex::max(0.0, amrex::min(h, zt) - amrex::max(zs, zb)) / (zt - zb);
    }
};

//! 1 where a horizontal face at height zf above the ground lies in a building or
//! on its roof or underside (tol: the rounding of zf): those block the face.
struct HorizontalBlock
{
    double zf, tol;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double operator() (double h, double zs) const noexcept { return (h > 0 && h >= zf - tol && zs <= zf + tol) ? 1.0 : 0.0; }
};

//! buildings.small_volfrac, at both ends: below small a fraction is 0, within
//! small of one it is 1. After it, a fluid cell has beta == 0 and a solid one beta == 1.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double clip_fraction (double beta, double small) noexcept
{
    if (beta >= 1 - small) { return 1; }
    return beta < small ? 0 : beta;
}

//! A wall on the low side of a cell along a normal, from the blocked fractions
//! of its low and high faces: the difference of two opposite faces is the wall
//! area in the cell facing that way. In a fluid cell any blocked face is a wall, so
//! a one-cell street has both. tol: the part of a face below which there is none.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool wall_on_low_side (double f_lo, double f_hi, double beta, double tol) noexcept
{
    return f_lo - f_hi > tol || (beta == 0 && f_lo > tol);
}

//! The fluid-centroid distance of a wall cell from its wall, along the normal.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double wall_d1 (double beta, double dn) noexcept { return 0.5 * (1 - beta) * dn; }

//! The distance of the reference cell (the next one out, its size dn_ref, its
//! fraction beta_ref) from the wall; 0 when there is none to use: outside the
//! domain (ok false), or solid. Equal layers give ERF's (3/2 - beta) dn.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double wall_d2 (double beta, double dn, double dn_ref, double beta_ref, bool ok) noexcept
{
    return (ok && beta_ref < 1) ? (1 - beta) * dn + 0.5 * dn_ref : 0.0;
}

} // namespace gust

#endif
