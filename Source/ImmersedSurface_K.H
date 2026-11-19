#ifndef GUST_IMMERSED_SURFACE_K_H_
#define GUST_IMMERSED_SURFACE_K_H_

// The kernels of an immersed ground that is a continuous surface (ImmersedSurface.H). The ground
// height h is given at the mesh's nodes. Over a cell's footprint it is piecewise linear on the four
// triangles that the footprint's centre (the mean of the four corner heights) makes with each pair
// of neighbouring corners; along a cell edge it is the line between the two nodes, so two cells
// that share a face see the same ground on it. Everything a cell or a face needs is a mean of
//     (t - h)_+^n,  n = 0, 1, 2,  t a level,
// over a triangle or a segment, which has a closed form (docs/immersed_terrain.rst, Geometry):
//     n = 0  the part of the footprint whose ground lies below t;
//     n = 1  the depth of fluid below t per unit footprint;
//     n = 2  twice the first moment, about the ground, of that fluid.
// Double precision in every build: the fractions are differences of these means.
// Every division is by a denominator made positive first (safe_den): a compiler may evaluate a
// branch that is not taken, and 0 / 0 there would raise the invalid flag (amrex.fpe_trap_invalid).

#include <AMReX_Algorithm.H>
#include <AMReX_GpuQualifiers.H>

namespace gust {

//! d where it is positive, else 1 (the quotient is then not used).
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double safe_den (double d) noexcept { return d > 0 ? d : 1.0; }

//! The means of (t - h)_+^n, n = 0, 1, 2.
struct SurfaceMoments
{
    double m0{0}, m1{0}, m2{0};
};

//! Over a triangle whose ground is linear with the corner heights v1, v2, v3.
//! With a <= b <= c the sorted heights, m their mean and
//! var = ((a - b)^2 + (b - c)^2 + (c - a)^2) / 36 the variance of h:
//!   t <= a      0, 0, 0;
//!   a < t <= b  the corner at a below t is a similar triangle of area p = (t - a)^2 / ((b - a)(c - a))
//!               on which t - h is linear with corner values (t - a, 0, 0): p, p (t - a) / 3, p (t - a)^2 / 6;
//!   b < t < c   the whole triangle less the corner at c above t, p = (c - t)^2 / ((c - b)(c - a)):
//!               1 - p, (t - m) + p (c - t) / 3, (t - m)^2 + var - p (c - t)^2 / 6;
//!   t >= c      1, t - m, (t - m)^2 + var.
//! A flat triangle (a = b = c) takes the first or the last branch. In the third branch m1 loses
//! digits as t -> b with a = b (an absolute error of about one ulp of the height, below the clip of
//! the fractions). m is formed from the differences, so a flat triangle's mean is its height exactly.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SurfaceMoments triangle_moments (double v1, double v2, double v3, double t) noexcept
{
    double a = v1, b = v2, c = v3;
    if (a > b) { amrex::Swap(a, b); }
    if (b > c) { amrex::Swap(b, c); }
    if (a > b) { amrex::Swap(a, b); }
    if (t <= a) { return {}; }
    if (t <= b) {
        const double x = t - a;
        const double p = x * x / safe_den((b - a) * (c - a));
        return {p, p * x / 3, p * x * x / 6};
    }
    const double m = a + ((b - a) + (c - a)) / 3;
    const double var = ((a - b) * (a - b) + (b - c) * (b - c) + (c - a) * (c - a)) / 36;
    const double full = (t - m) * (t - m) + var;
    if (t >= c) { return {1.0, t - m, full}; }
    const double y = c - t;
    const double p = y * y / safe_den((c - b) * (c - a));
    return {1 - p, (t - m) + p * y / 3, full - p * y * y / 6};
}

//! Over a cell's footprint with the corner heights h00, h10, h11, h01 (in order around it): the
//! mean of its four triangles, which have equal areas.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
SurfaceMoments footprint_moments (double h00, double h10, double h11, double h01, double t) noexcept
{
    const double hc = 0.25 * (h00 + h10 + h11 + h01);
    const SurfaceMoments s0 = triangle_moments(h00, h10, hc, t), s1 = triangle_moments(h10, h11, hc, t);
    const SurfaceMoments s2 = triangle_moments(h11, h01, hc, t), s3 = triangle_moments(h01, h00, hc, t);
    return {0.25 * (s0.m0 + s1.m0 + s2.m0 + s3.m0), 0.25 * (s0.m1 + s1.m1 + s2.m1 + s3.m1),
            0.25 * (s0.m2 + s1.m2 + s2.m2 + s3.m2)};
}

//! The mean of (t - h)_+ along a segment whose ground runs linearly from v1 to v2: with a <= b,
//! 0 for t <= a, (t - a)^2 / (2 (b - a)) for a < t < b, t - (a + b) / 2 for t >= b.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double segment_depth (double v1, double v2, double t) noexcept
{
    const double a = amrex::min(v1, v2), b = amrex::max(v1, v2);
    if (t <= a) { return 0; }
    if (t >= b) { return t - 0.5 * (a + b); }
    return (t - a) * (t - a) / safe_den(2 * (b - a));
}

//! The solid fraction of the layer [zb, zt] from the depths of fluid below its two levels,
//! 1 - (m1(zt) - m1(zb)) / (zt - zb), in [0, 1].
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double solid_fraction (double depth_b, double depth_t, double zb, double zt) noexcept
{
    return amrex::min(1.0, amrex::max(0.0, 1 - (depth_t - depth_b) / (zt - zb)));
}

//! The mean height above the ground of the fluid of a cell's layer, the fluid's first moment about
//! the ground over its volume: (m2(zt) - m2(zb)) / (2 (m1(zt) - m1(zb))); 0 in a cell with no fluid.
//! Above the ground it is the cell centre's height above the mean ground of the footprint.
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double fluid_height (const SurfaceMoments& b, const SurfaceMoments& t) noexcept
{
    const double depth = t.m1 - b.m1;
    return depth > 0 ? (t.m2 - b.m2) / safe_den(2 * depth) : 0.0;
}

} // namespace gust

#endif
