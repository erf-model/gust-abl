#ifndef GUST_TRANSFER_K_H_
#define GUST_TRANSFER_K_H_

// The horizontal bilinear interpolation per zeta-level shared by the grid
// hierarchy's prolongation (Transfer.H, prolong_fields) and the nested
// patch's parent reader (NestParent.H). The fine grid is r times finer than
// the coarse grid in x and y and has the same z-levels, so nothing is
// interpolated in the vertical.

#include "Stencil_K.H"

#include <AMReX_Box.H>
#include <AMReX_MultiFab.H>

namespace gust::transfer {

//! Where a fine point sits among the coarse points along one direction:
//! c = clamp(floor(m / den), imin, imax), t = (m - den c) / den and the
//! coarse index i_c = c + clo, with m = 2 (i - flo) + 1 - r, den = 2r for a
//! cell-centred direction (the fine centre at (i_c + 1/2 + t) coarse cells)
//! and m = i - flo, den = r with floor((m - 1) / den) for a nodal one. clo is
//! the coarse index of the fine grid's first cell (the coarse domain's low
//! end plus an offset: 0 for the grid hierarchy, the patch's position in the
//! parent for a nest). The clamp bounds give the linear extrapolation at a
//! non-periodic face of the coarse domain (t < 0 or t > 1); a periodic
//! direction is not clamped (bounds a period beyond the domain) and reads
//! the periodic images, which the caller fills as far as its points reach
//! (the grid hierarchy's valid points reach i_c = -1 and N, a nest's ghost
//! layers one coarse cell further). lin = 0 and
//! step = 0 for a non-periodic cell-centred direction one coarse cell wide
//! (constant there: no second point to extrapolate with).
struct Axis
{
    int flo{0}, clo{0};
    int r{2};
    int imin{0}, imax{0};
    int step{1};
    amrex::Real lin{1};
    amrex::Real inv{0};   // 1/den, a power of two
};

//! The axis of direction d: fine_lo is the fine domain's low index, cdomain
//! the coarse domain (cperiodic: periodic in d), r the ratio, nodal whether
//! the field is nodal in d, offset the coarse cells between the coarse
//! domain's low end and the fine grid's first cell (>= 0).
inline Axis make_axis (int fine_lo, const amrex::Box& cdomain, bool cperiodic, int d, int r, bool nodal, int offset)
{
    Axis a;
    a.flo = fine_lo;
    a.clo = cdomain.smallEnd(d) + offset;
    a.r = r;
    const int N = cdomain.length(d);
    if (cperiodic) {
        // No clamp: points within a period never reach these bounds.
        a.imin = -N;
        a.imax = 2 * N - 1;
        a.inv = amrex::Real(1) / amrex::Real(nodal ? r : 2 * r);
    } else if (nodal) {
        a.imin = 0;
        a.imax = N - 1;
        a.inv = amrex::Real(1) / amrex::Real(r);
    } else {
        a.inv = amrex::Real(1) / amrex::Real(2 * r);
        if (N >= 2) {
            a.imin = 0;
            a.imax = N - 2;
        } else {
            a.imin = 0;
            a.imax = 0;
            a.step = 0;
            a.lin = 0;
        }
    }
    // The bounds are coarse indices relative to the domain's low end; c is
    // relative to clo.
    a.imin -= offset;
    a.imax -= offset;
    return a;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int floor_div (int a, int b) noexcept
{
    return (a >= 0) ? a / b : -((-a + b - 1) / b);
}

template <bool Nodal>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void locate (const Axis& a, int i, int& ic, amrex::Real& t) noexcept
{
    const int den = Nodal ? a.r : 2 * a.r;
    const int m = Nodal ? (i - a.flo) : (2 * (i - a.flo) + 1 - a.r);
    int c = floor_div(Nodal ? m - 1 : m, den);
    c = amrex::min(amrex::max(c, a.imin), a.imax);
    t = a.lin * amrex::Real(m - den * c) * a.inv;
    ic = c + a.clo;
}

//! fine = the bilinear interpolation of the coarse image c (on fine's
//! BoxArray coarsened, fine's DistributionMapping, the coarse points around
//! every fine point filled) at fine's valid points of type T (stencil bits:
//! nodal in x, y), and at its ghost points within ng of them (0: valid
//! points only, as the grid hierarchy's prolongation; ng z must be 0).
template <int T>
void interpolate (const amrex::MultiFab& c, amrex::MultiFab& fine, const Axis& ax, const Axis& ay,
                  const amrex::IntVect& ng = amrex::IntVect(0))
{
    constexpr bool nx = stencil::nodal(T, 0);
    constexpr bool ny = stencil::nodal(T, 1);
    for (amrex::MFIter mfi(fine, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(ng);
        const auto f = fine.array(mfi);
        const auto q = c.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            int I = 0, J = 0;
            amrex::Real tx = 0, ty = 0;
            locate<nx>(ax, i, I, tx);
            locate<ny>(ay, j, J, ty);
            const int I1 = I + ax.step, J1 = J + ay.step;
            const amrex::Real one(1);
            f(i, j, k) = (one - ty) * ((one - tx) * q(I, J, k) + tx * q(I1, J, k))
                       + ty * ((one - tx) * q(I, J1, k) + tx * q(I1, J1, k));
        });
    }
}

} // namespace gust::transfer

#endif
