#include "CanopyDrag.H"

#include "CVFluxes.H"
#include "ImplicitOperator.H"
#include "Stencil_K.H"

#include <AMReX_GpuQualifiers.H>
#include <AMReX_MFIter.H>

#include <cmath>
#include <string>

namespace gust {

namespace {

//! The face's own component and |U|^2 at the c-face (i, j, k): the other two
//! components are the four-point means of their faces (Stencil_K.H).
template <int C>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void face_velocity (amrex::Array4<amrex::Real const> const& u, amrex::Array4<amrex::Real const> const& v,
                    amrex::Array4<amrex::Real const> const& w, int i, int j, int k, amrex::Real& own,
                    amrex::Real& speed2) noexcept
{
    using namespace stencil;
    constexpr int T = face(C);
    const amrex::Real uu = (C == 0) ? u(i, j, k) : interp<XF, T>(u, i, j, k);
    const amrex::Real vv = (C == 1) ? v(i, j, k) : interp<YF, T>(v, i, j, k);
    const amrex::Real ww = (C == 2) ? w(i, j, k) : interp<ZF, T>(w, i, j, k);
    own = (C == 0) ? uu : (C == 1) ? vv : ww;
    speed2 = uu * uu + vv * vv + ww * ww;
}

template <int C>
void add_sink (const amrex::MultiFab& fc, const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w,
               amrex::MultiFab& R)
{
    for (amrex::MFIter mfi(R, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto r = R.array(mfi);
        const auto f = fc.const_array(mfi);
        const auto ua = u.const_array(mfi), va = v.const_array(mfi), wa = w.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real own, s2;
            face_velocity<C>(ua, va, wa, i, j, k, own, s2);
            r(i, j, k) += f(i, j, k) * std::sqrt(s2) * own;
        });
    }
}

//! d(f |U| u_c)/d u_c = f (|U| + u_c^2 / |U|); 0 at rest.
template <int C>
void add_derivative (const amrex::MultiFab& fc, const amrex::MultiFab& u, const amrex::MultiFab& v,
                     const amrex::MultiFab& w, amrex::MultiFab& A, int comp)
{
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto a = A.array(mfi);
        const auto f = fc.const_array(mfi);
        const auto ua = u.const_array(mfi), va = v.const_array(mfi), wa = w.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real own, s2;
            face_velocity<C>(ua, va, wa, i, j, k, own, s2);
            const amrex::Real s = std::sqrt(s2);
            a(i, j, k, comp) += (s > 0) ? f(i, j, k) * (s + own * own / s) : amrex::Real(0);
        });
    }
}

//! f |U| with |U| frozen.
template <int C>
void add_picard (const amrex::MultiFab& fc, const amrex::MultiFab& u, const amrex::MultiFab& v,
                 const amrex::MultiFab& w, amrex::MultiFab& A, int comp)
{
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto a = A.array(mfi);
        const auto f = fc.const_array(mfi);
        const auto ua = u.const_array(mfi), va = v.const_array(mfi), wa = w.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            amrex::Real own, s2;
            face_velocity<C>(ua, va, wa, i, j, k, own, s2);
            a(i, j, k, comp) += f(i, j, k) * std::sqrt(s2);
        });
    }
}

template <int C>
void add_diagonal (CanopyDiagonal d, const amrex::MultiFab& fc, const amrex::MultiFab& u, const amrex::MultiFab& v,
                   const amrex::MultiFab& w, amrex::MultiFab& A, int comp)
{
    switch (d) {
    case CanopyDiagonal::derivative: add_derivative<C>(fc, u, v, w, A, comp); break;
    case CanopyDiagonal::picard: add_picard<C>(fc, u, v, w, A, comp); break;
    case CanopyDiagonal::none: break;
    }
}

//! The fields' layouts against the drag's (out on the c-faces).
void check_layout (const char* op, const CanopyFaceDrag& drag, int c, const amrex::MultiFab& u, const amrex::MultiFab& v,
                   const amrex::MultiFab& w, const amrex::MultiFab& out)
{
    const amrex::MultiFab* vel[3] = {&u, &v, &w};
    for (int d = 0; d < 3; ++d) {
        if (vel[d]->boxArray() != drag.f(d).boxArray() || vel[d]->DistributionMap() != drag.f(d).DistributionMap()
            || vel[d]->nGrow() < 1) {
            amrex::Abort(std::string("gust: ") + op + ": the velocity must be on the canopy's layout with a ghost cell");
        }
    }
    if (out.boxArray() != drag.f(c).boxArray() || out.DistributionMap() != drag.f(c).DistributionMap()) {
        amrex::Abort(std::string("gust: ") + op + ": the output must be on the canopy's face layout");
    }
}

} // namespace

std::string canopy_diagonal_name (CanopyDiagonal d)
{
    switch (d) {
    case CanopyDiagonal::derivative: return "derivative";
    case CanopyDiagonal::picard: return "picard";
    case CanopyDiagonal::none: return "none";
    }
    return "";
}

void CanopyFaceDrag::build (const Canopy& canopy, const amrex::Geometry& geom)
{
    if (!canopy.active()) { amrex::Abort("gust: CanopyFaceDrag::build needs an active canopy"); }
    m_diagonal = canopy.drag_diagonal();
    const amrex::MultiFab& cell = canopy.drag();
    const int klo = geom.Domain().smallEnd(2);
    const int khi = geom.Domain().bigEnd(2);
    for (int b = 0; b < static_cast<int>(cell.boxArray().size()); ++b) {
        const amrex::Box& bx = cell.boxArray()[b];
        if (bx.smallEnd(2) != klo || bx.bigEnd(2) != khi) {
            amrex::Abort("gust: CanopyFaceDrag::build: every box must span the whole z column");
        }
    }
    for (int c = 0; c < AMREX_SPACEDIM; ++c) {
        amrex::IntVect t(0);
        t[c] = 1;
        auto& fc = m_f[static_cast<std::size_t>(c)];
        fc.define(amrex::convert(cell.boxArray(), t), cell.DistributionMap(), 1, 0);
        for (amrex::MFIter mfi(fc, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto o = fc.array(mfi);
            const auto f = cell.const_array(mfi);
            if (c == 2) {
                // The boxes span the column (Canopy.H): the ground and top faces take the cell inside.
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const int kb = (k > klo) ? k - 1 : k;
                    const int kt = (k <= khi) ? k : k - 1;
                    o(i, j, k) = amrex::Real(0.5) * (f(i, j, kb) + f(i, j, kt));
                });
            } else {
                // The ghost column on the low side of the box (Canopy.H) gives the left cell.
                const int di = (c == 0) ? 1 : 0, dj = (c == 1) ? 1 : 0;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    o(i, j, k) = amrex::Real(0.5) * (f(i - di, j - dj, k) + f(i, j, k));
                });
            }
        }
    }
}

void add_canopy_drag (const CanopyFaceDrag& drag, const amrex::MultiFab& u, const amrex::MultiFab& v,
                      const amrex::MultiFab& w, amrex::MultiFab& Ru, amrex::MultiFab& Rv, amrex::MultiFab& Rw)
{
    if (drag.empty()) { return; }
    check_layout("add_canopy_drag", drag, 0, u, v, w, Ru);
    check_layout("add_canopy_drag", drag, 1, u, v, w, Rv);
    check_layout("add_canopy_drag", drag, 2, u, v, w, Rw);
    add_sink<0>(drag.f(0), u, v, w, Ru);
    add_sink<1>(drag.f(1), u, v, w, Rv);
    add_sink<2>(drag.f(2), u, v, w, Rw);
}

void add_canopy_drag_to_implicit (const CanopyFaceDrag& drag, int c, const amrex::MultiFab& u,
                                  const amrex::MultiFab& v, const amrex::MultiFab& w, amrex::MultiFab& A)
{
    if (drag.empty()) { return; }
    if (c < 0 || c > 2) { amrex::Abort("gust: add_canopy_drag_to_implicit: c must be 0, 1 or 2"); }
    check_layout("add_canopy_drag_to_implicit", drag, c, u, v, w, A);
    switch (c) {
    case 0: add_diagonal<0>(drag.diagonal(), drag.f(0), u, v, w, A, coef_P); break;
    case 1: add_diagonal<1>(drag.diagonal(), drag.f(1), u, v, w, A, coef_P); break;
    default: add_diagonal<2>(drag.diagonal(), drag.f(2), u, v, w, A, coef_P); break;
    }
}

void canopy_rows (const CanopyFaceDrag& drag, const amrex::MultiFab& u, const amrex::MultiFab& v,
                  const amrex::MultiFab& w, CanopyRows& rows)
{
    if (drag.empty()) { amrex::Abort("gust: canopy_rows: no canopy"); }
    for (int c = 0; c < 3; ++c) {
        auto& a = rows.ap[static_cast<std::size_t>(c)];
        if (a.empty()) { a.define(drag.f(c).boxArray(), drag.f(c).DistributionMap(), 1, 0); }
        check_layout("canopy_rows", drag, c, u, v, w, a);
        a.setVal(0.0);
    }
    add_diagonal<0>(drag.diagonal(), drag.f(0), u, v, w, rows.ap[0], 0);
    add_diagonal<1>(drag.diagonal(), drag.f(1), u, v, w, rows.ap[1], 0);
    add_diagonal<2>(drag.diagonal(), drag.f(2), u, v, w, rows.ap[2], 0);
}

} // namespace gust
