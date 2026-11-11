#include "BuildingForcing.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>
#include <AMReX_iMultiFab.H>

#include <cmath>
#include <limits>

namespace gust {

namespace {

//! The row the boundary conditions own among component c's: w on the ground, where they
//! set w from the terrain's slope. A fixed w = 0 there would fight them on a slope.
int ground_row (int c, const amrex::Geometry& geom)
{
    return c == 2 ? geom.Domain().smallEnd(2) : std::numeric_limits<int>::min();
}

} // namespace

void BuildingForcing::build (const Buildings& b, const BuildingsOptions& opt, const amrex::Geometry& geom, const amrex::MultiFab& detJ)
{
    if (!b.active()) { amrex::Abort("gust: BuildingForcing::build needs the buildings' geometry"); }
    m_b = &b;
    m_geom = geom;
    m_detJ = &detJ;
    m_z0 = static_cast<amrex::Real>(opt.z0);
    m_wall_form = opt.wall_form;
    m_rate.define(b.beta().boxArray(), b.beta().DistributionMap(), 4, 1);
    m_rate.setVal(0.0);
    // A side's wall cell has its centroid dx / 2 or dy / 2 from the side; below 10 z0 the law
    // reads the speed at 10 z0, so the stress no longer follows the cell size.
    if (side_law() && amrex::min(geom.CellSize(0), geom.CellSize(1)) < 20 * m_z0) {
        amrex::Print() << "gust: note: buildings.z0 = " << opt.z0 << " m: the sides' wall cells are within 10 z0 of their sides ("
                       << "cells of " << geom.CellSize(0) << " by " << geom.CellSize(1) << " m), where the wall law takes the "
                       << "speed as at 10 z0\n";
    }
}

void BuildingForcing::update (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w)
{
    if (!wall_law()) { return; }
    const amrex::Real z0 = m_z0, kappa = amrex::Real(0.41), umin = kBuildingMinSpeed;
    const amrex::Real dzeta = m_geom.CellSize(2);
    const amrex::MultiFab& bz = m_b->beta_face(2);
    m_rate.setVal(0.0);
    // Roofs and undersides: on u and v.
    for (amrex::MFIter mfi(m_rate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto r = m_rate.array(mfi);
        const auto uu = u.const_array(mfi);
        const auto vv = v.const_array(mfi);
        const auto beta = m_b->beta().const_array(mfi);
        const auto wl = m_b->wall().const_array(mfi);
        const auto fz = bz.const_array(mfi);
        const auto J = m_detJ->const_array(mfi);
        const int klo = m_geom.Domain().smallEnd(2), khi = m_geom.Domain().bigEnd(2);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int mask = static_cast<int>(wl(i, j, k, kWallMask));
            if (beta(i, j, k) >= amrex::Real(1) || (mask & (kWallRoof | kWallUnder)) == 0) { return; }
            auto speed = [&] (int kk) {
                const amrex::Real uc = amrex::Real(0.5) * (uu(i, j, kk) + uu(i + 1, j, kk));
                const amrex::Real vc = amrex::Real(0.5) * (vv(i, j, kk) + vv(i, j + 1, kk));
                return std::sqrt(uc * uc + vc * vc);
            };
            const amrex::Real own = speed(k);
            const amrex::Real inv = amrex::Real(1) / (amrex::max(own, umin) * J(i, j, k) * dzeta);
            const amrex::Real d1 = wl(i, j, k, kWallD1z), d2 = wl(i, j, k, kWallD2z);
            amrex::Real rate = 0, top = 0;
            if ((mask & kWallRoof) != 0) {
                const amrex::Real ref = k < khi ? speed(k + 1) : amrex::Real(0);
                const amrex::Real us2 = form7_ustar2(own, ref, beta(i, j, k), d1, d2, z0, kappa);
                rate += fz(i, j, k) * us2 * inv;
                top = amrex::max(top, us2);
            }
            if ((mask & kWallUnder) != 0) {
                const amrex::Real ref = k > klo ? speed(k - 1) : amrex::Real(0);
                const amrex::Real us2 = form7_ustar2(own, ref, beta(i, j, k), d1, d2, z0, kappa);
                rate += fz(i, j, k + 1) * us2 * inv;
                top = amrex::max(top, us2);
            }
            r(i, j, k, 0) = rate;
            r(i, j, k, 1) = rate;
            r(i, j, k, kBuildingUstar2) = top;
        });
    }
    if (side_law()) {
        // Sides normal to x: on v and w. The reference cell is the next one out along x. Only a
        // partial cell weighs it (beta > 0), and such a cell has one wall along x, facing away
        // from a cell that is not solid (Buildings_K.H: wall_on_low_side).
        const amrex::Real dx = m_geom.CellSize(0), dy = m_geom.CellSize(1);
        for (amrex::MFIter mfi(m_rate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto r = m_rate.array(mfi);
            const auto vv = v.const_array(mfi);
            const auto ww = w.const_array(mfi);
            const auto beta = m_b->beta().const_array(mfi);
            const auto wl = m_b->wall().const_array(mfi);
            const auto fx = m_b->beta_face(0).const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int mask = static_cast<int>(wl(i, j, k, kWallMask));
                if (beta(i, j, k) >= amrex::Real(1) || (mask & (kWallXlo | kWallXhi)) == 0) { return; }
                auto speed = [&] (int ii) {
                    const amrex::Real vc = amrex::Real(0.5) * (vv(ii, j, k) + vv(ii, j + 1, k));
                    const amrex::Real wc = amrex::Real(0.5) * (ww(ii, j, k) + ww(ii, j, k + 1));
                    return std::sqrt(vc * vc + wc * wc);
                };
                const amrex::Real own = speed(i);
                const amrex::Real inv = amrex::Real(1) / (amrex::max(own, umin) * dx);
                const amrex::Real d1 = wl(i, j, k, kWallD1x), d2 = wl(i, j, k, kWallD2x);
                amrex::Real rate = 0, top = r(i, j, k, kBuildingUstar2);
                if ((mask & kWallXlo) != 0) {
                    const amrex::Real us2 = form7_ustar2(own, speed(i + 1), beta(i, j, k), d1, d2, z0, kappa);
                    rate += fx(i, j, k) * us2 * inv;
                    top = amrex::max(top, us2);
                }
                if ((mask & kWallXhi) != 0) {
                    const amrex::Real us2 = form7_ustar2(own, speed(i - 1), beta(i, j, k), d1, d2, z0, kappa);
                    rate += fx(i + 1, j, k) * us2 * inv;
                    top = amrex::max(top, us2);
                }
                r(i, j, k, 1) += rate;
                r(i, j, k, 2) += rate;
                r(i, j, k, kBuildingUstar2) = top;
            });
        }
        // Sides normal to y: on u and w.
        for (amrex::MFIter mfi(m_rate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto r = m_rate.array(mfi);
            const auto uu = u.const_array(mfi);
            const auto ww = w.const_array(mfi);
            const auto beta = m_b->beta().const_array(mfi);
            const auto wl = m_b->wall().const_array(mfi);
            const auto fy = m_b->beta_face(1).const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int mask = static_cast<int>(wl(i, j, k, kWallMask));
                if (beta(i, j, k) >= amrex::Real(1) || (mask & (kWallYlo | kWallYhi)) == 0) { return; }
                auto speed = [&] (int jj) {
                    const amrex::Real uc = amrex::Real(0.5) * (uu(i, jj, k) + uu(i + 1, jj, k));
                    const amrex::Real wc = amrex::Real(0.5) * (ww(i, jj, k) + ww(i, jj, k + 1));
                    return std::sqrt(uc * uc + wc * wc);
                };
                const amrex::Real own = speed(j);
                const amrex::Real inv = amrex::Real(1) / (amrex::max(own, umin) * dy);
                const amrex::Real d1 = wl(i, j, k, kWallD1y), d2 = wl(i, j, k, kWallD2y);
                amrex::Real rate = 0, top = r(i, j, k, kBuildingUstar2);
                if ((mask & kWallYlo) != 0) {
                    const amrex::Real us2 = form7_ustar2(own, speed(j + 1), beta(i, j, k), d1, d2, z0, kappa);
                    rate += fy(i, j, k) * us2 * inv;
                    top = amrex::max(top, us2);
                }
                if ((mask & kWallYhi) != 0) {
                    const amrex::Real us2 = form7_ustar2(own, speed(j - 1), beta(i, j, k), d1, d2, z0, kappa);
                    rate += fy(i, j + 1, k) * us2 * inv;
                    top = amrex::max(top, us2);
                }
                r(i, j, k, 0) += rate;
                r(i, j, k, 2) += rate;
                r(i, j, k, kBuildingUstar2) = top;
            });
        }
    }
    m_rate.FillBoundary(m_geom.periodicity());
}

void BuildingForcing::add_sink (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w, amrex::MultiFab& Ru,
                                amrex::MultiFab& Rv, amrex::MultiFab& Rw) const
{
    if (!wall_law()) { return; }
    const amrex::MultiFab* q[3] = {&u, &v, &w};
    amrex::MultiFab* R[3] = {&Ru, &Rv, &Rw};
    // w has no wall parallel to it without the sides' law.
    const int nc = side_law() ? 3 : 2;
    for (int c = 0; c < nc; ++c) {
        const int di = (c == 0) ? 1 : 0, dj = (c == 1) ? 1 : 0, dk = (c == 2) ? 1 : 0;
        for (amrex::MFIter mfi(*R[c], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto r = m_rate.const_array(mfi, c);
            const auto x = q[c]->const_array(mfi);
            const auto a = R[c]->array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                a(i, j, k) += amrex::Real(0.5) * (r(i - di, j - dj, k - dk) + r(i, j, k)) * x(i, j, k);
            });
        }
    }
}

void BuildingForcing::close_viscous_fluxes (const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open,
                                            amrex::Array<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>, AMREX_SPACEDIM>& G) const
{
    if (!wall_law()) { return; }
    // Walls normal to d: z alone without the sides' law.
    for (int d = side_law() ? 0 : 2; d < 3; ++d) {
        for (int c = 0; c < 3; ++c) {
            if (c == d) { continue; }
            // The flux's point lies between the d-faces of the two cells of c's row.
            const int di = (c == 0) ? 1 : 0, dj = (c == 1) ? 1 : 0, dk = (c == 2) ? 1 : 0;
            amrex::MultiFab& Gcd = G[static_cast<std::size_t>(c)][static_cast<std::size_t>(d)];
            for (amrex::MFIter mfi(Gcd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto g = Gcd.array(mfi);
                const auto o = open[static_cast<std::size_t>(d)]->const_array(mfi);
                amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    g(i, j, k) *= amrex::Real(0.5) * (o(i - di, j - dj, k - dk) + o(i, j, k));
                });
            }
        }
    }
}

void BuildingForcing::apply_rows (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w, amrex::MultiFab& Ru,
                                  amrex::MultiFab& Rv, amrex::MultiFab& Rw) const
{
    const amrex::MultiFab* q[3] = {&u, &v, &w};
    amrex::MultiFab* R[3] = {&Ru, &Rv, &Rw};
    for (int c = 0; c < 3; ++c) {
        for (amrex::MFIter mfi(*R[c], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = R[c]->array(mfi);
            const auto x = q[c]->const_array(mfi);
            const auto b = m_b->beta_face(c).const_array(mfi);
            const int kg = ground_row(c, m_geom);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (b(i, j, k) == amrex::Real(1) && k != kg) { a(i, j, k) = x(i, j, k); }
            });
        }
    }
}

void BuildingForcing::add_to_implicit (int c, const amrex::Array<const amrex::MultiFab*, AMREX_SPACEDIM>& open, amrex::MultiFab& A) const
{
    const int di = (c == 0) ? 1 : 0, dj = (c == 1) ? 1 : 0, dk = (c == 2) ? 1 : 0;
    if (wall_law() && (c < 2 || side_law())) {
        for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = A.array(mfi);
            const auto r = m_rate.const_array(mfi, c);
            // Not on the ground's w: that row is the boundary conditions' identity.
            const int kg = ground_row(c, m_geom);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (k != kg) { a(i, j, k, coef_P) += amrex::Real(0.5) * (r(i - di, j - dj, k - dk) + r(i, j, k)); }
            });
        }
        // The couplings through the wall faces normal to d, as close_viscous_fluxes. The whole
        // coupling is scaled, its upwind weight too (the open mass flux has already reduced that
        // one): this changes the rate of convergence, not the fixed point. a_P keeps the closed
        // part, so the row stays diagonally dominant.
        for (int d = side_law() ? 0 : 2; d < 3; ++d) {
            if (d == c) { continue; }
            const int ei = (d == 0) ? 1 : 0, ej = (d == 1) ? 1 : 0, ek = (d == 2) ? 1 : 0;
            const int lo = coef_XM + 2 * d, hi = lo + 1;
            for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto a = A.array(mfi);
                const auto o = open[static_cast<std::size_t>(d)]->const_array(mfi);
                amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    a(i, j, k, lo) *= amrex::Real(0.5) * (o(i - di, j - dj, k - dk) + o(i, j, k));
                    a(i, j, k, hi) *= amrex::Real(0.5) * (o(i - di + ei, j - dj + ej, k - dk + ek) + o(i + ei, j + ej, k + ek));
                });
            }
        }
    }
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto a = A.array(mfi);
        const auto b = m_b->beta_face(c).const_array(mfi);
        const int kg = ground_row(c, m_geom);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (b(i, j, k) == amrex::Real(1) && k != kg) {
                for (int n = 0; n < n_implicit_coefs; ++n) { a(i, j, k, n) = 0; }
                a(i, j, k, coef_P) = 1;
            }
        });
    }
}

void BuildingForcing::add_fixed_rows (int c, amrex::MultiFab& mask) const
{
    for (amrex::MFIter mfi(mask, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto m = mask.array(mfi);
        const auto b = m_b->beta_face(c).const_array(mfi);
        const int kg = ground_row(c, m_geom);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            if (b(i, j, k) == amrex::Real(1) && k != kg) { m(i, j, k) = 1; }
        });
    }
}

amrex::Long count_sealed_fluid_cells (const Buildings& b, const amrex::Geometry& geom)
{
    // Each fluid cell takes the largest label among those it reaches through an open
    // face, until nothing changes: one label per connected region.
    const amrex::MultiFab& beta = b.beta();
    amrex::iMultiFab label(beta.boxArray(), beta.DistributionMap(), 1, 1);
    const amrex::Box& dom = geom.Domain();
    const auto lo = amrex::lbound(dom);
    const auto len = amrex::length(dom);
    label.setVal(-1);
    for (amrex::MFIter mfi(label, false); mfi.isValid(); ++mfi) {
        const auto l = label.array(mfi);
        const auto bc = beta.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            l(i, j, k) = bc(i, j, k) < amrex::Real(1) ? (i - lo.x) + len.x * ((j - lo.y) + len.y * (k - lo.z)) : -1;
        });
    }
    for (;;) {
        label.FillBoundary(geom.periodicity());
        amrex::ReduceOps<amrex::ReduceOpSum> op;
        amrex::ReduceData<amrex::Long> data(op);
        for (amrex::MFIter mfi(label, false); mfi.isValid(); ++mfi) {
            const auto l = label.array(mfi);
            const auto bx = b.beta_face(0).const_array(mfi);
            const auto by = b.beta_face(1).const_array(mfi);
            const auto bz = b.beta_face(2).const_array(mfi);
            const amrex::Real one = 1;
            op.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
            {
                if (l(i, j, k) < 0) { return {amrex::Long(0)}; }
                int m = l(i, j, k);
                if (bx(i, j, k) < one) { m = amrex::max(m, l(i - 1, j, k)); }
                if (bx(i + 1, j, k) < one) { m = amrex::max(m, l(i + 1, j, k)); }
                if (by(i, j, k) < one) { m = amrex::max(m, l(i, j - 1, k)); }
                if (by(i, j + 1, k) < one) { m = amrex::max(m, l(i, j + 1, k)); }
                if (bz(i, j, k) < one && k > lo.z) { m = amrex::max(m, l(i, j, k - 1)); }
                if (bz(i, j, k + 1) < one && k < lo.z + len.z - 1) { m = amrex::max(m, l(i, j, k + 1)); }
                const bool changed = m != l(i, j, k);
                l(i, j, k) = m;
                return {changed ? amrex::Long(1) : amrex::Long(0)};
            });
        }
        amrex::Long changed = amrex::get<0>(data.value(op));
        amrex::ParallelDescriptor::ReduceLongSum(changed);
        if (changed == 0) { break; }
    }
    const int top = label.max(0);
    amrex::ReduceOps<amrex::ReduceOpSum> op;
    amrex::ReduceData<amrex::Long> data(op);
    for (amrex::MFIter mfi(label, false); mfi.isValid(); ++mfi) {
        const auto l = label.const_array(mfi);
        op.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long>
        {
            return {(l(i, j, k) >= 0 && l(i, j, k) != top) ? amrex::Long(1) : amrex::Long(0)};
        });
    }
    amrex::Long sealed = amrex::get<0>(data.value(op));
    amrex::ParallelDescriptor::ReduceLongSum(sealed);
    return sealed;
}

} // namespace gust
