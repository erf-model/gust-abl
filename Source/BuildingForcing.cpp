#include "BuildingForcing.H"

#include "ImmersedSlope.H"

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
    m_rate.define(b.beta().boxArray(), b.beta().DistributionMap(), kNumBuildingRates, 1);
    m_rate.setVal(0.0);
    m_immersed_ground = opt.immersed_ground;
    m_slope_normal = opt.immersed_ground && opt.slope_normal;
    if (m_slope_normal && !b.from_surface()) {
        amrex::Abort("gust: BuildingForcing::build: terrain.immersed_slope = normal needs a ground built from a surface");
    }
    m_spacing.clear();   // a rebuild without an immersed ground keeps no spacing
    if (m_immersed_ground && wall_law()) { make_centroid_spacing(); }
    m_step_no_slip = opt.immersed_ground && opt.step_side_walls && opt.step_no_slip && side_law();
    for (auto& f : m_no_slip) { f.clear(); }
    if (m_step_no_slip) { make_no_slip_faces(); }
    // The weight of each face in a cell's speed: 0 on a fixed row with buildings.wall_speed =
    // open_faces, else 1. A fixed row is a blocked face (fraction 1), but for w on the ground,
    // which is the boundary conditions' row (its value follows the slope). One ghost layer:
    // periodic images, else 1 (outside a non-periodic side the ghosts are no fixed rows).
    for (int d = 0; d < 3; ++d) {
        amrex::MultiFab& wt = m_weight[static_cast<std::size_t>(d)];
        wt.define(b.beta_face(d).boxArray(), b.beta_face(d).DistributionMap(), 1, 1);
        wt.setVal(1.0);
        if (opt.wall_speed != BuildingWallSpeed::open_faces) { continue; }
        for (amrex::MFIter mfi(wt, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = wt.array(mfi);
            const auto f = b.beta_face(d).const_array(mfi);
            const int kg = ground_row(d, geom);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                a(i, j, k) = (f(i, j, k) == amrex::Real(1) && k != kg) ? amrex::Real(0) : amrex::Real(1);
            });
        }
        wt.FillBoundary(geom.periodicity());
    }
    // The wall cells whose fluid centroid is within 10 z0 of a wall with the law: there the
    // law reads the speed as at 10 z0, so the stress does not follow the distance. An immersed
    // ground has no floor (form7_ustar2_exact): none is counted.
    m_wall_cells = 0;
    m_below_floor = 0;
    if (wall_law() && !opt.immersed_ground) {
        const int bits = side_law() ? (kWallRoof | kWallUnder | kWallXlo | kWallXhi | kWallYlo | kWallYhi) : (kWallRoof | kWallUnder);
        const int n0 = side_law() ? kWallD1x : kWallD1z;
        const amrex::Real floor = amrex::Real(10) * m_z0;
        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpSum> op;
        amrex::ReduceData<amrex::Long, amrex::Long> data(op);
        for (amrex::MFIter mfi(b.beta(), false); mfi.isValid(); ++mfi) {
            const auto beta = b.beta().const_array(mfi);
            const auto wl = b.wall().const_array(mfi);
            op.eval(mfi.validbox(), data, [=] AMREX_GPU_DEVICE (int i, int j, int k) -> amrex::GpuTuple<amrex::Long, amrex::Long>
            {
                if (beta(i, j, k) >= amrex::Real(1) || (static_cast<int>(wl(i, j, k, kWallMask)) & bits) == 0) {
                    return {amrex::Long(0), amrex::Long(0)};
                }
                bool below = false;
                for (int n = n0; n <= kWallD1z; n += 2) { below = below || (wl(i, j, k, n) > amrex::Real(0) && wl(i, j, k, n) < floor); }
                return {amrex::Long(1), below ? amrex::Long(1) : amrex::Long(0)};
            });
        }
        m_wall_cells = amrex::get<0>(data.value(op));
        m_below_floor = amrex::get<1>(data.value(op));
        amrex::ParallelDescriptor::ReduceLongSum(m_wall_cells);
        amrex::ParallelDescriptor::ReduceLongSum(m_below_floor);
        if (m_below_floor > 0) {
            amrex::Print() << "gust: note: buildings.z0 = " << opt.z0 << " m: " << m_below_floor << " of " << m_wall_cells
                           << " wall cells have their fluid centroid within 10 z0 of a wall; the wall law reads the speed there "
                           << "as at 10 z0, so the stress does not follow the distance (docs/buildings.rst, Limits)\n";
        }
    }
}

namespace {

//! The law at the 10 z0 floor (buildings) and at the exact distances (an immersed ground).
struct FlooredLaw {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real ustar2 (amrex::Real a, amrex::Real b, amrex::Real beta, amrex::Real d1,
                                                                    amrex::Real d2, amrex::Real z0, amrex::Real kappa) noexcept
    { return form7_ustar2(a, b, beta, d1, d2, z0, kappa); }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void heat (amrex::Real us2, amrex::Real beta, amrex::Real d1, amrex::Real d2,
                                                           amrex::Real z0, amrex::Real kappa, amrex::Real& o, amrex::Real& r) noexcept
    { form7_heat_velocities(us2, beta, d1, d2, z0, kappa, o, r); }
};
struct ExactLaw {
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static amrex::Real ustar2 (amrex::Real a, amrex::Real b, amrex::Real beta, amrex::Real d1,
                                                                    amrex::Real d2, amrex::Real z0, amrex::Real kappa) noexcept
    { return form7_ustar2_exact(a, b, beta, d1, d2, z0, kappa); }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE static void heat (amrex::Real us2, amrex::Real beta, amrex::Real d1, amrex::Real d2,
                                                           amrex::Real z0, amrex::Real kappa, amrex::Real& o, amrex::Real& r) noexcept
    { form7_heat_velocities_exact(us2, beta, d1, d2, z0, kappa, o, r); }
};

//! The rates of roofs and undersides, on u and v, with the law LAW (each law its own kernel).
template <typename LAW>
void roof_rates (const Buildings& b, const std::array<amrex::MultiFab, 3>& weight, const amrex::MultiFab& detJ,
                 const amrex::Geometry& geom, const amrex::MultiFab& u, const amrex::MultiFab& v, amrex::Real z0,
                 amrex::Real kappa, amrex::Real umin, amrex::MultiFab& rate_mf)
{
    const amrex::Real dzeta = geom.CellSize(2);
    const amrex::MultiFab& bz = b.beta_face(2);
    for (amrex::MFIter mfi(rate_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto r = rate_mf.array(mfi);
        const auto uu = u.const_array(mfi);
        const auto vv = v.const_array(mfi);
        const auto beta = b.beta().const_array(mfi);
        const auto wx = weight[0].const_array(mfi);
        const auto wy = weight[1].const_array(mfi);
        const auto wl = b.wall().const_array(mfi);
        const auto fz = bz.const_array(mfi);
        const auto J = detJ.const_array(mfi);
        const int klo = geom.Domain().smallEnd(2), khi = geom.Domain().bigEnd(2);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int mask = static_cast<int>(wl(i, j, k, kWallMask));
            if (beta(i, j, k) >= amrex::Real(1) || (mask & (kWallRoof | kWallUnder)) == 0) { return; }
            auto speed = [&] (int kk) {
                const amrex::Real uc = open_face_mean(uu(i, j, kk), uu(i + 1, j, kk), wx(i, j, kk), wx(i + 1, j, kk));
                const amrex::Real vc = open_face_mean(vv(i, j, kk), vv(i, j + 1, kk), wy(i, j, kk), wy(i, j + 1, kk));
                return std::sqrt(uc * uc + vc * vc);
            };
            const amrex::Real own = speed(k);
            const amrex::Real inv = amrex::Real(1) / (amrex::max(own, umin) * J(i, j, k) * dzeta);
            const amrex::Real d1 = wl(i, j, k, kWallD1z), d2 = wl(i, j, k, kWallD2z);
            const amrex::Real inv_dn = amrex::Real(1) / (J(i, j, k) * dzeta);
            amrex::Real rate = 0, top = 0, heat_own = 0, heat_ref = 0, ho = 0, hr = 0;
            if ((mask & kWallRoof) != 0) {
                const amrex::Real ref = k < khi ? speed(k + 1) : amrex::Real(0);
                const amrex::Real us2 = LAW::ustar2(own, ref, beta(i, j, k), d1, d2, z0, kappa);
                const amrex::Real area = wall_area_fraction(fz(i, j, k), fz(i, j, k + 1), beta(i, j, k));
                rate += area * us2 * inv;
                top = amrex::max(top, us2);
                LAW::heat(us2, beta(i, j, k), d1, d2, z0, kappa, ho, hr);
                heat_own += area * ho * inv_dn;
                heat_ref += area * hr * inv_dn;
            }
            if ((mask & kWallUnder) != 0) {
                const amrex::Real ref = k > klo ? speed(k - 1) : amrex::Real(0);
                const amrex::Real us2 = LAW::ustar2(own, ref, beta(i, j, k), d1, d2, z0, kappa);
                const amrex::Real area = wall_area_fraction(fz(i, j, k + 1), fz(i, j, k), beta(i, j, k));
                rate += area * us2 * inv;
                top = amrex::max(top, us2);
                LAW::heat(us2, beta(i, j, k), d1, d2, z0, kappa, ho, hr);
                heat_own += area * ho * inv_dn;
                heat_ref += area * hr * inv_dn;
            }
            r(i, j, k, 0) = rate;
            r(i, j, k, 1) = rate;
            r(i, j, k, kBuildingUstar2) = top;
            r(i, j, k, kBuildingHeatOwn) = heat_own;
            r(i, j, k, kBuildingHeatRefZ) = heat_ref;
        });
    }
}

}   // namespace

void BuildingForcing::update (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::MultiFab& w)
{
    if (!wall_law()) { return; }
    const amrex::Real z0 = m_z0, kappa = amrex::Real(0.41), umin = kBuildingMinSpeed;
    m_rate.setVal(0.0);
    // Roofs and undersides: on u and v.
    if (m_slope_normal) {
        // terrain.immersed_slope = normal: its own kernel (ImmersedSlope.H), on u, v and w.
        immersed_slope_rates(*m_b, m_weight, *m_detJ, m_geom, u, v, w, z0, kappa, umin, m_rate);
    } else if (m_immersed_ground) {
        roof_rates<ExactLaw>(*m_b, m_weight, *m_detJ, m_geom, u, v, z0, kappa, umin, m_rate);
    } else {
        roof_rates<FlooredLaw>(*m_b, m_weight, *m_detJ, m_geom, u, v, z0, kappa, umin, m_rate);
    }
    if (side_law() && !m_step_no_slip) {
        // With terrain.immersed_side_walls = steps_no_slip the side bits take no law: their faces diffuse
        // to the solid cell (close_viscous_fluxes).
        // Sides normal to x: on v and w. The reference cell is the next one out along x. Only a
        // partial cell weighs it (beta > 0), and such a cell has one wall along x, facing away
        // from a cell that is not solid (Buildings_K.H: wall_on_low_side; an immersed ground has no
        // side walls by default, and its roofs follow ground_roof).
        const amrex::Real dx = m_geom.CellSize(0), dy = m_geom.CellSize(1);
        for (amrex::MFIter mfi(m_rate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto r = m_rate.array(mfi);
            const auto vv = v.const_array(mfi);
            const auto ww = w.const_array(mfi);
            const auto beta = m_b->beta().const_array(mfi);
            const auto wy = m_weight[1].const_array(mfi);
            const auto wz = m_weight[2].const_array(mfi);
            const auto wl = m_b->wall().const_array(mfi);
            const auto fx = m_b->beta_face(0).const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int mask = static_cast<int>(wl(i, j, k, kWallMask));
                if (beta(i, j, k) >= amrex::Real(1) || (mask & (kWallXlo | kWallXhi)) == 0) { return; }
                auto speed = [&] (int ii) {
                    const amrex::Real vc = open_face_mean(vv(ii, j, k), vv(ii, j + 1, k), wy(ii, j, k), wy(ii, j + 1, k));
                    const amrex::Real wc = open_face_mean(ww(ii, j, k), ww(ii, j, k + 1), wz(ii, j, k), wz(ii, j, k + 1));
                    return std::sqrt(vc * vc + wc * wc);
                };
                const amrex::Real own = speed(i);
                const amrex::Real inv = amrex::Real(1) / (amrex::max(own, umin) * dx);
                const amrex::Real d1 = wl(i, j, k, kWallD1x), d2 = wl(i, j, k, kWallD2x);
                amrex::Real rate = 0, top = r(i, j, k, kBuildingUstar2), heat_own = 0, heat_ref = 0, ho = 0, hr = 0;
                if ((mask & kWallXlo) != 0) {
                    const amrex::Real us2 = form7_ustar2(own, speed(i + 1), beta(i, j, k), d1, d2, z0, kappa);
                    const amrex::Real area = wall_area_fraction(fx(i, j, k), fx(i + 1, j, k), beta(i, j, k));
                    rate += area * us2 * inv;
                    top = amrex::max(top, us2);
                    form7_heat_velocities(us2, beta(i, j, k), d1, d2, z0, kappa, ho, hr);
                    heat_own += area * ho / dx;
                    heat_ref += area * hr / dx;
                }
                if ((mask & kWallXhi) != 0) {
                    const amrex::Real us2 = form7_ustar2(own, speed(i - 1), beta(i, j, k), d1, d2, z0, kappa);
                    const amrex::Real area = wall_area_fraction(fx(i + 1, j, k), fx(i, j, k), beta(i, j, k));
                    rate += area * us2 * inv;
                    top = amrex::max(top, us2);
                    form7_heat_velocities(us2, beta(i, j, k), d1, d2, z0, kappa, ho, hr);
                    heat_own += area * ho / dx;
                    heat_ref += area * hr / dx;
                }
                r(i, j, k, 1) += rate;
                r(i, j, k, 2) += rate;
                r(i, j, k, kBuildingUstar2) = top;
                r(i, j, k, kBuildingHeatOwn) += heat_own;
                r(i, j, k, kBuildingHeatRefX) = heat_ref;
            });
        }
        // Sides normal to y: on u and w.
        for (amrex::MFIter mfi(m_rate, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto r = m_rate.array(mfi);
            const auto uu = u.const_array(mfi);
            const auto ww = w.const_array(mfi);
            const auto beta = m_b->beta().const_array(mfi);
            const auto wx = m_weight[0].const_array(mfi);
            const auto wz = m_weight[2].const_array(mfi);
            const auto wl = m_b->wall().const_array(mfi);
            const auto fy = m_b->beta_face(1).const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int mask = static_cast<int>(wl(i, j, k, kWallMask));
                if (beta(i, j, k) >= amrex::Real(1) || (mask & (kWallYlo | kWallYhi)) == 0) { return; }
                auto speed = [&] (int jj) {
                    const amrex::Real uc = open_face_mean(uu(i, jj, k), uu(i + 1, jj, k), wx(i, jj, k), wx(i + 1, jj, k));
                    const amrex::Real wc = open_face_mean(ww(i, jj, k), ww(i, jj, k + 1), wz(i, jj, k), wz(i, jj, k + 1));
                    return std::sqrt(uc * uc + wc * wc);
                };
                const amrex::Real own = speed(j);
                const amrex::Real inv = amrex::Real(1) / (amrex::max(own, umin) * dy);
                const amrex::Real d1 = wl(i, j, k, kWallD1y), d2 = wl(i, j, k, kWallD2y);
                amrex::Real rate = 0, top = r(i, j, k, kBuildingUstar2), heat_own = 0, heat_ref = 0, ho = 0, hr = 0;
                if ((mask & kWallYlo) != 0) {
                    const amrex::Real us2 = form7_ustar2(own, speed(j + 1), beta(i, j, k), d1, d2, z0, kappa);
                    const amrex::Real area = wall_area_fraction(fy(i, j, k), fy(i, j + 1, k), beta(i, j, k));
                    rate += area * us2 * inv;
                    top = amrex::max(top, us2);
                    form7_heat_velocities(us2, beta(i, j, k), d1, d2, z0, kappa, ho, hr);
                    heat_own += area * ho / dy;
                    heat_ref += area * hr / dy;
                }
                if ((mask & kWallYhi) != 0) {
                    const amrex::Real us2 = form7_ustar2(own, speed(j - 1), beta(i, j, k), d1, d2, z0, kappa);
                    const amrex::Real area = wall_area_fraction(fy(i, j + 1, k), fy(i, j, k), beta(i, j, k));
                    rate += area * us2 * inv;
                    top = amrex::max(top, us2);
                    form7_heat_velocities(us2, beta(i, j, k), d1, d2, z0, kappa, ho, hr);
                    heat_own += area * ho / dy;
                    heat_ref += area * hr / dy;
                }
                r(i, j, k, 0) += rate;
                r(i, j, k, 2) += rate;
                r(i, j, k, kBuildingUstar2) = top;
                r(i, j, k, kBuildingHeatOwn) += heat_own;
                r(i, j, k, kBuildingHeatRefY) = heat_ref;
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
            if (m_step_no_slip && d < 2) {
                // No slip on the step faces: the face's factor is 2 there (the solid cell's zero half a
                // column away), the open part elsewhere.
                for (amrex::MFIter mfi(Gcd, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const auto g = Gcd.array(mfi);
                    const auto o = open[static_cast<std::size_t>(d)]->const_array(mfi);
                    const auto ns = m_no_slip[static_cast<std::size_t>(d)].const_array(mfi);
                    amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        const amrex::Real fa = ns(i - di, j - dj, k - dk) * amrex::Real(2) + (1 - ns(i - di, j - dj, k - dk)) * o(i - di, j - dj, k - dk);
                        const amrex::Real fb = ns(i, j, k) * amrex::Real(2) + (1 - ns(i, j, k)) * o(i, j, k);
                        g(i, j, k) *= amrex::Real(0.5) * (fa + fb);
                    });
                }
                continue;
            }
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
    if (!m_spacing.empty()) {
        // An immersed ground: u and v across the top face of a roof cell, over the distance between the
        // two cells' fluid centroids (make_centroid_spacing). The face at k lies on the roof cell k - 1.
        for (int c = 0; c < 2; ++c) {
            const int di = (c == 0) ? 1 : 0, dj = (c == 1) ? 1 : 0;
            amrex::MultiFab& Gcz = G[static_cast<std::size_t>(c)][2];
            for (amrex::MFIter mfi(Gcz, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto g = Gcz.array(mfi);
                const auto sp = m_spacing.const_array(mfi);
                amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    g(i, j, k) *= amrex::Real(0.5) * (sp(i - di, j - dj, k - 1) + sp(i, j, k - 1));
                });
            }
        }
    }
}

void BuildingForcing::make_no_slip_faces ()
{
    // The x face i is a step face when cell i has the low x bit or cell i - 1 the high one (the step kernel
    // sets the bit on the cell that is not solid); y likewise. On the valid faces, then the ghost layer by
    // FillBoundary (periodic images; outside a non-periodic side 0: no step there).
    const amrex::MultiFab& wall = m_b->wall();
    for (int d = 0; d < 2; ++d) {
        amrex::MultiFab& f = m_no_slip[static_cast<std::size_t>(d)];
        f.define(m_b->beta_face(d).boxArray(), m_b->beta_face(d).DistributionMap(), 1, 1);
        f.setVal(0.0);
        const int di = (d == 0) ? 1 : 0, dj = (d == 1) ? 1 : 0;
        const int blo = (d == 0) ? kWallXlo : kWallYlo, bhi = (d == 0) ? kWallXhi : kWallYhi;
        // A face on a non-periodic side is the domain's: the wall's ghost column there copies the nearest
        // domain column (image_ghost_columns), so its bits are not read.
        const bool per = m_geom.isPeriodic(d);
        const int nlo = m_geom.Domain().smallEnd(d), nhi = m_geom.Domain().bigEnd(d) + 1;
        for (amrex::MFIter mfi(f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const auto a = f.array(mfi);
            const auto wl = wall.const_array(mfi);
            amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const int n = (d == 0) ? i : j;
                const int own = (per || n < nhi) ? static_cast<int>(wl(i, j, k, kWallMask)) : 0;
                const int low = (per || n > nlo) ? static_cast<int>(wl(i - di, j - dj, k, kWallMask)) : 0;
                a(i, j, k) = ((own & blo) != 0 || (low & bhi) != 0) ? amrex::Real(1) : amrex::Real(0);
            });
        }
        f.FillBoundary(m_geom.periodicity());
    }
}

void BuildingForcing::make_centroid_spacing ()
{
    // s = (cell-centre spacing) / (d2 - d1) on a roof cell with a reference cell above, else 1:
    //   d2 - d1 = (1 - beta) Delta_n / 2 + Delta_ref / 2,  spacing = Delta_n / 2 + Delta_ref / 2,
    // so s >= 1, and s = 1 where beta = 0. On the cells and their ghost layer (wall() has one: periodic
    // images, else the nearest domain column; below and above the domain no bits, so s = 1).
    const amrex::MultiFab& wall = m_b->wall();
    m_spacing.define(wall.boxArray(), wall.DistributionMap(), 1, 1);
    m_spacing.setVal(1.0);
    const amrex::Real dzeta = m_geom.CellSize(2);
    const amrex::Box gdom = amrex::grow(m_geom.Domain(), 2, 1);
    for (amrex::MFIter mfi(m_spacing, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.growntilebox(1) & amrex::grow(gdom, amrex::IntVect(1, 1, -1));
        const auto sp = m_spacing.array(mfi);
        const auto wl = wall.const_array(mfi);
        const auto J = m_detJ->const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const int mask = static_cast<int>(wl(i, j, k, kWallMask));
            const amrex::Real d1 = wl(i, j, k, kWallD1z), d2 = wl(i, j, k, kWallD2z);
            if ((mask & kWallRoof) == 0 || !(d2 > d1)) { return; }
            const amrex::Real centres = amrex::Real(0.5) * (J(i, j, k) + J(i, j, k + 1)) * dzeta;
            sp(i, j, k) = centres / (d2 - d1);
        });
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
            if (m_step_no_slip && d < 2) {
                // The no-slip factors of close_viscous_fluxes. Where a factor exceeds 1, a_P grows by the
                // added coupling (the flux to the solid cell's zero is implicit in the own value).
                for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                    const auto a = A.array(mfi);
                    const auto o = open[static_cast<std::size_t>(d)]->const_array(mfi);
                    const auto ns = m_no_slip[static_cast<std::size_t>(d)].const_array(mfi);
                    amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                    {
                        auto f = [&] (int ii, int jj, int kk) { return ns(ii, jj, kk) * amrex::Real(2) + (1 - ns(ii, jj, kk)) * o(ii, jj, kk); };
                        const amrex::Real flo = amrex::Real(0.5) * (f(i - di, j - dj, k - dk) + f(i, j, k));
                        const amrex::Real fhi = amrex::Real(0.5) * (f(i - di + ei, j - dj + ej, k - dk + ek) + f(i + ei, j + ej, k + ek));
                        a(i, j, k, coef_P) += amrex::max(flo - 1, amrex::Real(0)) * a(i, j, k, lo)
                                              + amrex::max(fhi - 1, amrex::Real(0)) * a(i, j, k, hi);
                        a(i, j, k, lo) *= flo;
                        a(i, j, k, hi) *= fhi;
                    });
                }
                continue;
            }
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
        if (!m_spacing.empty() && c < 2) {
            // The centroid spacing of close_viscous_fluxes: the couplings across a roof cell's top face grow by
            // s >= 1, and a_P by the same amount, so a_P = sum a_nb still holds.
            for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
                const auto a = A.array(mfi);
                const auto sp = m_spacing.const_array(mfi);
                amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const amrex::Real sb = amrex::Real(0.5) * (sp(i - di, j - dj, k - 1) + sp(i, j, k - 1));
                    const amrex::Real st = amrex::Real(0.5) * (sp(i - di, j - dj, k) + sp(i, j, k));
                    a(i, j, k, coef_P) += (sb - 1) * a(i, j, k, coef_ZM) + (st - 1) * a(i, j, k, coef_ZP);
                    a(i, j, k, coef_ZM) *= sb;
                    a(i, j, k, coef_ZP) *= st;
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
