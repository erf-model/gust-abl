#include "Buildings.H"

#include "Buildings_K.H"
#include "ImmersedSurface_K.H"
#include "SurfaceRoughness.H"

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>

#include <algorithm>
#include <limits>

// Buildings::build_from_surface: the geometry of an immersed ground that is a continuous surface
// (Buildings.H, ImmersedSurface_K.H; docs/immersed_terrain.rst, Geometry).

namespace gust {

namespace {

//! The mesh's levels (node means over a cell's footprint or a face's edge) and the ground's height
//! at the nodes of a column, in double.
struct SurfaceNodes
{
    amrex::Array4<amrex::Real const> z;
    int k0, gc;
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double g (int i, int j) const noexcept { return static_cast<double>(z(i, j, k0, gc)); }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double h4 (int i, int j, int k) const noexcept
    {
        return 0.25 * (static_cast<double>(z(i, j, k)) + static_cast<double>(z(i + 1, j, k))
                       + static_cast<double>(z(i, j + 1, k)) + static_cast<double>(z(i + 1, j + 1, k)));
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double hx (int i, int j, int k) const noexcept
    {
        return 0.5 * (static_cast<double>(z(i, j, k)) + static_cast<double>(z(i, j + 1, k)));
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    double hy (int i, int j, int k) const noexcept
    {
        return 0.5 * (static_cast<double>(z(i, j, k)) + static_cast<double>(z(i + 1, j, k)));
    }
    //! The means of (t - h)_+^n over the footprint of cell (i, j).
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    SurfaceMoments moments (int i, int j, double t) const noexcept
    {
        return footprint_moments(g(i, j), g(i + 1, j), g(i + 1, j + 1), g(i, j + 1), t);
    }
};

} // namespace

void Buildings::build_from_surface (const BuildingsOptions& opt, const std::string& where, const amrex::Geometry& geom,
                                    const amrex::MultiFab& z_nd, const amrex::BoxArray& cell_ba,
                                    const amrex::DistributionMapping& dm)
{
    if (!opt.active() || !opt.immersed_ground) {
        amrex::Abort("gust: Buildings::build_from_surface needs the options of an immersed ground");
    }
    if (z_nd.nComp() < 2) {
        amrex::Abort("gust: Buildings::build_from_surface needs the ground's height as the node heights' second component");
    }
    m_active = true;
    m_hash.clear();
    for (int n = 0; n < cell_ba.size(); ++n) {
        if (cell_ba[n].length(2) != geom.Domain().length(2)) {
            amrex::Abort("gust: " + where + ": every box must span the column (the mesh's layout does)");
        }
    }
    const auto dx = geom.CellSizeArray();
    const double dxd = static_cast<double>(dx[0]), dyd = static_cast<double>(dx[1]);
    // At least the round-off of a difference of two depths, so small_volfrac = 0 still gives solid cells.
    const double small = std::max(opt.small_volfrac, 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()));
    const int k0 = geom.Domain().smallEnd(2), ktop = geom.Domain().bigEnd(2);
    const int gc = z_nd.nComp() - 1;

    // The cells: the solid fraction, and the surface's fields from the unclipped geometry. Ghost
    // layers below and above the domain stay 0.
    m_beta.clear();
    m_beta.define(cell_ba, dm, 1, 1);
    m_beta.setVal(0);
    m_surface.clear();
    m_surface.define(cell_ba, dm, kNumSurfaceComps, 1);
    m_surface.setVal(0);
    for (amrex::MFIter mfi(m_beta, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto b = m_beta.array(mfi);
        const auto s = m_surface.array(mfi);
        const SurfaceNodes z{z_nd.const_array(mfi), k0, gc};
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double zb = z.h4(i, j, k), zt = z.h4(i, j, k + 1), dz = zt - zb;
            const SurfaceMoments mb = z.moments(i, j, zb), mt = z.moments(i, j, zt);
            const double beta = clip_fraction(solid_fraction(mb.m1, mt.m1, zb, zt), small);
            b(i, j, k) = static_cast<amrex::Real>(beta);
            // The mean height of the cell's fluid above the ground; half the layer in a solid cell,
            // which has no fluid (a positive length for the closure's solid rows).
            const double depth = mt.m1 - mb.m1;
            const bool fluid = beta < 1 && depth > 0;
            s(i, j, k, kSurfaceHeight) = static_cast<amrex::Real>(fluid ? fluid_height(mb, mt) : 0.5 * dz);
            // The ground's mean slope over the footprint, from its four nodes.
            const double h00 = z.g(i, j), h10 = z.g(i + 1, j), h11 = z.g(i + 1, j + 1), h01 = z.g(i, j + 1);
            s(i, j, k, kSurfaceSlopeX) = static_cast<amrex::Real>(((h10 - h00) + (h11 - h01)) / (2 * dxd));
            s(i, j, k, kSurfaceSlopeY) = static_cast<amrex::Real>(((h01 - h00) + (h11 - h10)) / (2 * dyd));
        });
    }
    amrex::Gpu::streamSynchronize();
    image_ghost_columns(m_beta, geom);
    image_ghost_columns(m_surface, geom);

    // The ground through the top layer: no fluid above it for the wall law.
    amrex::ReduceOps<amrex::ReduceOpMax> rop;
    amrex::ReduceData<amrex::Real> rdata(rop);
    for (amrex::MFIter mfi(m_beta, false); mfi.isValid(); ++mfi) {
        const auto b = m_beta.const_array(mfi);
        rop.eval(mfi.validbox(), rdata, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real>
        {
            return {k == ktop ? b(i, j, k) : amrex::Real(0)};
        });
    }
    amrex::Real top = amrex::get<0>(rdata.value(rop));
    amrex::ParallelDescriptor::ReduceRealMax(top);
    if (top > 0) {
        amrex::Abort("gust: " + where + ": the ground reaches the top layer of the mesh; raise the top of the domain or lower the ground");
    }

    // The faces. The round-off of a height is several ulp of the top of the mesh.
    const double ztol = 64 * static_cast<double>(std::numeric_limits<amrex::Real>::epsilon()) * static_cast<double>(z_nd.max(0));
    for (int d = 0; d < 3; ++d) {
        amrex::IntVect nodal(0);
        nodal[d] = 1;
        auto& bf = m_beta_face[static_cast<std::size_t>(d)];
        bf.clear();
        bf.define(amrex::convert(cell_ba, nodal), dm, 1, 0);
        bf.setVal(0);
        const amrex::IntVect e = amrex::IntVect::TheDimensionVector(d);
        const int ei = e[0], ej = e[1], ek = e[2];
        for (amrex::MFIter mfi(bf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const amrex::Box bx = mfi.tilebox();
            const auto a = bf.array(mfi);
            const auto b = m_beta.const_array(mfi);
            const SurfaceNodes z{z_nd.const_array(mfi), k0, gc};
            // A face next to a solid cell is blocked and one between two fluid cells open, as the
            // buildings' faces: faces and cells agree after the clipping.
            auto finish = [=] AMREX_GPU_DEVICE (int i, int j, int k, double m) noexcept
            {
                const double lo = static_cast<double>(b(i - ei, j - ej, k - ek)), hi = static_cast<double>(b(i, j, k));
                if (lo == 1 || hi == 1) { return amrex::Real(1); }
                if (lo == 0 && hi == 0) { return amrex::Real(0); }
                return static_cast<amrex::Real>(clip_fraction(amrex::min(1.0, amrex::max(0.0, m)), small));
            };
            if (d == 0) {
                // The ground along the face's edge runs linearly between its two nodes.
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const double zb = z.hx(i, j, k), zt = z.hx(i, j, k + 1);
                    const double open = (segment_depth(z.g(i, j), z.g(i, j + 1), zt) - segment_depth(z.g(i, j), z.g(i, j + 1), zb)) / (zt - zb);
                    a(i, j, k) = finish(i, j, k, 1 - open);
                });
            } else if (d == 1) {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    const double zb = z.hy(i, j, k), zt = z.hy(i, j, k + 1);
                    const double open = (segment_depth(z.g(i, j), z.g(i + 1, j), zt) - segment_depth(z.g(i, j), z.g(i + 1, j), zb)) / (zt - zb);
                    a(i, j, k) = finish(i, j, k, 1 - open);
                });
            } else {
                // The part of the footprint whose ground lies at or above the face (ztol: a flat
                // ground on a face blocks it).
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    a(i, j, k) = finish(i, j, k, 1 - z.moments(i, j, z.h4(i, j, k) - ztol).m0);
                });
            }
        }
        amrex::Gpu::streamSynchronize();
        // Every copy of a face on a box edge or a periodic seam holds one value.
        bf.OverrideSync(geom.periodicity());
    }

    // The wall cells: the ground is the wall below a cell (the roof bit) wherever the blocked part
    // of its low face exceeds that of its high face, with no tolerance (ground_roof), so the roofs
    // of a column add to its footprint. The ground has no side walls: on a slope the
    // differences of the x and y faces are the same surface, whose area the low and high z faces
    // already count.
    m_wall.clear();
    m_wall.define(cell_ba, dm, kNumWallComps, 1);
    m_wall.setVal(0);
    for (amrex::MFIter mfi(m_wall, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto w = m_wall.array(mfi);
        const auto br = m_beta.const_array(mfi);
        const auto fz = m_beta_face[2].const_array(mfi);
        const auto s = m_surface.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const double own = static_cast<double>(br(i, j, k));
            if (own == 1) { return; }
            const double zl = static_cast<double>(fz(i, j, k)), zh = static_cast<double>(fz(i, j, k + 1));
            if (!ground_roof(zl, zh)) { return; }
            w(i, j, k, kWallMask) = static_cast<amrex::Real>(kWallRoof);
            // The distances of the wall law: the mean heights above the ground of the cell's fluid
            // and of the fluid of the cell above (0 when that one is solid or outside the domain).
            w(i, j, k, kWallD1z) = s(i, j, k, kSurfaceHeight);
            const bool ok = k < ktop && br(i, j, k + 1) < amrex::Real(1);
            w(i, j, k, kWallD2z) = ok ? s(i, j, k + 1, kSurfaceHeight) : amrex::Real(0);
        });
    }
    amrex::Gpu::streamSynchronize();
    image_ghost_columns(m_wall, geom);
}

} // namespace gust
