#include "NestDefect.H"

#include "MassFlux.H"
#include "Reductions.H"

#include <AMReX_ParallelDescriptor.H>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>

namespace gust {

namespace {

//! The parent cells under the nest's boxes, at type `like`'s, on the nest's DistributionMapping.
amrex::MultiFab image (const amrex::BoxArray& fine_ba, const amrex::DistributionMapping& dm, int r,
                       const amrex::IntVect& offset, const amrex::IndexType& type)
{
    amrex::BoxArray img = amrex::coarsen(fine_ba, amrex::IntVect(r, r, 1));
    img.shift(amrex::IntVect(offset[0], offset[1], 0));
    return amrex::MultiFab(amrex::convert(img, type), dm, 1, 0);
}

//! c = sum J q / sum J over the r x r nest cells of each parent image cell.
void restrict_cells (const amrex::MultiFab& J, const amrex::MultiFab& q, int r, const amrex::IntVect& offset,
                     amrex::MultiFab& c)
{
    const int ox = offset[0], oy = offset[1];
    for (amrex::MFIter mfi(c, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = c.array(mfi);
        const auto jj = J.const_array(mfi);
        const auto qq = q.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int I, int Jc, int k) noexcept
        {
            amrex::Real sjq = 0, sj = 0;
            for (int b = 0; b < r; ++b) {
                for (int a = 0; a < r; ++a) {
                    const int i = r * (I - ox) + a, j = r * (Jc - oy) + b;
                    sjq += jj(i, j, k) * qq(i, j, k);
                    sj += jj(i, j, k);
                }
            }
            o(I, Jc, k) = sjq / sj;
        });
    }
}

//! c = (1/r) sum a q over the r nest x-faces (D = 0) or y-faces (D = 1) of
//! each parent image face, over the parent's area ac there.
template <int D>
void restrict_faces (const amrex::MultiFab& area, const amrex::MultiFab& q, int r, const amrex::IntVect& offset,
                     const amrex::MultiFab& ac, amrex::MultiFab& c)
{
    const int ox = offset[0], oy = offset[1];
    const amrex::Real rinv = amrex::Real(1) / amrex::Real(r);
    for (amrex::MFIter mfi(c, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = c.array(mfi);
        const auto aa = area.const_array(mfi);
        const auto qq = q.const_array(mfi);
        const auto cc = ac.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int I, int J, int k) noexcept
        {
            amrex::Real s = 0;
            for (int b = 0; b < r; ++b) {
                const int i = (D == 0) ? r * (I - ox) : r * (I - ox) + b;
                const int j = (D == 0) ? r * (J - oy) + b : r * (J - oy);
                s += aa(i, j, k) * qq(i, j, k);
            }
            o(I, J, k) = s * rinv / cc(I, J, k);
        });
    }
}

//! c = (1/r^2) sum a_z Omega over the r x r nest z-faces of each parent image
//! z-face, over the parent's a_z there.
void restrict_zfaces (const amrex::MultiFab& az, const amrex::MultiFab& omega, int r, const amrex::IntVect& offset,
                      const amrex::MultiFab& ac, amrex::MultiFab& c)
{
    const int ox = offset[0], oy = offset[1];
    const amrex::Real r2inv = amrex::Real(1) / amrex::Real(r * r);
    for (amrex::MFIter mfi(c, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.tilebox();
        const auto o = c.array(mfi);
        const auto aa = az.const_array(mfi);
        const auto qq = omega.const_array(mfi);
        const auto cc = ac.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int I, int J, int k) noexcept
        {
            amrex::Real s = 0;
            for (int b = 0; b < r; ++b) {
                for (int a = 0; a < r; ++a) {
                    const int i = r * (I - ox) + a, j = r * (J - oy) + b;
                    s += aa(i, j, k) * qq(i, j, k);
                }
            }
            o(I, J, k) = s * r2inv / cc(I, J, k);
        });
    }
}

//! dst = src on the valid points of dst inside region (same layout).
void copy_in (const amrex::MultiFab& src, const amrex::Box& region, amrex::MultiFab& dst)
{
    for (amrex::MFIter mfi(dst, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & region;
        if (!bx.ok()) { continue; }
        const auto d = dst.array(mfi);
        const auto s = src.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { d(i, j, k) = s(i, j, k); });
    }
}

//! max |src - dst| over the points of region (dst's valid points; a NaN counts as infinite), over all ranks.
amrex::Real max_difference_in (const amrex::MultiFab& src, const amrex::Box& region, const amrex::MultiFab& dst)
{
    amrex::ReduceOps<amrex::ReduceOpMax> op;
    amrex::ReduceData<amrex::Real> data(op);
    for (amrex::MFIter mfi(dst, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & region;
        if (!bx.ok()) { continue; }
        const auto d = dst.const_array(mfi);
        const auto s = src.const_array(mfi);
        op.eval(bx, data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> amrex::GpuTuple<amrex::Real>
        {
            const amrex::Real e = std::abs(s(i, j, k) - d(i, j, k));
            return {e == e ? e : std::numeric_limits<amrex::Real>::infinity()};
        });
    }
    amrex::Real m = amrex::get<0>(data.value(op));
    amrex::ParallelDescriptor::ReduceRealMax(m);
    return m;
}

void restrict_horizontal (const State& fine, const Mesh& mesh_f, int r, const amrex::IntVect& offset,
                          const Mesh& mesh_c, const amrex::Geometry& geom_c, State& coarse)
{
    const amrex::BoxArray& fba = fine.p.boxArray();
    const amrex::DistributionMapping& fdm = fine.p.DistributionMap();
    if (!fba.coarsenable(amrex::IntVect(r, r, 1))) {
        amrex::Abort("gust: restrict_horizontal: the nest's boxes do not coarsen by " + std::to_string(r));
    }
    const amrex::IndexType CC = amrex::IndexType::TheCellType();
    const amrex::IndexType XF(amrex::IntVect(1, 0, 0)), YF(amrex::IntVect(0, 1, 0));

    // Cell centres: the detJ-weighted mean on the image, then copied to the parent (valid points; the
    // footprint lies inside the parent's domain, so no periodic shift).
    {
        amrex::MultiFab img = image(fba, fdm, r, offset, CC);
        const std::array<std::pair<const amrex::MultiFab*, amrex::MultiFab*>, 3> cells{{
            {&fine.theta, &coarse.theta}, {&fine.tke, &coarse.tke}, {&fine.p, &coarse.p}}};
        for (const auto& [src, dst] : cells) {
            restrict_cells(mesh_f.detJ, *src, r, offset, img);
            if (src == &fine.p) {
                amrex::MultiFab pc(img.boxArray(), img.DistributionMap(), 1, 0);
                pc.ParallelCopy(coarse.p, 0, 0, 1);
                // The level shift summed on one box over the footprint, in a fixed order: the same bits for
                // any nest layout or rank map (MultiFab::sum adds box by box, so its last bits follow the boxes).
                amrex::MultiFab diff(img.boxArray(), img.DistributionMap(), 1, 0);
                amrex::MultiFab::LinComb(diff, 1, pc, 0, -1, img, 0, 0, 1, 0);
                const amrex::BoxArray one(img.boxArray().minimalBox());
                amrex::MultiFab whole(one, amrex::DistributionMapping(one), 1, 0);
                whole.ParallelCopy(diff, 0, 0, 1);
                const amrex::Real d = deterministic_sum(whole, 0) / amrex::Real(one.numPts());
                img.plus(d, 0, 1, 0);
            }
            dst->ParallelCopy(img, 0, 0, 1);
        }
    }
    // u, v: the flux through the r nest faces over the parent face's area. ax, ay and az are the full face areas;
    // with partly open faces (buildings) the divergence uses flux_area, and this restriction would not be
    // divergence-free. App.cpp refuses buildings with nests.
    {
        amrex::MultiFab img = image(fba, fdm, r, offset, XF);
        amrex::MultiFab ac = image(fba, fdm, r, offset, XF);
        ac.ParallelCopy(mesh_c.ax, 0, 0, 1);
        restrict_faces<0>(mesh_f.ax, fine.u, r, offset, ac, img);
        // With the periodicity: a footprint face on a periodic seam (face 0 or N) also sets its image.
        coarse.u.ParallelCopy(img, 0, 0, 1, 0, 0, geom_c.periodicity());
    }
    {
        amrex::MultiFab img = image(fba, fdm, r, offset, YF);
        amrex::MultiFab ac = image(fba, fdm, r, offset, YF);
        ac.ParallelCopy(mesh_c.ay, 0, 0, 1);
        restrict_faces<1>(mesh_f.ay, fine.v, r, offset, ac, img);
        coarse.v.ParallelCopy(img, 0, 0, 1, 0, 0, geom_c.periodicity());
    }
}

void restrict_w (const State& fine, const Mesh& mesh_f, const amrex::Geometry& geom_f, int r,
                 const amrex::IntVect& offset, const Mesh& mesh_c, const amrex::Geometry& geom_c, State& coarse)
{
    const amrex::BoxArray& fba = fine.p.boxArray();
    const amrex::DistributionMapping& fdm = fine.p.DistributionMap();
    if (!fba.coarsenable(amrex::IntVect(r, r, 1))) {
        amrex::Abort("gust: restrict_w: the nest's boxes do not coarsen by " + std::to_string(r));
    }
    const amrex::IndexType ZF(amrex::IntVect(0, 0, 1));
    // w from the restricted vertical flux Omega: the restricted field is divergence-free on the parent cells
    // when the nest's field is on its own (the u, v fluxes above and Omega are conservative across the
    // different ground), except next to a face another nest's restriction set last (restrict_patches).
    {
        amrex::MultiFab om_f(fine.w.boxArray(), fdm, 1, 0);
        omega_from_w(fine.u, fine.v, fine.w, fine.z_phys_nd, mesh_f, geom_f, om_f, 0);
        amrex::MultiFab img = image(fba, fdm, r, offset, ZF);
        amrex::MultiFab ac = image(fba, fdm, r, offset, ZF);
        ac.ParallelCopy(mesh_c.az, 0, 0, 1);
        restrict_zfaces(mesh_f.az, om_f, r, offset, ac, img);
        amrex::MultiFab om_c(coarse.w.boxArray(), coarse.w.DistributionMap(), 1, 0);
        omega_from_w(coarse.u, coarse.v, coarse.w, coarse.z_phys_nd, mesh_c, geom_c, om_c, 0);
        om_c.ParallelCopy(img, 0, 0, 1);
        amrex::MultiFab wc(coarse.w.boxArray(), coarse.w.DistributionMap(), 1, 0);
        w_from_omega(coarse.u, coarse.v, om_c, coarse.z_phys_nd, mesh_c, geom_c, wc);
        amrex::Box fp = amrex::coarsen(geom_f.Domain(), amrex::IntVect(r, r, 1));
        fp.shift(amrex::IntVect(offset[0], offset[1], 0));
        copy_in(wc, amrex::convert(fp, ZF), coarse.w);
    }
}

} // namespace

amrex::Real relative_feedback_change (const FeedbackChange& c, const std::array<amrex::Real, 5>& floor)
{
    const amrex::Real inf = std::numeric_limits<amrex::Real>::infinity();
    if (c.first) { return inf; }
    amrex::Real r = 0;
    for (std::size_t n = 0; n < 5; ++n) {
        const amrex::Real d = amrex::max(c.size[n], floor[n]);
        const amrex::Real x = c.change[n];
        // The quotient's denominator is never 0: clang evaluates it before the branch, and 0 / 0 raises FE_INVALID
        // (trapped with amrex.fpe_trap_invalid). max(d, min) cannot fold back into d, as a guard d > 0 ? d : 1 does.
        const amrex::Real den = amrex::max(d, std::numeric_limits<amrex::Real>::min());
        const amrex::Real q = std::isnan(x) ? inf : d > 0 ? x / den : x > 0 ? inf : amrex::Real(0);
        r = amrex::max(r, q);
    }
    return r;
}

void FeedbackWindow::add (amrex::Real d)
{
    if (!(d >= 0 && d < std::numeric_limits<amrex::Real>::infinity())) {
        m_d.clear();
        return;
    }
    m_d.push_back(d);
    while (static_cast<int>(m_d.size()) > m_window) { m_d.pop_front(); }
}

amrex::Real FeedbackWindow::sum () const
{
    if (static_cast<int>(m_d.size()) < m_window) { return std::numeric_limits<amrex::Real>::infinity(); }
    amrex::Real s = 0;
    for (const amrex::Real d : m_d) { s += d; }
    return s;
}

bool FeedbackWindow::settled () const
{
    if (!(sum() <= m_tol)) { return false; }
    const std::size_t W = m_d.size();
    const std::size_t half = W / 2;
    amrex::Real d1 = 0, d2 = 0;
    for (std::size_t n = W - 2 * half; n < W - half; ++n) { d1 += m_d[n]; }
    for (std::size_t n = W - half; n < W; ++n) { d2 += m_d[n]; }
    if (std::max(d1, d2) <= amrex::Real(1.0e-3) * m_tol) { return true; }
    if (!(d2 < d1)) { return false; }
    const amrex::Real rho = d2 / d1;
    return d2 * rho / (1 - rho) <= m_tol;
}

void restrict_patches (const std::vector<NestImage>& nests, const Mesh& mesh_c, const amrex::Geometry& geom_c,
                       const std::function<void(State&)>& coarse_ghosts, State& coarse)
{
    for (const NestImage& n : nests) {
        restrict_horizontal(*n.fine, *n.mesh, n.ratio, n.offset, mesh_c, geom_c, coarse);
    }
    // The parent's ghost cells with the new u and v, for its Omega and for w_from_omega.
    coarse_ghosts(coarse);
    for (const NestImage& n : nests) {
        restrict_w(*n.fine, *n.mesh, *n.geom, n.ratio, n.offset, mesh_c, geom_c, coarse);
    }
    coarse_ghosts(coarse);
}

void restrict_patch (const State& fine, const Mesh& mesh_f, const amrex::Geometry& geom_f, int r,
                     const amrex::IntVect& offset, const Mesh& mesh_c, const amrex::Geometry& geom_c,
                     const std::function<void(State&)>& coarse_ghosts, State& coarse)
{
    restrict_patches({NestImage{&fine, &mesh_f, &geom_f, r, offset}}, mesh_c, geom_c, coarse_ghosts, coarse);
}

DefectRegion defect_region (const std::string& name, const amrex::IntVect& fine_cells, int r, const amrex::IntVect& offset,
                            const amrex::Box& dom_c, std::array<bool, 2> periodic)
{
    if (fine_cells[0] % r != 0 || fine_cells[1] % r != 0) {
        amrex::Abort("gust: " + name + ": the nest's cells do not coarsen by " + std::to_string(r));
    }
    DefectRegion g;
    g.name = name;
    const amrex::IntVect lo(offset[0], offset[1], dom_c.smallEnd(2));
    const amrex::IntVect hi(offset[0] + fine_cells[0] / r - 1, offset[1] + fine_cells[1] / r - 1, dom_c.bigEnd(2));
    g.footprint = amrex::Box(lo, hi);
    g.inner = g.footprint;
    for (int d = 0; d < 2; ++d) {
        if (!periodic[d]) { g.inner.grow(d, -kDefectMargin); }
    }
    if (!g.inner.ok()) {
        amrex::Abort("gust: " + name + ": with nest.feedback = defect or relax a nest spans at least "
                     + std::to_string(2 * kDefectMargin + 1) + " parent cells in x and y (the source keeps "
                     + std::to_string(kDefectMargin) + " cells from its faces)");
    }
    return g;
}

NestDefect::NestDefect (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    m_S.define(ba, dm);
    for (auto& f : m_S.fields()) { f.mf->setVal(0.0); }
}

void NestDefect::set (const ResidualFields& R, const std::vector<DefectRegion>& regions)
{
    amrex::MultiFab* S[5] = {&m_S.u, &m_S.v, &m_S.w, &m_S.theta, &m_S.tke};
    const amrex::MultiFab* r[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    for (int n = 0; n < 5; ++n) {
        if (r[n]->boxArray() != S[n]->boxArray() || r[n]->DistributionMap() != S[n]->DistributionMap()) {
            amrex::Abort("gust: NestDefect::set: the residual is not on the parent's layout");
        }
        // The change from the S of the last set (zero outside the inner boxes, before and after).
        amrex::Real ch = 0;
        for (const auto& g : regions) { ch = amrex::max(ch, max_difference_in(*r[n], amrex::convert(g.inner, S[n]->ixType()), *S[n])); }
        m_change.change[static_cast<std::size_t>(n)] = ch;
        S[n]->setVal(0.0);
        for (const auto& g : regions) { copy_in(*r[n], amrex::convert(g.inner, S[n]->ixType()), *S[n]); }
        m_change.size[static_cast<std::size_t>(n)] = S[n]->norm0();
    }
    m_change.first = !(m_active || m_have_previous);
    m_active = true;
}

void NestDefect::restore_previous (const std::array<const amrex::MultiFab*, 5>& f)
{
    amrex::MultiFab* dst[5] = {&m_S.u, &m_S.v, &m_S.w, &m_S.theta, &m_S.tke};
    for (std::size_t n = 0; n < 5; ++n) {
        amrex::MultiFab& d = *dst[n];
        if (f[n]->boxArray() != d.boxArray() || f[n]->DistributionMap() != d.DistributionMap()) {
            amrex::Abort("gust: NestDefect::restore_previous: the field is not on the parent's layout");
        }
        amrex::MultiFab::Copy(d, *f[n], 0, 0, 1, 0);
    }
    m_have_previous = true;
}

void NestDefect::add (const ResidualFields& C)
{
    if (!m_active) { amrex::Abort("gust: NestDefect::add: the source is not set"); }
    amrex::MultiFab* S[5] = {&m_S.u, &m_S.v, &m_S.w, &m_S.theta, &m_S.tke};
    const amrex::MultiFab* c[5] = {&C.u, &C.v, &C.w, &C.theta, &C.tke};
    for (int n = 0; n < 5; ++n) {
        if (c[n]->boxArray() != S[n]->boxArray() || c[n]->DistributionMap() != S[n]->DistributionMap()) {
            amrex::Abort("gust: NestDefect::add: the correction is not on the parent's layout");
        }
        amrex::MultiFab::Add(*S[n], *c[n], 0, 0, 1, 0);
    }
}

void NestDefect::add_to_residual (ResidualFields& R) const
{
    if (!m_active) { return; }
    const amrex::MultiFab* S[5] = {&m_S.u, &m_S.v, &m_S.w, &m_S.theta, &m_S.tke};
    amrex::MultiFab* r[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    for (int n = 0; n < 5; ++n) {
        if (r[n]->boxArray() != S[n]->boxArray() || r[n]->DistributionMap() != S[n]->DistributionMap()) {
            amrex::Abort("gust: NestDefect::add_to_residual: the residual is not on the parent's layout");
        }
        amrex::MultiFab::Subtract(*r[n], *S[n], 0, 0, 1, 0);
    }
}

std::vector<amrex::Real> NestDefect::max_norms () const
{
    return {m_S.u.norm0(), m_S.v.norm0(), m_S.w.norm0(), m_S.theta.norm0(), m_S.tke.norm0()};
}

} // namespace gust
