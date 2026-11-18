#include "NestInject.H"

#include "Reductions.H"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>

namespace gust {

namespace {

void require_layout (const char* op, const amrex::MultiFab& a, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
{
    if (amrex::convert(a.boxArray(), amrex::IndexType::TheCellType()) != ba || a.DistributionMap() != dm) {
        amrex::Abort(std::string("gust: NestInject::") + op + ": the field is not on the parent's layout");
    }
}

//! The cell box region converted to type t, less the ground and top z faces of the cell box domain: w there is
//! the boundary condition's (on the terrain it follows u), not data.
amrex::Box data_box (const amrex::Box& region, amrex::IndexType t, const amrex::Box& domain)
{
    amrex::Box pb = amrex::convert(region, t);
    if (t.nodeCentered(2)) {
        pb.setSmall(2, std::max(pb.smallEnd(2), domain.smallEnd(2) + 1));
        pb.setBig(2, std::min(pb.bigEnd(2), domain.bigEnd(2)));
    }
    return pb;
}

//! dst = src on the valid points of dst inside data_box(region).
void copy_region (const amrex::MultiFab& src, const amrex::Box& region, const amrex::Box& domain, amrex::MultiFab& dst)
{
    const amrex::Box pb = data_box(region, dst.ixType(), domain);
    for (amrex::MFIter mfi(dst, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & pb;
        if (!bx.ok()) { continue; }
        const auto d = dst.array(mfi);
        const auto s = src.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { d(i, j, k) = s(i, j, k); });
    }
}

//! dst = 0 on the valid points of dst inside the cell box region converted to dst's type.
void zero_region (const amrex::Box& region, amrex::MultiFab& dst)
{
    const amrex::Box pb = amrex::convert(region, dst.ixType());
    for (amrex::MFIter mfi(dst, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & pb;
        if (!bx.ok()) { continue; }
        const auto d = dst.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { d(i, j, k) = amrex::Real(0); });
    }
}

//! d = |a - b| and m = 1 on the valid points inside data_box(region) (d and m keep their values elsewhere).
void difference_region (const amrex::MultiFab& a, const amrex::MultiFab& b, const amrex::Box& region,
                        const amrex::Box& domain, amrex::MultiFab& d, amrex::MultiFab& m)
{
    const amrex::Box pb = data_box(region, d.ixType(), domain);
    for (amrex::MFIter mfi(d, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & pb;
        if (!bx.ok()) { continue; }
        const auto dd = d.array(mfi);
        const auto mm = m.array(mfi);
        const auto aa = a.const_array(mfi);
        const auto bb = b.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            dd(i, j, k) = std::abs(aa(i, j, k) - bb(i, j, k));
            mm(i, j, k) = amrex::Real(1);
        });
    }
}

//! a = |s| on the valid points inside data_box(region) (a keeps its values elsewhere).
void abs_region (const amrex::MultiFab& s, const amrex::Box& region, const amrex::Box& domain, amrex::MultiFab& a)
{
    const amrex::Box pb = data_box(region, a.ixType(), domain);
    for (amrex::MFIter mfi(a, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & pb;
        if (!bx.ok()) { continue; }
        const auto aa = a.array(mfi);
        const auto ss = s.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { aa(i, j, k) = std::abs(ss(i, j, k)); });
    }
}

//! t = q + omega (t - q) on the valid points of t inside data_box(region).
void blend_region (const amrex::MultiFab& q, amrex::Real omega, const amrex::Box& region, const amrex::Box& domain,
                   amrex::MultiFab& t)
{
    const amrex::Box pb = data_box(region, t.ixType(), domain);
    for (amrex::MFIter mfi(t, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & pb;
        if (!bx.ok()) { continue; }
        const auto ta = t.array(mfi);
        const auto qa = q.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ta(i, j, k) = qa(i, j, k) + omega * (ta(i, j, k) - qa(i, j, k));
        });
    }
}

} // namespace

std::string check_inject_factor (amrex::Real factor)
{
    if (std::isfinite(static_cast<double>(factor)) && factor > 0 && factor <= 1) { return {}; }
    std::ostringstream os;
    os << "nest.feedback_factor = " << factor << " with nest.feedback = inject: must be finite, > 0 and <= 1";
    return os.str();
}

NestInject::NestInject (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::Geometry& geom,
                        amrex::Real factor)
    : m_ba(amrex::convert(ba, amrex::IndexType::TheCellType())), m_dm(dm), m_geom(geom), m_factor(factor)
{
    const std::string err = check_inject_factor(factor);
    if (!err.empty()) { amrex::Abort("gust: " + err); }
    // u, v, w on their faces; p, theta, k on the cells.
    for (int n = 0; n < 6; ++n) {
        const amrex::BoxArray b = n < 3 ? amrex::convert(m_ba, amrex::IntVect::TheDimensionVector(n)) : m_ba;
        m_target[static_cast<std::size_t>(n)].define(b, m_dm, 1, 0);
        m_target[static_cast<std::size_t>(n)].setVal(0.0);
    }
    m_change.first = false;
}

void NestInject::restore_previous (const std::array<const amrex::MultiFab*, 5>& f)
{
    const int ti[5] = {0, 1, 2, 4, 5};
    for (int n = 0; n < 5; ++n) {
        amrex::MultiFab& t = m_target[static_cast<std::size_t>(ti[n])];
        if (f[n]->boxArray() != t.boxArray() || f[n]->DistributionMap() != t.DistributionMap()) {
            amrex::Abort("gust: NestInject::restore_previous: the field is not on the parent's layout");
        }
        amrex::MultiFab::Copy(t, *f[n], 0, 0, 1, 0);
    }
}

void NestInject::set (const std::vector<DefectRegion>& regions)
{
    for (std::size_t n = 0; n < regions.size(); ++n) {
        const auto& g = regions[n];
        if (!g.inner.ok() || !g.footprint.contains(g.inner) || !m_geom.Domain().contains(g.footprint)) {
            amrex::Abort("gust: NestInject::set: " + g.name + ": the inner region is not inside the footprint and the "
                         "parent's domain");
        }
        // Shared points would take omega twice (Nests.cpp refuses overlapping nests; this keeps the class honest).
        for (std::size_t m = 0; m < n; ++m) {
            if (g.inner.intersects(regions[m].inner)) {
                amrex::Abort("gust: NestInject::set: the inner regions of " + regions[m].name + " and " + g.name
                             + " overlap");
            }
        }
    }
    m_regions = regions;
    m_fluid.define(m_ba, m_dm, 1, 0);
    m_fluid.setVal(1.0);
    for (int d = 0; d < 3; ++d) {
        auto& o = m_open[static_cast<std::size_t>(d)];
        o.define(amrex::convert(m_ba, amrex::IntVect::TheDimensionVector(d)), m_dm, 1, 0);
        o.setVal(1.0);
    }
    for (const auto& g : m_regions) {
        zero_region(g.inner, m_fluid);
        for (auto& o : m_open) { zero_region(g.inner, o); }
    }
    m_active = true;
}

void NestInject::inject (const State& qstar, State& q, bool full)
{
    if (!m_active) { return; }
    const amrex::MultiFab* src[6] = {&qstar.u, &qstar.v, &qstar.w, &qstar.p, &qstar.theta, &qstar.tke};
    const amrex::MultiFab* dst[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    std::array<amrex::Real, 5> size{};
    for (int n = 0; n < 6; ++n) {
        require_layout("inject", *src[n], m_ba, m_dm);
        require_layout("inject", *dst[n], m_ba, m_dm);
        amrex::MultiFab& t = m_target[static_cast<std::size_t>(n)];
        if (src[n]->ixType() != t.ixType() || dst[n]->ixType() != t.ixType()) {
            amrex::Abort("gust: NestInject::inject: the state's fields differ in type");
        }
        for (const auto& g : m_regions) { copy_region(*src[n], g.inner, m_geom.Domain(), t); }
        if (!full) {
            // The target moves a fraction omega of the way from q to q* (q's chi_s is the last target, moved
            // by nothing in between: the stage holds it).
            for (const auto& g : m_regions) { blend_region(*dst[n], m_factor, g.inner, m_geom.Domain(), t); }
        }
        if (n == 3) { continue; }   // p: not a feedback field (no solved row reads it on chi_s)
        // The target's change (q on chi_s is the old target) and its size, over chi_s's data points.
        const std::size_t k = static_cast<std::size_t>(n < 3 ? n : n - 1);
        amrex::MultiFab d(t.boxArray(), m_dm, 1, 0);
        amrex::MultiFab a(t.boxArray(), m_dm, 1, 0);
        amrex::MultiFab m(t.boxArray(), m_dm, 1, 0);
        d.setVal(0.0);
        a.setVal(0.0);
        for (const auto& g : m_regions) {
            difference_region(t, *dst[n], g.inner, m_geom.Domain(), d, m);
            abs_region(t, g.inner, m_geom.Domain(), a);
        }
        m_change.change[k] = global_max(d, 0);
        size[k] = global_max(a, 0);
    }
    const amrex::Real U = amrex::max(size[0], amrex::max(size[1], size[2]));
    m_change.size = {U, U, U, amrex::Real(0), size[4]};
    m_held = true;
    hold(q, 6);
}

void NestInject::hold (State& q, int nfields) const
{
    if (!m_active) { return; }
    if (!m_held) { amrex::Abort("gust: NestInject::hold: no target yet (inject first)"); }
    // u, v, w first: hold(q, 3) is the velocity only.
    amrex::MultiFab* dst[6] = {&q.u, &q.v, &q.w, &q.p, &q.theta, &q.tke};
    for (int n = 0; n < nfields; ++n) {
        require_layout("hold", *dst[n], m_ba, m_dm);
        const amrex::MultiFab& t = m_target[static_cast<std::size_t>(n)];
        if (dst[n]->ixType() != t.ixType()) { amrex::Abort("gust: NestInject::hold: the state's fields differ in type"); }
        for (const auto& g : m_regions) { copy_region(t, g.inner, m_geom.Domain(), *dst[n]); }
    }
}

void NestInject::zero_rows (ResidualFields& R) const
{
    if (!m_active) { return; }
    for (amrex::MultiFab* r : {&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke}) { zero_points(*r); }
}

void NestInject::zero_points (amrex::MultiFab& f) const
{
    if (!m_active) { return; }
    require_layout("zero_points", f, m_ba, m_dm);
    for (const auto& g : m_regions) { zero_region(g.inner, f); }
}

std::vector<amrex::Real> NestInject::differences (const State& q, const State& qstar) const
{
    std::vector<amrex::Real> out(10, 0);
    if (!m_active) { return out; }
    const amrex::MultiFab* a[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    const amrex::MultiFab* b[5] = {&qstar.u, &qstar.v, &qstar.w, &qstar.theta, &qstar.tke};
    for (int n = 0; n < 5; ++n) {
        require_layout("differences", *a[n], m_ba, m_dm);
        require_layout("differences", *b[n], m_ba, m_dm);
        amrex::MultiFab d(a[n]->boxArray(), m_dm, 1, 0);
        amrex::MultiFab m(a[n]->boxArray(), m_dm, 1, 0);
        d.setVal(0.0);
        m.setVal(0.0);
        for (const auto& g : m_regions) { difference_region(*a[n], *b[n], g.inner, m_geom.Domain(), d, m); }
        const amrex::iMultiFab& owner = owner_mask(d, m_geom.periodicity());
        const double num = deterministic_dot(d, 0, d, 0, &owner);
        const double cnt = deterministic_dot(m, 0, m, 0, &owner);
        out[static_cast<std::size_t>(n)] = global_max(d, 0);
        out[static_cast<std::size_t>(5 + n)] = cnt > 0 ? static_cast<amrex::Real>(std::sqrt(num / cnt)) : amrex::Real(0);
    }
    return out;
}

} // namespace gust
