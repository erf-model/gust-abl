#include "NestRelax.H"

#include "ImplicitOperator.H"
#include "Reductions.H"

#include <AMReX_ParallelDescriptor.H>

#include <cmath>
#include <sstream>

namespace gust {

namespace {

int rate_index (Loc loc)
{
    return (loc == Loc::xface) ? 0 : (loc == Loc::yface) ? 1 : (loc == Loc::zface) ? 2 : 3;
}

int target_index (Equation eq)
{
    switch (eq) {
    case Equation::u: return 0;
    case Equation::v: return 1;
    case Equation::w: return 2;
    case Equation::theta: return 3;
    case Equation::tke: return 4;
    }
    return 4;
}

//! The rate of a point of index type t: c = max(c, c0 w(s)) on the points of
//! the footprint fp (cells) converted to t; face[d] says whether the
//! footprint has lateral faces across d.
void footprint_rate (amrex::MultiFab& rate, const amrex::Box& fp, const std::array<bool, 2>& face, amrex::Real c0)
{
    const amrex::IndexType t = rate.ixType();
    const amrex::Box pb = amrex::convert(fp, t);
    const int L0 = fp.smallEnd(0), H0 = fp.bigEnd(0), L1 = fp.smallEnd(1), H1 = fp.bigEnd(1);
    // A cell's distance from the low face is i - L + 1/2, a face point's i - L.
    const amrex::Real ox = t.nodeCentered(0) ? amrex::Real(0) : amrex::Real(0.5);
    const amrex::Real oy = t.nodeCentered(1) ? amrex::Real(0) : amrex::Real(0.5);
    const amrex::Real M = static_cast<amrex::Real>(kDefectMargin);
    // M added to a direction with no faces takes it out of the minimum (w <= 1) with no branch in the kernel.
    const amrex::Real nx = face[0] ? amrex::Real(0) : M;
    const amrex::Real ny = face[1] ? amrex::Real(0) : M;
    const amrex::Real invM = amrex::Real(1) / M;
    for (amrex::MFIter mfi(rate, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & pb;
        if (!bx.ok()) { continue; }
        const auto c = rate.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real sx = amrex::min(static_cast<amrex::Real>(i - L0) + ox,
                                              static_cast<amrex::Real>(H0 + 1 - i) - ox) + nx;
            const amrex::Real sy = amrex::min(static_cast<amrex::Real>(j - L1) + oy,
                                              static_cast<amrex::Real>(H1 + 1 - j) - oy) + ny;
            const amrex::Real w = amrex::min(amrex::Real(1), amrex::max(amrex::Real(0), amrex::min(sx, sy) * invM));
            c(i, j, k) = amrex::max(c(i, j, k), c0 * w);
        });
    }
}

//! The horizontal speed of the cell from its faces' u and v on the footprint's cells, 0 elsewhere.
void footprint_speed (const amrex::MultiFab& u, const amrex::MultiFab& v, const amrex::Box& fp, amrex::MultiFab& spd)
{
    spd.setVal(0.0);
    for (amrex::MFIter mfi(spd, false); mfi.isValid(); ++mfi) {
        const amrex::Box bx = mfi.validbox() & fp;
        if (!bx.ok()) { continue; }
        const auto s = spd.array(mfi);
        const auto uu = u.const_array(mfi);
        const auto vv = v.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            const amrex::Real uc = amrex::Real(0.5) * (uu(i, j, k) + uu(i + 1, j, k));
            const amrex::Real vc = amrex::Real(0.5) * (vv(i, j, k) + vv(i, j + 1, k));
            s(i, j, k) = std::sqrt(uc * uc + vc * vc);
        });
    }
}

void add_source (const amrex::MultiFab& rate, const amrex::MultiFab& q, const amrex::MultiFab& target, amrex::MultiFab& R)
{
    for (amrex::MFIter mfi(R, false); mfi.isValid(); ++mfi) {
        const auto r = R.array(mfi);
        const auto c = rate.const_array(mfi);
        const auto a = q.const_array(mfi);
        const auto p = target.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r(i, j, k) += c(i, j, k) * (a(i, j, k) - p(i, j, k));
        });
    }
}

void add_diagonal (const amrex::MultiFab& rate, amrex::MultiFab& A)
{
    for (amrex::MFIter mfi(A, false); mfi.isValid(); ++mfi) {
        const auto a = A.array(mfi);
        const auto c = rate.const_array(mfi);
        amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k, coef_P) += c(i, j, k);
        });
    }
}

void require_layout (const char* op, const amrex::MultiFab& a, const amrex::MultiFab& ref)
{
    if (a.boxArray() != ref.boxArray() || a.DistributionMap() != ref.DistributionMap()) {
        amrex::Abort(std::string("gust: NestRelax::") + op + ": the field is not on the parent's layout");
    }
}

static_assert(kDefectMargin > 0, "footprint_rate divides by the margin");

} // namespace

std::string check_feedback_factor (amrex::Real factor)
{
    if (std::isfinite(static_cast<double>(factor)) && factor > 0) { return {}; }
    std::ostringstream os;
    os << "nest.feedback_factor = " << factor << ": must be finite and > 0";
    return os.str();
}

NestRelax::NestRelax (const amrex::BoxArray& ba, const amrex::DistributionMapping& dm, const amrex::Geometry& geom,
                      amrex::Real factor)
    : m_geom(geom), m_factor(factor)
{
    const std::string err = check_feedback_factor(factor);
    if (!err.empty()) { amrex::Abort("gust: " + err); }
    const amrex::IntVect types[4] = {amrex::IntVect(1, 0, 0), amrex::IntVect(0, 1, 0), amrex::IntVect(0, 0, 1), amrex::IntVect(0)};
    for (int n = 0; n < 4; ++n) {
        m_rate[static_cast<std::size_t>(n)].define(amrex::convert(ba, types[n]), dm, 1, 0);
        m_rate[static_cast<std::size_t>(n)].setVal(0.0);
    }
    const int tt[5] = {0, 1, 2, 3, 3};
    for (int n = 0; n < 5; ++n) {
        m_target[static_cast<std::size_t>(n)].define(amrex::convert(ba, types[tt[n]]), dm, 1, 0);
        m_target[static_cast<std::size_t>(n)].setVal(0.0);
    }
}

void NestRelax::set (const State& qstar, const std::vector<DefectRegion>& regions)
{
    const amrex::MultiFab* src[5] = {&qstar.u, &qstar.v, &qstar.w, &qstar.theta, &qstar.tke};
    for (int n = 0; n < 5; ++n) {
        amrex::MultiFab& t = m_target[static_cast<std::size_t>(n)];
        require_layout("set", *src[n], t);
        amrex::MultiFab::Copy(t, *src[n], 0, 0, 1, 0);
    }
    for (auto& r : m_rate) { r.setVal(0.0); }
    m_nests.clear();
    const auto dx = m_geom.CellSizeArray();
    const amrex::Real dxh = amrex::min(dx[0], dx[1]);
    amrex::MultiFab spd(m_rate[3].boxArray(), m_rate[3].DistributionMap(), 1, 0);
    for (const auto& g : regions) {
        footprint_speed(qstar.u, qstar.v, g.footprint, spd);
        Nest n;
        n.name = g.name;
        n.speed = deterministic_sum(spd, 0) / static_cast<amrex::Real>(g.footprint.numPts());
        if (!std::isfinite(static_cast<double>(n.speed))) {
            amrex::Abort("gust: nest.feedback = relax: " + g.name + ": the restricted nest's speed is not finite");
        }
        n.rate = m_factor * n.speed / dxh;
        // A direction where the inner box keeps the footprint's whole span has no faces (a periodic nest).
        const std::array<bool, 2> face{g.inner.length(0) != g.footprint.length(0), g.inner.length(1) != g.footprint.length(1)};
        for (auto& r : m_rate) { footprint_rate(r, g.footprint, face, n.rate); }
        m_nests.push_back(n);
    }
    amrex::Gpu::streamSynchronize();
    m_active = true;
}

void NestRelax::add_to_residual (const State& q, ResidualFields& R) const
{
    if (!m_active) { return; }
    const amrex::MultiFab* qf[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    amrex::MultiFab* rf[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    const int ri[5] = {0, 1, 2, 3, 3};
    for (int n = 0; n < 5; ++n) {
        const amrex::MultiFab& t = m_target[static_cast<std::size_t>(n)];
        require_layout("add_to_residual", *rf[n], t);
        require_layout("add_to_residual", *qf[n], t);
        add_source(m_rate[static_cast<std::size_t>(ri[n])], *qf[n], t, *rf[n]);
    }
}

void NestRelax::add_to_implicit (Equation eq, amrex::MultiFab& A) const
{
    if (!m_active) { return; }
    const amrex::MultiFab& t = target(eq);
    require_layout("add_to_implicit", A, t);
    if (A.nComp() != n_implicit_coefs) { amrex::Abort("gust: NestRelax::add_to_implicit: A must hold the implicit coefficients"); }
    const Loc loc[5] = {Loc::xface, Loc::yface, Loc::zface, Loc::cell, Loc::cell};
    add_diagonal(rate(loc[target_index(eq)]), A);
}

const amrex::MultiFab& NestRelax::rate (Loc loc) const { return m_rate[static_cast<std::size_t>(rate_index(loc))]; }

const amrex::MultiFab& NestRelax::target (Equation eq) const { return m_target[static_cast<std::size_t>(target_index(eq))]; }

std::vector<amrex::Real> NestRelax::differences (const State& q) const
{
    std::vector<amrex::Real> out(10, 0);
    if (!m_active) { return out; }
    const amrex::MultiFab* qf[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    const int ri[5] = {0, 1, 2, 3, 3};
    for (int n = 0; n < 5; ++n) {
        const amrex::MultiFab& t = m_target[static_cast<std::size_t>(n)];
        const amrex::MultiFab& c = m_rate[static_cast<std::size_t>(ri[n])];
        require_layout("differences", *qf[n], t);
        amrex::MultiFab d(t.boxArray(), t.DistributionMap(), 1, 0);
        amrex::MultiFab m(t.boxArray(), t.DistributionMap(), 1, 0);
        for (amrex::MFIter mfi(d, false); mfi.isValid(); ++mfi) {
            const auto dd = d.array(mfi);
            const auto mm = m.array(mfi);
            const auto cc = c.const_array(mfi);
            const auto a = qf[n]->const_array(mfi);
            const auto p = t.const_array(mfi);
            amrex::ParallelFor(mfi.validbox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                const bool on = cc(i, j, k) > amrex::Real(0);
                dd(i, j, k) = on ? std::abs(a(i, j, k) - p(i, j, k)) : amrex::Real(0);
                mm(i, j, k) = on ? amrex::Real(1) : amrex::Real(0);
            });
        }
        const amrex::iMultiFab& owner = owner_mask(d, m_geom.periodicity());
        const double num = deterministic_dot(d, 0, d, 0, &owner);
        const double cnt = deterministic_dot(m, 0, m, 0, &owner);
        out[static_cast<std::size_t>(n)] = global_max(d, 0);
        out[static_cast<std::size_t>(5 + n)] = cnt > 0 ? static_cast<amrex::Real>(std::sqrt(num / cnt)) : amrex::Real(0);
    }
    return out;
}

} // namespace gust
