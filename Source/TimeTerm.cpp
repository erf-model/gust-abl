#include "TimeTerm.H"

#include <cmath>
#include <sstream>

namespace gust {

std::vector<std::string> check_time_term_options (const TimeTermOptions& opt)
{
    std::vector<std::string> err;
    const amrex::Real dt = opt.dt;
    if (std::isnan(dt) || !(dt > 0 && std::isfinite(dt) && std::isfinite(amrex::Real(1) / dt))) {
        std::ostringstream os;
        os << "unsteady.dt = " << dt << ": must be finite and > 0, with a finite 1/dt";
        err.push_back(os.str());
    }
    return err;
}

TimeTerm::TimeTerm (const TimeTermOptions& opt, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
    : m_opt(opt)
{
    const auto err = check_time_term_options(opt);
    if (!err.empty()) { amrex::Abort("gust: TimeTerm: " + err.front()); }
    m_inv_dt = amrex::Real(1) / opt.dt;
    using amrex::IntVect;
    m_u.define(amrex::convert(ba, IntVect(1, 0, 0)), dm, 1, ngrow_state);
    m_v.define(amrex::convert(ba, IntVect(0, 1, 0)), dm, 1, ngrow_state);
    m_w.define(amrex::convert(ba, IntVect(0, 0, 1)), dm, 1, ngrow_state);
    m_theta.define(ba, dm, 1, ngrow_state);
    m_tke.define(ba, dm, 1, ngrow_state);
}

void TimeTerm::set_dt (amrex::Real dt, const std::string& what)
{
    TimeTermOptions o;
    o.dt = dt;
    if (!check_time_term_options(o).empty()) {
        amrex::Abort("gust: TimeTerm: " + what + " = " + std::to_string(double(dt)) + ": must be finite and > 0, with a finite 1/dt");
    }
    m_opt.dt = dt;
    m_inv_dt = amrex::Real(1) / dt;
}

namespace {

void require_same_layout (const char* op, const char* what, const amrex::MultiFab& a, const amrex::MultiFab& b)
{
    if (a.boxArray() != b.boxArray() || a.DistributionMap() != b.DistributionMap()) {
        amrex::Abort(std::string("gust: TimeTerm::") + op + ": " + what + " is not on the time term's layout");
    }
}

//! r += (x - xn) * idt on the valid points.
void add_rate (const amrex::MultiFab& x, const amrex::MultiFab& xn, amrex::Real idt, amrex::MultiFab& r)
{
    for (amrex::MFIter mfi(r, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto ra = r.array(mfi);
        const auto xa = x.const_array(mfi);
        const auto na = xn.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            ra(i, j, k) += (xa(i, j, k) - na(i, j, k)) * idt;
        });
    }
}

} // namespace

void TimeTerm::store (const State& q)
{
    const amrex::MultiFab* src[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    amrex::MultiFab* dst[5] = {&m_u, &m_v, &m_w, &m_theta, &m_tke};
    const char* names[5] = {"q.u", "q.v", "q.w", "q.theta", "q.tke"};
    for (int n = 0; n < 5; ++n) {
        require_same_layout("store", names[n], *src[n], *dst[n]);
        if (src[n]->nGrow() < ngrow_state) {
            amrex::Abort(std::string("gust: TimeTerm::store: ") + names[n] + " has fewer than ngrow_state ghost cells");
        }
        amrex::MultiFab::Copy(*dst[n], *src[n], 0, 0, 1, ngrow_state);
    }
    m_stored = true;
}

void TimeTerm::add_to_residual (const State& q, ResidualFields& R) const
{
    if (!m_stored) { amrex::Abort("gust: TimeTerm::add_to_residual: q^n was never stored (TimeTerm::store)"); }
    const amrex::MultiFab* x[5] = {&q.u, &q.v, &q.w, &q.theta, &q.tke};
    const amrex::MultiFab* xn[5] = {&m_u, &m_v, &m_w, &m_theta, &m_tke};
    amrex::MultiFab* r[5] = {&R.u, &R.v, &R.w, &R.theta, &R.tke};
    const char* names[5] = {"u", "v", "w", "theta", "tke"};
    for (int n = 0; n < 5; ++n) {
        require_same_layout("add_to_residual", names[n], *x[n], *xn[n]);
        require_same_layout("add_to_residual", names[n], *r[n], *xn[n]);
        add_rate(*x[n], *xn[n], m_inv_dt, *r[n]);
    }
}

void TimeTerm::add_to_implicit (Equation eq, amrex::MultiFab& A) const
{
    const amrex::MultiFab& ref = previous(eq);
    require_same_layout("add_to_implicit", "A", A, ref);
    if (A.nComp() != n_implicit_coefs) { amrex::Abort("gust: TimeTerm::add_to_implicit: A must hold the implicit coefficients"); }
    A.plus(m_inv_dt, coef_P, 1, 0);
}

const amrex::MultiFab& TimeTerm::previous (Equation eq) const
{
    switch (eq) {
    case Equation::u: return m_u;
    case Equation::v: return m_v;
    case Equation::w: return m_w;
    case Equation::theta: return m_theta;
    case Equation::tke: return m_tke;
    }
    return m_tke;
}

} // namespace gust
