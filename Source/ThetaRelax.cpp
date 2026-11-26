#include "ThetaRelax.H"

#include "ImplicitOperator.H"

#include <AMReX.H>
#include <AMReX_GpuLaunch.H>

#include <cmath>
#include <cstdio>
#include <string>

namespace gust {

namespace {

void require_cells (const char* where, const char* name, const amrex::MultiFab& f, const amrex::MultiFab& ref)
{
    if (f.boxArray() != ref.boxArray() || f.DistributionMap() != ref.DistributionMap()) {
        amrex::Abort(std::string("gust: ") + where + ": " + name + " must have the layout of theta's target");
    }
    if (!f.ixType().cellCentered()) { amrex::Abort(std::string("gust: ") + where + ": " + name + " must be cell-centred"); }
}

void require_rate (const char* where, amrex::Real inv_tau)
{
    // isnan first: the ordered comparison traps on NaN.
    if (std::isnan(inv_tau) || !(inv_tau >= 0 && std::isfinite(inv_tau))) {
        amrex::Abort(std::string("gust: ") + where + ": the relaxation rate 1/tau must be finite and >= 0");
    }
}

} // namespace

void add_theta_relaxation (const amrex::MultiFab& theta, const amrex::MultiFab& target, amrex::Real inv_tau,
                           amrex::MultiFab& R)
{
    require_rate("add_theta_relaxation", inv_tau);
    require_cells("add_theta_relaxation", "R", R, target);
    require_cells("add_theta_relaxation", "theta", theta, target);
    for (amrex::MFIter mfi(R, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto r = R.array(mfi);
        const auto t = theta.const_array(mfi);
        const auto p = target.const_array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            r(i, j, k) += (t(i, j, k) - p(i, j, k)) * inv_tau;
        });
    }
}

void add_theta_relaxation_diagonal (amrex::Real inv_tau, amrex::MultiFab& A)
{
    require_rate("add_theta_relaxation_diagonal", inv_tau);
    if (A.nComp() != n_implicit_coefs) {
        amrex::Abort("gust: add_theta_relaxation_diagonal: A must hold the implicit coefficients");
    }
    for (amrex::MFIter mfi(A, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const auto a = A.array(mfi);
        amrex::ParallelFor(mfi.tilebox(), [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            a(i, j, k, coef_P) += inv_tau;
        });
    }
}

ThetaRelax::ThetaRelax (amrex::Real tau) : m_tau(tau)
{
    if (std::isnan(tau) || !(tau > 0 && std::isfinite(tau))) {
        char v[32];
        std::snprintf(v, sizeof(v), "%g", static_cast<double>(tau));
        amrex::Abort(std::string("gust: ThetaRelax: simplec.theta_relax_time = ") + v + ": must be finite and > 0 (s)");
    }
}

void ThetaRelax::set_target (const amrex::MultiFab& theta)
{
    if (!theta.ixType().cellCentered()) { amrex::Abort("gust: ThetaRelax::set_target: theta must be cell-centred"); }
    m_target.define(theta.boxArray(), theta.DistributionMap(), 1, 0);
    amrex::MultiFab::Copy(m_target, theta, 0, 0, 1, 0);
}

void ThetaRelax::add_to_residual (const State& q, ResidualFields& R) const
{
    if (!has_target()) { amrex::Abort("gust: ThetaRelax::add_to_residual: the target theta was never set (set_target)"); }
    add_theta_relaxation(q.theta, m_target, amrex::Real(1) / m_tau, R.theta);
}

void ThetaRelax::add_to_implicit (Equation eq, amrex::MultiFab& A) const
{
    if (eq != Equation::theta) { return; }
    if (!has_target()) { amrex::Abort("gust: ThetaRelax::add_to_implicit: the target theta was never set (set_target)"); }
    require_cells("ThetaRelax::add_to_implicit", "A", A, m_target);
    add_theta_relaxation_diagonal(amrex::Real(1) / m_tau, A);
}

} // namespace gust
