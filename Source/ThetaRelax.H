#ifndef GUST_THETA_RELAX_H_
#define GUST_THETA_RELAX_H_

// theta relaxed towards a parent run's theta (simplec.theta = relaxed_buoyant).
//
// A nested patch that starts from a parent run (init.first_guess = nest) can
// hold theta at the parent's (simplec.theta = frozen_buoyant) or solve it
// (simplec.theta = solve). This option sits between the two: theta is solved
// with a relaxation source towards the parent's theta theta_p,
//
//     R_theta(q) + (theta - theta_p) / tau = 0,
//
// with
//   R_theta  theta's residual (Residual.H), the steady form -d theta/dt;
//   theta_p  the patch's theta at the stage's first iteration: the parent's
//            theta interpolated by the first guess (fixed for the run);
//   tau      simplec.theta_relax_time (s, finite and > 0).
// Away from the boundary rows, theta - theta_p is O(tau R_theta(theta_p)):
// tau -> 0 holds theta at theta_p (frozen_buoyant), tau -> infinity solves it
// (solve). The buoyancy stays in the w equation, as with both.
//
// Implicit: 1/tau does not depend on theta, so it is the source's exact
// derivative; it is added to the diagonal a_P of theta's implicit rows
// (implicit_rows, Residual.cpp), so a line-sweep step solves
// (1/dtau + A_1 + 1/tau) dtheta = -R and stays a contraction at any dtau. The
// w rows' point-implicit buoyancy (Predictor.H) reads that diagonal too, so
// its predicted theta change shrinks with tau. Both are added right after the
// physical time term and the nests' relaxation, before the rows a boundary
// condition, the wall law or a building replaces (those rows lose the source,
// as they lose the others).
//
// theta_p is not written to the checkpoint: the option is refused with
// io.restart (App.cpp), as frozen_buoyant is.

#include "Residual.H"

#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

namespace gust {

//! R += (theta - target) * inv_tau on R's valid cells (R, theta and target cell-centred with the same layout;
//! inv_tau finite and >= 0). Its own kernel.
void add_theta_relaxation (const amrex::MultiFab& theta, const amrex::MultiFab& target, amrex::Real inv_tau,
                           amrex::MultiFab& R);

//! A's diagonal a_P (component coef_P of the implicit coefficients, ImplicitOperator.H) += inv_tau on A's valid
//! cells. Its own kernel.
void add_theta_relaxation_diagonal (amrex::Real inv_tau, amrex::MultiFab& A);

class ThetaRelax
{
public:
    ThetaRelax () = default;
    //! tau: simplec.theta_relax_time (s); aborts unless finite and > 0.
    explicit ThetaRelax (amrex::Real tau);

    //! theta_p := theta (valid cells; the layout is taken from theta). Called once, before the first residual.
    void set_target (const amrex::MultiFab& theta);
    [[nodiscard]] bool has_target () const { return !m_target.empty(); }
    [[nodiscard]] const amrex::MultiFab& target () const { return m_target; }
    [[nodiscard]] amrex::Real tau () const { return m_tau; }

    //! R.theta += (q.theta - theta_p) / tau. Aborts before set_target.
    void add_to_residual (const State& q, ResidualFields& R) const;
    //! The theta rows' a_P += 1/tau; the other equations' rows are left as they are.
    void add_to_implicit (Equation eq, amrex::MultiFab& A) const;

private:
    amrex::Real m_tau{0};
    amrex::MultiFab m_target;
};

} // namespace gust

#endif
