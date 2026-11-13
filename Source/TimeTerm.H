#ifndef GUST_TIME_TERM_H_
#define GUST_TIME_TERM_H_

// The physical time term of unsteady RANS by dual time stepping (backward
// Euler, BDF1). With q^n the solution at the last physical step and dt the
// physical step, every transported equation gains
//
//   R_q  <-  R_q + (q - q^n) / dt,       q = u, v, w, theta, k,
//
// per unit volume like the rest of the residual (Residual.H); continuity has
// no time term. The pseudo-time iteration then drives R_q, time term
// included, to zero, which is one backward-Euler step of the steady
// operator. The term's derivative with respect to q is 1/dt on the row's
// own diagonal, so the implicit rows (implicit_rows: the line solves, the
// momentum predictor's a_P and the SIMPLEC denominator d = a_P - sum a_nb
// built from them) gain 1/dt next to the local pseudo-time term 1/dtau, and
// the Newton core's J v, a difference of residuals, sees it too.
//
// The residual adds the term before the rows that replace their residual
// (the wall-k rows of SurfaceLayer.H with bc.surface_wall_k = dirichlet, the
// Dirichlet rows of the boundary conditions, the fixed rows of buildings),
// and implicit_rows adds 1/dt before the same rows are made identity rows,
// so a replaced row carries no time term in either. Rows that add to the
// residual (damping, the source form of the wall-k rows) keep it. The k rows
// of cells partly inside a building (BuildingTke.H) scale the residual and
// the row by (1 - f), and so keep (1 - f) of the term in both.
//
// q^n is a copy of u, v, w, theta and k with their ghost cells (ngrow_state),
// taken by store(). Only the valid points are read here. q^n is not yet
// written to a checkpoint (that comes with the unsteady outputs), and
// unsteady.dt is not yet read (the stage that steps q^n comes with the time
// driver; until then a deck that sets it stops at the unknown namespace).

#include "ImplicitOperator.H"
#include "Residual.H"
#include "State.H"

#include <AMReX_MultiFab.H>

#include <string>
#include <vector>

namespace gust {

//! The physical time term's options.
struct TimeTermOptions
{
    //! The physical step dt [s] (the input unsteady.dt): finite and > 0, with a finite 1/dt.
    amrex::Real dt{0};
};

//! The options' errors, each naming its input; empty when valid.
std::vector<std::string> check_time_term_options (const TimeTermOptions& opt);

class TimeTerm
{
public:
    //! Aborts, naming unsteady.dt, on invalid options. q^n is allocated on
    //! (ba, dm), the cell-centred BoxArray of the state, and set by store().
    TimeTerm (const TimeTermOptions& opt, const amrex::BoxArray& ba, const amrex::DistributionMapping& dm);

    //! q^n <- q (u, v, w, theta and k, ghost cells included). q must be on the term's layout.
    void store (const State& q);

    //! R_q += (q - q^n)/dt on the valid points of u, v, w, theta and k (not p). q and R
    //! on the term's layout; aborts before the first store().
    void add_to_residual (const State& q, ResidualFields& R) const;

    //! A(coef_P) += 1/dt: A holds the implicit coefficients (ImplicitOperator.H) of eq.
    void add_to_implicit (Equation eq, amrex::MultiFab& A) const;

    [[nodiscard]] amrex::Real dt () const { return m_opt.dt; }
    [[nodiscard]] amrex::Real inverse_dt () const { return m_inv_dt; }
    [[nodiscard]] bool stored () const { return m_stored; }
    //! q^n of one field (Equation::u ... tke).
    [[nodiscard]] const amrex::MultiFab& previous (Equation eq) const;

private:
    TimeTermOptions m_opt;
    amrex::Real m_inv_dt{0};
    amrex::MultiFab m_u, m_v, m_w, m_theta, m_tke;
    bool m_stored{false};
};

} // namespace gust

#endif
