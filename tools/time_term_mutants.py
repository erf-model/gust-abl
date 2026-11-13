#!/usr/bin/env python3
"""
The one-line mutants of the physical time term of dual time stepping
(Source/TimeTerm.cpp, Source/Residual.cpp).

Run against TimeTerm.* by the shared runner tools/mutants.py (two ranks, the
unit-test flags). The table is in docs/testing.rst.

    python3 tools/time_term_mutants.py BUILD_DIR [DT1 DT4 ...]
"""
import mutants

T = "Source/TimeTerm.cpp"
R = "Source/Residual.cpp"
M = [
 ("DT1", "the time term with the wrong sign", T,
  "ra(i, j, k) += (xa(i, j, k) - na(i, j, k)) * idt;", "ra(i, j, k) += (na(i, j, k) - xa(i, j, k)) * idt;"),
 ("DT2", "a factor 2 on 1/dt in the residual", T,
  "ra(i, j, k) += (xa(i, j, k) - na(i, j, k)) * idt;", "ra(i, j, k) += (xa(i, j, k) - na(i, j, k)) * idt * amrex::Real(2);"),
 ("DT3", "a factor 2 on 1/dt on the diagonal", T,
  "A.plus(m_inv_dt, coef_P, 1, 0);", "A.plus(amrex::Real(2) * m_inv_dt, coef_P, 1, 0);"),
 ("DT4", "q^n taken from the current iterate", T,
  "add_rate(*x[n], *xn[n], m_inv_dt, *r[n]);", "add_rate(*x[n], *x[n], m_inv_dt, *r[n]);"),
 ("DT5", "store() does not copy q", T,
  "amrex::MultiFab::Copy(*dst[n], *src[n], 0, 0, 1, ngrow_state);", "dst[n]->setVal(0);"),
 ("DT6", "the residual without the time term", R,
  "    if (in.time != nullptr) { in.time->add_to_residual(q, R); }\n", "\n"),
 ("DT7", "implicit_rows without 1/dt", R,
  "    if (in.time != nullptr) { in.time->add_to_implicit(eq, A); }\n", "\n"),
 ("DT8", "the time term added after the dirichlet wall-k rows", R,
  "    if (in.time != nullptr) { in.time->add_to_residual(q, R); }\n"
  "    if (in.surface != nullptr) { apply_wall_k(*ssp, *in.surface, q.tke, R.tke); }\n",
  "    if (in.surface != nullptr) { apply_wall_k(*ssp, *in.surface, q.tke, R.tke); }\n"
  "    if (in.time != nullptr) { in.time->add_to_residual(q, R); }\n"),
 ("DT9", "a denormal unsteady.dt (infinite 1/dt) accepted", T,
  "dt > 0 && std::isfinite(dt) && std::isfinite(amrex::Real(1) / dt)", "dt > 0 && std::isfinite(dt)"),
]

if __name__ == "__main__":
    mutants.run(M, "TimeTerm*", set(), __doc__)
