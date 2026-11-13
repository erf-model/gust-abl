#!/usr/bin/env python3
"""
The one-line mutants of the time driver of dual time stepping
(Source/TimeDriver.cpp, Source/TimeTerm.cpp, Source/Driver.cpp).

Run against TimeDriver.* by the shared runner tools/mutants.py (two ranks, the
unit-test flags). The table is in docs/testing.rst.

    python3 tools/time_driver_mutants.py BUILD_DIR [TD1 TD4 ...]
"""
import mutants

D = "Source/TimeDriver.cpp"
T = "Source/TimeTerm.cpp"
DR = "Source/Driver.cpp"
M = [
 ("TD1", "q^n stored at the first step only", D,
  "    m_term.store(q);\n", "    if (m_records.empty()) { m_term.store(q); }\n"),
 ("TD2", "the time term keeps the first step's dt", D,
  "    m_term.set_dt(dt, m_opt.cfl > 0 ? \"the step from unsteady.cfl\" : \"unsteady.dt\");\n", "\n"),
 ("TD3", "the last step not shortened", D,
  "    if (last) { dt = static_cast<amrex::Real>(left); }\n", "\n"),
 ("TD4", "a tolerance of 1e-6 on the last step (a sliver step is left)", D,
  "return left <= dt * (1 + 1.0e-3);", "return left <= dt * (1 + 1.0e-6);"),
 ("TD5", "the step from unsteady.cfl without the CFL", D,
  "return cfl * inv_rate.min(0);", "return inv_rate.min(0);"),
 ("TD6", "the inner stop without the floor at steady.residual_tol", D,
  "threshold[i] = amrex::max(drop * std::abs(res[i]), m_tol);", "threshold[i] = drop * std::abs(res[i]);"),
 ("TD7", "the inner stop's reference taken at every iteration", D,
  "        if (s.inner == 1) {\n", "        if (s.inner >= 1) {\n"),
 ("TD8", "one inner iteration more than unsteady.max_inner", D,
  "while (s.inner < m_opt.max_inner) {", "while (s.inner <= m_opt.max_inner) {"),
 ("TD9", "unsteady.inner_drop_orders above 15 accepted", D,
  "opt.inner_drop_orders > 0 && opt.inner_drop_orders <= 15", "opt.inner_drop_orders > 0"),
 ("TD10", "unsteady.dt with unsteady.cfl > 0 accepted", D,
  "if (use_cfl && given(\"dt\")) {", "if (false) {"),
 ("TD11", "the ghost cells not filled before dt", D,
  "    m_fill(q);\n", "\n"),
 ("TD12", "set_dt does not change 1/dt", T,
  "    m_inv_dt = amrex::Real(1) / dt;\n", "\n"),
 ("TD13", "the time a sum of rounded steps for a fixed step", D,
  "(m_opt.cfl > 0 ? m_time + double(dt) : fixed)", "(m_time + double(dt))"),
 ("TD14", "a target CFL over a field at rest accepted", D,
  "if (cfl > 0 && dt >= cfl / (2 * kMinRate)) {", "if (false) {"),
 ("TD15", "a step that is not finite accepted", D,
  "if (!std::isfinite(dt) || !(dt > 0)) {", "if (!(dt > 0)) {"),
 ("TD16", "the CFL ramp of an unsteady run stops at a steady verdict", DR,
  "if (v.status == SteadyStatus::NotSteady || (m_opt.stage == StageKind::unsteady && v.status == SteadyStatus::Steady)) {",
  "if (v.status == SteadyStatus::NotSteady) {"),
]

if __name__ == "__main__":
    mutants.run(M, "TimeDriver*", set(), __doc__)
