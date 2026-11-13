#!/usr/bin/env python3
"""
The one-line mutants of the unsteady outputs and the restart of the time
driver (Source/UnsteadyOutput.cpp, Source/TimeDriver.cpp).

Run against UnsteadyOutput.* and TimeDriver.Resume* by the shared runner
tools/mutants.py (two ranks, the unit-test flags). The table is in
docs/testing.rst.

    python3 tools/unsteady_output_mutants.py BUILD_DIR [UO1 UO4 ...]
"""
import mutants

U = "Source/UnsteadyOutput.cpp"
D = "Source/TimeDriver.cpp"
M = [
 ("UO1", "the weight t1 - t0 (the averaging start ignored)", U,
  "const double w = t1 - std::max(t0, double(m_start));", "const double w = t1 - t0;"),
 ("UO2", "West's a = w/W with W before the add", U,
  "const auto a = static_cast<amrex::Real>(w / m_weight);", "const auto a = static_cast<amrex::Real>(w / (m_weight - w + 1));"),
 ("UO3", "M2 from the old mean only", U,
  "m2(i, j, kk, n) += wt * d * (x[n] - mean(i, j, kk, n));", "m2(i, j, kk, n) += wt * d * d;"),
 ("UO4", "a tie takes the higher cell", U,
  "return {amrex::Math::abs(centre_of(z, i, j, k) - g - h) <= dmin ? k : std::numeric_limits<int>::max()};",
  "return {amrex::Math::abs(centre_of(z, i, j, k) - g - h) <= dmin ? -k : std::numeric_limits<int>::max()};"),
 ("UO5", "h measured from z = 0, not from the ground", U,
  "return {ground, top, amrex::Math::abs(centre_of(z, i, j, k) - g - h)};",
  "return {ground, top, amrex::Math::abs(centre_of(z, i, j, k) - h)};"),
 ("UO6", "a point on the high face not given the last cell", U,
  "+ dom.smallEnd(d), dom.bigEnd(d));", "+ dom.smallEnd(d), dom.bigEnd(d) + 1);"),
 ("UO7", "a restart drops the probe row at the checkpoint's time", U,
  "if (t <= keep_until) {", "if (t < keep_until) {"),
 ("UO8", "a resumed fixed step summed from the checkpoint's time (the anchor dropped)", D,
  "m_anchor_t = same_dt ? anchor_time : time;\n    m_anchor_n = same_dt ? anchor_steps : steps;",
  "m_anchor_t = time;\n    m_anchor_n = steps;"),
 ("UO9", "another step after a restart keeps the old anchor", D,
  "const bool same_dt = m_opt.cfl <= 0 && double(m_opt.dt) == dt_ref;", "const bool same_dt = true;"),
 ("UO10", "a resumed run numbers its steps from 1", D,
  "s.step = static_cast<int>(m_steps_ref + static_cast<long>(m_records.size())) + 1;",
  "s.step = static_cast<int>(m_records.size()) + 1;"),
 ("UO11", "a clipped last step keeps the anchor (the old (n_ref + n) dt)", D,
  "if (last && m_opt.cfl <= 0 && std::abs(m_time - fixed) > 1.0e-3 * double(m_opt.dt)) {", "if (false) {"),
 ("UO12", "every last step moves the anchor (a chain at full steps drifts)", D,
  "if (last && m_opt.cfl <= 0 && std::abs(m_time - fixed) > 1.0e-3 * double(m_opt.dt)) {", "if (last) {"),
 ("UO13", "the record's anchor not written", U,
  "<< \"anchor_time \" << hexfloat(r.anchor_time) << \"\\n\"", "<< \"anchor_time \" << hexfloat(0.0) << \"\\n\""),
 ("UO14", "a kept probe file reported as new (the checkpoint's row written twice)", U,
  "return keep_until >= 0 && exists;", "return false;"),
 ("UO15", "a restart without a probe file reported as kept (no first row)", U,
  "return keep_until >= 0 && exists;", "return keep_until >= 0;"),
]

if __name__ == "__main__":
    mutants.run(M, "UnsteadyOutput.*:TimeDriver.Resume*", set(), __doc__)
