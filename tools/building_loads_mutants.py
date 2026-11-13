#!/usr/bin/env python3
"""
The one-line mutants of the force on the buildings and the laminar closure
(Source/BuildingLoads.cpp, Source/Residual.cpp, Source/Closure.cpp).

Run against BuildingLoads.* and Closure.Laminar* by the shared runner
tools/mutants.py (two ranks, the unit-test flags). The table is in
docs/testing.rst.

    python3 tools/building_loads_mutants.py BUILD_DIR [BL1 BL4 ...]
"""
import mutants

L = "Source/BuildingLoads.cpp"
R = "Source/Residual.cpp"
C = "Source/Closure.cpp"
M = [
 ("BL1", "faces shared by two boxes counted twice (sum for sum_unique)", L,
  "f.sum_unique(0, false, geom.periodicity())", "f.sum(0, false)"),
 ("BL2", "the buildings' fixed rows kept in the residual of the force", R,
  "if (in.buildings != nullptr && in.building_rows) {", "if (in.buildings != nullptr) {"),
 ("BL3", "the faces of a non-periodic domain face counted", L,
  "const bool row = !periodic && (n == lo || n == hi);", "const bool row = false;"),
 ("BL4", "the face volume from half the mean J", L,
  "? -r(i, j, k) * amrex::Real(0.5) * (J(i - di, j - dj, k - dk) + J(i, j, k)) * cell",
  "? -r(i, j, k) * amrex::Real(0.25) * (J(i - di, j - dj, k - dk) + J(i, j, k)) * cell"),
 ("BL5", "turb.closure = none runs the k-L closure", C,
  "    if (laminar()) {\n        evaluate_laminar(q, f);\n        return;\n    }", "    if (false) {\n    }"),
 ("BL6", "turb.nu_molecular = 0 accepted with turb.closure = none", C,
  "if (!std::isnan(opt.nu_molecular) && !(opt.nu_molecular > 0)) {", "if (false) {"),
 ("BL7", "k-L keys accepted with turb.closure = none", C,
  "            if (!laminar_key(key)) {", "            if (false) {"),
 ("BL8", "k diffuses with alpha_molecular in place of nu_molecular", C,
  "f.gamma_tke.setVal(m_opt.nu_molecular, 0, 1, f.gamma_tke.nGrow());",
  "f.gamma_tke.setVal(m_opt.alpha_molecular, 0, 1, f.gamma_tke.nGrow());"),
 ("BL9", "the face volume from one cell's J", L,
  "? -r(i, j, k) * amrex::Real(0.5) * (J(i - di, j - dj, k - dk) + J(i, j, k)) * cell",
  "? -r(i, j, k) * amrex::Real(0.5) * (J(i, j, k) + J(i, j, k)) * cell"),
]

if __name__ == "__main__":
    mutants.run(M, "BuildingLoads.*:Closure.Laminar*", set(), __doc__)
