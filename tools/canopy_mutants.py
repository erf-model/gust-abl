#!/usr/bin/env python3
"""
The one-line mutants of the canopy's drag field
(Source/Canopy.cpp, Source/Canopy_K.H).

Run against Canopy.* by the shared runner tools/mutants.py (two ranks, the
unit-test flags). The table is in docs/testing.rst.

    python3 tools/canopy_mutants.py BUILD_DIR [C1 C4 ...]
"""
import mutants

F = "Source/Canopy.cpp"
K = "Source/Canopy_K.H"
M = [
 # profile (Canopy_K.H)
 ("C1", "lower-branch polynomial coefficient 12/6^3 -> 12/215", K,
  "* t + 12.0 / 216.0)", "* t + 12.0 / 215.0)"),
 ("C2", "upper-branch tail without the erfc term", K,
  "- sqrt_2pie * std::erfc(std::sqrt(0.5 * t)));", ");"),
 ("C3", "spanning cell misses the lower branch's top (E(1))", K,
  "(lm_lower_e(p.c / (1.0 - a)) - lm_lower_e(1.0)) + lm_upper_tail", "(lm_lower_e(p.c / (1.0 - a))) + lm_upper_tail"),
 ("C4", "uniform profile not clipped at the canopy top", K,
  "return eta <= 0.0 ? 0.0 : (eta >= 1.0 ? 1.0 : eta);", "return eta <= 0.0 ? 0.0 : eta;"),
 ("C5", "cell value not divided by the cell's height", K,
  "return cd * lai * fraction(z_b / h, z_t / h) / (z_t - z_b);", "return cd * lai * fraction(z_b / h, z_t / h);"),
 # geometry (Canopy.cpp)
 ("C6", "heights from sea level, not from the local ground", F,
  "zb = face(i, j, k) - g;", "zb = face(i, j, k);"),
 ("C7", "a patch across the periodic seam in x is not imaged", F,
  "if (lx > 0) { ex -= lx * std::round(ex / lx); }", ";"),
 ("C8", "overlapping patches replace instead of add", F,
  "sum += p.lm != 0 ?", "sum = p.lm != 0 ?"),
 ("C9", "ghost columns not imaged", F,
  "    image_ghost_columns(m_drag, geom);\n}\n\nvoid Canopy::restrict_from", "}\n\nvoid Canopy::restrict_from"),
 ("C10", "restriction without the fine detJ", F,
  "a += f(rx * i + c, ry * j + b, k) * J(rx * i + c, ry * j + b, k);", "a += f(rx * i + c, ry * j + b, k);"),
 ("C11", "restriction not divided by the coarse detJ", F,
  "d(i, j, k) /= J(i, j, k);", "d(i, j, k) /= amrex::Real(1);"),
 # inputs and restart
 ("C12", "restart hash of the file name, not its contents", F,
  's += "file " + file_hash(opt.file, "canopy.file") + "\\n";', 's += "file " + opt.file + "\\n";'),
 ("C13", "restart with another canopy accepted", F,
  "if (active && rec != hash) {", "if (false) {"),
 ("C14", "patch height 0 accepted", F,
  "} else if (!(v[3] > 0)) {", "} else if (!(v[3] >= 0)) {"),
 ("C15", "laimax with the uniform profile accepted", F,
  'no_effect("laimax", "with canopy.profile = uniform");', ";"),
 ("C16", "negative map values accepted", F,
  "for (const double v : pts.z) { if (!(v >= 0)) { ++bad; } }", ";"),
 ("C17", "patch wider than a periodic direction accepted", F,
  "if (geom.isPeriodic(d) && p.diameter > static_cast<double>(geom.ProbLength(d))) {", "if (false) {"),
 ("C18", "cells above every patch top skip with a wrong bound (half the tallest)", F,
  "hmax = std::max(hmax, p.height);", "hmax = std::max(hmax, 0.5 * p.height);"),
]

if __name__ == "__main__":
    mutants.run(M, "Canopy.*", set(), __doc__)
