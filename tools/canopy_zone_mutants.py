#!/usr/bin/env python3
"""
The one-line mutants of the count behind the canopy's inflow-face and
damping-layer warnings (BoundaryConditions::canopy_in_boundary_zones,
Source/BoundaryDamping.cpp).

Run against CanopyTke.* by the shared runner tools/mutants.py (two ranks, the
unit-test flags). The table is in docs/testing.rst.

    python3 tools/canopy_zone_mutants.py BUILD_DIR [Z1 Z3 ...]
"""
import mutants

B = "Source/BoundaryDamping.cpp"
M = [
 ("Z1", "every cell counted, not only canopy cells", B,
  "            if (!(fa(i, j, k) > 0)) { return {0, 0}; }\n", "\n"),
 ("Z2", "the inflow face's second column counted in place of its first", B,
  "(in_xlo && i == x0)", "(in_xlo && i == x0 + 1)"),
 ("Z3", "the sponge left out", B,
  "(rayleigh_to_profile && coef.rayleigh_part(above_ground ? h : zc) > 0) || coef.sponge_part(x, y) > 0;",
  "(rayleigh_to_profile && coef.rayleigh_part(above_ground ? h : zc) > 0);"),
 ("Z4", "the Rayleigh layer left out", B,
  "(rayleigh_to_profile && coef.rayleigh_part(above_ground ? h : zc) > 0) || coef.sponge_part(x, y) > 0;",
  "coef.sponge_part(x, y) > 0;"),
 ("Z6", "the x high face's second-last column counted", B,
  "(in_xhi && i == x1)", "(in_xhi && i == x1 - 1)"),
 ("Z7", "the y high face's second-last row counted", B,
  "(in_yhi && j == y1)", "(in_yhi && j == y1 - 1)"),
 ("Z8", "the Rayleigh layer counted when it damps only w", B,
  "(rayleigh_to_profile && coef.rayleigh_part", "(coef.rayleigh_part"),
 ("Z5", "an outflow face counted as an inflow face", B,
  "== FaceType::inflow;\n    const bool in_xhi", "!= FaceType::periodic;\n    const bool in_xhi"),
]

if __name__ == "__main__":
    mutants.run(M, "CanopyTke.*", set(), __doc__)
