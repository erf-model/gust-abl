#!/usr/bin/env python3
"""
The one-line mutants of the canopy terms of the k-L closure
(Source/CanopyTke.cpp, Source/Closure.cpp, Source/Closure_K.H,
Source/Residual.cpp, Source/Canopy.cpp).

Run against Canopy.*, CanopyDrag.* and CanopyTke.* by the shared runner
tools/mutants.py (two ranks, the unit-test flags). The table is in
docs/testing.rst.

    python3 tools/canopy_tke_mutants.py BUILD_DIR [T1 T4 ...]
"""
import mutants

T = "Source/CanopyTke.cpp"
C = "Source/Closure.cpp"
K = "Source/Closure_K.H"
R = "Source/Residual.cpp"
O = "Source/Canopy.cpp"
M = [
 ("T1", "wake production without beta_p", T,
  "pc(i, j, k) = beta_p * fd * wind * wind * wind;", "pc(i, j, k) = fd * wind * wind * wind;"),
 ("T2", "wake production with |U|^2 in place of |U|^3", T,
  "pc(i, j, k) = beta_p * fd * wind * wind * wind;", "pc(i, j, k) = beta_p * fd * wind * wind;"),
 ("T3", "the sink not on the Patankar diagonal", T,
  "            ds(i, j, k) += c;\n", "\n"),
 ("T4", "the sink added to the exact derivative with the wrong sign", T,
  "de(i, j, k) -= above ? c : amrex::Real(0);", "de(i, j, k) += above ? c : amrex::Real(0);"),
 ("T5", "the residual without the wake production", R,
  "            amrex::MultiFab::Subtract(R.tke, cf->Pc, 0, 0, 1, 0);\n", "\n"),
 ("T6", "the residual without the short-circuit sink", R,
  "            amrex::MultiFab::Add(R.tke, cf->Dc, 0, 0, 1, 0);\n", "\n"),
 ("T7", "the capped length is not alpha / f", C,
  "const amrex::Real L_cap = alpha / (amrex::max(fa(i, j, k), amrex::Real(0)) + std::numeric_limits<amrex::Real>::min());",
  "const amrex::Real L_cap = 2 * alpha / (amrex::max(fa(i, j, k), amrex::Real(0)) + std::numeric_limits<amrex::Real>::min());"),
 ("T8", "the live decisions without the canopy", C,
  "            o(i, j, k) = encode(closure::switches_for(kk(i, j, k), n2(i, j, k), l_g, c, f(i, j, k), alpha));\n        });\n    }\n}\n\n//! The decisions with the floor held",
  "            o(i, j, k) = encode(closure::switches_for(kk(i, j, k), n2(i, j, k), l_g, c));\n        });\n    }\n}\n\n//! The decisions with the floor held"),
 ("T9", "the cap's switch not encoded (not frozen)", C,
  " | (s.canopy_capped ? 32 : 0);", ";"),
 ("T10", "the capped decision taken against alpha / (2 f)", K,
  "s.canopy_capped = (L < L_free);", "s.canopy_capped = (L < amrex::Real(0.5) * L_free);"),
 ("T11", "P_b's derivative without d(cmu'/cmu)/dk (the bug fixed here)", C,
  "(rd.d * nt + ratio * dnut) * n2r(i, j, k)", "(ratio * dnut) * n2r(i, j, k)"),
 ("T12", "canopy.beta_p outside [0, 1] accepted", O,
  "        if (!(t.beta_p >= 0 && t.beta_p <= 1)) {", "        if (false) {"),
 ("T13", "canopy.length_alpha with canopy.tke = off accepted", O,
  '        for (const char* k : {"beta_p", "beta_d", "length_alpha"}) {', '        for (const char* k : {"beta_p", "beta_d"}) {'),
]

if __name__ == "__main__":
    mutants.run(M, "Canopy*", set(), __doc__)
