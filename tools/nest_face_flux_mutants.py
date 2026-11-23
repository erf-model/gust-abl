#!/usr/bin/env python3
"""
The one-line mutants of the conservative nest faces (nest.feedback_flux = nest; Source/NestFaceFlux.cpp and
NestDefect::add in Source/NestDefect.cpp).

Run against the NestFaceFlux.* unit tests by the shared runner tools/mutants.py
(two ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/nest_face_flux_mutants.py BUILD_DIR [F1 F4 ...]
"""
import mutants

F = "Source/NestFaceFlux.cpp"
M = [
 ("F1", "the ring sums with the wrong sign", F,
  "    const amrex::Real f = amrex::Real(sgn) * dxinv;\n    const amrex::Real rinv",
  "    const amrex::Real f = -amrex::Real(sgn) * dxinv;\n    const amrex::Real rinv"),
 ("F2", "the transverse trapezoid's end nodes at full weight", F, "? amrex::Real(0.5) * rinv : rinv", "? rinv : rinv"),
 ("F3", "a cell-centred flux's normal average shifted by one nest cell", F, "r * (pD + off - oD) + h - 1",
  "r * (pD + off - oD) + h"),
 ("F4", "the nest's fluxes without the ghost exchange", F, "                Gf.FillBoundary(m_geom->periodicity());\n",
  ""),
 ("F5", "the upwind blend's term dropped", F,
  "+ w * ((au(i, j, k) - ac(i, j, k)) - (fu(i, j, k) - fc(i, j, k))))",
  "+ 0 * w * ((au(i, j, k) - ac(i, j, k)) - (fu(i, j, k) - fc(i, j, k))))"),
 ("F6", "the diffusive flux added, not subtracted", F,
  "c(i, j, k) = ((ac(i, j, k) - ad(i, j, k)) - (fc(i, j, k) - fd(i, j, k))) /",
  "c(i, j, k) = ((ac(i, j, k) + ad(i, j, k)) - (fc(i, j, k) + fd(i, j, k))) /"),
 ("F7", "the ring's face on the wrong side", F,
  "const int off = (sgn > 0 ? 1 : 0) + flo;\n    const amrex::Real f = amrex::Real(sgn) * dxinv;\n    for",
  "const int off = (sgn > 0 ? 0 : 1) + flo;\n    const amrex::Real f = amrex::Real(sgn) * dxinv;\n    for"),
 ("F8", "the correction subtracted from the source", "Source/NestDefect.cpp",
  "amrex::MultiFab::Add(*S[n], *c[n], 0, 0, 1, 0);", "amrex::MultiFab::Subtract(*S[n], *c[n], 0, 0, 1, 0);"),
]
EQUIVALENT = set()


if __name__ == "__main__":
    mutants.run(M, "NestFaceFlux.*", EQUIVALENT, __doc__)
