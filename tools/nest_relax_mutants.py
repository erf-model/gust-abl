#!/usr/bin/env python3
"""
The one-line mutants of the two-way nests' relaxation source
(Source/NestRelax.cpp, nest.feedback = relax): the rate, the target and the
implicit diagonal.

Run against the NestRelax.* unit tests by the shared runner tools/mutants.py
(two ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/nest_relax_mutants.py BUILD_DIR [X1 X4 ...]
"""
import mutants

F = "Source/NestRelax.cpp"
M = [
 # the source
 ("X1", "the source with the wrong sign", F, "r(i, j, k) += c(i, j, k) * (a(i, j, k) - p(i, j, k));",
  "r(i, j, k) -= c(i, j, k) * (a(i, j, k) - p(i, j, k));"),
 ("X2", "k left out of the source", F, "for (int n = 0; n < 5; ++n) {\n        const amrex::MultiFab& t = m_target[static_cast<std::size_t>(n)];\n        require_layout(\"add_to_residual\"",
  "for (int n = 0; n < 4; ++n) {\n        const amrex::MultiFab& t = m_target[static_cast<std::size_t>(n)];\n        require_layout(\"add_to_residual\""),
 # the rate
 ("X3", "a margin of three parent cells", F, "const amrex::Real M = static_cast<amrex::Real>(kDefectMargin);",
  "const amrex::Real M = static_cast<amrex::Real>(kDefectMargin + 1);"),
 ("X4", "cells measured from the face like face points", F,
  "const amrex::Real ox = t.nodeCentered(0) ? amrex::Real(0) : amrex::Real(0.5);",
  "const amrex::Real ox = amrex::Real(0);"),
 ("X5", "faces also across a periodic span", F, "const amrex::Real ny = face[1] ? amrex::Real(0) : M;",
  "const amrex::Real ny = amrex::Real(0);"),
 ("X6", "overlapping footprints: the last one, not the larger rate", F,
  "c(i, j, k) = amrex::max(c(i, j, k), c0 * w);", "c(i, j, k) = c0 * w;"),
 ("X7", "the rate floored at 0.1 of its full value", F,
  "amrex::min(amrex::Real(1), amrex::max(amrex::Real(0), amrex::min(sx, sy) * invM));",
  "amrex::min(amrex::Real(1), amrex::max(amrex::Real(0.1), amrex::min(sx, sy) * invM));"),
 ("X8", "U from |u| alone", F, "s(i, j, k) = std::sqrt(uc * uc + vc * vc);", "s(i, j, k) = std::abs(uc);"),
 ("X9", "dx the larger horizontal cell size", F, "const amrex::Real dxh = amrex::min(dx[0], dx[1]);",
  "const amrex::Real dxh = amrex::max(dx[0], dx[1]);"),
 # the implicit diagonal
 ("X10", "no diagonal", F, "a(i, j, k, coef_P) += c(i, j, k);", "a(i, j, k, coef_P) += 0;"),
]
EQUIVALENT = set()


if __name__ == "__main__":
    mutants.run(M, "NestRelax.*", EQUIVALENT, __doc__)
