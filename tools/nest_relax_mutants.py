#!/usr/bin/env python3
"""
The one-line mutants of the two-way nests' relaxation source
(Source/NestRelax.cpp, nest.feedback = relax and relax_consistent): the rate,
the target, the implicit diagonal and relax_consistent's added term.

Run against the NestRelax.* and NestRelaxConsistent.* unit tests by the shared runner tools/mutants.py
(two ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/nest_relax_mutants.py BUILD_DIR [X1 X4 C1 ...]
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
 # the boxes outside every footprint, skipped (tracker N-114, item 14)
 ("X11", "only the first footprint's boxes kept", F,
  "if (amrex::convert(fp, vb.ixType()).intersects(vb)) { return true; }",
  "return amrex::convert(fp, vb.ixType()).intersects(vb);"),
 ("X12", "boxes kept by the inner box, not the footprint", F,
  "m_footprints.push_back(g.footprint);", "m_footprints.push_back(g.inner);"),
 ("X13", "the boxes that touch a footprint skipped", F,
  "if (!touches_footprint(fps, mfi.validbox())) { continue; }\n        const auto r = R.array(mfi);",
  "if (touches_footprint(fps, mfi.validbox())) { continue; }\n        const auto r = R.array(mfi);"),
 # relax_consistent's added term s = w R_c(q*) (tracker N-114, item 12)
 ("C1", "the added term with the wrong sign", F, "(a(i, j, k) - p(i, j, k)) - s(i, j, k);",
  "(a(i, j, k) - p(i, j, k)) + s(i, j, k);"),
 ("C2", "the added term left out", F, "(a(i, j, k) - p(i, j, k)) - s(i, j, k);", "(a(i, j, k) - p(i, j, k));"),
 ("C3", "the residual not weighted", F, "? ww(i, j, k) * rr(i, j, k) :", "? rr(i, j, k) :"),
 ("C4", "the residual weighted by w squared", F, "? ww(i, j, k) * rr(i, j, k) :",
  "? ww(i, j, k) * ww(i, j, k) * rr(i, j, k) :"),
 ("C5", "the weight from a footprint one cell wider", F,
  "for (auto& w : wt) { footprint_rate(w, g.footprint, face, amrex::Real(1)); }",
  "for (auto& w : wt) { footprint_rate(w, amrex::grow(g.footprint, amrex::IntVect(AMREX_D_DECL(1, 1, 0))), face, "
  "amrex::Real(1)); }"),
]
EQUIVALENT = set()


if __name__ == "__main__":
    mutants.run(M, "NestRelax.*:NestRelaxConsistent.*", EQUIVALENT, __doc__)
