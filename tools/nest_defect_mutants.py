#!/usr/bin/env python3
"""
The one-line mutants of the two-way nests' restriction and defect source
(Source/NestDefect.cpp): the restriction I q_f and the source's region.

Run against the NestDefect.* unit tests by the shared runner tools/mutants.py
(two ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/nest_defect_mutants.py BUILD_DIR [D1 D4 ...]
"""
import mutants

F = "Source/NestDefect.cpp"
M = [
 # the restriction
 ("D1", "cell values a plain mean (no detJ weight)", F, "sjq += jj(i, j, k) * qq(i, j, k);", "sjq += qq(i, j, k);"),
 ("D2", "u and v not divided by the parent face area", F, "o(I, J, k) = s * rinv / cc(I, J, k);",
  "o(I, J, k) = s * rinv;"),
 ("D3", "Omega summed over one nest face row only", F,
  "                for (int a = 0; a < r; ++a) {\n                    const int i = r * (I - ox) + a, j = r * (J - oy) + b;\n                    s += aa",
  "                for (int a = 0; a < 1; ++a) {\n                    const int i = r * (I - ox) + a, j = r * (J - oy) + b;\n                    s += aa"),
 ("D4", "the pressure level not shifted to the parent's", F, "img.plus(d, 0, 1, 0);", ""),
 ("D5", "the pressure shift with the wrong sign", F, "(pc.sum(0) - img.sum(0))", "(img.sum(0) - pc.sum(0))"),
 # the region
 ("D6", "a margin of one parent cell", F, "g.inner.grow(d, -kDefectMargin);", "g.inner.grow(d, -kDefectMargin + 1);"),
 ("D7", "the margin also across a periodic direction", F, "if (!periodic[d]) { g.inner", "if (true) { g.inner"),
 # the source
 ("D8", "the source over the whole footprint", F, "copy_in(*r[n], amrex::convert(g.inner",
  "copy_in(*r[n], amrex::convert(g.footprint"),
 ("D9", "the source added, not subtracted", F, "amrex::MultiFab::Subtract(*r[n], *S[n], 0, 0, 1, 0);",
  "amrex::MultiFab::Add(*r[n], *S[n], 0, 0, 1, 0);"),
 ("D10", "u copied without the periodicity", F, "coarse.u.ParallelCopy(img, 0, 0, 1, 0, 0, geom_c.periodicity());",
  "coarse.u.ParallelCopy(img, 0, 0, 1);"),
]
EQUIVALENT = set()


if __name__ == "__main__":
    mutants.run(M, "NestDefect.*", EQUIVALENT, __doc__)
