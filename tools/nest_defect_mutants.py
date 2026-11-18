#!/usr/bin/env python3
"""
The mutants (one line each, D11 one block) of the two-way nests' restriction and defect source
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
 ("D2", "u and v not divided by the parent face area", F, "o(I, J, k) = ao > 0 ? s * rinv / den : amrex::Real(0);",
  "o(I, J, k) = s * rinv;"),
 ("D3", "Omega summed over one nest face row only", F,
  "                for (int a = 0; a < r; ++a) {\n                    const int i = r * (I - ox) + a, j = r * (J - oy) + b;\n                    s += aa",
  "                for (int a = 0; a < 1; ++a) {\n                    const int i = r * (I - ox) + a, j = r * (J - oy) + b;\n                    s += aa"),
 ("D4", "the pressure level not shifted to the parent's", F, "img.plus(d, 0, 1, 0);", ""),
 ("D5", "the pressure shift with the wrong sign", F, "LinComb(diff, 1, pc, 0, -1, img", "LinComb(diff, -1, pc, 0, 1, img"),
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
 # the order for several nests
 ("D11", "each nest restricted in full before the next (the order before N-114 item 8)", F,
  "        restrict_horizontal(*n.fine, *n.mesh, n.ratio, n.offset, mesh_c, geom_c, coarse);\n    }\n"
  "    // The parent's ghost cells with the new u and v, for its Omega and for w_from_omega.\n"
  "    coarse_ghosts(coarse);\n    for (const NestImage& n : nests) {\n",
  "        restrict_horizontal(*n.fine, *n.mesh, n.ratio, n.offset, mesh_c, geom_c, coarse);\n"
  "        coarse_ghosts(coarse);\n"),
 # the open face areas with buildings (N-114 item 7)
 ("D12", "u: the nest's full face areas in the flux", F, "restrict_faces<0>(mesh_f.flux_area(0), fine.u",
  "restrict_faces<0>(mesh_f.ax, fine.u"),
 ("D13", "u: the parent's full face area as the divisor", F, "ac.ParallelCopy(mesh_c.flux_area(0), 0, 0, 1);",
  "ac.ParallelCopy(mesh_c.ax, 0, 0, 1);"),
 ("D14", "v: the nest's full face areas in the flux", F, "restrict_faces<1>(mesh_f.flux_area(1), fine.v",
  "restrict_faces<1>(mesh_f.ay, fine.v"),
 ("D15", "v: the parent's full face area as the divisor", F, "ac.ParallelCopy(mesh_c.flux_area(1), 0, 0, 1);",
  "ac.ParallelCopy(mesh_c.ay, 0, 0, 1);"),
 ("D16", "Omega: the nest's full face areas in the flux", F, "restrict_zfaces(mesh_f.flux_area(2), om_f",
  "restrict_zfaces(mesh_f.az, om_f"),
 ("D17", "Omega: the parent's full face area as the divisor", F, "ac.ParallelCopy(mesh_c.flux_area(2), 0, 0, 1);",
  "ac.ParallelCopy(mesh_c.az, 0, 0, 1);"),
 ("D18", "u, v: a closed parent face divided by its zero area", F,
  "o(I, J, k) = ao > 0 ? s * rinv / den : amrex::Real(0);", "o(I, J, k) = s * rinv / ao;"),
 ("D19", "Omega: a closed parent face divided by its zero area", F,
  "o(I, J, k) = ao > 0 ? s * r2inv / den : amrex::Real(0);", "o(I, J, k) = s * r2inv / ao;"),
]
EQUIVALENT = set()


if __name__ == "__main__":
    mutants.run(M, "NestDefect.*", EQUIVALENT, __doc__)
