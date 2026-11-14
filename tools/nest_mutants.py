#!/usr/bin/env python3
"""
The one-line mutants of the nested patch's parent reader
(Source/NestParent.cpp, Source/NestParent.H, the interpolation kernel in
Source/Transfer_K.H) and of the checkpoint Header's z-levels
(Source/StateIO.cpp).

Run against NestParent.* and IO.CheckpointHeader* by the shared runner
tools/mutants.py (two ranks, the unit-test flags). The table is in
docs/testing.rst.

    python3 tools/nest_mutants.py BUILD_DIR [N1 N4 ...]
"""
import mutants

F = "Source/NestParent.cpp"
H = "Source/NestParent.H"
K = "Source/Transfer_K.H"
S = "Source/StateIO.cpp"
M = [
 # interpolation
 ("N1", "u interpolated with the cell-centred x axis (wrong stagger)", F,
  "(m_u, patch.u, m_nx, m_cy, ng)", "(m_u, patch.u, m_cx, m_cy, ng)"),
 ("N2", "v interpolated with the cell-centred y axis (wrong stagger)", F,
  "(m_v, patch.v, m_cx, m_ny, ng)", "(m_v, patch.v, m_cx, m_cy, ng)"),
 ("N3", "w interpolated with the nodal x axis (wrong stagger)", F,
  "(m_w, patch.w, m_cx, m_cy, ng)", "(m_w, patch.w, m_nx, m_cy, ng)"),
 ("N4", "one term of the kernel reads the level above (wrong k)", K,
  "f(i, j, k) = (one - ty) * ((one - tx) * q(I, J, k)", "f(i, j, k) = (one - ty) * ((one - tx) * q(I, J, k + 1)"),
 ("N5", "offset subtracted from the bounds with the wrong sign", K,
  "    a.imin -= offset;", "    a.imin += offset;"),
 ("N6", "the patch's position in the parent ignored", K,
  "a.clo = cdomain.smallEnd(d) + offset;", "a.clo = cdomain.smallEnd(d);"),
 ("N7", "periodic direction clamped as before (extrapolates past i_c = -1)", K,
  "        a.imin = -N;", "        a.imin = -1;"),
 # ghost layers
 ("N8", "no horizontal ghost layers", F,
  "const amrex::IntVect ng(ngrow_state, ngrow_state, 0);", "const amrex::IntVect ng(0, 0, 0);"),
 ("N9", "no z ghost rows", F,
  "        zero_gradient_z(*mf);", "        (void)mf;"),
 ("N10", "z ghost rows below the ground copy the second row (wrong k)", F,
  "a(i, j, k) = a(i, j, klo); });", "a(i, j, k) = a(i, j, klo + 1); });"),
 ("N11", "parent images filled without the periodic images", F,
  "amrex::IntVect(kNestImageGhost, kNestImageGhost, 0), pgeom.periodicity());",
  "amrex::IntVect(kNestImageGhost, kNestImageGhost, 0), amrex::Periodicity::NonPeriodic());"),
 # placement
 ("N12", "margin of one parent cell", H,
  "kNestMargin = 2;", "kNestMargin = 1;"),
 ("N13", "ratio not checked", F,
  "if (patch.n_cell[d] != kNestRatio * n) {", "if (false) {"),
 ("N14", "z-levels not compared", F,
  "if (patch_zlevels[k] != h.zlevels[k]) {", "if (false) {"),
 ("N15", "alignment tolerance of half a cell", F,
  "return std::abs(c - n) <= 1.0e-4;", "return std::abs(c - n) <= 0.5;"),
 ("N16", "a version-1 parent accepted", F,
  "if (h.version < 2 || h.zlevels.empty()) {", "if (false) {"),
 # checkpoint Header
 ("N17", "z-levels written in single-precision decimal", S,
  'for (const amrex::Real z : zl) { os << " " << hexfloat(z); }',
  'for (const amrex::Real z : zl) { os << " " << static_cast<float>(z); }'),
 ("N18", "z-level count not checked", S,
  "if (nz1 != n_cell[2] + 1) {", "if (false) {"),
 ("N19", "a missing Header key names the token read, not the caller's key", S,
  "header_error(key, dir, std::string(\"corrupt Header: expected '\") + want + \"'\");",
  "header_error(tok, dir, std::string(\"corrupt Header: expected '\") + want + \"'\");"),
]

# The call of zero_ground_flux after the interpolation (Source/App.cpp) is a
# regression-level mutant: rt_nest_patch_hill fails without it (ground
# Omega 0.80 m/s against 1e-12 x max |w|); run by hand, see docs/testing.rst.

# Equivalent mutants found while writing the tests, not in the table:
#  shifting the parent image in z (img.shift(..., 1)): ParallelCopy keeps
#  each parent level at its own index, so nothing changes;
#  w through the cell-centred stencil: w is cell-centred in x and y;
#  reading the parent checkpoint's ghost cells: they hold the neighbours'
#  values (CheckpointGhostImagesAreTheNeighboursValues), and the images never
#  reach past a non-periodic parent face (kNestMargin).
EQUIVALENT = set()


if __name__ == "__main__":
    mutants.run(M, "NestParent.*:IO.CheckpointHeader*", EQUIVALENT, __doc__)
