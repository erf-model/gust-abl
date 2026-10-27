#!/usr/bin/env python3
"""
The one-line mutants of the canopy's drag on the momentum equations
(Source/CanopyDrag.cpp, Source/Residual.cpp, Source/Canopy.cpp).

Run against Canopy.* and CanopyDrag.* by the shared runner tools/mutants.py
(two ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/canopy_drag_mutants.py BUILD_DIR [D1 D4 ...]
"""
import mutants

D = "Source/CanopyDrag.cpp"
R = "Source/Residual.cpp"
C = "Source/Canopy.cpp"
M = [
 ("D1", "the sink with the wrong sign", D,
  "r(i, j, k) += f(i, j, k) * std::sqrt(s2) * own;", "r(i, j, k) -= f(i, j, k) * std::sqrt(s2) * own;"),
 ("D2", "x and y faces take one cell's f, not the mean of two", D,
  "o(i, j, k) = amrex::Real(0.5) * (f(i - di, j - dj, k) + f(i, j, k));", "o(i, j, k) = f(i, j, k);"),
 ("D3", "z faces take the cell above only", D,
  "const int kb = (k > klo) ? k - 1 : k;", "const int kb = k;"),
 ("D4", "the derivative diagonal without the u_c^2 / |U| term", D,
  "f(i, j, k) * (s + own * own / s)", "f(i, j, k) * s"),
 ("D5", "|U| without w", D,
  "speed2 = uu * uu + vv * vv + ww * ww;", "speed2 = uu * uu + vv * vv;"),
 ("D6", "v at an x- or z-face not averaged over its four faces", D,
  "(C == 1) ? v(i, j, k) : interp<YF, T>(v, i, j, k);", "v(i, j, k);"),
 ("D7", "picard diagonal with the derivative's term", D,
  "a(i, j, k, comp) += f(i, j, k) * std::sqrt(s2);", "a(i, j, k, comp) += f(i, j, k) * (std::sqrt(s2) + own * own / std::sqrt(s2));"),
 ("D8", "the residual without the drag", R,
  "if (in.canopy != nullptr) { add_canopy_drag(*in.canopy, q.u, q.v, q.w, R.u, R.v, R.w); }", ";"),
 ("D9", "implicit_rows without the drag diagonal", R,
  "        amrex::MultiFab::Add(A, ap, 0, coef_P, 1, 0);\n", "\n"),
 ("D11", "the drag rows on the off-diagonal of the rows", R,
  "amrex::MultiFab::Add(A, ap, 0, coef_P, 1, 0);", "amrex::MultiFab::Add(A, ap, 0, coef_XM, 1, 0);"),
 ("D12", "canopy_rows from the derivative whatever the option", D,
  "    add_diagonal<1>(drag.diagonal(), drag.f(1), u, v, w, rows.ap[1], 0);", "    add_diagonal<1>(CanopyDiagonal::derivative, drag.f(1), u, v, w, rows.ap[1], 0);"),
 ("D10", "a drag_diagonal that is not one of the three accepted", C,
  '        err.emplace_back("canopy.drag_diagonal must be derivative, picard or none (got " + opt.drag_diagonal_name + ")");', ";"),
]

if __name__ == "__main__":
    mutants.run(M, "Canopy*", set(), __doc__)
