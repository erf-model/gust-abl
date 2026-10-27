#!/usr/bin/env python3
"""The canopy's drag field in a GUST plotfile against its closed form.

    check_canopy.py PLOTFILE PATCHES [--rtol R]

For every column the integral of canopy_drag over the column,
sum_k canopy_drag * detJ * dzeta, must be the sum of cd * lai over the patches
of PATCHES (ERF's format, "type x y height diameter cd lai laimax") whose
disc holds the column centre, with the nearest periodic image in x (the
deck is periodic in x only). That is the leaf area index times the drag
coefficient, whatever the profile and wherever the canopy top falls in a
cell (docs/canopy.rst). Columns outside every patch must hold 0 in every
cell. rtol (default 5e-6) is relative to the largest cd * lai: it covers a
single-precision plotfile.
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402


def read_patches(path):
    out = []
    with open(path) as fh:
        for line in fh:
            tok = line.split("#")[0].split()
            if tok:
                t, x, y, h, d, cd, lai, _ = (float(v) for v in tok)
                out.append((x, y, 0.5 * d, cd * lai))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("patches")
    ap.add_argument("--rtol", type=float, default=5e-6)
    a = ap.parse_args()
    hdr, data = erf_plotfile.read_fields(a.plotfile, ["canopy_drag", "detJ"])
    f, J = data["canopy_drag"], data["detJ"]
    nx, ny, nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]
    dx, dy, dzeta = hdr["dx"]
    lx = hdr["prob_hi"][0] - hdr["prob_lo"][0]
    patches = read_patches(a.patches)
    scale = max(p[3] for p in patches)
    worst, n_in, bad = 0.0, 0, []
    for i in range(nx):
        x = hdr["prob_lo"][0] + (i + 0.5) * dx
        for j in range(ny):
            y = hdr["prob_lo"][1] + (j + 0.5) * dy
            want = 0.0
            for px, py, r, cdlai in patches:
                ex = x - px
                ex -= lx * round(ex / lx)
                if ex * ex + (y - py) ** 2 <= r * r:
                    want += cdlai
            col = [f[i][j][k] for k in range(nz)]
            got = sum(col[k] * J[i][j][k] * dzeta for k in range(nz))
            if not all(math.isfinite(v) for v in col):
                bad.append("column (%d, %d) holds a value that is not finite" % (i, j))
                continue
            if want == 0.0:
                if any(v != 0.0 for v in col):
                    bad.append("column (%d, %d) is outside every patch but holds drag" % (i, j))
                continue
            n_in += 1
            err = abs(got - want) / scale
            worst = max(worst, err)
            if err > a.rtol:
                bad.append("column (%d, %d): integral %.12g, cd lai %.12g" % (i, j, got, want))
    print("check_canopy: %d columns with canopy; largest |integral - cd lai| / %.4g = %.3e (rtol %.1e)"
          % (n_in, scale, worst, a.rtol))
    if n_in == 0:
        bad.append("no column lies in a patch")
    for b in bad[:20]:
        print("  " + b)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
