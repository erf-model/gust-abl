#!/usr/bin/env python3
"""The tables of the forested escarpment (Tests/Canonical_RANS/Forested_Escarpment).

    make_case.py OUT_DIR [--forest-at-inflow]

Writes, on 25 m points that cover the domain and two cells beyond it:
  terrain.xyz        z = H/2 (1 + tanh((x - xc)/L)), H = 200 m, L = 250 m, xc = 3000 m
                     (steepest slope 22 deg, at x = xc), constant in y;
  canopy_height.xyz  25 m on the forest, 0 elsewhere;
  canopy_lai.xyz     5 on the forest, 0 elsewhere;
the forest covering x in [2200, 3600] m (the slope and the plateau's edge; the
inflow at x = 0 is open ground). --forest-at-inflow extends it to x = 0 (the
deck the inflow-face warning must flag).
"""
import math
import os
import sys

H, L, XC = 200.0, 250.0, 3000.0
DX = 25.0
X0, X1, Y0, Y1 = -2 * DX, 6000.0 + 2 * DX, -2 * DX, 200.0 + 2 * DX


def main(argv):
    if len(argv) < 2 or argv[1].startswith("-"):
        raise SystemExit(__doc__)
    out = argv[1]
    at_inflow = "--forest-at-inflow" in argv[2:]
    os.makedirs(out, exist_ok=True)
    xs = [X0 + i * DX for i in range(int(round((X1 - X0) / DX)) + 1)]
    ys = [Y0 + j * DX for j in range(int(round((Y1 - Y0) / DX)) + 1)]
    lo = X0 if at_inflow else 2200.0
    with open(os.path.join(out, "terrain.xyz"), "w") as t, \
         open(os.path.join(out, "canopy_height.xyz"), "w") as h, \
         open(os.path.join(out, "canopy_lai.xyz"), "w") as a:
        for y in ys:
            for x in xs:
                z = 0.5 * H * (1 + math.tanh((x - XC) / L))
                forest = lo <= x <= 3600.0
                t.write("%.3f %.3f %.6f\n" % (x, y, z))
                h.write("%.3f %.3f %.3f\n" % (x, y, 25.0 if forest else 0.0))
                a.write("%.3f %.3f %.3f\n" % (x, y, 5.0 if forest else 0.0))


if __name__ == "__main__":
    main(sys.argv)
