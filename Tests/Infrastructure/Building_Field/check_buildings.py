#!/usr/bin/env python3
"""The building fractions in a GUST plotfile against the height map.

    check_buildings.py PLOTFILE MAP [--rtol R]

For every column, sum_k building_beta * detJ * dzeta must be the mean over the
column's footprint of the tiles' heights (ERF-format MAP, Source/Buildings.H),
periodic in x. The deck has buildings.small_volfrac = 0, so nothing is clipped.
The mask of two columns is pinned too. rtol (default 1e-5) is relative to the tallest building: it covers a
single-precision plotfile.
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402


def edges(c):
    e = [0.5 * (a + b) for a, b in zip(c, c[1:])]
    return [2 * c[0] - e[0]] + e + [2 * c[-1] - e[-1]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("map")
    ap.add_argument("--rtol", type=float, default=1e-5)
    a = ap.parse_args()
    v = [float(t) for t in open(a.map).read().split()]
    mx, my = int(v[0]), int(v[1])
    xe, ye = edges(v[2:2 + mx]), edges(v[2 + mx:2 + mx + my])
    h = v[2 + mx + my:]
    hdr, data = erf_plotfile.read_fields(a.plotfile, ["building_beta", "building_wall", "detJ"])
    b, J = data["building_beta"], data["detJ"]
    nx, ny, nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]
    dx, dy, dz = hdr["dx"]
    lx = hdr["prob_hi"][0] - hdr["prob_lo"][0]
    worst = 0.0
    for i in range(nx):
        for j in range(ny):
            want = 0.0
            for ti in range(mx):
                for s in (-lx, 0.0, lx):
                    wx = min((i + 1) * dx, xe[ti + 1] + s) - max(i * dx, xe[ti] + s)
                    if wx <= 0:
                        continue
                    for tj in range(my):
                        wy = min((j + 1) * dy, ye[tj + 1]) - max(j * dy, ye[tj])
                        if wy > 0:
                            want += wx * wy * h[ti * my + tj]
            want /= dx * dy
            got = sum(b[i][j][k] * J[i][j][k] * dz for k in range(nz))
            worst = max(worst, abs(got - want))
    scale = max(h)
    print(f"largest column volume error {worst:.3e} m against a height of {scale} m")
    # The wall mask in the plotfile: the 30 m block's tiles reach from x = 190 m, inside cell 9 of
    # 20 m, so that cell is a rim with its wall at high x (bit 4), the ground cell included.
    w = data["building_wall"]
    interior = [int(w[12][10][k]) for k in range(nz)]
    ok = int(w[9][10][0]) == 4 and interior.count(1) == 1 and all(v in (0, 1) for v in interior)
    print(f"mask of rim ground cell {int(w[9][10][0])}, interior column {interior}")
    return 0 if worst <= a.rtol * scale and ok else 1


if __name__ == "__main__":
    sys.exit(main())
