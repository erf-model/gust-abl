#!/usr/bin/env python3
"""Terrain files of the Rampart_2D case: a smooth cliff h(x) = H/2 (1 + tanh((x - XC) / W)), H = 900 m, half a cell past the middle
of a 16 km domain with 200 m cells. Two widths: W = 170 m (the ground rises 476 m between two neighbouring nodes, slope
69 deg) and W = 750 m (119 m between two nodes, slope 31 deg).

    python3 make_terrain.py OUTDIR [DX]  (writes rampart_steep.xyz and rampart_gentle.xyz, and prints the layer ratios)

DX (default 200 m) is the spacing of the terrain nodes in x: 25 for the mesh sweep of README.md, whose meshes of
200, 100, 50 and 25 m all have their nodes on it.
"""
import math
import os
import sys

H, XC, DX, TOP, NZ = 900.0, 8100.0, 200.0, 10000.0, 60
XS = [DX * i for i in range(81)]             # 0 ... 16 km, the mesh nodes
YS = [-200.0, 0.0, 200.0, 400.0, 600.0]      # covers the 400 m periodic width with a margin
WIDTHS = (("rampart_steep.xyz", 170.0), ("rampart_gentle.xyz", 750.0))


def cliff(x, w):
    return 0.5 * H * (1.0 + math.tanh((x - XC) / w))


def ratio(dz0):
    """mesh.stretch_ratio that puts the top of NZ layers starting at dz0 on TOP (bisection)."""
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        lo, hi = (r, hi) if dz0 * (r ** NZ - 1.0) / (r - 1.0) < TOP else (lo, r)
    return 0.5 * (lo + hi)


def main_to(out, dx=DX):
    os.makedirs(out, exist_ok=True)
    xs = [dx * i for i in range(int(round(XS[-1] / dx)) + 1)]
    for name, w in WIDTHS:
        with open(os.path.join(out, name), "w") as f:
            f.write(f"# tanh cliff, {H:g} m high, centre x = {XC:g} m, width W = {w:g} m\n")
            for y in YS:
                for x in xs:
                    f.write(f"{x} {y} {cliff(x, w):.4f}\n")
        step = max(abs(cliff(xs[i + 1], w) - cliff(xs[i], w)) for i in range(len(xs) - 1))
        print(f"{os.path.join(out, name)}: largest step between two nodes {step:.1f} m, "
              f"largest slope {math.degrees(math.atan(0.5 * H / w)):.1f} deg")
    for dz0 in (2.0, 8.0):
        print(f"mesh.dz0 = {dz0:g}: mesh.stretch_ratio = {ratio(dz0):.6f} ({NZ} layers to {TOP:g} m)")


if __name__ == "__main__":
    main_to(sys.argv[1] if len(sys.argv) > 1 else ".", float(sys.argv[2]) if len(sys.argv) > 2 else DX)
