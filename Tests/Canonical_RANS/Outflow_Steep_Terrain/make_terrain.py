#!/usr/bin/env python3
"""Terrain files of the Outflow_Steep_Terrain case: the cos^2 ridge of Ishihara and Zhou (40 m high, 100 m half-width,
maximum slope 32 deg) with its crest 50 m beyond the outflow face at x = 2400 m, so that the domain ends on the steepest
part of the rising slope; and the same terrain levelled beyond x = 2380 m (8.2 m high, slope 26 deg just inside), so
that the terrain is steep inside the domain but has no slope normal to the outflow face.

    python3 make_terrain.py OUTDIR      (writes outflow_slope.xyz and outflow_slope_flat.xyz)
"""
import math
import os
import sys

H, HALF, XC, XEDGE = 40.0, 100.0, 2450.0, 2380.0
XS = [2.0 * i for i in range(1201)]          # 0 ... 2400 m
YS = [-16.0, 0.0, 16.0, 32.0, 48.0]          # covers the 32 m periodic width with a margin


def ridge(x):
    return H * math.cos(math.pi * (x - XC) / (2 * HALF)) ** 2 if abs(x - XC) <= HALF else 0.0


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else "."
    os.makedirs(out, exist_ok=True)
    for name, flat in (("outflow_slope.xyz", False), ("outflow_slope_flat.xyz", True)):
        with open(os.path.join(out, name), "w") as f:
            f.write(f"# cos^2 ridge, crest at x = {XC:g} m, outflow at 2400 m"
                    f"{f'; no normal slope beyond x = {XEDGE:g} m' if flat else ''}\n")
            for y in YS:
                for x in XS:
                    f.write(f"{x} {y} {ridge(min(x, XEDGE) if flat else x):.4f}\n")
        print(os.path.join(out, name))


if __name__ == "__main__":
    main()
