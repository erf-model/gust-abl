#!/usr/bin/env python3
"""The distance the closure's length reads beside the cube, from a GUST turbulence plotfile.

    check_cube_length.py PLOTFILE --z0-ground Z0G --z0-wall Z0B
    check_cube_length.py --latest PREFIX --z0-ground Z0G --z0-wall Z0B

With buildings.wall_tke = on the closure's d + z0 is the smaller of the ground's d + z0 and the nearest
wall's d_b + buildings.z0. The plotfile's walldist is the distance the closure is given, to which it
adds the ground's z0: d_b + z0_b - z0_g where a wall is the nearer. On this deck (5 m cells) a wall cell's
fluid centroid is 2.5 m from its wall. Checked in a side's wall cell of the second layer, in the roof's
wall cell, in a side's wall cell of the first layer (2.5 m above the ground and 2.5 m from the side: the
smoother of the two sets the length) and far upstream, where it is the height above the ground. The exit
code is the verdict.
"""
import argparse
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("--latest", action="store_true", help="plotfile is a prefix: the last plotfile written with it")
    ap.add_argument("--z0-ground", type=float, required=True, help="bc.surface_z0 of the run (m)")
    ap.add_argument("--z0-wall", type=float, required=True, help="buildings.z0 of the run (m)")
    ap.add_argument("--tol", type=float, default=1e-5, help="absolute tolerance (m); the plotfile may be single precision")
    a = ap.parse_args()
    name = sorted(glob.glob(a.plotfile + "[0-9]*"))[-1] if a.latest else a.plotfile
    hdr, d = erf_plotfile.read_fields(name, ["walldist"])
    wd = d["walldist"]
    ny = hdr["hi"][1] - hdr["lo"][1] + 1
    jc = ny // 2
    wall = 2.5 + a.z0_wall - a.z0_ground
    first = min(2.5, wall)
    want = [("a side's wall cell, second layer", wd[14][8][1], wall),
            ("the roof's wall cell", wd[14][jc][4], wall),
            ("a side's wall cell, first layer", wd[14][8][0], first),
            ("far upstream, ninth layer", wd[1][jc][8], 42.5)]
    bad = []
    for name, got, ref in want:
        print(f"walldist in {name}: {got:.6f} m (expected {ref:.6f})")
        if abs(got - ref) > a.tol:
            bad.append(name)
    for b in bad:
        print(f"FAIL: walldist in {b} is not the expected one")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
