#!/usr/bin/env python3
"""The flow around one cube in a GUST plotfile.

    check_cube.py PLOTFILE
    check_cube.py --latest PREFIX

Cell values, 5 m cells: the velocity is 0 in every solid cell; the sum of u over each x-plane is the
first plane's within --flux-tol (cell means, so not exact where a plane cuts the block's faces); the
air reverses behind the block and speeds up over its roof; u in a wall cell of the roof and of a side
is pinned to this deck's (--roof, --side: the defaults are those of buildings.wall_form =
fraction_stress; roof_stress has 3.888 and 1.505). The exit code is the verdict.
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
    ap.add_argument("--flux-tol", type=float, default=0.02)
    ap.add_argument("--roof", type=float, default=3.833, help="u in the roof's wall cell (m/s)")
    ap.add_argument("--side", type=float, default=2.387, help="u in a side's wall cell (m/s)")
    a = ap.parse_args()
    name = sorted(glob.glob(a.plotfile + "[0-9]*"))[-1] if a.latest else a.plotfile
    hdr, d = erf_plotfile.read_fields(name, ["u", "v", "w", "building_beta"])
    u, v, w, beta = d["u"], d["v"], d["w"], d["building_beta"]
    nx, ny, nz = [hdr["hi"][n] - hdr["lo"][n] + 1 for n in range(3)]
    cells = [(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)]
    bad = []
    top = max(abs(f[i][j][k]) for f in (u, v, w) for i, j, k in cells if beta[i][j][k] >= 1)
    print(f"largest speed in a solid cell: {top:.3e} m/s")
    if top > 1e-10:
        bad.append("the solid cells move")
    flux = [sum(u[i][j][k] for j in range(ny) for k in range(nz)) for i in range(nx)]
    print(f"plane flux / first plane: min {min(flux) / flux[0]:.4f}, max {max(flux) / flux[0]:.4f}")
    if max(abs(f / flux[0] - 1) for f in flux) > a.flux_tol:
        bad.append("the volume flux changes along x")
    jc = ny // 2
    behind = u[17][jc][0]
    over = u[14][jc][8]
    far = u[1][jc][8]
    print(f"u behind the cube {behind:.3f}, 20 m over its roof {over:.3f}, upstream at that height {far:.3f} m/s")
    # The roof's wall cell against this deck's value: moves with the wall law and the closed face.
    wall = u[14][jc][4]
    print(f"u in the roof's wall cell {wall:.3f} m/s (pinned at {a.roof})")
    if abs(wall - a.roof) > 0.02:
        bad.append("the roof's wall cell moved")
    # A wall cell of the side at y = 45 m, second layer: moves with the sides' wall law and their closed faces.
    side = u[14][8][1]
    print(f"u in a side's wall cell {side:.3f} m/s (pinned at {a.side})")
    if abs(side - a.side) > 0.02:
        bad.append("the side's wall cell moved")
    if behind >= 0:
        bad.append("no reversal behind the cube")
    if over <= far:
        bad.append("no speed-up over the roof")
    for b in bad:
        print("FAIL:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
