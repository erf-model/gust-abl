#!/usr/bin/env python3
"""The flow in a courtyard of buildings in a GUST plotfile.

    check_courtyard.py PREFIX [--k-corner K] [--k-top K]

The last plotfile with PREFIX. Cell values, 5 m cells: the velocity is 0 in every solid cell and k there is the
closure's floor; the sum of u over each x-plane is the first plane's within --flux-tol; k in an inner corner
cell of the courtyard (two sides) is pinned in the second layer and in the fourth, within 5 %, and the two
inner corners agree (the deck is symmetric in y). The defaults are those of buildings.wall_speed = open_faces;
all_faces has 0.0081 and 0.432. The exit code is the verdict.
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
    ap.add_argument("prefix")
    ap.add_argument("--flux-tol", type=float, default=0.02)
    ap.add_argument("--k-corner", type=float, default=0.0227, help="k in an inner corner cell, second layer (m2/s2)")
    ap.add_argument("--k-top", type=float, default=0.633, help="k in the same corner, fourth layer (m2/s2)")
    a = ap.parse_args()
    name = sorted(glob.glob(a.prefix + "[0-9]*"))[-1]
    hdr, d = erf_plotfile.read_fields(name, ["u", "v", "w", "tke", "building_beta", "building_wall"])
    u, v, w, k, beta, wall = d["u"], d["v"], d["w"], d["tke"], d["building_beta"], d["building_wall"]
    nx, ny, nz = [hdr["hi"][n] - hdr["lo"][n] + 1 for n in range(3)]
    cells = [(i, j, kk) for i in range(nx) for j in range(ny) for kk in range(nz)]
    bad = []
    top = max(abs(f[i][j][kk]) for f in (u, v, w) for i, j, kk in cells if beta[i][j][kk] >= 1)
    k_solid = max(k[i][j][kk] for i, j, kk in cells if beta[i][j][kk] >= 1)
    print(f"largest speed in a solid cell: {top:.3e} m/s; largest k there {k_solid:.3e} m2/s2")
    if top > 1e-10:
        bad.append("the solid cells move")
    if k_solid > 1.001e-6:
        bad.append("k in the solid cells is above its floor")
    flux = [sum(u[i][j][kk] for j in range(ny) for kk in range(nz)) for i in range(nx)]
    print(f"plane flux / first plane: min {min(flux) / flux[0]:.4f}, max {max(flux) / flux[0]:.4f}")
    if max(abs(f / flux[0] - 1) for f in flux) > a.flux_tol:
        bad.append("the volume flux changes along x")
    # The inner corners: x 70..75 m, y 55..60 m (a side at high x and one at low y) and its mirror, y 60..65 m.
    lo, hi = (14, 11), (14, 12)
    if int(wall[lo[0]][lo[1]][1]) != 4 + 8 or int(wall[hi[0]][hi[1]][1]) != 4 + 16:
        bad.append("the pinned cells are not the courtyard's inner corners")
    for layer, pin, what in ((1, a.k_corner, "second"), (3, a.k_top, "fourth")):
        got, mirror = k[lo[0]][lo[1]][layer], k[hi[0]][hi[1]][layer]
        print(f"k in the inner corner, {what} layer: {got:.4f} m2/s2 (pinned at {pin}); its mirror {mirror:.4f}")
        if abs(got - pin) > 0.05 * pin:
            bad.append(f"k in the inner corner moved ({what} layer)")
        if abs(got - mirror) > 0.01 * pin:
            bad.append(f"the two inner corners differ ({what} layer)")
    for b in bad:
        print("FAIL:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
