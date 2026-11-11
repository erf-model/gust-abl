#!/usr/bin/env python3
"""The flow around a block on a ridge in a GUST plotfile.

    check_block_hill.py PREFIX

The last plotfile with PREFIX. Horizontal speeds are 0 in every solid cell, the volume flux (u times the layer's thickness,
summed over each x-plane) varies by less than --flux-tol, and the air reverses behind the block. u two cells
above the roof's wall cell is pinned to this deck's (the wall cell itself lies in a thin reversed layer, one cell
deep, whose strength follows the cell's solid fraction), as is k in that wall cell, within 5 %; k in the solid
cells is the closure's floor (buildings.wall_tke = on). The exit code is the verdict.
"""
import argparse
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402

ROOF = 7.423


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix")
    ap.add_argument("--flux-tol", type=float, default=0.05)
    ap.add_argument("--roof", type=float, default=ROOF)
    ap.add_argument("--k-roof", type=float, default=1.015, help="k in the roof's wall cell (m2/s2)")
    ap.add_argument("--k-solid", type=float, default=1.001e-6, help="largest k allowed in a solid cell (m2/s2)")
    a = ap.parse_args()
    name = sorted(glob.glob(a.prefix + "[0-9]*"))[-1]
    hdr, d = erf_plotfile.read_fields(name, ["u", "v", "w", "tke", "detJ", "building_beta"])
    u, v, J, beta = d["u"], d["v"], d["detJ"], d["building_beta"]
    nx, ny, nz = [hdr["hi"][n] - hdr["lo"][n] + 1 for n in range(3)]
    cells = [(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)]
    bad = []
    top = max(abs(f[i][j][k]) for f in (u, v) for i, j, k in cells if beta[i][j][k] >= 1)
    print(f"largest horizontal speed in a solid cell: {top:.3e} m/s")
    if top > 1e-10:
        bad.append("the solid cells move")
    flux = [sum(u[i][j][k] * J[i][j][k] for j in range(ny) for k in range(nz)) for i in range(nx)]
    print(f"plane flux / first plane: min {min(flux) / flux[0]:.4f}, max {max(flux) / flux[0]:.4f}")
    if max(abs(f / flux[0] - 1) for f in flux) > a.flux_tol:
        bad.append("the volume flux changes along x")
    jc = 14
    behind = min(u[i][jc][k] for i in (41, 42, 43) for k in (0, 1))
    k_roof = min(k for k in range(nz) if beta[37][jc][k] < 1)
    # Two cells above the roof's wall cell: above the thin reversed layer on the roof, where the value is steady.
    roof = u[37][jc][k_roof + 2]
    print(f"lowest u behind the block {behind:.3f} m/s; u in the roof's wall cell {u[37][jc][k_roof]:.3f} m/s; "
          f"two cells above it {roof:.3f} m/s (pinned at {a.roof})")
    k_wall = d["tke"][37][jc][k_roof]
    k_solid = max(d["tke"][i][j][k] for i, j, k in cells if beta[i][j][k] >= 1)
    print(f"k in the roof's wall cell {k_wall:.3f} m2/s2 (pinned at {a.k_roof}); largest k in a solid cell {k_solid:.3e} m2/s2")
    if abs(k_wall - a.k_roof) > 0.05 * a.k_roof:
        bad.append("k in the roof's wall cell moved")
    if k_solid > a.k_solid:
        bad.append("k in the solid cells is above its floor")
    if behind >= 0:
        bad.append("no reversal behind the block")
    if abs(roof - a.roof) > 0.02:
        bad.append("u above the roof moved")
    for b in bad:
        print("FAIL:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
