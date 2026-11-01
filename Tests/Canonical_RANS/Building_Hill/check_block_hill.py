#!/usr/bin/env python3
"""The flow around a block on a ridge in a GUST plotfile.

    check_block_hill.py PREFIX

The last plotfile with PREFIX. Horizontal speeds are 0 in every solid cell, the volume flux (u times the layer's thickness,
summed over each x-plane) varies by less than --flux-tol, and the air reverses behind the block. The value in
the roof's wall cell is pinned to this deck's. The exit code is the verdict.
"""
import argparse
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402

ROOF = 3.994


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix")
    ap.add_argument("--flux-tol", type=float, default=0.05)
    ap.add_argument("--roof", type=float, default=ROOF)
    a = ap.parse_args()
    name = sorted(glob.glob(a.prefix + "[0-9]*"))[-1]
    hdr, d = erf_plotfile.read_fields(name, ["u", "v", "w", "detJ", "building_beta"])
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
    roof = u[37][jc][k_roof]
    print(f"lowest u behind the block {behind:.3f} m/s; u in the roof's wall cell {roof:.3f} m/s (pinned at {a.roof})")
    if behind >= 0:
        bad.append("no reversal behind the block")
    if abs(roof - a.roof) > 0.02:
        bad.append("the roof's wall cell moved")
    for b in bad:
        print("FAIL:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
