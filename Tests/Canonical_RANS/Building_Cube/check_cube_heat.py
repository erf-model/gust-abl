#!/usr/bin/env python3
"""The heat balance of the cube with warm walls in a GUST plotfile.

    check_cube_heat.py --latest PREFIX --log LOG --theta-wall 305 [--heat H] [--balance-tol 0.02]

Cell values, 5 m cells, air entering at 300 K. Checked: theta in every solid cell is the walls' temperature;
no fluid cell is colder than the inflow or warmer than the walls (beyond --overshoot); the air in a side's
wall cell and behind the block is warmer than the inflow; the heat the walls give the air, from the log line
"the walls take ... K m^3/s", equals the heat carried through the x-planes downstream of the block, the sum
of u (theta - 300) dy dz, within --balance-tol (the ground and the top exchange none in this deck; measured
0.4 % apart, the planes between 577.7 and 585.2 K m^3/s: the planes nearest the block carry part of the heat by
turbulent diffusion along x, which the sum leaves out); and it is pinned to this deck's value (--heat, within
2 %). The indices and the 25 m^2 of a cell's face are those of inputs_cube_heated. The exit code is the verdict.
"""
import argparse
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("--latest", action="store_true", help="plotfile is a prefix: the last plotfile written with it")
    ap.add_argument("--log", required=True, help="the run's log")
    ap.add_argument("--theta-wall", type=float, required=True, help="buildings.theta_surface (K)")
    ap.add_argument("--theta-in", type=float, default=300.0, help="theta of the inflow (K)")
    ap.add_argument("--heat", type=float, default=None, help="the heat the walls give the air (K m^3/s, > 0), within 2 %%")
    ap.add_argument("--balance-tol", type=float, default=0.02)
    ap.add_argument("--overshoot", type=float, default=1e-3, help="K outside [theta-in, theta-wall] allowed in the air")
    a = ap.parse_args()
    name = sorted(glob.glob(a.plotfile + "[0-9]*"))[-1] if a.latest else a.plotfile
    hdr, d = erf_plotfile.read_fields(name, ["u", "theta", "building_beta"])
    u, th, beta = d["u"], d["theta"], d["building_beta"]
    nx, ny, nz = [hdr["hi"][n] - hdr["lo"][n] + 1 for n in range(3)]
    cells = [(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)]
    bad = []
    solid = [th[i][j][k] for i, j, k in cells if beta[i][j][k] >= 1]
    print(f"theta in the solid cells: {min(solid):.6f} to {max(solid):.6f} K (the walls' {a.theta_wall})")
    if max(abs(t - a.theta_wall) for t in solid) > 1e-4:
        bad.append("the solid cells are not at the walls' temperature")
    air = [th[i][j][k] for i, j, k in cells if beta[i][j][k] < 1]
    print(f"theta in the air: {min(air):.6f} to {max(air):.6f} K")
    if min(air) < a.theta_in - a.overshoot or max(air) > a.theta_wall + a.overshoot:
        bad.append("the air is outside the range of the inflow and the walls")
    jc = ny // 2
    side, wake = th[14][8][1], th[17][jc][0]
    print(f"theta in a side's wall cell {side:.4f} K, behind the block near the ground {wake:.4f} K")
    if not (side > a.theta_in + 0.01 and wake > a.theta_in + 0.01):
        bad.append("the air beside and behind the block is not warmed")
    m = re.search(r"the walls take (\S+) K m\^3/s", open(a.log).read())
    if m is None:
        print("FAIL: the log has no line with the walls' heat")
        return 1
    given = -float(m.group(1))
    # Downstream of the block (it ends at x = 80 m, i = 16), short of the outflow face's first columns.
    planes = [sum(u[i][j][k] * (th[i][j][k] - a.theta_in) for j in range(ny) for k in range(nz)) * 25.0
              for i in range(24, nx - 4)]
    carried = sum(planes) / len(planes)
    print(f"heat the walls give the air {given:.2f} K m^3/s; carried downstream {carried:.2f} "
          f"(planes {min(planes):.2f} to {max(planes):.2f})")
    if abs(carried - given) > a.balance_tol * given:
        bad.append("the heat carried downstream is not the heat the walls give")
    if a.heat is not None:
        print(f"pinned at {a.heat}")
        if abs(given - a.heat) > 0.02 * a.heat:
            bad.append("the walls' heat moved")
    for b in bad:
        print("FAIL:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
