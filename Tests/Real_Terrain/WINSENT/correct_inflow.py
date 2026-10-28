#!/usr/bin/env python3
"""One pass of the inflow correction: make the run's wind at the valley lidar L140 match the lidar.

The lidar profile measured at L140 is imposed on the inflow face as height above the face's ground, 5.8 km upstream
and 83 m higher than L140, so the wind that reaches L140 is slower than the lidar. This pass compares the run at
L140 with the target (the lidar profile, ERF's inflow.txt) at each height above the ground, and writes a new inflow
file whose u and v at each height are the current ones scaled by target/run speed and turned by the direction
difference. Below the lowest and above the highest compared height the end corrections are held. k and theta are
kept. Repeat (run, correct) until the bias is small.
"""
import argparse
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sample_stations as ss  # noqa: E402


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plotfile", help="the converged run's plotfile, or its prefix")
    p.add_argument("--case", required=True, help="folder with terrain.txt and stations.txt")
    p.add_argument("--target", required=True, help="the lidar profile: ``z u v theta k`` (prepare_gust.py's inflow.txt)")
    p.add_argument("--current", required=True, help="the inflow file the run used")
    p.add_argument("--out", required=True, help="the corrected inflow file to write")
    p.add_argument("--zmin", type=float, default=20.0, help="lowest compared height (m)")
    p.add_argument("--zmax", type=float, default=1500.0, help="highest compared height (m)")
    p.add_argument("--dz", type=float, default=10.0)
    p.add_argument("--z0", type=float, default=0.1)
    p.add_argument("--relax", type=float, default=1.0,
                   help="fraction of the correction applied, in (0, 1]: factor**relax and relax * turn (default 1)")
    a = p.parse_args()
    if not 0.0 < a.relax <= 1.0:
        sys.exit("--relax must be in (0, 1]")

    cols = ss.al.Columns(a.plotfile)
    terrain = ss.al.read_terrain(os.path.join(a.case, "terrain.txt"))
    x, y = ss.read_stations(os.path.join(a.case, "stations.txt"))["L140"]
    tgt = np.loadtxt(a.target)
    cur = np.loadtxt(a.current)
    zs = np.arange(a.zmin, a.zmax + 0.5 * a.dz, a.dz)
    fac, dth = [], []
    for z in zs:
        u, v = ss.uv_at(cols, terrain, x, y, z, a.z0)
        ut, vt = np.interp(z, tgt[:, 0], tgt[:, 1]), np.interp(z, tgt[:, 0], tgt[:, 2])
        fac.append(math.hypot(ut, vt) / max(math.hypot(u, v), 1e-3))
        dth.append(math.atan2(vt, ut) - math.atan2(v, u))
    fac, dth = np.array(fac), np.array(dth)
    f = np.interp(cur[:, 0], zs, fac) ** a.relax  # end values held outside [zmin, zmax]
    t = a.relax * np.interp(cur[:, 0], zs, dth)
    c, s = np.cos(t), np.sin(t)
    new = cur.copy()
    new[:, 1] = f * (c * cur[:, 1] - s * cur[:, 2])
    new[:, 2] = f * (s * cur[:, 1] + c * cur[:, 2])
    np.savetxt(a.out, new, fmt="%10.3f %9.4f %9.4f %8.2f %8.4f",
               header=f"z u v theta k: {a.current} corrected at L140 against {a.target} using {cols.path}")
    for z in (20, 50, 100, 200, 300, 400, 800, 1500):
        i = int(np.argmin(abs(zs - z)))
        print(f"z {zs[i]:6.0f} m: speed factor {fac[i]:.3f}, turn {math.degrees(dth[i]):+5.1f} deg")
    sel = (zs >= 100) & (zs <= 400)
    print(f"L140 bias 100-400 m before this pass: mean factor {fac[sel].mean():.3f}")


if __name__ == "__main__":
    sys.exit(main())
