#!/usr/bin/env python3
"""Write the vector-mean lidar profile of one WINSENT class as ``z u v theta k`` in this case's frame.

Reads the class's 10-minute period files (*_lidar_data.csv: time, height in m NHN, wind_speed, wind_direction) and
averages u and v over the periods at each height kept by at least half of them (as ERF's prepare_winsent.py and
lidar_spread_plot.py do). Heights become metres above the lidar ground (468.8 m NHN). theta and k are placeholders
(--theta, --k): column_inflow.py takes only the direction from this file.
"""
import argparse
import csv
import glob
import math
import os

import numpy as np

LIDAR_GROUND_NHN = 468.8


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", help="the class's lidar folder, e.g. .../lidar_data/high_LAI/unstable")
    p.add_argument("--out", required=True)
    p.add_argument("--wd", type=float, default=293.0, help="the frame's wind direction (deg from north)")
    p.add_argument("--theta", type=float, default=300.0)
    p.add_argument("--k", type=float, default=0.2)
    a = p.parse_args()

    files = sorted(glob.glob(os.path.join(a.folder, "*_lidar_data.csv")))
    acc = {}
    for fn in files:
        with open(fn) as fh:
            for r in csv.DictReader(fh):
                h = round(float(r["height"]) - LIDAR_GROUND_NHN, 2)
                s, d = float(r["wind_speed"]), math.radians(float(r["wind_direction"]))
                acc.setdefault(h, []).append((-s * math.sin(d), -s * math.cos(d)))  # east, north
    az = math.radians(a.wd + 180.0)
    ex = (math.sin(az), math.cos(az))
    ey = (-ex[1], ex[0])
    rows = []
    for h in sorted(acc):
        if len(acc[h]) < 0.5 * len(files):
            continue
        e, n = np.mean(acc[h], axis=0)
        rows.append((h, e * ex[0] + n * ex[1], e * ey[0] + n * ey[1], a.theta, a.k))
    np.savetxt(a.out, np.array(rows), fmt="%10.3f %9.4f %9.4f %8.2f %8.4f",
               header=f"z u v theta k: vector mean of {len(files)} lidar periods ({a.folder}), frame +x toward "
                      f"{(a.wd + 180) % 360:g} deg; theta and k are placeholders")
    print(f"{a.out}: {len(files)} periods, {len(rows)} heights {rows[0][0]:.0f}-{rows[-1][0]:.0f} m")


if __name__ == "__main__":
    main()
