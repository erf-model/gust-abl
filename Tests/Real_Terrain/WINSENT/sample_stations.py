#!/usr/bin/env python3
"""Wind speed and direction at the WINSENT stations, from a GUST or an ERF plotfile.

Reads the last plotfile of a prefix (or a named one), the case's terrain.txt and stations.txt, and prints, for
MMNW and L140 at the heights that ERF's stations sample, the frame's u and v, the speed and the direction the wind
blows from (degrees from north). Between cell centres u and v are linear in ln(z + z0); across columns bilinear in
x and y (tools/askervein_lines.py, which reads both codes' plotfiles). The frame's +x points toward
--wd + 180 degrees (prepare_winsent.py's rotation).
"""
import argparse
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
import askervein_lines as al  # noqa: E402

HEIGHTS = {"MMNW": [11, 25, 35, 46, 47, 59, 72, 74, 87, 98], "L140": [100, 150, 200, 250, 300, 350, 400]}


def read_stations(path):
    out = {}
    for line in open(path):
        t = line.split()
        if t:
            kv = dict(s.split("=") for s in t[1:])
            out[t[0]] = (float(kv["x"]), float(kv["y"]))
    return out


def uv_at(cols, terrain, x, y, height, z0):
    h = cols.hdr
    lo, dx = h["prob_lo"], h["dx"]
    fx, fy = (x - lo[0]) / dx[0] - 0.5, (y - lo[1]) / dx[1] - 0.5
    i0, j0 = int(math.floor(fx)), int(math.floor(fy))
    a, b = fx - i0, fy - j0
    U = V = 0.0
    for di, dj, w in ((0, 0, (1 - a) * (1 - b)), (1, 0, a * (1 - b)), (0, 1, (1 - a) * b), (1, 1, a * b)):
        i, j = i0 + di, j0 + dj
        u, v, z = cols.column(i, j)
        zg = al.ground(terrain, lo[0] + (i + 0.5) * dx[0], lo[1] + (j + 0.5) * dx[1])
        agl = [zz - zg for zz in z]
        if height < agl[0]:
            f = math.log((height + z0) / z0) / math.log((agl[0] + z0) / z0)
            U += w * f * u[0]
            V += w * f * v[0]
            continue
        k, _ = al.bracket(agl, height)
        t = (math.log(height + z0) - math.log(agl[k] + z0)) / (math.log(agl[k + 1] + z0) - math.log(agl[k] + z0))
        U += w * ((1 - t) * u[k] + t * u[k + 1])
        V += w * ((1 - t) * v[k] + t * v[k + 1])
    return U, V


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plotfile", help="plotfile, or a prefix (its highest step is read)")
    p.add_argument("--case", default=".", help="folder with terrain.txt and stations.txt")
    p.add_argument("--wd", type=float, default=293.0, help="the frame's wind direction, deg from north (293)")
    p.add_argument("--z0", type=float, default=0.1)
    a = p.parse_args()
    cols = al.Columns(a.plotfile)
    terrain = al.read_terrain(os.path.join(a.case, "terrain.txt"))
    st = read_stations(os.path.join(a.case, "stations.txt"))
    az = math.radians(a.wd + 180.0)
    ex = (math.sin(az), math.cos(az))
    ey = (-ex[1], ex[0])
    print(f"# {cols.path} ({cols.kind})")
    print("station,agl,u_frame,v_frame,speed,dir_from")
    for name, hs in HEIGHTS.items():
        x, y = st[name]
        for hgt in hs:
            u, v = uv_at(cols, terrain, x, y, hgt, a.z0)
            east, north = u * ex[0] + v * ey[0], u * ex[1] + v * ey[1]
            d = math.degrees(math.atan2(-east, -north)) % 360.0
            print(f"{name},{hgt},{u:.3f},{v:.3f},{math.hypot(u, v):.3f},{d:.1f}")


if __name__ == "__main__":
    sys.exit(main())
