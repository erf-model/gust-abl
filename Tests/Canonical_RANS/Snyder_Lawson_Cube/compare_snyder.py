#!/usr/bin/env python3
"""Compare a GUST run of the Snyder and Lawson (1994) cube with the measurements.

    compare_snyder.py PLOTFILE [--latest] [--csv OUT.csv] [--pin NAME=VALUE ...]

u is interpolated to each measured point from the cell centres: bilinearly in x and y, along each of the four
columns in z (their own heights, so stretched layers are read as they are): linearly, but in ln((z + z0) / z0)
between the floor's first two cells, where the profile is the log law (z0 = 1 mm, the floor's): linear
interpolation there reads a log profile low (0.07 m/s at z = 0.02 m between centres at 0.01 and 0.03 m). A solid cell's u is 0,
the wall's value. The cube's upwind face is at x = 1.0 m and the domain centre at y = 0.6 m in the decks.

Reported, for each plane (horizontal z = 0.02 m, vertical y = 0) and in all:
- the bias, mean(u_GUST - u_meas), and the rms difference (m/s);
- the hit rate q: the share of points with |u_GUST - u_meas| <= 0.25 |u_meas| or <= D (VDI 3783 Part 9's
  form; here D = 0.1 m/s, about 0.03 of the approach speed at the cube's height);
- the reattachment length behind the cube on the centreline at z = 0.02 m, the zero crossing of u between the
  stations at x = 0.4 and 0.5 m (measured: 1.36 H from the downwind face on the horizontal plane, 1.33 H on the
  vertical one), and u at x = 0.3 m on the
  centreline at z = 0.02 m (measured: -0.80 m/s).
--pin NAME=VALUE fails the run when NAME (rms_all, q_all, rms_horizontal, rms_vertical, reattach) differs from
VALUE by more than 2 % (a regression pin, not a target). --probe X,Y,Z=U fails it when u at (x, y, z) (the
coordinates of the measurements) differs from U by more than 0.01 m/s: a point away from the measurements, such as
the wake 1.5 m behind the cube, where a closure option acts beyond the zone around the cube. The exit code is the
verdict.
"""
import argparse
import csv
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402

X_FACE, Y_MID, H = 1.0, 0.6, 0.2
Z0_FLOOR = 0.001
W_REL, D_ABS = 0.25, 0.1


def read_meas(path):
    rows = []
    with open(path) as fh:
        for r in csv.DictReader(ln for ln in fh if not ln.startswith('#')):
            rows.append((r['plane'], float(r['x']), float(r['y']), float(r['z']), float(r['u'])))
    return rows


class Field:
    def __init__(self, plotfile):
        hdr, d = erf_plotfile.read_fields(plotfile, ["u", "z", "building_beta"])
        self.u, self.z, self.beta = d["u"], d["z"], d["building_beta"]
        self.lo, self.dx = hdr["prob_lo"], hdr["dx"]
        self.n = [hdr["hi"][k] - hdr["lo"][k] + 1 for k in range(3)]

    def column(self, i, j, z):
        zs, us = self.z[i][j], self.u[i][j]
        bs = self.beta[i][j]
        vals = [0.0 if b >= 1 else v for v, b in zip(us, bs)]
        if z <= zs[0]:
            return vals[0]
        for k in range(1, len(zs)):
            if z <= zs[k]:
                if k == 1:
                    def lg(v):
                        return math.log((v + Z0_FLOOR) / Z0_FLOOR)
                    t = (lg(z) - lg(zs[0])) / (lg(zs[1]) - lg(zs[0]))
                else:
                    t = (z - zs[k - 1]) / (zs[k] - zs[k - 1])
                return (1 - t) * vals[k - 1] + t * vals[k]
        return vals[-1]

    def at(self, X, Y, z):
        fx = (X - self.lo[0]) / self.dx[0] - 0.5
        fy = (Y - self.lo[1]) / self.dx[1] - 0.5
        i0, j0 = int(fx // 1), int(fy // 1)
        tx, ty = fx - i0, fy - j0
        out = 0.0
        for di, wx in ((0, 1 - tx), (1, tx)):
            for dj, wy in ((0, 1 - ty), (1, ty)):
                i = min(max(i0 + di, 0), self.n[0] - 1)
                j = min(max(j0 + dj, 0), self.n[1] - 1)
                out += wx * wy * self.column(i, j, z)
        return out


def stats(pairs):
    n = len(pairs)
    bias = sum(g - m for g, m in pairs) / n
    rms = (sum((g - m) ** 2 for g, m in pairs) / n) ** 0.5
    q = sum(1 for g, m in pairs if abs(g - m) <= W_REL * abs(m) or abs(g - m) <= D_ABS) / n
    return n, bias, rms, q


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("--latest", action="store_true", help="plotfile is a prefix: the last plotfile written with it")
    ap.add_argument("--meas", default=os.path.join(HERE, "snyder_lawson_u.csv"))
    ap.add_argument("--csv", help="write each point's measured and computed u here")
    ap.add_argument("--pin", action="append", default=[], help="NAME=VALUE, within 2 %%")
    ap.add_argument("--probe", action="append", default=[], help="X,Y,Z=U: u at that point within 0.01 m/s")
    a = ap.parse_args()
    name = sorted(glob.glob(a.plotfile + "[0-9]*"))[-1] if a.latest else a.plotfile
    f = Field(name)
    meas = read_meas(a.meas)
    res = []
    for plane, x, y, z, um in meas:
        res.append((plane, x, y, z, um, f.at(X_FACE + x, Y_MID + y, z)))
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh, lineterminator="\n")
            w.writerow(["plane", "x", "y", "z", "u_meas", "u_gust"])
            for r in res:
                w.writerow([r[0], r[1], r[2], r[3], f"{r[4]:.3f}", f"{r[5]:.4f}"])
    got = {}
    print(f"{name}: {len(res)} points")
    for label, sel in (("horizontal", lambda r: r[0] == "horizontal"), ("vertical", lambda r: r[0] == "vertical"),
                       ("all", lambda r: True)):
        n, bias, rms, q = stats([(r[5], r[4]) for r in res if sel(r)])
        print(f"  {label:10s} n {n:3d}  bias {bias:+.3f}  rms {rms:.3f} m/s  hit rate q {q:.3f}")
        got[f"rms_{label}"], got[f"q_{label}"] = rms, q
    for label, sel in (("upwind (x < 0)", lambda r: r[1] < 0), ("over and beside (0 <= x <= 0.2)", lambda r: 0 <= r[1] <= H),
                       ("wake (x > 0.2)", lambda r: r[1] > H)):
        n, bias, rms, q = stats([(r[5], r[4]) for r in res if sel(r)])
        print(f"  {label:32s} n {n:3d}  bias {bias:+.3f}  rms {rms:.3f}  q {q:.3f}")
    # The wake on the centreline at z = 0.02 m.
    xs = [0.2 + 0.005 * k for k in range(161)]
    us = [f.at(X_FACE + x, Y_MID, 0.02) for x in xs]
    reattach = None
    for k in range(1, len(xs)):
        if us[k - 1] < 0 <= us[k] and xs[k] > 0.25:
            reattach = xs[k - 1] + (xs[k] - xs[k - 1]) * (-us[k - 1]) / (us[k] - us[k - 1]) - H
            break
    u03 = f.at(X_FACE + 0.3, Y_MID, 0.02)
    print(f"  centreline z = 0.02 m: u at x = 0.3 m {u03:+.3f} m/s (measured -0.80); reattachment "
          + (f"{reattach:.3f} m = {reattach / H:.2f} H" if reattach is not None else "none in 0.8 m")
          + " behind the cube (measured 1.33 to 1.36 H)")
    got["reattach"] = reattach if reattach is not None else float("nan")
    bad = []
    for p in a.pin:
        k, v = p.split("=")
        v = float(v)
        if not (abs(got[k] - v) <= 0.02 * abs(v)):
            bad.append(f"{k} = {got[k]:.4f}, pinned at {v}")
    for p in a.probe:
        pos, v = p.split("=")
        x, y, z = (float(c) for c in pos.split(","))
        u = f.at(X_FACE + x, Y_MID + y, z)
        print(f"  u at x = {x:g}, y = {y:g}, z = {z:g} m: {u:.4f} m/s (pinned at {v})")
        if not abs(u - float(v)) <= 0.01:
            bad.append(f"u at ({x:g}, {y:g}, {z:g}) = {u:.4f}, pinned at {v}")
    for b in bad:
        print("FAIL:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
