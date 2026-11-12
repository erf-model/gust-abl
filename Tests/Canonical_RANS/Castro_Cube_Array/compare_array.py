#!/usr/bin/env python3
"""Compare a GUST run of the Castro et al. (2006) staggered cube array with the measurements.

    compare_array.py PLOTFILE TURBPLOT MEAS.csv [--latest] [--csv OUT.csv] [--pin NAME=VALUE ...]

MEAS.csv has the columns station, quantity, z_over_h, value: quantity u is u / u_r and uw is -<u'w'> / u_r^2 in
units of 1e-3 (as Fig. 3c plots it), at the
stations of Munoz-Esparza et al. (2020, Fig. 2): P0 above the cube's centre, P1 one h behind it, P2 one h beside it
and P3 one h behind and beside (half an h in front of the next cube). README.md says where the values come from.
u_r = 10 m/s, the tunnel's free-stream speed.

In the deck's tile (h = 0.02 m) the cube of P0 is centred at (0.02, 0.06) m, so P0 = (0.02, 0.06), P1 = (0.04, 0.06),
P2 = (0.02, 0.04) and P3 = (0.04, 0.04). Every station lies on a cell face, so each value is the bilinear mean of the
four nearest cell-centre columns, interpolated linearly in z between their own centres. A solid cell's u is 0.

GUST's -<u'w'> is the closure's: nu_t (du/dz + dw/dx), from cell-centre differences of u and w, with nu_t (Kmv)
from the turbulence plotfile (io.turb_plot_file).

Reported for u: bias, rms and the hit rate q (|du| <= 0.25 |u_meas| or <= 0.01, VDI 3783 Part 9's form with
D = 0.01 u_r = 0.1 m/s, our choice) per station and in all; and the rms when GUST's u is scaled to the measured u at P0's highest point (how
Munoz-Esparza et al. set their forcing). For -<u'w'>: rms per station and the peak near the roof (0.8 h < z < 1.3 h).
--latest reads PLOTFILE and TURBPLOT as prefixes: the last plotfiles written with them.
--pin NAME=VALUE fails the run when NAME (rms_u, q_u, rms_uw, peak_uw_P1, peak_uw_P3) differs from VALUE by more than 2 %
(a regression pin, not a target). The exit code is the verdict.
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

H, U_R = 0.02, 10.0
STATIONS = {"P0": (0.02, 0.06), "P1": (0.04, 0.06), "P2": (0.02, 0.04), "P3": (0.04, 0.04)}
W_REL, D_ABS = 0.25, 0.01


def read_meas(path):
    with open(path) as fh:
        # uw is given in units of 1e-3 u_r^2, as plotted.
        return [(r["station"], r["quantity"], float(r["z_over_h"]),
                 float(r["value"]) * (1e-3 if r["quantity"] == "uw" else 1.0))
                for r in csv.DictReader(ln for ln in fh if not ln.startswith("#"))]


class Field:
    def __init__(self, plotfile, turbplot):
        hdr, d = erf_plotfile.read_fields(plotfile, ["u", "w", "z", "building_beta"])
        d["nut"] = erf_plotfile.read_fields(turbplot, ["Kmv"])[1]["Kmv"]
        self.d = d
        self.lo, self.dx = hdr["prob_lo"], hdr["dx"]
        self.n = [hdr["hi"][k] - hdr["lo"][k] + 1 for k in range(3)]

    def wcol(self, i, j):
        return [0.0 if b >= 1 else v for v, b in zip(self.d["w"][i][j], self.d["building_beta"][i][j])]

    def col(self, name, i, j):
        i %= self.n[0]
        j %= self.n[1]
        if name == "u":
            return [0.0 if b >= 1 else v for v, b in zip(self.d["u"][i][j], self.d["building_beta"][i][j])]
        if name == "uw":
            nu, zs = self.d["nut"][i][j], self.d["z"][i][j]
            u = self.col("u", i, j)
            # w of a solid neighbour is 0, as u is. Next to a roof the solid cell's u = 0 sits at its centre, half a
            # cell inside the roof, so du/dz there is too small; no measured point is that close to a roof.
            wp, wm = (self.wcol((i + di) % self.n[0], j) for di in (1, -1))
            out = []
            for k in range(len(zs)):
                if self.d["building_beta"][i][j][k] >= 1:
                    out.append(0.0)
                    continue
                k0, k1 = max(k - 1, 0), min(k + 1, len(zs) - 1)
                dudz = (u[k1] - u[k0]) / (zs[k1] - zs[k0])
                dwdx = (wp[k] - wm[k]) / (2 * self.dx[0])
                out.append(nu[k] * (dudz + dwdx))
            return out
        raise KeyError(name)

    def at(self, name, x, y, z):
        fx = (x - self.lo[0]) / self.dx[0] - 0.5
        fy = (y - self.lo[1]) / self.dx[1] - 0.5
        i0, j0 = int(fx // 1), int(fy // 1)
        tx, ty = fx - i0, fy - j0
        out = 0.0
        for di, wx in ((0, 1 - tx), (1, tx)):
            for dj, wy in ((0, 1 - ty), (1, ty)):
                if wx * wy == 0:
                    continue
                zs = self.d["z"][(i0 + di) % self.n[0]][(j0 + dj) % self.n[1]]
                vs = self.col(name, i0 + di, j0 + dj)
                out += wx * wy * interp(zs, vs, z)
        return out


def interp(zs, vs, z):
    if z <= zs[0]:
        return vs[0]
    for k in range(1, len(zs)):
        if z <= zs[k]:
            t = (z - zs[k - 1]) / (zs[k] - zs[k - 1])
            return (1 - t) * vs[k - 1] + t * vs[k]
    return vs[-1]


def stats(pairs, hit=True):
    n = len(pairs)
    d = [g - m for g, m in pairs]
    bias = sum(d) / n
    rms = (sum(v * v for v in d) / n) ** 0.5
    q = sum(1 for g, m in pairs if abs(g - m) <= max(W_REL * abs(m), D_ABS)) / n if hit else None
    return n, bias, rms, q


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("turbplot")
    ap.add_argument("meas")
    ap.add_argument("--latest", action="store_true")
    ap.add_argument("--csv")
    ap.add_argument("--pin", action="append", default=[])
    a = ap.parse_args()
    pf, tp = a.plotfile, a.turbplot
    if a.latest:
        found = [sorted(glob.glob(pre + "[0-9]*")) for pre in (pf, tp)]
        for pre, names in zip((pf, tp), found):
            if not names:
                print(f"FAIL: no plotfile {pre}*")
                return 1
        pf, tp = found[0][-1], found[1][-1]
        if pf[-5:] != tp[-5:]:
            print(f"FAIL: {pf} and {tp} are from different steps")
            return 1
    f = Field(pf, tp)
    meas = read_meas(a.meas)
    rows = []
    for st, qn, zh, v in meas:
        x, y = STATIONS[st]
        g = f.at("u" if qn == "u" else "uw", x, y, zh * H)
        g = g / U_R if qn == "u" else g / U_R ** 2
        rows.append((st, qn, zh, v, g))
    print(f"{pf}: {len(rows)} points")
    res = {}
    for qn, unit, scale in (("u", "u/u_r", 1.0), ("uw", "-<u'w'>/u_r^2 x 1e3", 1e3)):
        sel = [r for r in rows if r[1] == qn]
        for st in sorted({r[0] for r in sel}):
            n, b, rms, q = stats([(r[4], r[3]) for r in sel if r[0] == st], hit=qn == "u")
            print(f"  {qn:2s} {st}  n {n:3d}  bias {b * scale:+.3f}  rms {rms * scale:.3f} {unit}"
                  + (f"  q {q:.3f}" if q is not None else ""))
        n, b, rms, q = stats([(r[4], r[3]) for r in sel], hit=qn == "u")
        print(f"  {qn:2s} all n {n:3d}  bias {b * scale:+.3f}  rms {rms * scale:.3f} {unit}"
              + (f"  q {q:.3f}" if q is not None else ""))
        res["rms_" + qn] = rms * scale
        if q is not None:
            res["q_" + qn] = q
    top = max((r for r in rows if r[0] == "P0" and r[1] == "u"), key=lambda r: r[2])
    s = top[3] / top[4]
    n, b, rms, q = stats([(r[4] * s, r[3]) for r in rows if r[1] == "u"])
    print(f"  u scaled by {s:.3f} to the measured u at P0, z = {top[2]:.2f} h: rms {rms:.3f}, q {q:.3f}")
    for st in ("P1", "P3"):
        band = [r for r in rows if r[0] == st and r[1] == "uw" and 0.8 < r[2] < 1.3]
        if band:
            zz = [z * H for z in [0.8 + 0.01 * k for k in range(51)]]
            gpk = max(f.at("uw", *STATIONS[st], z) for z in zz) / U_R ** 2 * 1e3
            mpk = max(r[3] for r in band) * 1e3
            print(f"  -<u'w'> peak at {st}, 0.8 h < z < 1.3 h: GUST {gpk:.2f}, measured {mpk:.2f} (x 1e-3 u_r^2)")
            res["peak_uw_" + st] = gpk
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["station", "quantity", "z_over_h", "measured", "gust"])
            w.writerows(rows)
    bad = 0
    for p in a.pin:
        name, val = p.split("=")
        if name not in res:
            print(f"FAIL: unknown pin {name} (known: {', '.join(sorted(res))})")
            bad += 1
            continue
        got, want = res[name], float(val)
        if not math.isfinite(got) or abs(got - want) > 0.02 * abs(want):
            print(f"FAIL: {name} = {got:.4f}, pinned at {want}")
            bad += 1
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
