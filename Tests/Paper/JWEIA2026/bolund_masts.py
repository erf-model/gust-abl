#!/usr/bin/env python3
"""Bolund, mast by mast: speed and turbulent kinetic energy at 2 m and 5 m above the ground at every mast with a
sonic anemometer, against the measurements of the blind comparison (Bechmann et al. 2011, Bolund_Measurements.zip
of windninja-tools).

Both the model and the measurements are normalised by the speed at 5 m at the upstream reference mast (M0 for the
westerly cases 1 and 3, M9 for the easterly case 4): s = S / S_ref and k_n = k / S_ref^2, with S the horizontal
speed. The errors are R_S = mean |s_model - s_obs| and R_k = mean |k_n,model - k_n,obs| over the masts (in %).

    python3 bolund_masts.py --runs RUNS --wntools ~/windninja-tools [--out-table bolund_masts.csv]

RUNS is the --out directory of paper_runs.py (runs M_nosensor_bolund_* and M_sensor_bolund_*).
"""
import argparse
import csv
import glob
import io
import math
import os
import sys
import zipfile

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

ZIP = None   # set from --wntools
CASES = {"270": ("case 1", "M0"), "239": ("case 3", "M0"), "90": ("case 4", "M9")}
RUNS = {"without": "M_nosensor", "with": "M_sensor"}


def measurements(case):
    txt = zipfile.ZipFile(ZIP).read(f"Measurements/Dir_{case}.dat").decode()
    rows = list(csv.reader(io.StringIO(txt), delimiter=" ", skipinitialspace=True))
    head = rows[0]
    out = []
    for r in rows[1:]:
        if not r or not r[0].endswith("S"):        # sonic anemometers only (they carry the TKE)
            continue
        d = dict(zip(head, r))
        ident = d["ID"]
        h = int(ident[3:5])
        if h not in (2, 5):
            continue
        us = float(d["u*[m/s]"])
        out.append(dict(mast=ident[:2], h=h, x=float(d["x[m]"]), y=float(d["y[m]"]),
                        S=float(d["vel/u*"]) * us, k=float(d["tke/u*^2"]) * us * us,
                        u=float(d["u/u*"]) * us, v=float(d["v/u*"]) * us))
    return out


class Model:
    def __init__(self, turb):
        self.h, self.d = ep.read_fields(turb, ["x_velocity", "y_velocity", "KE", "z_phys"])
        lo, hi = self.h["lo"], self.h["hi"]
        self.nx, self.ny, self.nz = hi[0] - lo[0] + 1, hi[1] - lo[1] + 1, hi[2] - lo[2] + 1
        self.x0, self.y0 = self.h["prob_lo"][0], self.h["prob_lo"][1]
        self.dx, self.dy = self.h["dx"][0], self.h["dx"][1]

    def at(self, x, y, h):
        """Horizontal speed and k at height h above the ground, bilinear in x, y and linear in ln(height)."""
        fx = (x - self.x0) / self.dx - 0.5; fy = (y - self.y0) / self.dy - 0.5
        i0, j0 = int(math.floor(fx)), int(math.floor(fy)); tx, ty = fx - i0, fy - j0
        S = K = 0.0
        for (i, j, w) in ((i0, j0, (1 - tx) * (1 - ty)), (i0 + 1, j0, tx * (1 - ty)), (i0, j0 + 1, (1 - tx) * ty), (i0 + 1, j0 + 1, tx * ty)):
            z = [self.d["z_phys"][i][j][k] for k in range(self.nz)]
            g = 1.5 * z[0] - 0.5 * z[1]           # the ground below the first cell centre
            hz = [zz - g for zz in z]
            sp = [math.hypot(self.d["x_velocity"][i][j][k], self.d["y_velocity"][i][j][k]) for k in range(self.nz)]
            ke = [self.d["KE"][i][j][k] for k in range(self.nz)]
            k = next(k for k in range(1, self.nz) if hz[k] >= h)
            t = (math.log(h) - math.log(hz[k - 1])) / (math.log(hz[k]) - math.log(hz[k - 1])) if h > hz[0] else 0.0
            S += w * (sp[k - 1] + t * (sp[k] - sp[k - 1]))
            K += w * (ke[k - 1] + t * (ke[k] - ke[k - 1]))
        return S, K


def lastturb(d):
    p = sorted(x for x in glob.glob(f"{d}/turb[0-9]*"))
    return p[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True, help="the --out directory of paper_runs.py")
    ap.add_argument("--wntools", required=True, help="checkout of github.com/nwagenbrenner/windninja-tools")
    ap.add_argument("--out-table", default="bolund_masts.csv")
    a = ap.parse_args()
    global ZIP
    ZIP = os.path.join(a.wntools, "data", "bolund", "Bolund_Measurements.zip")
    rows, summary = [], []
    for case, (name, refm) in CASES.items():
        obs = measurements(case)
        oref = next(o for o in obs if o["mast"] == refm and o["h"] == 5)
        for lab, run in RUNS.items():
            m = Model(lastturb(f"{os.path.abspath(a.runs)}/runs/{run}_bolund_{case}"))
            Sref, _ = m.at(oref["x"], oref["y"], 5.0)
            err = {2: ([], []), 5: ([], [])}
            for o in obs:
                if o["mast"] == refm:
                    continue
                S, K = m.at(o["x"], o["y"], float(o["h"]))
                s_m, k_m = S / Sref, K / Sref ** 2
                s_o, k_o = o["S"] / oref["S"], o["k"] / oref["S"] ** 2
                err[o["h"]][0].append(abs(s_m - s_o)); err[o["h"]][1].append(abs(k_m - k_o))
                rows.append(dict(case=name, run=lab, mast=o["mast"], h=o["h"], s_obs=round(s_o, 3), s_model=round(s_m, 3),
                                 k_obs=round(k_o, 4), k_model=round(k_m, 4)))
            for h in (2, 5):
                es, ek = err[h]
                summary.append((name, lab, h, len(es), 100 * sum(es) / len(es), 100 * sum(ek) / len(ek)))
    with open(a.out_table, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print("case | run | h (m) | n | R_S (%) | R_k (%)")
    for s in summary:
        print(f"{s[0]} | {s[1]} | {s[2]} | {s[3]} | {s[4]:.1f} | {s[5]:.2f}")


if __name__ == "__main__":
    main()
