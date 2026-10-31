#!/usr/bin/env python3
"""Mesh runs of the forest edge (not a CTest; README, "Expected Results"):

    mesh_study.py --gust "LAUNCHER... GUST" --deck inputs_forest_edge [--max-iterations 20000]

Builds its decks from inputs_forest_edge and runs, each in its own folder:

  edge_3m/         the deck's pattern on 3 m columns, LAI 2;
  winsent_coarse/  WINSENT's canopy (15 m, LAI 5, C_d 0.15; 704 m period, 320 m forest) on its
                   32 m columns and its layers (an 8 m first layer, top 2400 m);
  winsent_fine/    the same canopy on 5.5 m columns (319 m forest), a 1 m first layer stretched
                   by 1.08 (top 2330 m).

Prints, for each: iterations, the last r, u_ref, the peak w/u_ref and the adjustment length, and winds
112 m behind the edge. The exit code is 0.
"""
import argparse
import glob
import os
import shlex
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

ZLEVELS = [0.00, 8.00, 16.95, 26.98, 38.19, 50.75, 64.80, 80.52, 98.12, 117.82, 139.87, 164.55, 192.17, 223.09,
           257.69, 296.41, 339.75, 388.27, 442.56, 503.33, 571.35, 647.47, 732.68, 828.04, 934.77, 1054.23,
           1187.93, 1337.58, 1505.07, 1692.53, 1902.34, 2137.17, 2400.00]
CASES = {
    # name: nx, ny, dx, h, lai, cd, L, forest, top, F, mesh lines replaced
    "edge_3m": dict(nx=230, ny=4, dx=3.0, h=18.0, lai=2.0, cd=0.2, L=690.0, forest=306.0, F=1.25e-3,
                    nz=49, mgs="58 4", hi="690. 12. 200.", layers=None),
    "winsent_coarse": dict(nx=22, ny=4, dx=32.0, h=15.0, lai=5.0, cd=0.15, L=704.0, forest=320.0, F=1.0417e-4,
                           nz=32, mgs="11 4", hi="704. 128. 2400.", layers="zlevels.txt"),
    "winsent_fine": dict(nx=128, ny=4, dx=5.5, h=15.0, lai=5.0, cd=0.15, L=704.0, forest=319.0, F=1.0729e-4,
                         nz=68, mgs="64 4", hi="704. 22. 2330.", layers=1.08),
}


def write_case(name, c, base):
    os.makedirs(name, exist_ok=True)
    out = []
    for ln in base.splitlines():
        key = ln.split("=")[0].strip()
        if key == "mesh.n_cell":
            ln = "mesh.n_cell = %d %d %d" % (c["nx"], c["ny"], c["nz"])
        elif key == "mesh.prob_hi":
            ln = "mesh.prob_hi = " + c["hi"]
        elif key == "mesh.max_grid_size":
            ln = "mesh.max_grid_size = " + c["mgs"]
        elif key == "mesh.dz0" and c["layers"] == "zlevels.txt":
            ln = "mesh.zlevels_file = zlevels.txt"
        elif key == "mesh.stretch_ratio" and c["layers"] == "zlevels.txt":
            continue
        elif key == "mesh.stretch_ratio" and c["layers"]:
            ln = "mesh.stretch_ratio = %g" % c["layers"]
        elif key == "canopy.cd":
            ln = "canopy.cd = %g" % c["cd"]
        elif key == "forcing.pressure_gradient":
            ln = "forcing.pressure_gradient = %g 0. 0." % c["F"]
        out.append(ln)
    with open(os.path.join(name, "deck"), "w") as f:
        f.write("\n".join(out) + "\n")
    if c["layers"] == "zlevels.txt":
        with open(os.path.join(name, "zlevels.txt"), "w") as f:
            f.write("\n".join("%.2f" % z for z in ZLEVELS) + "\n")
    dx, dy = c["dx"], (float(c["hi"].split()[1]) / c["ny"])
    with open(os.path.join(name, "h.xyz"), "w") as h, open(os.path.join(name, "lai.xyz"), "w") as a:
        for j in range(-2, c["ny"] + 2):
            for i in range(-4, c["nx"] + 4):
                x, y = (i + 0.5) * dx, (j + 0.5) * dy   # the columns' centres: an exact hit, sharp and 2-D
                forest = (x % c["L"]) < c["forest"]
                h.write("%g %g %g\n" % (x, y, c["h"] if forest else 0.0))
                a.write("%g %g %g\n" % (x, y, c["lai"] if forest else 0.0))


def score(name, c):
    plt = sorted(glob.glob(os.path.join(name, "edge_plt*")))[-1]
    _, a = ep.read_fields(plt, ["u", "w", "z"])
    u, w = a["u"], a["w"]
    nx, nz = len(u), len(u[0][0])
    zc = [a["z"][0][0][k] for k in range(nz)]
    dx = c["L"] / nx
    xs = [(i + 0.5) * dx for i in range(nx)]

    def col(f, i, q):
        for k in range(nz - 1):
            if zc[k] <= q <= zc[k + 1]:
                t = (q - zc[k]) / (zc[k + 1] - zc[k])
                return f[i][0][k] * (1.0 - t) + f[i][0][k + 1] * t
        raise ValueError(q)

    h = c["h"]
    iref = min(range(nx), key=lambda i: abs(xs[i] - (c["L"] - 8.5 * h)))
    uref = col(u, iref, 2.0 * h)
    wt = [col(w, i, h) / uref for i in range(nx)]
    peak = max((i for i in range(nx) if xs[i] < c["forest"]), key=lambda i: wt[i])
    adj = None
    for i in range(peak, nx - 1):
        if wt[i] >= 0.01 > wt[i + 1]:
            adj = xs[i] + dx * (wt[i] - 0.01) / (wt[i] - wt[i + 1])
            break
    i112 = min(range(nx), key=lambda i: abs(xs[i] - 112.0))
    log = open(os.path.join(name, "run.log")).read().splitlines()
    rows = [ln.split() for ln in log if ln[:1] == " " and ln.split()[:1] and ln.split()[0].isdigit()]
    steady = any(ln.startswith("gust: steady after") for ln in log)
    print("%s: %s iterations, r %s, %s; u_ref %.3f; peak w/u_ref %.3f at %.0f m; adjustment %s m (%s h)"
          % (name, rows[-1][0], rows[-1][-1], "verdict" if steady else "no verdict", uref, wt[peak], xs[peak],
             "%.1f" % adj if adj else "none", "%.2f" % (adj / h) if adj else "-"))
    print("  at %.0f m: u/u_ref at h/2 %.3f, h %.3f, 2 h %.3f, 60 m %.3f"
          % (xs[i112], col(u, i112, h / 2) / uref, col(u, i112, h) / uref, col(u, i112, 2 * h) / uref,
             col(u, i112, 60.0) / uref))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gust", required=True)
    ap.add_argument("--deck", required=True)
    ap.add_argument("--max-iterations", type=int, default=20000)
    ap.add_argument("--no-run", action="store_true")
    args = ap.parse_args()
    base = open(args.deck).read()
    procs = {}
    for name, c in CASES.items():
        if not args.no_run:
            write_case(name, c, base)
            with open(os.path.join(name, "run.log"), "w") as log:
                procs[name] = subprocess.Popen(shlex.split(args.gust) + ["deck", "steady.max_iterations=%d" % args.max_iterations],
                                               cwd=name, stdout=log, stderr=subprocess.STDOUT)
    for p in procs.values():
        p.wait()
    for name, c in CASES.items():
        score(name, c)
    return 0


if __name__ == "__main__":
    sys.exit(main())
