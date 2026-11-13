#!/usr/bin/env python3
"""Drag, lift, shedding frequency and wake length of the laminar cylinder (README.md) against the
published values.

    compare_cylinder.py --run DIR --re 40 [--cells-per-d 10] [--t-avg 90] [--check ...]

DIR holds the run's buildings.force_file (cylinder_force.csv) and its plotfiles (cylinder_plt*).
U = 1 m/s, D = 1 m and the span Ly = 0.1 m (the decks), so
  Cd = 2 Fx / (U^2 D Ly),   Cl = 2 Fz / (U^2 D Ly),   St = f D / U.
Steady runs (one force row after the solve): Cd and Cl of the last row. Unsteady runs: the mean Cd
and the lift amplitude (max - min) / 2 over t >= --t-avg, and f from the mean period between the
upward zero crossings of Cl - mean(Cl) there.
The wake length Lw is measured along the line through the axis (zc = H/2 - d, d = 1/n, a cell face:
the mean of the two cell rows beside it): from the rear of the circle (x = D/2) to the first point behind it where u
(u_mean in an unsteady plotfile with averages) turns positive, linearly between cell centres.
--check compares with the bands given (each "lo,hi") and sets the exit code.
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

LY = 0.1      # the span (m): mesh.prob_hi y
X0 = -6.0     # mesh.prob_lo x (m)

# Published values (README.md, Data and acknowledgements), checked against the papers:
# Dennis & Chang (1970) table 2 (L there is in radii: L/D = L/2), Fornberg (1980) p. 845, Park, Kwon &
# Choi (1998) table 1 (time means and amplitudes), Williamson (1996) eq. for laminar parallel shedding
# St = -3.3265/Re + 0.1816 + 1.6e-4 Re. Fornberg's and Park's wake lengths are read from their figures by
# Li, Peng & Kouser (2022) table 1.
PUBLISHED = {
    20: {"Cd": [("Dennis & Chang 1970", 2.045), ("Fornberg 1980", 2.000), ("Park et al. 1998", 2.01)],
         "Lw/D": [("Dennis & Chang 1970", 0.94), ("Fornberg 1980 (via Li et al. 2022)", 0.91)]},
    40: {"Cd": [("Dennis & Chang 1970", 1.522), ("Fornberg 1980", 1.498), ("Park et al. 1998", 1.51)],
         "Lw/D": [("Dennis & Chang 1970", 2.345), ("Fornberg 1980 (via Li et al. 2022)", 2.24)]},
    100: {"Cd steady": [("Dennis & Chang 1970", 1.056), ("Fornberg 1980", 1.058)],
          "Lw/D steady": [("Dennis & Chang 1970", 6.555)],
          "Cd mean": [("Park et al. 1998", 1.33)],
          "Cl amplitude": [("Park et al. 1998", 0.3321)],
          "Lw/D mean": [("Park et al. 1998 (via Li et al. 2022)", 1.405)],
          "St": [("Williamson 1996", 0.1643)]},
}


def read_forces(path):
    """Rows (t, Fx, Fy, Fz) of a force file."""
    with open(path) as f:
        r = csv.DictReader(f)
        return [(float(x["t"]), float(x["Fx"]), float(x["Fy"]), float(x["Fz"])) for x in r]


def coefficients(rows, t_avg):
    """Cd, Cl (steady: last row) or mean Cd, lift amplitude, St (unsteady)."""
    q = 0.5 * LY  # 1/2 U^2 D Ly
    if len(rows) == 1 or t_avg is None:
        if len(rows) > 1:
            print(f"compare_cylinder: {len(rows)} force rows (an unsteady run?) and no --t-avg: the last row only")
        t, fx, _, fz = rows[-1]
        return {"Cd": fx / q, "Cl": fz / q}
    win = [r for r in rows if r[0] >= t_avg]
    if len(win) < 3:
        sys.exit(f"compare_cylinder: fewer than 3 force rows after t = {t_avg}")
    cd = [r[1] / q for r in win]
    cl = [r[3] / q for r in win]
    t = [r[0] for r in win]
    # Time-weighted means (trapezoid), so an uneven dt does not bias them.
    span = t[-1] - t[0]
    mean = lambda y: sum(0.5 * (y[i] + y[i + 1]) * (t[i + 1] - t[i]) for i in range(len(y) - 1)) / span  # noqa: E731
    cdm, clm = mean(cd), mean(cl)
    ups = [t[i] - (cl[i] - clm) * (t[i + 1] - t[i]) / (cl[i + 1] - cl[i])
           for i in range(len(cl) - 1) if cl[i] - clm < 0 <= cl[i + 1] - clm]
    st = (len(ups) - 1) / (ups[-1] - ups[0]) if len(ups) >= 3 else float("nan")
    return {"Cd mean": cdm, "Cl amplitude": 0.5 * (max(cl) - min(cl)), "St": st, "periods": max(len(ups) - 1, 0)}


def wake_length(plotfile, n, H):
    """Lw/D on the axis row, from u_mean when the plotfile has it, else u. H: the domain height."""
    hdr = erf_plotfile._read_header(plotfile)
    name = "u_mean" if "u_mean" in hdr["names"] else "u"
    hdr, d = erf_plotfile.read_fields(plotfile, [name])
    u = d[name]
    dx = 1.0 / n
    zc = 0.5 * H - 1.0 / n
    k = int(round(zc / dx))   # the row above the axis; k - 1 is the row below
    nx = hdr["hi"][0] - hdr["lo"][0] + 1
    ny = hdr["hi"][1] - hdr["lo"][1] + 1
    line = [sum(u[i][j][kk] for j in range(ny) for kk in (k - 1, k)) / (2 * ny) for i in range(nx)]
    xs = [X0 + (i + 0.5) * dx for i in range(nx)]
    for i in range(nx - 1):
        if xs[i] > 0.5 and line[i] < 0 <= line[i + 1]:
            x = xs[i] - line[i] * (xs[i + 1] - xs[i]) / (line[i + 1] - line[i])
            return x - 0.5, name
    return 0.0, name   # no reversed flow on the axis behind the cylinder


def band(text):
    lo, hi = (float(v) for v in text.split(","))
    return lo, hi


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run", required=True, help="the run directory")
    ap.add_argument("--re", type=int, required=True, choices=sorted(PUBLISHED), help="Reynolds number U D / nu")
    ap.add_argument("--cells-per-d", type=int, default=10, help="D / d of the run")
    ap.add_argument("--force-file", default="cylinder_force.csv")
    ap.add_argument("--plot-prefix", default="cylinder_plt")
    ap.add_argument("--height", type=float, default=12.0, help="the domain height H (m), mesh.prob_hi z")
    ap.add_argument("--t-avg", type=float, default=None, help="unsteady: start of the averaging window (s)")
    ap.add_argument("--check-cd", type=band, help="band for Cd (steady) or the mean Cd (unsteady)")
    ap.add_argument("--check-cl", type=band, help="band for Cl (steady) or the lift amplitude (unsteady)")
    ap.add_argument("--check-lw", type=band, help="band for Lw/D")
    ap.add_argument("--check-st", type=band, help="band for St (unsteady)")
    a = ap.parse_args()

    rows = read_forces(os.path.join(a.run, a.force_file))
    c = coefficients(rows, a.t_avg)
    plots = sorted(glob.glob(os.path.join(a.run, a.plot_prefix + "[0-9]*")))
    lw, field = wake_length(plots[-1], a.cells_per_d, a.height) if plots else (float("nan"), "-")
    print(f"Re = {a.re}, D/d = {a.cells_per_d}, {len(rows)} force rows, plotfile {os.path.basename(plots[-1]) if plots else '-'}")
    for key, val in c.items():
        print(f"  {key:13s} {val:.4f}" if key != "periods" else f"  {key:13s} {val}")
    print(f"  {'Lw/D':13s} {lw:.4f} (from {field})")
    print("published:")
    for key, refs in PUBLISHED[a.re].items():
        print("  " + key + ": " + ", ".join(f"{who} {v}" for who, v in refs))

    bad = []
    unsteady = "Cd mean" in c
    checks = [("Cd", a.check_cd, c["Cd mean"] if unsteady else c["Cd"]),
              ("Cl", a.check_cl, c["Cl amplitude"] if unsteady else c["Cl"]),
              ("Lw/D", a.check_lw, lw),
              ("St", a.check_st, c.get("St", float("nan")))]
    for key, b, val in checks:
        if b is None:
            continue
        ok = math.isfinite(val) and b[0] <= val <= b[1]
        print(f"check {key}: {val:.4f} in [{b[0]}, {b[1]}]: {'yes' if ok else 'NO'}")
        if not ok:
            bad.append(key)
    if bad:
        print("FAIL: " + ", ".join(bad))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
