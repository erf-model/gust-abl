#!/usr/bin/env python3
"""Profiles of U and k at the inflow, mid-domain and outflow of the empty-domain ABL runs, against the inflow law.

    python3 plot_hhabl.py <reader dir> <out.png> <label>=<turb plotfile> ...   [--ustar U --z0 Z --cmu0 C --kappa K]
"""
import argparse, math, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument("readerdir"); ap.add_argument("out"); ap.add_argument("runs", nargs="+")
ap.add_argument("--ustar", type=float, default=0.6278); ap.add_argument("--z0", type=float, default=0.03)
ap.add_argument("--cmu0", type=float, default=0.5562); ap.add_argument("--kappa", type=float, default=0.41)
ap.add_argument("--stations", type=float, nargs="+", default=[0.0, None])
ap.add_argument("--zmax", type=float, default=300.0)
ap.add_argument("--klim", type=float, nargs=2, default=None, help="k axis range (m2/s2); default 0 to twice the equilibrium k")
a = ap.parse_args()
sys.path.insert(0, a.readerdir)
import erf_plotfile as ep  # noqa: E402

def column(d, hdr, i):
    lo, hi = hdr["lo"], hdr["hi"]; ny, nz = hi[1]-lo[1]+1, hi[2]-lo[2]+1
    m = lambda f, k: sum(d[f][i][j][k] for j in range(ny)) / ny  # noqa: E731
    z = [m("z_phys", k) for k in range(nz)]
    u = [math.hypot(m("x_velocity", k), m("y_velocity", k)) for k in range(nz)]
    return z, u, [m("KE", k) for k in range(nz)]

fig, ax = plt.subplots(1, 2, figsize=(9, 4.2), dpi=150)
zz = [a.z0 + (a.zmax - a.z0) * t / 200 for t in range(201)]
ax[0].plot([a.ustar / a.kappa * math.log((z + a.z0) / a.z0) for z in zz], zz, "k--", lw=1.2, label="log law")
ax[1].plot([a.ustar**2 / a.cmu0**2] * len(zz), zz, "k--", lw=1.2, label="equilibrium $u_*^2/C_{\\mu0}^2$")
colors = ["#d95f02", "#1b6ca8", "tab:red", "tab:green", "tab:purple"]  # first two: without / with the sensor (plot_paper.py)
for c, run in zip(colors, a.runs):
    label, pf = run.split("=", 1)
    hdr, d = ep.read_fields(pf, ["x_velocity", "y_velocity", "KE", "z_phys"])
    lo, hi = hdr["lo"], hdr["hi"]; nx = hi[0]-lo[0]+1
    xc = [hdr["prob_lo"][0] + (i + 0.5) * hdr["dx"][0] for i in range(nx)]
    for st, ls in zip(a.stations, ["-", ":"]):
        i = nx - 1 if st is None else min(range(nx), key=lambda i: abs(xc[i] - st))
        z, u, k = column(d, hdr, i)
        name = ("%s, hill position" % label if abs(st) < 1e-9 else "%s, x = %.0f m" % (label, xc[i])) if st is not None else "%s, outflow" % label
        ax[0].plot(u, z, ls, color=c, lw=1.4, label=name)
        ax[1].plot(k, z, ls, color=c, lw=1.4, label=name)
for x, t in zip(ax, ["$U$ (m/s)", "$k$ (m$^2$/s$^2$)"]):
    x.set_xlabel(t); x.set_ylabel("z (m)"); x.set_yscale("log"); x.set_ylim(0.5, a.zmax); x.grid(alpha=0.3)
keq = a.ustar**2 / a.cmu0**2
ax[1].set_xlim(*(a.klim or (0.0, 2.0 * keq)))   # a wide range: the drift reads as a fraction of k
ax[0].set_xlim(0.0, None)
ax[0].legend(fontsize=7); ax[1].legend(fontsize=7, loc="upper right")
fig.tight_layout(); fig.savefig(a.out)
print("wrote", a.out)
