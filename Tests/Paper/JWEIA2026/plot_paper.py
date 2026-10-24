#!/usr/bin/env python3
"""The paper's comparison figures: observations, GUST - CFD without the lee sensor and GUST - CFD with it.

    python3 plot_paper.py --runs RUNS [all|sites|ridge|hill] [--outdir DIR] [--manning CSV]

Field sites (Askervein, Bolund cases 1, 3, 4, Big Southern Butte): the sampled csv files of the two runs
(windninja2019.py sample) and the masts. Ishihara-Zhou ridge and RUSHIL hills: the two plotfiles, with the
tunnel data (and the tunnel/LES reattachment ranges for the ridge). Only these three sources are drawn.
RUNS is the --out directory of paper_runs.py (runs M_nosensor_* and M_sensor_*).
"""
import argparse
import csv
import glob
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
for d in ("tools", "Tests/Real_Terrain/IshiharaZhou", "Tests/Real_Terrain/RUSHIL"):
    sys.path.insert(0, os.path.join(REPO, d))
R = "."          # set by --runs
MANNING = None   # set by --manning
import askervein_lines as al  # noqa: E402
import ishihara_zhou as iz  # noqa: E402
import rushil  # noqa: E402

OBS = "black"
RUNS = [  # (key, label, style): the order is the legend order
    ("nosensor", "GUST - CFD without sensor", dict(color="#d95f02", ls="--", lw=1.5)),
    ("sensor", "GUST - CFD with sensor", dict(color="#1b6ca8", ls="-", lw=1.6)),
]
SITES = {  # name: (site, case, title); samples R/samples/M_sensor_<name>.csv and M_nosensor_<name>.csv
    "askervein": ("askervein", "210", "Askervein"),
    "bolund_270": ("bolund", "270", "Bolund case 1 (270°)"),
    "bolund_239": ("bolund", "239", "Bolund case 3 (239°)"),
    "bolund_90": ("bolund", "90", "Bolund case 4 (90°)"),
    "bsb": ("bsb", "222", "Big Southern Butte"),
}


def converged(tag):
    log = f"{R}/runs/{tag}/log.txt"
    return os.path.exists(log) and "steady after" in open(log).read()


def last_plt(d):
    p = sorted(x for x in glob.glob(f"{d}/plt[0-9]*") if "level" not in x)
    if not p:
        raise SystemExit(f"no plotfile in {d}")
    return p[-1]


def have(path):
    if not os.path.exists(path):
        print(f"  skipped: {path} is missing")
        return False
    return True


def sites(outdir):
    masts = list(csv.DictReader(open(f"{R}/wn/masts.csv")))
    for name, (site, case, title) in SITES.items():
        s_csv, n_csv = f"{R}/samples/M_sensor_{name}.csv", f"{R}/samples/M_nosensor_{name}.csv"
        if not have(s_csv):
            continue
        # The Butte without the sensor does not converge: observations and the sensor run only.
        pairs = [(r, path) for r, path in zip(RUNS, (n_csv, s_csv))
                 if os.path.exists(path) and (name != "bsb" or converged(f"M_{r[0]}_{name}"))]
        obs = [r for r in masts if r["site"] == site and r["case"] == case]
        lines = [n for n in dict.fromkeys(r["line"] for r in obs) if n != "-"]
        fig, axes = plt.subplots(1, len(lines), figsize=(4.2 * len(lines), 3.6), squeeze=False)
        for ax, line in zip(axes[0], lines):
            for (key, label, style), path in pairs:
                rows = sorted((r for r in csv.DictReader(open(path)) if r["line"] == line and not r["mast"]),
                              key=lambda r: float(r["d"]))
                ax.plot([float(r["d"]) for r in rows], [float(r["speedup"]) for r in rows], label=label, **style)
            o = [r for r in obs if r["line"] == line]
            ax.plot([float(r["d"]) for r in o], [float(r["speedup"]) for r in o], "o", color=OBS, ms=4,
                    label="observed")
            ax.axhline(0, color="0.8", lw=0.8)
            ax.set_title(f"{title}: line {line}", fontsize=10)
            ax.set_xlabel("distance (m)")
            ax.set_ylabel("speed-up")
            ax.grid(alpha=0.3)
        handles, labels = axes[0][0].get_legend_handles_labels()
        one = len(lines) == 1   # a single panel is too narrow for the legend in one row
        fig.legend(handles, labels, loc="lower center", ncol=1 if one else 3, fontsize=8, frameon=False)
        fig.tight_layout(rect=(0, 0.2 if one else 0.07, 1, 1))
        out = f"{outdir}/{name}.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(out)


def ridge_line(cols, s):
    xs, us = [], []
    for i in range(0, 551):
        xr = -150.0 + i
        u, zs = iz.column_u(cols, iz.XC + xr)
        zrel = [z - iz.ground(iz.XC + xr) for z in zs]
        zt = s["zline"]
        us.append(u[0] * math.log((zt + s["z0"]) / s["z0"]) / math.log((zrel[0] + s["z0"]) / s["z0"])
                  if zt < zrel[0] else iz.at_height(u, zrel, zt))
        xs.append(xr)
    return xs, us


def ridge(outdir):
    for surf in ("smooth", "rough"):
        dirs = {k: f"{R}/runs/M_{k}_ridge_{surf}" for k in ("nosensor", "sensor")}
        if not all(have(d) and glob.glob(f"{d}/plt[0-9]*") for d in dirs.values()):
            continue
        cols = {k: al.Columns(last_plt(d)) for k, d in dirs.items()}
        s = iz.SURF[surf]
        st = list(csv.DictReader(open(f"{R}/iz/stations_{surf}.csv")))
        fig = plt.figure(figsize=(13, 5.5))
        ax = fig.add_subplot(1, 3, 1)
        for key, label, style in RUNS:
            xs, us = ridge_line(cols[key], s)
            ax.plot(xs, us, label=label, **style)
        line = [(float(r["x"]), float(r["U"])) for r in st if abs(float(r["z"]) - s["zline"]) < 1e-6]
        ax.plot([p[0] for p in line], [p[1] for p in line], "o", color=OBS, ms=4,
                label=f"tunnel (z = {s['zline']:g} m)")
        ref = iz.REF[surf]
        ax.axvspan(*ref["tunnel_reattach"], color=OBS, alpha=0.12, label="tunnel reattachment")
        ax.axvline(ref["les"][1], color="grey", ls=":", label="LES reattachment")
        ax.axhline(0, color="0.7", lw=0.8)
        ax.set_xlabel("x from the crest (m)")
        ax.set_ylabel(f"u at {s['zline']:g} m above the ground (m/s)")
        ax.set_title(f"Ridge, {surf}: near-surface line", fontsize=10)
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ax2 = fig.add_subplot(1, 3, (2, 3))
        xst = sorted({float(r["x"]) for r in st})
        sc = 50 / 10.0   # 50 m per 10 m/s
        for x in xst:
            g = iz.ground(iz.XC + x)
            for key, label, style in RUNS:
                u, zs = iz.column_u(cols[key], iz.XC + x)
                zz = [z for z in zs if z - g <= 200]
                ax2.plot([x + u[k] * sc for k in range(len(zz))], zz, label=label if x == xst[0] else None,
                         **dict(style, lw=1.1))
            pts = sorted((float(r["z"]), float(r["U"])) for r in st if float(r["x"]) == x)
            ax2.plot([x + p[1] * sc for p in pts], [p[0] for p in pts], "o", ms=3, color=OBS,
                     label="tunnel" if x == xst[0] else None)
            ax2.axvline(x, color="0.85", lw=0.6)
        gx = [iz.XC - 200 + 2 * i for i in range(301)]
        ax2.fill_between([x - iz.XC for x in gx], [iz.ground(x) for x in gx], 0, color="0.75")
        x0, y0 = -215, 182
        ax2.annotate("", xy=(x0 + 10 * sc, y0), xytext=(x0, y0), arrowprops=dict(arrowstyle="|-|", color="k", lw=1.2, mutation_scale=4))
        ax2.text(x0 + 5 * sc, y0 + 5, "10 m/s", ha="center", va="bottom", fontsize=9)
        ax2.set_xlabel("x from the crest (m); profiles of U drawn from each station (scale bar: 10 m/s)")
        ax2.set_ylabel("z (m)")
        ax2.set_ylim(0, 200)
        ax2.set_title("U profiles at the tunnel stations", fontsize=10)
        ax2.legend(fontsize=8, loc="upper right")
        fig.tight_layout()
        out = f"{outdir}/ridge_{surf}.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(out)


def manning_measured(h, zf):
    """The digitized tunnel points of Manning et al. (2010), Figs. 4.3/4.4 (Hills 3/5); [] for other hills.
    Their Meteodyn WT curves are in the same file and are not drawn."""
    fig = {3: "fig4_3", 5: "fig4_4"}.get(h)
    path = MANNING
    if not fig or not path or not os.path.exists(path):
        return []
    return [(float(r["x"]), float(r["y"])) for r in csv.DictReader(l for l in open(path) if l[0] != "#")
            if r["figure"] == fig and r["panel"] == f"Z{zf:g}h" and r["series"] == "measured"]


def hill(outdir, hills=(3, 5, 8)):
    for h in hills:
        dirs = {k: f"{R}/runs/M_{k}_hill_{h}" for k in ("nosensor", "sensor")}
        if not all(have(d) and glob.glob(f"{d}/plt[0-9]*") for d in dirs.values()):
            continue
        hl = rushil.HILLS[h] * rushil.H
        g = rushil.ground_fn(hl)
        cols = {k: al.Columns(last_plt(d)) for k, d in dirs.items()}
        uref = {}
        for k, c in cols.items():
            us, zs = rushil.column(c, 200.0)
            uref[k] = rushil.u_at(us, zs, 0.8 * rushil.H)
        st = list(csv.DictReader(open(f"{R}/rushil/stations_h{h}.csv")))
        uref_t = uref["sensor"]   # the tunnel is scaled as in benchmark_plots.py (by the run's inflow)
        fig, axs = plt.subplots(3, 4, figsize=(15, 10))
        for n, zf in enumerate((0.1, 0.2, 0.3, 0.4, 0.5, 0.6)):
            ax = axs.flat[n]
            for key, label, style in RUNS:
                xa, sp = [], []
                for i in range(-500, 601):
                    x = rushil.XC + i / 100.0 * hl
                    if not 50 < x < 9300:
                        continue
                    u, z = rushil.column(cols[key], x)
                    xa.append(i / 100.0)
                    sp.append(rushil.u_at(u, [zz - g(x) for zz in z], zf * rushil.H) / uref[key])
                ax.plot(xa, sp, label=label, **style)
            meas = manning_measured(h, zf)
            if meas:   # Hills 3 and 5: the tunnel points as digitized by Manning et al. (2010), as before
                ax.plot([q[0] for q in meas], [q[1] for q in meas], "o", ms=4, color=OBS, label="tunnel")
                ax.axhline(0, color="0.7", lw=0.8)
                ax.set_xlim(-5, 6)
                ax.set_xlabel("X / a")
                ax.set_title(f"Z = {zf:g} H", fontsize=10)
                ax.grid(alpha=0.3)
                if n % 4 == 0:
                    ax.set_ylabel("U / U_in(0.8 H)")
                continue
            pts = []
            for x in sorted({float(r["x"]) for r in st}):
                prof = sorted((float(r["z"]), float(r["U"])) for r in st if float(r["x"]) == x)
                zt = zf * rushil.H
                for k in range(len(prof) - 1):
                    if prof[k][0] <= zt <= prof[k + 1][0]:
                        t = (zt - prof[k][0]) / (prof[k + 1][0] - prof[k][0])
                        pts.append((x / hl, ((1 - t) * prof[k][1] + t * prof[k + 1][1]) / uref_t))
                        break
            ax.plot([p[0] for p in pts], [p[1] for p in pts], "o", ms=4, color=OBS, label="tunnel")
            ax.axhline(0, color="0.7", lw=0.8)
            ax.set_xlim(-5, 6)
            ax.set_xlabel("X / a")
            ax.set_title(f"Z = {zf:g} H", fontsize=10)
            ax.grid(alpha=0.3)
            if n % 4 == 0:
                ax.set_ylabel("U / U_in(0.8 H)")
        stations = [x for x in sorted({float(r["x"]) for r in st}) if -1.0 * hl <= x <= 3.2 * hl][:6]
        for n, x in enumerate(stations):
            ax = axs.flat[6 + n]
            for key, label, style in RUNS:
                u, z = rushil.column(cols[key], rushil.XC + x)
                zr = [zz - g(rushil.XC + x) for zz in z]
                k = [i for i, v in enumerate(zr) if v <= 1.2 * rushil.H]
                ax.plot([u[i] / uref[key] for i in k], [zr[i] / rushil.H for i in k], label=label, **style)
            prof = sorted((float(r["z"]), float(r["U"])) for r in st if float(r["x"]) == x)
            ax.plot([p[1] / uref_t for p in prof], [p[0] / rushil.H for p in prof], "o", ms=3, color=OBS,
                    label="tunnel")
            ax.axvline(0, color="0.7", lw=0.8)
            ax.set_ylim(0, 1.2)
            ax.set_title(f"X = {x / hl:.2f} a", fontsize=10)
            ax.grid(alpha=0.3)
            if (6 + n) % 4 == 0:
                ax.set_ylabel("Z / H")
            ax.set_xlabel("U / U_in(0.8 H)")
        for n in range(6 + len(stations), 12):
            axs.flat[n].axis("off")
        fig.suptitle(f"RUSHIL Hill {h} (a / H = {h})", fontsize=12)
        hd, lb = axs.flat[0].get_legend_handles_labels()
        fig.legend(hd, lb, loc="lower center", ncol=3, fontsize=10, frameon=False)
        fig.tight_layout(rect=(0, 0.04, 1, 0.97))
        out = f"{outdir}/hill_{h}.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        print(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", nargs="?", default="all", choices=["all", "sites", "ridge", "hill"])
    ap.add_argument("--runs", required=True, help="the --out directory of paper_runs.py")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--manning", help="tunnel points of Manning et al. (2010) digitized by the user (not distributed)")
    a = ap.parse_args()
    global R, MANNING
    R, MANNING = os.path.abspath(a.runs), a.manning
    if a.what in ("all", "sites"):
        sites(a.outdir)
    if a.what in ("all", "ridge"):
        ridge(a.outdir)
    if a.what in ("all", "hill"):
        hill(a.outdir)


if __name__ == "__main__":
    main()
