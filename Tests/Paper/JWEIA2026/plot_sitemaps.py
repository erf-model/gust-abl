#!/usr/bin/env python3
"""Site maps of the field cases: terrain height (filled contours), measurement lines and masts, and the wind
direction of each case, in the coordinates of the computational domain.

    python3 plot_sitemaps.py --runs RUNS [--out sitemaps.png]

Terrain: RUNS/wn/bolund.xyz and bsb.xyz (x y z points on a regular grid, from paper_runs.py prepare) and the ERF
terrain file of Askervein (extern_erf/askervein.txt: nx, ny, the x and y node lists, then z). Masts: RUNS/wn/masts.csv. Wind directions: bc.inflow_direction of the decks
(the direction the wind blows towards, counter-clockwise from the x axis).
"""
import argparse
import csv
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
R = "."   # set by --runs
LINE_COLORS = {"A": "#1b6ca8", "AA": "#d95f02", "B": "#1a9850", "TSW": "#1b6ca8", "TWSW": "#d95f02"}


def read_xyz(path, stride=1):
    x, y, z = [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            a, b, c = line.split()[:3]
            x.append(float(a)); y.append(float(b)); z.append(float(c))
    x, y, z = np.array(x), np.array(y), np.array(z)
    xs, ys = np.unique(x), np.unique(y)
    Z = np.full((len(ys), len(xs)), np.nan)
    Z[np.searchsorted(ys, y), np.searchsorted(xs, x)] = z
    return xs[::stride], ys[::stride], Z[::stride, ::stride]


def read_erf(path):
    v = [float(t) for t in open(path).read().split()]
    nx, ny = int(v[0]), int(v[1])
    xs = np.array(v[2:2 + nx]); ys = np.array(v[2 + nx:2 + nx + ny])
    z = np.array(v[2 + nx + ny:2 + nx + ny + nx * ny])
    return xs, ys, z.reshape(ny, nx)


def masts(site):
    return [r for r in csv.DictReader(open(f"{R}/wn/masts.csv")) if r["site"] == site]


def arrow(ax, deg, label, xy, length, color="k"):
    dx, dy = length * math.cos(math.radians(deg)), length * math.sin(math.radians(deg))
    ax.annotate("", xy=(xy[0] + dx, xy[1] + dy), xytext=xy,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=1.6))
    ax.text(xy[0] - 0.15 * dx, xy[1] - 0.15 * dy, label, fontsize=8, ha="center", va="center", color=color)


def panel(ax, xs, ys, Z, title, levels):
    cf = ax.contourf(xs, ys, Z, levels=levels, cmap="Greys", extend="both")
    ax.contour(xs, ys, Z, levels=levels[::2], colors="0.35", linewidths=0.4)
    ax.set_aspect("equal")
    ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(ys[0], ys[-1])
    ax.set_title(title, fontsize=10)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    return cf


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs", required=True, help="the --out directory of paper_runs.py")
    ap.add_argument("--out", default="sitemaps.png")
    a = ap.parse_args()
    global R
    R = os.path.abspath(a.runs)
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.4))
    figb, axb = plt.subplots(1, 1, figsize=(10, 5.6))

    # Askervein: lines A, AA, B; wind 210 deg (towards 60 deg in the domain)
    xs, ys, Z = read_erf(os.path.join(REPO, "extern_erf", "askervein.txt"))
    k = (np.abs(xs) <= 1500); m = (np.abs(ys) <= 1500)
    cf = panel(axs[0], xs[k], ys[m], Z[np.ix_(m, k)], "Askervein (210$^\\circ$)", np.linspace(0, 120, 13))
    for line in ("A", "AA", "B"):
        pts = sorted(((float(r["x"]), float(r["y"])) for r in masts("askervein") if r["line"] == line))
        axs[0].plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=3, lw=1, color=LINE_COLORS[line], label=f"line {line}")
    arrow(axs[0], 60, "wind", (-1300, -1300), 450)
    axs[0].legend(fontsize=8, loc="upper left")
    fig.colorbar(cf, ax=axs[0], shrink=0.75, label="terrain height (m)")

    # Bolund (own figure): lines A and B with the masts; winds of cases 1, 3, 4
    # Big Southern Butte: transects TSW, TWSW and the other stations; wind 222 deg (towards 48 deg)
    xs, ys, Z = read_xyz(f"{R}/wn/bsb.xyz", stride=2)
    k = (xs >= 2000) & (xs <= 15000); m = (ys >= 5000) & (ys <= 17000)
    cf = panel(axs[1], xs[k], ys[m], Z[np.ix_(m, k)], "Big Southern Butte (222$^\\circ$)", np.linspace(0, 800, 17))
    st = masts("bsb")
    oth = [r for r in st if r["line"] == "-"]
    axs[1].plot([float(r["x"]) for r in oth], [float(r["y"]) for r in oth], "o", ms=3, color="#7f7f7f", label="stations")
    for line in ("TSW", "TWSW"):
        pts = sorted(((float(r["x"]), float(r["y"])) for r in st if r["line"] == line), key=lambda p: p[0] + p[1])
        axs[1].plot([p[0] for p in pts], [p[1] for p in pts], "o-", ms=3, lw=1, color=LINE_COLORS[line], label=line)
    arrow(axs[1], 48, "wind", (3000, 7000), 1500)
    axs[1].legend(fontsize=8, loc="upper left")
    fig.colorbar(cf, ax=axs[1], shrink=0.75, label="terrain height (m)")

    xs, ys, Z = read_xyz(f"{R}/wn/bolund.xyz")
    k = (xs >= -230) & (xs <= 170); m = (ys >= -130) & (ys <= 100)
    cfb = panel(axb, xs[k], ys[m], Z[np.ix_(m, k)], "Bolund", np.linspace(0, 12, 13))
    seen = {}
    for r in masts("bolund"):
        seen.setdefault(r["mast"][:2], (float(r["x"]), float(r["y"]), r["line"]))
    for name, (x, y, line) in sorted(seen.items()):
        axb.plot(x, y, "o", ms=5, color=LINE_COLORS[line], mec="k", mew=0.5)
        off = {"M0": (6, -12), "M1": (-22, -12), "M2": (-6, 8), "M3": (4, 8), "M4": (6, 4), "M6": (-4, 9), "M7": (-24, 6), "M8": (6, 6)}.get(name, (5, 5))
        axb.annotate(name, (x, y), textcoords="offset points", xytext=off, fontsize=10, weight="bold", bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85))
    for line in ("A", "B"):
        pts = sorted({(float(r["x"]), float(r["y"])) for r in masts("bolund") if r["line"] == line})
        axb.plot([pts[0][0] - 20, pts[-1][0] + 20],
                    [pts[0][1] - 20 * (pts[-1][1] - pts[0][1]) / max(pts[-1][0] - pts[0][0], 1e-9),
                     pts[-1][1] + 20 * (pts[-1][1] - pts[0][1]) / max(pts[-1][0] - pts[0][0], 1e-9)],
                    "--", lw=1, color=LINE_COLORS[line], label=f"line {line}")
    for deg, lab, xy, txt in ((0, "case 1 (270$^\\circ$)", (-215, 60), (-170, 75)), (31, "case 3 (239$^\\circ$)", (-225, -120), (-150, -120)),
                              (180, "case 4 (90$^\\circ$)", (160, 60), (120, 75))):
        dx, dy = 45 * math.cos(math.radians(deg)), 45 * math.sin(math.radians(deg))
        axb.annotate("", xy=(xy[0] + dx, xy[1] + dy), xytext=xy, arrowprops=dict(arrowstyle="-|>", color="k", lw=1.8))
        axb.text(txt[0], txt[1], lab, fontsize=9, ha="center")
    axb.legend(fontsize=8, loc="lower right")
    figb.colorbar(cfb, ax=axb, shrink=0.75, label="terrain height (m)")

    figb.tight_layout()
    figb.savefig(a.out.replace("sitemaps", "bolund_map"), dpi=150)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(a.out, a.out.replace("sitemaps", "bolund_map"))


if __name__ == "__main__":
    main()
