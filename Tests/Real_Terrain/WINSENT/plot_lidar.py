#!/usr/bin/env python3
"""Plot one GUST run against the valley lidar L140, plus its MMNW profile (no measurement there).

Reads a sample_stations.py CSV and the lidar profile (prepare_gust.py's inflow.txt, the case mean in height above
the ground) and writes a three-panel PNG: L140 speed (with the 100-400 m bias and RMS error), L140 direction, and
the MMNW speed.
"""
import argparse
import csv
import glob
import math
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

LIDAR_GROUND_NHN = 468.8


def lidar_spread(folder):
    """Heights above the lidar ground, the periods' mean speed and its standard deviation, the vector-mean
    direction and its circular standard deviation, and the number of periods (ERF's lidar_spread_plot.py)."""
    files = sorted(glob.glob(os.path.join(folder, "*_lidar_data.csv")))
    acc = {}
    for fn in files:
        with open(fn) as fh:
            for r in csv.DictReader(fh):
                h = round(float(r["height"]) - LIDAR_GROUND_NHN, 2)
                s, d = float(r["wind_speed"]), math.radians(float(r["wind_direction"]))
                acc.setdefault(h, []).append((s, -math.sin(d) * s, -math.cos(d) * s))
    hs = sorted(h for h, v in acc.items() if len(v) >= 0.5 * len(files))
    s_mean = np.array([np.mean([v[0] for v in acc[h]]) for h in hs])
    s_sd = np.array([np.std([v[0] for v in acc[h]]) for h in hs])
    d_mean, d_sd = [], []
    for h in hs:
        fe, fn_ = np.mean([v[1] for v in acc[h]]), np.mean([v[2] for v in acc[h]])
        dm = (math.degrees(math.atan2(fe, fn_)) + 180.0) % 360.0
        d_mean.append(dm)
        d_sd.append(np.std([((math.degrees(math.atan2(v[1], v[2])) + 180.0 - dm + 180.0) % 360.0) - 180.0
                            for v in acc[h]]))
    return np.array(hs), s_mean, s_sd, np.array(d_mean), np.array(d_sd), len(files)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv", help="sample_stations.py output")
    p.add_argument("--lidar", required=True, help="the lidar profile: z u v theta k")
    p.add_argument("--wd", type=float, default=293.0, help="the frame's wind direction (deg from north)")
    p.add_argument("--label", default="GUST")
    p.add_argument("--title", default="WINSENT high LAI, neutral: GUST against the lidar")
    p.add_argument("--lidar-folder", help="the case's lidar period files (*_lidar_data.csv): shade +-1 sd")
    p.add_argument("--powerlaw", type=float, nargs=2, metavar=("C", "a"),
                   help="the benchmark power law U = C h^a (high LAI neutral: 1.033 0.323)")
    p.add_argument("--canopy-height", type=float, default=25.0, help="canopy height shaded on the MMNW panel (m)")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    g = pd.read_csv(a.csv, comment="#")
    r = np.loadtxt(a.lidar)
    h = np.arange(100, 401, 10.0)
    lu, lv = np.interp(h, r[:, 0], r[:, 1]), np.interp(h, r[:, 0], r[:, 2])
    lw = np.hypot(lu, lv)
    ld = (a.wd - np.degrees(np.arctan2(lv, lu))) % 360
    d, m = g[g.station == "L140"], g[g.station == "MMNW"]
    b = d.speed.values - np.interp(d.agl, h, lw)

    lidar_c, run_c = "#2a2a28", "#2a78d6"
    mk = dict(color=run_c, lw=2, marker="o", ms=8, mec="white")
    fig, ax = plt.subplots(1, 3, figsize=(13, 5.6))
    if a.lidar_folder:
        hs, sm, ssd, dm, dsd, n = lidar_spread(a.lidar_folder)
        ax[0].fill_betweenx(hs, sm - ssd, sm + ssd, color="#d9d8d2", lw=0, label=f"lidar ±1 sd ({n} periods)")
        ax[1].fill_betweenx(hs, dm - dsd, dm + dsd, color="#d9d8d2", lw=0)
    if a.powerlaw:
        hp = np.linspace(100, 400, 60)
        ax[0].plot(a.powerlaw[0] * hp ** a.powerlaw[1], hp, color="#55554f", lw=1.6, ls="--",
                   label=f"benchmark power law {a.powerlaw[0]:g} h^{a.powerlaw[1]:g}")
    ax[0].plot(lw, h, color=lidar_c, lw=2, label="lidar L140 (case mean)")
    ax[0].plot(d.speed, d.agl, label=a.label, **mk)
    ax[0].text(0.03, 0.97, f"100-400 m: bias {b.mean():+.2f} m/s\nrms {np.sqrt((b ** 2).mean()):.2f} m/s",
               transform=ax[0].transAxes, va="top", fontsize=9, color="#55554f")
    ax[0].set_title("L140 (valley lidar): wind speed")
    ax[0].set_xlabel("speed (m/s)")
    ax[0].set_ylabel("height above ground (m)")
    ax[1].plot(ld, h, color=lidar_c, lw=2)
    ax[1].plot(d.dir_from, d.agl, **mk)
    ax[1].set_title("L140: wind direction")
    ax[1].set_xlabel("direction from (deg)")
    ax[1].set_xlim(285, 305)
    ax[2].plot(m.speed, m.agl, **mk)
    ax[2].axhspan(0, a.canopy_height, color="#1baf7a", alpha=0.12, lw=0)
    ax[2].text(0.05, 0.03, f"canopy ({a.canopy_height:g} m)", transform=ax[2].transAxes, color="#55554f", fontsize=9)
    ax[2].set_title("MMNW (plateau mast): GUST only, no data")
    ax[2].set_xlabel("speed (m/s)")
    for x in ax:
        x.grid(color="#e4e3dd", lw=0.8)
        x.set_axisbelow(True)
        for sp in ("top", "right"):
            x.spines[sp].set_visible(False)
    fig.suptitle(a.title, fontsize=13, y=0.995)
    fig.legend(*ax[0].get_legend_handles_labels(), loc="upper center", bbox_to_anchor=(0.5, 0.945), ncol=4,
               frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    fig.savefig(a.out, dpi=140, bbox_inches="tight")
    print(f"{a.out}: L140 100-400 m bias {b.mean():+.2f} m/s, rms {np.sqrt((b ** 2).mean()):.2f} m/s")


if __name__ == "__main__":
    main()
