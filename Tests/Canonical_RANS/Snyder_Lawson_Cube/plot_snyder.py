#!/usr/bin/env python3
"""Profiles of the Snyder and Lawson (1994) cube in the layout of Wise et al. (2026), Fig. A3.

    plot_snyder.py OUT.png PLOTFILE_OR_PREFIX:LABEL [PLOTFILE_OR_PREFIX:LABEL ...]

(a) the horizontal plane z = 0.02 m and (b) the vertical plane y = 0. Each profile is drawn from its station's
grey line as x = x_s + 0.025 u (the arrows give the scale), the measurements as black dots, each run as a line.
Unlike the paper's panel (a), y is the true lateral position (its labels there read twice y). A prefix takes the
last plotfile written with it. Not a test.
"""
import glob
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from compare_snyder import Field, read_meas, X_FACE, Y_MID, H  # noqa: E402

SCALE = 0.025
ST = {"horizontal": [-0.3, -0.2, -0.1, 0.04, 0.16, 0.3, 0.4, 0.5, 0.7],
      "vertical": [-0.3, -0.2, -0.1, 0.06, 0.14, 0.3, 0.4, 0.5, 0.7]}
COLOURS = ["tab:blue", "tab:red", "tab:green", "tab:orange", "tab:purple"]
STYLES = ["-", "-", "--", ":", "-"]


def main():
    out = sys.argv[1]
    runs = []
    for arg in sys.argv[2:]:
        path, label = arg.rsplit(":", 1)
        if not os.path.isdir(path) or not os.path.exists(os.path.join(path, "Header")):
            path = sorted(glob.glob(path + "[0-9]*"))[-1]
        runs.append((Field(path), label))
    meas = read_meas(os.path.join(HERE, "snyder_lawson_u.csv"))
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.0, 7.6), gridspec_kw={"height_ratios": [0.8, 0.4]})
    for ax, plane, top in ((ax1, "horizontal", 0.4), (ax2, "vertical", 0.4)):
        for xs in ST[plane]:
            ax.axvline(xs, color="0.8", lw=2.5, zorder=0)
        if plane == "horizontal":
            ax.add_patch(plt.Rectangle((0, -H / 2), H, H, color="k", zorder=3))
        else:
            ax.add_patch(plt.Rectangle((0, 0), H, H, color="k", zorder=3))
        for n, (f, label) in enumerate(runs):
            for xs in ST[plane]:
                if plane == "horizontal":
                    # Beside the cube, its two sides apart.
                    parts = [(-0.4, -H / 2), (H / 2, 0.4)] if 0 <= xs <= H else [(-0.4, 0.4)]
                    for k, (c0, c1) in enumerate(parts):
                        cs = [c0 + (c1 - c0) * m / 100 for m in range(101)]
                        ax.plot([xs + SCALE * f.at(X_FACE + xs, Y_MID + c, 0.02) for c in cs], cs, color=COLOURS[n],
                                ls=STYLES[n], lw=1.6, label=label if (xs == ST[plane][0] and k == 0) else None, zorder=2)
                else:
                    z0 = H if 0 <= xs <= H else 0.0
                    cs = [z0 + 0.002 * k for k in range(int((top - z0) / 0.002) + 1)]
                    ax.plot([xs + SCALE * f.at(X_FACE + xs, Y_MID, c) for c in cs], cs, color=COLOURS[n], ls=STYLES[n], lw=1.6,
                            label=label if xs == ST[plane][0] else None, zorder=2)
        m = [r for r in meas if r[0] == plane]
        ax.plot([r[1] + SCALE * r[4] for r in m], [r[2] if plane == "horizontal" else r[3] for r in m], "o",
                color="k", ms=4.5, zorder=4, label="Snyder and Lawson (1994)" if plane == "horizontal" else None)
        ax.set_xlim(-0.32, 0.8)
        ax.set_aspect("equal")
        ax.set_xticks([-0.3 + 0.1 * k for k in range(12)])
        if plane == "horizontal":
            ax.set_ylim(-0.4, 0.4)
            ax.set_ylabel("y [m]")
            ax.set_title("(a) z = 0.02 m", loc="left", fontsize=10)
            ax.annotate("", xy=(-0.3 + 2 * SCALE, -0.35), xytext=(-0.3, -0.35), arrowprops=dict(arrowstyle="->", lw=1.5))
            ax.text(-0.3 + 2 * SCALE + 0.01, -0.35, r"2 m s$^{-1}$", va="center", fontsize=8)
            fig.legend(loc="lower center", fontsize=8, ncol=2, frameon=False)
        else:
            ax.set_ylim(0, 0.4)
            ax.set_ylabel("z [m]")
            ax.set_xlabel("x [m] from the cube's upwind face")
            ax.set_title("(b) y = 0", loc="left", fontsize=10)
            ax.annotate("", xy=(0.6 + 3 * SCALE, 0.37), xytext=(0.6, 0.37), arrowprops=dict(arrowstyle="->", lw=1.5))
            ax.text(0.6 + 3 * SCALE + 0.01, 0.37, r"3 m s$^{-1}$", va="center", fontsize=8)
    fig.tight_layout(rect=(0, 0.04 + 0.02 * ((len(runs) + 2) // 2), 1, 1))
    fig.savefig(out, dpi=200)
    print("wrote", out)


if __name__ == "__main__":
    main()
