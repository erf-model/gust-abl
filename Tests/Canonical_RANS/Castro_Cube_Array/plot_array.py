#!/usr/bin/env python3
"""Plot GUST's profiles in the cube array against the measurements, in the layout of Munoz-Esparza et al. (2020, Fig. 3).

    plot_array.py MEAS.csv OUT.png LABEL:PLOTFILE:TURBPLOT [LABEL:PLOTFILE:TURBPLOT ...]

(a) u / u_r at P0 and P3, (b) u / u_r at P1 and P2, (c) -<u'w'> / u_r^2 at P1 and P3. Lines are GUST runs (one line
style per run), symbols the measurements. Not a test.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_array import Field, H, STATIONS, U_R, read_meas  # noqa: E402

PANELS = [("a", "u", ("P0", "P3")), ("b", "u", ("P1", "P2")), ("c", "uw", ("P1", "P3"))]
COLOUR = {"P0": "#1b1b1b", "P1": "#1b1b1b", "P2": "#1f4fd1", "P3": "#8a8a8a"}
MARKER = {("a", "P0"): "s", ("a", "P3"): "^", ("b", "P1"): "o", ("b", "P2"): "x", ("c", "P1"): "o", ("c", "P3"): "^"}
STYLES = ["-", "--", ":", "-."]


def main():
    meas, out, runs = sys.argv[1], sys.argv[2], sys.argv[3:]
    m = read_meas(meas)
    fields = []
    for r in runs:
        label, pf, tp = r.split(":")
        fields.append((label, Field(pf, tp)))
    zs = [0.01 * k * H for k in range(1, 301)]
    fig, axs = plt.subplots(1, 3, figsize=(11, 4.2), sharey=True)
    for ax, (tag, qn, sts) in zip(axs, PANELS):
        for st in sts:
            for (label, f), ls in zip(fields, STYLES):
                x, y = STATIONS[st]
                v = [f.at("u" if qn == "u" else "uw", x, y, z) for z in zs]
                v = [w / U_R for w in v] if qn == "u" else [w / U_R ** 2 * 1e3 for w in v]
                ax.plot(v, [z / H for z in zs], ls, color=COLOUR[st], lw=1.4,
                        label=f"{st} GUST {label}" if len(fields) > 1 else f"{st} GUST")
            pts = [(z, v) for s, q, z, v in m if s == st and q == qn]
            scale = 1e3 if qn == "uw" else 1.0
            mk = MARKER[(tag, st)]
            ax.plot([v * scale for z, v in pts], [z for z, v in pts], mk, color=COLOUR[st],
                    mfc="none" if mk != "x" else COLOUR[st], ms=5, label=f"{st} exp.")
        ax.set_ylim(0, 3)
        ax.grid(alpha=0.3)
        ax.set_xlabel("u / u_r" if qn == "u" else "-<u'w'> / u_r^2  (x 1e-3)")
        ax.set_title(f"({tag})", loc="left")
        ax.legend(fontsize=7, loc="upper left")
    axs[0].set_ylabel("z / h")
    axs[0].set_xlim(-0.2, 0.8)
    axs[1].set_xlim(-0.2, 0.8)
    axs[2].set_xlim(-2, 12)
    fig.tight_layout()
    fig.savefig(out, dpi=150)


if __name__ == "__main__":
    main()
