#!/usr/bin/env python3
"""The kinds of building the geometry holds, on flat ground and on a hill.

    plot_building_kinds.py GUST OUTDIR

Writes an "x y height zstart" map, runs GUST (no stage) on a flat and on a
hilly mesh in the current folder, and plots building_beta on a slice along x
from each plotfile, wall cells outlined. The two PNGs go in OUTDIR
(docs/figures/buildings).
"""
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402

# x from, x to, top, start
PARTS = [(60, 100, 30, 0), (140, 160, 20, 0), (160, 190, 35, 0), (230, 280, 30, 10), (320, 345, 25, 0), (350, 375, 25, 0),
         (420, 422.5, 20, 0), (462.5, 500, 17.5, 0), (540, 620, 22, 14)]
LABELS = [(80, "block"), (165, "stepped roof"), (255, "void deck"), (347.5, "narrow street"), (421, "thin wall"),
          (481, "off-mesh wall,\nroof in a layer"), (580, "lifted slab")]
DECK = """mesh.n_cell = 140 8 32
mesh.prob_lo = 0. 0. 0.
mesh.prob_hi = 700. 40. 160.
mesh.is_periodic = 0 1 0
mesh.max_grid_size = 64
buildings.type = immersed_forcing
buildings.format = xyz
buildings.file = kinds.xyz
io.write_checkpoint = 0
"""
HILL = """terrain.type = agnesi
terrain.hill_shape = ridge_x
terrain.hill_height = 30.
terrain.hill_half_width = 150.
terrain.smoothing = 0
"""


def height(x):
    for a, b, h, z in PARTS:
        if a < x < b:
            return h, z
    return 0.0, 0.0


def plot(plt_dir, title, out):
    hdr, d = erf_plotfile.read_fields(plt_dir, ["building_beta", "building_wall", "z", "detJ"])
    nx, _, nz = [hdr["hi"][n] - hdr["lo"][n] + 1 for n in range(3)]
    dx, _, dzeta = hdr["dx"]
    j = 4
    polys, vals, walls = [], [], []
    for i in range(nx):
        for k in range(nz):
            zc, dz = d["z"][i][j][k], d["detJ"][i][j][k] * dzeta
            q = [(i * dx, zc - dz / 2), ((i + 1) * dx, zc - dz / 2), ((i + 1) * dx, zc + dz / 2), (i * dx, zc + dz / 2)]
            polys.append(q)
            vals.append(d["building_beta"][i][j][k])
            if d["building_wall"][i][j][k] > 0:
                walls.append(q)
    ground = [d["z"][i][j][0] - d["detJ"][i][j][0] * dzeta / 2 for i in range(nx)]
    fig, ax = plt.subplots(figsize=(11, 3.6), dpi=150)
    pc = PolyCollection(polys, array=vals, cmap="Greys", clim=(0, 1), edgecolors="#e4e4e4", linewidths=0.2)
    ax.add_collection(pc)
    ax.add_collection(PolyCollection(walls, facecolors="none", edgecolors="#c8501e", linewidths=0.8, label="wall cell"))
    xs = [(i + 0.5) * dx for i in range(nx)]
    ax.fill_between(xs, min(ground) - 5, ground, color="#b9a78c", linewidth=0)
    for x, text in LABELS:
        i = int(x / dx)
        ax.text(x, ground[i] + 40, text, ha="center", va="bottom", fontsize=8, color="#333333")
    ax.set_xlim(0, 700)
    ax.set_ylim(min(ground) - 5, max(ground) + 62)
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    ax.set_title(title, loc="left", fontsize=10)
    ax.legend(loc="upper right", frameon=False, fontsize=8)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.colorbar(pc, ax=ax, label="solid fraction", fraction=0.025, pad=0.01)
    fig.tight_layout()
    fig.savefig(out)


def main():
    gust, outdir = sys.argv[1:3]
    with open("kinds.xyz", "w") as f:
        for ti in range(280):
            for tj in range(16):
                x, y = 1.25 + 2.5 * ti, 1.25 + 2.5 * tj
                h, z = height(x)
                f.write(f"{x} {y} {h} {z}\n")
    for name, extra, title in (("flat", "", "Flat ground"), ("hill", HILL, "On a hill: heights from the local ground")):
        with open(f"inputs_{name}", "w") as f:
            f.write(DECK + extra + f"io.plot_file = {name}_plt\n")
        subprocess.run([gust, f"inputs_{name}"], check=True)
        plot(f"{name}_plt00000", title, os.path.join(outdir, f"kinds_{name}.png"))


if __name__ == "__main__":
    main()
