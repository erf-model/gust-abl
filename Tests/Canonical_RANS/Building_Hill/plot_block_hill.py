#!/usr/bin/env python3
"""Figure of the flow around the block on the ridge, for the docs and the PR. Not a test.

    plot_block_hill.py PLOTFILE OUT.png

u on the block's centre plane at true heights, with velocity arrows; the block and the ground are grey.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402


def main():
    hdr, d = erf_plotfile.read_fields(sys.argv[1], ["u", "w", "z", "building_beta"])
    u, w, z, b = [np.array(d[n], dtype=float) for n in ("u", "w", "z", "building_beta")]
    dx = hdr["dx"][0]
    nx = u.shape[0]
    jc, kmax = 14, 18
    x = np.tile(((np.arange(nx) + 0.5) * dx)[:, None], (1, kmax))
    zz, uu, ww, bb = z[:, jc, :kmax], u[:, jc, :kmax], w[:, jc, :kmax], b[:, jc, :kmax]
    fig, ax = plt.subplots(figsize=(10, 3.9), constrained_layout=True)
    norm = TwoSlopeNorm(vmin=min(-1.0, uu.min()), vcenter=0.0, vmax=uu.max())
    # Cell corners from the centres, so the cells show as the mesh has them.
    xe = np.tile((np.arange(nx + 1) * dx)[:, None], (1, kmax + 1))
    ze = np.zeros((nx + 1, kmax + 1))
    zc = np.vstack([zz[0], 0.5 * (zz[1:] + zz[:-1]), zz[-1]])
    ze[:, 1:-1] = 0.5 * (zc[:, 1:] + zc[:, :-1])
    ze[:, 0] = 1.5 * zc[:, 0] - 0.5 * zc[:, 1]
    ze[:, -1] = 1.5 * zc[:, -1] - 0.5 * zc[:, -2]
    m = ax.pcolormesh(xe, ze, np.ma.masked_where(bb >= 1, uu), cmap="RdBu_r", norm=norm, shading="flat")
    ax.pcolormesh(xe, ze, np.ma.masked_where(bb < 1, bb), cmap="Greys", vmin=0, vmax=2, shading="flat")
    s = 2
    p, q = np.ma.masked_where(bb >= 1, uu), np.ma.masked_where(bb >= 1, ww)
    ax.quiver(x[::s, ::s], zz[::s, ::s], p[::s, ::s], q[::s, ::s], color="0.15", width=0.0016, scale=170)
    ax.fill_between(xe[:, 0], 0, ze[:, 0], color="0.45")
    ax.set_aspect("equal")
    ax.set_ylim(0, ze[:, -1].min())
    ax.set_xlim(0, xe.max())
    ax.set_xlabel("x (m)")
    ax.set_ylabel("z (m)")
    fig.colorbar(m, ax=ax, label="u (m/s)", shrink=0.8)
    fig.savefig(sys.argv[2], dpi=140)


if __name__ == "__main__":
    main()
