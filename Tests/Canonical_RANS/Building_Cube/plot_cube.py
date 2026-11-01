#!/usr/bin/env python3
"""Figure of the flow around the block, for the docs and the PR. Not a test.

    plot_cube.py PLOTFILE OUT.png

u on the centre plane (x-z) and 2.5 m above the ground (x-y), with velocity arrows; the block is grey.
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
    hdr, d = erf_plotfile.read_fields(sys.argv[1], ["u", "v", "w", "building_beta"])
    u, v, w, b = [np.array(d[n], dtype=float) for n in ("u", "v", "w", "building_beta")]
    dx, dy, dz = hdr["dx"]
    nx, ny, nz = u.shape
    x, y, z = [(np.arange(n) + 0.5) * s for n, s in zip((nx, ny, nz), (dx, dy, dz))]
    jc, kmax = ny // 2, 14
    norm = TwoSlopeNorm(vmin=min(-1.0, u.min()), vcenter=0.0, vmax=u.max())
    fig, ax = plt.subplots(2, 1, figsize=(9, 6.4), constrained_layout=True)

    def panel(a, h, f, g, s, xlab, title, every):
        m = a.pcolormesh(h, f, np.ma.masked_where(s.T >= 1, g[0].T), cmap="RdBu_r", norm=norm, shading="nearest")
        a.pcolormesh(h, f, np.ma.masked_where(s.T < 1, s.T), cmap="Greys", vmin=0, vmax=2, shading="nearest")
        p, q = np.ma.masked_where(s.T >= 1, g[0].T), np.ma.masked_where(s.T >= 1, g[1].T)
        a.quiver(h[::every], f[::every], p[::every, ::every], q[::every, ::every], color="0.15", width=0.0022, scale=130)
        a.set_aspect("equal")
        a.set_xlabel("x (m)")
        a.set_ylabel(xlab)
        a.set_title(title, loc="left", fontsize=11)
        return m

    panel(ax[0], x, z[:kmax], (u[:, jc, :kmax], w[:, jc, :kmax]), b[:, jc, :kmax], "z (m)", "Centre plane", 1)
    m = panel(ax[1], x, y, (u[:, :, 0], v[:, :, 0]), b[:, :, 0], "y (m)", "2.5 m above the ground", 1)
    fig.colorbar(m, ax=ax, label="u (m/s)", shrink=0.8)
    fig.savefig(sys.argv[2], dpi=140)


if __name__ == "__main__":
    main()
