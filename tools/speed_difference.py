"""Wind-speed differences against a reference run, in m/s: counts above thresholds and a map on the terrain.

Used by the nest cases (Tests/Real_Terrain/Askervein_Nest/askervein_nest.py two_way,
Tests/Real_Terrain/Columbia_Gorge/columbia_gorge.py compare) to answer "how many points differ by
more than 0.25, 0.5 or 1 m/s, and where".

    at_height(zn, zr, fr)     the reference field fr (columns at heights zr) at the heights zn, column by column
    coarsen(a, r)             the r x r mean of a in x and y (a fine field onto the coarse cells)
    near_ground(z, agl)       the cells whose centres are less than agl above their column's first centre
    counts(d, thresholds)     cells of |d| above each threshold
    print_table(rows, ...)    one line per (run, region): cells, largest |d|, counts above the thresholds
    column_max(d, mask)       the largest |d| of each column over the masked cells (0 where none)
    plot_map(ax, ...)         terrain contours with each column coloured by its largest |d| (matplotlib)

The speed is the horizontal one, sqrt(u^2 + v^2), at the cell centres. A joint nest's levels differ from the
reference's, so a nest is compared at its own cell heights (at_height), never by index.
"""
import numpy as np

THRESHOLDS = (0.25, 0.5, 1.0)   # m/s
AGL = 200.0                     # m above the first cell centre: the layer the tables and maps cover
COLOURS = ("#f2c94c", "#e07b00", "#c0392b")


def speed(f):
    return np.hypot(f["u"], f["v"])


def at_height(zn, zr, fr):
    out = np.empty_like(zn, dtype=float)
    for a in range(zn.shape[0]):
        for b in range(zn.shape[1]):
            out[a, b] = np.interp(zn[a, b], zr[a, b], fr[a, b])
    return out


def coarsen(a, r=2):
    nx, ny = a.shape[0] // r, a.shape[1] // r
    return a[:nx * r, :ny * r].reshape(nx, r, ny, r, *a.shape[2:]).mean(axis=(1, 3))


def near_ground(z, agl=AGL):
    return (z - z[:, :, :1]) < agl


def counts(d, thresholds=THRESHOLDS):
    d = np.abs(np.asarray(d, dtype=float))
    return [int(np.count_nonzero(d > t)) for t in thresholds]


def print_table(rows, thresholds=THRESHOLDS, out=print):
    """rows: (run, region, differences); the differences are the cells to count (already masked)."""
    head = f"{'run':16s} {'region':8s} {'cells':>8s} {'max |dS|':>9s} " + " ".join(f"{'> %g m/s' % t:>16s}" for t in thresholds)
    out(head)
    for run, region, d in rows:
        d = np.abs(np.asarray(d, dtype=float)).ravel()
        n = max(d.size, 1)
        cols = " ".join(f"{c:7d} ({100.0 * c / n:5.2f}%)" for c in counts(d, thresholds))
        out(f"{run:16s} {region:8s} {d.size:8d} {d.max() if d.size else 0.0:9.3f} {cols}")


def column_max(d, mask):
    return np.where(mask, np.abs(d), 0.0).max(axis=2)


def plot_map(ax, terrain, layers, box, title, thresholds=THRESHOLDS, colours=COLOURS):
    """terrain: (x, y, h) with h[i, j] at (x[i], y[j]); layers: (X, Y, dmax) arrays of column centres and
    their largest |d| (NaN: not plotted); box: (x0, y0, width, height) of the nest. Returns the counts of
    columns in each band [t_n, t_n+1)."""
    import matplotlib.pyplot as plt
    x, y, h = terrain
    cs = ax.contour(x, y, h.T, levels=15, colors="#777777", linewidths=0.6)
    ax.clabel(cs, fontsize=6, fmt="%.0f")
    ax.add_patch(plt.Rectangle(box[:2], box[2], box[3], fill=False, ec="#2a6fdb", lw=1.4, ls="--"))
    bands = list(thresholds) + [np.inf]
    n = []
    for b in range(len(thresholds)):
        k = 0
        for X, Y, dmax in layers:
            m = (dmax > bands[b]) & (dmax <= bands[b + 1])
            k += int(np.count_nonzero(m))
            ax.scatter(X[m], Y[m], s=10 + 6 * b, c=colours[b], edgecolors="none", zorder=3)
        n.append(k)
    ax.set_title(title, fontsize=9)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("y (m)")
    return n


def legend(fig, thresholds=THRESHOLDS, colours=COLOURS):
    import matplotlib.pyplot as plt
    labels = [f"{a:g}-{b:g} m/s" for a, b in zip(thresholds, thresholds[1:])] + [f"> {thresholds[-1]:g} m/s"]
    h = [plt.Line2D([], [], marker="o", ls="", color=c, label=lab) for c, lab in zip(colours, labels)]
    h.append(plt.Line2D([], [], color="#2a6fdb", ls="--", label="nest"))
    fig.legend(handles=h, loc="lower center", ncol=len(h), frameon=False)
