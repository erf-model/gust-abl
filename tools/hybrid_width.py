#!/usr/bin/env python3
"""Recommend a range for terrain.hybrid_half_width from the terrain and the mesh, before running GUST.

    hybrid_width.py --deck inputs [--terrain-file terrain.xyz] [--max-cells 16] [--plot widths.png]
    hybrid_width.py --terrain terrain.xyz --domain X0 Y0 X1 Y1 --dx 50 [--dy 50] [--periodic x|y|xy]
    hybrid_width.py --agnesi H A ridge_x|ridge_y|round --domain X0 Y0 X1 Y1 --dx 64

The hybrid mesh (terrain.representation = hybrid, docs/immersed_terrain.rst) follows the envelope h_s of the
terrain h_t over a square of half-width w and immerses the rest, the gap h_t - h_s, as on the flat mesh. The tool
builds h_s for every width from 0 to --max-cells cells, as Source/HybridTerrain.cpp does (the minimum over the
square, then the mean of the minima over the square, then the smaller of that and h_t; the square is cut at a
non-periodic face and wraps across a periodic one), and reports for each width:

  - the largest slope of h_s (the mesh's lower face), and
  - the largest gap and the gap as a share of the terrain's relief (max - min of h_t after a least-squares plane
    is removed, so that a regional tilt does not count as relief).

Two rules turn these into a range. Both come from a sweep of w on two 2-D Agnesi ridges (inflow/outflow, dz 2 m,
dx 64 and 32 m, ground one layer above the lower face), each run against the terrain-following run on the same
mesh, in m/s at 5-200 m above the ground:

  1. Lower end (convergence). On the ridge with slope 1.0, widths of 0 and 1 cell did not converge (Poisson stop or
     no drop in 20000 iterations); 2 cells and more converged in about half the terrain-following iterations. On
     the ridge with slope 0.26 every width converged. So when the terrain's largest slope exceeds --steep-slope
     (0.5) the lower end is --min-cells (2) cells; otherwise it is 0 (the terrain-following mesh, which then needs
     no hybrid mesh).
  2. Upper end (close to terrain-following). With the largest gap at 0.16 of the relief the speeds stayed within
     about 0.25 m/s of terrain-following on the gentle ridge (crest at 10 m: -0.11 and -0.05 m/s) and the crest was
     0.66 m/s low on the steep one; with 0.38 the crest was 1 to 5 m/s low. So the upper end is the widest width whose
     largest gap is at most --gap-fraction (0.2, between the two; no width between them was run) of the relief.
     This measures agreement with terrain-following, not with measurements: where terrain-following itself is in
     error (a cliff, a separated lee) a larger difference need not be worse.

The suggested value is the lower end when it lies within the range. When the lower end exceeds the upper one no
width meets both rules: the tool says so, gives the mesh spacing at which --min-cells cells fit under the upper
end, and exits with status 2. These rules are calibrated on two analytic ridges only; treat the range as a starting
point and check a run against a terrain-following one where that converges.

Terrain from a file is interpolated linearly to the mesh nodes (GUST interpolates by inverse distance:
terrain.idw_*), and terrain.smoothing is not applied, so the numbers are close to, not equal to, GUST's. A deck
with terrain.type = erf_grid needs --terrain with the same points as x y z.
Exit status: 0 a range found, 2 no width meets both rules, 3 an input error.
"""
import argparse
import os
import sys

import numpy as np


def fail(msg):
    print(f"hybrid_width.py: {msg}", file=sys.stderr)
    sys.exit(3)


# --- the envelope (Source/HybridTerrain.cpp: lower_envelope) ------------------------------------------------


def _line_pass(a, axis, p, periodic, mean):
    """The running minimum or mean over [i - p, i + p] along axis: wrapped when periodic, else cut at the ends."""
    a = np.moveaxis(a, axis, 0)
    n = a.shape[0]
    if periodic and 2 * p + 1 > n:
        fail(f"the square ({2 * p + 1} nodes) is longer than the periodic line ({n} nodes)")
    out = np.empty_like(a)
    if mean:
        tot = np.zeros_like(a)
        cnt = np.zeros(n)
    else:
        out[...] = np.inf
    for s in range(-p, p + 1):
        if periodic:
            b = np.roll(a, -s, axis=0)
            ok = np.ones(n, bool)
        else:
            idx = np.arange(n) + s
            ok = (idx >= 0) & (idx < n)
            b = a[np.clip(idx, 0, n - 1)]
        if mean:
            tot[ok] += b[ok]
            cnt[ok] += 1
        else:
            out[ok] = np.minimum(out[ok], b[ok])
    if mean:
        out = tot / cnt.reshape((n,) + (1,) * (a.ndim - 1))
    return np.moveaxis(out, 0, axis)


def lower_envelope(h, px, py, perx=False, pery=False):
    """h_s on the nodes h[i, j]: px, py nodes on each side in x and y."""
    if px == 0 and py == 0:
        return h.copy()
    m = _line_pass(_line_pass(h, 0, px, perx, False), 1, py, pery, False)
    e = _line_pass(_line_pass(m, 0, px, perx, True), 1, py, pery, True)
    return np.minimum(e, h)


def max_slope(h, dx, dy):
    """The largest |grad h| from one-sided node differences (the steepest cell edge or its diagonal)."""
    s = 0.0
    if h.shape[0] > 1:
        s = max(s, float(np.abs(np.diff(h, axis=0)).max()) / dx)
    if h.shape[1] > 1:
        s = max(s, float(np.abs(np.diff(h, axis=1)).max()) / dy)
    if h.shape[0] > 1 and h.shape[1] > 1:
        gx = np.diff(h, axis=0)[:, :-1] / dx
        gy = np.diff(h, axis=1)[:-1, :] / dy
        s = max(s, float(np.sqrt(gx ** 2 + gy ** 2).max()))
    return s


def relief(xs, ys, h):
    """max - min of h after a least-squares plane is removed."""
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    A = np.column_stack([np.ones(h.size), X.ravel(), Y.ravel()])
    c, *_ = np.linalg.lstsq(A, h.ravel(), rcond=None)
    r = h - (A @ c).reshape(h.shape)
    return float(r.max() - r.min())


# --- the terrain at the nodes --------------------------------------------------------------------------------


def nodes(lo, hi, d, periodic):
    n = int(round((hi - lo) / d))
    if n < 1 or abs(lo + n * d - hi) > 1e-6 * max(1.0, abs(hi)):
        fail(f"the domain {lo:g} to {hi:g} is not a whole number of cells of {d:g}")
    return lo + d * np.arange(n if periodic else n + 1)


def agnesi(xs, ys, H, A, shape, centre):
    X, Y = np.meshgrid(xs - centre[0], ys - centre[1], indexing="ij")
    r2 = {"ridge_x": X ** 2, "ridge_y": Y ** 2, "round": X ** 2 + Y ** 2}[shape]
    return H / (1 + r2 / A ** 2)


def from_file(path, xs, ys):
    try:
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    except ImportError:
        fail("a terrain file needs scipy")
    if not os.path.isfile(path):
        fail(f"terrain file {path} not found")
    p = np.loadtxt(path, comments="#", ndmin=2)
    if p.shape[1] < 3 or len(p) < 3:
        fail(f"terrain file {path}: needs at least three x y z points")
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    h = LinearNDInterpolator(p[:, :2], p[:, 2])(X, Y)
    bad = np.isnan(h)
    if bad.any():
        h[bad] = NearestNDInterpolator(p[:, :2], p[:, 2])(X[bad], Y[bad])
    return h


def read_deck(path):
    keys = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            keys[k.strip()] = v.split()
    return keys


def site_from_args(a):
    """xs, ys, h, dx, dy, periodic flags and a label from the arguments."""
    if a.deck:
        if not os.path.isfile(a.deck):
            fail(f"--deck {a.deck} not found")
        k = read_deck(a.deck)
        try:
            n = [int(v) for v in k["mesh.n_cell"]]
            lo = [float(v) for v in k["mesh.prob_lo"]]
            hi = [float(v) for v in k["mesh.prob_hi"]]
        except (KeyError, ValueError):
            fail(f"--deck {a.deck}: needs mesh.n_cell, mesh.prob_lo and mesh.prob_hi")
        per = [int(v) for v in k.get("mesh.is_periodic", ["0", "0", "0"])]
        dx, dy = (hi[0] - lo[0]) / n[0], (hi[1] - lo[1]) / n[1]
        xs, ys = nodes(lo[0], hi[0], dx, per[0]), nodes(lo[1], hi[1], dy, per[1])
        ttype = k.get("terrain.type", ["flat"])[0]
        if ttype == "agnesi":
            h = agnesi(xs, ys, float(k["terrain.hill_height"][0]), float(k["terrain.hill_half_width"][0]),
                       k.get("terrain.hill_shape", ["round"])[0], (0.5 * (lo[0] + hi[0]), 0.5 * (lo[1] + hi[1])))
        elif ttype in ("xyz", "erf_grid"):
            f = a.terrain_file or (k.get("terrain.file", [None])[0])
            if f is None:
                fail(f"--deck {a.deck}: terrain.file is not in the deck; give --terrain-file")
            if ttype == "erf_grid" and not a.terrain_file:
                fail("terrain.type = erf_grid: give the same points as x y z with --terrain-file")
            if not a.terrain_file and not os.path.isabs(f):
                f = os.path.join(os.path.dirname(os.path.abspath(a.deck)), f)
            h = from_file(f, xs, ys)
        else:
            fail(f"--deck {a.deck}: terrain.type = {ttype} has no terrain to immerse")
        return xs, ys, h, dx, dy, bool(per[0]), bool(per[1]), os.path.basename(a.deck)
    if not a.domain or not a.dx:
        fail("give --deck, or --terrain or --agnesi with --domain and --dx")
    x0, y0, x1, y1 = a.domain
    dx, dy = a.dx, (a.dy or a.dx)
    perx, pery = "x" in a.periodic, "y" in a.periodic
    xs, ys = nodes(x0, x1, dx, perx), nodes(y0, y1, dy, pery)
    if a.agnesi:
        H, A, shape = a.agnesi
        if shape not in ("ridge_x", "ridge_y", "round"):
            fail(f"--agnesi shape must be ridge_x, ridge_y or round (got {shape})")
        h = agnesi(xs, ys, float(H), float(A), shape, (0.5 * (x0 + x1), 0.5 * (y0 + y1)))
        label = f"Agnesi {shape} H {H} A {A}"
    elif a.terrain:
        h = from_file(a.terrain, xs, ys)
        label = os.path.basename(a.terrain)
    else:
        fail("give --terrain or --agnesi with --domain")
    return xs, ys, h, dx, dy, perx, pery, label


# --- the recommendation --------------------------------------------------------------------------------------


def sweep(xs, ys, h, dx, dy, perx, pery, max_cells):
    """One row per width w = n d, n = 0..max_cells: (n, w, slope of h_s, largest gap)."""
    d = cell(xs, ys, dx, dy, perx, pery)
    rows = []
    for n in range(max_cells + 1):
        w = n * d
        px, py = int(np.floor(w / dx + 1e-9)), int(np.floor(w / dy + 1e-9))
        if len(xs) == 1:
            px = 0
        if len(ys) == 1:
            py = 0
        if (perx and 2 * px + 1 > len(xs)) or (pery and 2 * py + 1 > len(ys)):
            break
        e = lower_envelope(h, px, py, perx, pery)
        rows.append((n, w, max_slope(e, dx, dy), float((h - e).max())))
    return rows


def cell(xs, ys, dx, dy, perx, pery):
    """The width step, as GUST's default counts a cell: the larger spacing of the directions with more than one
    cell (a periodic direction has one node per cell, the others one more)."""
    nx, ny = len(xs) - (0 if perx else 1), len(ys) - (0 if pery else 1)
    return dx if ny <= 1 else dy if nx <= 1 else max(dx, dy)


def recommend(rows, slope_t, rel, steep_slope, gap_fraction, min_cells):
    """(lower n, upper n or None, suggested n or None)."""
    lo = min_cells if slope_t > steep_slope else 0
    ok = [r[0] for r in rows if r[3] <= gap_fraction * rel + 1e-9]
    hi = max(ok) if ok else None
    return lo, hi, (lo if hi is not None and lo <= hi else None)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--deck", help="a GUST deck (mesh.*, terrain.type and its keys)")
    p.add_argument("--terrain-file", help="the terrain file, when it is not in the deck or overrides it")
    p.add_argument("--terrain", help="x y z points (without --deck)")
    p.add_argument("--agnesi", nargs=3, metavar=("H", "A", "SHAPE"), help="an Agnesi hill (without --deck)")
    p.add_argument("--domain", nargs=4, type=float, metavar=("X0", "Y0", "X1", "Y1"))
    p.add_argument("--dx", type=float)
    p.add_argument("--dy", type=float)
    p.add_argument("--periodic", default="", help="x, y or xy")
    p.add_argument("--max-cells", type=int, default=16)
    p.add_argument("--steep-slope", type=float, default=0.5)
    p.add_argument("--gap-fraction", type=float, default=0.2)
    p.add_argument("--min-cells", type=int, default=2)
    p.add_argument("--plot", help="a PNG of h_t and h_s for the suggested width")
    a = p.parse_args(argv)
    for name in ("max_cells", "min_cells"):
        if getattr(a, name) < 0:
            fail(f"--{name.replace('_', '-')} must be >= 0")
    for name in ("steep_slope", "gap_fraction"):
        if not np.isfinite(getattr(a, name)) or getattr(a, name) <= 0:
            fail(f"--{name.replace('_', '-')} must be > 0")
    if a.dx is not None and not a.dx > 0:
        fail("--dx must be > 0")
    if a.dy is not None and not a.dy > 0:
        fail("--dy must be > 0")
    if set(a.periodic) - set("xy"):
        fail(f"--periodic must be x, y or xy (got {a.periodic})")

    xs, ys, h, dx, dy, perx, pery, label = site_from_args(a)
    st, rel = max_slope(h, dx, dy), relief(xs, ys, h)
    if rel <= 0:
        fail("the terrain has no relief after the plane is removed: nothing to immerse")
    rows = sweep(xs, ys, h, dx, dy, perx, pery, a.max_cells)
    lo, hi, sug = recommend(rows, st, rel, a.steep_slope, a.gap_fraction, a.min_cells)
    d = cell(xs, ys, dx, dy, perx, pery)

    print(f"{label}: dx {dx:g} m, dy {dy:g} m; terrain slope up to {st:.2f}, relief {rel:.1f} m (plane removed)")
    print(f"{'cells':>5} {'w (m)':>8} {'h_s slope':>9} {'gap max (m)':>11} {'gap/relief':>10}  note")
    for n, w, s, g in rows:
        note = []
        if n < lo:
            note.append("below the lower end")
        if g > a.gap_fraction * rel + 1e-9:
            note.append("gap above the upper end")
        print(f"{n:5d} {w:8.1f} {s:9.2f} {g:11.1f} {g / rel:10.3f}  {', '.join(note)}")
    print()
    print(f"Lower end: {lo} cells ({lo * d:g} m)" + (f": the terrain's slope exceeds {a.steep_slope:g}" if lo else
          f": the slope is at most {a.steep_slope:g}, so terrain-following (terrain.representation = body_fitted, "
          f"or hybrid with 0) should converge"))
    if hi is None:
        print(f"Upper end: none: even w = 0 has a gap above {a.gap_fraction:g} of the relief")
    else:
        print(f"Upper end: {hi} cells ({hi * d:g} m): the largest gap is at most {a.gap_fraction:g} of the relief"
              + (" (the widest width tried; raise --max-cells to look further)" if hi == rows[-1][0] else ""))
    if sug is not None:
        print(f"Suggested: terrain.hybrid_half_width = {sug * d:g}   (range {lo * d:g} to {hi * d:g} m)")
        status = 0
    else:
        # The width in metres that the upper end allows, found on this mesh; min_cells cells fit under it at
        # dx = w_hi / min_cells.
        w_hi = hi * d if hi is not None else 0.0
        need = w_hi / a.min_cells if a.min_cells else float("nan")
        print(f"No width meets both rules: the lower end ({lo * d:g} m) is above the upper end ({w_hi:g} m).")
        if w_hi > 0:
            print(f"  {a.min_cells} cells fit under the upper end at about dx = {need:g} m (here {d:g} m). On this mesh")
        else:
            print("  On this mesh")
        print(f"  terrain.hybrid_half_width = {lo * d:g} is the narrowest width that converged in the calibration; expect"
              f" speeds near crests and edges more than 0.25 m/s from a terrain-following run.")
        status = 2
    if a.plot:
        plot(a.plot, xs, ys, h, dx, dy, perx, pery, (sug if sug is not None else lo) * d)
    return status


def plot(path, xs, ys, h, dx, dy, perx, pery, w):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    e = lower_envelope(h, int(np.floor(w / dx + 1e-9)) if len(xs) > 1 else 0,
                       int(np.floor(w / dy + 1e-9)) if len(ys) > 1 else 0, perx, pery)
    j = int(np.argmax(h.max(axis=0)))
    fig, ax = plt.subplots(figsize=(7, 3), constrained_layout=True)
    ax.plot(xs, h[:, j], color="0.2", lw=1.5, label="terrain $h_t$")
    ax.plot(xs, e[:, j], color="#1f6fb4", lw=1.5, label=f"envelope $h_s$, w = {w:g} m")
    ax.fill_between(xs, e[:, j], h[:, j], color="#1f6fb4", alpha=0.2, lw=0, label="immersed gap")
    ax.set_xlabel("x (m)")
    ax.set_ylabel("height (m)")
    ax.set_title(f"Section through the highest row, y = {ys[j]:g} m", fontsize=10)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(path, dpi=150)


if __name__ == "__main__":
    sys.exit(main())
