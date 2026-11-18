#!/usr/bin/env python3
"""Check a nest's placement from the terrain alone, before running GUST.

    nest_placement.py --terrain terrain.xyz --domain X0 Y0 X1 Y1 --dx 200 \\
                      --nest X0 Y0 X1 Y1 [--nest ...] [--ratio 2] [--direction 270 | D1 D2 ... | all] [--plot map.png]
    nest_placement.py --terrain terrain.xyz --domain X0 Y0 X1 Y1 --dx 200 --recommend WIDTH HEIGHT [--plot map.png]

--nest checks given nests; --recommend slides a nest of the given size over the outer mesh, drops the
positions with a flagged inflow or outflow face for any of the wind directions and lists the rest by the
share of high-misfit columns they cover (corners on outer nodes, ready for nest.<name>.lo / hi). The
domain, the nest corners and the --recommend sizes are whole outer cells, as GUST needs them.

Two questions, answered from the terrain file and the outer mesh only:

1. Is a nest worth adding here?  The outer mesh sees the ground only at its nodes, dx apart. The misfit
   |h - h_c| between the terrain h and the bilinear surface h_c through the outer nodes is where the outer
   cells are wrong: on the Columbia Gorge the outer error follows it (Spearman +0.55; the column curvature
   |lap h| +0.69). The tool reports the share of high-misfit columns (misfit > --misfit) that each nest
   covers.

2. Are the nest's faces in a good place?  A nest's error enters through its faces, and in a joint run
   (nest.feedback = defect) the fed-back flow leaves through the downstream faces into outer cells that
   cannot carry it. For each face the tool reports:
     - misfit along the face (95th percentile, m);
     - each channel or ridge the face cuts (a local minimum or maximum of h along the face with
       prominence > --prominence) and each step (a change of more than --step-rise across --step-width
       along the face, outside the channels and ridges: an escarpment, a plateau edge, a steep valley side),
       with its size, the angle between its axis and the face normal (0 = it runs straight through the
       face; for a step the axis is the contour line) and the misfit within one outer cell of it;
     - how far a channel's or ridge's axis turns within --bend-length of the face (a note, not a flag);
     - the misfit of the outer cells beyond the face (up to --downstream m), which carry the fed-back
       flow to the outflow when the face is an outflow face.
   An inflow or outflow face is flagged when it cuts a channel, ridge or step that the outer cells cannot
   follow (misfit within one outer cell above --misfit) at more than --angle degrees to the face normal.
   For each wind direction (--direction; "all" = every 22.5 degrees) a face is inflow or outflow when
   the wind is within 60 degrees of its normal, and a side face otherwise; side faces are reported, not
   flagged, and a nest is clear only when it is clear for every direction. The angle is measured against
   the face normal, not the wind.

   What the runs support (PR-70, a channel site with 200 m outer cells, the wind normal to the faces):
   the 20 degree angle and the 5 m misfit gate were fitted on four nests (25 and 28 degrees with misfit
   above 5 m raised the outer cells above 0.25 m/s from 12 to 663 and 817 in a joint run; 0 and 14 degrees
   stayed at 12 and 26); they have not been tested on an independent case. A nest next to an 89 degree
   channel bend stayed at 12 cells, so a bend is a note. The step rule was tested on an escarpment with
   50 m cells (Tests/Canonical_RANS/Nest_Placement_Runs, N-114 item 11): for winds along an axis the two
   flagged nests made the joint run worse and the two clear ones did not; for an oblique wind blowing
   within 10 degrees of the step's contours both nests were flagged but neither run was worse (the
   angle is taken against the face normal, not the wind). Other cell sizes have not been tested by a
   run. A nest closer than --margin to an outer face gets a note: on the Gorge a
   nest 4 outer cells from the inflow face did not converge in a joint run.

The thresholds scale with the outer cell size dx; the factors are the values fitted at dx = 200 m:
   --misfit      0.025 dx (5 m)              --prominence  0.125 dx (25 m)
   --step-rise   0.125 dx (25 m)             --step-width  2 dx
   --bend-length 5 dx                        --downstream  15 dx
   --margin      10 dx
Each can be given in m instead. The misfit is the outer mesh's error in the ground height, so it is
compared with a length tied to dx, not to the terrain's relief (a term in the relief would make the verdict
depend on how far the domain reaches).

The terrain file holds x y z points (any order, '#' comments), as GUST's terrain.type = xyz reads them;
it is interpolated (linear) onto a regular grid of spacing dx / (2 ratio).

Exit status: 0 when every nest checked is clear (or --recommend found a clear position), 1 when a checked
nest has a flagged inflow or outflow face, 2 when --recommend finds no clear position, 3 for an input
error (the message names the input) or when no position fits.
"""
import argparse
import math
import sys

import numpy as np

ALL_DIRECTIONS = [22.5 * k for k in range(16)]
PER_DX = {"misfit": 0.025, "prominence": 0.125, "step_rise": 0.125, "step_width": 2.0, "bend_length": 5.0,
          "downstream": 15.0, "margin": 10.0}
LENGTHS = ("misfit", "prominence", "step_rise", "step_width", "bend_length", "downstream", "margin")
INPUT_ERROR = 3


def on_cells(v, o, dx):
    """True when v - o is a whole number of outer cells."""
    q = (v - o) / dx
    return abs(q - round(q)) < 1e-6


def load_terrain(path, domain, step):
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    p = np.loadtxt(path, comments="#", ndmin=2)
    x0, y0, x1, y1 = domain
    if p.shape[1] < 3 or len(p) < 3:
        fail(f"--terrain {path}: needs at least three x y z points")
    tol = 0.5 * step
    if p[:, 0].min() > x0 + tol or p[:, 0].max() < x1 - tol or p[:, 1].min() > y0 + tol or p[:, 1].max() < y1 - tol:
        fail(f"--domain {x0:g} {y0:g} {x1:g} {y1:g} reaches beyond the points of --terrain {path} "
             f"(x {p[:, 0].min():g} to {p[:, 0].max():g}, y {p[:, 1].min():g} to {p[:, 1].max():g})")
    xs = np.arange(x0, x1 + 0.5 * step, step)
    ys = np.arange(y0, y1 + 0.5 * step, step)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    h = LinearNDInterpolator(p[:, :2], p[:, 2])(X, Y)
    if np.isnan(h).any():
        h[np.isnan(h)] = NearestNDInterpolator(p[:, :2], p[:, 2])(X[np.isnan(h)], Y[np.isnan(h)])
    return xs, ys, h


def coarse_surface(xs, ys, h, dx):
    """h_c: the bilinear surface through the terrain at the outer mesh's nodes, on the fine grid."""
    from scipy.interpolate import RegularGridInterpolator
    step = xs[1] - xs[0]
    m = int(round(dx / step))
    f = RegularGridInterpolator((xs[::m], ys[::m]), h[::m, ::m], bounds_error=False, fill_value=None)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    return f(np.stack([X, Y], axis=-1))


def derivatives(h, step, smooth):
    """Gradient and Hessian of h smoothed over `smooth` m."""
    from scipy.ndimage import gaussian_filter
    hs = gaussian_filter(h, smooth / step) if smooth > 0 else h
    hx, hy = np.gradient(hs, step)
    hxx, hxy = np.gradient(hx, step)
    _, hyy = np.gradient(hy, step)
    return (hx, hy), (hxx, hxy, hyy)


def axis_at(H, i, j):
    """Unit vector along a channel or ridge at (i, j): the eigenvector of the Hessian with the
    smaller |eigenvalue| (the direction the terrain curves least)."""
    hxx, hxy, hyy = (a[i, j] for a in H)
    w, v = np.linalg.eigh(np.array([[hxx, hxy], [hxy, hyy]]))
    return v[:, int(np.argmin(np.abs(w)))]


def contour_at(G, i, j):
    """Unit vector along the contour line at (i, j) (normal to the gradient): a step's axis."""
    gx, gy = G[0][i, j], G[1][i, j]
    g = math.hypot(gx, gy)
    return np.array([-gy / g, gx / g]) if g > 0 else np.array([1.0, 0.0])


def trace_turn(H, xs, ys, i, j, length):
    """Largest change of the feature axis direction (degrees) over +-length along the axis from (i, j)."""
    step = xs[1] - xs[0]
    a0 = axis_at(H, i, j)
    worst = 0.0
    for sgn in (1.0, -1.0):
        x, y = xs[i], ys[j]
        a = sgn * a0
        for _ in range(int(length / step)):
            x, y = x + step * a[0], y + step * a[1]
            ii, jj = int(round((x - xs[0]) / step)), int(round((y - ys[0]) / step))
            if not (0 <= ii < len(xs) and 0 <= jj < len(ys)):
                break
            b = axis_at(H, ii, jj)
            if np.dot(a, b) < 0:
                b = -b
            worst = max(worst, math.degrees(math.acos(min(1.0, abs(float(np.dot(a0, b)))))))
            a = b
    return worst


FACES = (("west", 0, 1, 0, 3, (-1, 0)), ("east", 2, 1, 2, 3, (1, 0)),
         ("south", 0, 1, 2, 1, (0, -1)), ("north", 0, 3, 2, 3, (0, 1)))


def face_kind(normal, direction):
    """inflow, outflow or side for the wind from `direction` (degrees from north the wind comes from)."""
    th = math.radians(direction)
    c = -math.sin(th) * normal[0] - math.cos(th) * normal[1]   # flow (towards) . outward normal
    return "outflow" if c > 0.5 else ("inflow" if c < -0.5 else "side")


def crossings(prof, prominence):
    """Channels and ridges along a face profile: [(index, type, prominence m)] and their bases [(lo, hi)]
    (the stretch of the profile down to the prominence's reference level)."""
    from scipy.signal import find_peaks, peak_widths
    found, bases = [], []
    for sgn, what in ((1, "ridge"), (-1, "channel")):
        idx, props = find_peaks(sgn * prof, prominence=prominence)
        found += [(int(k), what, float(p)) for k, p in zip(idx, props["prominences"])]
        if len(idx):
            _, _, lo, hi = peak_widths(sgn * prof, idx, rel_height=1.0, prominence_data=(
                props["prominences"], props["left_bases"], props["right_bases"]))
            bases += list(zip(lo, hi))
    return sorted(found), bases


def steps(prof, w, rise, bases):
    """Steps along a face profile: (index, 'step', change m) where h changes by more than `rise` across
    +-w points, outside the base of every channel and ridge (their flanks are reported with them)."""
    from scipy.signal import find_peaks
    n = len(prof)
    if n < 2 * w + 1 or w < 1:
        return []
    d = np.zeros(n)
    d[w:n - w] = np.abs(prof[2 * w:] - prof[:n - 2 * w])
    idx, _ = find_peaks(d, height=rise, distance=2 * w)
    return [(int(k), "step", float(d[k])) for k in idx if not any(lo <= k <= hi for lo, hi in bases)]


class Site:
    """The terrain on a fine grid, the outer mesh's misfit, the thresholds and the derivatives."""

    def __init__(self, xs, ys, h, args):
        self.xs, self.ys, self.h = xs, ys, h
        self.step = xs[1] - xs[0]
        self.relief = float(h.max() - h.min())
        resolve_thresholds(args)
        self.args = args
        self.misfit_map = np.abs(h - coarse_surface(xs, ys, h, args.dx))
        self.G, self.H = derivatives(h, self.step, args.dx / 2)
        self.X, self.Y = np.meshgrid(xs, ys, indexing="ij")
        self.high = self.misfit_map > args.misfit

    @classmethod
    def from_file(cls, args):
        xs, ys, h = load_terrain(args.terrain, args.domain, args.dx / (2 * args.ratio))
        return cls(xs, ys, h, args)

    def faces(self, nest):
        """One dict per face: name, normal, misfit95, beyond95, features (text), faults (text: the cuts
        that flag the face when it is an inflow or outflow face) and notes (bends)."""
        a_, xs, ys, h, misfit_map, step = self.args, self.xs, self.ys, self.h, self.misfit_map, self.step
        X, Y = self.X, self.Y
        r = int(round(a_.dx / step))
        w = max(1, int(round(0.5 * a_.step_width / step)))
        out = []
        for name, ia, ja, ib, jb, n in FACES:
            a, b, n = np.array((nest[ia], nest[ja]), float), np.array((nest[ib], nest[jb]), float), np.array(n, float)
            ln = float(np.linalg.norm(b - a))
            s = np.arange(0.0, ln + 0.5 * step, step)
            px, py = a[0] + (b[0] - a[0]) * s / ln, a[1] + (b[1] - a[1]) * s / ln
            ii = np.clip(np.rint((px - xs[0]) / step).astype(int), 0, len(xs) - 1)
            jj = np.clip(np.rint((py - ys[0]) / step).astype(int), 0, len(ys) - 1)
            d = (X - a[0]) * n[0] + (Y - a[1]) * n[1]
            tt = (X - a[0]) * (b[0] - a[0]) / ln + (Y - a[1]) * (b[1] - a[1]) / ln
            band = (d > 0) & (d <= a_.downstream) & (tt >= 0) & (tt <= ln)
            prof = h[ii, jj]
            peaks, bases = crossings(prof, a_.prominence)
            found = peaks + steps(prof, w, a_.step_rise, bases)
            feats, faults, notes = [], [], []
            for k, what, size in sorted(found):
                ax = contour_at(self.G, ii[k], jj[k]) if what == "step" else axis_at(self.H, ii[k], jj[k])
                ang = math.degrees(math.acos(min(1.0, abs(float(np.dot(ax, n))))))
                local = float(misfit_map[max(ii[k] - r, 0):ii[k] + r + 1, max(jj[k] - r, 0):jj[k] + r + 1].max())
                turn = None if what == "step" else trace_turn(self.H, xs, ys, ii[k], jj[k], a_.bend_length)
                feats.append(f"{what} {size:.0f} {ang:.0f} {'-' if turn is None else f'{turn:.0f}'} {local:.1f}")
                # Only a feature the outer cells cannot follow (local misfit above --misfit) can be
                # carried wrongly across the face.
                if local <= a_.misfit:
                    continue
                if ang > a_.angle:
                    faults.append(f"cuts a {what} at {ang:.0f} deg (misfit {local:.0f} m)")
                if turn is not None and turn > a_.turn:
                    notes.append(f"{what} bends {turn:.0f} deg near the face (misfit {local:.0f} m)")
            out.append(dict(name=name, normal=n, m95=float(np.percentile(misfit_map[ii, jj], 95)),
                            beyond=float(np.percentile(misfit_map[band], 95)) if band.any() else None,
                            feats=feats, faults=faults, notes=notes))
        return out

    def evaluate(self, nest):
        """cover, misfit inside / outside, faces, flagged {direction: [(face, kind, faults)]}, edge notes."""
        a_ = self.args
        x0, y0, x1, y1 = nest
        inside = (self.X >= x0) & (self.X <= x1) & (self.Y >= y0) & (self.Y <= y1)
        cover = self.high[inside].sum() / max(self.high.sum(), 1)
        faces = self.faces(nest)
        flagged = {}
        for dr in a_.direction:
            bad = [(f["name"], k, f["faults"]) for f in faces
                   for k in [face_kind(f["normal"], dr)] if k != "side" and f["faults"]]
            if bad:
                flagged[dr] = bad
        X0, Y0, X1, Y1 = a_.domain
        gaps = {"west": x0 - X0, "east": X1 - x1, "south": y0 - Y0, "north": Y1 - y1}
        edge = [f"{nm} face {g:g} m ({g / a_.dx:.0f} outer cells) from the outer face, under --margin {a_.margin:g} m"
                for nm, g in gaps.items() if g < a_.margin - 1e-6]
        m_in = float(self.misfit_map[inside].mean())
        m_out = float(self.misfit_map[~inside].mean()) if (~inside).any() else 0.0
        return dict(nest=nest, cover=float(cover), m_in=m_in, m_out=m_out, faces=faces, flagged=flagged, edge=edge)


def resolve_thresholds(a):
    """Fill the thresholds not given from dx (PER_DX); record which were scaled."""
    a.scaled = []
    for name, f in PER_DX.items():
        if getattr(a, name) is None:
            setattr(a, name, f * a.dx)
            a.scaled.append(name)


def summary(site):
    a, h, misfit_map = site.args, site.h, site.misfit_map
    print(f"terrain: relief {site.relief:.0f} m; outer cells {a.dx:g} m see it with misfit |h - h_c| mean "
          f"{misfit_map.mean():.1f} m, 95th percentile {np.percentile(misfit_map, 95):.1f} m, max {misfit_map.max():.1f} m")
    print("thresholds (* from dx): " + ", ".join(
        f"{nm.replace('_', '-')} {getattr(a, nm):g} m{'*' if nm in a.scaled else ''}" for nm in LENGTHS)
        + f", angle {a.angle:g} deg, turn {a.turn:g} deg")
    print(f"wind from: {directions_text(a.direction)}")
    print(f"columns with misfit > {a.misfit:g} m: {100 * site.high.mean():.1f}% of the domain")
    if not site.high.any():
        print("  no column exceeds the misfit threshold: the outer mesh follows this terrain, a nest adds little here")


def directions_text(ds):
    return "all 16 (every 22.5 deg)" if list(ds) == ALL_DIRECTIONS else " ".join(f"{d:g}" for d in ds)


def verdict(ev, ndir):
    if not ev["flagged"]:
        return "faces look clear" + (f" for all {ndir} wind directions" if ndir > 1 else "")
    parts = []
    for dr, bad in ev["flagged"].items():
        parts.append(f"wind from {dr:g}: " + "; ".join(f"{nm} ({kd}) {', '.join(fl)}" for nm, kd, fl in bad))
    return f"move the nest: flagged for {len(ev['flagged'])} of {ndir} wind directions: " + " | ".join(parts)


def report(ev, args):
    x0, y0, x1, y1 = ev["nest"]
    one = len(args.direction) == 1
    print(f"\nnest {x0:g} {y0:g} {x1:g} {y1:g}: covers {100 * ev['cover']:.0f}% of the high-misfit columns "
          f"(misfit inside: mean {ev['m_in']:.1f} m, outside {ev['m_out']:.1f} m)")
    print(f"  {'face':6s} {'kind' if one else '':8s} {'misfit95':>8s} {'beyond95':>8s}  "
          f"features cut (type, size m, angle deg, turn deg, misfit m)  faults / notes")
    for f in ev["faces"]:
        kind = face_kind(f["normal"], args.direction[0]) if one else ""
        b = f"{f['beyond']:8.1f}" if f["beyond"] is not None and (not one or kind == "outflow") else ""
        text = "; ".join(f["faults"] + [f"note: {t}" for t in f["notes"]]) or "ok"
        if one and kind == "side" and f["faults"]:
            text += " (side face: not flagged)"
        print(f"  {f['name']:6s} {kind:8s} {f['m95']:8.1f} {b:>8s}  {', '.join(f['feats']) if f['feats'] else '-':52s}  {text}")
    for e in ev["edge"]:
        print(f"  note: {e}")
    print("  verdict: " + verdict(ev, len(args.direction)))


def check(args, site=None):
    site = site or Site.from_file(args)
    summary(site)
    results = [site.evaluate(nest) for nest in args.nest]
    for ev in results:
        report(ev, args)
    if args.plot:
        plot(args, site, results)
    return 1 if any(ev["flagged"] for ev in results) else 0


def recommend(args, site=None):
    """Slide a nest of the given size over the outer mesh (corners on outer nodes, --search-step apart,
    at least --margin from the outer faces); keep the positions with no flagged inflow or outflow face for
    any wind direction and rank them by the share of high-misfit columns covered, then by the largest face
    misfit (smaller first). Returns (exit status, clear positions)."""
    site = site or Site.from_file(args)
    summary(site)
    w, hgt = args.recommend
    X0, Y0, X1, Y1 = args.domain
    st = args.search_step if args.search_step is not None else 4 * args.dx
    snap = lambda v, o: o + args.dx * math.ceil((v - o) / args.dx - 1e-9)
    cands, seen = [], set()
    for x0 in np.arange(X0 + args.margin, X1 - args.margin - w + 0.5, st):
        for y0 in np.arange(Y0 + args.margin, Y1 - args.margin - hgt + 0.5, st):
            x0s, y0s = snap(x0, X0), snap(y0, Y0)
            if x0s + w > X1 - args.margin + 1e-6 or y0s + hgt > Y1 - args.margin + 1e-6 or (x0s, y0s) in seen:
                continue
            seen.add((x0s, y0s))
            ev = site.evaluate((x0s, y0s, x0s + w, y0s + hgt))
            face = max(f["m95"] for f in ev["faces"])
            cands.append((bool(ev["flagged"]), -round(ev["cover"], 2), face, ev))
    if not cands:
        fail("--recommend: no position fits: the nest plus --margin is larger than --domain")
    cands.sort(key=lambda c: (c[0], c[1], c[2]))
    ok = [c for c in cands if not c[0]]
    print(f"\n{len(cands)} positions of a {w:g} x {hgt:g} m nest tried ({st:g} m apart, {args.margin:g} m from the outer faces); "
          f"{len(ok)} have clear inflow and outflow faces for every wind direction given")
    flagged = [c for c in cands if c[0]]
    if ok:
        print(f"{'rank':4s} {'nest.<name>.lo':>16s} {'nest.<name>.hi':>16s} {'covers':>7s} {'face misfit95':>13s}")
        for k, c in enumerate(ok[:args.top]):
            x0, y0, x1, y1 = c[3]["nest"]
            print(f"{k + 1:4d} {x0:8g} {y0:7g} {x1:8g} {y1:7g} {100 * -c[1]:6.0f}% {c[2]:13.1f}")
    best_flagged = min(flagged, key=lambda c: (c[1], c[2])) if flagged else None
    if flagged:
        c = best_flagged
        x0, y0, x1, y1 = c[3]["nest"]
        print(f"best-covering flagged position, for comparison: {x0:g} {y0:g} {x1:g} {y1:g} covers {100 * -c[1]:.0f}%")
    if ok:
        report(ok[0][3], args)
        n = ok[0][3]["nest"]
        print(f"\nrecommended: nest.<name>.lo = {n[0]:g} {n[1]:g}   nest.<name>.hi = {n[2]:g} {n[3]:g}")
    else:
        msg = (f"WARNING: no position of a {w:g} x {hgt:g} m nest has clear inflow and outflow faces for wind from "
               f"{directions_text(args.direction)}. Nothing is recommended. Try another size, fewer wind directions, a "
               f"smaller --search-step or --margin, or run nest.feedback = none against defect on the best-covering "
               f"flagged position above.")
        print("\n" + msg)
        print(msg, file=sys.stderr)
    if args.plot:
        plot(args, site, [c[3] for c in ok[:3]] if ok else [best_flagged[3]])
    return (0 if ok else 2), [c[3] for c in ok]


def plot(args, site, results):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xs, ys, h, misfit_map = site.xs, site.ys, site.h, site.misfit_map
    fig, ax = plt.subplots(figsize=(8, 6.5))
    im = ax.pcolormesh(xs, ys, misfit_map.T, cmap="magma_r", shading="auto",
                       vmin=0, vmax=max(2 * args.misfit, np.percentile(misfit_map, 99)))
    fig.colorbar(im, ax=ax, label="outer-mesh misfit |h - h_c| (m)")
    cs = ax.contour(xs, ys, h.T, levels=15, colors="#4a7f4a", linewidths=0.5)
    ax.clabel(cs, fontsize=6, fmt="%.0f")
    for ev in results:
        bad = {nm for lst in ev["flagged"].values() for nm, _, _ in lst}
        for name, ia, ja, ib, jb, _ in FACES:
            n = ev["nest"]
            ax.plot([n[ia], n[ib]], [n[ja], n[jb]], color="#c0392b" if name in bad else "#2a6fdb", lw=2.2)
    ax.set_aspect("equal")
    ax.set_title(f"nest faces: blue clear, red flagged (outer dx {args.dx:g} m, wind from "
                 f"{directions_text(args.direction)})", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=110)
    print(f"\nwrote {args.plot}")


def fail(msg):
    print(f"nest_placement: error: {msg}", file=sys.stderr)
    raise SystemExit(INPUT_ERROR)


class Parser(argparse.ArgumentParser):
    def error(self, message):
        self.print_usage(sys.stderr)
        fail(message)


def parse(argv):
    ap = Parser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--terrain", required=True, help="x y z points (GUST terrain.type = xyz)")
    ap.add_argument("--domain", type=float, nargs=4, required=True, metavar=("X0", "Y0", "X1", "Y1"), help="outer domain, m")
    ap.add_argument("--dx", type=float, required=True, help="outer cell size, m")
    ap.add_argument("--ratio", type=int, default=2, help="nest refinement ratio (default 2)")
    ap.add_argument("--nest", type=float, nargs=4, action="append", metavar=("X0", "Y0", "X1", "Y1"),
                    help="check a nest's corners (nest.<name>.lo / hi), m; repeatable")
    ap.add_argument("--recommend", type=float, nargs=2, metavar=("WIDTH", "HEIGHT"),
                    help="find places for a nest of this size (m, x and y) instead of checking --nest")
    ap.add_argument("--search-step", type=float, help="spacing of the positions tried, m (default 4 dx)")
    ap.add_argument("--margin", type=float, help="smallest distance from a nest to the outer faces, m (default 10 dx)")
    ap.add_argument("--top", type=int, default=5, help="positions listed (default 5)")
    ap.add_argument("--direction", nargs="+", default=["270"],
                    help="wind directions, degrees from north the wind comes from, or 'all' (every 22.5 deg); default 270")
    ap.add_argument("--misfit", type=float, help="high-misfit threshold, m (default 0.025 dx)")
    ap.add_argument("--prominence", type=float,
                    help="smallest channel depth or ridge height counted at a face, m (default 0.125 dx)")
    ap.add_argument("--step-rise", type=float,
                    help="smallest change of height across --step-width counted as a step, m (default 0.125 dx)")
    ap.add_argument("--step-width", type=float, help="distance along a face over which a step's rise is measured, m (default 2 dx)")
    ap.add_argument("--angle", type=float, default=20.0, help="flag a feature cut at more than this angle to the face normal, degrees (default 20)")
    ap.add_argument("--turn", type=float, default=20.0, help="note a channel or ridge whose axis turns more than this near the face, degrees (default 20)")
    ap.add_argument("--bend-length", type=float, help="distance along a feature checked for a bend, m (default 5 dx)")
    ap.add_argument("--downstream", type=float, help="depth of the outer band beyond a face, m (default 15 dx)")
    ap.add_argument("--plot", help="write a map (PNG)")
    a = ap.parse_args(argv)
    if not a.dx > 0:
        ap.error("--dx must be > 0")
    for name in LENGTHS[:-1]:
        v = getattr(a, name)
        if v is not None and not v > 0:
            ap.error(f"--{name.replace('_', '-')} must be > 0")
    if a.ratio < 2:
        ap.error("--ratio must be >= 2")
    if not (0 < a.angle < 90 and 0 < a.turn < 180):
        ap.error("--angle must be in (0, 90) and --turn in (0, 180)")
    if a.direction == ["all"]:
        a.direction = list(ALL_DIRECTIONS)
    else:
        try:
            a.direction = [float(d) for d in a.direction]
        except ValueError:
            ap.error("--direction takes numbers (degrees) or the single word 'all'")
        if not all(0 <= d < 360 for d in a.direction):
            ap.error("--direction must be in [0, 360)")
    if not (a.domain[0] < a.domain[2] and a.domain[1] < a.domain[3]):
        ap.error("--domain needs X0 < X1 and Y0 < Y1")
    if (a.nest is None) == (a.recommend is None):
        ap.error("give --nest (check) or --recommend WIDTH HEIGHT (search), not both")
    if a.margin is not None and a.margin < 0:
        ap.error("--margin must be >= 0")
    if a.search_step is not None and not a.search_step > 0:
        ap.error("--search-step must be > 0")
    if not a.top > 0:
        ap.error("--top must be > 0")
    if a.nest is not None and (a.search_step is not None or a.top != 5):
        ap.error("--search-step and --top go with --recommend, not --nest")
    X0, Y0, X1, Y1 = a.domain
    if not (on_cells(X1, X0, a.dx) and on_cells(Y1, Y0, a.dx)):
        ap.error(f"--domain {X0:g} {Y0:g} {X1:g} {Y1:g} must be a whole number of --dx {a.dx:g} cells in x and y")
    if a.recommend is not None:
        if not min(a.recommend) > 0:
            ap.error("--recommend sizes must be > 0")
        if not all(on_cells(v, 0.0, a.dx) for v in a.recommend):
            ap.error(f"--recommend {a.recommend[0]:g} {a.recommend[1]:g} must be whole --dx {a.dx:g} cells")
    for n in a.nest or []:
        if not (X0 < n[0] < n[2] < X1 and Y0 < n[1] < n[3] < Y1):
            ap.error(f"--nest {n[0]:g} {n[1]:g} {n[2]:g} {n[3]:g} must lie inside --domain with X0 < X1 and Y0 < Y1")
        if not all(on_cells(v, o, a.dx) for v, o in zip(n, (X0, Y0, X0, Y0))):
            ap.error(f"--nest {n[0]:g} {n[1]:g} {n[2]:g} {n[3]:g}: the corners must lie on outer nodes "
                     f"(--domain lo plus whole --dx {a.dx:g} cells), as nest.<name>.lo / hi must")
    return a


def main(argv):
    a = parse(argv)
    return check(a) if a.nest else recommend(a)[0]


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
