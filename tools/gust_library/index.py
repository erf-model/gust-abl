"""Build a wind library from a sweep's case folder.

    build(case_dir, heights, reference_point=None, directions=None)   writes <case_dir>/wind_library/

A case folder written by tools/make_gust_case.py with ``wind.directions``
holds one ``dir_<DDD>/`` run per direction. For each, the library reads the
run's log (GUST's verdict, the iterations and the last relative residual r)
and its final plotfile, samples u, v, w and k at each height above the ground
(sample.py), and stores them divided by the sector's reference speed (k by
its square): the horizontal wind speed at the reference point and the first
height. Fields are in the case frame (GUST's mesh coordinates); u is towards
+x (east), v towards +y (north).

Written in <case_dir>/wind_library/ (built in wind_library.tmp/ and swapped in
whole, so a refused build leaves the previous library as it was):

    index.json             the heights, the reference point, the grid and one entry per sector
    grid.npz               the column centres x, y and the ground under them
    sector_<DDD>.npz       u, v, w, k (len(heights), nx, ny) divided by the reference speed (k by its square)

A sector whose run did not reach GUST's steady verdict is written with
``"converged": false``; the query (query.py) refuses it. A run that wrote no
verdict, no plotfile after it, an unsteady run and a two-way joint run are
refused here.
"""
import json
import math
import os
import re
import shutil
import sys

import numpy as np
import yaml

from . import plotfile, sample

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from gust_case import wind  # noqa: E402

FORMAT_VERSION = 1
FOLDER = "wind_library"
INDEX = "index.json"
VERDICT = re.compile(r"^gust: (steady|not steady|diverged) after (\d+) iterations")
GRID_TOL = 1.0e-6     # m: the sectors' column centres and ground must agree
R_COLUMN = 9          # a table row: Iter, CFL, line solver, six residuals, r, then optional notes


class LibraryError(ValueError):
    """A case folder, run or option the library refuses; the message says which."""


def deck_value(path, key, default=None):
    """The value of ``key`` in a GUST deck (a string), or default."""
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split("#")[0].split()
            if len(parts) >= 3 and parts[0] == key and parts[1] == "=":
                return " ".join(parts[2:])
    return default


def sweep(case_dir):
    """(case name, [(direction from, folder name)]) from case_resolved.yaml."""
    path = os.path.join(case_dir, "case_resolved.yaml")
    if not os.path.isfile(path):
        raise LibraryError("%s: not found; the library needs a case folder written by make_gust_case.py" % path)
    try:
        with open(path, encoding="utf-8") as fh:
            case = yaml.safe_load(fh)
        name = case["case"]["name"]
        spec = case["wind"].get("directions")
    except (yaml.YAMLError, TypeError, KeyError, AttributeError) as err:
        raise LibraryError("%s: not a case description written by make_gust_case.py (%s)" % (path, err))
    if spec is None:
        raise LibraryError("%s: no wind.directions; a library needs a sweep (one run per direction)" % path)
    return name, [(d, "dir_" + wind.direction_label(d)) for d in wind.sweep_directions(spec)]


def table_r(line):
    """r of an iteration-table row of the level-0 run (None for any other line, or a row whose r is '-')."""
    tok = line.split()
    if len(tok) <= R_COLUMN or not tok[0].isdigit():
        return None
    try:
        if tok[1] not in ("-", "inf"):   # CFL; a Newton row may print "-" or "inf"
            float(tok[1])
        return float(tok[R_COLUMN])
    except ValueError:
        return None


def read_log(folder):
    """GUST's verdict in folder/log.txt: (status, iterations, last r before it, plotfile written after it).

    The plotfile is the last <io.plot_file><step> written after the verdict (the deck's io.plot_file; one given
    on the command line is not seen); r is None when the table has no r before the verdict."""
    log = os.path.join(folder, "log.txt")
    if not os.path.isfile(log):
        raise LibraryError("%s: not found; run run.sh first" % log)
    deck = os.path.join(folder, "inputs")
    if not os.path.isfile(deck):
        raise LibraryError("%s: not found; the run folder holds the deck GUST ran" % deck)
    prefix = deck_value(deck, "io.plot_file", "plt")
    status, iterations, last_r, plt = None, None, None, None
    with open(log, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            m = VERDICT.match(line)
            if m:
                status, iterations, plt = m.group(1), int(m.group(2)), None
            elif line.startswith("gust: joint run:"):
                raise LibraryError("%s: a two-way joint run (nest.names); the library reads single-domain runs" % log)
            elif line.startswith("gust: unsteady run"):
                raise LibraryError("%s: an unsteady run; the library reads steady runs" % log)
            elif status is None and table_r(line) is not None:
                last_r = table_r(line)
            elif status is not None and line.startswith("gust: wrote "):
                name = line.split()[2]
                if re.fullmatch(re.escape(prefix) + r"\d+", name):
                    plt = os.path.join(folder, name)
    if status is None:
        raise LibraryError("%s: no verdict line (gust: steady/not steady/diverged after N iterations); "
                           "the run did not finish" % log)
    if plt is None or not os.path.isfile(os.path.join(plt, "Header")):
        raise LibraryError("%s: no plotfile %s<step> written after the verdict (io.write_plotfile = 1 writes one; "
                           "an io.plot_file given on the command line is not seen)" % (log, prefix))
    return status, iterations, last_r, plt


def sector(folder, heights, reference_point):
    """One sector's sampled fields and its entry (before the index is written)."""
    status, iterations, last_r, plt = read_log(folder)
    cells, hdr = plotfile.read_cells(plt, ["u", "v", "w", "tke", "z"])
    z_nodes = plotfile.read_nodes(plt)
    agl = sample.heights_above_ground(cells["z"], z_nodes)
    dx = hdr["dx"]
    nx, ny, nz = cells["u"].shape
    if min(nx, ny, nz) < 2:
        raise LibraryError("%s: %d x %d x %d cells; the library needs at least 2 in each direction" % (plt, nx, ny, nz))
    xc = hdr["prob_lo"][0] + (np.arange(nx) + 0.5) * dx[0]
    yc = hdr["prob_lo"][1] + (np.arange(ny) + 0.5) * dx[1]
    if reference_point is None:
        reference_point = (0.5 * (hdr["prob_lo"][0] + hdr["prob_hi"][0]), 0.5 * (hdr["prob_lo"][1] + hdr["prob_hi"][1]))
    fields = {f: sample.at_heights(cells[f], agl, heights) for f in ("u", "v", "w", "tke")}
    if not all(np.isfinite(a).all() for a in fields.values()):
        raise LibraryError("%s: a sampled field is not finite" % plt)
    u_ref = float(sample.at_point(xc, yc, fields["u"][0], *reference_point))
    v_ref = float(sample.at_point(xc, yc, fields["v"][0], *reference_point))
    s_ref = math.hypot(u_ref, v_ref)
    if not s_ref > 0.0:
        raise LibraryError("%s: no wind at the reference point (%g, %g) m, %g m above the ground"
                           % (plt, reference_point[0], reference_point[1], heights[0]))
    deck = os.path.join(folder, "inputs")
    z0, z0_in = deck_value(deck, "bc.surface_z0"), deck_value(deck, "bc.inflow_z0")
    entry = {"status": status, "converged": status == "steady", "iterations": iterations, "last_r": last_r,
             "plotfile": os.path.relpath(plt, os.path.dirname(folder)),
             "reference_speed": s_ref,
             "reference_direction": wind.meteorological_from_towards(math.degrees(math.atan2(v_ref, u_ref))),
             "surface_z0": float(z0) if z0 is not None else None,
             "inflow_z0": float(z0_in) if z0_in is not None else None}
    stored = {"u": fields["u"] / s_ref, "v": fields["v"] / s_ref, "w": fields["w"] / s_ref,
              "k": fields["tke"] / s_ref ** 2}
    grid = {"x": xc, "y": yc, "ground": sample.ground(z_nodes)}
    return entry, stored, grid, reference_point


def build(case_dir, heights=(10.0,), reference_point=None, directions=None, quiet=True):
    """Write <case_dir>/wind_library/; returns the index (a dict). Every run is read before anything is written."""
    heights = [float(h) for h in heights]
    if not heights:
        raise LibraryError("--heights: give at least one height")
    if reference_point is not None:
        reference_point = tuple(float(v) for v in reference_point)
        if len(reference_point) != 2 or not all(math.isfinite(v) for v in reference_point):
            raise LibraryError("--reference-point: give two finite coordinates x y (m, case frame)")
    name, runs = sweep(case_dir)
    if directions is not None:
        labels = {wind.direction_label(d): (d, f) for d, f in runs}
        chosen = []
        for d in directions:
            lab = wind.direction_label(d)
            if lab not in labels:
                raise LibraryError("--directions: %g is not a direction of the sweep (%s)"
                                   % (d, ", ".join("%g" % r[0] for r in runs)))
            if labels[lab] in chosen:
                raise LibraryError("--directions: %g is given twice" % d)
            chosen.append(labels[lab])
        runs = chosen
    sectors, grid0 = [], None
    for direction, folder in sorted(runs, key=lambda r: r[0] % 360.0):
        try:
            entry, stored, grid, reference_point = sector(os.path.join(case_dir, folder), heights, reference_point)
        except sample.SampleError as err:   # an option out of range: the same for every run (one mesh)
            raise LibraryError(str(err))
        except (ValueError, OSError) as err:   # LibraryError and the plotfile reader's
            raise LibraryError("%s: %s" % (folder, err))
        if grid0 is None:
            grid0 = grid
        elif any(grid[k].shape != grid0[k].shape or np.abs(grid[k] - grid0[k]).max() > GRID_TOL for k in grid):
            raise LibraryError("%s: the mesh differs from %s's; every sector needs the same mesh"
                               % (folder, sectors[0][0]["folder"]))
        label = wind.direction_label(direction)
        entry = dict({"direction": direction % 360.0, "folder": folder, "file": "sector_%s.npz" % label}, **entry)
        sectors.append((entry, stored))
        if not quiet:
            print("library: %s  %s after %d iterations, reference wind %.4g m/s from %.3f degrees"
                  % (folder, entry["status"], entry["iterations"], entry["reference_speed"], entry["reference_direction"]))
    index = {"format_version": FORMAT_VERSION, "case": name, "heights": heights,
             "reference_point": list(reference_point), "reference_height": heights[0], "grid": "grid.npz",
             "frame": "case frame (GUST mesh coordinates, m); u towards +x (east), v towards +y (north)",
             "stored": "u, v, w divided by reference_speed, k by its square",
             "sectors": [e for e, _ in sectors]}
    out_dir = os.path.join(case_dir, FOLDER)
    tmp_dir, old_dir = out_dir + ".tmp", out_dir + ".old"
    for d in (tmp_dir, old_dir):
        shutil.rmtree(d, ignore_errors=True)
    os.makedirs(tmp_dir)
    for entry, stored in sectors:
        np.savez(os.path.join(tmp_dir, entry["file"]), **stored)
    np.savez(os.path.join(tmp_dir, "grid.npz"), **grid0)
    with open(os.path.join(tmp_dir, INDEX), "w", encoding="utf-8") as fh:
        json.dump(index, fh, indent=1, allow_nan=False)
        fh.write("\n")
    if os.path.isdir(out_dir):
        os.replace(out_dir, old_dir)
    os.replace(tmp_dir, out_dir)
    shutil.rmtree(old_dir, ignore_errors=True)
    return index
