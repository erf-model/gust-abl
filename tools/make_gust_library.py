#!/usr/bin/env python3
"""Build and query a wind library from a sweep written by make_gust_case.py (gust_library/).

    tools/make_gust_library.py build CASE_DIR [--heights 10 50] [--reference-point X Y] [--directions D ...]
    tools/make_gust_library.py query CASE_DIR --direction D --speed S --height H [--point X Y] [--out FILE.npz]

build reads each dir_<DDD>/ run of the case folder (after run.sh) and writes
CASE_DIR/wind_library/: index.json, grid.npz and one sector_<DDD>.npz per
direction. A build replaces the case's library only when every run is read.
--heights: heights above the ground (m), each at or above the first cell
centre and at or below the last in every column (default 10).
--reference-point: where the reference speed is taken, at the first height
(m, case frame: x and y from domain.center; default the domain centre).
--directions: a subset of the sweep (default all of it).

query blends the two sectors either side of --direction (degrees the wind
comes from, the sweep's inflow direction) and scales by --speed, the local
wind speed at the reference point and first height (not the inflow speed).
It prints u, v, w, the speed, the direction and k at --point (default the
reference point), and with --out writes the fields at the height to a
.npz file (u, v, w, k, x, y, ground).

query and ensemble also write WindNinja-format files at the height:
--grids DIR: <prefix>_vel.asc and <prefix>_ang.asc (speed, and the direction
the wind comes from), ESRI ASCII grids with a .prj (written when the case
gave domain.center_lat_lon; --prj FILE copies a projection file instead);
--points FILE: a CSV of the values at points ('x y [name]' per line, m, case
frame) in DIR (or the current folder). --speed-units mps|mph|kph for the
speed grids (default mps; the CSV stays in m/s). --prefix: the start of the
file names. --prj, --speed-units and --prefix are refused without the
options they apply to. Columns that are not square are resampled onto square
cells (an ESRI ASCII grid has one cell size); the CSV is not resampled.

    tools/make_gust_library.py ensemble CASE_DIR --height H (--members FILE | --directions D ... --speeds S ...)
        [--direction-weights W ...] [--speed-weights W ...] [--percentiles 10 50 90] [--thresholds T ...]
        [--grids DIR] [--points FILE] [--out FILE.npz]

ensemble takes the wind of every member (a direction, a speed at the
reference point and first height, a weight) from the library, with no new
runs, and gives per column the mean speed, its standard deviation, the vector
mean and its direction, weighted percentiles of the speed, the probability
of a speed above each threshold and the probability of a reversal (a local
wind against the member's inflow direction). Members come from --members
(lines 'direction speed [weight]') or from every --directions with every
--speeds, weighted by the product of --direction-weights and --speed-weights
(uniform when left out). The weights are divided by their sum. --grids
writes <prefix>_vel.asc (mean speed) and <prefix>_ang.asc (direction of the
vector mean) and one grid per other statistic (<prefix>_<name>.asc: std_speed,
p<P>, pexceed_<T>, preversal; the vector mean's components are in --out and
the CSV only).

Exit status: 0 done; 2 refused (the message names the option or the run).
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gust_case.wind import meteorological_from_towards  # noqa: E402
from gust_library import Library, LibraryError, build, ensemble, grids  # noqa: E402
from gust_library.index import FOLDER  # noqa: E402
from gust_library.sample import SampleError, at_point  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="command", required=True)
    b = sub.add_parser("build", help="write CASE_DIR/wind_library/ from the sweep's runs")
    b.add_argument("case_dir")
    b.add_argument("--heights", type=float, nargs="+", default=[10.0], metavar="H",
                   help="heights above the ground, m (default 10)")
    b.add_argument("--reference-point", type=float, nargs=2, metavar=("X", "Y"),
                   help="where the reference speed is taken, m, case frame: from domain.center "
                        "(default the domain centre)")
    b.add_argument("--directions", type=float, nargs="+", metavar="D", help="a subset of the sweep's directions")
    b.add_argument("--quiet", action="store_true", help="print nothing when the build succeeds")
    q = sub.add_parser("query", help="the wind for a direction and speed at a stored height")
    q.add_argument("case_dir", help="the case folder (or its wind_library/ folder)")
    q.add_argument("--direction", type=float, required=True, help="degrees the wind comes from, clockwise from north")
    q.add_argument("--speed", type=float, required=True,
                   help="m/s at the reference point and first height (the local speed there, not the inflow speed)")
    q.add_argument("--height", type=float, required=True, help="m above the ground, one of the library's heights")
    q.add_argument("--point", type=float, nargs=2, metavar=("X", "Y"),
                   help="m, case frame: from domain.center (default the reference point)")
    q.add_argument("--out", help="write the fields at the height to this .npz file")
    _output_options(q)
    e = sub.add_parser("ensemble", help="statistics over a distribution of directions and speeds")
    e.add_argument("case_dir", help="the case folder (or its wind_library/ folder)")
    e.add_argument("--height", type=float, required=True, help="m above the ground, one of the library's heights")
    e.add_argument("--members", help="file of members, 'direction speed [weight]' per line")
    e.add_argument("--directions", type=float, nargs="+", metavar="D", help="degrees the wind comes from")
    e.add_argument("--direction-weights", type=float, nargs="+", metavar="W", help="one per direction (default equal)")
    e.add_argument("--speeds", type=float, nargs="+", metavar="S",
                   help="m/s at the reference point and first height")
    e.add_argument("--speed-weights", type=float, nargs="+", metavar="W", help="one per speed (default equal)")
    e.add_argument("--percentiles", type=float, nargs="*", default=[10.0, 50.0, 90.0], metavar="P",
                   help="weighted percentiles of the speed, each in (0, 100] (default 10 50 90)")
    e.add_argument("--thresholds", type=float, nargs="*", default=[], metavar="T",
                   help="m/s: the probability of a speed above each (default none)")
    e.add_argument("--out", help="write the statistics at the height to this .npz file")
    _output_options(e)
    a = ap.parse_args(argv)
    try:
        if a.command == "build":
            index = build(a.case_dir, a.heights, a.reference_point, a.directions, quiet=a.quiet)
            if not a.quiet:
                print("library: wrote %s (%d sectors, heights %s m)"
                      % (os.path.join(a.case_dir, FOLDER), len(index["sectors"]),
                         ", ".join("%g" % h for h in index["heights"])))
            return 0
        _check_output_options(a)
        lib = Library(a.case_dir)
        if a.command == "ensemble":
            return _ensemble(lib, a)
        w = lib.wind(a.direction, a.speed, a.height)
        prepared = _prepare_outputs(lib, a, lambda pts: ensemble.wind_at_points(lib, a.direction, a.speed,
                                                                                 a.height, pts))
        x, y = a.point if a.point else lib.index["reference_point"]
        vals = {f: float(at_point(lib.x, lib.y, w[f], x, y, "--point")) for f in ("u", "v", "w", "k")}
        sp = math.hypot(vals["u"], vals["v"])
        frm = meteorological_from_towards(math.degrees(math.atan2(vals["v"], vals["u"])))
        print("blend: " + ", ".join("%g degrees x %.6g" % dw for dw in w["weights"]))
        print("at (%g, %g) m, %g m above the ground: u %.6g  v %.6g  w %.6g m/s; speed %.6g m/s from %.4f degrees; "
              "k %.6g m^2/s^2" % (x, y, a.height, vals["u"], vals["v"], vals["w"], sp, frm, vals["k"]))
        if a.out:
            out = a.out if a.out.endswith(".npz") else a.out + ".npz"   # numpy adds the suffix
            np.savez(out, **{f: w[f] for f in ("u", "v", "w", "k", "x", "y", "ground")})
            print("wrote %s" % out)
        if a.grids or a.points:
            prefix = a.prefix or "%s_%sdeg_%smps_%sm_agl" % (_case_name(lib), _g(a.direction % 360.0), _g(a.speed),
                                                            _g(a.height))
            fields = dict(speed=np.hypot(w["u"], w["v"]), u=w["u"], v=w["v"])
            _write_outputs(lib, a, prefix, fields, ("speed", "u", "v"), (), (), prepared)
        return 0
    except (LibraryError, SampleError, OSError) as err:
        print("make_gust_library.py: %s" % err, file=sys.stderr)
        return 2


def _output_options(p):
    p.add_argument("--grids", metavar="DIR", help="write WindNinja-format ESRI ASCII grids (and .prj) to DIR")
    p.add_argument("--points", metavar="FILE", help="write a CSV of the values at the points of FILE "
                                                    "('x y [name]' per line, m, case frame)")
    p.add_argument("--prefix", help="the start of the output file names (default from the case and the options)")
    p.add_argument("--speed-units", choices=list(grids.SPEED_UNITS), default="mps",
                   help="units of the speed grids: mps, mph or kph (default mps; the CSV is in m/s)")
    p.add_argument("--prj", metavar="FILE", help="copy this projection file beside each grid (in place of the "
                                                 "WGS 84 / UTM one written for a centre given in latitude and "
                                                 "longitude)")


def _g(x):
    return ("%g" % x).replace(".", "p")


def _case_name(lib):
    return lib.index.get("case") or "library"


def _check_output_options(a):
    """Refuses output options that would be ignored, before anything is read or written."""
    if not a.grids:
        if a.prj:
            raise LibraryError("--prj: only used with --grids")
        if a.speed_units != "mps":
            raise LibraryError("--speed-units: only used with --grids (the CSV is in m/s)")
    if a.prefix is not None and not (a.grids or a.points):
        raise LibraryError("--prefix: only used with --grids or --points")
    if a.prj and not os.path.isfile(a.prj):
        raise LibraryError("--prj: %s not found" % a.prj)
    if a.grids and os.path.exists(a.grids) and not os.path.isdir(a.grids):
        raise LibraryError("--grids: %s is a file, not a folder" % a.grids)


def _prepare_outputs(lib, a, point_values):
    """Everything --grids and --points can refuse, checked before any file is written: the grid's cells and
    georeference, the points file and the values at its points (point_values(points) at [(x, y)])."""
    prepared = {}
    if a.grids:
        grids.grid_of(lib.x, lib.y)
        grids.georeference(lib)
    if a.points:
        points = grids.read_points(a.points)
        prepared["points"] = points
        prepared["values"] = point_values([(x, y) for _, x, y in points])
    return prepared


def _write_outputs(lib, a, prefix, fields, wind, extra, speed_names, prepared):
    """--grids and --points from the fields and what _prepare_outputs read."""
    written = []
    if a.grids:
        files, geo, note = grids.export_grids(lib, a.grids, prefix, fields, wind, extra, speed_names,
                                              a.speed_units, a.prj)
        written += files
        print("grids: corner and .prj from %s%s%s" % (geo.source, "; .prj copied from %s" % a.prj if a.prj else "",
                                                    "; " + note if note else ""))
    if a.points:
        out_dir = a.grids or "."
        written += grids.export_points(lib, os.path.join(out_dir, prefix + "_points.csv"), prepared["points"],
                                       prepared["values"], a.height)
    for f in written:
        print("wrote %s" % f)


def _ensemble(lib, a):
    if a.members and (a.directions or a.speeds or a.direction_weights or a.speed_weights):
        raise LibraryError("--members: give the members either in a file or as --directions and --speeds, not both")
    if a.members:
        members = ensemble.from_file(a.members)
        source = "--members %s" % a.members
    elif a.directions or a.speeds:
        members = ensemble.product(a.directions, a.direction_weights, a.speeds, a.speed_weights)
        source = "--directions and --speeds"
    else:
        raise LibraryError("ensemble: give --members FILE, or --directions and --speeds")
    stats = ensemble.ensemble(lib, members, a.height, a.percentiles, a.thresholds, source)
    prepared = _prepare_outputs(lib, a, lambda pts: ensemble.at_points(lib, members, a.height, pts, a.percentiles,
                                                                        a.thresholds, source))
    x, y = lib.index["reference_point"]
    ref = ensemble.at_points(lib, members, a.height, [(x, y)], (), (), source)
    print("ensemble: %d members, %d directions; at the reference point (%g, %g) m, %g m above the ground: mean speed "
          "%.6g m/s, standard deviation %.6g m/s, probability of a reversal %.4g"
          % (len(members), len({m[0] for m in members}), x, y, a.height, ref["mean_speed"][0],
             ref["std_speed"][0], ref["preversal"][0]))
    if a.out:
        out = a.out if a.out.endswith(".npz") else a.out + ".npz"
        np.savez(out, **stats)
        print("wrote %s" % out)
    if a.grids or a.points:
        prefix = a.prefix or "%s_ensemble_%sm_agl" % (_case_name(lib), _g(a.height))
        extra = [n for n in stats if n not in ("mean_speed", "mean_direction", "mean_u", "mean_v", "x", "y")]
        _write_outputs(lib, a, prefix, stats, ("mean_speed", "mean_u", "mean_v"), extra,
                       ensemble.speed_names(a.percentiles), prepared)
    return 0


if __name__ == "__main__":
    sys.exit(main())
