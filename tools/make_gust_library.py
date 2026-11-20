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

Exit status: 0 done; 2 refused (the message names the option or the run).
"""
import argparse
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gust_case.wind import meteorological_from_towards  # noqa: E402
from gust_library import Library, LibraryError, build  # noqa: E402
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
    a = ap.parse_args(argv)
    try:
        if a.command == "build":
            index = build(a.case_dir, a.heights, a.reference_point, a.directions, quiet=a.quiet)
            if not a.quiet:
                print("library: wrote %s (%d sectors, heights %s m)"
                      % (os.path.join(a.case_dir, FOLDER), len(index["sectors"]),
                         ", ".join("%g" % h for h in index["heights"])))
            return 0
        lib = Library(a.case_dir)
        w = lib.wind(a.direction, a.speed, a.height)
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
        return 0
    except (LibraryError, SampleError, OSError) as err:
        print("make_gust_library.py: %s" % err, file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
