#!/usr/bin/env python3
"""Download terrain for a GUST case: SRTM-format tiles from the AWS Terrain Tiles.

    tools/fetch_terrain.py case.yaml                       the domain of a case file, written to dem.asc beside it
    tools/fetch_terrain.py case.yaml --out terrain/dem.asc --sea-floor keep
    tools/fetch_terrain.py --box 57.15 57.23 -7.45 -7.30 --out askervein.asc

The output is an ESRI ASCII grid in longitude and latitude; give it to the
case file as ``terrain.file`` (with ``domain.center_lat_lon``). Beside it,
``<name>.sources.txt`` holds the attribution the data need. Tiles are kept in
a cache (``--cache``, GUST_DEM_CACHE or ~/.cache/gust/terrain_tiles) and
fetched again only when missing. No key or account is needed; the data are
those of https://registry.opendata.aws/terrain-tiles/ (docs/case_setup.rst).

Exit status: 0 written; 2 refused (nothing written).
"""
import argparse
import datetime
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gust_case import CaseError, load_case  # noqa: E402
from gust_case import dem  # noqa: E402


def main(argv=None):
    parser = argparse.ArgumentParser(prog="fetch_terrain.py", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_file", nargs="?", help="a YAML case description with domain.center_lat_lon")
    parser.add_argument("--box", type=float, nargs=4, metavar=("SOUTH", "NORTH", "WEST", "EAST"),
                        help="the box in degrees, in place of a case file")
    parser.add_argument("--out", help="the grid file (.asc); default dem.asc beside the case file, or in the "
                                      "current folder with --box")
    parser.add_argument("--sea-floor", choices=dem.SEA_FLOOR, default="zero",
                        help="zero: heights below 0 m (the sea floor over the sea) become 0 m (default); "
                             "keep: left as they are (land below sea level)")
    parser.add_argument("--cache", help="folder of downloaded tiles (default GUST_DEM_CACHE or "
                                        "~/.cache/gust/terrain_tiles)")
    args = parser.parse_args(sys.argv[1:] if argv is None else argv)
    if (args.case_file is None) == (args.box is None):
        parser.error("give a case file or --box, one of the two")
    try:
        if args.case_file is not None:
            case, given = load_case(args.case_file)
            box = dem.box_for_case(case, given)
            out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.case_file)), "dem.asc")
        else:
            box = tuple(args.box)
            out = args.out or "dem.asc"

        def progress(done, total, name):
            sys.stdout.write("tile %d of %d: %s\n" % (done, total, name))
            sys.stdout.flush()

        result = dem.download(box, out, args.cache or dem.default_cache_dir(), sea_floor=args.sea_floor,
                              date=datetime.date.today().isoformat(), progress=progress)
    except CaseError as err:
        for problem in err.problems:
            sys.stderr.write("fetch_terrain.py: %s\n" % problem)
        return 2
    sys.stdout.write("wrote %s: %d x %d points, %s to %s m\n" % (
        result.path, result.shape[1], result.shape[0],
        "-" if result.lowest is None else "%.6g" % result.lowest,
        "-" if result.highest is None else "%.6g" % result.highest))
    for note in result.notes:
        sys.stdout.write("note: %s\n" % note)
    sys.stdout.write("attribution: %s\n" % result.sources)
    return 0


if __name__ == "__main__":
    sys.exit(main())
