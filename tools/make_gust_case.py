#!/usr/bin/env python3
"""Write a GUST case folder from a YAML case description.

    tools/make_gust_case.py case.yaml                    write the case
    tools/make_gust_case.py case.yaml --dry-run          check it and print the deck, write nothing
    tools/make_gust_case.py case.yaml --direction 240 --speed 8 --name from_240
    tools/make_gust_case.py case.yaml --set mesh.layers=48 --set solver.cfl_max=8
    tools/make_gust_case.py --list-keys                  every key, its default and its range
    tools/make_gust_case.py --example                    a starter case file

The case file holds sections (domain, terrain, mesh, wind, atmosphere,
surface, model, solver, output, run, nests, gust_inputs); docs/case_setup.rst
describes them. The options below change single keys without editing the
file; ``--set section.key=value`` reaches every key. A file named on the
command line is relative to the current folder, one named in the case file to
the case file.

The tool writes the deck(s), terrain.xyz and roughness.xyz when asked for, a
run script and a summary, then runs the deck checker (check_gust_inputs.py)
on every deck it wrote.

Exit status: 0 written; 1 written, but the deck checker found an error (or a
warning with --strict); 2 the case description was refused, nothing written.
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gust_case import CaseError, apply_override, build, key_table_rst, key_table_text, load_yaml, validate  # noqa: E402

EXAMPLE = """\
# A starter case for tools/make_gust_case.py: wind from the west over flat ground.
# tools/make_gust_case.py --list-keys lists every key.
case:
  name: flat_west
domain:
  half_width: 2000.      # m from the centre to each face
  top: 1000.             # m
mesh:
  cell_size: 50.         # m
  layers: 32
  first_layer: 4.        # m; GUST stretches the layers to reach the top
wind:
  speed: 8.              # m/s at the height below
  height: 10.            # m above the ground
  direction: 270.        # degrees the wind comes from, clockwise from north
surface:
  roughness: 0.03        # m
# terrain:
#   file: dem.tif        # GeoTIFF, ESRI ASCII grid or x y z points
# domain:
#   center_lat_lon: [57.19, -7.37]
"""

# option -> (case key, help, type); the value is given to the reader as text
SHORTCUTS = [
    ("--name", "case.name", "name of the case folder", str),
    ("--speed", "wind.speed", "wind speed at wind.height (m/s)", float),
    ("--height", "wind.height", "reference height of the wind speed (m)", float),
    ("--direction", "wind.direction", "direction the wind comes from (degrees clockwise from north)", float),
    ("--drive", "wind.drive", "inflow_only, couette, geostrophic_column or geostrophic_inflow", str),
    ("--cell-size", "mesh.cell_size", "horizontal cell size (m)", float),
    ("--layers", "mesh.layers", "number of layers", int),
    ("--top", "domain.top", "height of the domain (m)", float),
    ("--roughness", "surface.roughness", "uniform roughness length (m)", float),
    ("--max-iterations", "solver.max_iterations", "steady.max_iterations", int),
    ("--ranks", "run.ranks", "MPI ranks of the run script", int),
    ("--executable", "run.executable", "the GUST executable of the run script", str),
]


def parse_args(argv):
    parser = argparse.ArgumentParser(prog="make_gust_case.py", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_file", nargs="?", help="the YAML case description")
    keys = parser.add_argument_group("single keys (each replaces the case file's value)")
    for option, key, text, kind in SHORTCUTS:
        keys.add_argument(option, type=kind, metavar="VALUE", help="%s: %s" % (key, text))
    keys.add_argument("--directions", type=float, nargs="+", metavar="DEG",
                      help="wind.directions: a sweep, one run folder per direction (replaces wind.direction)")
    keys.add_argument("--terrain", metavar="FILE", help="terrain.file, relative to the current folder")
    keys.add_argument("--output-dir", metavar="DIR", help="case.output_dir, relative to the current folder")
    keys.add_argument("--set", action="append", default=[], metavar="section.key=value", dest="assignments",
                      help="any key of the case description, the value read as YAML; repeatable. "
                           "gust_inputs.<namespace.key>=value sets one GUST input")
    actions = parser.add_argument_group("what to do")
    actions.add_argument("--dry-run", action="store_true", help="check the case and print the summary and the decks; write nothing")
    actions.add_argument("--overwrite", action="store_true", help="write into a case folder that exists")
    actions.add_argument("--no-check", action="store_true", help="do not run the deck checker on the decks written")
    actions.add_argument("--strict", action="store_true", help="exit 1 on a warning of the deck checker, not only on an error")
    actions.add_argument("--quiet", action="store_true", help="print only problems")
    actions.add_argument("--list-keys", action="store_true", help="print every key of the case description and exit")
    actions.add_argument("--example", action="store_true", help="print a starter case file and exit")
    actions.add_argument("--write-key-table", metavar="FILE", help="write the key table as reStructuredText "
                         "(docs/case_setup_keys.inc) and exit")
    return parser, parser.parse_args(argv)


def overrides_from_args(args):
    """The section.key=value assignments of the options, in a fixed order; --set comes last."""
    out = []
    for option, key, _, _ in SHORTCUTS:
        value = getattr(args, option[2:].replace("-", "_"))
        if value is not None:
            out.append((key, value))
    if args.directions is not None:
        out.append(("wind.directions", list(args.directions)))
        out.append(("wind.direction", None))
    if args.terrain is not None:
        out.append(("terrain.file", os.path.abspath(args.terrain)))
    if args.output_dir is not None:
        out.append(("case.output_dir", os.path.abspath(args.output_dir)))
    return out


def main(argv=None):
    parser, args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.direction is not None and args.directions is not None:
        parser.error("--direction and --directions: give one (--direction is one run, --directions a sweep)")
    if args.list_keys:
        sys.stdout.write(key_table_text())
        return 0
    if args.example:
        sys.stdout.write(EXAMPLE)
        return 0
    if args.write_key_table:
        with open(args.write_key_table, "w", encoding="utf-8") as out:
            out.write(key_table_rst())
        return 0
    if args.case_file is None:
        parser.print_usage(sys.stderr)
        sys.stderr.write("make_gust_case.py: a case file is required (--example prints a starter)\n")
        return 2
    try:
        raw = load_yaml(args.case_file)
        for key, value in overrides_from_args(args):
            node = raw
            parts = key.split(".")
            for part in parts[:-1]:
                if not isinstance(node.get(part), dict):
                    node[part] = {}
                node = node[part]
            node[parts[-1]] = value
        if args.direction is not None and args.directions is None and isinstance(raw.get("wind"), dict):
            raw["wind"].pop("directions", None)     # --direction replaces a sweep of the file
        for assignment in args.assignments:
            apply_override(raw, assignment)
        case, given = validate(raw)
        result = build(case, given, os.path.dirname(os.path.abspath(args.case_file)), dry_run=args.dry_run,
                       overwrite=args.overwrite, check=not args.no_check)
    except CaseError as err:
        for problem in err.problems:
            sys.stderr.write("make_gust_case.py: %s\n" % problem)
        sys.stderr.write("make_gust_case.py: the case was refused (%d problem%s); nothing was written\n"
                         % (len(err.problems), "" if len(err.problems) == 1 else "s"))
        return 2
    if args.dry_run:
        sys.stdout.write(result.summary)
        for rel, deck, _, _ in result.decks:
            sys.stdout.write("\n==== %s (not written)\n%s" % (rel, deck.text()))
        return 0
    if not args.quiet:
        sys.stdout.write(result.summary)
    elif result.findings:
        for rel, text in result.findings:
            sys.stderr.write("%s\n%s\n" % (rel, text.rstrip()))
    if result.errors or (args.strict and result.warnings):
        sys.stderr.write("make_gust_case.py: the deck checker found %d error%s and %d warning%s; see %s\n"
                         % (result.errors, "" if result.errors == 1 else "s", result.warnings,
                            "" if result.warnings == 1 else "s", os.path.join(result.case_dir, "case_summary.txt")))
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
