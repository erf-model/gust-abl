#!/usr/bin/env python3
"""Check GUST input files before running them.

    python3 tools/check_gust_inputs.py inputs
    python3 tools/check_gust_inputs.py inputs steady.max_iterations=500 io.plot_file=test_plt
    python3 tools/check_gust_inputs.py --ranks 16 --min-severity WARNING --strict Tests/*/*/inputs*
    python3 tools/check_gust_inputs.py --self-test

An argument that contains '=' and does not name an existing file is a
command-line definition, as on gust's own command line; it applies to every
deck. For each deck the script prints, from most to least severe:

    ERROR    GUST would stop at start-up
    WARNING  GUST starts, but the setting is very likely not what was meant
    ADVICE   a setting that could be chosen better
    INFO     derived numbers (cell sizes, layers, boxes, wind direction)

Files named in a deck are looked up from the current directory, as gust does.
The exit code is 1 when any deck has an ERROR, or with --strict when any deck
has a WARNING. See gust_input_checker/__init__.py for the design and
gust_input_checker/checks/ for what is checked.
"""
import argparse
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from gust_input_checker import Options, Severity, check_deck, registry  # noqa: E402


def split_arguments(items):
    """The decks and the command-line definitions among the positional arguments."""
    decks, overrides = [], []
    for item in items:
        if "=" in item and not os.path.exists(item):
            overrides.append(item)
            continue
        hits = sorted(glob.glob(item))
        decks.extend(hits if hits else [item])
    return decks, overrides


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("inputs", nargs="*", help="input files, and key=value definitions applied to each of them")
    ap.add_argument("--ranks", type=int, help="MPI ranks the run will use; enables the boxes-per-rank advice")
    ap.add_argument("--gpu", action="store_true", help="the run uses a GPU build (default sweep and box-size advice)")
    ap.add_argument("--min-severity", default="INFO", choices=[s.name for s in Severity],
                    help="hide findings below this level")
    ap.add_argument("--strict", action="store_true", help="exit 1 on WARNING as well as ERROR")
    ap.add_argument("--json", action="store_true", help="JSON output")
    ap.add_argument("--only", action="append", metavar="NAME",
                    help="run only the named checker(s), e.g. --only mesh --only faces")
    ap.add_argument("--list-checkers", action="store_true", help="list the checkers and exit")
    ap.add_argument("--reference", metavar="PATH",
                    help="path of Tests/inputs_reference (found automatically inside the repository)")
    ap.add_argument("--made-by-run", action="append", metavar="PATH", default=[],
                    help="a file or checkpoint directory (or a glob pattern) that an earlier step of the run writes: "
                         "a key naming it is not reported missing (repeatable)")
    ap.add_argument("--self-test", action="store_true", help="run the built-in tests and exit")
    args = ap.parse_args(argv)

    if args.self_test:
        from gust_input_checker import selftest
        return selftest.main()
    if args.list_checkers:
        for cls in registry():
            print("%-24s %s" % (cls.name, (cls.__doc__ or "").strip().splitlines()[0] if cls.__doc__ else ""))
        return 0
    if args.ranks is not None and args.ranks < 1:
        ap.error("--ranks must be at least 1")
    known = [cls.name for cls in registry()]
    for name in args.only or []:
        if name not in known:
            ap.error("--only %s: no such checker (see --list-checkers)" % name)

    decks, overrides = split_arguments(args.inputs)
    if not decks:
        ap.error("give at least one input file")

    opts = Options(ranks=args.ranks, gpu=args.gpu, reference=args.reference, made_by_run=args.made_by_run)
    minsev = Severity[args.min_severity]
    worst = Severity.INFO
    reports = []
    for path in decks:
        if os.path.isdir(path):
            continue
        rep = check_deck(path, overrides=overrides, options=opts, only=args.only)
        reports.append(rep)
        if rep.worst is not None and rep.worst > worst:
            worst = rep.worst
        if args.json:
            print(rep.as_json())
        else:
            print(rep.format(minsev))
            print()
    if len(reports) > 1 and not args.json:
        print("== %d decks: %d with errors, %d with warnings" % (
            len(reports), sum(1 for r in reports if r.count(Severity.ERROR)),
            sum(1 for r in reports if r.count(Severity.WARNING))))
    if worst >= Severity.ERROR or (args.strict and worst >= Severity.WARNING):
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
