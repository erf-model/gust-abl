#!/usr/bin/env python3
"""When a steady run settles at its probes, in m/s.

    probe_settle.py [--tol A,B,...] [--tail N] RUN[:REFERENCE] ...

RUN is a steady.probes_file output (<io.plot_file>_steady_probes.csv) or a directory holding one. For each run the
reference state is the mean of the last N rows (default 100) of REFERENCE (default: the run itself). The error of a
row is the length of the velocity difference (u, v, w) to the reference at each probe. Printed for each run:
  - the first iteration after which the largest error over the probes stays within each tolerance (default 0.5 and
    0.25 m/s), or "never" when a row of the reference's last N is outside it;
  - the same for the 99th percentile of the probes;
  - the largest error of the last row, and the largest distance of one of the last 200 rows from the last row
    (a run that still drifts settles against its own end state: this number shows it; "-" below 201 rows);
  - the number of probes outside each tolerance at some row of the tail.
A run that repeats a cycle is "never" within a tolerance smaller than half its swing. A run shorter than 2 N rows
is measured against a mean that holds its start, and the tool says so. Needs numpy.
"""
import glob
import os
import sys

import numpy as np


def load(path):
    if os.path.isdir(path):
        found = sorted(glob.glob(os.path.join(path, "*_steady_probes.csv")))
        if len(found) != 1:
            raise SystemExit("probe_settle: %s holds %d files *_steady_probes.csv (one expected)" % (path, len(found)))
        path = found[0]
    with open(path) as fh:
        skip = 0
        for ln in fh:
            skip += 1
            if not ln.startswith("#"):
                break                                   # the column line
    a = np.loadtxt(path, delimiter=",", skiprows=skip, dtype=np.float32, ndmin=2)
    if a.shape[1] < 6 or (a.shape[1] - 1) % 5:
        raise SystemExit("probe_settle: %s is not a probe file (iteration, then u, v, w, theta, tke per probe)" % path)
    return a[:, 0], a[:, 1:].reshape(len(a), -1, 5)[:, :, :3]


def settle(its, err, tol, tail):
    """The first iteration after which err stays within tol; None when a row of the last `tail` is outside."""
    bad = np.nonzero(~(err <= tol))[0]                  # a NaN row is outside
    if bad.size == 0:
        return int(its[0])
    return None if bad[-1] >= len(err) - min(tail, len(err)) else int(its[bad[-1] + 1])


def analyse(run, ref, tols, tail):
    its, v = load(run)
    rv = load(ref)[1] if ref else v
    if rv.shape[1] != v.shape[1]:
        raise SystemExit("probe_settle: %s and %s have different numbers of probes" % (run, ref))
    state = rv[-tail:].mean(0)
    e = np.sqrt(((v - state) ** 2).sum(-1))
    emax, e99 = e.max(1), np.percentile(e, 99, axis=1)
    worst = e[-tail:].max(0)
    drift = float(np.sqrt(((v[-201:] - v[-1]) ** 2).sum(-1)).max()) if len(v) >= 201 else None
    if len(v) < 2 * tail:
        print("probe_settle: %s has %d rows, fewer than twice the tail of %d: the reference holds the run's start"
              % (run, len(v), tail))
    return {"iterations": int(its[-1]), "probes": v.shape[1],
            "max": [settle(its, emax, t, tail) for t in tols], "p99": [settle(its, e99, t, tail) for t in tols],
            "end": float(emax[-1]), "drift": drift, "outside": [int((~(worst <= t)).sum()) for t in tols]}


def main(argv):
    tols, tail, runs = [0.5, 0.25], 100, []
    it = iter(argv[1:])
    try:
        for a in it:
            if a == "--tol":
                tols = [float(x) for x in next(it, "").split(",")]
            elif a == "--tail":
                tail = int(next(it, ""))
            elif a in ("-h", "--help"):
                print(__doc__)
                return 0
            else:
                runs.append(a)
    except ValueError:
        runs = []                                       # --tol or --tail without a usable value
    if not runs or tail < 1 or any(t <= 0.0 for t in tols):
        print(__doc__)
        return 2
    show = lambda x: "never" if x is None else "%d" % x  # noqa: E731
    head = ["all<%g" % t for t in tols] + ["99%%<%g" % t for t in tols]
    print("%-28s %6s | %s | %8s %8s | probes outside %s in the tail"
          % ("run", "its", " ".join("%9s" % h for h in head), "end", "last200", " / ".join("%g" % t for t in tols)))
    for arg in runs:
        run, _, ref = arg.partition(":")
        r = analyse(run, ref, tols, tail)
        name = os.path.basename(os.path.normpath(run)) + ("*" if ref else "")
        print("%-28s %6d | %s | %8.3f %8s | %s of %d"
              % (name, r["iterations"], " ".join("%9s" % show(x) for x in r["max"] + r["p99"]), r["end"],
                 "-" if r["drift"] is None else "%.3f" % r["drift"],
                 " / ".join("%d" % n for n in r["outside"]), r["probes"]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
