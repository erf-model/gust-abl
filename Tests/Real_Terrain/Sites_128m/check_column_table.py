#!/usr/bin/env python3
"""Compare a column table that GUST just wrote with a committed one (Flatirons_InflowFile, Perdigao_TargetWind).

    python3 check_column_table.py NEW --committed <case>/<case>_inflow.txt [--tol 1e-4] [--theta-tol 1e-3]

NEW is the io.column_guess_file of a run of the case's column deck (inputs_flatirons_column,
inputs_perdigao_column). The two tables must have the same heights
(to 1e-9 m), u, v and k must agree to TOL relative to the largest |value| of each column (u and v together, by
the speed), and theta to THETA_TOL in K. The friction velocity in the header must agree to TOL too. A column that is
zero everywhere in the committed table is compared absolutely. Exit status 0 when they agree.
"""
import argparse
import os
import re
import sys

def read(path):
    rows, ustar = [], None
    with open(path) as f:
        for line in f:
            s = line.strip()
            if s.startswith("#"):
                m = re.search(r"u\* ([0-9.eE+-]+) m/s", s)
                if m:
                    ustar = float(m.group(1))
                continue
            if s:
                rows.append([float(t) for t in s.split()])
    if not rows or any(len(r) != 5 for r in rows):
        raise SystemExit("%s: expected rows of z u v theta k" % path)
    if ustar is None:
        raise SystemExit("%s: no 'u* <value> m/s' in the header" % path)
    return rows, ustar


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("new")
    p.add_argument("--committed", required=True, help="the committed table")
    p.add_argument("--tol", type=float, default=1e-4)
    p.add_argument("--theta-tol", type=float, default=1e-3, help="K")
    a = p.parse_args(argv)
    if not a.tol > 0:
        p.error("--tol must be > 0")
    if not a.theta_tol > 0:
        p.error("--theta-tol must be > 0")
    new, us_new = read(a.new)
    old, us_old = read(a.committed)
    if len(new) != len(old):
        print("FAIL: %d rows against %d committed" % (len(new), len(old)))
        return 1
    bad = []
    dz = max(abs(n[0] - o[0]) for n, o in zip(new, old))
    if dz > 1e-9:
        bad.append("heights differ by %.3g m" % dz)
    speed = max((o[1] ** 2 + o[2] ** 2) ** 0.5 for o in old)
    scales = {"u": speed, "v": speed, "theta": 1.0, "k": max(abs(o[4]) for o in old)}
    tols = {"u": a.tol, "v": a.tol, "theta": a.theta_tol, "k": a.tol}
    for j, name in enumerate(("u", "v", "theta", "k"), start=1):
        err = max(abs(n[j] - o[j]) for n, o in zip(new, old)) / (scales[name] or 1.0)
        print("%-5s max difference %.3e%s" % (name, err, " K" if name == "theta" else " of %.6g" % scales[name]))
        if not err <= tols[name]:
            bad.append("%s differs by %.3e (tolerance %.1e)" % (name, err, tols[name]))
    eu = abs(us_new - us_old) / (abs(us_old) or 1.0)
    print("u*    %.6g m/s against %.6g m/s, relative difference %.3e" % (us_new, us_old, eu))
    if not eu <= a.tol:
        bad.append("u* differs by %.3e (tolerance %.1e)" % (eu, a.tol))
    if bad:
        print("FAIL: " + "; ".join(bad))
        return 1
    print("PASS: the table matches %s to %.1e" % (os.path.basename(a.committed), a.tol))
    return 0


if __name__ == "__main__":
    sys.exit(main())
