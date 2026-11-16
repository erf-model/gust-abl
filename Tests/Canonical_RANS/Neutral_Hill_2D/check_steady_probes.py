#!/usr/bin/env python3
"""RT-6.4b.1: steady.probes_file on the steady 2-D hill, with a restart.

    check_steady_probes.py -- LAUNCHER... GUST DECK [FLAGS ...]

DECK is Tests/Canonical_RANS/Neutral_Hill_2D/inputs_steady_hill2d_restart. The
script runs it from the log-law first guess with the four probes of
../Unsteady_Hill_2D/probes_hill2d.txt: ten iterations in straight/, and five
iterations then a restart from their checkpoint to ten in restart/ (each
emptied first). It checks:
  - the probe file is <io.plot_file>_steady_probes.csv (no <io.plot_file>_probes.csv, the unsteady
    run's name), its column line starts with "iteration" and it has the rows 0, 1, ..., 10;
  - row 0 is the first guess, before any iteration: v = 0 and theta = 300 K exactly at every probe
    (the deck's first guess has no v and a uniform theta below the inversion); row 1 has neither;
  - row 10 is the final plotfile's u, v, w, theta and tke at each probe's cell (the header's i j k),
    and differs from row 0 (the flow moved);
  - every row between differs from the one before it (a row per iteration, not the same state twice);
  - the restart: the log keeps the rows up to iteration 5, the first six rows are byte for byte the
    five-iteration run's, the file has the rows 0 ... 10 once each, and rows 6 ... 10 are the
    straight run's within 1e-9 of the largest value of each column.
The exit code is the verdict.
"""
import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

PROBES = os.path.join(HERE, "..", "Unsteady_Hill_2D", "probes_hill2d.txt")
FIELDS = ["u", "v", "w", "theta", "tke"]


def run(cmd, where, extra, clean=True):
    if clean:
        shutil.rmtree(where, ignore_errors=True)
        os.makedirs(where)
    full = cmd + extra
    print("check_steady_probes: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_steady_probes: the command failed with code %d in %s/" % (out.returncode, where))
    return out.stdout


def wrote(stdout, what):
    names = [ln.split()[-1] for ln in stdout.splitlines()
             if ln.startswith("gust: wrote ") and what in ln and "_level" not in ln]
    if not names:
        raise SystemExit("check_steady_probes: the run wrote no %s" % what)
    return names[-1]


def read_probes(path):
    """The probe file: per probe (i, j, k, h_cell) from the header, the column line, and the rows as text."""
    cells, cols, rows = [], "", []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("# probe "):
                arrow = ln.split("->")[1].split()
                cells.append((int(arrow[0]), int(arrow[1]), int(arrow[2]), float(arrow[3])))
            elif ln.startswith("#"):
                continue
            elif not cols:
                cols = ln.strip()
            elif ln.strip():
                rows.append(ln)
    return cells, cols, rows


def numbers(row):
    return [float(x) for x in row.split(",")]


def main(argv):
    if len(argv) < 4 or argv[1] != "--":
        print(__doc__)
        return 2
    cmd = argv[2:]
    bad = []

    def check(good, what):
        print("  %s: %s" % ("ok" if good else "FAIL", what))
        if not good:
            bad.append(what)

    fresh = ["init.first_guess=profile", "steady.probes_file=" + PROBES]
    out = run(cmd, "straight", fresh + ["steady.max_iterations=10"])
    files = sorted(glob.glob(os.path.join("straight", "*probes.csv")))
    check(len(files) == 1 and files[0].endswith("_steady_probes.csv"),
          "one probe file, <io.plot_file>_steady_probes.csv (found %s)" % ", ".join(files))
    if len(files) != 1:
        return 1
    check("gust: steady: 4 probes to " + os.path.basename(files[0]) + ", a row per iteration" in out,
          "the log names the probe file and its four probes")
    cells, cols, rows = read_probes(files[0])
    check(len(cells) == 4 and cols.startswith("iteration,p1_u,p1_v,p1_w,p1_theta,p1_tke,p2_u")
          and len(cols.split(",")) == 21, "the column line: iteration, then five fields for each of the four probes")
    its = [numbers(r)[0] for r in rows]
    check(its == [float(n) for n in range(11)], "the rows are the iterations 0, 1, ..., 10 (found %s)" % its)
    if its != [float(n) for n in range(11)] or len(cells) != 4:
        return 1
    first, last = numbers(rows[0]), numbers(rows[-1])

    second = numbers(rows[1])
    guess = all(first[1 + 5 * n + 1] == 0.0 and first[1 + 5 * n + 3] == 300.0 for n in range(4))
    stepped = all(second[1 + 5 * n + 1] != 0.0 and second[1 + 5 * n + 3] != 300.0 for n in range(4))
    check(guess and stepped, "row 0 is the first guess (v = 0 and theta = 300 K exactly at every probe) and row 1 is not")

    _, a = ep.read_fields(os.path.join("straight", wrote(out, "_plt")), FIELDS)
    worst = 0.0
    for n, (i, j, k, _) in enumerate(cells):
        for m, f in enumerate(FIELDS):
            ref = a[f][i][j][k]
            worst = max(worst, abs(last[1 + 5 * n + m] - ref) / max(abs(ref), 1e-30))
    check(worst <= 1e-12, "row 10 is the final plotfile at the probes' cells (largest relative difference %.2e, "
          "<= 1e-12)" % worst)
    moved = max(abs(x - y) for x, y in zip(first[1:], last[1:]))
    check(moved > 1e-3, "row 10 differs from row 0 (largest change %.3g, > 1e-3)" % moved)
    same = [n for n in range(1, 11) if rows[n].split(",", 1)[1] == rows[n - 1].split(",", 1)[1]]
    check(not same, "every row differs from the one before it (the same at %s)" % same)

    out5 = run(cmd, "restart", fresh + ["steady.max_iterations=5"])
    five = os.path.join("restart", os.path.basename(files[0]))
    with open(five) as fh:
        text5 = fh.read()
    _, _, rows5 = read_probes(five)
    check([numbers(r)[0] for r in rows5] == [float(n) for n in range(6)], "five iterations: the rows 0 ... 5")
    chk = wrote(out5, "chk")
    out10 = run(cmd, "restart", ["steady.probes_file=" + PROBES, "steady.max_iterations=10", "io.restart=" + chk], clean=False)
    check("keeps the probe rows of " + os.path.basename(five) + " up to iteration = 5 (the checkpoint is at 5)" in out10,
          "the restart's log keeps the rows up to iteration 5")
    with open(five) as fh:
        text10 = fh.read()
    check(text10.startswith(text5), "the restart keeps the five-iteration run's file byte for byte and appends to it")
    _, _, rows10 = read_probes(five)
    its = [numbers(r)[0] for r in rows10]
    check(its == [float(n) for n in range(11)], "the restart's rows are the iterations 0 ... 10, once each (found %s)" % its)
    if its == [float(n) for n in range(11)]:
        big = [max(abs(numbers(r)[c]) for r in rows) for c in range(21)]
        worst = max(abs(x - y) / max(b, 1e-30) for ra, rb in zip(rows[6:], rows10[6:])
                    for x, y, b in zip(numbers(ra)[1:], numbers(rb)[1:], big[1:]))
        check(worst <= 1e-9, "rows 6 ... 10 after the restart are the straight run's (largest relative difference "
              "%.2e, <= 1e-9)" % worst)

    if bad:
        print("check_steady_probes: %d check(s) failed" % len(bad))
        return 1
    print("check_steady_probes: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
