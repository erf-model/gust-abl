#!/usr/bin/env python3
"""The Terrain_Shapes_3D case on 200 m cells: the residual verdict against the probes.

    check_terrain_shapes.py -- LAUNCHER... GUST DECK [FLAGS ...]

DECK is inputs_terrain_shapes. The script writes the round gentle and the round steep butte on 200 m cells with
make_terrain.py (80 x 64 x 40 cells, 8 m first layers, 2560 probes) and runs, each in its own directory (emptied
first):
  gentle  the 32 deg butte   to the steady verdict (at most 300 iterations)
  steep   the 65 deg butte   200 iterations
It checks:
  - gentle reaches the steady verdict within 250 iterations (136 measured on 2 ranks) with no ground-step warning;
  - steep prints the ground-step warning (94 first layers, above the tested 64) and has no verdict at 200;
  - steep's probes are settled all the same: against the mean of rows 151 ... 200, the velocity at every probe
    stays within 0.5 m/s from iteration 40 at the latest and within 0.25 m/s from iteration 100 at the latest
    (26 and 59 measured), and no probe moves more than 0.2 m/s between iterations 100 and 200 (this last number
    does not use the end state as its reference);
  - gentle's probes stay within 0.25 m/s of its last row from iteration 40 at the latest (29 measured);
  - the rows are the run's: steep's first guess is more than 1 m/s from row 200 somewhere, and gentle's first
    probe (300 m from the inflow face) ends within 0.02 m/s of the inflow log law u*/kappa ln(z/z0) at its cell
    centre's height, read from the file's header (9.466 m/s at 72.65 m; 9.4662 measured);
  - a value that is not finite in a probe file fails the check.
Mutant: steady.cfl_max = 2 in the two runs gives 48 for steep's 0.5 m/s and 53 for gentle, and fails the two
limits of 40. The 0.25 m/s limit of 100 and the 0.2 m/s between iterations 100 and 200 are guards against a run
that drifts; that mutant does not fail them.
The exit code is the verdict.
"""
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_terrain as mt  # noqa: E402

DX = 200.0
RUNS = (("gentle", "round_gentle", 300), ("steep", "round_steep", 200))


def rows_of(path):
    """The probe file as (iteration, [(u, v, w) per probe]) rows."""
    out, cols = [], False
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("#"):
                continue
            if not cols:
                cols = True
                continue
            x = ln.split(",")
            x = [float(v) for v in x]
            if not all(math.isfinite(v) for v in x):
                raise SystemExit("check_terrain_shapes: %s has a value that is not finite at iteration %g"
                                 % (path, x[0]))
            out.append((int(x[0]), [(x[n], x[n + 1], x[n + 2]) for n in range(1, len(x), 5)]))
    return out


def error(row, ref):
    return max(math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2) for a, b in zip(row, ref))


def settled_from(rows, ref, tol):
    """The first iteration after which every probe stays within tol of ref; None when the last row is outside."""
    last_bad = None
    for n, (_, row) in enumerate(rows):
        if not error(row, ref) <= tol:                  # a NaN row is outside
            last_bad = n
    if last_bad is None:
        return rows[0][0]
    return None if last_bad == len(rows) - 1 else rows[last_bad + 1][0]


def main(argv):
    if len(argv) < 4 or argv[1] != "--":
        print(__doc__)
        return 2
    cmd = argv[2:]
    deck_at = [n for n, a in enumerate(cmd) if os.path.basename(a) == "inputs_terrain_shapes"]
    if not deck_at:
        raise SystemExit("check_terrain_shapes: the command does not name inputs_terrain_shapes")
    bad = []

    def check(good, what):
        print("  %s: %s" % ("ok" if good else "FAIL", what))
        if not good:
            bad.append(what)

    shutil.rmtree("terrain", ignore_errors=True)
    os.makedirs("terrain")
    mt.write("terrain", DX, [r[1] for r in RUNS])
    # The deck is for 100 m cells and keeps a checkpoint: the check's copy has 200 m cells and no checkpoint
    # (a key given both in the deck and on the command line aborts).
    with open(cmd[deck_at[0]]) as fh:
        deck = fh.read()
    deck, n1 = re.subn(r"(?m)^mesh\.n_cell\s*=.*$", "mesh.n_cell = 80 64 40", deck)
    deck, n2 = re.subn(r"(?m)^io\.write_checkpoint\s*=.*$", "io.write_checkpoint = 0", deck)
    if n1 != 1 or n2 != 1:
        raise SystemExit("check_terrain_shapes: the deck has no single mesh.n_cell and io.write_checkpoint line")
    with open(os.path.join("terrain", "deck"), "w") as fh:
        fh.write(deck)
    run_cmd = list(cmd)
    run_cmd[deck_at[0]] = os.path.join("..", "terrain", "deck")

    logs, probes = {}, {}
    for name, shape, its in RUNS:
        shutil.rmtree(name, ignore_errors=True)
        os.makedirs(name)
        shutil.copy(os.path.join("terrain", "probes.txt"), name)
        full = run_cmd + ["terrain.file=../terrain/%s.xyz" % shape, "mesh.dz0=8",
                          "mesh.stretch_ratio=%.6f" % mt.ratio(8.0), "steady.cfl_max=16",
                          "steady.max_iterations=%d" % its]
        print("check_terrain_shapes: in %s/: %s" % (name, " ".join(full)))
        out = subprocess.run(full, cwd=name, check=False, capture_output=True, text=True)
        if out.returncode != 0:
            sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
            raise SystemExit("check_terrain_shapes: the run failed with code %d in %s/" % (out.returncode, name))
        logs[name] = out.stdout
        probes[name] = rows_of(os.path.join(name, "plt_steady_probes.csv"))

    m = re.search(r"gust: steady after (\d+) iterations", logs["gentle"])
    check(m is not None and int(m.group(1)) <= 250,
          "gentle: the steady verdict within 250 iterations (%s)" % (m.group(1) if m else "no verdict"))
    check("the ground steps by up to" not in logs["gentle"], "gentle: no ground-step warning")
    check("the ground steps by up to" in logs["steep"], "steep: the ground-step warning")
    check("gust: not steady after 200 iterations" in logs["steep"], "steep: no steady verdict at 200 iterations")

    rows = probes["steep"]
    check(len(rows) == 201 and len(rows[0][1]) == 2560, "steep: 201 probe rows of 2560 probes (%d of %d)"
          % (len(rows), len(rows[0][1]) if rows else 0))
    if len(rows) == 201:
        tail = [r[1] for r in rows[151:]]
        ref = [tuple(sum(row[p][c] for row in tail) / len(tail) for c in range(3)) for p in range(len(tail[0]))]
        for tol, limit in ((0.5, 40), (0.25, 100)):
            at = settled_from(rows, ref, tol)
            check(at is not None and at <= limit, "steep: every probe within %g m/s of the mean of rows 151 ... 200 "
                  "from iteration %s (<= %d)" % (tol, "never" if at is None else at, limit))
        began = error(rows[0][1], rows[200][1])
        check(began > 1.0, "steep: the rows are the run's: the first guess is %.2f m/s from row 200 somewhere (> 1)"
              % began)
        moved = error(rows[100][1], rows[200][1])
        check(moved < 0.2, "steep: the largest change at a probe between iterations 100 and 200 is %.4f m/s "
              "(< 0.2; 0.087 measured)" % moved)

    rows = probes["gentle"]
    at = settled_from(rows, rows[-1][1], 0.25) if rows else None
    check(at is not None and at <= 40, "gentle: every probe within 0.25 m/s of the last row from iteration %s (<= 40)"
          % ("never" if at is None else at))

    # The first probe is 300 m from the inflow face, upstream of the butte: the inflow's log law at its cell centre.
    with open(os.path.join("gentle", "plt_steady_probes.csv")) as fh:
        m = re.search(r"(?m)^# probe 1: .* -> \d+ \d+ \d+ (\S+)$", fh.read(20000))
    if m and rows:
        law = mt.USTAR / mt.KAPPA * math.log(float(m.group(1)) / mt.Z0)
        u = rows[-1][1][0][0]
        check(abs(u - law) < 0.02, "gentle: u at the first probe is %.4f m/s, the inflow log law at %.2f m is %.4f "
              "(within 0.02)" % (u, float(m.group(1)), law))
    else:
        check(False, "gentle: the probe file's header gives the first probe's cell-centre height")

    if bad:
        print("check_terrain_shapes: %d check(s) failed" % len(bad))
        return 1
    print("check_terrain_shapes: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
