#!/usr/bin/env python3
"""Runs the hill sweep of library_hill.yaml, builds a wind library and checks its blend (CTest rt_library_hill).

    check_library_hill.py --gust EXE [--mpi "mpiexec -n 2"] [-- key=value ...]

The run happens in the current directory: library_hill.yaml is copied there,
hill.asc is written there, the case folder is written with
tools/make_gust_case.py --strict and run.sh runs the five directions (260 to
280 degrees, 5 apart) with GUST=EXE and MPI=<--mpi>. The arguments after --
go to run.sh.

Checks:
- each run reached GUST's steady verdict;
- a library of 260 and 280 queried at 260 gives the 260 run's wind at 10 m,
  scaled to the speed asked (1e-12 of the speed): a stored sector comes back
  unchanged;
- the blend's wind at the reference point and 10 m has the speed asked
  (1e-12 relative);
- the blend at 270 of the 260 and 280 runs against the 270 run, at 10 m and
  50 m: the largest horizontal-wind difference more than 4 cells from the
  lateral faces under 1 % of the reference speed (0.67 % at 10 m in double
  precision; a blend without the renormalisation is 1.5 % off), and under
  0.4 % for the library of 265 and 275 (0.23 %). The faces are left out: a
  lateral face is an inflow face from 260 and 280 and a slip wall from 270,
  which gives 1.2 to 1.6 % in the outflow corners;
- the same library of 260 and 280 queried at 265 and 275 against those runs
  (1 %);
- k of the blend at 270 against the 270 run: the interior mean of
  k / k_270 - 1 within 0.5 % (+0.1 % at 10 m, +0.3 % at 50 m; scaling k by
  the renormalised speed squared gives +3 %).
The vertical wind w is printed, not checked (blend differences up to 0.02 m/s;
|w| reaches 1 m/s over the hill).
Measured in double precision; single precision is not measured here.
"""
import argparse
import math
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))
from gust_library import Library, build, index  # noqa: E402
from gust_library.sample import at_point  # noqa: E402

CASE = "library_hill"
SPEED = 8.0                                   # m/s asked of the blend at the reference point (the hilltop) and 10 m
HEIGHTS = (10.0, 50.0)                        # m above the ground
EDGE = 4                                      # cells left out next to each lateral face
TOL = {(260.0, 280.0): 0.01, (265.0, 275.0): 0.004}   # of the reference speed
TOL_K = 0.005                                 # interior mean of k / k_270 - 1

# The hill: a Gaussian on flat ground at 0 m, in the grid's frame (m).
HILL_X, HILL_Y, HILL_H, HILL_SIGMA = 2000.0, 2000.0, 60.0, 200.0
GRID_X0, GRID_Y0, GRID_CELL, GRID_N = 0.0, 0.0, 20.0, 201
FAILED = []


def hill_height(x, y):
    return HILL_H * math.exp(-((x - HILL_X) ** 2 + (y - HILL_Y) ** 2) / (2.0 * HILL_SIGMA ** 2))


def write_hill(path):
    """hill.asc: an ESRI ASCII grid of the hill, cell corners from (GRID_X0, GRID_Y0), north row first."""
    with open(path, "w") as fh:
        fh.write("ncols %d\nnrows %d\nxllcorner %.1f\nyllcorner %.1f\ncellsize %.1f\nNODATA_value -9999\n"
                 % (GRID_N, GRID_N, GRID_X0, GRID_Y0, GRID_CELL))
        for j in reversed(range(GRID_N)):
            y = GRID_Y0 + (j + 0.5) * GRID_CELL
            fh.write(" ".join("%.4f" % hill_height(GRID_X0 + (i + 0.5) * GRID_CELL, y) for i in range(GRID_N)))
            fh.write("\n")


def check(name, value, target, tol, kind):
    err = abs(value - target) / abs(target) if kind == "rel" else abs(value - target)
    ok = math.isfinite(value) and err <= tol
    print("%s %-66s %.6g (target %.6g, %s error %.3g, tolerance %.3g)"
          % ("PASS" if ok else "FAIL", name, value, target, kind, err, tol))
    if not ok:
        FAILED.append(name)


def run_fields(lib, folder):
    """The run's wind at the library's heights, scaled to SPEED at the library's reference point (k to SPEED**2)."""
    _, stored, _, _ = index.sector(os.path.join(CASE, folder), list(HEIGHTS), tuple(lib.index["reference_point"]))
    out = {f: SPEED * stored[f] for f in ("u", "v", "w")}
    out["k"] = SPEED ** 2 * stored["k"]
    return out


def compare(tag, lib, direction, folder, tol, k_check=False):
    """The blend at direction against the run in folder, at every height."""
    run = run_fields(lib, folder)
    for n, h in enumerate(HEIGHTS):
        w = lib.wind(direction, SPEED, h)
        du = np.hypot(w["u"] - run["u"][n], w["v"] - run["v"][n])
        inner = du[EDGE:-EDGE, EDGE:-EDGE] / SPEED
        dw = np.abs(w["w"] - run["w"][n])[EDGE:-EDGE, EDGE:-EDGE]
        print("     %s at %g, %g m: interior RMS %.3g %%, all cells max %.3g %%; interior max |dw| %.3g m/s "
              "(run max |w| %.3g m/s)" % (tag, direction, h, 100.0 * math.sqrt((inner ** 2).mean()),
                                         100.0 * du.max() / SPEED, dw.max(), np.abs(run["w"][n]).max()))
        check("%s at %g, %g m: interior max |dU| / reference speed" % (tag, direction, h), float(inner.max()),
              0.0, tol, "abs")
        if k_check:
            bias = float((w["k"] / run["k"][n])[EDGE:-EDGE, EDGE:-EDGE].mean() - 1.0)
            check("%s at %g, %g m: interior mean k / k_run - 1" % (tag, direction, h), bias, 0.0, TOL_K, "abs")


def main():
    argv = sys.argv[1:]
    extra = []
    if "--" in argv:
        n = argv.index("--")
        argv, extra = argv[:n], argv[n + 1:]
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--gust", required=True, help="the GUST executable")
    ap.add_argument("--mpi", default="", help="the launcher run.sh puts before GUST (empty: none)")
    opts = ap.parse_args(argv)

    shutil.copyfile(os.path.join(HERE, CASE + ".yaml"), CASE + ".yaml")
    write_hill("hill.asc")
    if os.path.isdir(CASE):
        shutil.rmtree(CASE)
    made = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "make_gust_case.py"), CASE + ".yaml",
                           "--quiet", "--strict"])
    check("make_gust_case.py --strict exit status", float(made.returncode), 0.0, 0.0, "abs")
    if made.returncode != 0:
        return 1
    env = dict(os.environ, GUST=os.path.abspath(opts.gust), MPI=opts.mpi)
    ran = subprocess.run(["bash", os.path.join(CASE, "run.sh")] + extra, env=env)
    check("run.sh exit status", float(ran.returncode), 0.0, 0.0, "abs")
    if ran.returncode != 0:
        return 1
    for d in (260, 265, 270, 275, 280):
        try:
            status = index.read_log(os.path.join(CASE, "dir_%03d" % d))[0]
        except index.LibraryError as err:
            print(err)
            status = None
        check("from %d: GUST's steady verdict" % d, float(status == "steady"), 1.0, 0.0, "abs")
    if FAILED:
        return 1

    for pair, tol in TOL.items():
        build(CASE, HEIGHTS, None, pair)
        lib = Library(CASE)
        tag = "library %g/%g" % pair
        if pair == (260.0, 280.0):
            w = lib.wind(260.0, SPEED, HEIGHTS[0])
            run = run_fields(lib, "dir_260")
            diff = max(np.abs(w["u"] - run["u"][0]).max(), np.abs(w["v"] - run["v"][0]).max())
            check(tag + " at 260: the 260 run's wind (m/s)", diff, 0.0, 1e-12 * SPEED, "abs")
        x, y = lib.index["reference_point"]
        w = lib.wind(270.0, SPEED, HEIGHTS[0])
        s_ref = math.hypot(float(at_point(lib.x, lib.y, w["u"], x, y)), float(at_point(lib.x, lib.y, w["v"], x, y)))
        check(tag + " at 270: speed at the reference point (m/s)", s_ref, SPEED, 1e-12, "rel")
        compare(tag, lib, 270.0, "dir_270", tol, k_check=True)
        if pair == (260.0, 280.0):
            for d in (265.0, 275.0):
                compare(tag, lib, d, "dir_%03d" % d, 0.01)

    if FAILED:
        print("check_library_hill: %d check(s) failed: %s" % (len(FAILED), "; ".join(FAILED)))
        return 1
    print("check_library_hill: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
