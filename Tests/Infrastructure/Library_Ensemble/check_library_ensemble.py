#!/usr/bin/env python3
"""Runs the hill sweep of library_ensemble.yaml, builds a wind library and checks its ensemble statistics and
WindNinja-format grids (CTest rt_library_ensemble).

    check_library_ensemble.py --gust EXE [--mpi "mpiexec -n 4"] [-- key=value ...]

The run happens in the current directory: library_ensemble.yaml is copied there, hill.asc is written there (UTM
metres of domain.center_lat_lon's zone), the case folder is written with tools/make_gust_case.py --strict and
run.sh runs the four directions (0, 90, 180 and 270 degrees) with GUST=EXE and MPI=<--mpi>. The arguments after
-- go to run.sh.

Checks:
- each run reached GUST's steady verdict;
- make_gust_library.py ensemble over 10 members (the directions 0, 45, 90, 180 and 270 degrees with weights 0.1,
  0.15, 0.2, 0.25, 0.3, each at 6 and 10 m/s with weights 0.5, 0.5) at 10 m and 50 m: the mean speed, its standard
  deviation, the 10th, 50th and 90th percentiles, the probability of a speed above 8 m/s and the probability of a
  reversal in every column equal the same statistics computed here, cell by cell, from the runs' own plotfiles
  (each run's wind scaled to the member's speed at the reference point; 45 degrees, between two runs, blended
  here by hand), within 1e-12 of the mean speed (probabilities within 1e-12);
- the _vel and _ang grids it writes: header (40 columns and rows, cell size 50 m, NODATA_value -9999), corner
  at the domain's south-west corner in UTM (the easting and northing of domain.center_lat_lon less 1000 m,
  within 1e-6 m), north row first, values equal to the mean speed and the direction of the vector mean (1e-7
  relative and 1e-6 degrees);
- the .prj names the UTM zone of domain.center_lat_lon (13N, central meridian -105);
- the library's ground written as a grid with the same corner has its centroid (weights: the ground above its
  lowest value) at the hill's top in UTM within 1 m: the corner puts the fields on the map where the terrain
  came from;
- the CSV at points: at the hilltop the values equal the statistics of the members' winds there.
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
from gust_case import geo  # noqa: E402
from gust_library import Library, build, grids, index  # noqa: E402
from gust_library.sample import at_point  # noqa: E402

CASE = "library_ensemble"
LAT, LON = 40.0, -105.5
HALF_WIDTH, CELL = 1000.0, 50.0
HEIGHTS = (10.0, 50.0)
DIRECTIONS = (0.0, 90.0, 180.0, 270.0)                 # the runs
MEMBER_DIRECTIONS = (0.0, 45.0, 90.0, 180.0, 270.0)      # 45: between two runs, blended
DIRECTION_WEIGHTS = (0.1, 0.15, 0.2, 0.25, 0.3)
SPEEDS = (6.0, 10.0)                        # m/s at the reference point and 10 m
PERCENTILES = (10.0, 50.0, 90.0)
THRESHOLD = 8.0                             # m/s
TOL = 1.0e-12

# The hill: a Gaussian 60 m high, 200 m wide, at the domain centre; the grid in UTM metres around it.
HILL_H, HILL_SIGMA = 60.0, 200.0
GRID_CELL, GRID_HALF = 20.0, 1600.0
FAILED = []


def check(name, value, target, tol, kind="abs"):
    err = abs(value - target) / abs(target) if kind == "rel" else abs(value - target)
    ok = math.isfinite(value) and err <= tol
    print("%s %-74s %.9g (target %.9g, %s error %.3g, tolerance %.3g)"
          % ("PASS" if ok else "FAIL", name, value, target, kind, err, tol))
    if not ok:
        FAILED.append(name)


def write_hill(path, east, north):
    """hill.asc: the hill at (east, north), UTM metres, north row first."""
    n = int(round(2 * GRID_HALF / GRID_CELL))
    x0, y0 = east - GRID_HALF, north - GRID_HALF
    with open(path, "w") as fh:
        fh.write("ncols %d\nnrows %d\nxllcorner %.3f\nyllcorner %.3f\ncellsize %.1f\nNODATA_value -9999\n"
                 % (n, n, x0, y0, GRID_CELL))
        for j in reversed(range(n)):
            y = y0 + (j + 0.5) * GRID_CELL - north
            x = x0 + (np.arange(n) + 0.5) * GRID_CELL - east
            fh.write(" ".join("%.4f" % h for h in HILL_H * np.exp(-(x ** 2 + y ** 2) / (2.0 * HILL_SIGMA ** 2))))
            fh.write("\n")


def read_asc(path):
    """Header dict and the values as written (rows north first), read here without gust_library."""
    with open(path) as fh:
        header = {}
        for _ in range(6):
            k, v = fh.readline().split()
            header[k] = float(v)
        rows = np.array([[float(a) for a in line.split()] for line in fh if line.strip()])
    return header, rows


def members():
    return [(d, s, wd * ws) for d, wd in zip(MEMBER_DIRECTIONS, DIRECTION_WEIGHTS) for s, ws in zip(SPEEDS, (0.5, 0.5))]


def unit_wind(stored, d, height_index, xc, yc, ref):
    """u, v (nx, ny) for a wind from d of 1 m/s at the reference point and the first height. A run's direction is
    its own fields (stored per unit reference speed). 45 degrees is the half-and-half mix of the runs from 0 and 90,
    divided by the length of the mix of their unit wind vectors at the reference point (the vectors are read here
    from the fields, not from the library's index)."""
    if d in stored:
        return stored[d]["u"][height_index], stored[d]["v"][height_index]
    assert d == 45.0
    a, b = stored[0.0], stored[90.0]
    ex = 0.5 * at_point(xc, yc, a["u"][0], *ref) + 0.5 * at_point(xc, yc, b["u"][0], *ref)
    ey = 0.5 * at_point(xc, yc, a["v"][0], *ref) + 0.5 * at_point(xc, yc, b["v"][0], *ref)
    length = math.hypot(ex, ey)
    return ((0.5 * a["u"][height_index] + 0.5 * b["u"][height_index]) / length,
            (0.5 * a["v"][height_index] + 0.5 * b["v"][height_index]) / length)


def member_winds(stored, height_index, xc, yc, ref):
    """u, v (members, nx, ny): each member's unit wind scaled by its speed."""
    unit = {d: unit_wind(stored, d, height_index, xc, yc, ref) for d in MEMBER_DIRECTIONS}
    u = np.array([s * unit[d][0] for d, s, _ in members()])
    v = np.array([s * unit[d][1] for d, s, _ in members()])
    return u, v


def reference_statistics(u, v):
    """The statistics, cell by cell with plain loops (an implementation separate from gust_library.ensemble)."""
    ms = members()
    w = [m[2] for m in ms]
    out = {k: np.zeros(u.shape[1:]) for k in ("mean_speed", "std_speed", "pexceed", "preversal", "mean_direction")}
    for p in PERCENTILES:
        out["p%g" % p] = np.zeros(u.shape[1:])
    for i in range(u.shape[1]):
        for j in range(u.shape[2]):
            s = [math.hypot(u[m, i, j], v[m, i, j]) for m in range(len(ms))]
            mean = sum(wm * sm for wm, sm in zip(w, s))
            out["mean_speed"][i, j] = mean
            out["std_speed"][i, j] = math.sqrt(sum(wm * (sm - mean) ** 2 for wm, sm in zip(w, s)))
            out["pexceed"][i, j] = sum(wm for wm, sm in zip(w, s) if sm > THRESHOLD)
            rev = 0.0
            for m, (d, _, wm) in enumerate(ms):
                ex, ey = -math.sin(math.radians(d)), -math.cos(math.radians(d))   # towards which the inflow blows
                if u[m, i, j] * ex + v[m, i, j] * ey < 0.0:
                    rev += wm
            out["preversal"][i, j] = rev
            mu = sum(wm * u[m, i, j] for m, wm in enumerate(w))
            mv = sum(wm * v[m, i, j] for m, wm in enumerate(w))
            out["mean_direction"][i, j] = math.degrees(math.atan2(-mu, -mv)) % 360.0
            ranked = sorted(zip(s, w))
            tol = max(1e-12, 4.0 * len(ms) * sys.float_info.epsilon)  # the cumulative weight rounds by about n eps
            for p in PERCENTILES:
                cum = 0.0
                for sm, wm in ranked:
                    cum += wm
                    if cum >= p / 100.0 - tol:
                        out["p%g" % p][i, j] = sm
                        break
    return out


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

    east, north, zone, south = geo.latlon_to_utm(LAT, LON)
    east, north = float(east), float(north)
    shutil.copyfile(os.path.join(HERE, CASE + ".yaml"), CASE + ".yaml")
    write_hill("hill.asc", east, north)
    if os.path.isdir(CASE):
        shutil.rmtree(CASE)
    made = subprocess.run([sys.executable, os.path.join(ROOT, "tools", "make_gust_case.py"), CASE + ".yaml",
                           "--quiet", "--strict"])
    check("make_gust_case.py --strict exit status", float(made.returncode), 0.0, 0.0)
    if made.returncode != 0:
        return 1
    env = dict(os.environ, GUST=os.path.abspath(opts.gust), MPI=opts.mpi)
    ran = subprocess.run(["bash", os.path.join(CASE, "run.sh")] + extra, env=env)
    check("run.sh exit status", float(ran.returncode), 0.0, 0.0)
    if ran.returncode != 0:
        return 1
    for d in DIRECTIONS:
        try:
            status = index.read_log(os.path.join(CASE, "dir_%03d" % d))[0]
        except index.LibraryError as err:
            print(err)
            status = None
        check("from %d: GUST's steady verdict" % d, float(status == "steady"), 1.0, 0.0)
    if FAILED:
        return 1

    build(CASE, HEIGHTS, None, None, quiet=True)
    lib = Library(CASE)
    ref = tuple(lib.index["reference_point"])
    stored = {}
    for d in DIRECTIONS:
        _, fields, grid, _ = index.sector(os.path.join(CASE, "dir_%03d" % d), list(HEIGHTS), ref)
        stored[d] = fields
    xc, yc = grid["x"], grid["y"]
    tool = os.path.join(ROOT, "tools", "make_gust_library.py")
    with open("points.txt", "w") as fh:
        fh.write("%r %r hilltop\n" % ref)
    for n, h in enumerate(HEIGHTS):
        out = "grids_%gm" % h
        cmd = [sys.executable, tool, "ensemble", CASE, "--height", "%g" % h,
               "--directions"] + ["%g" % d for d in MEMBER_DIRECTIONS] + ["--direction-weights"] + \
              ["%g" % w for w in DIRECTION_WEIGHTS] + ["--speeds"] + ["%g" % s for s in SPEEDS] + \
              ["--thresholds", "%g" % THRESHOLD, "--grids", out, "--points", "points.txt", "--out", out + ".npz"]
        rc = subprocess.run(cmd).returncode
        check("ensemble at %g m: exit status" % h, float(rc), 0.0, 0.0)
        if rc != 0:
            continue
        st = dict(np.load(out + ".npz"))
        u, v = member_winds(stored, n, xc, yc, ref)
        want = reference_statistics(u, v)
        scale = float(want["mean_speed"].max())
        for name, key, tol in (("mean speed", "mean_speed", TOL * scale), ("standard deviation", "std_speed", TOL * scale),
                               ("p10", "p10", TOL * scale), ("p50", "p50", TOL * scale), ("p90", "p90", TOL * scale),
                               ("probability above 8 m/s", "pexceed_8", TOL), ("probability of a reversal", "preversal", TOL)):
            got = st[key]
            ref_ = want["pexceed" if key.startswith("pexceed") else key]
            check("ensemble at %g m: %s, max |tool - runs|" % (h, name), float(np.abs(got - ref_).max()), 0.0, tol)
        print("     at %g m: mean speed %.4g to %.4g m/s, standard deviation up to %.4g m/s, probability above 8 m/s "
              "%.3g to %.3g, probability of a reversal up to %.3g"
              % (h, want["mean_speed"].min(), want["mean_speed"].max(), want["std_speed"].max(),
                 want["pexceed"].min(), want["pexceed"].max(), want["preversal"].max()))
        stem = os.path.join(out, "%s_ensemble_%sm_agl" % (CASE, ("%g" % h).replace(".", "p")))
        for kind, field, tol, how in (("vel", want["mean_speed"], 1e-7, "rel"), ("ang", want["mean_direction"], 1e-6, "deg")):
            header, rows = read_asc(stem + "_%s.asc" % kind)
            check("_%s at %g m: ncols, nrows" % (kind, h), header["ncols"] + 1000 * header["nrows"], 40 + 40000, 0.0)
            check("_%s at %g m: cellsize (m)" % (kind, h), header["cellsize"], CELL, 0.0)
            check("_%s at %g m: NODATA_value" % (kind, h), header["NODATA_value"], -9999.0, 0.0)
            check("_%s at %g m: xllcorner (m, UTM)" % (kind, h), header["xllcorner"], east - HALF_WIDTH, 1e-6)
            check("_%s at %g m: yllcorner (m, UTM)" % (kind, h), header["yllcorner"], north - HALF_WIDTH, 1e-6)
            values = rows[::-1].T                                  # (east, north)
            if how == "rel":
                err = float(np.abs(values / field - 1.0).max())
            else:
                err = float(np.abs((values - field + 180.0) % 360.0 - 180.0).max())
            check("_%s at %g m: values, north row first (max %s error)" % (kind, h, how), err, 0.0, tol)
        with open(stem + "_vel.prj") as fh:
            prj = fh.read()
        check("_vel.prj at %g m: names WGS_1984_UTM_Zone_13N" % h,
              float('PROJCS["WGS_1984_UTM_Zone_%d%s"' % (zone, "S" if south else "N") in prj and zone == 13), 1.0, 0.0)
        check("_vel.prj at %g m: central meridian -105" % h, float('"Central_Meridian",-105.0]' in prj), 1.0, 0.0)
        with open(stem + "_points.csv") as fh:
            lines = [line.rstrip("\n").split(",") for line in fh]
        col = lines[0].index("mean_speed")
        pu = np.array([at_point(lib.x, lib.y, u[m], *ref) for m in range(u.shape[0])])
        pv = np.array([at_point(lib.x, lib.y, v[m], *ref) for m in range(u.shape[0])])
        at_top = reference_statistics(pu[:, None, None], pv[:, None, None])["mean_speed"][0, 0]
        check("CSV at %g m: mean speed at the hilltop (m/s)" % h, float(lines[1][col]), at_top, 1e-8, "rel")
        check("CSV at %g m: easting of the hilltop (m)" % h, float(lines[1][lines[0].index("easting")]),
              east + ref[0], 1e-6)

    geo_ = grids.georeference(lib)
    grids.write_grid("ground.asc", lib.x, lib.y, lib.ground, geo_)
    header, rows = read_asc("ground.asc")
    # the hill is symmetric about its top: the centroid of the ground above its lowest value is the top
    weight = rows - rows.min()
    ce = header["xllcorner"] + (np.arange(rows.shape[1]) + 0.5) * header["cellsize"]
    cn = header["yllcorner"] + (rows.shape[0] - np.arange(rows.shape[0]) - 0.5) * header["cellsize"]
    top_e = float((weight * ce[None, :]).sum() / weight.sum())
    top_n = float((weight * cn[:, None]).sum() / weight.sum())
    check("ground grid: easting of the hill's centroid vs its top (m)", top_e, east, 1.0)
    check("ground grid: northing of the hill's centroid vs its top (m)", top_n, north, 1.0)

    if FAILED:
        print("check_library_ensemble: %d check(s) failed: %s" % (len(FAILED), "; ".join(FAILED)))
        return 1
    print("check_library_ensemble: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
