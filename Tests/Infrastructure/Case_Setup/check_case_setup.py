#!/usr/bin/env python3
"""Writes a case with tools/make_gust_case.py, runs its run.sh and checks the result.

    check_case_setup.py CASE --gust EXE [--mpi "mpiexec -n 2"] [-- key=value ...]

CASE is flat_sweep, flat_nest, column_target, column_axis, hill or upstream (the YAML file
of that name beside this script). The run happens in the current directory: the YAML
file is copied there (hill.asc or regional.xyz is written there), the case folder
is written with --overwrite, and run.sh runs with GUST=EXE and MPI=<--mpi>
(empty: no launcher). The arguments after -- go to run.sh, which passes them
to the 3-D runs.

The checks pin closed forms, not the code's own output:
- flat_sweep, flat_nest: the planar-mean wind at wind.height is the log law's
  wind.speed (3e-4) along the asked direction (0.05 degrees), GUST reached
  its steady verdict, and the nested patch has half the parent's cell size;
- column_target: the column's table holds wind.speed at wind.height (2e-4,
  the column's stopping rule |R - 1| < 1e-4) along the asked direction, and
  the 3-D run that took the column's geostrophic wind meets it to 1 % and
  0.3 degrees, and its own column first guess, solved from the forcing
  run.sh passed, has the table's wind at the top; library.json (run.sh's
  --collect) holds wind.speed at wind.height (2e-4) along the asked direction
  (0.01 degrees: the stopping rule allows 0.0057), G equal to the column log's
  last forcing.ug and forcing.vg (the values run.sh passes on), and G turned
  clockwise of the 10 m wind by 10 to 25 degrees (between half the neutral
  drag law's angle and all of it, 23.5 degrees);
- column_axis: with the 10 m wind along y and along x, every cell's wind of the
  3-D run is the column table's at its height to 4 % of the top wind (a face
  parallel to the 10 m wind made a slip wall gives 73 % from 0, 72 % from 270), and no lateral face
  of the decks is a slip wall; library.json holds wind.speed at wind.height along the asked
  direction from 0, 270 and 240 (a direction off the axes: a quarter turn is exact in floating
  point), and the same |G| and turning angle for the three (the neutral column is
  rotation-invariant: 1e-6 relative and 1e-4 degrees in double precision; 2e-4 and 0.0115
  degrees, the bound of the stopping rule, with --precision SINGLE);
- hill: GUST's ground follows the hill's formula (crest in the right cell,
  relief to 2 %, RMS difference below 0.3 % of the relief);
- upstream: wind.upstream_roughness = auto over a half-plane regional map that
  also sets the ground: the wind from 270 degrees takes the map's 0.03 m west of
  the domain (the decks' bc.inflow_z0 and the column's bc.surface_z0), and the
  10 m wind more than 4 cells from the x faces varies by less than 1 % along x.
"""
import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "extern_erf"))
import erf_plotfile  # noqa: E402

SPEED, HEIGHT, Z0 = 8.0, 10.0, 0.03          # wind.speed, wind.height, surface.roughness of every YAML here
FAILED = []

# The hill of hill.yaml, in the grid's projected frame (m).
HILL_X, HILL_Y = 601700.0, 4101500.0         # 100 m east and 100 m south of the domain centre
HILL_H, HILL_SIGMA, HILL_SLOPE = 120.0, 300.0, 0.01
GRID_X0, GRID_Y0, GRID_CELL, GRID_N = 600000.0, 4100000.0, 20.0, 161
CENTER_X, CENTER_Y, HALF_WIDTH = 601600.0, 4101600.0, 1200.0   # domain.center, domain.half_width


def hill_height(x, y):
    """The hill's elevation (m): a Gaussian on a plane rising to the east."""
    r2 = (x - HILL_X) ** 2 + (y - HILL_Y) ** 2
    return 1000.0 + HILL_H * math.exp(-r2 / (2.0 * HILL_SIGMA ** 2)) + HILL_SLOPE * (x - GRID_X0)


def write_hill(path):
    """hill.asc: an ESRI ASCII grid of the hill, cell corners from (GRID_X0, GRID_Y0), north row first."""
    with open(path, "w") as fh:
        fh.write("ncols %d\nnrows %d\nxllcorner %.1f\nyllcorner %.1f\ncellsize %.1f\nNODATA_value -9999\n"
                 % (GRID_N, GRID_N, GRID_X0, GRID_Y0, GRID_CELL))
        for j in reversed(range(GRID_N)):
            y = GRID_Y0 + (j + 0.5) * GRID_CELL
            fh.write(" ".join("%.4f" % hill_height(GRID_X0 + (i + 0.5) * GRID_CELL, y) for i in range(GRID_N)))
            fh.write("\n")


# The regional map of upstream.yaml: z0 = Z0 west of x = UPSTREAM_STEP, UPSTREAM_EAST east (m).
UPSTREAM_STEP, UPSTREAM_EAST, UPSTREAM_REACH, UPSTREAM_SPACING = 1200.0, 1.0, 8000.0, 200.0


def write_regional(path):
    """regional.xyz: x y z0 points every UPSTREAM_SPACING out to UPSTREAM_REACH from the domain centre."""
    n = int(round(UPSTREAM_REACH / UPSTREAM_SPACING))
    with open(path, "w") as fh:
        for j in range(-n, n + 1):
            for i in range(-n, n + 1):
                x = i * UPSTREAM_SPACING
                fh.write("%.1f %.1f %g\n" % (x, j * UPSTREAM_SPACING, Z0 if x < UPSTREAM_STEP else UPSTREAM_EAST))


def deck_value(path, key):
    with open(path) as fh:
        for ln in fh:
            parts = ln.split()
            if len(parts) >= 3 and parts[0] == key and parts[1] == "=":
                return float(parts[2])
    return float("nan")


def check_upstream(folder):
    """The decks take the map's z0 upwind of the inflow face; the 10 m wind (linear in ln(z + z0) between the cell
    centres around it, mean over y) away from the inflow and outflow faces holds along x."""
    check("bc.inflow_z0 of the 3-D deck (m)", deck_value(os.path.join(folder, "inputs"), "bc.inflow_z0"), Z0,
          1e-9, "rel")
    check("bc.surface_z0 of the column deck (m)",
          deck_value(os.path.join(folder, "column", "inputs"), "bc.surface_z0"), Z0, 1e-9, "rel")
    check("GUST's steady verdict", steady(folder), 1.0, 0.0, "abs")
    hdr, data = erf_plotfile.read_fields(last_plotfile(folder), ["u", "v", "z"])
    along_x = []
    for i, plane in enumerate(data["z"]):
        winds = [log_interp(column, [math.hypot(a, b) for a, b in zip(data["u"][i][j], data["v"][i][j])], HEIGHT)
                 for j, column in enumerate(plane)]
        along_x.append(sum(winds) / len(winds))
    inner = along_x[4:-4]
    check("(max - min) / max of the 10 m wind along x, 4 cells from the x faces", (max(inner) - min(inner)) / max(inner),
          0.0, 1e-2, "abs")


def check(name, value, target, tol, kind):
    err = abs(value - target) / abs(target) if kind == "rel" else abs(value - target)
    ok = math.isfinite(value) and err <= tol
    print("%s %-62s %.6g (target %.6g, %s error %.3g, tolerance %.3g)"
          % ("PASS" if ok else "FAIL", name, value, target, kind, err, tol))
    if not ok:
        FAILED.append(name)


def towards(direction_from):
    """GUST's flow direction (degrees counter-clockwise from x) of a meteorological direction."""
    return (270.0 - direction_from) % 360.0


def angle_error(a, b):
    return abs((a - b + 180.0) % 360.0 - 180.0)


def log_interp(z, f, zt):
    """f at height zt, linear in ln(z + z0) between the two levels around it."""
    k = max(j for j in range(len(z) - 1) if z[j] <= zt)
    t = (math.log(zt + Z0) - math.log(z[k] + Z0)) / (math.log(z[k + 1] + Z0) - math.log(z[k] + Z0))
    return (1.0 - t) * f[k] + t * f[k + 1]


def last_plotfile(folder):
    plts = sorted(d for d in os.listdir(folder) if d.startswith("plt") and d[3:].isdigit())
    if not plts:
        raise SystemExit("check_case_setup: no plotfile in %s" % folder)
    return os.path.join(folder, plts[-1])


def wind_at_height(folder):
    """(speed, direction) of the planar-mean wind at HEIGHT on flat ground; the dx of the plotfile."""
    z, avg, hdr = erf_plotfile.planar_averages(last_plotfile(folder), ["u", "v", "z"])
    u = log_interp(avg["z"], avg["u"], HEIGHT)
    v = log_interp(avg["z"], avg["v"], HEIGHT)
    return math.hypot(u, v), math.degrees(math.atan2(v, u)) % 360.0, hdr["dx"]


def steady(folder):
    with open(os.path.join(folder, "log.txt")) as fh:
        return 1.0 if any(ln.startswith("gust: steady after") for ln in fh) else 0.0


def check_flat(folder, direction_from, label):
    speed, direction, dx = wind_at_height(folder)
    check("%s: speed at %g m (m/s)" % (label, HEIGHT), speed, SPEED, 3e-4, "rel")
    check("%s: direction error (degrees)" % label, angle_error(direction, towards(direction_from)), 0.0, 0.05, "abs")
    check("%s: GUST's steady verdict" % label, steady(folder), 1.0, 0.0, "abs")
    return dx


def column_table(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            ln = ln.split("#")[0].split()
            if ln:
                rows.append([float(x) for x in ln])
    return rows


def check_column_table(path, direction_from):
    rows = column_table(path)
    z = [r[0] for r in rows]
    u = log_interp(z[1:], [r[1] for r in rows[1:]], HEIGHT)
    v = log_interp(z[1:], [r[2] for r in rows[1:]], HEIGHT)
    check("column table: speed at %g m (m/s)" % HEIGHT, math.hypot(u, v), SPEED, 2e-4, "rel")
    check("column table: direction error (degrees)",
          angle_error(math.degrees(math.atan2(v, u)), towards(direction_from)), 0.0, 0.05, "abs")
    return rows[-1][1], rows[-1][2]


def check_column_handover(log, top_wind):
    """The 3-D run solves its own column first guess from the forcing.ug/vg run.sh passed it. The same
    1-D problem with the same geostrophic wind gives the column table's wind at the top (the log prints
    4 decimals). The 10 m wind of a small 3-D run is set by the inflow table, so only this check sees
    what run.sh passed."""
    pattern = re.compile(r"column guess: converged .* u (\S+) v (\S+) m/s at the top")
    with open(log) as fh:
        found = [pattern.search(ln) for ln in fh]
    found = [m for m in found if m]
    u, v = (float(found[-1].group(1)), float(found[-1].group(2))) if found else (float("nan"),) * 2
    error = math.hypot(u - top_wind[0], v - top_wind[1]) / math.hypot(*top_wind)
    check("3-D column guess: top wind against the column table (relative)", error, 0.0, 1e-4, "abs")


def check_column_axis(folder, label):
    """Flat ground: every cell's (u, v) against the column table at the cell's height, relative to the table's
    top wind. A face that blocks the flow aloft (a slip wall) distorts the whole field."""
    with open(os.path.join(folder, "inputs")) as fh:
        faces = {ln.split()[0]: ln.split()[2] for ln in fh if re.match(r"bc\.[xy](lo|hi) ", ln)}
    check("%s: lateral faces that are slip walls" % label, float(sum(v == "slip_wall" for v in faces.values())),
          0.0, 0.0, "abs")
    rows = column_table(os.path.join(folder, "column", "column_guess.dat"))[1:]
    zc = [r[0] for r in rows]
    top = math.hypot(rows[-1][1], rows[-1][2])
    hdr, data = erf_plotfile.read_fields(last_plotfile(folder), ["u", "v", "z"])
    worst = 0.0
    for i, plane in enumerate(data["z"]):
        for j, column in enumerate(plane):
            for k, z in enumerate(column):
                n = max([0] + [m for m in range(len(zc) - 1) if zc[m] <= z])
                n = min(n, len(zc) - 2)
                t = min(max((z - zc[n]) / (zc[n + 1] - zc[n]), 0.0), 1.0)
                uc = (1.0 - t) * rows[n][1] + t * rows[n + 1][1]
                vc = (1.0 - t) * rows[n][2] + t * rows[n + 1][2]
                worst = max(worst, math.hypot(data["u"][i][j][k] - uc, data["v"][i][j][k] - vc) / top)
    check("%s: max |wind - column table| / top wind" % label, worst, 0.0, 4e-2, "abs")


def log_forcing(log):
    """The last forcing.ug and forcing.vg lines of a column log, as run.sh's awk reads them."""
    found = {}
    with open(log) as fh:
        for line in fh:
            parts = line.split()
            if len(parts) >= 3 and parts[0] in ("forcing.ug", "forcing.vg") and parts[1] == "=":
                found[parts[0]] = float(parts[2])
    return found["forcing.ug"], found["forcing.vg"]


def check_library(case_dir, directions, single=False):
    """library.json: the reference wind is wind.speed from each asked direction, G is turned clockwise of it, and
    |G| and the turning angle do not depend on the direction. GUST stops adjusting each column at |R - 1| < 1e-4
    (R: the target over the column's wind), so each G is within 1e-4 (relative) and 1e-4 rad of the exact one:
    two directions differ by at most 2e-4 and 0.0115 degrees. In double precision the columns agree far closer
    (1e-9 measured), so the checks there are 1e-6 and 1e-4 degrees; in single precision the bound."""
    path = os.path.join(case_dir, "library.json")
    if not os.path.isfile(path):
        check("library.json written by run.sh", 1.0, 0.0, 0.0, "abs")
        return
    g_tol, turn_tol = (2e-4, 0.0115) if single else (1e-6, 1e-4)
    with open(path) as fh:
        entries = json.load(fh)["directions"]
    check("library.json: directions", float(len(entries)), float(len(directions)), 0.0, "abs")
    for e, d in zip(entries, directions):
        try:
            check("library from %g: run folder" % d, float(e["folder"] in ("dir_%03d" % d, ".")), 1.0, 0.0, "abs")
            check("library from %g: speed at %g m (m/s)" % (d, HEIGHT), e["wind_at_height"]["speed"], SPEED, 2e-4,
                  "rel")
            check("library from %g: direction error at %g m (degrees)" % (d, HEIGHT),
                  angle_error(e["wind_at_height"]["direction"], d), 0.0, 0.01, "abs")
            # G turned clockwise by between half the neutral drag law's angle and all of it (alpha = 23.5 degrees
            # for this u* at 45 N; the faces of the decks are chosen from alpha / 2): a first guess along the
            # asked wind, never adjusted, turns 0
            check("library from %g: G turned clockwise of the reference wind by 10 to 25 degrees" % d,
                  float(-25.0 < e["turning_deg"] < -10.0), 1.0, 0.0, "abs")
            ug, vg = log_forcing(os.path.join(case_dir, e["folder"], "column", "log.txt"))
            check("library from %g: G = the column log's last forcing.ug, forcing.vg (run.sh's)" % d,
                  math.hypot(e["geostrophic_wind"]["ug"] - ug, e["geostrophic_wind"]["vg"] - vg), 0.0, 0.0, "abs")
            check("library from %g: |G| against the first direction's" % d, e["geostrophic_wind"]["speed"],
                  entries[0]["geostrophic_wind"]["speed"], g_tol, "rel")
            check("library from %g: turning against the first direction's (degrees)" % d,
                  abs(e["turning_deg"] - entries[0]["turning_deg"]), 0.0, turn_tol, "abs")
        except (KeyError, TypeError) as err:
            check("library from %g: entry complete (%s: %s)" % (d, type(err).__name__, err), 1.0, 0.0, 0.0, "abs")


def check_hill(folder, top, first_layer):
    """The ground from the first cell centre: z = zeta + h (1 - zeta / top) (terrain.smoothing = 0)."""
    hdr, data = erf_plotfile.read_fields(last_plotfile(folder), ["z"])
    zeta0 = 0.5 * first_layer
    x0, y0 = CENTER_X - HALF_WIDTH, CENTER_Y - HALF_WIDTH          # the domain's south-west corner
    dx, dy = hdr["dx"][0], hdr["dx"][1]
    nx, ny = len(data["z"]), len(data["z"][0])
    ground, exact = [], []
    for i in range(nx):
        for j in range(ny):
            ground.append((data["z"][i][j][0] - zeta0) / (1.0 - zeta0 / top))
            exact.append(hill_height(x0 + (i + 0.5) * dx, y0 + (j + 0.5) * dy))
    gmean, emean = sum(ground) / len(ground), sum(exact) / len(exact)
    diff = [(g - gmean) - (e - emean) for g, e in zip(ground, exact)]
    relief = max(exact) - min(exact)
    top_g = max(range(len(ground)), key=lambda n: ground[n] - HILL_SLOPE * (x0 + (n // ny + 0.5) * dx))
    top_e = max(range(len(exact)), key=lambda n: exact[n] - HILL_SLOPE * (x0 + (n // ny + 0.5) * dx))
    shift = math.hypot((top_g // ny - top_e // ny) * dx, (top_g % ny - top_e % ny) * dy)
    check("hill: crest position error (m)", shift, 0.0, 0.0, "abs")
    check("hill: relief (m)", max(ground) - min(ground), relief, 2e-2, "rel")
    check("hill: RMS of ground - formula / relief",
          math.sqrt(sum(d * d for d in diff) / len(diff)) / relief, 0.0, 3e-3, "abs")


def main():
    argv = sys.argv[1:]
    extra = []
    if "--" in argv:
        n = argv.index("--")
        argv, extra = argv[:n], argv[n + 1:]
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("case", choices=("flat_sweep", "flat_nest", "column_target", "column_axis", "hill", "upstream"))
    ap.add_argument("--gust", required=True, help="the GUST executable")
    ap.add_argument("--mpi", default="", help="the launcher run.sh puts before GUST (empty: none)")
    ap.add_argument("--precision", choices=("DOUBLE", "SINGLE"), default="DOUBLE", help="GUST's precision "
                    "(GUST_PRECISION): the library's direction-invariance tolerances")
    opts = ap.parse_args(argv)

    yaml_file = opts.case + ".yaml"
    shutil.copyfile(os.path.join(HERE, yaml_file), yaml_file)
    if opts.case == "hill":
        write_hill("hill.asc")
    if opts.case == "upstream":
        write_regional("regional.xyz")
    if os.path.isdir(opts.case):
        shutil.rmtree(opts.case)
    tool = os.path.join(ROOT, "tools", "make_gust_case.py")
    made = subprocess.run([sys.executable, tool, yaml_file, "--quiet", "--strict"])   # exit 1 on a warning too
    check("make_gust_case.py --strict exit status (0: written, no error or warning)", float(made.returncode), 0.0, 0.0, "abs")
    if made.returncode not in (0, 1):
        return 1
    env = dict(os.environ, GUST=os.path.abspath(opts.gust), MPI=opts.mpi)
    ran = subprocess.run(["bash", os.path.join(opts.case, "run.sh")] + extra, env=env)
    check("run.sh exit status", float(ran.returncode), 0.0, 0.0, "abs")
    if ran.returncode != 0:
        return 1

    if opts.case == "flat_sweep":
        for d in (240.0, 0.0):
            check_flat(os.path.join(opts.case, "dir_%03d" % d), d, "from %g" % d)
    elif opts.case == "flat_nest":
        dx_parent = check_flat(opts.case, 240.0, "parent")
        dx_patch = check_flat(os.path.join(opts.case, "nest_centre"), 240.0, "patch")
        check("patch cell size / parent cell size (x)", dx_patch[0] / dx_parent[0], 0.5, 1e-12, "abs")
        check("patch cell size / parent cell size (y)", dx_patch[1] / dx_parent[1], 0.5, 1e-12, "abs")
    elif opts.case == "column_target":
        top = check_column_table(os.path.join(opts.case, "column", "column_guess.dat"), 240.0)
        check_column_handover(os.path.join(opts.case, "log.txt"), top)
        speed, direction, _ = wind_at_height(opts.case)
        check("3-D run: speed at %g m (m/s)" % HEIGHT, speed, SPEED, 1e-2, "rel")
        check("3-D run: direction error (degrees)", angle_error(direction, towards(240.0)), 0.0, 0.3, "abs")
        check_library(opts.case, (240.0,), opts.precision == "SINGLE")
    elif opts.case == "column_axis":
        for d in (0.0, 270.0):
            check_column_axis(os.path.join(opts.case, "dir_%03d" % d), "from %g" % d)
        check_library(opts.case, (0.0, 270.0, 240.0), opts.precision == "SINGLE")
    elif opts.case == "upstream":
        check_upstream(opts.case)
    else:
        check_hill(opts.case, 1500.0, 4.0)

    if FAILED:
        print("check_case_setup: %d check(s) failed: %s" % (len(FAILED), "; ".join(FAILED)))
        return 1
    print("check_case_setup: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
