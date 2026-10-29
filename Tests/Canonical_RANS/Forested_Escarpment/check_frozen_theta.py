#!/usr/bin/env python3
"""The forested escarpment with a prescribed stable theta (simplec.theta = frozen).

    check_frozen_theta.py --max-iterations N -- LAUNCHER... GUST [ARGS ...]

LAUNCHER... GUST is the application with its launcher, ending with the
application (the script replaces the number after "-n" or "-np" for the
4-rank run). ARGS are the deck (inputs_frozen_theta: the forested escarpment with a stable
sounding, 300 K to 1000 m, an 8 K inversion to 1100 m, 3 K/km above; its
values are SOUNDING_Z, SOUNDING_T here) and the flags. The script writes the
tables with make_case.py in each run directory (emptied first) and runs:

  frozen/   2 ranks, simplec.theta = frozen, to the verdict;
  np4/      the same on 4 ranks;
  file/     2 ranks, frozen, 3 iterations, with inputs_frozen_theta_file (next to
            the deck): the inflow from a table this script writes, whose theta
            (TABLE_Z, TABLE_T) differs from the deck's sounding, and
            init.z_agl_min = 5;
  solve/    2 ranks, theta solved (the default), 5 iterations.

and checks:
  - the frozen run is steady within N iterations, with a theta residual of 0
    in every row of the iteration table;
  - theta is the sounding at each cell centre's height above the ground,
    h = (mean of the cell's eight node heights) - (mean of its four ground
    nodes), the ground nodes the tanh terrain at the node x as make_case.py
    writes it (level 0) or the level plotfile's z_phys_nd (the coarse levels,
    whose ground is smoothed) and the cell-centre heights the plotfile's z: at
    every cell of
    the final level-0 plotfile within 1e-9 K, and on every sequenced level
    the logged minimum, maximum and mean of the frozen theta equal the
    sounding's over that level's cells (its level plotfile's z) within 1e-9 K;
  - the turbulence plotfile's R_theta is 0 (theta is not an unknown);
  - with the table, theta is the table's theta at h (not floored at
    init.z_agl_min) on level 0 and by the logged min, max and mean on every
    level, within 1e-9 K, and every first cell has Rt > 0 (the table is stable
    from the ground; theta's ground ghosts extend the first two cells, where
    the neutral wall's mirror would give Rt = 0);
  - the 4-rank run reaches the same answer (1e-9 of the largest |u|, p and k
    of their own largest values);
  - the solved run's theta residual is not 0 (the old path, theta moves).
The exit code is the verdict.
"""
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

SOUNDING_Z = [0.0, 1000.0, 1100.0, 2400.0]
SOUNDING_T = [300.0, 300.0, 308.0, 308.0 + 0.003 * 1300.0]
H, L, XC = 200.0, 250.0, 3000.0   # make_case.py's terrain
XLEN = 6000.0                      # the deck's domain length in x (from x = 0)


# The table of the file-profile deck (inputs_frozen_theta_file): its own theta, 299 K + 5 K/km to 900 m, 6 K more
# by 1000 m, then 2 K/km; the wind the log law with u* 0.35 m/s and z0 0.1 m; k 0.3.
TABLE_Z = [0.0, 900.0, 1000.0, 2400.0]
TABLE_T = [299.0, 303.5, 309.5, 312.3]
TABLE_DZ = 10.0
Z_AGL_MIN = 5.0   # given to the file-profile run: the frozen theta must ignore it


def interp(h, zs, ts):
    if h <= zs[0]:
        return ts[0]
    for a in range(len(zs) - 1):
        if h <= zs[a + 1]:
            f = (h - zs[a]) / (zs[a + 1] - zs[a])
            return ts[a] + f * (ts[a + 1] - ts[a])
    return ts[-1]


def sounding(h):
    return interp(h, SOUNDING_Z, SOUNDING_T)


def table_theta(h):
    return interp(h, TABLE_Z, TABLE_T)


def write_table(where):
    """frozen_inflow.txt: z u v theta k every TABLE_DZ m (the breakpoints among them, so the table's own linear
    interpolation is table_theta)."""
    n = int(round(TABLE_Z[-1] / TABLE_DZ))
    with open(os.path.join(where, "frozen_inflow.txt"), "w") as f:
        for i in range(n + 1):
            z = i * TABLE_DZ
            u = 0.35 / 0.41 * math.log((z + 0.1) / 0.1)
            f.write("%.1f %.6f 0.0 %.12f 0.3\n" % (z, u, table_theta(z)))


def ground(x):
    """The terrain at a node, as make_case.py writes it (rounded to 1e-6 m)."""
    return float("%.6f" % (0.5 * H * (1 + math.tanh((x - XC) / L))))


def with_ranks(cmd, n):
    out = list(cmd)
    for i, a in enumerate(out[:-1]):
        if a in ("-n", "-np"):
            out[i + 1] = str(n)
            return out
    if n != 2:
        raise SystemExit("check_frozen_theta: no -n in the launcher, cannot run on %d ranks" % n)
    return out


def run(cmd, where, extra, table=False):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    subprocess.run([sys.executable, os.path.join(HERE, "make_case.py"), where], check=True)
    if table:
        write_table(where)
    full = cmd + extra
    print("check_frozen_theta: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_frozen_theta: the command failed with code %d in %s/" % (out.returncode, where))
    return out.stdout


def verdict(stdout):
    for ln in stdout.splitlines():
        if ln.startswith("gust: steady after "):
            return int(ln.split()[3])
    return None


ROW = re.compile(r"^(?:level \d+:)?\s*\d+\s+\S+e[+-]\d+\s+[a-z_]+\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$")


def theta_residuals(stdout):
    """The theta-res column of every row of the iteration tables (every level): u, v, w, tke, theta, p, r."""
    return [float(m.group(5)) for m in map(ROW.match, stdout.splitlines()) if m]


def plotfiles(stdout, where):
    """Every plotfile the run wrote (the coarse levels' and the final one), last of each prefix."""
    wrote = [ln.split()[-1] for ln in stdout.splitlines() if ln.startswith("gust: wrote ") and "_plt" in ln]
    last = {}
    for w in wrote:
        last[re.sub(r"\d+$", "", w)] = w
    if not last:
        raise SystemExit("check_frozen_theta: no plotfile written in %s/" % where)
    return [os.path.join(where, w) for w in last.values()]


def values(a, f):
    return [x for plane in a[f] for col in plane for x in col]


def level_ground(plt):
    """{(i, j): node height} of the ground nodes (k = 0) of a level plotfile's z_phys_nd (MeshHierarchy writes it next
    to Cell_H, in the same VisMF layout; the coarse levels' ground is level 0's injected and smoothed)."""
    level_dir = os.path.join(plt, "Level_0")
    with open(os.path.join(level_dir, "z_phys_nd_H")) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    nbox = int(re.match(r"\((\d+)\s+(\d+)", lines[4]).group(1))
    boxes = []
    for b in range(nbox):
        lo, hi = re.findall(r"\(([-\d,]+)\)", lines[5 + b])[:2]
        boxes.append(([int(v) for v in lo.split(",")], [int(v) for v in hi.split(",")]))
    idx = 5 + nbox + 1
    nfab = int(lines[idx])
    fabs = [re.match(r"FabOnDisk:\s+(\S+)\s+(\d+)", lines[idx + 1 + f]).groups() for f in range(nfab)]
    out = {}
    for (blo, bhi), (fname, off) in zip(boxes, fabs):
        flo, fhi, data = ep._read_fab(os.path.join(level_dir, fname), int(off), 1)
        fnx, fny = fhi[0] - flo[0] + 1, fhi[1] - flo[1] + 1
        k = 0
        for jj in range(blo[1], bhi[1] + 1):
            for ii in range(blo[0], bhi[0] + 1):
                out[(ii, jj)] = data[((k - flo[2]) * fny + (jj - flo[1])) * fnx + (ii - flo[0])]
    return out


def sounding_at_cells(plt, theta_of=sounding):
    """theta(h) at every cell of a level: h from the plotfile's cell-centre z and the cell's four ground nodes (the
    level plotfile's z_phys_nd on the coarse levels, make_case.py's terrain on level 0)."""
    _, a = ep.read_fields(plt, ["z"])
    nx = len(a["z"])
    zg_node = level_ground(plt) if "_level" in plt else None
    out = []
    for i in range(nx):
        plane = []
        for j, col in enumerate(a["z"][i]):
            if zg_node is None:
                dx = XLEN / nx
                zg = 0.5 * (ground(i * dx) + ground((i + 1) * dx))
            else:
                zg = 0.25 * (zg_node[(i, j)] + zg_node[(i + 1, j)] + zg_node[(i, j + 1)] + zg_node[(i + 1, j + 1)])
            plane.append([theta_of(z - zg) for z in col])
        out.append(plane)
    return out


LOGGED = re.compile(r"^gust: level (\d+): theta frozen at the inflow profile's theta by height above the ground: "
                    r"min (\S+), max (\S+), mean (\S+) K$", re.M)


def main(argv):
    nmax = None
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] == "--max-iterations" and i + 1 < len(argv):
            nmax = int(argv[i + 1])
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if nmax is None or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    fresh = ["init.first_guess=profile"]
    ok = True

    frozen = run(with_ranks(cmd, 2), "frozen", fresh + ["simplec.theta=frozen"])
    n2 = verdict(frozen)
    good = n2 is not None and n2 <= nmax
    ok = ok and good
    print("  frozen, 2 ranks: steady after %s iterations (<= %d): %s" % (n2, nmax, "pass" if good else "FAIL"))
    res = theta_residuals(frozen)
    good = len(res) > 0 and max(res) == 0
    ok = ok and good
    print("  frozen: theta residual 0 in all %d table rows: %s" % (len(res), "pass" if good else "FAIL"))
    plts = plotfiles(frozen, "frozen")
    logged = {int(m.group(1)): (float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in LOGGED.finditer(frozen)}
    for plt in plts:
        m = re.search(r"_level(\d+)_", plt)
        lev = int(m.group(1)) if m else 0
        exact = sounding_at_cells(plt)
        flat = [x for plane in exact for col in plane for x in col]
        want = (min(flat), max(flat), sum(flat) / len(flat))
        got = logged.get(lev)
        good = got is not None and all(abs(g - w) <= 1e-9 for g, w in zip(got, want))
        ok = ok and good
        print("  frozen: level %d logged theta min, max, mean %s; sounding(h) over its %d cells %.9f, %.9f, %.9f K "
              "(1e-9): %s" % (lev, "%.9f, %.9f, %.9f" % got if got else "none", len(flat), *want, "pass" if good else "FAIL"))
        if lev == 0:
            _, a = ep.read_fields(plt, ["theta"])
            err = max(abs(t - e) for pa, pe in zip(a["theta"], exact) for ca, ce in zip(pa, pe) for t, e in zip(ca, ce))
            good = err <= 1e-9
            ok = ok and good
            print("  frozen: level 0 final theta, max |theta - sounding(h)| = %.2e K (<= 1e-9): %s"
                  % (err, "pass" if good else "FAIL"))
    final = [p for p in plts if "_level" not in p]
    _, a = ep.read_fields(final[0], ["u", "v", "w", "p", "tke"])

    four = run(with_ranks(cmd, 4), "np4", fresh + ["simplec.theta=frozen"])
    _, b = ep.read_fields([p for p in plotfiles(four, "np4") if "_level" not in p][0], ["u", "v", "w", "p", "tke"])
    umax = max(abs(x) for x in values(a, "u"))
    worst = 0.0
    for f in ("u", "v", "w", "p", "tke"):
        va, vb = values(a, f), values(b, f)
        big = umax if f in ("u", "v", "w") else max(abs(x) for x in va)
        worst = max(worst, max(abs(x - y) for x, y in zip(va, vb)) / big)
    good = verdict(four) is not None and worst <= 1e-9
    ok = ok and good
    print("  frozen, 4 ranks: steady after %s iterations, max relative difference to 2 ranks %.2e (<= 1e-9): %s"
          % (verdict(four), worst, "pass" if good else "FAIL"))

    turb = [ln.split()[-1] for ln in frozen.splitlines() if ln.startswith("gust: wrote ") and "_turb" in ln]
    if turb:
        _, t = ep.read_fields(os.path.join("frozen", turb[-1]), ["R_theta"])
        rmax = max(abs(x) for x in values(t, "R_theta"))
    good = bool(turb) and rmax == 0
    ok = ok and good
    print("  frozen: turbulence plotfile max |R_theta| = %s (0): %s"
          % ("%.3e" % rmax if turb else "no plotfile", "pass" if good else "FAIL"))

    # The WINSENT path: theta from an inflow table (not the deck's sounding), with init.z_agl_min given.
    file_cmd = list(cmd)
    file_cmd[-1 - next(k for k, a in enumerate(reversed(file_cmd)) if a.endswith("inputs_frozen_theta"))] = \
        os.path.join(HERE, "inputs_frozen_theta_file")
    table = run(with_ranks(file_cmd, 2), "file", fresh + ["simplec.theta=frozen", "steady.max_iterations=3",
                                                          "init.z_agl_min=%g" % Z_AGL_MIN], table=True)
    logged = {int(m.group(1)): (float(m.group(2)), float(m.group(3)), float(m.group(4))) for m in LOGGED.finditer(table)}
    for plt in plotfiles(table, "file"):
        m = re.search(r"_level(\d+)_", plt)
        lev = int(m.group(1)) if m else 0
        exact = sounding_at_cells(plt, table_theta)
        flat = [x for plane in exact for col in plane for x in col]
        want = (min(flat), max(flat), sum(flat) / len(flat))
        got = logged.get(lev)
        good = got is not None and all(abs(g - w) <= 1e-9 for g, w in zip(got, want))
        ok = ok and good
        print("  file profile: level %d logged theta min, max, mean %s; table(h) %.9f, %.9f, %.9f K (1e-9): %s"
              % (lev, "%.9f, %.9f, %.9f" % got if got else "none", *want, "pass" if good else "FAIL"))
        if lev == 0:
            _, a0 = ep.read_fields(plt, ["theta"])
            err = max(abs(t - e) for pa, pe in zip(a0["theta"], exact) for ca, ce in zip(pa, pe) for t, e in zip(ca, ce))
            good = err <= 1e-9
            ok = ok and good
            print("  file profile: level 0 theta, max |theta - table(h)| = %.2e K (<= 1e-9, init.z_agl_min = %g "
                  "not applied): %s" % (err, Z_AGL_MIN, "pass" if good else "FAIL"))

    # The first cell sees the table's stable gradient (theta's ground ghosts extend the first two cells; the
    # neutral wall's mirror would give d theta/dz = 0 there, so Rt = 0): Rt > 0 in every first cell.
    tf = [ln.split()[-1] for ln in table.splitlines() if ln.startswith("gust: wrote ") and "_turb" in ln]
    if tf:
        _, rt = ep.read_fields(os.path.join("file", tf[-1]), ["Rt"])
        first = [col[0] for plane in rt["Rt"] for col in plane]
    good = bool(tf) and min(first) > 0
    ok = ok and good
    print("  file profile: first-cell Rt min %s over %d columns (> 0): %s"
          % ("%.3e" % min(first) if tf else "no plotfile", len(first) if tf else 0, "pass" if good else "FAIL"))

    solve = run(with_ranks(cmd, 2), "solve", fresh + ["steady.max_iterations=5"])
    res = theta_residuals(solve)
    good = len(res) > 0 and min(res) > 0
    ok = ok and good
    print("  solved theta: theta residual > 0 in all %d table rows: %s" % (len(res), "pass" if good else "FAIL"))
    if not ok:
        raise SystemExit("check_frozen_theta: FAIL")
    print("check_frozen_theta: pass")


if __name__ == "__main__":
    main(sys.argv)
