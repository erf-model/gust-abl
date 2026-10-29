#!/usr/bin/env python3
"""The column first guess with a uniform canopy (UT-8.5.1, ColumnGuess.H):

    check_canopy_column.py --tol T --max-ratio R -- LAUNCHER... GUST [ARGS ...]

ARGS are the deck (Tests/Canonical_RANS/Canopy_Drag/inputs_canopy_drag) and the
flags, without canopy.type (the script gives the map canopy). The script
writes a uniform map canopy (20 m, LAI 3, Cd 0.25 on every column, the
tables in the run folder) and runs:

  column/   init.first_guess = column, with io.column_guess_file;
  profile/  init.first_guess = profile;
  patches/  init.first_guess = column with the two-patch canopy of the case
            (canopy.file, a canopy that differs between columns), 1 iteration.

and checks:
  - the column run's first u and v residuals are below T times the profile
    run's (the column is the 3-D steady column: the drag, the k terms and the
    length cap are the stage's on a flat column);
  - the converged 3-D u and v (the column mean over the periodic domain) are
    within T of the column table at every level, relative to the table's
    largest speed;
  - the column run reaches the verdict in at most R times the profile run's
    iterations;
  - the two-patch canopy is refused, naming init.first_guess.
The exit code is the verdict.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

N, L = 16, 1280.0   # the deck's columns


def write_map(where, h, lai):
    xs = [-2 * L / N + i * L / N for i in range(N + 5)]
    with open(os.path.join(where, "h.xyz"), "w") as fh, open(os.path.join(where, "lai.xyz"), "w") as fl:
        for y in xs:
            for x in xs:
                fh.write("%.6f %.6f %.6f\n" % (x, y, h))
                fl.write("%.6f %.6f %.6f\n" % (x, y, lai))


def run(cmd, where, extra, expect_fail=False):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    write_map(where, 20.0, 3.0)
    full = cmd + extra
    print("check_canopy_column: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if (out.returncode != 0) != expect_fail:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_canopy_column: unexpected exit code %d in %s/" % (out.returncode, where))
    return out.stdout + out.stderr


def verdict(stdout):
    for ln in stdout.splitlines():
        if ln.startswith("gust: steady after "):
            return int(ln.split()[3])
    return None


def first_residuals(stdout):
    """(u, v) of the first level-0 table row."""
    rows = [ln for ln in stdout.splitlines() if ln.strip() and ln.split()[0] == "1" and "lusgs" in ln]
    cols = rows[0].split()
    return float(cols[3]), float(cols[4])


def last_plotfile(stdout, where):
    wrote = [ln.split()[-1] for ln in stdout.splitlines()
             if ln.startswith("gust: wrote ") and "_plt" in ln and "_level" not in ln]
    return os.path.join(where, wrote[-1])


def main(argv):
    tol, ratio = None, None
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] == "--tol":
            tol = float(argv[i + 1]); i += 2
        elif argv[i] == "--max-ratio":
            ratio = float(argv[i + 1]); i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if tol is None or ratio is None or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    canopy = ["canopy.type=map", "canopy.height_file=h.xyz", "canopy.lai_file=lai.xyz", "canopy.cd=0.25"]
    ok = True

    col = run(cmd, "column", canopy + ["init.first_guess=column", "io.column_guess_file=column.txt"])
    pro = run(cmd, "profile", canopy + ["init.first_guess=profile"])
    ncol, npro = verdict(col), verdict(pro)
    good = ncol is not None and npro is not None
    ok = ok and good
    print("  steady: column after %s, profile after %s iterations: %s" % (ncol, npro, "pass" if good else "FAIL"))
    if not good:
        raise SystemExit("check_canopy_column: FAIL")
    (uc, vc), (up, vp) = first_residuals(col), first_residuals(pro)
    good = uc <= tol * up and vc <= tol * vp
    ok = ok and good
    print("  first residuals, column / profile: u %.2e / %.2e, v %.2e / %.2e (<= %g): %s"
          % (uc, up, vc, vp, tol, "pass" if good else "FAIL"))
    good = ncol <= ratio * npro
    ok = ok and good
    print("  iterations: column %d <= %g x profile %d: %s" % (ncol, ratio, npro, "pass" if good else "FAIL"))

    # The converged 3-D column against the table.
    table = [[float(x) for x in ln.split()] for ln in open(os.path.join("column", "column.txt")) if not ln.startswith("#")]
    _, a = ep.read_fields(last_plotfile(col, "column"), ["u", "v", "z"])
    nx, ny, nz = len(a["u"]), len(a["u"][0]), len(a["u"][0][0])
    umax = max(abs(r[1]) for r in table) + max(abs(r[2]) for r in table)
    worst = 0.0
    for k in range(nz):
        um = sum(a["u"][i][j][k] for i in range(nx) for j in range(ny)) / (nx * ny)
        vm = sum(a["v"][i][j][k] for i in range(nx) for j in range(ny)) / (nx * ny)
        zc = a["z"][0][0][k]
        row = min(table[1:], key=lambda r: abs(r[0] - zc))
        worst = max(worst, abs(um - row[1]) / umax, abs(vm - row[2]) / umax)
    good = worst <= tol
    ok = ok and good
    print("  converged 3-D u, v against the column table: largest difference %.2e of the table's speed (<= %g): %s"
          % (worst, tol, "pass" if good else "FAIL"))

    bad = run(cmd, "patches", ["canopy.type=patches", "canopy.file=" + os.path.join(HERE, "canopy_patches.txt"),
                               "init.first_guess=column", "steady.max_iterations=1"], expect_fail=True)
    good = "init.first_guess = column needs a horizontally uniform canopy" in bad
    ok = ok and good
    print("  a canopy that differs between columns is refused: %s" % ("pass" if good else "FAIL"))
    print("check_canopy_column: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
