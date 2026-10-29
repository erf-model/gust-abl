#!/usr/bin/env python3
"""The forested escarpment: a dense canopy on steep terrain with an inflow face.

    check_forested_escarpment.py --max-iterations N -- LAUNCHER... GUST [ARGS ...]

LAUNCHER... GUST is the application with its launcher, ending with the
application (the script sets the rank count: it replaces the number after
"-n" or "-np" for the 1- and 4-rank runs). ARGS are the deck
(Tests/Canonical_RANS/Forested_Escarpment/inputs_forested_escarpment) and the
flags. The script writes the tables with make_case.py in each run directory
(emptied first) and runs:

  np2/      2 ranks, to the verdict;
  np1/      1 rank, to the verdict;
  np4/      4 ranks, to the verdict;
  first/    2 ranks, 40 iterations; then continued from its checkpoint to 80;
  straight/ 2 ranks, 80 iterations;
  inflow/   the forest extended to the inflow face, 1 iteration.

and checks:
  - the 2-rank run is steady within N iterations at the default CFL cap (8),
    and logs no inflow-face or damping-layer warning;
  - the 1- and 4-rank runs reach the same answer: u, v and w within 1e-9 of
    the largest |u|, p and k within 1e-9 of their own largest values (v is
    rounding noise in this 2-D case, so it is measured against |u|);
  - the canopy acts: 12 m above the ground, u at the forest's centre on the
    slope (x = 2900 m) is below half of u on open ground upstream (x = 1500 m);
  - the continued and straight plotfiles are identical in every field;
  - with the forest at the inflow face, the inflow-face warning names a
    positive count of canopy cells in the first column.
The exit code is the verdict.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

FIELDS = ["u", "v", "w", "p", "theta", "tke"]


def with_ranks(cmd, n):
    out = list(cmd)
    for i, a in enumerate(out[:-1]):
        if a in ("-n", "-np"):
            out[i + 1] = str(n)
            return out
    if n != 2:
        raise SystemExit("check_forested_escarpment: no -n in the launcher, cannot run on %d ranks" % n)
    return out


def run(cmd, where, extra, forest_at_inflow=False):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    make = [sys.executable, os.path.join(HERE, "make_case.py"), where] + (["--forest-at-inflow"] if forest_at_inflow else [])
    subprocess.run(make, check=True)
    full = cmd + extra
    print("check_forested_escarpment: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_forested_escarpment: the command failed with code %d in %s/" % (out.returncode, where))
    return out.stdout


def verdict(stdout):
    for ln in stdout.splitlines():
        if ln.startswith("gust: steady after "):
            return int(ln.split()[3])
    return None


def last_plotfile(stdout, where):
    wrote = [ln.split()[-1] for ln in stdout.splitlines()
             if ln.startswith("gust: wrote ") and "_plt" in ln and "_level" not in ln]
    if not wrote:
        raise SystemExit("check_forested_escarpment: no plotfile written in %s/" % where)
    return os.path.join(where, wrote[-1])


def values(a, f):
    return [x for plane in a[f] for col in plane for x in col]


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

    two = run(with_ranks(cmd, 2), "np2", fresh)
    n2 = verdict(two)
    good = n2 is not None and n2 <= nmax
    ok = ok and good
    print("  2 ranks: steady after %s iterations (<= %d): %s" % (n2, nmax, "pass" if good else "FAIL"))
    if n2 is None:
        raise SystemExit("check_forested_escarpment: FAIL")
    good = "canopy cells lie" not in two
    ok = ok and good
    print("  no inflow-face or damping-layer warning: %s" % ("pass" if good else "FAIL"))
    _, a = ep.read_fields(last_plotfile(two, "np2"), FIELDS + ["z"])

    for n in (1, 4):
        out = run(with_ranks(cmd, n), "np%d" % n, fresh)
        good = verdict(out) is not None
        ok = ok and good
        _, b = ep.read_fields(last_plotfile(out, "np%d" % n), FIELDS)
        worst = 0.0
        umax = max(abs(x) for x in values(a, "u"))
        for f in ("u", "v", "w", "p", "tke"):
            va, vb = values(a, f), values(b, f)
            big = umax if f in ("u", "v", "w") else max(abs(x) for x in va)
            worst = max(worst, max(abs(x - y) for x, y in zip(va, vb)) / big)
        good = good and worst <= 1e-9
        ok = ok and good
        print("  %d rank(s): steady after %s iterations, max relative difference to 2 ranks %.2e (<= 1e-9): %s"
              % (n, verdict(out), worst, "pass" if good else "FAIL"))

    # 12 m above the ground: the layer whose centre is nearest 12 m (8 m first layer: the second, 12.5 m). The ground
    # is taken as the first cell centre minus 4 m, which is approximate on the slope (the first layer thins with the
    # terrain by a few tenths of a metre) and does not change the layer chosen.
    def u_at(x):
        ii = int(round(x / 25.0 - 0.5))   # the column whose centre is nearest x
        zg = a["z"][ii][3][0] - 4.0
        k = min(range(len(a["z"][ii][3])), key=lambda kk: abs(a["z"][ii][3][kk] - zg - 12.0))
        return a["u"][ii][3][k]
    u_in, u_open = u_at(2900.0), u_at(1500.0)
    good = u_in < 0.5 * u_open
    ok = ok and good
    print("  canopy: u 12 m above the ground %.3f m/s in the forest, %.3f m/s upstream (ratio %.2f < 0.5): %s"
          % (u_in, u_open, u_in / u_open, "pass" if good else "FAIL"))

    first = run(with_ranks(cmd, 2), "first", fresh + ["steady.max_iterations=40"])
    chk = [ln.split()[-1] for ln in first.splitlines() if ln.startswith("gust: wrote ") and "_chk" in ln]
    if not chk:
        raise SystemExit("check_forested_escarpment: no checkpoint written in first/")
    cont = subprocess.run(with_ranks(cmd, 2) + ["steady.max_iterations=80", "io.restart=" + chk[-1]], cwd="first",
                          check=False, capture_output=True, text=True)
    if cont.returncode != 0:
        sys.stdout.write(cont.stdout[-4000:] + cont.stderr[-4000:])
        raise SystemExit("check_forested_escarpment: the restart failed with code %d" % cont.returncode)
    straight = run(with_ranks(cmd, 2), "straight", fresh + ["steady.max_iterations=80"])
    _, c = ep.read_fields(last_plotfile(cont.stdout, "first"), FIELDS)
    _, s = ep.read_fields(last_plotfile(straight, "straight"), FIELDS)
    for f in FIELDS:
        d = max(abs(x - y) for x, y in zip(values(c, f), values(s, f)))
        good = d == 0
        ok = ok and good
        print("  restart: max |%s(continued) - %s(straight)| = %.3e (bit for bit): %s"
              % (f, f, d, "pass" if good else "FAIL"))

    inflow = run(with_ranks(cmd, 2), "inflow", fresh + ["steady.max_iterations=1"], forest_at_inflow=True)
    line = [ln for ln in inflow.splitlines() if "canopy cells lie in the first cell column of an inflow face" in ln]
    count = int(line[0].split()[2]) if line else 0
    good = count > 0
    ok = ok and good
    print("  forest at the inflow face: warning with %d canopy cells: %s" % (count, "pass" if good else "FAIL"))
    if not ok:
        raise SystemExit("check_forested_escarpment: FAIL")
    print("check_forested_escarpment: pass")


if __name__ == "__main__":
    main(sys.argv)
