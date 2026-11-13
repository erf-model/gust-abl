#!/usr/bin/env python3
"""RT-6.2.1: the 2-D hill as unsteady RANS on 1, 2 and 4 ranks.

    check_unsteady_hill2d.py -- LAUNCHER... GUST [ARGS ...]

LAUNCHER... GUST is the application with its launcher, ending with the
application (the script sets the rank count: it replaces the number after
"-n" or "-np" for the 1- and 4-rank runs). ARGS are the deck
(Tests/Canonical_RANS/Unsteady_Hill_2D/inputs_unsteady_hill2d) and the flags.
The script runs the deck from the log-law first guess in np2/, np1/ and np4/
(each emptied first) and checks:
  - the 2-rank run logs ten physical steps of 60 s, ending at t = 600 s, every
    step converged (none stopped at unsteady.max_inner), and the summary line;
  - the 1- and 4-rank runs log the same steps with the same inner iteration
    counts, and their final u, v, w within 1e-9 of the largest |u| of the
    2-rank run, p, theta and k within 1e-9 of their own largest values;
  - one step of 1e-3 s (in short/) converges within 20 inner iterations
    (measured 12): the physical time term (1/dt = 1000 1/s on every diagonal)
    reaches the stage through the application's level solver. Without it
    (res_in.time not set) this step and every 60 s step stop at
    unsteady.max_inner = 60 with their residual still above 0.93 of its start;
  - the 2-rank plotfile's header carries the physical time, 600 s;
  - the flow has moved: u at the hill top 50 m above the ground differs from
    the short run's (nearly the first guess) by more than 0.1 m/s.
The exit code is the verdict.
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

FIELDS = ["u", "v", "w", "p", "theta", "tke"]
STEP = re.compile(r"^gust: step (\d+): t = (\S+) s, dt = (\S+) s, (\d+) inner iterations, inner residual (\S+)(.*)$")


def with_ranks(cmd, n):
    out = list(cmd)
    for i, a in enumerate(out[:-1]):
        if a in ("-n", "-np"):
            out[i + 1] = str(n)
            return out
    if n != 2:
        raise SystemExit("check_unsteady_hill2d: no -n in the launcher, cannot run on %d ranks" % n)
    return out


def run(cmd, where, extra):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    full = cmd + extra
    print("check_unsteady_hill2d: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_unsteady_hill2d: the command failed with code %d in %s/" % (out.returncode, where))
    return out.stdout


def steps(stdout):
    return [STEP.match(ln) for ln in stdout.splitlines() if STEP.match(ln)]


def last_plotfile(stdout, where):
    wrote = [ln.split()[-1] for ln in stdout.splitlines()
             if ln.startswith("gust: wrote ") and "_plt" in ln and "_level" not in ln]
    if not wrote:
        raise SystemExit("check_unsteady_hill2d: no plotfile written in %s/" % where)
    return os.path.join(where, wrote[-1])


def plotfile_time(path):
    # AMReX plotfile Header: version, the number of variables, their names, the dimension, then the time.
    with open(os.path.join(path, "Header")) as fh:
        lines = fh.read().splitlines()
    nvar = int(lines[1])
    return float(lines[2 + nvar + 1])


def values(a, f):
    return [x for plane in a[f] for col in plane for x in col]


def main(argv):
    if len(argv) < 3 or argv[1] != "--":
        raise SystemExit(__doc__)
    cmd = argv[2:]
    fresh = ["init.first_guess=profile", "io.write_checkpoint=0"]
    ten = fresh + ["unsteady.stop_time=600", "unsteady.dt=60"]
    ok = True

    def check(good, what):
        nonlocal ok
        ok = ok and good
        print("  %s: %s" % (what, "pass" if good else "FAIL"))

    two = run(with_ranks(cmd, 2), "np2", ten)
    s2 = steps(two)
    check(len(s2) == 10, "2 ranks: 10 physical steps (got %d)" % len(s2))
    check(bool(s2) and all(abs(float(m.group(3)) - 60.0) < 1e-9 for m in s2), "every dt is 60 s")
    check(bool(s2) and abs(float(s2[-1].group(2)) - 600.0) < 1e-9, "the last step ends at t = 600 s")
    check(bool(s2) and all(m.group(6) == "" for m in s2), "every step converged (none at unsteady.max_inner, none diverged)")
    check("gust: unsteady: 10 physical steps to t = 600 s" in two and "0 steps stopped at unsteady.max_inner" in two,
          "the summary line")
    print("  inner iterations per step: %s" % " ".join(m.group(4) for m in s2))
    plt2 = last_plotfile(two, "np2")
    tp = plotfile_time(plt2)
    check(abs(tp - 600.0) < 1e-9, "the plotfile's time is %g s (600)" % tp)
    _, a = ep.read_fields(plt2, FIELDS + ["z"])
    umax = max(abs(x) for x in values(a, "u"))

    for n in (1, 4):
        out = run(with_ranks(cmd, n), "np%d" % n, ten)
        sn = steps(out)
        check(bool(sn) and [m.group(4) for m in sn] == [m.group(4) for m in s2],
              "%d rank(s): the same inner iteration counts as 2 ranks" % n)
        _, b = ep.read_fields(last_plotfile(out, "np%d" % n), FIELDS)
        worst = 0.0
        for f in FIELDS:
            va, vb = values(a, f), values(b, f)
            big = umax if f in ("u", "v", "w") else max(abs(x) for x in va)
            worst = max(worst, max(abs(x - y) for x, y in zip(va, vb)) / big)
        check(worst <= 1e-9, "%d rank(s): max relative difference to 2 ranks %.2e (<= 1e-9)" % (n, worst))

    short = run(with_ranks(cmd, 2), "short", fresh + ["unsteady.stop_time=1e-3", "unsteady.dt=1e-3"])
    ss = steps(short)
    check(len(ss) == 1 and ss[0].group(6) == "" and int(ss[0].group(4)) <= 20,
          "one step of 1e-3 s converged in %s inner iterations (<= 20)" % (ss[0].group(4) if ss else "?"))
    _, c = ep.read_fields(last_plotfile(short, "short"), ["u", "z"])
    # The hill top is at x = 0, the column nearest the domain's centre (128 columns over -2560..2560 m).
    ii = 63
    zg = a["z"][ii][0][0] - 15.625 / 2
    k = min(range(len(a["z"][ii][0])), key=lambda kk: abs(a["z"][ii][0][kk] - zg - 50.0))
    du = abs(a["u"][ii][0][k] - c["u"][ii][0][k])
    check(du > 0.1, "u at the hill top 50 m above the ground moved by %.3f m/s (> 0.1) over 600 s" % du)

    print("check_unsteady_hill2d: %s" % ("pass" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
