#!/usr/bin/env python3
"""RT-6.2.1 and RT-6.3.1: the 2-D hill as unsteady RANS on 1, 2 and 4 ranks,
its outputs and a restart.

    check_unsteady_hill2d.py -- LAUNCHER... GUST [ARGS ...]

LAUNCHER... GUST is the application with its launcher, ending with the
application (the script sets the rank count: it replaces the number after
"-n" or "-np" for the 1- and 4-rank runs). ARGS are the deck
(Tests/Canonical_RANS/Unsteady_Hill_2D/inputs_unsteady_hill2d) and the flags.
The script runs the deck from the log-law first guess in np2/, np1/ and np4/
(each emptied first), with the probes of probes_hill2d.txt, the averages from
t = 90 s (mid-step 2) and a plotfile and a checkpoint every five steps, and
checks:
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
    the short run's (nearly the first guess) by more than 0.1 m/s;
  - RT-6.3.1 (sub-phase 6.3, beyond the plan), the outputs: the probe file has
    a row at t = 0 and one per step; the plotfiles of step 5 (t = 300 s) and
    step 10 carry u_mean ... tke_mean and u_var ... tke_var, which at every
    probe's cell equal the weighted mean and variance of that probe's rows
    computed here in two passes, with the weights 30 s for step 2 (90..120 s)
    and 60 s for the later steps;
  - the restart (in restart/): from the 2-rank checkpoint of step 5 to 600 s,
    with the 2-rank probe file copied in. The log continues from t = 300 s at
    step 6 with the 2-rank run's inner counts; the final plotfile (state and
    averages) and the probe file are byte for byte the 2-rank run's, and so are
    the final checkpoint's fields, its Unsteady record and averages.
The exit code is the verdict.
"""
import filecmp
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

FIELDS = ["u", "v", "w", "p", "theta", "tke"]
AVG = ["u", "v", "w", "theta", "tke"]
PROBES = os.path.join(HERE, "probes_hill2d.txt")
EPS = 2.0 ** -52
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


def plotfiles(stdout):
    return [ln.split()[-1] for ln in stdout.splitlines()
            if ln.startswith("gust: wrote ") and "_plt" in ln and "_level" not in ln]


def checkpoints(stdout):
    return [ln.split()[-1] for ln in stdout.splitlines() if re.match(r"^gust: wrote chk\d+$", ln)]


def read_probes(path):
    """The probe file: the cells (i, j, k) from the header and the rows (t, values...)."""
    cells, rows = [], []
    with open(path) as fh:
        for ln in fh:
            if ln.startswith("# probe "):
                arrow = ln.split("->")[1].split()
                cells.append(tuple(int(x) for x in arrow[:3]))
            elif ln[0].isdigit():
                rows.append([float(x) for x in ln.split(",")])
    return cells, rows


def hand_averages(rows, start):
    """Two-pass weighted mean and variance of every column after the time start: the weight of the
    row at t^{n+1} is t^{n+1} - max(t^n, start), when > 0."""
    weights = []
    for prev, row in zip(rows, rows[1:]):
        w = row[0] - max(prev[0], start)
        if w > 0:
            weights.append((w, row[1:]))
    total = sum(w for w, _ in weights)
    ncol = len(rows[0]) - 1
    mean = [sum(w * x[c] for w, x in weights) / total for c in range(ncol)]
    var = [sum(w * (x[c] - mean[c]) ** 2 for w, x in weights) / total for c in range(ncol)]
    return mean, var, total


def same_tree(a, b, skip=()):
    """The files of directory a that differ from b's (or are missing in b)."""
    bad = []
    for root, _, files in os.walk(a):
        for f in files:
            pa = os.path.join(root, f)
            rel = os.path.relpath(pa, a)
            if rel in skip:
                continue
            pb = os.path.join(b, rel)
            if not os.path.exists(pb) or not filecmp.cmp(pa, pb, shallow=False):
                bad.append(rel)
    return bad


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
    outputs = ["unsteady.stop_time=600", "unsteady.dt=60", "unsteady.average_start=90", "unsteady.plot_interval=5",
               "unsteady.checkpoint_interval=5", "unsteady.probes_file=" + PROBES]
    ten = ["init.first_guess=profile"] + outputs
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

    # RT-6.3.1: the averages against the probes' rows.
    cells, rows = read_probes(os.path.join("np2", "unsteady_hill2d_plt_probes.csv"))
    check(len(rows) == 11 and abs(rows[0][0]) == 0 and abs(rows[-1][0] - 600) < 1e-9,
          "the probe file has 11 rows, t = 0 to 600 s (got %d)" % len(rows))
    plts = plotfiles(two)
    check(len(plts) == 2, "two plotfiles, steps 5 and 10 (got %d)" % len(plts))
    for plt, upto in zip(plts, (300.0, 600.0)):
        path = os.path.join("np2", plt)
        check(abs(plotfile_time(path) - upto) < 1e-9, "%s carries t = %g s" % (plt, upto))
        _, avg = ep.read_fields(path, [f + "_mean" for f in AVG] + [f + "_var" for f in AVG])
        mean, var, total = hand_averages([r for r in rows if r[0] <= upto + 1e-9], 90.0)
        worst_m = worst_v = 0.0
        for n, (i, j, k) in enumerate(cells):
            for c, f in enumerate(AVG):
                col = n * len(AVG) + c
                series = [r[1 + col] for r in rows]
                scale = max(abs(x) for x in series) or 1.0
                worst_m = max(worst_m, abs(avg[f + "_mean"][i][j][k] - mean[col]) / scale)
                # Round-off of the deviations d (about eps scale each) bounds the variance's error by
                # about eps scale sqrt(var): theta's variance is 1e-11 on a scale of 300.
                tol = 1e-12 * var[col] + 64 * EPS * scale * math.sqrt(var[col]) + 1e-300
                worst_v = max(worst_v, abs(avg[f + "_var"][i][j][k] - var[col]) / tol)
        check(worst_m <= 1e-13 and worst_v <= 1,
              "%s: the averages at the probes are the rows' weighted mean (%.1e of the scale) and variance (%.2f of "
              "its round-off bound) over %g s" % (plt, worst_m, worst_v, total))

    # RT-6.3.1: the restart from the checkpoint of step 5.
    chks = checkpoints(two)
    check(len(chks) == 2, "two checkpoints, steps 5 and 10 (got %d)" % len(chks))
    shutil.rmtree("restart", ignore_errors=True)
    os.makedirs("restart")
    shutil.copy(os.path.join("np2", "unsteady_hill2d_plt_probes.csv"), "restart")
    full = with_ranks(cmd, 2) + ["io.restart=../np2/" + chks[0]] + outputs
    print("check_unsteady_hill2d: in restart/: %s" % " ".join(full))
    out = subprocess.run(full, cwd="restart", check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_unsteady_hill2d: the restart failed with code %d" % out.returncode)
    rs = steps(out.stdout)
    check("gust: unsteady: continues from t = 300 s after 5 physical steps" in out.stdout, "the restart continues at t = 300 s")
    check([m.group(1) for m in rs] == [str(n) for n in range(6, 11)]
          and [m.group(4) for m in rs] == [m.group(4) for m in s2[5:]],
          "the restart runs steps 6 to 10 with the 2-rank inner counts")
    rp = plotfiles(out.stdout)
    check(bool(rp) and rp[-1] == plts[-1] and not same_tree(os.path.join("np2", plts[-1]), os.path.join("restart", rp[-1])),
          "the restart's final plotfile is byte for byte the 2-rank run's")
    check(filecmp.cmp(os.path.join("np2", "unsteady_hill2d_plt_probes.csv"),
                      os.path.join("restart", "unsteady_hill2d_plt_probes.csv"), shallow=False),
          "the restart's probe file is byte for byte the 2-rank run's")
    rc = checkpoints(out.stdout)
    # The driver history records wall-clock times; everything else is compared.
    diff = same_tree(os.path.join("np2", chks[-1]), os.path.join("restart", rc[-1]), skip=("Driver",)) if rc else ["(none)"]
    check(not diff, "the restart's final checkpoint is byte for byte the 2-rank run's but Driver (differs: %s)"
          % (" ".join(diff) or "nothing"))

    print("check_unsteady_hill2d: %s" % ("pass" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
