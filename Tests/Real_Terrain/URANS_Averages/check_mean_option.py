#!/usr/bin/env python3
"""The time-mean option of the plotfile column reader (tools/askervein_lines.py,
Columns(mean=True); `askervein_lines.py sample --mean`, `rushil.py compare --mean`).

    check_mean_option.py --prepare RUNDIR DECK
    check_mean_option.py RUNDIR PREFIX

--prepare empties RUNDIR (a plotfile of an earlier run would enter the mean) and
writes RUNDIR/inputs_two_in_y, DECK with mesh.n_cell = 128 2 64: the mast
sampler interpolates between four columns, so the run needs two cells in y.

RUNDIR holds the plotfiles PREFIX<step> of an unsteady run with
unsteady.average_start = 0, a constant unsteady.dt and unsteady.plot_interval = 1.
With these the running mean of the last plotfile is the arithmetic mean of the
fields of every step, so the check needs no stored number:
  - Columns(last, mean=True) gives, in every cell of three columns, the mean
    over the plotfiles of Columns(plotfile) (u and v; relative 1e-10);
  - the mean differs from the last step's field by more than 1e-3 of the
    largest |u| (the run is still moving: a reader that returned u for u_mean
    would pass the first check on a field at rest);
  - z is the same with and without the option;
  - the command lines forward the option: `askervein_lines.py sample --mean`
    at two masts (a terrain file and a lines file written here; the run has
    two cells in y) gives the speed of the mean columns checked above, and
    differs from plain `sample` of the last plotfile; `check --mean` passes
    against those speeds and plain `check` fails; `rushil.py compare
    --mean` (a station file written here) prints the U_ref of the mean
    columns, which differs from the one without the option.
A NaN in a sampled mean fails. Pure Python (the CTest rule of this repository).
"""
import csv
import glob
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
import askervein_lines as al  # noqa: E402


def prepare(rundir, deck):
    shutil.rmtree(rundir, ignore_errors=True)
    os.makedirs(rundir)
    text, n = re.subn(r"(?m)^mesh\.n_cell\s*=.*$", "mesh.n_cell = 128 2 64", open(deck).read())
    if n != 1:
        sys.exit(f"FAIL: {deck}: {n} lines mesh.n_cell (one expected)")
    with open(os.path.join(rundir, "inputs_two_in_y"), "w") as f:
        f.write(text)


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "--prepare":
        return prepare(*sys.argv[2:])
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    rundir, prefix = sys.argv[1:]
    stem = os.path.join(rundir, prefix)
    plts = sorted((p for p in glob.glob(stem + "[0-9]*") if os.path.isdir(p) and p[len(stem):].isdigit()),
                  key=lambda p: int(p[len(stem):]))
    if len(plts) < 3:
        sys.exit(f"FAIL: {len(plts)} plotfiles {prefix}* in {rundir}; at least 3 steps are needed")
    inst = [al.Columns(p) for p in plts]
    mean = al.Columns(plts[-1], mean=True)
    h = mean.hdr
    nx, ny = h["hi"][0] - h["lo"][0] + 1, h["hi"][1] - h["lo"][1] + 1
    worst, moved, scale = 0.0, 0.0, 0.0
    for i in (nx // 4, nx // 2, (3 * nx) // 4):
        j = ny // 2
        um, vm, zm = mean.column(i, j)
        cols = [c.column(i, j) for c in inst]
        if zm != cols[-1][2]:
            sys.exit(f"FAIL: column {i}: z differs with the option")
        for k in range(len(um)):
            ua = sum(c[0][k] for c in cols) / len(cols)
            va = sum(c[1][k] for c in cols) / len(cols)
            if not (math.isfinite(um[k]) and math.isfinite(vm[k])):
                sys.exit(f"FAIL: column {i}, cell {k}: the mean is not finite")
            scale = max(scale, abs(ua))
            worst = max(worst, abs(um[k] - ua), abs(vm[k] - va))
            moved = max(moved, abs(um[k] - cols[-1][0][k]))
    print(f"{len(plts)} plotfiles: |mean - average of the steps| <= {worst:.3e} m/s, "
          f"|mean - last step| up to {moved:.3e} m/s, largest |u| {scale:.3f} m/s")
    if not worst <= 1e-10 * scale:
        sys.exit("FAIL: the mean read with the option is not the average of the steps' fields")
    if not moved > 1e-3 * scale:
        sys.exit("FAIL: the mean equals the last step's field; the check cannot tell u_mean from u")
    command_lines(rundir, plts, mean)
    print("PASS")


def run(args, expect=0):
    r = subprocess.run([sys.executable] + args, capture_output=True, text=True)
    if r.returncode != expect:
        sys.exit(f"FAIL: exit code {r.returncode} (expected {expect}) of {' '.join(args)}\n{r.stdout}{r.stderr}")
    return r.stdout


def speeds(path):
    with open(path) as fh:
        return [float(r["speed"]) for r in csv.DictReader(ln for ln in fh if not ln.startswith("#"))]


def command_lines(rundir, plts, mean):
    h = mean.hdr
    lo, dx = h["prob_lo"], h["dx"]
    nx, ny = h["hi"][0] - h["lo"][0] + 1, h["hi"][1] - h["lo"][1] + 1
    if ny < 2:
        sys.exit("FAIL: the sampler needs two cells in y (mesh.n_cell of the run)")
    # ERF's terrain format, here the ground of the lowest cell faces is not known to the script: the
    # sampler only needs a ground below the first cell centre, so the hill of the deck (Agnesi, 100 m, 500 m).
    xs = [lo[0] + i * dx[0] for i in range(nx + 1)]
    ys = [lo[1] + j * dx[1] for j in range(ny + 1)]
    terrain, lines = os.path.join(rundir, "mean_terrain.txt"), os.path.join(rundir, "mean_lines.csv")
    with open(terrain, "w") as f:
        f.write(f"{len(xs)} {len(ys)}\n" + " ".join(f"{x:.3f}" for x in xs) + "\n" + " ".join(f"{y:.3f}" for y in ys) + "\n")
        f.write(" ".join(f"{100.0 / (1.0 + (x / 500.0) ** 2):.4f}" for x in xs for _ in ys) + "\n")
    ymid = lo[1] + 0.5 * ny * dx[1]
    with open(lines, "w") as f:
        f.write("line,d,x,y\nA,0,0.0,%.3f\nA,600,600.0,%.3f\n" % (ymid, ymid))
    tool = os.path.join(HERE, "..", "..", "..", "tools", "askervein_lines.py")
    base = [tool, "sample", "--terrain", terrain, "--lines", lines, "--height", "50"]
    out = os.path.join(rundir, "mean_last.csv")
    run(base[:2] + [plts[-1]] + base[2:] + ["--out", out])
    last = speeds(out)
    out = os.path.join(rundir, "mean_sample.csv")
    run(base[:2] + [plts[-1]] + base[2:] + ["--out", out, "--mean"])
    got = speeds(out)
    if "time mean" not in open(out).readline():
        sys.exit("FAIL: the CSV of sample --mean does not say it holds the time mean")
    # The speed is not linear in u, v (v is not zero: Coriolis), so the expected value is the sampler's on
    # the mean columns checked above, not the average of the steps' samples.
    grid = al.read_terrain(terrain)
    avg = [al.speed_at(mean, grid, x, ymid, 50.0, 0.03) for x in (0.0, 600.0)]
    d_avg = max(abs(g - a) for g, a in zip(got, avg))
    d_last = max(abs(g - s) for g, s in zip(got, last))
    print(f"sample --mean at {len(got)} masts: |speed - speed of the mean columns| <= {d_avg:.2e} m/s, "
          f"|speed - last step's| up to {d_last:.3e} m/s")
    if not (len(got) == 2 and d_avg <= 2e-6 and d_last > 1e-3):      # the CSV carries six decimals
        sys.exit("FAIL: sample --mean does not sample the mean columns, or equals the last step's sample")
    ref = os.path.join(rundir, "mean_reference.csv")
    with open(ref, "w") as f:
        f.write("line,d,x,y,speed\nA,0,0,0,%.6f\nA,600,0,0,%.6f\n" % tuple(avg))
    chk = [tool, "check", plts[-1], ref, "--terrain", terrain, "--lines", lines, "--height", "50",
           "--compare_lines", "A", "--rel_tol", "1e-5", "--out", os.path.join(rundir, "mean_check.csv")]
    run(chk + ["--mean"])
    run(chk, expect=1)
    # rushil.py compare: the hill of this run is not a RUSHIL hill; only U_ref of the inflow column is used.
    with open(os.path.join(rundir, "stations_h5.csv"), "w") as f:
        f.write("x,z,U\n-4480.0,50.0,5.0\n")
    rushil = os.path.join(HERE, "..", "RUSHIL", "rushil.py")
    cmp_ = [rushil, "compare", plts[-1], "--hill", "5", "--outdir", rundir]
    uref = lambda text: float(text.rsplit("U_ref(0.8H)", 1)[1].split()[0])
    with_mean, without = run(cmp_ + ["--mean"]), run(cmp_)
    sys.path.insert(0, os.path.dirname(rushil))
    import rushil as ru
    us, zs = ru.column(mean, 200.0)
    want = ru.u_at(us, zs, 0.8 * ru.H)
    print(f"rushil.py compare: U_ref {uref(with_mean):.3f} m/s with --mean ({want:.3f} from the mean columns), "
          f"{uref(without):.3f} m/s without")
    if "time mean" not in with_mean or abs(uref(with_mean) - want) > 6e-4 or abs(uref(with_mean) - uref(without)) < 2e-3:
        sys.exit("FAIL: rushil.py compare --mean does not read the mean")


if __name__ == "__main__":
    main()
