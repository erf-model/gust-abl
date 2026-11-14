#!/usr/bin/env python3
"""One-way nesting validation on Askervein: 20 m patches against the uniform 20 m run.

    askervein_nest.py run    NAME N          -- COMMAND [ARGS ...]
    askervein_nest.py nest   NAME PARENT D   -- COMMAND [ARGS ...]
    askervein_nest.py compare PATCH REFERENCE [PARENT] [--ratio R] [--rms X] [--below PATCH2 ...]
    askervein_nest.py study  -- MPIEXEC ... GUST

run: the whole 6 km domain (inputs_askervein_domain) on N x N x 18 cells (150:
40 m, the parent; 300: 20 m, the reference) to the steady verdict. The last
plotfile and checkpoint are copied to NAME_plt_final and NAME_chk_final, and
the 10 m speeds at the masts of lines A and AA are written to NAME_lines.csv
(tools/askervein_lines.py sample).

nest: a square patch of half-width D m centred on the hilltop HT (0, 0)
(inputs_askervein_patch), at 2 2 1 of the finished run PARENT; outputs as for
run. D must be a multiple of the parent's cell size and at least 900 m (the
masts reach 870 m from HT).

COMMAND is the solver with its deck; the terrain, the inflow table (run), the
mesh, the parent and the output names are appended.

compare: the 10 m speeds of PATCH and of REFERENCE at the 30 masts
(NAME_lines.csv of finished runs in the working directory). Prints each
mast's relative difference, the RMS and the largest over the masts, and with
PARENT the same for the parent. --rms X: the patch's RMS must be at most X.
--ratio R (needs PARENT): the patch's RMS must be at most R times the
parent's. --below names patches with nest faces nearer the hill, whose RMS
must not be smaller than this patch's (to 5 % of it). The exit code is the
verdict.

study: the 40 m parent, the 20 m reference and patches of D = 1000, 1400 and
2000 m, in the working directory; the table is printed and written to
askervein_nest.csv.
"""
import csv
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(ROOT, "Tests", "Canonical_RANS", "Nest_Validation"))
import nest_validation as nv  # noqa: E402  (finish, geometry: the runs of the hill study)

TERRAIN = os.path.join(ROOT, "extern_erf", "askervein.txt")
INFLOW = os.path.join(ROOT, "Tests", "Real_Terrain", "Askervein_210", "askervein_210_inflow.txt")
LINES = os.path.join(ROOT, "Tests", "Real_Terrain", "Askervein_210", "askervein_lines.csv")
SAMPLER = os.path.join(ROOT, "tools", "askervein_lines.py")
L = 3000.0       # the domain is (-L, L) in x and y
D_MIN = 900.0    # the farthest mast is 870 m from HT


def fail(msg):
    raise SystemExit("askervein_nest: " + msg)


def sample(name):
    proc = subprocess.run([sys.executable, SAMPLER, "sample", name + "_plt_final", "--terrain", TERRAIN,
                           "--lines", LINES, "--out", name + "_lines.csv"],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        print(proc.stdout)
        fail(f"{name}: the masts could not be sampled")


def run_domain(name, n, cmd):
    its = nv.finish(name, cmd + [f"mesh.n_cell={n} {n} 18", f"terrain.file={TERRAIN}", f"bc.inflow_file={INFLOW}"])
    sample(name)
    return its


def run_nest(name, parent, d, cmd, extra=()):
    pdx = nv.geometry(parent)[0]
    m = d / pdx
    if abs(m - round(m)) > 1e-9:
        fail(f"a patch of half-width {d} m does not lie on the {pdx} m cells of {parent}")
    if d < D_MIN or d >= L:
        fail(f"the half-width must be from {D_MIN:g} m to below {L:g} m")
    n = 4 * int(round(m))
    its = nv.finish(name, cmd + [f"mesh.n_cell={n} {n} 18", f"mesh.prob_lo={-d} {-d} 0.", f"mesh.prob_hi={d} {d} 440.",
                                 f"terrain.file={TERRAIN}", f"bc.nest_parent={parent}_chk_final"] + list(extra))
    sample(name)
    return its


def speeds(name):
    with open(name + "_lines.csv") as fh:
        return {(r["line"], int(r["d"])): float(r["speed"])
                for r in csv.DictReader(ln for ln in fh if not ln.startswith("#"))}


def differences(name, ref):
    """(RMS, largest, per-mast) relative difference of the 10 m speed from the reference's."""
    a, b = speeds(name), speeds(ref)
    if set(a) != set(b) or not a:
        fail(f"{name}_lines.csv and {ref}_lines.csv do not hold the same masts")
    rel = {k: (a[k] - b[k]) / b[k] for k in b}
    vals = list(rel.values())
    rms = math.sqrt(sum(v * v for v in vals) / len(vals))
    mx = 0.0
    for v in vals:
        if not abs(v) <= mx:   # a NaN is kept
            mx = abs(v)
    return rms, mx, rel


def compare(args):
    names, ratio, rms_max, below = [], None, None, []
    it = iter(args)
    for a in it:
        if a == "--ratio":
            ratio = float(next(it))
        elif a == "--rms":
            rms_max = float(next(it))
        elif a == "--below":
            below = list(it)
        else:
            names.append(a)
    if len(names) not in (2, 3) or (ratio is not None and len(names) != 3):
        print(__doc__)
        return 2
    patch, ref = names[0], names[1]
    rms, mx, rel = differences(patch, ref)
    prel = None
    if len(names) == 3:
        prms, pmx, prel = differences(names[2], ref)
    for k in sorted(rel):
        extra = "" if prel is None else f"; {names[2]} {100 * prel[k]:+7.3f} %"
        print(f"askervein_nest: line {k[0]:>2} d = {k[1]:5d} m: {patch} {100 * rel[k]:+7.3f} %{extra}")
    print(f"askervein_nest: {patch} vs {ref}: {len(rel)} masts, RMS {100 * rms:.3f} %, largest {100 * mx:.3f} %")
    bad = []
    if rms_max is not None and not rms <= rms_max:
        bad.append(f"the RMS {rms:.4e} is above {rms_max:g}")
    if prel is not None:
        print(f"askervein_nest: {names[2]} vs {ref}: RMS {100 * prms:.3f} %, largest {100 * pmx:.3f} %")
        if ratio is not None and not rms <= ratio * prms:
            bad.append(f"the patch's RMS {rms:.4e} is above {ratio:g} x its parent's {prms:.4e}")
    for other in below:
        orms = differences(other, ref)[0]
        print(f"askervein_nest: {other} vs {ref}: RMS {100 * orms:.3f} %")
        if not rms <= 1.05 * orms:
            bad.append(f"the RMS {rms:.4e} is above that of {other} ({orms:.4e}), whose nest faces are nearer the hill")
    for b in bad:
        print("askervein_nest: FAIL: " + b)
    if not bad:
        print("askervein_nest: PASS")
    return 1 if bad else 0


def study(solver):
    dom = solver + [os.path.join(HERE, "inputs_askervein_domain")]
    nst = solver + [os.path.join(HERE, "inputs_askervein_patch")]
    rows = [("par40", 150 * 150 * 18, run_domain("par40", 150, dom)),
            ("ref20", 300 * 300 * 18, run_domain("ref20", 300, dom))]
    for d in (1000, 1400, 2000):
        n = d // 10
        rows.append((f"nest_d{d}", n * n * 18, run_nest(f"nest_d{d}", "par40", float(d), nst)))
    with open("askervein_nest.csv", "w") as f:
        f.write("run,cells,iterations,rms,largest\n")
        for name, ncell, its in rows:
            rms, mx, _ = differences(name, "ref20")
            print(f"askervein_nest: {name}: {ncell} cells, {its} iterations, RMS {100 * rms:.3f} %, "
                  f"largest {100 * mx:.3f} % from ref20")
            f.write(f"{name},{ncell},{its},{rms:.6e},{mx:.6e}\n")
    print("askervein_nest: wrote askervein_nest.csv")
    return 0


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    if argv[0] == "compare":
        return compare(argv[1:])
    if "--" not in argv:
        print(__doc__)
        return 2
    sep = argv.index("--")
    args, cmd = argv[:sep], argv[sep + 1:]
    if not cmd:
        print(__doc__)
        return 2
    if args[0] == "run" and len(args) == 3:
        run_domain(args[1], int(args[2]), cmd)
    elif args[0] == "nest" and len(args) == 4:
        run_nest(args[1], args[2], float(args[3]), cmd)
    elif args[0] == "study" and len(args) == 1:
        return study(cmd)
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
