#!/usr/bin/env python3
"""One-way nesting validation on Askervein: 20 m patches against the uniform 20 m run.

    askervein_nest.py run    NAME N          -- COMMAND [ARGS ...]
    askervein_nest.py nest   NAME PARENT D   -- COMMAND [ARGS ...]
    askervein_nest.py compare PATCH REFERENCE [PARENT] [--ratio R] [--rms X] [--below PATCH2 ...]
    askervein_nest.py study  -- MPIEXEC ... GUST [-- ARGS ...]
    askervein_nest.py two_way [relax] [inject] [relax_consistent] [keep] -- COMMAND [ARGS ...]
    askervein_nest.py two_way paper [relax] [inject] [relax_consistent] [keep] -- MPIEXEC ... GUST [ARGS ...]

run: the whole 6 km domain (inputs_askervein_domain) on N x N x 18 cells (150:
40 m, the parent; 300: 20 m, the reference) to the steady verdict. The last
plotfile and checkpoint are copied to NAME_plt_final and NAME_chk_final, and
the 10 m speeds at the masts of lines A and AA are written to NAME_lines.csv
(tools/askervein_lines.py sample).

nest: a square patch of half-width D m centred on the hilltop HT (0, 0)
(inputs_askervein_patch), at 2 2 1 of the finished run PARENT; outputs as for
run. D must be a multiple of the parent's cell size and at least 900 m (the
masts reach 870 m from HT in x).

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

two_way: two-way nesting (nest.feedback = defect, PR-70) in the working
directory: the 20 m reference (ref_*), and the 40 m outer domain with the
nest ht (-1000, -1000) to (1000, 1000) m, one way (ow_*) and joint (jt_*),
one after another, each with the recommended setting (TWO_WAY). Prints the
RMS of u, v, w and k against the reference over the nest (the reference
interpolated to the nest's cell heights: a joint nest's levels follow the
outer domain's) and over the outer cells outside its footprint, the cells of the lowest 200 m whose
horizontal speed differs from the reference's by more than 0.25, 0.5 and
1 m/s (tools/speed_difference.py), the RMS of the masts' 10 m speeds against the
reference's, and the RMS of the speed-up on lines A and AA against the
observations (digitized from Wagenbrenner et al. 2019, Fig. 8, in
Tests/Real_Terrain/WindNinja2019/wagenbrenner2019_digitized.csv; speed-up
relative to the inflow's 10 m speed). No verdict: the numbers are those of
docs/boundary_conditions.rst. "relax", "inject" and "relax_consistent" add
the joint runs with nest.feedback = relax (rx_*) at nest.feedback_factor = 1
(the default before tracker N-114 item 13; now 0.5), inject (ij_*) and
relax_consistent (rc_*) at their default nest.feedback_factor (0.5 and 1),
compared with the others; "keep" skips a
run whose plotfile is already in the working directory. With "paper"
COMMAND names no deck: the runs
take a copy of Tests/Real_Terrain/WindNinja2019/inputs_askervein without the
keys the runs set (PAPER_SET; inputs_askervein_paper; the vertical mesh of the
Askervein comparison: a 1 m first layer, 64 levels; the log-law inflow): the
reference on 256 x 256 cells (23.4 m), the outer domain on 128 x 128 and the
nest ht (-937.5, -937.5) to (937.5, 937.5) m, on outer cell faces.

study: the 40 m parent, the 20 m reference and patches of D = 1000, 1400 and
2000 m, in the working directory (ARGS are given to every run after its
deck); the table is printed and written to
askervein_nest.csv. Targets A1 to A3 of the README are judged; the exit code
is the verdict.
"""
import csv
import math
import os
import re
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
PAPER_DECK = os.path.join(ROOT, "Tests", "Real_Terrain", "WindNinja2019", "inputs_askervein")
PAPER_SET = ("mesh.n_cell", "io.write_checkpoint", "io.plot_file", "io.chk_file")   # keys the runs set
OBSERVED = os.path.join(ROOT, "Tests", "Real_Terrain", "WindNinja2019", "wagenbrenner2019_digitized.csv")
L = 3000.0       # the domain is (-L, L) in x and y
D_MIN = 900.0    # the masts reach 870 m from HT in x


def fail(msg):
    raise SystemExit("askervein_nest: " + msg)


def unsample(name):
    """Remove the mast table of an earlier run of this name before the run starts."""
    if os.path.exists(name + "_lines.csv"):
        os.remove(name + "_lines.csv")


def sample(name):
    proc = subprocess.run([sys.executable, SAMPLER, "sample", name + "_plt_final", "--terrain", TERRAIN,
                           "--lines", LINES, "--out", name + "_lines.csv"],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        print(proc.stdout)
        fail(f"{name}: the masts could not be sampled")


def run_domain(name, n, cmd):
    unsample(name)
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
    unsample(name)
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
        mx = nv.larger(mx, abs(v))   # a NaN stays
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
    flags = []
    if "--" in solver:
        flags = solver[solver.index("--") + 1:]
        solver = solver[:solver.index("--")]
    dom = solver + [os.path.join(HERE, "inputs_askervein_domain")] + flags
    nst = solver + [os.path.join(HERE, "inputs_askervein_patch")] + flags
    rows = [("par40", 150 * 150 * 18, run_domain("par40", 150, dom)),
            ("ref20", 300 * 300 * 18, run_domain("ref20", 300, dom))]
    unsteady = []
    for d in (1000, 1400, 2000):
        n = d // 10   # 20 m cells
        try:
            rows.append((f"nest_d{d}", n * n * 18, run_nest(f"nest_d{d}", "par40", float(d), nst)))
        except SystemExit as err:
            print(err)
            unsteady.append(d)
    rms_of = {}
    with open("askervein_nest.csv", "w") as f:
        f.write("run,cells,iterations,rms,largest\n")
        for name, ncell, its in rows:
            rms, mx, _ = differences(name, "ref20")
            rms_of[name] = rms
            print(f"askervein_nest: {name}: {ncell} cells, {its} iterations, RMS {100 * rms:.3f} %, "
                  f"largest {100 * mx:.3f} % from ref20")
            f.write(f"{name},{ncell},{its},{rms:.6e},{mx:.6e}\n")
    print("askervein_nest: wrote askervein_nest.csv")
    bad = []

    def judge(target, ok, text):
        print(f"askervein_nest: {target}: {'met' if ok else 'MISSED'}: {text}")
        if not ok:
            bad.append(target)

    par = rms_of["par40"]
    nan = float("nan")   # a patch that was not finished misses its targets
    for d in (1000, 1400, 2000):
        rms_of.setdefault(f"nest_d{d}", nan)
    for d in (1400, 2000):
        r = rms_of[f"nest_d{d}"]
        judge("A1", r <= 0.5 * par, f"D = {d} m: RMS {100 * r:.3f} %, the parent's {100 * par:.3f} %")
    for d, e in ((1400, 1000), (2000, 1400)):
        r, q = rms_of[f"nest_d{d}"], rms_of[f"nest_d{e}"]
        judge("A2", r <= 1.05 * q, f"D = {d} m: RMS {100 * r:.3f} %, {100 * q:.3f} % at D = {e} m")
    warned = []
    for d in (1000, 1400, 2000):
        if d in unsteady:
            continue
        with open(f"nest_d{d}.log") as f:
            if re.search(r"WARNING: nest faces: the net flux", f.read()):
                warned.append(d)
    judge("A3", not warned and not unsteady,
          "patches short of the steady verdict: " + (", ".join(f"D = {d} m" for d in unsteady) if unsteady else "none")
          + "; net-flux warnings: " + (", ".join(f"D = {d} m" for d in warned) if warned else "none"))
    print("askervein_nest: " + ("FAIL: missed " + ", ".join(sorted(set(bad))) if bad else "PASS"))
    return 1 if bad else 0


# The two-way runs: the recommended setting (docs: Recommended settings) and the nest's checkpoint.
TWO_WAY = ["turb.realizability=0.309", "turb.nut_relax=0.5", "convection.hybrid_sigma_max=0.1",
           "convection.hybrid_height=290.", "turb.max_geom_lscale=100.", "io.write_checkpoint=1"]
TWO_WAY_NEST = ["nest.names=ht", "nest.ht.lo=-1000 -1000", "nest.ht.hi=1000 1000"]


def deck_without_cap(deck, copy):
    """A copy of the deck without turb.max_geom_lscale, which TWO_WAY gives (a key given twice aborts)."""
    with open(deck) as fh, open(copy, "w") as out:
        out.writelines(line for line in fh if line.split("=")[0].strip() != "turb.max_geom_lscale")
    return copy


def two_way(cmd, paper=False, extra=(), keep=False):
    import glob
    import numpy as np
    sys.path.insert(0, os.path.join(ROOT, "extern_erf"))
    import erf_plotfile
    if paper:   # the deck sets the inflow and the recommended setting; its copy leaves out what the runs set
        with open(PAPER_DECK) as fh, open("inputs_askervein_paper", "w") as out:
            out.writelines(line for line in fh if line.split("=")[0].strip() not in PAPER_SET)
        at = next((i for i, a in enumerate(cmd) if "=" in a and not a.startswith("-")), len(cmd))  # deck before ARGS
        base = cmd[:at] + ["inputs_askervein_paper"] + cmd[at:] + [f"terrain.file={TERRAIN}", "io.write_checkpoint=1"]
        n_ref, n_out, half, nz = 256, 128, 937.5, 64
    else:
        cmd = [deck_without_cap(a, "inputs_askervein_two_way") if os.path.basename(a) == "inputs_askervein_domain"
               else a for a in cmd]   # the deck sets the 30 m cap of the one-way runs
        base = cmd + [f"terrain.file={TERRAIN}", f"bc.inflow_file={INFLOW}"] + TWO_WAY
        n_ref, n_out, half, nz = 300, 150, 1000.0, 18
    nest = ["nest.names=ht", f"nest.ht.lo={-half:g} {-half:g}", f"nest.ht.hi={half:g} {half:g}"]
    runs = [("ref", [f"mesh.n_cell={n_ref} {n_ref} {nz}"]),
            ("ow", [f"mesh.n_cell={n_out} {n_out} {nz}"] + nest + ["nest.feedback=none"]),
            ("jt", [f"mesh.n_cell={n_out} {n_out} {nz}"] + nest + ["nest.feedback=defect"])]
    # relax at factor 1, the default before tracker N-114 item 13 (now 0.5), so the README's table stays
    # reproducible; inject and relax_consistent at their defaults (0.5 and 1).
    joint = {"rx": ("relax", ["nest.feedback_factor=1"]), "ij": ("inject", []), "rc": ("relax_consistent", [])}
    runs += [(k, [f"mesh.n_cell={n_out} {n_out} {nz}"] + nest + [f"nest.feedback={v}"] + fa)
             for k, (v, fa) in joint.items() if v in extra]
    for name, args in runs:
        if keep and glob.glob(name + "_plt[0-9]*"):
            continue
        with open(name + ".log", "w") as fh:
            rc = subprocess.run(base + args + [f"io.plot_file={name}_plt", f"io.chk_file={name}_chk"],
                                stdout=fh, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            fail(f"{name}: the run exited {rc}; see {name}.log")

    def last(prefix):
        found = [(int(m.group(1)), p) for p in glob.glob(prefix + "[0-9]*")
                 for m in [re.fullmatch(re.escape(prefix) + r"(\d+)", p)] if m and os.path.isdir(p)]
        if not found:
            fail(f"no plotfile {prefix}*")
        return max(found)[1]

    fields = ["u", "v", "w", "tke"]
    rd = lambda p: {k: np.array(v, dtype=float) for k, v in erf_plotfile.read_fields(p, fields + ["z"])[1].items()}
    tags = [name for name, _ in runs if name != "ref"]   # ow, jt and the extra joint runs
    plt = {"ref": last("ref_plt")}
    for t in tags:
        plt[t], plt[t + "_ht"] = last(t + "_plt"), last(t + "_plt_ht")
    R = rd(plt["ref"])
    O = {t: rd(plt[t]) for t in tags}
    N = {t: rd(plt[t + "_ht"]) for t in tags}
    dx = 6000.0 / n_out                  # the footprint: outer cells i0 .. i0 + ni - 1 in x and y
    i0, ni = round((3000.0 - half) / dx), round(2 * half / dx)
    co = lambda a: 0.25 * (a[0::2, 0::2] + a[1::2, 0::2] + a[0::2, 1::2] + a[1::2, 1::2])
    rms = lambda a: float(np.sqrt(np.mean(a * a)))
    print(f"{'field':5s} " + " ".join(f"{'nest ' + t:>10s}" for t in tags) + " "
          + " ".join(f"{'out ' + t:>10s}" for t in tags))
    for f in fields:
        rc, rn = co(R[f]), R[f][2 * i0:2 * (i0 + ni), 2 * i0:2 * (i0 + ni)]
        zr = R["z"][2 * i0:2 * (i0 + ni), 2 * i0:2 * (i0 + ni)]

        def at(n):   # the reference at the nest's cell heights, column by column (a joint nest's levels move)
            out_ = np.empty_like(n["z"])
            for a in range(out_.shape[0]):
                for b in range(out_.shape[1]):
                    out_[a, b] = np.interp(n["z"][a, b], zr[a, b], rn[a, b])
            return out_
        out = np.ones(rc.shape[:2], bool)
        out[i0:i0 + ni, i0:i0 + ni] = False
        out = np.broadcast_to(out[:, :, None], rc.shape)
        print(f"{f:5s} " + " ".join(f"{rms(N[t][f] - at(N[t])):10.4e}" for t in tags) + " "
              + " ".join(f"{rms((O[t][f] - rc)[out]):10.4e}" for t in tags))
    # Horizontal speed differences in m/s in the lowest 200 m (tools/speed_difference.py), nest at its own heights.
    sys.path.insert(0, os.path.join(ROOT, "tools"))
    import speed_difference as sd
    sl = (slice(2 * i0, 2 * (i0 + ni)), slice(2 * i0, 2 * (i0 + ni)))
    Rs = sd.speed(R)
    Rc, mc = sd.coarsen(Rs), sd.near_ground(sd.coarsen(R["z"]))
    outer = np.ones(Rc.shape[:2], bool)
    outer[i0:i0 + ni, i0:i0 + ni] = False
    rows = []
    label = {"ow": "one-way", "jt": "joint", "rx": "relax, f = 1", "ij": "inject", "rc": "relax_consistent"}
    for t in tags:
        n_, o_ = N[t], O[t]
        rows.append((label[t], "nest",
                     (sd.speed(n_) - sd.at_height(n_["z"], R["z"][sl], Rs[sl]))[sd.near_ground(n_["z"])]))
        rows.append((label[t], "outer", (sd.speed(o_) - Rc)[mc & outer[:, :, None]]))
    print(f"horizontal speed |S - S_ref| below {sd.AGL:g} m above the ground "
          f"(reference mean {Rs[sd.near_ground(R['z'])].mean():.2f} m/s)")
    sd.print_table(rows)
    sampled = set(plt)
    for k, p in plt.items():
        proc = subprocess.run([sys.executable, SAMPLER, "sample", p, "--terrain", TERRAIN, "--lines", LINES,
                               "--out", f"tw_{k}_lines.csv"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                              text=True)
        if proc.returncode != 0:   # the outer domain's ground can lie above the terrain file's at a mast
            print(f"{p}: the masts are not sampled ({proc.stdout.strip().splitlines()[-1]})")
            sampled.discard(k)
    for k in [k for k in tags + [t + "_ht" for t in tags] if k in sampled]:
        r, mx, _ = differences("tw_" + k, "tw_ref")
        print(f"masts, {k}: 10 m speed relative to the reference: RMS {r:.4f}, largest {mx:.4f}")
    uref = 8.9 if paper else inflow_speed(10.0)
    obs = {}
    with open(OBSERVED) as fh:
        for row in csv.reader(line for line in fh if line[0] != "#"):
            if row[0] == "fig8" and row[2] == "observed" and row[1] in ("LineA", "LineAA"):
                obs.setdefault(row[1][4:], []).append((float(row[3]), float(row[4])))
    print(f"speed-up on the masts against the observations (U_ref = {uref:.3f} m/s, the inflow at 10 m):")
    for k in [k for k in ["ref"] + tags + [t + "_ht" for t in tags] if k in sampled]:
        with open(f"tw_{k}_lines.csv") as fh:
            rows = list(csv.DictReader(line for line in fh if line[0] != "#"))
        for line in sorted(obs):
            ds = []
            for d, so in obs[line]:   # the nearest mast within 15 m of each observation
                near = [r for r in rows if r["line"] == line and abs(float(r["d"]) - d) < 15.0]
                if near:
                    r = min(near, key=lambda r: abs(float(r["d"]) - d))
                    ds.append(float(r["speed"]) / uref - 1.0 - so)
            print(f"  {k:5s} line {line:2s}: {len(ds)} masts, RMS {rms(np.array(ds)):.4f}, "
                  f"largest {max(ds, key=abs):+.4f}")
    return 0


def inflow_speed(z):
    """The speed of the inflow table INFLOW at height z (linear in ln z between the rows around it)."""
    rows = [[float(v) for v in line.split()] for line in open(INFLOW) if line.strip() and line[0] != "#"]
    for a, b in zip(rows, rows[1:]):
        if 0.0 < a[0] < z <= b[0]:
            t = math.log(z / a[0]) / math.log(b[0] / a[0])
            return math.hypot(a[1], a[2]) + t * (math.hypot(b[1], b[2]) - math.hypot(a[1], a[2]))
    fail(f"{INFLOW}: no rows around z = {z} m")


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
    elif args[0] == "two_way" and set(args[1:]) <= {"paper", "relax", "inject", "relax_consistent", "keep"} \
            and len(set(args[1:])) == len(args) - 1:
        return two_way(cmd, paper="paper" in args, extra=[a for a in args if a in ("relax", "inject", "relax_consistent")],
                       keep="keep" in args)
    elif args[0] == "study" and len(args) == 1:
        return study(cmd)
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
