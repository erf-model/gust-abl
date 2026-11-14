#!/usr/bin/env python3
"""One-way nesting validation on the hill: patches against a uniform fine run.

    nest_validation.py run    NAME DX            -- COMMAND [ARGS ...]
    nest_validation.py nest   NAME PARENT D      -- COMMAND [ARGS ...]
    nest_validation.py compare PATCH REFERENCE [PARENT] [--rms-hill X] [--max-hill X] [--ratio R]
                               [--below PATCH2 ...]
    nest_validation.py study  [--stages one,chain,zone,ranks] -- MPIEXEC ... GUST [-- ARGS ...]

run: the whole domain (inputs_hill_domain) on DX m cells (100, 50, 25 or
12.5) to the steady verdict; the last plotfile and checkpoint are copied to
NAME_plt_final and NAME_chk_final.

nest: a square patch of half-width D m centred on the hill
(inputs_hill_nest), at 2 2 1 of the run PARENT (a whole-domain run or another
patch, by its NAME), to the steady verdict; outputs as for run.

COMMAND is the solver with its deck, e.g. mpiexec -n 2 gust inputs_hill_nest;
the mesh, the parent and the output names are appended (and
mesh.max_grid_size = 64 for a whole domain or 48 for a patch, unless COMMAND
gives it).

compare: PATCH against the whole-domain run REFERENCE on the same cells
(names of finished runs in the working directory). Prints, for the regions
below, the RMS and the largest difference of the velocity vector as a
fraction of U_REF, and those of k as a fraction of U_REF^2. With PARENT (the
patch's parent, or for a chain the whole-domain run it starts from) it prints
the same for the parent against the reference, on the parent's cells (the
reference averaged over each parent cell), and the patch's hill RMS must be
below the parent's. --rms-hill and --max-hill are limits on the patch's hill
numbers; --ratio R (needs PARENT): the patch's hill RMS must be at most R
times the parent's; --below names patches with nest faces nearer the hill,
whose hill RMS must not be smaller than this patch's (to 5 % of it). The exit
code is the verdict.

Regions (cell centres; r is the horizontal distance from the hill's centre,
a = 250 m the hill's half-width):
  hill   r <= a, layers centred below 100 m over flat ground
  lee    the square of half-width a centred 2 a downwind of the centre
         (wind from 20 degrees), the same layers; patches of D >= 3 a hold it
  faces  the two rings of patch cells beside the nest faces, all layers (for
         a parent: its cells centred within that strip, or within one of
         its own cells of the faces when that is wider)
  patch  every patch cell

study: the runs of the README, in the working directory (ARGS are given to
every run after its deck), with the measured
table printed and written to nest_validation.csv. Stages: one (one nest: the
50 m parent, 25 m patches of D = 400, 600, 800 and 1000 m, the 25 m
reference), chain (100 m whole domain, patches of 50, 25 and 12.5 m, the
12.5 m reference), zone (D = 600 m with bc.nest_relax_width = 0, 3 and 5), ranks
(D = 600 m in boxes of 24 x 24 cells on 1, 2 and 4 ranks; the number after -n
or -np in the launch command is replaced). Targets T1 to T5 of the README are
judged for the stages run; a run that does not reach the steady verdict
misses its target. The exit code is the verdict.
"""
import glob
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402

LX, LY = 3200.0, 2400.0          # inputs_hill_domain
XC, YC = 0.5 * LX, 0.5 * LY      # the hill's centre
A = 250.0                        # terrain.hill_half_width
WIND = math.radians(20.0)        # bc.inflow_direction
Z_TOP = 100.0                    # the hill and lee regions end here (m above flat ground)
# The speed scale (u* / 0.4) ln(h / z0) with u* = 0.4 m/s, h = 150 m, z0 = 0.1 m: 7.313 m/s. The inflow speed at h is
# 7.14 m/s (the solver's kappa is 0.41); the scale is kept at the value the README's targets were set with.
U_REF = 0.4 / 0.4 * math.log(150.0 / 0.1)
FIELDS = ["u", "v", "w", "tke", "z"]


def fail(msg):
    raise SystemExit("nest_validation: " + msg)


def run_solver(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=14400)
    if proc.returncode != 0:
        print(proc.stdout)
        fail(f"the run exited {proc.returncode}")
    return proc.stdout


def finish(name, cmd):
    """Run cmd with its outputs named; keep the last plotfile and checkpoint and the log."""
    for d in (glob.glob(name + "_plt[0-9]*") + glob.glob(name + "_chk[0-9]*") + glob.glob(name + "_plt_level*")
              + [name + "_plt_final", name + "_chk_final"]):
        shutil.rmtree(d, ignore_errors=True)   # nothing of an earlier run of this name is left to be read
    log = run_solver(cmd + [f"io.plot_file={name}_plt", f"io.chk_file={name}_chk"])
    with open(name + ".log", "w") as f:
        f.write(log)
    m = re.search(r"^gust: steady after (\d+) iterations", log, re.M)
    if not m:
        print(log[-3000:])
        fail(f"{name}: the run did not reach the steady verdict")
    for kind in ("plt", "chk"):
        dirs = sorted(d for d in glob.glob(f"{name}_{kind}[0-9]*") if os.path.isdir(d))
        if not dirs:
            fail(f"no {name}_{kind}NNNNN was written")
        shutil.rmtree(f"{name}_{kind}_final", ignore_errors=True)
        shutil.copytree(dirs[-1], f"{name}_{kind}_final")
        for d in dirs:
            shutil.rmtree(d)
    print(f"nest_validation: {name}: steady after {m.group(1)} iterations")
    return int(m.group(1))


def cells(length, dx):
    if not dx > 0.0:
        fail(f"the cell size {dx} m is not positive")
    n = length / dx
    if abs(n - round(n)) > 1e-9:
        fail(f"{length} m is not a whole number of {dx} m cells")
    return int(round(n))


def boxes(cmd, default):
    """mesh.max_grid_size, unless the command gives it (a key given twice aborts the solver)."""
    return [] if any(a.startswith("mesh.max_grid_size") for a in cmd) else [f"mesh.max_grid_size={default}"]


def run_domain(name, dx, cmd):
    return finish(name, cmd + [f"mesh.n_cell={cells(LX, dx)} {cells(LY, dx)} 16"] + boxes(cmd, 64))


def geometry(name):
    """(dx, x0, y0, nx, ny) of a finished run, from its plotfile header."""
    hdr = erf_plotfile._read_header(name + "_plt_final")
    nx, ny = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(2)]
    return hdr["dx"][0], hdr["prob_lo"][0], hdr["prob_lo"][1], nx, ny


def run_nest(name, parent, d, cmd, extra=()):
    pdx, px0, py0, pnx, pny = geometry(parent)
    dx = 0.5 * pdx
    n = 2 * cells(d, dx)
    lo = (XC - d, YC - d)
    for v, p0 in zip(lo, (px0, py0)):
        if abs((v - p0) / pdx - round((v - p0) / pdx)) > 1e-9:
            fail(f"a patch of half-width {d} m does not lie on the {pdx} m cells of {parent}")
    if lo[0] < px0 or lo[1] < py0 or XC + d > px0 + pnx * pdx or YC + d > py0 + pny * pdx:
        fail(f"a patch of half-width {d} m does not fit in {parent}")
    extra = list(extra) + boxes(list(cmd) + list(extra), 48)
    return finish(name, cmd + [f"mesh.n_cell={n} {n} 16", f"mesh.prob_lo={lo[0]} {lo[1]} 0.",
                               f"mesh.prob_hi={XC + d} {YC + d} 800.", f"bc.nest_parent={parent}_chk_final"]
                  + list(extra))


class Run:
    def __init__(self, name):
        self.name = name
        plt = name + "_plt_final"
        if not os.path.isdir(plt):
            fail(f"{plt} is not here: {name} is not a finished run in this directory")
        hdr = erf_plotfile._read_header(plt)
        ncomp, _, bxs, fabs = erf_plotfile._read_cell_h(os.path.join(plt, "Level_0"))
        want = 1
        for d in range(3):
            want *= hdr["hi"][d] - hdr["lo"][d] + 1
        have = sum((hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1) for lo, hi in bxs)
        if have != want or len(fabs) != len(bxs):
            fail(f"{plt}: the boxes hold {have} cells, the domain has {want}")
        for (lo, hi), (fname, off) in zip(bxs, fabs):
            npts = (hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1)
            path = os.path.join(plt, "Level_0", fname)
            if not os.path.isfile(path) or os.path.getsize(path) < off + 8 * ncomp * npts:
                fail(f"{path} is missing or shorter than its box")
        hdr, self.f = erf_plotfile.read_fields(plt, FIELDS)
        self.dx = hdr["dx"][0]
        self.x0, self.y0 = hdr["prob_lo"][0], hdr["prob_lo"][1]
        self.nx, self.ny, self.nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]

    def xy(self, i, j):
        return self.x0 + (i + 0.5) * self.dx, self.y0 + (j + 0.5) * self.dx


def offset(patch, ref):
    """Index of the reference cell under patch cell (0, 0); the two have the same cells."""
    if abs(patch.dx - ref.dx) > 1e-9 * ref.dx or patch.nz != ref.nz:
        fail(f"{patch.name} and {ref.name} do not have the same cells")
    ox, oy = (patch.x0 - ref.x0) / ref.dx, (patch.y0 - ref.y0) / ref.dx
    if abs(ox - round(ox)) > 1e-6 or abs(oy - round(oy)) > 1e-6:
        fail(f"{patch.name} does not lie on the cells of {ref.name}")
    return int(round(ox)), int(round(oy))


def low_layers(ref):
    """Layers centred below Z_TOP over flat ground (the reference's corner column, 1700 m from the hill)."""
    return [k for k in range(ref.nz) if ref.f["z"][0][0][k] <= Z_TOP]


def in_hill(x, y):
    return math.hypot(x - XC, y - YC) <= A


def in_lee(x, y):
    cx, cy = XC + 2.0 * A * math.cos(WIND), YC + 2.0 * A * math.sin(WIND)
    return abs(x - cx) <= A and abs(y - cy) <= A


def larger(mx, v):
    """max(mx, v), in which a NaN in either stays."""
    if math.isnan(mx) or math.isnan(v):
        return math.nan
    return max(mx, v)


class Norm:
    def __init__(self):
        self.n, self.s2, self.mx, self.k2, self.kmx = 0, 0.0, 0.0, 0.0, 0.0

    def add(self, du, dv, dw, dk):
        e2 = du * du + dv * dv + dw * dw
        self.n += 1
        self.s2 += e2
        self.k2 += dk * dk
        self.mx = larger(self.mx, math.sqrt(e2))
        self.kmx = larger(self.kmx, abs(dk))

    def out(self):
        """RMS and max of the velocity difference / U_REF, RMS and max of the k difference / U_REF^2."""
        if self.n == 0:
            return None
        return (math.sqrt(self.s2 / self.n) / U_REF, self.mx / U_REF,
                math.sqrt(self.k2 / self.n) / U_REF ** 2, self.kmx / U_REF ** 2)


REGIONS = ("hill", "lee", "faces", "patch")


def patch_norms(patch, ref):
    ox, oy = offset(patch, ref)
    low = set(low_layers(ref))
    norms = {r: Norm() for r in REGIONS}
    P, R = patch.f, ref.f
    for i in range(patch.nx):
        ring_i = i < 2 or i >= patch.nx - 2
        for j in range(patch.ny):
            x, y = patch.xy(i, j)
            hill, lee = in_hill(x, y), in_lee(x, y)
            ring = ring_i or j < 2 or j >= patch.ny - 2
            for k in range(patch.nz):
                d = [P[f][i][j][k] - R[f][i + ox][j + oy][k] for f in ("u", "v", "w", "tke")]
                norms["patch"].add(*d)
                if ring:
                    norms["faces"].add(*d)
                if k in low:
                    if hill:
                        norms["hill"].add(*d)
                    if lee:
                        norms["lee"].add(*d)
    return {r: norms[r].out() for r in REGIONS}


def parent_norms(parent, ref, patch):
    """The parent against the reference averaged over each parent cell, in the regions of patch."""
    ratio = parent.dx / ref.dx
    r = int(round(ratio))
    if abs(ratio - r) > 1e-9 or parent.nz != ref.nz:
        fail(f"{parent.name} is not a coarsening of {ref.name}")
    ox, oy = (parent.x0 - ref.x0) / ref.dx, (parent.y0 - ref.y0) / ref.dx
    if abs(ox - round(ox)) > 1e-6 or abs(oy - round(oy)) > 1e-6:
        fail(f"{parent.name} does not lie on the cells of {ref.name}")
    ox, oy = int(round(ox)), int(round(oy))
    low = set(low_layers(ref))
    xlo, ylo = patch.x0, patch.y0
    xhi, yhi = patch.x0 + patch.nx * patch.dx, patch.y0 + patch.ny * patch.dx
    w = max(2.0 * patch.dx, parent.dx)   # the two rings of patch cells, or one parent cell when that is wider
    norms = {reg: Norm() for reg in REGIONS}
    P, R = parent.f, ref.f
    for i in range(parent.nx):
        for j in range(parent.ny):
            x, y = parent.xy(i, j)
            if not (xlo < x < xhi and ylo < y < yhi):
                continue
            hill, lee = in_hill(x, y), in_lee(x, y)
            ring = x < xlo + w or x > xhi - w or y < ylo + w or y > yhi - w
            for k in range(parent.nz):
                d = []
                for f in ("u", "v", "w", "tke"):
                    s = 0.0
                    for a in range(r):
                        for b in range(r):
                            s += R[f][r * i + ox + a][r * j + oy + b][k]
                    d.append(P[f][i][j][k] - s / (r * r))
                norms["patch"].add(*d)
                if ring:
                    norms["faces"].add(*d)
                if k in low:
                    if hill:
                        norms["hill"].add(*d)
                    if lee:
                        norms["lee"].add(*d)
    return {reg: norms[reg].out() for reg in REGIONS}


def show(label, norms):
    for reg in REGIONS:
        o = norms[reg]
        if o is None:
            print(f"nest_validation: {label}: {reg:5s}: no cells")
        else:
            print(f"nest_validation: {label}: {reg:5s}: velocity RMS {100 * o[0]:.3f} %, max {100 * o[1]:.3f} % of "
                  f"{U_REF:.3f} m/s; k RMS {100 * o[2]:.4f} %, max {100 * o[3]:.4f} % of U_REF^2")


def compare(args):
    names, rms_hill, max_hill, ratio, below = [], None, None, None, []
    it = iter(args)
    try:
        for a in it:
            if a == "--rms-hill":
                rms_hill = float(next(it))
            elif a == "--max-hill":
                max_hill = float(next(it))
            elif a == "--ratio":
                ratio = float(next(it))
            elif a == "--below":
                below = list(it)
            else:
                names.append(a)
    except (StopIteration, ValueError):
        names = []
    if len(names) not in (2, 3) or (ratio is not None and len(names) != 3):
        print(__doc__)
        return 2
    patch, ref = Run(names[0]), Run(names[1])
    mine = patch_norms(patch, ref)
    show(f"{patch.name} vs {ref.name}", mine)
    bad = []
    hill = mine["hill"]
    if hill is None:
        bad.append("the patch has no cell in the hill region")
    else:
        if rms_hill is not None and not hill[0] <= rms_hill:
            bad.append(f"hill: velocity RMS {hill[0]:.4e} of U_REF is above {rms_hill:g}")
        if max_hill is not None and not hill[1] <= max_hill:
            bad.append(f"hill: largest velocity difference {hill[1]:.4e} of U_REF is above {max_hill:g}")
    if len(names) == 3:
        parent = Run(names[2])
        theirs = parent_norms(parent, ref, patch)
        show(f"{parent.name} vs {ref.name}", theirs)
        if hill is not None and theirs["hill"] is not None and not hill[0] < theirs["hill"][0]:
            bad.append(f"hill: the patch's velocity RMS {hill[0]:.4e} is not below its parent's {theirs['hill'][0]:.4e}")
        if hill is not None and ratio is not None and not (theirs["hill"] is not None
                                                           and hill[0] <= ratio * theirs["hill"][0]):
            bad.append(f"hill: the patch's velocity RMS {hill[0]:.4e} is above {ratio:g} x its parent's")
    for other in below:
        o = patch_norms(Run(other), ref)["hill"]
        if o is None:
            bad.append(f"{other} has no cell in the hill region")
            continue
        print(f"nest_validation: {other} vs {ref.name}: hill : velocity RMS {100 * o[0]:.3f} %")
        if hill is not None and not hill[0] <= 1.05 * o[0]:
            bad.append(f"hill: velocity RMS {hill[0]:.4e} is above that of {other} ({o[0]:.4e}), whose nest faces "
                       "are nearer the hill")
    for b in bad:
        print("nest_validation: FAIL: " + b)
    if not bad:
        print("nest_validation: PASS")
    return 1 if bad else 0


def study(stages, solver):
    flags = []
    if "--" in solver:
        flags = solver[solver.index("--") + 1:]
        solver = solver[:solver.index("--")]
    dom = solver + [os.path.join(HERE, "inputs_hill_domain")] + flags
    nst = solver + [os.path.join(HERE, "inputs_hill_nest")] + flags
    rows, bad = [], []

    def judge(target, ok, text):
        print(f"nest_validation: {target}: {'met' if ok else 'MISSED'}: {text}")
        if not ok:
            bad.append(target)

    def row(case, patch, ref, parent, its):
        p, r, q = Run(patch), Run(ref), Run(parent)
        mine, theirs = patch_norms(p, r), parent_norms(q, r, p)
        show(f"{patch} vs {ref}", mine)
        show(f"{parent} vs {ref}", theirs)
        rows.append((case, patch, p.nx * p.ny * p.nz, its, mine, theirs))
        return mine["hill"][0], theirs["hill"][0]

    if "one" in stages or "zone" in stages or "ranks" in stages:
        it_ref = run_domain("ref25", 25.0, dom)
        it_par = run_domain("par50", 50.0, dom)
        print(f"nest_validation: whole domain: 25 m {it_ref} iterations, 50 m {it_par} iterations")
    if "one" in stages:
        last = None
        for d in (400, 600, 800, 1000):
            its = run_nest(f"one_d{d}", "par50", float(d), nst)
            h, ph = row(f"one nest, D = {d} m", f"one_d{d}", "ref25", "par50", its)
            if d >= 600:
                judge("T1", h <= 1e-3 and h <= 0.25 * ph,
                      f"D = {d} m: hill RMS {100 * h:.3f} % of U_ref, the parent's {100 * ph:.3f} %")
            if last is not None:
                judge("T2", h <= 1.05 * last, f"D = {d} m: hill RMS {100 * h:.3f} %, {100 * last:.3f} % one size down")
            last = h
    if "zone" in stages:
        for n in (0, 3, 5):
            its = run_nest(f"zone_n{n}", "par50", 600.0, nst, [f"bc.nest_relax_width={n}"] if n else [])
            row(f"D = 600 m, zone of {n} cells", f"zone_n{n}", "ref25", "par50", its)
    if "chain" in stages:
        it_ref = run_domain("ref12", 12.5, dom)
        it_par = run_domain("par100", 100.0, dom)
        print(f"nest_validation: whole domain: 12.5 m {it_ref} iterations, 100 m {it_par} iterations")
        try:
            i1 = run_nest("chain_50", "par100", 1000.0, nst)
            i2 = run_nest("chain_25", "chain_50", 800.0, nst)
            i3 = run_nest("chain_12", "chain_25", 600.0, nst)
        except SystemExit as err:
            print(err)
            judge("T3", False, "the chain was not finished")
            judge("T4", False, "a patch of the chain did not reach the steady verdict")
        else:
            h, ph = row(f"chain of three (iterations {i1} + {i2} + {i3}), D = 600 m", "chain_12", "ref12", "par100",
                        i3)
            judge("T3", h <= 0.5 * ph and h <= 5e-3,
                  f"chain: hill RMS {100 * h:.3f} % of U_ref, the 100 m run's {100 * ph:.3f} %")
            judge("T4", True, f"the patches of the chain reached the steady verdict ({i1}, {i2}, {i3} iterations)")
    if "ranks" in stages:
        at = [i for i, a in enumerate(solver) if a in ("-n", "-np")]
        if not at:
            fail("ranks: the launch command has no -n or -np to replace")
        hs, ph = [], None
        for np in (1, 2, 4):
            cmd = list(nst)
            cmd[at[0] + 1] = str(np)
            try:
                its = run_nest(f"rank_np{np}", "par50", 600.0, cmd, ["mesh.max_grid_size=24"])
            except SystemExit as err:
                print(err)
                judge("T5", False, f"{np} ranks: the patch did not reach the steady verdict")
                continue
            h, ph = row(f"D = 600 m, four boxes, {np} ranks", f"rank_np{np}", "ref25", "par50", its)
            hs.append(h)
            judge("T5", h <= 1e-3 and h <= 0.25 * ph, f"{np} ranks: hill RMS {100 * h:.3f} % of U_ref, {its} iterations")
        if len(hs) == 3:
            judge("T5", max(hs) - min(hs) <= 0.1 * min(hs),
                  f"hill RMS on 1, 2 and 4 ranks: {', '.join(f'{100 * h:.4f}' for h in hs)} %")
            a, b = Run("rank_np1"), Run("rank_np4")
            gap = patch_norms(b, a)["patch"]
            print(f"nest_validation: rank_np4 vs rank_np1: largest velocity difference {gap[1]:.3e} of U_ref")
    if not rows:
        fail("no stage was run")
    with open("nest_validation.csv", "w") as f:
        f.write("case,patch,cells,iterations," + ",".join(
            f"{who}_{reg}_{q}" for who in ("patch", "parent") for reg in REGIONS
            for q in ("vel_rms", "vel_max", "k_rms", "k_max")) + "\n")
        for case, patch, ncell, its, mine, theirs in rows:
            vals = []
            for norms in (mine, theirs):
                for reg in REGIONS:
                    vals += ["" if norms[reg] is None else f"{v:.6e}" for v in (norms[reg] or (None,) * 4)]
            f.write(f"\"{case}\",{patch},{ncell},{its}," + ",".join(vals) + "\n")
    print("nest_validation: wrote nest_validation.csv")
    print("nest_validation: " + ("FAIL: missed " + ", ".join(sorted(set(bad))) if bad else "PASS"))
    return 1 if bad else 0


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
        run_domain(args[1], float(args[2]), cmd)
    elif args[0] == "nest" and len(args) == 4:
        run_nest(args[1], args[2], float(args[3]), cmd)
    elif args[0] == "study":
        stages = ["one", "chain", "zone", "ranks"]
        if len(args) == 3 and args[1] == "--stages":
            stages = args[2].split(",")
        elif len(args) != 1:
            print(__doc__)
            return 2
        if not stages or any(s not in ("one", "chain", "zone", "ranks") for s in stages):
            print(__doc__)
            return 2
        return study(stages, cmd)
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
