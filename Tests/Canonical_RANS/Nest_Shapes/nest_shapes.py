#!/usr/bin/env python3
"""One-way patches around canonical shapes against uniform fine runs (tracker row PR-76).

    nest_shapes.py files   SHAPE
    nest_shapes.py run     SHAPE NAME DX         -- COMMAND [ARGS ...]
    nest_shapes.py nest    SHAPE NAME PARENT D   -- COMMAND [ARGS ...]
    nest_shapes.py compare SHAPE PATCH REFERENCE [PARENT]
    nest_shapes.py study   [--shapes a,b,...] [--ranks a,b,...] -- MPIEXEC ... GUST [-- ARGS ...]

SHAPE is ridge, hill, escarpment (terrain, inputs_terrain_*), cube or plate
(buildings, inputs_building_*); SHAPES below holds each one's domain, cells,
parent cell size, distances and regions.

files: write the shape's input files (escarpment terrain, building map,
plate z-levels) into the working directory. run and nest write them too.

run: the whole domain on DX m cells to the steady verdict; the last plotfile
and checkpoint are copied to NAME_plt_final and NAME_chk_final.

nest: a patch at 2 2 1 of the finished run PARENT, its nest faces D m from
the shape's centre (a square of half-width D, or across the whole periodic
width where the shape is 2-D), to the steady verdict; outputs as for run.

COMMAND is the solver with its deck, e.g. mpiexec -n 2 gust
inputs_terrain_nest; the shape's inputs, the mesh, the parent and the output
names are appended (mesh.max_grid_size too, unless COMMAND gives it).

compare: PATCH against the whole-domain run REFERENCE on the same cells. For
the regions below it prints the RMS and the largest difference of the
velocity vector as a fraction of the shape's U_ref, and those of k as a
fraction of U_ref^2; with PARENT, the same for the parent against the
reference averaged over each parent cell.

Regions (cell centres; layers by the reference's corner column, over flat
ground):
  roi    around the shape (SHAPES[...]["roi"]), layers centred below ztop
  lee    downwind of the shape, the same layers (patches that hold it)
  faces  the two rings of patch cells beside the nest faces, all layers
  patch  every patch cell

study: for each shape the reference and parent runs and a patch at each
distance; with --ranks, the middle distance again in four boxes on 1, 2 and 4
ranks (the number after -n or -np of the launch command is replaced). The
targets of the README are judged; the table is printed and written to
nest_shapes.csv. The exit code is the verdict.
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

FIELDS = ["u", "v", "w", "tke", "z"]
KAPPA_REF = 0.4   # the speed scale (u* / 0.4) ln(h / z0); the solver's kappa is 0.41


def u_ref(height, z0):
    return 0.4 / KAPPA_REF * math.log(height / z0)


# Plate z-levels: 0.025 m layers to 1.5 m, then each layer 1.15 times the last up to at least 20 m.
def plate_zlevels():
    z = [0.025 * k for k in range(61)]
    dz = 0.025
    while z[-1] < 20.0:
        dz *= 1.15
        z.append(z[-1] + dz)
    return [round(v, 9) for v in z]   # the file holds these values, and prob_hi z is the last of them


ESC_H, ESC_L = 100.0, 100.0   # escarpment rise and tanh length (m)


def write_files(shape):
    s = SHAPES[shape]
    if shape == "escarpment":
        xs = [-100.0 + 10.0 * i for i in range(int((s["lx"] + 200.0) / 10.0) + 1)]
        ys = [-100.0, s["ly"] + 100.0]
        with open("escarpment_terrain.txt", "w") as f:
            f.write(f"{len(xs)}\n{len(ys)}\n")
            f.writelines(f"{x:.6f}\n" for x in xs)
            f.writelines(f"{y:.6f}\n" for y in ys)
            for x in xs:
                z = 0.5 * ESC_H * (1.0 + math.tanh((x - s["xc"]) / ESC_L))
                f.write(f"{z:.9f}\n{z:.9f}\n")
    elif shape == "cube":
        # ERF format: points own tiles halfway to their neighbours; the 20 m tiles are [115, 145] x [105, 135].
        with open("cube_nest.txt", "w") as f:
            f.write("4 4\n110. 120. 140. 150.\n100. 110. 130. 140.\n"
                    "0. 0. 0. 0.\n0. 20. 20. 0.\n0. 20. 20. 0.\n0. 0. 0. 0.\n")
    elif shape == "plate":
        # xyz with zstart: the plate x in [20, 24] m, z in [1.0, 1.1] m, across the 4 m periodic width.
        with open("plate_nest.txt", "w") as f:
            for x, h, zs in ((19.0, 0.0, 0.0), (21.0, 1.1, 1.0), (23.0, 1.1, 1.0), (25.0, 0.0, 0.0)):
                for y in (1.0, 3.0):
                    f.write(f"{x:g} {y:g} {h:g} {zs:g}\n")
        with open("plate_zlevels.txt", "w") as f:
            f.writelines(f"{z!r}\n" for z in plate_zlevels())


def agnesi(shape_kind):
    return ["terrain.type=agnesi", "terrain.hill_height=150.", "terrain.hill_half_width=150.",
            f"terrain.hill_shape={shape_kind}"]


# Each shape: domain lx x ly (periodic in y), the vertical cells, the parent cell size (the reference and the patches
# are on half of it), the centre and the distances D of the nest faces, the inputs given to every run, U_ref, the
# regions and whether the accuracy targets apply.
SHAPES = {
    "ridge": dict(kind="terrain", lx=3200.0, ly=200.0, z=["mesh.n_cell_z=16"], top=800.0, pdx=50.0,
                  xc=1600.0, yc=100.0, two_d=True, ds=(300.0, 600.0, 900.0),
                  args=agnesi("ridge_x") + ["convection.hybrid_height=375."],
                  uref=u_ref(150.0, 0.1), ztop=150.0,
                  roi=lambda x, y: abs(x - 1600.0) <= 150.0,
                  lee=lambda x, y: 1750.0 < x <= 2050.0, targets=True),
    "hill": dict(kind="terrain", lx=3200.0, ly=2400.0, z=["mesh.n_cell_z=16"], top=800.0, pdx=50.0,
                 xc=1600.0, yc=1200.0, two_d=False, ds=(300.0, 600.0, 900.0),
                 args=agnesi("round") + ["convection.hybrid_height=375."],
                 uref=u_ref(150.0, 0.1), ztop=150.0,
                 roi=lambda x, y: math.hypot(x - 1600.0, y - 1200.0) <= 150.0,
                 lee=lambda x, y: 1750.0 < x <= 2050.0 and abs(y - 1200.0) <= 150.0, targets=True),
    "escarpment": dict(kind="terrain", lx=3200.0, ly=200.0, z=["mesh.n_cell_z=16"], top=800.0, pdx=50.0,
                       xc=1600.0, yc=100.0, two_d=True, ds=(300.0, 600.0, 900.0),
                       args=["terrain.type=erf_grid", "terrain.file=escarpment_terrain.txt",
                             "convection.hybrid_height=250."],
                       uref=u_ref(100.0, 0.1), ztop=100.0,
                       roi=lambda x, y: abs(x - 1600.0) <= 200.0,
                       lee=lambda x, y: 1800.0 < x <= 2200.0, targets=True),
    "cube": dict(kind="building", lx=480.0, ly=240.0, z=["mesh.n_cell_z=24"], top=120.0, pdx=10.0,
                 xc=130.0, yc=120.0, two_d=False, ds=(40.0, 70.0, 100.0),
                 args=["buildings.file=cube_nest.txt", "buildings.format=erf", "bc.surface_z0=0.1",
                       "buildings.z0=0.1", "turb.max_geom_lscale=40.", "convection.hybrid_height=50."],
                 domain_args=["bc.inflow_z0=0.1"],
                 uref=u_ref(20.0, 0.1), ztop=40.0,
                 roi=lambda x, y: 100.0 <= x <= 160.0 and 90.0 <= y <= 150.0,
                 lee=lambda x, y: 145.0 < x <= 165.0 and 105.0 <= y <= 135.0, targets=True),
    "plate": dict(kind="building", lx=64.0, ly=4.0, z=["zlevels"], top=None, pdx=1.0,
                  xc=22.0, yc=2.0, two_d=True, ds=(4.0, 8.0, 12.0),
                  args=["buildings.file=plate_nest.txt", "buildings.format=xyz", "bc.surface_z0=0.01",
                        "buildings.z0=0.01", "turb.max_geom_lscale=40.", "convection.hybrid_height=2.75",
                        "mesh.zlevels_file=plate_zlevels.txt"],
                  domain_args=["bc.inflow_z0=0.01"],
                  uref=u_ref(1.1, 0.01), ztop=3.0,
                  roi=lambda x, y: 19.0 <= x <= 25.0,
                  lee=lambda x, y: 24.0 < x <= 28.0, targets=False),
}


def fail(msg):
    raise SystemExit("nest_shapes: " + msg)


def run_solver(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=14400)
    if proc.returncode != 0:
        print(proc.stdout[-6000:])
        fail(f"the run exited {proc.returncode}")
    return proc.stdout


def finish(name, cmd):
    """Run cmd with its outputs named; keep the last plotfile and checkpoint and the log."""
    for d in (glob.glob(name + "_plt[0-9]*") + glob.glob(name + "_chk[0-9]*")
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
        shutil.copytree(dirs[-1], f"{name}_{kind}_final")
        for d in dirs:
            shutil.rmtree(d)
    print(f"nest_shapes: {name}: steady after {m.group(1)} iterations", flush=True)
    return int(m.group(1))


def cells(length, dx):
    n = length / dx
    if not dx > 0.0 or abs(n - round(n)) > 1e-9:
        fail(f"{length} m is not a whole number of {dx} m cells")
    return int(round(n))


def boxes(cmd, default):
    """mesh.max_grid_size, unless the command gives it (a key given twice aborts the solver)."""
    return [] if any(a.startswith("mesh.max_grid_size") for a in cmd) else [f"mesh.max_grid_size={default}"]


def vertical(s):
    """n_cell z and prob_hi z of the shape."""
    if s["z"] == ["zlevels"]:
        z = plate_zlevels()
        return len(z) - 1, z[-1]
    return int(s["z"][0].split("=")[1]), s["top"]


def run_domain(shape, name, dx, cmd):
    s = SHAPES[shape]
    write_files(shape)
    nz, top = vertical(s)
    return finish(name, cmd + s["args"] + s.get("domain_args", [])
                  + [f"mesh.n_cell={cells(s['lx'], dx)} {cells(s['ly'], dx)} {nz}",
                     f"mesh.prob_hi={s['lx']} {s['ly']} {top!r}"] + boxes(cmd, 64))


def geometry(name):
    """(dx, x0, y0, nx, ny) of a finished run, from its plotfile header."""
    hdr = erf_plotfile._read_header(name + "_plt_final")
    nx, ny = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(2)]
    return hdr["dx"][0], hdr["prob_lo"][0], hdr["prob_lo"][1], nx, ny


def run_nest(shape, name, parent, d, cmd, extra=()):
    s = SHAPES[shape]
    write_files(shape)
    nz, top = vertical(s)
    pdx, px0, py0, pnx, pny = geometry(parent)
    dx = 0.5 * pdx
    xlo, xhi = s["xc"] - d, s["xc"] + d
    if s["two_d"]:
        ylo, yhi = py0, py0 + pny * pdx
        faces = ["mesh.is_periodic=0 1 0"]
    else:
        ylo, yhi = s["yc"] - d, s["yc"] + d
        faces = ["mesh.is_periodic=0 0 0", "bc.ylo=nest", "bc.yhi=nest"]
    for v, p0 in ((xlo, px0), (ylo, py0)):
        if abs((v - p0) / pdx - round((v - p0) / pdx)) > 1e-9:
            fail(f"a patch at D = {d} m does not lie on the {pdx} m cells of {parent}")
    if xlo < px0 or ylo < py0 or xhi > px0 + pnx * pdx or yhi > py0 + pny * pdx:
        fail(f"a patch at D = {d} m does not fit in {parent}")
    extra = list(extra) + boxes(list(cmd) + list(extra), 48)
    return finish(name, cmd + s["args"] + faces
                  + [f"mesh.n_cell={cells(xhi - xlo, dx)} {cells(yhi - ylo, dx)} {nz}",
                     f"mesh.prob_lo={xlo!r} {ylo!r} 0.", f"mesh.prob_hi={xhi!r} {yhi!r} {top!r}",
                     f"bc.nest_parent={parent}_chk_final"] + extra)


class Run:
    def __init__(self, name):
        self.name = name
        plt = name + "_plt_final"
        if not os.path.isdir(plt):
            fail(f"{plt} is not here: {name} is not a finished run in this directory")
        hdr, self.f = erf_plotfile.read_fields(plt, FIELDS)
        self.dx = hdr["dx"][0]
        self.x0, self.y0 = hdr["prob_lo"][0], hdr["prob_lo"][1]
        self.nx, self.ny, self.nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]

    def xy(self, i, j):
        return self.x0 + (i + 0.5) * self.dx, self.y0 + (j + 0.5) * self.dx


def larger(mx, v):
    """max(mx, v), in which a NaN in either stays."""
    return math.nan if math.isnan(mx) or math.isnan(v) else max(mx, v)


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

    def out(self, uref):
        """RMS and max of the velocity difference / U_ref, RMS and max of the k difference / U_ref^2."""
        if self.n == 0:
            return None
        return (math.sqrt(self.s2 / self.n) / uref, self.mx / uref,
                math.sqrt(self.k2 / self.n) / uref ** 2, self.kmx / uref ** 2)


REGIONS = ("roi", "lee", "faces", "patch")


def low_layers(s, ref):
    """Layers centred below ztop over flat ground (the reference's corner column)."""
    return set(k for k in range(ref.nz) if ref.f["z"][0][0][k] <= s["ztop"])


def offset(a, ref, ratio):
    ox, oy = (a.x0 - ref.x0) / ref.dx, (a.y0 - ref.y0) / ref.dx
    if (abs(a.dx - ratio * ref.dx) > 1e-9 * ref.dx or a.nz != ref.nz
            or abs(ox - round(ox)) > 1e-6 or abs(oy - round(oy)) > 1e-6):
        fail(f"{a.name} does not lie on the cells of {ref.name} (ratio {ratio})")
    return int(round(ox)), int(round(oy))


def norms(shape, run, ref, patch, ratio):
    """run (the patch, ratio 1, or its parent, ratio 2) against ref, in the regions of patch."""
    s = SHAPES[shape]
    ox, oy = offset(run, ref, ratio)
    low = low_layers(s, ref)
    xlo, ylo = patch.x0, patch.y0
    xhi, yhi = patch.x0 + patch.nx * patch.dx, patch.y0 + patch.ny * patch.dx
    w = max(2.0 * patch.dx, run.dx)   # two rings of patch cells, or one parent cell when that is wider
    out = {r: Norm() for r in REGIONS}
    P, R = run.f, ref.f
    for i in range(run.nx):
        for j in range(run.ny):
            x, y = run.xy(i, j)
            if not (xlo < x < xhi and ylo < y < yhi):
                continue
            roi, lee = s["roi"](x, y), s["lee"](x, y)
            ring = x < xlo + w or x > xhi - w or (not s["two_d"] and (y < ylo + w or y > yhi - w))
            for k in range(run.nz):
                d = []
                for f in ("u", "v", "w", "tke"):
                    acc = 0.0
                    for a in range(ratio):
                        for b in range(ratio):
                            acc += R[f][ratio * i + ox + a][ratio * j + oy + b][k]
                    d.append(P[f][i][j][k] - acc / (ratio * ratio))
                out["patch"].add(*d)
                if ring:
                    out["faces"].add(*d)
                if k in low:
                    if roi:
                        out["roi"].add(*d)
                    if lee:
                        out["lee"].add(*d)
    return {r: out[r].out(s["uref"]) for r in REGIONS}


def show(shape, label, nm):
    for reg in REGIONS:
        o = nm[reg]
        if o is None:
            print(f"nest_shapes: {label}: {reg:5s}: no cells")
        else:
            print(f"nest_shapes: {label}: {reg:5s}: velocity RMS {100 * o[0]:.3f} %, max {100 * o[1]:.3f} % of "
                  f"{SHAPES[shape]['uref']:.3f} m/s; k RMS {100 * o[2]:.4f} %, max {100 * o[3]:.4f} % of U_ref^2")


def compare(args):
    if len(args) not in (3, 4) or args[0] not in SHAPES:
        print(__doc__)
        return 2
    shape = args[0]
    patch, ref = Run(args[1]), Run(args[2])
    show(shape, f"{patch.name} vs {ref.name}", norms(shape, patch, ref, patch, 1))
    if len(args) == 4:
        parent = Run(args[3])
        show(shape, f"{parent.name} vs {ref.name}", norms(shape, parent, ref, patch, 2))
    return 0


def study(shapes, ranks, solver):
    flags = []
    if "--" in solver:
        flags = solver[solver.index("--") + 1:]
        solver = solver[:solver.index("--")]
    rows, bad = [], []

    def judge(target, ok, text):
        print(f"nest_shapes: {target}: {'met' if ok else 'MISSED'}: {text}", flush=True)
        if not ok:
            bad.append(target)

    def deck(shape, which):
        kind = SHAPES[shape]["kind"]
        return solver + [os.path.join(HERE, f"inputs_{kind}_{which}")] + flags

    def row(shape, case, patch, ref, parent, its, its_ref):
        p, r, q = Run(patch), Run(ref), Run(parent)
        mine, theirs = norms(shape, p, r, p, 1), norms(shape, q, r, p, 2)
        show(shape, f"{patch} vs {ref}", mine)
        show(shape, f"{parent} vs {ref}", theirs)
        rows.append((shape, case, patch, p.nx * p.ny * p.nz, its, its_ref, mine, theirs))
        return mine["roi"][0], theirs["roi"][0]

    for shape in shapes:
        s = SHAPES[shape]
        fdx = 0.5 * s["pdx"]
        try:
            it_ref = run_domain(shape, f"{shape}_ref", fdx, deck(shape, "domain"))
            it_par = run_domain(shape, f"{shape}_par", s["pdx"], deck(shape, "domain"))
        except SystemExit as err:
            print(err)
            judge("T4", False, f"{shape}: a whole-domain run did not reach the steady verdict")
            continue
        print(f"nest_shapes: {shape}: whole domain {fdx:g} m {it_ref} iterations, {s['pdx']:g} m {it_par}")
        hs = []
        for d in s["ds"]:
            name = f"{shape}_d{d:g}"
            try:
                its = run_nest(shape, name, f"{shape}_par", d, deck(shape, "nest"))
            except SystemExit as err:
                print(err)
                judge("T4", False, f"{name} did not reach the steady verdict")
                continue
            h, ph = row(shape, f"{shape}, D = {d:g} m", name, f"{shape}_ref", f"{shape}_par", its, it_ref)
            hs.append((d, h, ph))
        if s["targets"] and hs:
            d, h, ph = hs[-1]
            judge("T1", d == s["ds"][-1] and h <= 5e-3 and h <= 0.25 * ph,
                  f"{shape}, D = {d:g} m: ROI RMS {100 * h:.3f} % of U_ref, the parent's {100 * ph:.3f} %")
            for (d0, h0, _), (d1, h1, _) in zip(hs, hs[1:]):
                judge("T2", h1 <= 1.05 * h0,
                      f"{shape}: ROI RMS {100 * h1:.3f} % at D = {d1:g} m, {100 * h0:.3f} % at D = {d0:g} m")
        if shape in ranks:
            at = [i for i, a in enumerate(solver) if a in ("-n", "-np")]
            if not at:
                fail("ranks: the launch command has no -n or -np to replace")
            d = s["ds"][1]
            n_half = cells(d, fdx)   # four boxes: half the patch's cells in x, and in y
            ny_half = n_half if not s["two_d"] else cells(0.5 * s["ly"], fdx)
            vals = []
            for np in (1, 2, 4):
                cmd = deck(shape, "nest")
                cmd[at[0] + 1] = str(np)
                name = f"{shape}_np{np}"
                mgs = f"mesh.max_grid_size={n_half} {ny_half}"
                try:
                    its = run_nest(shape, name, f"{shape}_par", d, cmd, [mgs])
                except SystemExit as err:
                    print(err)
                    judge("T3", False, f"{name} did not reach the steady verdict")
                    continue
                h, _ = row(shape, f"{shape}, D = {d:g} m, four boxes, {np} ranks", name, f"{shape}_ref",
                           f"{shape}_par", its, it_ref)
                vals.append(h)
            if len(vals) == 3:
                judge("T3", max(vals) - min(vals) <= 0.1 * min(vals),
                      f"{shape}: ROI RMS on 1, 2 and 4 ranks {', '.join(f'{100 * v:.4f}' for v in vals)} %")
                gap = norms(shape, Run(f"{shape}_np4"), Run(f"{shape}_np1"), Run(f"{shape}_np1"), 1)["patch"]
                print(f"nest_shapes: {shape}_np4 vs {shape}_np1: largest velocity difference {gap[1]:.3e} of U_ref")
    if not rows:
        fail("no patch was run")
    with open("nest_shapes.csv", "w") as f:
        f.write("shape,case,patch,cells,iterations,reference_iterations," + ",".join(
            f"{who}_{reg}_{q}" for who in ("patch", "parent") for reg in REGIONS
            for q in ("vel_rms", "vel_max", "k_rms", "k_max")) + "\n")
        for shape, case, patch, ncell, its, its_ref, mine, theirs in rows:
            vals = []
            for nm in (mine, theirs):
                for reg in REGIONS:
                    vals += ["" if nm[reg] is None else f"{v:.6e}" for v in (nm[reg] or (None,) * 4)]
            f.write(f"{shape},\"{case}\",{patch},{ncell},{its},{its_ref}," + ",".join(vals) + "\n")
    print("nest_shapes: wrote nest_shapes.csv")
    print("nest_shapes: " + ("FAIL: missed " + ", ".join(sorted(set(bad))) if bad else "PASS"))
    return 1 if bad else 0


def shape_list(text):
    out = [t for t in text.split(",") if t]
    if any(t not in SHAPES for t in out):
        fail(f"unknown shape in {text!r}; the shapes are {', '.join(SHAPES)}")
    return out


def main(argv):
    if not argv:
        print(__doc__)
        return 2
    if argv[0] == "files" and len(argv) == 2 and argv[1] in SHAPES:
        write_files(argv[1])
        return 0
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
    if args[0] == "run" and len(args) == 4 and args[1] in SHAPES:
        run_domain(args[1], args[2], float(args[3]), cmd)
    elif args[0] == "nest" and len(args) == 5 and args[1] in SHAPES:
        run_nest(args[1], args[2], args[3], float(args[4]), cmd)
    elif args[0] == "study":
        shapes, ranks = list(SHAPES), []
        it = iter(args[1:])
        for a in it:
            if a == "--shapes":
                shapes = shape_list(next(it, ""))
            elif a == "--ranks":
                ranks = shape_list(next(it, ""))
            else:
                print(__doc__)
                return 2
        return study(shapes, ranks, cmd)
    else:
        print(__doc__)
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
