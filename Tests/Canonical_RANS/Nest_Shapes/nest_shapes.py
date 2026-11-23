#!/usr/bin/env python3
"""One-way patches around canonical shapes against uniform fine runs.

    nest_shapes.py files   SHAPE
    nest_shapes.py run     SHAPE NAME DX         -- COMMAND [ARGS ...]
    nest_shapes.py nest    SHAPE NAME PARENT D [--nest-ratio R] -- COMMAND [ARGS ...]
    nest_shapes.py compare SHAPE PATCH REFERENCE [PARENT] [--rms-roi X] [--ratio R] [--max-patch X]
    nest_shapes.py study   [--shapes a,b,...] [--ranks a,b,...] -- MPIEXEC ... GUST [-- ARGS ...]

SHAPE is ridge, hill, escarpment (terrain, inputs_terrain_*), cube, plate or
solar (buildings, inputs_building_*); SHAPES below holds each one's domain,
cells, parent cell size, distances and regions.

files: write the shape's input files (escarpment terrain, building map,
plate and solar z-levels) into the working directory. run and nest write
them too.

run: the whole domain on DX m cells to the steady verdict; the last plotfile
and checkpoint are copied to NAME_plt_final and NAME_chk_final. The solver's
output goes to NAME.log as it runs; there is no time limit (steady.max_iterations
ends a run that does not converge).

nest: a patch at R R 1 (default 2 2 1) of the finished run PARENT, its nest faces D m from
the shape's centre (a square of half-width D, or across the whole periodic
width where the shape is 2-D), to the steady verdict; outputs as for run.
D may be DX,DY (a rectangle of half-widths DX and DY). PARENT may itself be a
patch (a chain of nests). --nest-ratio R: refinement R R 1 (bc.nest_ratio =
R; 2, 4, 8 or 16), the patch's cells 1/R of the parent's.

COMMAND is the solver with its deck, e.g. mpiexec -n 2 gust
inputs_terrain_nest; the shape's inputs, the mesh, the parent and the output
names are appended (mesh.max_grid_size too, unless COMMAND gives it).

compare: PATCH against the whole-domain run REFERENCE on the same cells. For
the regions below it prints the RMS and the largest difference of the
velocity vector as a fraction of the shape's U_ref, and those of k as a
fraction of U_ref^2; with PARENT, the same for the parent against the
reference averaged over each parent cell. --rms-roi: a limit on the patch's
ROI velocity RMS (a fraction of U_ref); --ratio R (needs PARENT): the patch's
ROI RMS must be at most R times the parent's; --max-patch X: a limit on the
largest velocity difference over the whole patch (a fraction of U_ref). The
cells of PATCH and of PARENT are each a whole number of REFERENCE's (one for
a patch on the reference's cells; REFERENCE need only cover PATCH, so a patch
can be compared with a patch), and REFERENCE is averaged over each. The exit code is the verdict.

Regions (cell centres; layers by the reference's corner column, over flat
ground):
  roi    around the shape (SHAPES[...]["roi"]), layers centred below ztop
  lee    downwind of the shape, the same layers (patches that hold it)
  faces  the two rings of patch cells beside the nest faces, all layers
  patch  every patch cell

study: for each shape the reference and parent runs and a patch at each
distance; with --ranks, the middle distance again in four boxes on 1, 2 and 4
ranks (the number after -n or -np of the launch command is replaced). The
solar farm is a chain: the 1 m parent, a 0.5 m middle patch at D = 64 m, and
0.25 m patches of the middle one at each distance, judged against the middle
patch; then, for the record, 0.25 m patches of the 1 m parent at ratio 4 at
each distance. Without --shapes the study runs every shape but the solar farm,
whose 0.25 m reference takes hours; name it with --shapes solar. The targets
of the README are judged; the table is printed and written to
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
RANK_TOL = 1e-8   # T3: the largest velocity difference between rank counts, as a fraction of U_ref
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


# Solar farm z-levels: 0.02 m layers to 3.2 m (twice the panel top), then each layer 1.1 times the last up to at
# least 60 m. The panel (z 1.5 to 1.6 m) is five layers thick.
def solar_zlevels():
    z = [0.02 * k for k in range(161)]
    dz = 0.02
    while z[-1] < 60.0:
        dz *= 1.1
        z.append(z[-1] + dz)
    return [round(v, 9) for v in z]


# Solar farm: SOLAR_ROWS flat panels SOLAR_CHORD m long in x, at a pitch of SOLAR_PITCH m, the first centred at
# SOLAR_X0 m; each panel spans z in [SOLAR_BOTTOM, SOLAR_TOP] m, across the periodic width (rows without ends).
SOLAR_ROWS, SOLAR_CHORD, SOLAR_PITCH, SOLAR_X0 = 6, 4.0, 10.0, 55.0
SOLAR_BOTTOM, SOLAR_TOP = 1.5, 1.6
SOLAR_END = SOLAR_X0 + (SOLAR_ROWS - 1) * SOLAR_PITCH + 0.5 * SOLAR_CHORD   # the last row's trailing edge (107 m)


def solar_panel(x):
    """Whether the 2 m tile centred at x lies in a panel."""
    return any(abs(x - (SOLAR_X0 + r * SOLAR_PITCH)) < 0.5 * SOLAR_CHORD for r in range(SOLAR_ROWS))


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
    elif shape == "solar":
        # xyz with zstart: points every 2 m own the tiles [x - 1, x + 1]; a panel is two tiles.
        with open("solar_nest.txt", "w") as f:
            for i in range(31):
                x = 50.0 + 2.0 * i
                h, zs = (SOLAR_TOP, SOLAR_BOTTOM) if solar_panel(x) else (0.0, 0.0)
                for y in (0.5, 1.5):
                    f.write(f"{x:g} {y:g} {h:g} {zs:g}\n")
        with open("solar_zlevels.txt", "w") as f:
            f.writelines(f"{z!r}\n" for z in solar_zlevels())


def agnesi(shape_kind):
    return ["terrain.type=agnesi", "terrain.hill_height=150.", "terrain.hill_half_width=150.",
            f"terrain.hill_shape={shape_kind}"]


# Each shape: domain lx x ly (periodic in y), the vertical cells (n_cell_z, or a z-levels writer `zlevels`), the
# parent cell size pdx, the centre and the distances D of the nest faces, the inputs given to every run, U_ref, the
# regions and whether the accuracy targets apply. The reference and the patches are on pdx / 2**levels: levels = 1
# (default) is one nest at 2 2 1; levels = 2 is a chain through a middle patch of half-width mid. direct = 2**levels
# adds patches straight from the parent at that ratio (no target). manual = True leaves the shape out of a study
# without --shapes.
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
                       "buildings.z0=0.1", "turb.max_geom_lscale=40.", "convection.hybrid_height=50.",
                       "steady.cfl_max=8", "steady.max_iterations=4000"],
                 domain_args=["bc.inflow_z0=0.1"],
                 uref=u_ref(20.0, 0.1), ztop=40.0,
                 roi=lambda x, y: 100.0 <= x <= 160.0 and 90.0 <= y <= 150.0,
                 lee=lambda x, y: 145.0 < x <= 165.0 and 105.0 <= y <= 135.0, targets=True),
    "plate": dict(kind="building", lx=64.0, ly=4.0, z=["zlevels"], zlevels=plate_zlevels, top=None, pdx=1.0,
                  xc=22.0, yc=2.0, two_d=True, ds=(4.0, 8.0, 12.0),
                  args=["buildings.file=plate_nest.txt", "buildings.format=xyz", "bc.surface_z0=0.01",
                        "buildings.z0=0.01", "turb.max_geom_lscale=40.", "convection.hybrid_height=2.75",
                        "mesh.zlevels_file=plate_zlevels.txt", "steady.cfl_max=4", "steady.max_iterations=20000",
                        "steady.confidence=0.95"],
                  domain_args=["bc.inflow_z0=0.01"],
                  uref=u_ref(1.1, 0.01), ztop=3.0,
                  roi=lambda x, y: 19.0 <= x <= 25.0,
                  lee=lambda x, y: 24.0 < x <= 28.0, targets=False),
    # Six rows of flat panels on 0.02 m layers; D from the array's centre (80 m), the array 27 m either side.
    # A chain of two nests: 1 m parent, 0.5 m middle patch (D = mid), 0.25 m patches (D = ds) of the middle one;
    # direct=4: 0.25 m patches of the 1 m parent at ratio 4 too (for the record, no target). ly = 2 m gives the
    # 1 m parent two cells across the periodic width, as many as the state's ghost cells. The pressure solve to
    # 1e-6: at the default 1e-3 its error shows as jumps in the k residual through the ground's k rows, whose
    # weight u*/d0 is about 35 1/s on the 0.02 m first layer (fields within 4e-8 of U_ref of the 1e-6 run; README).
    "solar": dict(kind="building", lx=160.0, ly=2.0, z=["zlevels"], zlevels=solar_zlevels, top=None, pdx=1.0,
                  xc=80.0, yc=1.0, two_d=True, ds=(32.0, 40.0, 48.0), levels=2, mid=64.0, direct=4,
                  manual=True,
                  args=["buildings.file=solar_nest.txt", "buildings.format=xyz", "bc.surface_z0=0.01",
                        "buildings.z0=0.01", "turb.max_geom_lscale=40.", "convection.hybrid_height=4.",
                        "mesh.zlevels_file=solar_zlevels.txt", "steady.cfl_max=4", "steady.max_iterations=20000",
                        "steady.confidence=0.95", "simplec.pressure_reltol=1e-6"],
                  domain_args=["bc.inflow_z0=0.01"],
                  uref=u_ref(SOLAR_TOP, 0.01), ztop=3.2,
                  roi=lambda x, y: SOLAR_X0 - 0.5 * SOLAR_CHORD <= x <= SOLAR_END,
                  lee=lambda x, y: SOLAR_END < x <= SOLAR_END + 4.0, targets=True),
}


def fail(msg):
    raise SystemExit("nest_shapes: " + msg)


def run_solver(cmd, log_name):
    """Run cmd with its output written to log_name as it runs (a long run can be followed); return the output."""
    with open(log_name, "w") as f:
        proc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
    with open(log_name) as f:
        log = f.read()
    if proc.returncode != 0:
        print(log[-6000:])
        fail(f"the run exited {proc.returncode}; its output is in {log_name}")
    return log


def finish(name, cmd):
    """Run cmd with its outputs named; keep the last plotfile and checkpoint and the log."""
    for d in (glob.glob(name + "_plt[0-9]*") + glob.glob(name + "_chk[0-9]*")
              + [name + "_plt_final", name + "_chk_final"]):
        shutil.rmtree(d, ignore_errors=True)   # nothing of an earlier run of this name is left to be read
    log = run_solver(cmd + [f"io.plot_file={name}_plt", f"io.chk_file={name}_chk"], name + ".log")
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
        z = s["zlevels"]()
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


def half_widths(text):
    """D or DX,DY (m) -> (DX, DY)."""
    try:
        v = [float(t) for t in text.split(",")]
    except ValueError:
        v = []
    if len(v) not in (1, 2) or not all(x > 0.0 for x in v):
        fail(f"D = {text!r}: one or two half-widths > 0 (m)")
    return (v[0], v[-1])


def run_nest(shape, name, parent, d, cmd, extra=(), ratio=2):
    """d: the half-width, or (DX, DY)."""
    s = SHAPES[shape]
    write_files(shape)
    nz, top = vertical(s)
    pdx, px0, py0, pnx, pny = geometry(parent)
    dx = pdx / ratio
    dx_, dy_ = d if isinstance(d, tuple) else (d, d)
    xlo, xhi = s["xc"] - dx_, s["xc"] + dx_
    if s["two_d"]:
        ylo, yhi = py0, py0 + pny * pdx
        faces = ["mesh.is_periodic=0 1 0"]
    else:
        ylo, yhi = s["yc"] - dy_, s["yc"] + dy_
        faces = ["mesh.is_periodic=0 0 0", "bc.ylo=nest", "bc.yhi=nest"]
    for v, p0 in ((xlo, px0), (ylo, py0)):
        if abs((v - p0) / pdx - round((v - p0) / pdx)) > 1e-9:
            fail(f"a patch at D = {d} m does not lie on the {pdx} m cells of {parent}")
    if xlo < px0 or ylo < py0 or xhi > px0 + pnx * pdx or yhi > py0 + pny * pdx:
        fail(f"a patch at D = {d} m does not fit in {parent}")
    if ratio != 2:   # the default adds no key: the ratio-2 commands stay as they were
        extra = list(extra) + [f"bc.nest_ratio={ratio}"]
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


def cell_ratio(a, ref):
    """The whole number of ref's cells in one of a's (in x and y)."""
    r = a.dx / ref.dx
    if abs(r - round(r)) > 1e-9 * r or round(r) < 1:
        fail(f"the cells of {a.name} ({a.dx} m) are not a whole number of {ref.name}'s ({ref.dx} m)")
    return int(round(r))


def norms(shape, run, ref, patch, ratio):
    """run (the patch, ratio 1, or its parent, ratio r) against ref, in the regions of patch."""
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
    names, rms_roi, ratio, max_patch = [], None, None, None
    it = iter(args)
    try:
        for a in it:
            if a == "--rms-roi":
                rms_roi = float(next(it))
            elif a == "--max-patch":
                max_patch = float(next(it))
            elif a == "--ratio":
                ratio = float(next(it))
            else:
                names.append(a)
    except (StopIteration, ValueError):
        names = []
    if len(names) not in (3, 4) or names[0] not in SHAPES or (ratio is not None and len(names) != 4):
        print(__doc__)
        return 2
    shape = names[0]
    patch, ref = Run(names[1]), Run(names[2])
    mine = norms(shape, patch, ref, patch, cell_ratio(patch, ref))
    show(shape, f"{patch.name} vs {ref.name}", mine)
    bad = []
    roi = mine["roi"]
    if roi is None:
        bad.append("the patch has no cell in the ROI")
    elif rms_roi is not None and not roi[0] <= rms_roi:
        bad.append(f"ROI: velocity RMS {roi[0]:.4e} of U_ref is above {rms_roi:g}")
    if max_patch is not None and not mine["patch"][1] <= max_patch:
        bad.append(f"patch: largest velocity difference {mine['patch'][1]:.4e} of U_ref is above {max_patch:g}")
    if len(names) == 4:
        parent = Run(names[3])
        theirs = norms(shape, parent, ref, patch, cell_ratio(parent, ref))
        show(shape, f"{parent.name} vs {ref.name}", theirs)
        if ratio is not None and not (roi is not None and theirs["roi"] is not None
                                      and roi[0] <= ratio * theirs["roi"][0]):
            bad.append(f"ROI: the patch's velocity RMS is above {ratio:g} x its parent's")
    for b in bad:
        print("nest_shapes: FAIL: " + b)
    if not bad:
        print("nest_shapes: PASS")
    return 1 if bad else 0


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
        mine, theirs = norms(shape, p, r, p, cell_ratio(p, r)), norms(shape, q, r, p, cell_ratio(q, r))
        show(shape, f"{patch} vs {ref}", mine)
        show(shape, f"{parent} vs {ref}", theirs)
        if parent != f"{shape}_par":   # a chain: the whole-domain parent too, for the record
            w = Run(f"{shape}_par")
            show(shape, f"{shape}_par vs {ref}", norms(shape, w, r, p, cell_ratio(w, r)))
        rows.append((shape, case, patch, p.nx * p.ny * p.nz, its, its_ref, mine, theirs))
        return mine["roi"][0], theirs["roi"][0]

    for shape in shapes:
        s = SHAPES[shape]
        levels = s.get("levels", 1)
        if levels not in (1, 2):
            fail(f"{shape}: levels = {levels}; the study runs one or two nests")
        if s.get("direct") and s["direct"] != 2 ** levels:
            fail(f"{shape}: direct = {s['direct']}; it must be 2**levels = {2 ** levels}, the reference's refinement")
        if levels == 2 and not s["mid"] >= max(s["ds"]) + 2.0 * s["pdx"]:
            fail(f"{shape}: mid = {s['mid']:g} m; it must be at least the largest D plus two parent cells")
        fdx = s["pdx"] / 2 ** levels
        parent = f"{shape}_par" if levels == 1 else f"{shape}_mid"
        try:
            it_ref = run_domain(shape, f"{shape}_ref", fdx, deck(shape, "domain"))
            it_par = run_domain(shape, f"{shape}_par", s["pdx"], deck(shape, "domain"))
        except SystemExit as err:
            print(err)
            judge("T4", False, f"{shape}: a whole-domain run did not reach the steady verdict")
            continue
        print(f"nest_shapes: {shape}: whole domain {fdx:g} m {it_ref} iterations, {s['pdx']:g} m {it_par}")
        if levels == 2:
            try:
                it_mid = run_nest(shape, parent, f"{shape}_par", s["mid"], deck(shape, "nest"))
            except SystemExit as err:
                print(err)
                judge("T4", False, f"{parent} did not reach the steady verdict")
                continue
            print(f"nest_shapes: {shape}: middle patch {0.5 * s['pdx']:g} m, D = {s['mid']:g} m, {it_mid} iterations")
            row(shape, f"{shape}, middle patch, D = {s['mid']:g} m, from {shape}_par", parent, f"{shape}_ref",
                f"{shape}_par", it_mid, it_ref)
        hs = []
        for d in s["ds"]:
            name = f"{shape}_d{d:g}"
            try:
                its = run_nest(shape, name, parent, d, deck(shape, "nest"))
            except SystemExit as err:
                print(err)
                judge("T4", False, f"{name} did not reach the steady verdict")
                continue
            case = f"{shape}, D = {d:g} m" + ("" if levels == 1 else f", from {parent}")
            h, ph = row(shape, case, name, f"{shape}_ref", parent, its, it_ref)
            hs.append((d, h, ph))
        if s.get("direct"):   # the same patches straight from the whole-domain parent (no target)
            for d in s["ds"]:
                name = f"{shape}_r{s['direct']}_d{d:g}"
                try:
                    its = run_nest(shape, name, f"{shape}_par", d, deck(shape, "nest"), ratio=s["direct"])
                except SystemExit as err:
                    print(err)
                    judge("T4", False, f"{name} did not reach the steady verdict")
                    continue
                row(shape, f"{shape}, D = {d:g} m, ratio {s['direct']} from {shape}_par", name, f"{shape}_ref",
                    f"{shape}_par", its, it_ref)
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
                    its = run_nest(shape, name, parent, d, cmd, [mgs])
                except SystemExit as err:
                    print(err)
                    judge("T3", False, f"{name} did not reach the steady verdict")
                    continue
                h, _ = row(shape, f"{shape}, D = {d:g} m, four boxes, {np} ranks", name, f"{shape}_ref",
                           parent, its, it_ref)
                vals.append(h)
            if len(vals) == 3:
                # The same four boxes on every rank count: the fields differ by round-off only.
                gaps = [norms(shape, Run(f"{shape}_np{np}"), Run(f"{shape}_np1"), Run(f"{shape}_np1"), 1)["patch"][1]
                        for np in (2, 4)]
                judge("T3", max(gaps) <= RANK_TOL,
                      f"{shape}: largest velocity difference from 1 rank {gaps[0]:.3e} (2 ranks), {gaps[1]:.3e} "
                      f"(4 ranks) of U_ref, limit {RANK_TOL:g}; ROI RMS {', '.join(f'{100 * v:.4f}' for v in vals)} %")
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
    elif args[0] == "nest" and len(args) in (5, 7) and args[1] in SHAPES:
        ratio = 2
        if len(args) == 7:
            if args[5] != "--nest-ratio" or not args[6].isdigit():
                print(__doc__)
                return 2
            ratio = int(args[6])
            if ratio not in (2, 4, 8, 16):
                fail(f"--nest-ratio {ratio}: 2, 4, 8 or 16")
        run_nest(args[1], args[2], args[3], half_widths(args[4]), cmd, ratio=ratio)
    elif args[0] == "study":
        shapes, ranks = [k for k in SHAPES if not SHAPES[k].get("manual")], []
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
