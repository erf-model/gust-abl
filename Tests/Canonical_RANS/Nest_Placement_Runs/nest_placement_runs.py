#!/usr/bin/env python3
"""The step rule and oblique wind of tools/nest_placement.py against joint runs.

    nest_placement_runs.py [--only WIND ...] -- COMMAND [ARGS ...]

COMMAND runs GUST on a deck (an argument ending in inputs_outer); ARGS go to every run after the deck. In the
working directory the script

1. writes the terrain file escarpment.xyz: z = (H/2)(1 + tanh(s/L)), H = 60 m, L = 60 m, s the distance from the
   line through the middle of the domain at 35 degrees to the x axis, on points 12.5 m apart;
2. runs tools/nest_placement.py on the two nests (NESTS) for the three winds (WINDS): xcut, whose x faces cut the
   step, and ycut, whose y faces cut it. The contour lines meet an x face at 35 degrees to its normal and a y face
   at 55 degrees, so the tool flags xcut for the wind from 270, ycut for the wind from 180, and both for the oblique
   wind from 225 (every face is then inflow or outflow);
3. for each wind: the uniform 25 m reference (ref_<wind>/) and, for each nest, the joint run
   (nest.feedback = defect) in <nest>_<wind>/ (the first plotfiles are the one-way pass, the last the joint run)
   and again with the conservative faces (nest.feedback_flux = nest) in <nest>_<wind>_flux/;
4. prints, for each nest and wind, the tool's verdict and the run's numbers against the reference: RMS of u, v, w
   and k over the nest (the reference at the nest's cell heights) and over the outer cells outside the footprint
   (the reference's 2 x 2 means), one way, joint and joint with the conservative faces, and the cells of the lowest 200 m whose horizontal speed
   differs from the reference's by more than 0.25, 0.5 and 1 m/s (tools/speed_difference.py). The table is also
   written to nest_placement_runs.csv.

--only runs the winds named (270, 180, 225). A run folder whose log already has the steady verdict is reused, not
run again (delete it to rerun). No verdict: the numbers are those of README.md. The exit code is 0
when every run exits 0 and reaches the steady verdict.
"""
import csv
import glob
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..", "..")
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(ROOT, "extern_erf"))
sys.path.insert(0, os.path.join(ROOT, "tools"))
import erf_plotfile  # noqa: E402
import speed_difference as sd  # noqa: E402

DOMAIN = (0.0, 0.0, 3000.0, 2400.0)
DX = 50.0
N_OUT, N_REF, NZ = (60, 48), (120, 96), 16
H, LW, BETA = 60.0, 60.0, 35.0           # escarpment height and width (m), contour angle to the x axis (degrees)
CENTRE = (1500.0, 1200.0)
NESTS = {"xcut": ((1250.0, 800.0), (1750.0, 1600.0)),    # the step crosses the xlo and xhi faces
         "ycut": ((750.0, 1000.0), (2250.0, 1400.0))}    # the step crosses the ylo and yhi faces
# Wind from (degrees from north, as nest_placement.py reads it): GUST's bc.inflow_direction (towards, counter-
# clockwise from x) and the side faces.
WINDS = {270: (0.0, {"xlo": "inflow", "xhi": "outflow", "ylo": "slip_wall", "yhi": "slip_wall"}),
         180: (90.0, {"xlo": "slip_wall", "xhi": "slip_wall", "ylo": "inflow", "yhi": "outflow"}),
         225: (45.0, {"xlo": "inflow", "xhi": "outflow", "ylo": "inflow", "yhi": "outflow"})}
FIELDS = ["u", "v", "w", "tke"]


def fail(msg):
    raise SystemExit("nest_placement_runs: " + msg)


def write_terrain(path):
    step = 12.5
    xs = np.arange(DOMAIN[0], DOMAIN[2] + 0.5 * step, step)
    ys = np.arange(DOMAIN[1], DOMAIN[3] + 0.5 * step, step)
    X, Y = np.meshgrid(xs, ys, indexing="ij")
    b = math.radians(BETA)
    s = (X - CENTRE[0]) * math.sin(b) - (Y - CENTRE[1]) * math.cos(b)   # > 0 south-east of the line
    Z = 0.5 * H * (1.0 + np.tanh(s / LW))
    with open(path, "w") as fh:
        fh.write(f"# escarpment H = {H:g} m, L = {LW:g} m, contours at {BETA:g} degrees to x through "
                 f"({CENTRE[0]:g}, {CENTRE[1]:g}) (nest_placement_runs.py)\n")
        np.savetxt(fh, np.column_stack([X.ravel(), Y.ravel(), Z.ravel()]), fmt="%.3f")


def placement(name, wind):
    """The tool's verdict for one nest and wind: (exit status, flagged faces)."""
    (x0, y0), (x1, y1) = NESTS[name]
    cmd = [sys.executable, os.path.join(ROOT, "tools", "nest_placement.py"), "--terrain", "escarpment.xyz",
           "--domain", *[f"{v:g}" for v in DOMAIN], "--dx", f"{DX:g}", "--nest", f"{x0:g}", f"{y0:g}", f"{x1:g}",
           f"{y1:g}", "--direction", f"{wind:g}"]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    with open(f"placement_{name}_{wind}.txt", "w") as fh:
        fh.write(proc.stdout)
    verdict = next((ln for ln in proc.stdout.splitlines() if ln.strip().startswith("verdict:")), "")
    flagged = sorted(set(re.findall(r"(\w+) \((?:inflow|outflow)\) cuts", verdict)))
    if proc.returncode not in (0, 1):
        fail(f"nest_placement.py exited {proc.returncode} for {name} at {wind}:\n{proc.stdout}")
    return proc.returncode, flagged


def run(cmd, deck_at, args, where):
    logfile = os.path.join(where, "log")
    if os.path.exists(logfile):
        log = open(logfile).read()
        if re.search(r"^gust: (joint run: )?steady after \d+", log, re.M):
            return log
    os.makedirs(where, exist_ok=True)
    full = cmd[:deck_at + 1] + args + cmd[deck_at + 1:]
    with open(os.path.join(where, "log"), "w") as fh:
        rc = subprocess.run(full, cwd=where, stdout=fh, stderr=subprocess.STDOUT).returncode
    log = open(os.path.join(where, "log")).read()
    if rc != 0:
        fail(f"{where}: the run exited {rc}; see {where}/log")
    if not re.search(r"^gust: (joint run: )?steady after \d+", log, re.M):
        fail(f"{where}: no steady verdict; see {where}/log")
    return log


def plotfiles(where, prefix):
    found = [(int(m.group(1)), p) for p in glob.glob(os.path.join(where, prefix + "[0-9]*"))
             for m in [re.fullmatch(re.escape(prefix) + r"(\d+)", os.path.basename(p))] if m and os.path.isdir(p)]
    return [p for _, p in sorted(found)]


def read(p):
    return {k: np.array(v, dtype=float) for k, v in erf_plotfile.read_fields(p, FIELDS + ["z"])[1].items()}


def rms(a):
    return float(np.sqrt(np.mean(a * a)))


def compare(ref, where, name):
    """One way and joint against the reference: {(region, field, pass): RMS} and the speed-difference counts."""
    (x0, y0), (x1, y1) = NESTS[name]
    i0, j0 = round((x0 - DOMAIN[0]) / DX), round((y0 - DOMAIN[1]) / DX)
    ni, nj = round((x1 - x0) / DX), round((y1 - y0) / DX)
    sl = (slice(2 * i0, 2 * (i0 + ni)), slice(2 * j0, 2 * (j0 + nj)))
    rc = {f: sd.coarsen(ref[f]) for f in FIELDS + ["z"]}
    outside = np.ones(rc["u"].shape[:2], bool)
    outside[i0:i0 + ni, j0:j0 + nj] = False
    m3 = np.broadcast_to(outside[:, :, None], rc["u"].shape)
    outer, nest = plotfiles(where, "outer_plt"), plotfiles(where, "outer_plt_" + name)
    if len(outer) != 2 or len(nest) != 2:
        fail(f"{where}: expected two plotfiles each (one way, joint): {outer}, {nest}")
    out = {}
    rs, rcs = sd.speed(ref), sd.speed(rc)
    for k, (po, pn) in zip(("ow", "jt"), zip(outer, nest)):
        O, N = read(po), read(pn)
        for f in FIELDS:
            out[("nest", f, k)] = rms(N[f] - sd.at_height(N["z"], ref["z"][sl], ref[f][sl]))
            out[("outer", f, k)] = rms((O[f] - rc[f])[m3])
        dn = sd.speed(N) - sd.at_height(N["z"], ref["z"][sl], rs[sl])
        do = sd.speed(O) - rcs
        out[("nest", "cells", k)] = sd.counts(dn[sd.near_ground(N["z"])])
        out[("outer", "cells", k)] = sd.counts(do[sd.near_ground(rc["z"]) & m3])
    return out


def main():
    if "--" not in sys.argv:
        print(__doc__)
        return 2
    opts = sys.argv[1:sys.argv.index("--")]
    cmd = sys.argv[sys.argv.index("--") + 1:]
    winds = [int(w) for w in opts[opts.index("--only") + 1:]] if "--only" in opts else list(WINDS)
    if any(w not in WINDS for w in winds):
        fail(f"--only takes winds from {sorted(WINDS)}")
    deck_at = next((i for i, a in enumerate(cmd) if a.endswith("inputs_outer")), None)
    if deck_at is None:
        fail("COMMAND has no inputs_outer")
    cmd[deck_at] = os.path.abspath(cmd[deck_at])
    write_terrain("escarpment.xyz")
    terrain = [f"terrain.file={os.path.abspath('escarpment.xyz')}"]
    rows = []
    for wind in winds:
        direction, faces = WINDS[wind]
        base = terrain + [f"bc.inflow_direction={direction:g}"] + [f"bc.{k}={v}" for k, v in faces.items()]
        where = f"ref_{wind}"
        run(cmd, deck_at, base + [f"mesh.n_cell={N_REF[0]} {N_REF[1]} {NZ}", "io.plot_file=ref_plt"], where)
        ref = read(plotfiles(where, "ref_plt")[-1])
        for name, ((x0, y0), (x1, y1)) in NESTS.items():
            status, flagged = placement(name, wind)
            where = f"{name}_{wind}"
            nest = base + [f"mesh.n_cell={N_OUT[0]} {N_OUT[1]} {NZ}", f"nest.names={name}",
                           f"nest.{name}.lo={x0:g} {y0:g}", f"nest.{name}.hi={x1:g} {y1:g}",
                           "nest.feedback=defect", "io.plot_file=outer_plt"]
            n = []
            for suffix, extra in (("", []), ("_flux", ["nest.feedback_flux=nest"])):
                log = run(cmd, deck_at, nest + extra, where + suffix)
                m = re.search(r"^gust: joint run: steady after (\d+) exchanges", log, re.M)
                n.append(int(m.group(1)) if m else -1)
            res = compare(ref, where, name)
            for key, val in compare(ref, where + "_flux", name).items():
                if key[2] == "jt":
                    res[(key[0], key[1], "flux")] = val
            verdict = "flagged " + "+".join(flagged) if status == 1 else "clear"
            rows.append((wind, name, verdict, n, res))
            line = f"{wind:3d} {name:4s} {verdict:22s} {n[0]:4d}/{n[1]:4d} exch |"
            for region in ("nest", "outer"):
                ow, jt, fl = (res[(region, "u", k)] for k in ("ow", "jt", "flux"))
                line += f" {region} u {ow:.3e} -> {jt:.3e} ({jt / ow:.2f}), flux {fl:.3e} ({fl / ow:.2f})"
                line += " cells>0.25 " + " -> ".join(str(res[(region, "cells", k)][0]) for k in ("ow", "jt", "flux"))
                line += " |"
            print("nest_placement_runs: " + line, flush=True)
    with open("nest_placement_runs.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["wind", "nest", "tool", "exchanges", "exchanges_flux", "region", "field", "one_way", "joint",
                    "joint_flux"])
        for wind, name, verdict, n, res in rows:
            for region in ("nest", "outer"):
                for f in FIELDS:
                    w.writerow([wind, name, verdict, *n, region, f,
                                *(f"{res[(region, f, k)]:.4e}" for k in ("ow", "jt", "flux"))])
                for t, a, b, c in zip(sd.THRESHOLDS, *(res[(region, "cells", k)] for k in ("ow", "jt", "flux"))):
                    w.writerow([wind, name, verdict, *n, region, f"cells>{t:g}", a, b, c])
    print("nest_placement_runs: wrote nest_placement_runs.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
