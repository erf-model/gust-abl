#!/usr/bin/env python3
"""The canopy cases of the body-fitted mesh over an immersed ground and on a hybrid mesh (README.md).

    check_immersed_canopy.py GUST OUTDIR [--cases NAME ...] [--immersed F ...] [--hybrid H ...]
                             [--limits DU DK DKTOP] [--mpi "EXE FLAG"] [--fresh] [key=value ...]

Each case is a deck of another case folder, run as it is (body-fitted) and with the ground offset:
  homogeneous  ../Homogeneous_Canopy/inputs_homogeneous_canopy, LAI 5 (a uniform 20 m canopy, flat ground);
  edge         ../Forest_Edge/inputs_forest_edge, LAI 2 (forest, clearing, forest; flat ground);
  hill         ../Forested_Hill/inputs_deep_hill (a 20 m canopy, LAI 10, on cosine hills 4 m high).
The script writes the cases' tables (canopy height, LAI, terrain) as their own check scripts do, and runs:
  * the deck (body-fitted), with init.first_guess = profile;
  * terrain.representation = immersed_forcing for every F after --immersed (default 0 0.3 0.7): the deck's
    layers with two more layers of the first layer's thickness below them (mesh.zlevels_file) and
    terrain.flat_height = (2 + F) dz0, so that the ground cuts the fraction F off the third layer;
  * terrain.representation = hybrid for every H after --hybrid (default 0 0.5): the deck's mesh and
    terrain.flat_height = H dz0 (H = 0: the lower face is the ground where the mesh follows it).
The runs take 2, 4, 1, 2, ... ranks in this order when --mpi gives the launcher (its executable and the flag
of the rank count, "mpiexec -n"); one rank each without it.

Printed for each run: the iterations and the largest differences from the body-fitted run over every column
and every fluid cell above the wall cell, at equal height above the ground (the body-fitted column
interpolated linearly in height), in percent of the body-fitted run's largest value of the field:
  du        u at every height;
  dk in     tke below the canopy top (of the tallest canopy, also over the clearing);
  dk above  tke above the canopy top;
  dk top    tke within 1.5 layers of the canopy top (the body-fitted layer there), which the first two leave
            out: tke rises across one layer there, and a value between two layers' centres is not resolved.
Also printed (drag): the largest difference of a column's integral of the canopy's drag field (canopy_drag in
the plotfile) over the cells with a part above the ground from the body-fitted run's, in percent of C_d LAI. A canopy that
stood on the mesh's lower face would leave a part of its leaf area in the cells below the ground.

Every run has steady.sweep = zjacobi and solver.poisson_backend = gmres_fft, which the immersed ground needs.
key=value arguments replace or add keys in every run's deck (key= with no value removes the key); terrain.*
keys go to the immersed and hybrid runs only. A run
folder with a `done` file is read, not run again, unless --fresh is given. With --limits the script is a
test: it exits with 1 when a run has no steady verdict, when u or tke differs from the body-fitted run by
more than DU (u), DK (tke away from the canopy top) or DKTOP (tke at the top) percent, or when a column's drag integral departs from
the body-fitted run's by more than 1e-6 percent of C_d LAI.
"""
import math
import os
import re
import shlex
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402

CD = 0.2
FRESH = False
MPI = None
RANKS = (2, 4, 1)
COUNT = [0]


def hill_ground(x):
    return 0.5 * 4.0 * (1.0 + math.cos(math.pi * x / (2.0 * 40.0) + math.pi))


def flat_ground(x):
    return 0.0


# deck: the other case's deck; prefix: its io.plot_file; nz, dz0, ratio: its layers; lx, nx: its length and
# columns in x; hc: the tallest canopy; lai: the leaf area index; canopy(x): whether column x has canopy.
CASES = {
    "homogeneous": dict(deck="../Homogeneous_Canopy/inputs_homogeneous_canopy", prefix="homogeneous", nz=143, dz0=0.25,
                        ratio=1.02, lx=800.0, nx=8, hc=20.0, lai=5.0, ground=flat_ground, canopy=lambda x: True),
    "edge": dict(deck="../Forest_Edge/inputs_forest_edge", prefix="edge", nz=49, dz0=1.0, ratio=1.05, lx=690.0, nx=115,
                 hc=18.0, lai=2.0, ground=flat_ground, canopy=lambda x: (x % 690.0) <= 308.4),
    "hill": dict(deck="../Forested_Hill/inputs_deep_hill", prefix="hill", nz=39, dz0=1.0, ratio=1.07, lx=160.0, nx=40,
                 hc=20.0, lai=10.0, ground=hill_ground, canopy=lambda x: True),
}


def write_tables(name, c, cwd):
    """The case's tables, as its own check script writes them."""
    os.makedirs(cwd, exist_ok=True)
    if name == "homogeneous":
        pts = [(x, y) for y in [-200.0 + 200.0 * i for i in range(8)] for x in [-200.0 + 200.0 * i for i in range(8)]]
    elif name == "edge":
        pts = [(3.0 * i, y) for y in (-9.0, -3.0, 3.0, 9.0, 15.0, 21.0, 27.0, 33.0) for i in range(-8, 240)]
    else:
        pts = [(i * 1.0, j * 4.0) for j in range(-2, 7) for i in range(-8, 169)]
    with open(os.path.join(cwd, "h.xyz"), "w") as h, open(os.path.join(cwd, "lai.xyz"), "w") as a:
        for x, y in pts:
            on = c["canopy"](x)
            h.write("%g %g %g\n" % (x, y, c["hc"] if on else 0.0))
            a.write("%g %g %g\n" % (x, y, c["lai"] if on else 0.0))
    if name == "hill":
        with open(os.path.join(cwd, "t.xyz"), "w") as t:
            for x, y in pts:
                t.write("%g %g %.6f\n" % (x, y, hill_ground(x)))


def levels(c, below=0):
    """The deck's levels above `below` more layers of the first layer's thickness."""
    z = [k * c["dz0"] for k in range(below + 1)]
    for k in range(c["nz"]):
        z.append(z[-1] + c["dz0"] * c["ratio"] ** k)
    return z


def run(gust, name, c, cwd, extra):
    """Run the case's deck with the keys of extra in cwd; (iterations or None, plotfile or None, ranks, log)."""
    write_tables(name, c, cwd)
    log = os.path.join(cwd, "log")
    ranks = RANKS[COUNT[0] % len(RANKS)] if MPI else 1
    COUNT[0] += 1
    if FRESH or not os.path.exists(os.path.join(cwd, "done")):
        # The run's deck: a copy with each key=value replacing the deck's line of that key, or added (the
        # reader refuses a key given twice); key= with no value removes the key.
        lines = open(os.path.join(HERE, c["deck"])).read().splitlines()
        for kv in extra:
            key, value = (x.strip() for x in kv.split("=", 1))
            at = [n for n, ln in enumerate(lines) if ln.split("=", 1)[0].strip() == key and not ln.lstrip().startswith("#")]
            if not value:
                lines = [ln for n, ln in enumerate(lines) if n not in at]
            elif at:
                lines[at[0]] = kv
            else:
                lines.append(kv)
        open(os.path.join(cwd, "deck"), "w").write("\n".join(lines) + "\n")
        cmd = (MPI + [str(ranks)] if MPI else []) + [gust, "deck"]
        with open(log, "w") as f:
            code = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT).returncode
        open(os.path.join(cwd, "done"), "w").write("%d %d\n" % (code, ranks))
    else:
        ranks = int(open(os.path.join(cwd, "done")).read().split()[1])   # those the folder's run took
    code = int(open(os.path.join(cwd, "done")).read().split()[0])
    text = open(log).read()
    m = re.search(r"gust: steady after (\d+)", text)
    its = int(m.group(1)) if (m and code == 0) else None   # a run that ends with an error has no verdict
    # The plotfile this run wrote: an earlier run's may lie in the folder.
    wrote = [ln.split()[-1] for ln in text.splitlines() if ln.startswith("gust: wrote " + c["prefix"] + "_plt")]
    return its, (os.path.join(cwd, wrote[-1]) if wrote else None), ranks, text


def columns(plt, c, ztop, offset):
    """Per column i (the fields are uniform in y): the fluid centroids' heights above the ground (-1 in a
    solid cell), u, tke, the wall cell's layer and the integral of canopy_drag over the column's cells that hold fluid."""
    _, d = erf_plotfile.read_fields(plt, ["u", "tke", "z", "detJ", "canopy_drag"])
    a = {k: np.asarray(v)[:, 0, :] for k, v in d.items()}
    nx, nz = a["u"].shape
    dz = a["detJ"] * ztop / nz           # the cells' heights
    lo, hi = a["z"] - 0.5 * dz, a["z"] + 0.5 * dz
    dx = c["lx"] / nx
    out = []
    for i in range(nx):
        g = offset + 0.5 * (c["ground"](i * dx) + c["ground"]((i + 1) * dx))   # the mean of the column's nodes
        open_h = np.clip(hi[i] - np.maximum(lo[i], g), 0, None)
        fluid = open_h > 1e-6 * dz[i]
        zc = np.where(fluid, hi[i] - 0.5 * open_h - g, -1.0)
        out.append(dict(z=zc, u=a["u"][i], tke=a["tke"][i], wall=int(np.argmax(fluid)),
                        drag=float((a["canopy_drag"][i] * dz[i])[fluid].sum())))
    return out


def compare(run_cols, ref_cols, c):
    """The largest differences from the body-fitted run in percent of its largest u and tke: (du, dk in,
    dk above, dk top), and the largest difference of the drag integrals in percent of C_d LAI."""
    umax = max(np.abs(r["u"]).max() for r in ref_cols)
    kmax = max(np.abs(r["tke"]).max() for r in ref_cols)
    out = np.zeros(4)
    drag = 0.0
    for d, b in zip(run_cols, ref_cols):
        sel = np.arange(d["wall"] + 1, len(d["z"]))
        sel = sel[(d["z"][sel] > b["z"][0]) & (d["z"][sel] < b["z"][-1])]
        z = d["z"][sel]
        if len(z) < len(b["z"]) // 2:
            raise SystemExit("check_immersed_canopy: a column has %d cells to compare (of %d)" % (len(z), len(b["z"])))
        du = np.abs(d["u"][sel] - np.interp(z, b["z"], b["u"])) / umax * 100
        dk = np.abs(d["tke"][sel] - np.interp(z, b["z"], b["tke"])) / kmax * 100
        # The band at the canopy top: 1.5 of the body-fitted layer there (the spacing of its centres).
        band = 1.5 * np.interp(c["hc"], 0.5 * (b["z"][1:] + b["z"][:-1]), np.diff(b["z"]))
        for n, (v, m) in enumerate(((du, z > -1), (dk, z < c["hc"] - band), (dk, z > c["hc"] + band),
                                    (dk, np.abs(z - c["hc"]) <= band))):
            if m.any():
                out[n] = max(out[n], v[m].max())
        drag = max(drag, abs(d["drag"] - b["drag"]) / (CD * c["lai"]) * 100)
    return out, drag


def main(argv):
    global FRESH, MPI
    gust, out = os.path.abspath(argv[1]), os.path.abspath(argv[2])
    rest, limits, keys, bad = argv[3:], None, [], []
    lists = {"--immersed": None, "--hybrid": None, "--cases": None}
    into = None
    while rest:
        a = rest.pop(0)
        if a == "--limits":
            limits, into = [float(rest.pop(0)) for _ in range(3)], None
        elif a == "--mpi":
            MPI, into = shlex.split(rest.pop(0)), None
        elif a == "--fresh":
            FRESH, into = True, None
        elif a in lists:
            lists[a] = []
            into = lists[a]
        elif "=" in a:
            keys.append(a)
            into = None
        elif into is not None:
            into.append(a)
        else:
            raise SystemExit("check_immersed_canopy: unknown argument " + a)
    immersed = [0.0, 0.3, 0.7] if lists["--immersed"] is None else [float(x) for x in lists["--immersed"]]
    hybrid = [0.0, 0.5] if lists["--hybrid"] is None else [float(x) for x in lists["--hybrid"]]
    names = list(CASES) if lists["--cases"] is None else lists["--cases"]
    for name in names:
        if name not in CASES:
            raise SystemExit("check_immersed_canopy: unknown case " + name)
        c = CASES[name]
        deck = dict((x.strip() for x in ln.split("#", 1)[0].split("=", 1)) for ln in open(os.path.join(HERE, c["deck"])).read().splitlines()
                    if "=" in ln.split("#", 1)[0])
        deck_n, hi = deck["mesh.n_cell"].split(), deck["mesh.prob_hi"].split()
        ztop = float(hi[2])               # detJ is a cell's height over prob_hi / n_cell
        where = os.path.join(out, name)
        common = ["init.first_guess=profile", "canopy.height_file=h.xyz", "canopy.lai_file=lai.xyz",
                  "steady.sweep=zjacobi", "steady.zline_sweeps=16", "steady.cfl_max=4", "solver.poisson_backend=gmres_fft"] + [k for k in keys if not k.startswith("terrain.")]
        ground_keys = [k for k in keys if k.startswith("terrain.")]
        its, plt, ranks, text = run(gust, name, c, os.path.join(where, "body_fitted"), common)
        if its is None or plt is None:
            print("FAIL %s body-fitted: no steady verdict" % name)
            bad.append(name + " body-fitted: no steady verdict")
            continue
        ref = columns(plt, c, ztop, 0.0)
        dx = c["lx"] / len(ref)
        whole = max(abs(r["drag"] / (CD * c["lai"]) - 1) * 100 for i, r in enumerate(ref)
                    if c["canopy"](i * dx) and c["canopy"]((i + 1) * dx))
        print("%s, body-fitted: %d iterations on %d rank(s); drag integral within %.1e %% of C_d LAI" % (name, its, ranks, whole))
        print("  %-16s %5s %6s %8s %8s %9s %9s %10s" % ("run", "ranks", "iter", "du", "dk in", "dk above", "dk top", "drag"))
        runs = []
        for f in immersed:
            z = levels(c, 2)
            folder = os.path.join(where, "f%04.2f" % f)
            os.makedirs(folder, exist_ok=True)
            open(os.path.join(folder, "zlevels.txt"), "w").write("".join("%.12g\n" % v for v in z))
            extra = ["terrain.representation=immersed_forcing", "terrain.flat_height=%.10g" % ((2 + f) * c["dz0"]),
                     "terrain.smoothing=", "mesh.dz0=", "mesh.stretch_ratio=", "mesh.zlevels_file=zlevels.txt"]
            runs.append(("immersed F=%g" % f, folder, extra, z[-1], (2 + f) * c["dz0"]))
        for h in hybrid:
            extra = ["terrain.representation=hybrid", "terrain.flat_height=%.10g" % (h * c["dz0"])]
            runs.append(("hybrid H=%g" % h, os.path.join(where, "h%04.2f" % h), extra, ztop, h * c["dz0"]))
        for label, folder, extra, top, offset in runs:
            # The immersed mesh has two more layers: its n_cell from the deck's.
            if label.startswith("immersed"):
                extra += ["mesh.n_cell=%s %s %d" % (deck_n[0], deck_n[1], int(deck_n[2]) + 2),
                          "mesh.prob_hi=%s %s %.12g" % (hi[0], hi[1], top)]
            its, plt, ranks, text = run(gust, name, c, folder, common + extra + ground_keys)
            if plt is None:
                print("  %-16s  no plotfile: %s" % (label, text.strip().splitlines()[-1][:150] if text.strip() else "no output"))
                bad.append("%s %s: no plotfile" % (name, label))
                continue
            e, drag = compare(columns(plt, c, top, offset), ref, c)
            print("  %-16s %5d %6s %8.3f %8.3f %9.3f %9.3f %10.1e" % (label, ranks, its, e[0], e[1], e[2], e[3], drag))
            if limits:
                if its is None:
                    bad.append("%s %s: no steady verdict" % (name, label))
                for n, what in enumerate(("u", "tke in the canopy", "tke above the canopy", "tke at the canopy top")):
                    lim = limits[(0, 1, 1, 2)[n]]
                    if e[n] > lim:
                        bad.append("%s %s: %s differs from the body-fitted run by %.3f %% (limit %g)" % (name, label, what, e[n], lim))
                if drag > 1e-6:
                    bad.append("%s %s: a column's drag integral departs from the body-fitted run's by %.2e %% of C_d LAI" % (name, label, drag))
    for m in bad:
        print("FAIL " + m)
    if limits and not bad:
        print("PASS")
    return 1 if (limits and bad) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
