#!/usr/bin/env python3
"""A step of the roughness length under a flat ground: immersed and hybrid against body-fitted (README.md).

    check_roughness_step.py GUST OUTDIR [--limits DU DK DKW] [--mpi "EXE FLAG"] [--fresh] [key=value ...]
                            [--immersed F ...] [--hybrid H ...]

Runs, each with bc.surface_z0_file = z0_step.xyz (z0 = 0.03 m for x < 1120 m, 0.3 m beyond):
  * inputs_body_fitted_step, and once more with the uniform bc.surface_z0 = 0.03 (no map: the control);
  * inputs_immersed_step with terrain.flat_height = (2 + F) dz for every F after --immersed (default 0 0.3 0.7;
    F is the cut cell's solid fraction);
  * inputs_hybrid_step with terrain.flat_height = H dz for every H after --hybrid (default 0 0.5; H = 0: the
    lower face is the ground and the surface layer acts on it; else the ground cuts the first layer).
The runs take 2, 4, 1, 2, ... ranks in this order when --mpi gives the launcher (its executable and the flag
of the rank count, "mpiexec -n"); one rank each without it.

Printed for each run: the iterations and the largest differences from the body-fitted run, in percent:
  * u over every column and every fluid cell above the wall cell, at equal height above the ground (the
    body-fitted column interpolated linearly in ln(z + z0), z0 the column's);
  * tke the same way, over the columns away from the step ("dk far") and over the NEAR columns downstream of it
    ("dk near": the first 400 m, where the new internal layer grows from less than one layer to a few and tke
    changes across it within one or two cells, so that a value between two cell centres is not resolved:
    printed, no limit);
  * tke in the wall cell of every column (the cut cell, or the first cell of a ground on a face) against the
    body-fitted first cell ("dk wall"): the wall law's u*^2 with the column's z0, the step included.
The control gives the size of the step's effect: the largest difference of u between the body-fitted runs with
the map and with the uniform z0.

key=value arguments replace or add keys in every run's deck (key= with no value removes the key), the
terrain.* ones in the immersed and hybrid runs only. A run folder with a `done` file is read, not run again,
unless --fresh is given. With --limits the script is a test: it exits with 1 when a run has no steady verdict,
when its columns differ along y by more than 1e-8 of the field's largest value, when u, tke away from the step or the wall cell's tke differs from the
body-fitted run by more than DU, DK or DKW percent, or when the step's effect on u is below 5 DU percent (the
map would not be seen by the limits).
"""
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

DZ = 15.625
DX, NX = 40.0, 64
X_STEP = 1120.0
Z0_UP, Z0_DOWN = 0.03, 0.3
NEAR = 10   # columns downstream of the step with no limit on tke above the wall cell
MAP = os.path.join(HERE, "z0_step.xyz")

FRESH = False
MPI = None
RANKS = (2, 4, 1)
COUNT = [0]


def run(gust, deck, cwd, extra, prefix):
    """Run deck with the keys of extra in cwd; (iterations or None, plotfile or None, ranks, log)."""
    os.makedirs(cwd, exist_ok=True)
    log = os.path.join(cwd, "log")
    ranks = RANKS[COUNT[0] % len(RANKS)] if MPI else 1
    COUNT[0] += 1
    if FRESH or not os.path.exists(os.path.join(cwd, "done")):
        # The run's deck: a copy with each key=value replacing the deck's line of that key, or added (the
        # reader refuses a key given twice); key= with no value removes the key.
        lines = open(os.path.join(HERE, deck)).read().splitlines()
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
    text = open(log).read()
    m = re.search(r"gust: steady after (\d+)", text)
    its = int(m.group(1)) if m else None
    # The plotfile this run wrote: an earlier run's may lie in the folder.
    wrote = [ln.split()[-1] for ln in text.splitlines() if ln.startswith("gust: wrote " + prefix + "_plt")]
    return its, (os.path.join(cwd, wrote[-1]) if wrote else None), ranks, text


def fields(plt):
    _, d = erf_plotfile.read_fields(plt, ["u", "tke"])
    return {k: np.asarray(v) for k, v in d.items()}


def column_z0(i):
    return Z0_UP if (i + 0.5) * DX < X_STEP else Z0_DOWN


def heights(n, h):
    """The fluid centroid of each of n layers above a ground at h (-1 in a solid layer); the wall cell's layer."""
    lo, hi = np.arange(n) * DZ, (np.arange(n) + 1) * DZ
    open_h = np.clip(hi - np.maximum(lo, h), 0, DZ)
    zc = np.where(open_h > 1e-9, hi - 0.5 * open_h - h, -1.0)
    return zc, int(np.argmax(open_h > 1e-9))


def compare(d, h, b, zb):
    """Largest differences from the body-fitted run in percent: u, tke far from the step, tke near it, wall-cell tke."""
    zc, kc = heights(d["u"].shape[2], h)
    # Above the wall cell (the cut cell, or the first cell of a ground on a face), inside the reference's heights.
    sel = np.arange(kc + 1, d["u"].shape[2])
    sel = sel[(zc[sel] > zb[0]) & (zc[sel] < zb[-1])]
    out = np.zeros(4)
    first = int(round(X_STEP / DX))
    for i in range(NX):
        z0 = column_z0(i)
        xs, xb = np.log(zc[sel] + z0), np.log(zb + z0)
        du = np.abs(d["u"][i, 0, sel] / np.interp(xs, xb, b["u"][i, 0, :]) - 1).max() * 100
        dk = np.abs(d["tke"][i, 0, sel] / np.interp(xs, xb, b["tke"][i, 0, :]) - 1).max() * 100
        near = first <= i < first + NEAR
        out[0] = max(out[0], du)
        out[2 if near else 1] = max(out[2 if near else 1], dk)
        out[3] = max(out[3], abs(d["tke"][i, 0, kc] / b["tke"][i, 0, 0] - 1) * 100)
    return out


def main(argv):
    global FRESH, MPI
    gust, out = os.path.abspath(argv[1]), os.path.abspath(argv[2])
    rest, limits, keys, bad = argv[3:], None, [], []
    lists = {"--immersed": None, "--hybrid": None}
    into = None
    while rest:
        a = rest.pop(0)
        if a == "--limits":
            limits, into = (float(rest.pop(0)), float(rest.pop(0)), float(rest.pop(0))), None
        elif a == "--mpi":
            MPI, into = shlex.split(rest.pop(0)), None
        elif a == "--fresh":
            FRESH, into = True, None
        elif a in lists:
            lists[a] = []
            into = lists[a]
        elif "=" in a:
            keys.append(a)
        elif into is not None:
            into.append(float(a))
        else:
            raise SystemExit("check_roughness_step: unknown argument " + a)
    immersed = [0, 0.3, 0.7] if lists["--immersed"] is None else lists["--immersed"]
    hybrid = [0, 0.5] if lists["--hybrid"] is None else lists["--hybrid"]
    use_map = ["bc.surface_z0_file=" + MAP]
    plain = [k for k in keys if not k.startswith("terrain.") or k.startswith("terrain.idw")]

    its, plt, ranks, _ = run(gust, "inputs_body_fitted_step", os.path.join(out, "body_fitted"), use_map + plain, "body_fitted_step")
    if its is None or plt is None:
        print("FAIL body-fitted: no steady verdict")
        return 1
    b = fields(plt)
    zb = (np.arange(b["u"].shape[2]) + 0.5) * DZ
    print("body-fitted with the map: %d iterations on %d rank(s)" % (its, ranks))
    its, plt, ranks, _ = run(gust, "inputs_body_fitted_step", os.path.join(out, "uniform"), ["bc.surface_z0=%g" % Z0_UP] + plain + ["terrain.idw_neighbors=", "terrain.idw_power="], "body_fitted_step")
    if its is None or plt is None:
        print("FAIL body-fitted with the uniform z0: no steady verdict")
        return 1
    c = fields(plt)
    effect = np.abs(c["u"][:, 0, :] / b["u"][:, 0, :] - 1).max() * 100
    print("body-fitted with the uniform z0 = %g m: %d iterations on %d rank(s); the step's effect on u: %.2f %%" % (Z0_UP, its, ranks, effect))
    if limits and effect < 5 * limits[0]:
        bad.append("the step's effect on u (%.2f %%) is below 5 times the limit on u" % effect)

    print("%-16s %5s %6s %9s %9s %9s %9s" % ("run", "ranks", "iter", "du", "dk far", "dk near", "dk wall"))
    cases = [("immersed F=%g" % f, "inputs_immersed_step", "immersed_step", "f%04.2f" % f, (2 + f) * DZ) for f in immersed]
    cases += [("hybrid H=%g" % h, "inputs_hybrid_step", "hybrid_step", "h%04.2f" % h, h * DZ) for h in hybrid]
    for name, deck, prefix, folder, h in cases:
        its, plt, ranks, text = run(gust, deck, os.path.join(out, folder), use_map + ["terrain.flat_height=%.10g" % h] + keys, prefix)
        if plt is None:
            print("%-16s  no plotfile: %s" % (name, text.strip().splitlines()[-1][:150] if text.strip() else "no output"))
            bad.append(name + ": no plotfile")
            continue
        d = fields(plt)
        spread = max(np.abs(d[k] - d[k][:, :1, :]).max() / np.abs(d[k]).max() for k in ("u", "tke"))
        e = compare(d, h, b, zb)
        print("%-16s %5d %6s %9.3f %9.3f %9.3f %9.3f" % (name, ranks, its, e[0], e[1], e[2], e[3]))
        if limits:
            if its is None:
                bad.append(name + ": no steady verdict")
            if spread > 1e-8:
                bad.append("%s: the columns differ along y by %.2e" % (name, spread))
            for n, what in ((0, "u"), (1, "tke away from the step"), (3, "the wall cell's tke")):
                if e[n] > limits[n if n < 3 else 2]:
                    bad.append("%s: %s differs from the body-fitted run by %.3f %% (limit %g)" % (name, what, e[n], limits[n if n < 3 else 2]))
    for m in bad:
        print("FAIL " + m)
    if limits and not bad:
        print("PASS")
    return 1 if (limits and bad) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
