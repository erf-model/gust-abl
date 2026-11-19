#!/usr/bin/env python3
"""Run GUST on the immersed ridge and on its body-fitted reference and compare them (rt_immersed_ridge_*).

    check_immersed_ridge.py --slope vertical|normal --body-fitted DECK [--wind across|along|hill]
                            [--crest PCT] [--lee PCT] [--rms PCT] [--max PCT] [--reuse] -- COMMAND [ARGS ...]

COMMAND runs the application on inputs_immersed_ridge (--wind across, the default),
inputs_immersed_along (--wind along) or inputs_immersed_hill (--wind hill, the round hill, read along
the row of columns just north of its centre); terrain.immersed_slope = SLOPE is added to it.
The same COMMAND with DECK in place of the immersed deck runs the body-fitted ridge in the sub-folder
body_fitted/ (--reuse: a body-fitted plotfile already there is read instead). Both must reach the
steady verdict. The wind (u across the ridge, v along it) is then read at fixed heights above the ground, d = 20, 40, 80 and 160 m,
in every column: in a column the fluid cells (building_beta = 0 in the immersed run) are taken at
their centres' heights above the column's mean ground and u is interpolated linearly in
ln((d + z0)/z0). At each d:
  * the crest speed (the largest u over the columns) of the two runs differs by at most --crest
    percent of the body-fitted value;
  * the lee minimum (the smallest u downstream of the crest) by at most --lee percent of the
    body-fitted far-field speed at that height;
  * the root mean square of the difference over the columns is at most --rms percent of that speed.
With --wind along there is no crest or lee: the largest difference over the columns is at most --max
percent of the body-fitted far-field speed, and the rms at most --rms percent.
Prints one row per check and the measured table; the exit code is the verdict.
"""
import glob
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402

Z0 = 0.1
HILL = 100.0
# Per pair of decks: the decks' tag, the wind's component, the columns along x (number, first node, width),
# terrain.hill_half_width, and for the round hill the row of columns read (the first north of the centre)
# with the nodes' y on its two sides.
CASES = {
    "across": {"tag": "ridge", "comp": "u", "nx": 128, "x0": -2560.0, "dx": 40.0, "half": 500.0, "row": 0, "y": None},
    "along": {"tag": "along", "comp": "v", "nx": 128, "x0": -2560.0, "dx": 40.0, "half": 200.0, "row": 0, "y": None},
    "hill": {"tag": "hill", "comp": "u", "nx": 32, "x0": -1280.0, "dx": 80.0, "half": 300.0, "row": 16, "y": (0.0, 80.0)},
}
BASE = 31.25          # terrain.flat_height of the immersed deck
HEIGHTS = (20.0, 40.0, 80.0, 160.0)


def hill(x, y, half):
    return HILL / (1 + (x * x + y * y) / half ** 2)


def ground(i, base, case):
    """The column's mean ground: the mean of its nodes' heights."""
    xs = (case["x0"] + i * case["dx"], case["x0"] + (i + 1) * case["dx"])
    ys = case["y"] or (0.0,)
    return base + sum(hill(x, y, case["half"]) for x in xs for y in ys) / (len(xs) * len(ys))


def run(cmd, cwd):
    os.makedirs(cwd, exist_ok=True)
    out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    sys.stdout.write(out.stdout[-1500:])
    sys.stderr.write(out.stderr[-1500:])
    return out.returncode, out.stdout + out.stderr


def written_plotfile(cwd, log, prefix):
    """The plotfile this run wrote ("gust: wrote <name>"): an earlier run's may lie in the folder."""
    wrote = [ln.split()[-1] for ln in log.splitlines() if ln.startswith("gust: wrote " + prefix)]
    if not wrote:
        raise SystemExit("check_immersed_ridge: the run wrote no " + prefix + "*")
    return os.path.join(cwd, wrote[-1])


def iterations(log):
    for ln in log.splitlines():
        if ln.startswith("gust: steady after"):
            return ln.split()[3]
    return "?"


def columns(plt, base, immersed, case):
    """Per column, the (height above the ground, wind) of its fluid cells, lowest first."""
    comp, NX, row = case["comp"], case["nx"], case["row"]
    names = [comp, "z"] + (["building_beta"] if immersed else [])
    _, f = erf_plotfile.read_fields(plt, names)
    if len(f[comp]) != NX:
        raise SystemExit("check_immersed_ridge: %s has %d columns, not %d" % (plt, len(f[comp]), NX))
    cols = []
    for i in range(NX):
        g = ground(i, base, case)
        pts = []
        for k in range(len(f[comp][i][row])):
            if immersed and f["building_beta"][i][row][k] != 0:
                continue
            pts.append((f["z"][i][row][k] - g, f[comp][i][row][k]))
        cols.append(pts)
    return cols


def at(pts, d):
    """u at the height d above the ground, linear in ln((d + z0)/z0); below the first fluid cell, the log law through it."""
    def lg(x):
        return math.log((x + Z0) / Z0)
    for a, b in zip(pts, pts[1:]):
        if a[0] <= d <= b[0]:
            return a[1] + (lg(d) - lg(a[0])) / (lg(b[0]) - lg(a[0])) * (b[1] - a[1])
    return pts[0][1] * lg(d) / lg(pts[0][0])


def main(argv):
    if "--" not in argv:
        raise SystemExit(__doc__)
    sep = argv.index("--")
    opts, cmd = argv[:sep], argv[sep + 1:]
    slope = body_fitted = None
    lim = {"--crest": 5.0, "--lee": 10.0, "--rms": 5.0, "--max": 5.0}
    wind = "across"
    reuse = False
    i = 0
    while i < len(opts):
        if opts[i] == "--reuse":
            reuse = True
            i += 1
            continue
        if opts[i] == "--slope":
            slope = opts[i + 1]
        elif opts[i] == "--wind":
            wind = opts[i + 1]
        elif opts[i] == "--body-fitted":
            body_fitted = opts[i + 1]
        elif opts[i] in lim:
            lim[opts[i]] = float(opts[i + 1])
        else:
            raise SystemExit("check_immersed_ridge: unknown option " + opts[i])
        i += 2
    if slope not in ("vertical", "normal") or body_fitted is None or wind not in CASES:
        raise SystemExit("check_immersed_ridge: --slope vertical|normal, --body-fitted DECK and --wind across|along|hill")
    case = CASES[wind]
    tag, NX = case["tag"], case["nx"]
    deck = [n for n, a in enumerate(cmd) if a.endswith("inputs_immersed_" + tag)]
    if len(deck) != 1:
        raise SystemExit("check_immersed_ridge: COMMAND names inputs_immersed_" + tag + " once")

    ok = True

    def check(name, good, detail):
        nonlocal ok
        ok = ok and good
        print(("PASS " if good else "FAIL ") + name + ": " + detail)

    code, log = run(cmd + ["terrain.immersed_slope=" + slope], os.getcwd())
    check("immersed steady verdict", code == 0 and "gust: steady after" in log, "exit %d, %s iterations" % (code, iterations(log)))
    # The run names the option it used: a key that were ignored would pass the comparisons below.
    check("terrain.immersed_slope in the log", "(terrain.immersed_slope = " + slope + ")" in log, slope)
    if not ok:
        return 1
    im = columns(written_plotfile(os.getcwd(), log, "immersed_" + tag + "_plt"), BASE, True, case)

    sub = os.path.join(os.getcwd(), "body_fitted")
    old = sorted(glob.glob(os.path.join(sub, "body_fitted_" + tag + "_plt*"))) if reuse else []
    if old:
        bf_plt = old[-1]
        print("reusing " + bf_plt)
    else:
        # No plotfile of an earlier body-fitted run may stay: --reuse of a later test reads the last by name.
        for stale in glob.glob(os.path.join(sub, "body_fitted_" + tag + "_*")):
            shutil.rmtree(stale)
        bf_cmd = list(cmd)
        bf_cmd[deck[0]] = os.path.abspath(body_fitted)
        code, log = run(bf_cmd, sub)
        check("body-fitted steady verdict", code == 0 and "gust: steady after" in log,
              "exit %d, %s iterations" % (code, iterations(log)))
        if not ok:
            return 1
        bf_plt = written_plotfile(sub, log, "body_fitted_" + tag + "_plt")
    bf = columns(bf_plt, 0.0, False, case)

    print("  d (m)   far field bf/im    crest bf/im (column)      lee minimum bf/im (column)    rms    max")
    for d in HEIGHTS:
        ub = [at(c, d) for c in bf]
        ui = [at(c, d) for c in im]
        nfar = NX // 16
        far = sum(ub[:nfar]) / nfar
        fari = sum(ui[:nfar]) / nfar
        cb = max(range(NX), key=lambda n: ub[n])
        ci = max(range(NX), key=lambda n: ui[n])
        lb = min(range(NX // 2, NX), key=lambda n: ub[n])
        li = min(range(NX // 2, NX), key=lambda n: ui[n])
        rms = math.sqrt(sum((a - b) ** 2 for a, b in zip(ub, ui)) / NX)
        mx = max(abs(a - b) for a, b in zip(ub, ui))
        print("  %5.0f   %6.3f / %6.3f    %6.3f / %6.3f (%d / %d)    %6.3f / %6.3f (%d / %d)    %.3f  %.3f"
              % (d, far, fari, ub[cb], ui[ci], cb, ci, ub[lb], ui[li], lb, li, rms, mx))
        if wind == "along":
            e = 100 * mx / far
            check("largest difference at %g m" % d, e <= lim["--max"], "%.2f %% of the far-field speed (max %.1f %%)" % (e, lim["--max"]))
            e = 100 * rms / far
            check("rms at %g m" % d, e <= lim["--rms"], "%.2f %% of the far-field speed (max %.1f %%)" % (e, lim["--rms"]))
            continue
        e = 100 * (ui[ci] / ub[cb] - 1)
        check("crest speed at %g m" % d, abs(e) <= lim["--crest"], "%+.2f %% (max %.1f %%)" % (e, lim["--crest"]))
        e = 100 * (ui[li] - ub[lb]) / far
        check("lee minimum at %g m" % d, abs(e) <= lim["--lee"], "%+.2f %% of the far-field speed (max %.1f %%)" % (e, lim["--lee"]))
        e = 100 * rms / far
        check("rms at %g m" % d, e <= lim["--rms"], "%.2f %% of the far-field speed (max %.1f %%)" % (e, lim["--rms"]))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
