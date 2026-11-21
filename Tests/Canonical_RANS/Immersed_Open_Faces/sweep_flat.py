#!/usr/bin/env python3
"""Flat ground with an inflow and an outflow face, immersed against body-fitted (README.md).

    sweep_flat.py GUST OUTDIR [--limits DU DK] [--fresh] [key=value ...] [F ...]

Runs inputs_body_fitted_flat once and inputs_immersed_flat with terrain.flat_height = (2 + F) dz for every
solid fraction F of the cut cell (default 0, 0.1, ..., 0.9, 0.99, 1), one rank each, and prints for each run:
iterations, the mass flux through the inflow and outflow faces (per column, the mean over y), max |w|, the pressure in the outflow
column, the speed of the three cells above the cut cell against the log law (and the cut cell's), and the
largest difference of u and tke from the body-fitted outflow column at equal height above the ground.

key=value arguments replace or add keys in every run's deck (turb.shear_stencil=log_law, ...; key= with no value
removes the key), the terrain.* ones (terrain.small_volfrac=0.1, ...) in the immersed runs only. A run folder
with a `done` file is read, not run again, unless --fresh is given. With --limits the script is a test: it exits
with 1 when a run has no steady verdict, when the outflow flux differs from the inflow flux by more than 1e-6 of
it, or when the largest difference of u (tke) from the body-fitted column exceeds DU (DK) percent.
"""
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402

DZ = 15.625
Z0, KAPPA, USTAR = 0.1, 0.41, 0.35


FRESH = False


def run(gust, deck, cwd, extra, prefix):
    os.makedirs(cwd, exist_ok=True)
    log = os.path.join(cwd, "log")
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
        with open(log, "w") as f:
            code = subprocess.run([gust, "deck"], cwd=cwd, stdout=f, stderr=subprocess.STDOUT).returncode
        open(os.path.join(cwd, "done"), "w").write(str(code))
    text = open(log).read()
    m = re.search(r"gust: steady after (\d+)", text)
    its = int(m.group(1)) if m else None
    wrote = [ln.split()[-1] for ln in text.splitlines() if ln.startswith("gust: wrote " + prefix + "_plt")]
    return its, (os.path.join(cwd, wrote[-1]) if wrote else None), text


def fields(plt):
    _, d = erf_plotfile.read_fields(plt, ["u", "w", "p", "tke"])
    return {k: np.asarray(v) for k, v in d.items()}


def main(argv):
    global FRESH
    gust, out = os.path.abspath(argv[1]), os.path.abspath(argv[2])
    rest, limits, keys, fs, bad = argv[3:], None, [], [], []
    while rest:
        a = rest.pop(0)
        if a == "--limits":
            limits = (float(rest.pop(0)), float(rest.pop(0)))
        elif a == "--fresh":
            FRESH = True
        elif "=" in a:
            keys.append(a)
        else:
            fs.append(float(a))
    fs = fs or [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0]
    its, plt, _ = run(gust, "inputs_body_fitted_flat", os.path.join(out, "body_fitted"),
                      [k for k in keys if not k.startswith("terrain.")], "body_fitted_flat")
    if its is None or plt is None:
        print("body-fitted: no steady verdict")
        return 1
    b = fields(plt)
    nz = b["u"].shape[2]
    zb = (np.arange(nz) + 0.5) * DZ
    ub, kb = b["u"][-1, 0, :], b["tke"][-1, 0, :]
    print("body-fitted: iterations %s  flux in %.4f out %.4f  max|w| %.2e  p out [%.2e, %.2e]" % (
        its, b["u"][0].sum() * DZ / b["u"].shape[1], b["u"][-1].sum() * DZ / b["u"].shape[1], np.abs(b["w"]).max(),
        b["p"][-1].min(), b["p"][-1].max()))
    for j in range(3):
        print("   cell %d: u/log %.4f" % (j, ub[j] / (USTAR / KAPPA * math.log((zb[j] + Z0) / Z0))))
    print("%5s %6s %10s %10s %9s %9s %10s %7s %7s %7s %7s %8s %8s" % (
        "F", "iter", "flux in", "flux out", "out-in %", "max|w|", "max|p| out", "cut", "c+1", "c+2", "c+3", "du max%", "dk max%"))
    for f in fs:
        h = (2 + f) * DZ
        its, plt, text = run(gust, "inputs_immersed_flat", os.path.join(out, "f%04.2f" % f), ["terrain.flat_height=%.10g" % h] + keys, "immersed_flat")
        if plt is None:
            print("%5.2f  no plotfile: %s" % (f, text.strip().splitlines()[-1][:150] if text.strip() else "no output"))
            bad.append("F = %g: no plotfile" % f)
            continue
        d = fields(plt)
        n = d["u"].shape[2]
        lo, hi = np.arange(n) * DZ, (np.arange(n) + 1) * DZ
        open_h = np.clip(hi - np.maximum(lo, h), 0, DZ)                   # fluid height of each layer
        zc = np.where(open_h > 0, hi - 0.5 * open_h - h, -1.0)             # fluid centroid above the ground
        # The flux per column, the mean over y (columns differ with side walls).
        fin, fout = (d["u"][0] * open_h).sum() / d["u"].shape[1], (d["u"][-1] * open_h).sum() / d["u"].shape[1]
        kc = int(np.argmax(open_h > 1e-9))
        u = d["u"][-1, 0, :]
        ratio = [u[kc + j] / (USTAR / KAPPA * math.log((zc[kc + j] + Z0) / Z0)) for j in range(4)]
        sel = np.arange(kc + 1, n)
        sel = sel[(zc[sel] > zb[0]) & (zc[sel] < zb[-1])]
        # In ln(z + z0), u as its ratio to the log law: exact between two layers for a logarithmic profile.
        law = lambda z: USTAR / KAPPA * np.log((z + Z0) / Z0)                      # noqa: E731
        xs, xb = np.log(zc[sel] + Z0), np.log(zb + Z0)
        du = np.abs(u[sel] / law(zc[sel]) / np.interp(xs, xb, ub / law(zb)) - 1).max() * 100
        dk = np.abs(d["tke"][-1, 0, sel] / np.interp(xs, xb, kb) - 1).max() * 100
        print("%5.2f %6s %10.4f %10.4f %9.4f %9.2e %10.2e %7.4f %7.4f %7.4f %7.4f %8.3f %8.3f" % (
            f, its, fin, fout, (fout / fin - 1) * 100, np.abs(d["w"]).max(), np.abs(d["p"][-1]).max(), *ratio, du, dk))
        if limits:
            if its is None:
                bad.append("F = %g: no steady verdict" % f)
            if abs(fout / fin - 1) > 1e-6:
                bad.append("F = %g: outflow flux differs from the inflow flux by %.2e of it (limit 1e-6)" % (f, abs(fout / fin - 1)))
            if du > limits[0]:
                bad.append("F = %g: u differs from the body-fitted column by %.3f %% (limit %g)" % (f, du, limits[0]))
            if dk > limits[1]:
                bad.append("F = %g: tke differs from the body-fitted column by %.3f %% (limit %g)" % (f, dk, limits[1]))
    for m in bad:
        print("FAIL " + m)
    if limits and not bad:
        print("PASS")
    return 1 if (limits and bad) else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
