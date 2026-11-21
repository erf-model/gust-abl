#!/usr/bin/env python3
"""Flat ground with an inflow and an outflow face on layers with a fine band, immersed against body-fitted.

    sweep_band.py GUST OUTDIR [key=value ...]

The immersed mesh has coarse layers at z = 0, a fine band and stretched layers above it (band_levels). The
ground is placed in one layer of a zone at solid fraction F of that layer (0, 0.1, ..., 0.9, 0.99, 1):

    uniform     67 layers of 15.625 m (the mesh of sweep_flat.py), ground in the third
    coarse      ground in the second 16 m layer, the band 40 m and more above it
    transition  ground in the 12 m layer between the coarse layers and the band
    band        ground in a 4 m layer inside the band
    band_top    ground in the last 4 m layer of the band, the stretched layers above it

The body-fitted reference of a zone has the immersed mesh's layers from the bottom face of that layer upward.
key=value arguments replace or add deck keys in every run (turb.shear_stencil=log_law, ...). Printed for each
run: iterations, the mass flux through the inflow and outflow columns, max |w|, max |p| in the outflow column,
the cut cell's speed over the log law, and the largest difference of u and tke from the body-fitted outflow
column at equal height above the ground. The reference is interpolated in ln(z + z0), u as its ratio to the
log law, so that the interpolation between two layers is exact for a logarithmic profile.
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

Z0, KAPPA, USTAR = 0.1, 0.41, 0.35
FS = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99, 1.0]


def band_levels():
    """Coarse 16 m layers to 48 m, 12, 8 and 6 m, ten 4 m layers (74 to 114 m), then x 1.15 up to 40 m."""
    lev = [0.0, 16.0, 32.0, 48.0, 60.0, 68.0, 74.0]
    lev += [74.0 + 4.0 * (i + 1) for i in range(10)]
    dz = 4.0
    while lev[-1] < 1114.0:
        dz = min(dz * 1.15, 40.0)
        lev.append(lev[-1] + dz)
    return np.array(lev)


ZONES = {  # name: (levels, index of the layer that holds the ground)
    "uniform": (np.arange(68) * 15.625, 2),
    "coarse": (band_levels(), 1),
    "transition": (band_levels(), 3),
    "band": (band_levels(), 8),
    "band_top": (band_levels(), 15),
}


def loglaw(z):
    return USTAR / KAPPA * np.log((z + Z0) / Z0)


def run(gust, cwd, base, lev, keys, prefix):
    os.makedirs(cwd, exist_ok=True)
    keys = dict(keys)
    keys.update({"mesh.n_cell": "64 4 %d" % (len(lev) - 1), "mesh.prob_hi": "2560. 160. %.10g" % lev[-1],
                 "mesh.zlevels_file": "zlevels.txt"})
    with open(os.path.join(cwd, "zlevels.txt"), "w") as f:
        f.write("".join("%.10g\n" % z for z in lev))
    with open(os.path.join(cwd, "deck"), "w") as f:
        f.write("".join(ln for ln in open(os.path.join(HERE, base)) if ln.split("=")[0].strip() not in keys))
        f.write("".join("%s = %s\n" % kv for kv in keys.items()))
    log = os.path.join(cwd, "log")
    if not os.path.exists(os.path.join(cwd, "done")):
        with open(log, "w") as f:
            code = subprocess.run([gust, "deck"], cwd=cwd, stdout=f, stderr=subprocess.STDOUT).returncode
        open(os.path.join(cwd, "done"), "w").write(str(code))
    text = open(log).read()
    m = re.search(r"gust: steady after (\d+)", text)
    wrote = [ln.split()[-1] for ln in text.splitlines() if ln.startswith("gust: wrote " + prefix + "_plt")]
    if not wrote:
        return None, None, (text.strip().splitlines() or ["no output"])[-1][:160]
    _, d = erf_plotfile.read_fields(os.path.join(cwd, wrote[-1]), ["u", "w", "p", "tke"])
    return (int(m.group(1)) if m else None), {k: np.asarray(v) for k, v in d.items()}, ""


def main(argv):
    gust, out = os.path.abspath(argv[1]), os.path.abspath(argv[2])
    keys = dict(a.split("=", 1) for a in argv[3:])
    for name, (lev, c) in ZONES.items():
        blev = lev[c:] - lev[c]
        its, b, msg = run(gust, os.path.join(out, name, "body_fitted"), "inputs_body_fitted_flat", blev, keys, "body_fitted_flat")
        if b is None:
            print("%s: body-fitted: %s" % (name, msg))
            continue
        zb, dzb = 0.5 * (blev[1:] + blev[:-1]), np.diff(blev)
        ub, kb = b["u"][-1, 0, :], b["tke"][-1, 0, :]
        rb = ub / loglaw(zb)
        print("%s (ground layer %.4g m, layers above %s m)" % (name, lev[c + 1] - lev[c], " ".join("%.3g" % v for v in np.diff(lev)[c + 1:c + 5])))
        print("  body-fitted: iter %s  flux in %.3f out %.3f  max|w| %.2e  max|p| out %.2e  u/log first cells %s" % (
            its, (b["u"][0, 0, :] * dzb).sum(), (ub * dzb).sum(), np.abs(b["w"]).max(), np.abs(b["p"][-1]).max(), np.round(rb[:3], 4)))
        print("  %5s %5s %10s %10s %9s %9s %10s %8s %8s %8s  %s" % (
            "F", "iter", "flux in", "flux out", "out-in %", "max|w|", "max|p| out", "cut/log", "du max%", "dk max%", "du % in the four cells above the cut cell"))
        for f in FS:
            h = lev[c] + f * (lev[c + 1] - lev[c])
            k2 = dict(keys)
            k2["terrain.flat_height"] = "%.10g" % h
            its, d, msg = run(gust, os.path.join(out, name, "f%04.2f" % f), "inputs_immersed_flat", lev, k2, "immersed_flat")
            if d is None:
                print("  %5.2f  no plotfile: %s" % (f, msg))
                continue
            lo, hi = lev[:-1], lev[1:]
            open_h = np.clip(hi - np.maximum(lo, h), 0, None)                          # fluid height of each layer
            zc = np.where(open_h > 0, hi - 0.5 * open_h - h, -1.0)                      # fluid centroid above the ground
            kc = int(np.argmax(open_h > 1e-9 * (hi - lo)))
            u, tke = d["u"][-1, 0, :], d["tke"][-1, 0, :]
            fin, fout = (d["u"][0, 0, :] * open_h).sum(), (u * open_h).sum()
            sel = np.arange(kc + 1, len(lo))
            sel = sel[(zc[sel] > zb[0]) & (zc[sel] < zb[-1])]
            xs, xb = np.log(zc[sel] + Z0), np.log(zb + Z0)
            du = (u[sel] / loglaw(zc[sel]) / np.interp(xs, xb, rb) - 1) * 100
            dk = (tke[sel] / np.interp(xs, xb, kb) - 1) * 100
            print("  %5.2f %5s %10.4f %10.4f %9.4f %9.2e %10.2e %8.4f %8.3f %8.3f  %s" % (
                f, its, fin, fout, (fout / fin - 1) * 100, np.abs(d["w"]).max(), np.abs(d["p"][-1]).max(),
                u[kc] / loglaw(zc[kc]), np.abs(du).max(), np.abs(dk).max(), np.round(du[:4], 2)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
