#!/usr/bin/env python3
"""A coarse pair of one case, body-fitted against immersed, as a test (README.md, CTest).

    check_coarse_pair.py --body-fitted DECK --n-cell "NX NY" --agl H --z0 Z0 --border B
                         --rms R --top T [--key KEY=VALUE ...] -- [launcher ...] GUST IMMERSED_DECK [key=value ...]

Runs the immersed deck in ./immersed and the body-fitted deck in ./body_fitted with NX x NY columns (each deck
keeps its own layers) and the --key arguments (terrain.file=..., mesh.zlevels_file=... for the immersed deck
only), then reads the speed at H above the ground (m) in every column. A column's ground is the lowest node
of the terrain-following mesh, the solid height of the immersed one. Speed at H: linear between cell values,
the log law with Z0 below the first. Exit 1 unless both runs reach the steady verdict, the rms difference of the
two speeds over the columns more than B from the side faces is below R, and the difference in the column with the
highest ground is below T (both in percent of the mean body-fitted speed at H in those columns).
"""
import glob
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


def deck_value(deck, key):
    """The last value of key in a deck."""
    val = None
    for line in open(deck):
        left, _, right = line.partition("=")
        if left.strip() == key:
            val = right.split("#")[0].strip()
    return val


def run(cmd, cwd):
    os.makedirs(cwd, exist_ok=True)
    for old in glob.glob(os.path.join(cwd, "plt?????")):
        subprocess.run(["rm", "-rf", old], check=True)
    out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    text = out.stdout + out.stderr
    open(os.path.join(cwd, "log"), "w").write(text)
    sys.stdout.write(text[-800:])
    m = re.search(r"gust: steady after (\d+)", text)
    plts = sorted(glob.glob(os.path.join(cwd, "plt?????")))
    return (int(m.group(1)) if m else None), (plts[-1] if plts else None)


def speed_at(plt, agl, z0, levels, flat):
    """Speed at agl above the ground and the ground's height, per column."""
    names = ["u", "v"] + (["building_beta"] if levels is not None else ["z"])
    h, f = erf_plotfile.read_fields(plt, names)
    f = {k: np.asarray(v) for k, v in f.items()}
    s = np.hypot(f["u"], f["v"])
    n = s.shape
    ground, out = np.zeros(n[:2]), np.zeros(n[:2])
    for i in range(n[0]):
        for j in range(n[1]):
            if levels is not None:
                lo, hi = levels[:-1], levels[1:]
                b = f["building_beta"][i, j, :]
                g = float(np.sum(b * (hi - lo)))
                live = b < 1
                op = np.clip(hi - np.maximum(lo, g), 0, None)
                z, sc = (hi - 0.5 * op - g)[live], s[i, j, :][live]      # fluid centroids above the ground
                ground[i, j] = g - flat
            else:
                zc = f["z"][i, j, :]
                g = zc[0] - (zc[1] - zc[0]) / (1 + (zc[2] - zc[1]) / (zc[1] - zc[0]))   # half the first layer
                z, sc = zc - g, s[i, j, :]
                ground[i, j] = g
            if agl < z[0]:
                out[i, j] = sc[0] * math.log((agl + z0) / z0) / math.log((z[0] + z0) / z0)
            else:
                out[i, j] = np.interp(agl, z, sc)
    return h, ground, out


def main(argv):
    if "--" not in argv:
        print(__doc__)
        return 2
    sep = argv.index("--")
    opts, cmd = argv[:sep], argv[sep + 1:]
    o, keys = {}, []
    while opts:
        a = opts.pop(0)
        if a == "--key":
            keys.append(opts.pop(0))
        else:
            o[a.lstrip("-")] = opts.pop(0)
    imm_deck = next(a for a in cmd if os.path.isfile(a) and os.path.basename(a).startswith("inputs_"))
    at = cmd.index(imm_deck)
    bf_deck = o["body-fitted"]
    agl, z0, border = float(o["agl"]), float(o["z0"]), float(o["border"])
    cwd = os.getcwd()

    def one(deck, sub, immersed):
        """Run a copy of the deck with the changed keys replaced (the reader refuses a key given twice)."""
        nz = deck_value(deck, "mesh.n_cell").split()[2]
        extra = ["mesh.n_cell=%s %s" % (o["n-cell"], nz)] + [k for k in keys if immersed or not k.startswith("mesh.zlevels_file")]
        rest = []
        for a in cmd[at + 1:]:
            (extra if "=" in a and deck_value(deck, a.split("=", 1)[0]) is not None else rest).append(a)
        new = dict(e.split("=", 1) for e in extra)
        lines = [ln for ln in open(deck) if ln.partition("=")[0].strip() not in new]
        os.makedirs(os.path.join(cwd, sub), exist_ok=True)
        with open(os.path.join(cwd, sub, "deck"), "w") as fh:
            fh.writelines(lines)
            fh.write("".join("%s = %s\n" % kv for kv in new.items()))
        return run(cmd[:at] + ["deck"] + rest, os.path.join(cwd, sub))

    its_i, plt_i = one(imm_deck, "immersed", True)
    its_b, plt_b = one(bf_deck, "body_fitted", False)
    print("iterations: body-fitted %s, immersed %s" % (its_b, its_i))
    if its_i is None or its_b is None or plt_i is None or plt_b is None:
        print("FAIL a run did not reach the steady verdict")
        return 1
    zfile = next(k.split("=", 1)[1] for k in keys if k.startswith("mesh.zlevels_file="))
    levels = np.array([float(x) for x in open(zfile).read().split()])
    hb, gb, sb = speed_at(plt_b, agl, z0, None, 0.0)
    _, gi, si = speed_at(plt_i, agl, z0, levels, float(deck_value(imm_deck, "terrain.flat_height")))
    n = sb.shape
    x = hb["prob_lo"][0] + (np.arange(n[0]) + 0.5) * hb["dx"][0]
    y = hb["prob_lo"][1] + (np.arange(n[1]) + 0.5) * hb["dx"][1]
    xx, yy = np.meshgrid(x, y, indexing="ij")
    inner = np.minimum.reduce([xx - x[0], x[-1] - xx, yy - y[0], y[-1] - yy]) > border
    ref = float(sb[inner].mean())
    rms = 100 * float(np.sqrt(np.mean((si - sb)[inner] ** 2))) / ref
    top = np.unravel_index(int(np.argmax(np.where(inner, gb, -1e30))), n)
    dtop = 100 * float(si[top] - sb[top]) / ref
    print("ground: body-fitted %.1f to %.1f m, immersed %.1f to %.1f m" % (gb.min(), gb.max(), gi.min(), gi.max()))
    print("speed at %g m, %d columns more than %g m from the side faces: mean body-fitted %.3f m/s" % (agl, int(inner.sum()), border, ref))
    print("rms difference %.2f %% (limit %s); highest ground (%.0f m): body-fitted %.3f, immersed %.3f m/s, %+.2f %% (limit %s)" % (
        rms, o["rms"], gb[top], sb[top], si[top], dtop, o["top"]))
    bad = []
    if not rms < float(o["rms"]):
        bad.append("the rms difference")
    if not abs(dtop) < float(o["top"]):
        bad.append("the difference on the highest ground")
    print("FAIL " + ", ".join(bad) if bad else "PASS")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
