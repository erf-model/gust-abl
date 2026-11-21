"""compare_dx.py: crest speed-up of the grid-resolution ridge, terrain-following (TF) against immersed, per dx.

    python3 compare_dx.py DX=TF_DIR,IMMERSED_DIR [DX=...] [--plot FILE]

Each folder holds one run's plotfiles (*_plt*, *_turb*; the last of each is read). In every column the wind speed S
and k are read at d = 10, 20, 40 and 80 m above that column's ground: TF, the lowest node; immersed, the solid
depth (the sum of building_beta times each layer's thickness, the layers taken from the cell centres) above the
mesh bottom, with the cut cell at the centroid of its fluid and solid cells left out. Linear interpolation in d;
below a run's first fluid centre the value is that centre's (np.interp holds it). The speed-up is
dS = S / S_far - 1, S_far being the run's own mean over the first sixteenth of the columns. The crest column is
the one centred at x = -dx/2 (the crest lies on a face).

Speeds at 10 m on the steep slopes at coarse dx come from cells in or near the cut band and are less certain than
the crest values.

Printed per run: the peak dS at 10 and 40 m and its x, the lee minimum at 10 m, k at 10 m in the crest column;
then the immersed minus TF dS in the crest column per height and the observed order p of that gap between
successive dx (gap ~ dx^p); then per dx the immersed / TF far-field speed and crest speed at each height (the gap
above leaves out the far-field difference, as each run is divided by its own far field).
"""
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402

AGL = (10.0, 20.0, 40.0, 80.0)
X0, LX = -4096.0, 8192.0


def read(d, imm):
    last = lambda pat: sorted(glob.glob(os.path.join(d, pat)))[-1]  # noqa: E731
    _, f = erf_plotfile.read_fields(last("*_plt?????"), ["u", "v", "z"] + (["building_beta"] if imm else []))
    f = {k: np.asarray(v) for k, v in f.items()}
    _, t = erf_plotfile.read_fields(last("*_turb?????"), ["KE"])
    f["tke"] = np.asarray(t["KE"])
    return f


def column(f, imm, i):
    s = np.hypot(f["u"][i, 0], f["v"][i, 0])
    k = f["tke"][i, 0]
    z = f["z"][i, 0]
    if imm:
        b = f["building_beta"][i, 0]
        faces = [0.0]
        for zc_ in z:
            faces.append(2 * zc_ - faces[-1])
        faces = np.array(faces)
        g = float(np.sum(b * np.diff(faces)))
        zc = np.where(faces[:-1] < g, 0.5 * (g + faces[1:]), z)
        live = b < 1
        return zc[live] - g, s[live], k[live]
    return z - (z[0] - 0.5 * (z[1] - z[0])), s, k


def speedup(d, imm):
    f = read(d, imm)
    n = f["u"].shape[0]
    x = X0 + (np.arange(n) + 0.5) * LX / n
    S = np.zeros((n, len(AGL)))
    K = np.zeros_like(S)
    for i in range(n):
        z, s, k = column(f, imm, i)
        S[i], K[i] = np.interp(AGL, z, s), np.interp(AGL, z, k)
    far = S[: max(1, n // 16)].mean(0)
    return x, S / far - 1, K, far, S


def main(argv):
    plot = None
    if "--plot" in argv:
        i = argv.index("--plot")
        if i + 1 >= len(argv):
            sys.exit("compare_dx.py: --plot needs a file name")
        plot = argv[i + 1]
        argv = argv[:i] + argv[i + 2:]
    if not argv or any(a.startswith("-") for a in argv):
        sys.exit(__doc__)
    runs = []
    for a in argv:
        dx, _, dirs = a.partition("=")
        tf, _, im = dirs.partition(",")
        try:
            float(dx)
        except ValueError:
            sys.exit(f"compare_dx.py: {a!r}: DX must be a number\n{__doc__}")
        if not tf or not im:
            sys.exit(f"compare_dx.py: {a!r} is not DX=TF_DIR,IMMERSED_DIR")
        for d in (tf, im):
            if not glob.glob(os.path.join(d, "*_plt?????")) or not glob.glob(os.path.join(d, "*_turb?????")):
                sys.exit(f"compare_dx.py: no *_plt?????/*_turb????? plotfiles in {d}")
        runs.append((dx, tf, im))
    res, raw = {}, {}
    print(f"{'dx':>4} {'run':9} | {'peak dS10':>9} {'x':>6} | {'peak dS40':>9} {'x':>6} | {'lee min dS10':>12} {'x':>5} | k10 crest")
    for dx, tf, im in runs:
        for lab, d, imm in (("TF", tf, False), ("immersed", im, True)):
            x, s, k, far, S = speedup(d, imm)
            res[dx, lab] = (x, s)
            raw[dx, lab] = (far, S[np.argmin(abs(x))])
            i10, i40, ic = np.argmax(s[:, 0]), np.argmax(s[:, 2]), np.argmin(abs(x))
            il = np.argmin(np.where(x > 0, s[:, 0], 9))
            print(f"{dx:>4} {lab:9} | {s[i10, 0]:9.3f} {x[i10]:6.0f} | {s[i40, 2]:9.3f} {x[i40]:6.0f} | "
                  f"{s[il, 0]:12.3f} {x[il]:5.0f} | {k[ic, 0]:.2f}")
    print("\nimmersed minus TF dS in the crest column; p: observed order between successive dx")
    for a, z in enumerate(AGL):
        gap = []
        for dx, *_ in runs:
            (xt, st), (xi, si) = res[dx, "TF"], res[dx, "immersed"]
            gap.append(si[np.argmin(abs(xi)), a] - st[np.argmin(abs(xt)), a])
        dxs = [float(r[0]) for r in runs]
        p = [np.log(gap[i + 1] / gap[i]) / np.log(dxs[i + 1] / dxs[i]) for i in range(len(gap) - 1)
             if gap[i] * gap[i + 1] > 0]
        print(f"  d {int(z):3d} m: " + "  ".join(f"dx {r[0]}: {g:+.3f}" for r, g in zip(runs, gap))
              + ("   p = " + ", ".join(f"{q:.2f}" for q in p) if p else ""))
    print("\nimmersed / TF speed: far field | crest column, at " + ", ".join(f"{int(z)}" for z in AGL) + " m")
    for dx, *_ in runs:
        (ft, ct), (fi, ci) = raw[dx, "TF"], raw[dx, "immersed"]
        print(f"  dx {dx:>4}: " + " ".join(f"{q:.3f}" for q in fi / ft) + " | " + " ".join(f"{q:.3f}" for q in ci / ct))
    if plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
        cols = plt.rcParams["axes.prop_cycle"].by_key()["color"]
        for a in range(3):
            for c, (dx, *_) in enumerate(runs):
                for lab, sty, lw in (("TF", "--", 1.2), ("immersed", "-", 2.0)):
                    x, s = res[dx, lab]
                    ax[a].plot(x, s[:, a], sty, color=cols[c % len(cols)], lw=lw, label=f"{lab} dx {dx}")
            ax[a].axhline(0, color="0.7", lw=0.6)
            ax[a].set_ylabel(f"dS at {int(AGL[a])} m")
            ax[a].set_xlim(-1600, 1600)
        ax[0].legend(ncol=len(runs), fontsize=8)
        ax[-1].set_xlabel("x (m), crest at 0")
        fig.savefig(plot, dpi=140)
        print("wrote", plot)


if __name__ == "__main__":
    main(sys.argv[1:])
