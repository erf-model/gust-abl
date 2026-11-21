#!/usr/bin/env python3
"""Big Southern Butte, 222 degrees: runs against the 3 m sensors and against each other (README.md).

    compare_butte.py OUT.png DATA LABEL=DIR [LABEL=DIR ...]

A run folder holds `deck`, `log`, `plt?????`; one with `zlevels.txt` is an immersed run. DATA is the folder
`windninja2019.py prepare` wrote (masts.csv, bsb.xyz).

Each column's ground comes from the run itself: the lowest node of a terrain-following mesh, the solid height of an
immersed one. Speed at 3 m above that ground: linear between points, the log law (z0 = 0.01 m) below the first."""
import sys, os, csv, glob, re, math
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(HERE + "/../../..")
sys.path.insert(0, REPO + "/extern_erf"); sys.path.insert(0, REPO + "/Tests/Real_Terrain/WindNinja2019"); sys.dont_write_bytecode = True
import numpy as np, erf_plotfile, windninja2019 as wn
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
DATA = ""
AGL, Z0, BORDER = 3.0, 0.01, 1800.0     # BORDER: the flat strip and its ramp in bsb_flat_border.xyz
SIZE = (18690.0, 20160.0)


def iters(d):
    m = re.search(r"steady after (\d+)", open(d + "/log").read())
    return m.group(1) if m else "not steady"


class Run:
    def __init__(self, label, d):
        self.label, self.d = label, d
        self.imm = os.path.exists(d + "/zlevels.txt")
        h, f = erf_plotfile.read_fields(sorted(glob.glob(d + "/plt?????"))[-1], ["u", "v"] + (["building_beta"] if self.imm else ["z"]))
        f = {k: np.asarray(v) for k, v in f.items()}
        self.lo, self.dx = h["prob_lo"], h["dx"]; s = np.hypot(f["u"], f["v"]); n = s.shape; self.n = n
        self.ground = np.zeros(n[:2]); self.s3 = np.zeros(n[:2])
        if self.imm:
            lev = np.array([float(x) for x in open(d + "/zlevels.txt").read().split()]); lo, hi = lev[:-1], lev[1:]
            flat = float(re.search(r"terrain.flat_height\s*=\s*(\S+)", open(d + "/deck").read()).group(1))
        for i in range(n[0]):
            for j in range(n[1]):
                if self.imm:
                    b = f["building_beta"][i, j, :]; g = float(np.sum(b * (hi - lo))); live = b < 1
                    op = np.clip(hi - np.maximum(lo, g), 0, None); z = (hi - 0.5 * op - g)[live]; sc = s[i, j, :][live]; self.ground[i, j] = g - flat
                else:
                    zc = f["z"][i, j, :]
                    g = zc[0] - (zc[1] - zc[0]) / (1 + (zc[2] - zc[1]) / (zc[1] - zc[0]))     # half the first layer, from the stretch ratio
                    z = zc - g; sc = s[i, j, :]; self.ground[i, j] = g
                self.s3[i, j] = sc[0] * math.log((AGL + Z0) / Z0) / math.log((z[0] + Z0) / Z0) if AGL < z[0] else np.interp(AGL, z, sc)

    def at(self, x, y):
        i = min(max(int((x - self.lo[0]) / self.dx[0]), 0), self.n[0] - 1); j = min(max(int((y - self.lo[1]) / self.dx[1]), 0), self.n[1] - 1)
        return self.s3[i, j]


_G = {}


def mast_class(r):
    """windninja2019.mast_class's rule for the Butte, on bsb.xyz (the same ground, shifted by its minimum)."""
    if not _G:
        a = np.loadtxt(DATA + "/bsb.xyz"); k = int(np.argmax(a[:, 2])); _G.update(a=a, zs=a[k, 2], xs=a[k, 0], ys=a[k, 1])
    a = _G["a"]; x, y = float(r["x"]), float(r["y"]); zg = a[int(np.argmin((a[:, 0] - x) ** 2 + (a[:, 1] - y) ** 2)), 2]
    if zg >= _G["zs"] - 150:
        return "r"
    return "w" if (x - _G["xs"]) * math.cos(math.radians(48)) + (y - _G["ys"]) * math.sin(math.radians(48)) < 0 else "l"


def rms(a):
    return float(np.sqrt(np.mean(np.square(a)))) if len(a) else float("nan")


def main(out, data, specs):
    global DATA
    DATA = data
    runs = [Run(*s.split("=", 1)) for s in specs]
    obs = wn.read_masts(DATA + "/masts.csv", "bsb", "222")
    xo = np.array([float(r["x"]) for r in obs]); yo = np.array([float(r["y"]) for r in obs]); so = np.array([float(r["speed"]) for r in obs])
    cls = np.array([mast_class(r) for r in obs])
    edge = np.minimum.reduce([xo, SIZE[0] - xo, yo, SIZE[1] - yo])
    print("%d sensors at %g m; %d within %.0f m of a side face (the flattened strip and ramp)" % (len(obs), AGL, int(np.sum(edge < BORDER)), BORDER))
    print("%-34s %10s %8s %8s | %8s %8s %8s   (speed error, m/s)" % ("run", "iterations", "rms", "bias", "windward", "ridge", "lee"))
    mod = {}
    for r in runs:
        m = np.array([r.at(x, y) for x, y in zip(xo, yo)]); mod[r.label] = m; e = m - so
        print("%-34s %10s %8.2f %8.2f | %8.2f %8.2f %8.2f" % (r.label, iters(r.d), rms(e), e.mean(), rms(e[cls == "w"]), rms(e[cls == "r"]), rms(e[cls == "l"])))
    print("sensors by class: windward %d, ridge %d, lee %d" % tuple(int(np.sum(cls == c)) for c in "wrl"))
    print("\nrun against run: speed at %g m, all columns more than %.0f m from the side faces; and at the sensors" % (AGL, BORDER))
    r0 = runs[0]; X = r0.lo[0] + (np.arange(r0.n[0]) + 0.5) * r0.dx[0]; Y = r0.lo[1] + (np.arange(r0.n[1]) + 0.5) * r0.dx[1]
    XX, YY = np.meshgrid(X, Y, indexing="ij"); inner = np.minimum.reduce([XX, SIZE[0] - XX, YY, SIZE[1] - YY]) > BORDER
    pairs = [(a, b) for n, a in enumerate(runs) for b in runs[n + 1:]]
    for a, b in pairs:
        dS = (b.s3 - a.s3)[inner]
        print("%-30s minus %-30s: rms %.2f, mean %+.2f m/s; |dS| > 0.25 / 0.5 / 1 m/s in %.0f / %.0f / %.0f %% of columns; at the sensors rms %.2f" % (
            b.label, a.label, rms(dS), dS.mean(), 100 * np.mean(np.abs(dS) > 0.25), 100 * np.mean(np.abs(dS) > 0.5), 100 * np.mean(np.abs(dS) > 1), rms(mod[b.label] - mod[a.label])))
    n = len(runs); fig, axes = plt.subplots(2, max(n, len(pairs)), figsize=(5.2 * max(n, len(pairs)), 10.5), squeeze=False)
    for ax, r in zip(axes[0], runs):
        ax.plot([0, 16], [0, 16], color="0.6", lw=1)
        for c, col, lab in (("w", "#1f6f8b", "windward"), ("r", "k", "ridge top"), ("l", "#c0562b", "lee")):
            ax.plot(so[cls == c], mod[r.label][cls == c], "o", color=col, ms=5, label=lab)
        ax.set_xlim(0, 16); ax.set_ylim(0, 16); ax.set_aspect(1); ax.grid(alpha=0.3); ax.set_xlabel("measured speed at 3 m (m/s)"); ax.set_ylabel("run (m/s)")
        ax.set_title("%s\n%s iterations, rms %.2f m/s" % (r.label, iters(r.d), rms(mod[r.label] - so)), loc="left", fontsize=10); ax.legend(fontsize=8, frameon=False)
    for ax in axes[0][n:]: ax.axis("off")
    for ax, (a, b) in zip(axes[1], pairs):
        m = ax.pcolormesh(X / 1000, Y / 1000, (b.s3 - a.s3).T, vmin=-3, vmax=3, cmap="RdBu_r", shading="auto")
        ax.contour(X / 1000, Y / 1000, a.ground.T, levels=np.arange(100, 800, 100), colors="k", linewidths=0.4)
        ax.plot(xo / 1000, yo / 1000, "k.", ms=3); ax.set_aspect(1); ax.set_xlabel("x (km)"); ax.set_ylabel("y (km)")
        ax.set_title("speed at 3 m: %s\nminus %s (m/s)" % (b.label, a.label), loc="left", fontsize=9); fig.colorbar(m, ax=ax, shrink=0.8)
    for ax in axes[1][len(pairs):]: ax.axis("off")
    fig.tight_layout(); fig.savefig(out, dpi=120, bbox_inches="tight"); print(out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
