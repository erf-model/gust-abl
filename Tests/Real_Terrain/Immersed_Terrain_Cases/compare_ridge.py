#!/usr/bin/env python3
"""Ishihara-Zhou ridge, body-fitted against immersed: the near-surface line, U and k at the tunnel stations
(README.md; the measures of ../IshiharaZhou/ishihara_zhou.py compare).

    compare_ridge.py SURFACE BODY_FITTED_DIR IMMERSED_DIR OUT.png DATA      (SURFACE = smooth | rough)

A run folder holds `log`, `plt?????` and, for the immersed run, `zlevels.txt`. DATA is the folder
`ishihara_zhou.py prepare` wrote (stations_smooth.csv, stations_rough.csv)."""
import sys, os, csv, glob, re, math
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(HERE + "/../../..")
sys.path.insert(0, REPO + "/extern_erf"); sys.dont_write_bytecode = True
import numpy as np, erf_plotfile
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
H, L, XC, FLAT = 40.0, 100.0, 800.0, 0.75
SURF = {"smooth": dict(z0=0.01, zline=5.0, reattach=(150.0, 200.0), les=180.0), "rough": dict(z0=0.3, zline=7.0, reattach=(250.0, 300.0), les=239.0)}
CB, CI = "#1f6f8b", "#c0562b"


def ground(x):
    s = abs(x - XC)
    return H * math.cos(0.5 * math.pi * s / L) ** 2 if s <= L else 0.0


def iters(d):
    m = re.search(r"steady after (\d+)", open(d + "/log").read())
    return m.group(1) if m else "not steady"


class Run:
    def __init__(self, d, imm):
        self.imm = imm
        h, f = erf_plotfile.read_fields(sorted(glob.glob(d + "/plt?????"))[-1], ["u", "tke"] + (["building_beta"] if imm else ["z"]))
        self.f = {k: np.asarray(v) for k, v in f.items()}
        self.lo, self.dx = h["prob_lo"], h["dx"]; self.n = self.f["u"].shape; self.j = self.n[1] // 2
        if imm: self.lev = np.array([float(x) for x in open(d + "/zlevels.txt").read().split()])

    def col(self, i):
        """Heights above the local ground, u and k of column i (fluid cells; the cut cell at its open part's centre)."""
        g = ground(self.lo[0] + (i + 0.5) * self.dx[0]); u, k = self.f["u"][i, self.j, :], self.f["tke"][i, self.j, :]
        if not self.imm:
            return g, self.f["z"][i, self.j, :] - g, u, k
        gi = g + FLAT; lo, hi = self.lev[:-1], self.lev[1:]
        op = np.clip(hi - np.maximum(lo, gi), 0, None); live = self.f["building_beta"][i, self.j, :] < 1
        return g, (hi - 0.5 * op - gi)[live], u[live], k[live]

    def at(self, x, z, z0, above_ground):
        """u and k at x: z above the ground, or z absolute (the tunnel stations), linear in x between two columns."""
        fx = (x - self.lo[0]) / self.dx[0] - 0.5; i0 = int(math.floor(fx)); t = fx - i0; out = np.zeros(2)
        for i, w in ((i0, 1 - t), (i0 + 1, t)):
            g, zr, u, k = self.col(i); zt = z if above_ground else z - g
            if zt < zr[0]:
                f = math.log((max(zt, 0) + z0) / z0) / math.log((zr[0] + z0) / z0) if above_ground else 1.0
                out += w * np.array([u[0] * f, k[0]])
            else:
                out += w * np.array([np.interp(zt, zr, u), np.interp(zt, zr, k)])
        return out

    def profile(self, x, ztop=200.0):
        i = int((x - self.lo[0]) / self.dx[0]); g, zr, u, k = self.col(i); m = zr <= ztop
        return zr[m] + g, u[m], k[m]


def main(surface, dbf, dim, out, IZ):
    s = SURF[surface]; runs = [(Run(dbf, False), CB, "body-fitted, same settings", dbf), (Run(dim, True), CI, "immersed", dim)]
    st = list(csv.DictReader(open(IZ + "/stations_%s.csv" % surface))); xst = sorted({float(r["x"]) for r in st})
    fig = plt.figure(figsize=(13, 10)); gs = fig.add_gridspec(2, 3, hspace=0.3)
    a0 = fig.add_subplot(gs[0, 0]); a1 = fig.add_subplot(gs[0, 1:]); a2 = fig.add_subplot(gs[1, 1:]); a3 = fig.add_subplot(gs[1, 0])
    xs = np.arange(-150.0, 401.0, 1.0)
    for run, c, lab, d in runs:
        ul = np.array([run.at(XC + x, s["zline"], s["z0"], True)[0] for x in xs])
        kl = np.array([run.at(XC + x, s["zline"], s["z0"], True)[1] for x in xs])
        rev = xs[(ul < 0) & (xs >= -100) & (xs <= 350)]
        m = np.array([run.at(XC + float(r["x"]), float(r["z"]), s["z0"], False) for r in st])
        eu = m[:, 0] - np.array([float(r["U"]) for r in st]); ek = m[:, 1] - np.array([float(r["k"]) for r in st])
        print("%-28s %s: %s iterations; reversed flow at %g m: %s; U rms %.3f m/s, k rms %.3f m2/s2 over %d points" % (
            lab, surface, iters(d), s["zline"], ("%.0f to %.0f m" % (rev.min(), rev.max())) if rev.size else "none",
            np.sqrt(np.mean(eu ** 2)), np.sqrt(np.mean(ek ** 2)), len(st)))
        a0.plot(xs, ul, color=c, lw=2, label="%s (%s iterations)" % (lab, iters(d))); a3.plot(xs, kl, color=c, lw=2)
        for n, x in enumerate(xst):
            z, u, k = run.profile(XC + x)
            a1.plot(x + 5 * u, z, color=c, lw=1.3, label=lab if n == 0 else None); a2.plot(x + 25 * k, z, color=c, lw=1.3)
    ln = [r for r in st if abs(float(r["z"]) - s["zline"]) < 1e-6]
    a0.plot([float(r["x"]) for r in ln], [float(r["U"]) for r in ln], "ko", ms=5, label="tunnel, z = %g m" % s["zline"])
    a3.plot([float(r["x"]) for r in ln], [float(r["k"]) for r in ln], "ko", ms=5)
    a0.axvspan(*s["reattach"], color="k", alpha=0.12, label="tunnel reattachment"); a0.axvline(s["les"], color="grey", ls=":", label="LES reattachment")
    a0.axhline(0, color="0.7", lw=0.8); a0.set_xlabel("x from the crest (m)"); a0.set_ylabel("u at %g m above the ground (m/s)" % s["zline"])
    a0.set_title("Ishihara-Zhou ridge, %s: near-surface line" % surface, loc="left", fontsize=10); a0.legend(fontsize=7, frameon=False); a0.grid(alpha=0.3)
    a3.set_xlabel("x from the crest (m)"); a3.set_ylabel("k at %g m (m$^2$/s$^2$)" % s["zline"]); a3.set_title("k on the same line (station heights are absolute)", loc="left", fontsize=10); a3.grid(alpha=0.3)
    gx = np.arange(-200.0, 401.0, 2.0)
    for a, key, sc, ttl, xl in ((a1, "U", 5, "U at the tunnel stations", "x from the crest (m); profiles: 50 m per 10 m/s"),
                                (a2, "k", 25, "k at the tunnel stations", "x from the crest (m); profiles: 50 m per 2 m$^2$/s$^2$")):
        for n, x in enumerate(xst):
            pts = sorted((float(r["z"]), float(r[key])) for r in st if float(r["x"]) == x)
            a.plot([x + sc * p[1] for p in pts], [p[0] for p in pts], "ko", ms=3, label="tunnel" if n == 0 else None); a.axvline(x, color="0.85", lw=0.6)
        a.fill_between(gx, [ground(XC + x) for x in gx], 0, color="0.75"); a.set_ylim(0, 200); a.set_xlabel(xl); a.set_ylabel("z (m)"); a.set_title(ttl, loc="left", fontsize=10)
    a1.legend(fontsize=8, frameon=False, loc="upper right")
    fig.savefig(out, dpi=120, bbox_inches="tight"); print(out)


if __name__ == "__main__":
    main(*sys.argv[1:6])
