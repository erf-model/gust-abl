#!/usr/bin/env python3
"""Body-fitted against immersed, one mast case: speed-up along the mast lines, sections of speed and k, mast
profiles (README.md).

    compare_masts.py CASE BODY_FITTED_DIR IMMERSED_DIR OUT.png DATA [BOLUND_ZIP]

CASE = bolund_270, bolund_239, bolund_90 or askervein_210. A run folder holds the run's deck as `deck`, its
output as `log`, its plotfiles `plt?????` and, for the immersed run, the levels as `zlevels.txt`. DATA is the
folder `windninja2019.py prepare` wrote (masts.csv, bolund.xyz). BOLUND_ZIP (optional) is Bolund_Measurements.zip
of the Bolund blind comparison: with it the mast panels show the measured profiles of speed and k."""
import sys, os, csv, glob, re, subprocess, io, zipfile, math
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(HERE + "/../../..")
sys.path.insert(0, REPO + "/extern_erf"); sys.path.insert(0, REPO + "/Tests/Real_Terrain/WindNinja2019"); sys.path.insert(0, REPO + "/tools")
sys.dont_write_bytecode = True
import numpy as np, erf_plotfile, windninja2019 as wn
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
W = ""       # DATA
ZIP = None   # BOLUND_ZIP (Bechmann et al. 2011)
REF = {"270": "M0", "239": "M0", "90": "M9"}
CB, CI = "#1f6f8b", "#c0562b"
# flat: terrain.flat_height of the immersed deck; section: the line of the sections and its d range; zshow: height
# shown above the highest ground; ztop: top of the profiles sampled; prof: top of the mast profile panels.
SITES = {
    "bolund": dict(name="Bolund", xyz="bolund.xyz", flat=0.75, z0=0.01, agl=5, drange=(-150, 200), zshow=33, ztop=40, prof=20,
                   smax=16, kmax=7, dmax=4, kn=14),
    "askervein": dict(name="Askervein", erf=REPO + "/extern_erf/askervein.txt", flat=1.5, z0=0.03, agl=10, section="A", drange=(-1100, 1000),
                      zshow=150, ztop=200, prof=60, smax=18, kmax=4, dmax=3, kn=6, uref=8.9),
}


def measured(case):
    """Sonic masts of the Bolund blind comparison: mast, height, x, y, speed and k (as bolund_masts.py reads them)."""
    rows = list(csv.reader(io.StringIO(zipfile.ZipFile(ZIP).read("Measurements/Dir_%s.dat" % case).decode()), delimiter=" ", skipinitialspace=True))
    out = []
    for r in rows[1:]:
        if not r or not r[0].endswith("S"): continue
        d = dict(zip(rows[0], r)); us = float(d["u*[m/s]"])
        out.append(dict(mast=d["ID"][:2], h=int(d["ID"][3:5]), x=float(d["x[m]"]), y=float(d["y[m]"]), S=float(d["vel/u*"]) * us, k=float(d["tke/u*^2"]) * us * us))
    return out


def last_plt(d):
    return sorted(glob.glob(d + "/plt?????"))[-1]


def iters(d):
    m = re.search(r"steady after (\d+)", open(d + "/log").read())
    return m.group(1) if m else "not steady"


def sample(d, site, case, S, imm):
    out = d + "/s.csv"
    cmd = [sys.executable, HERE + "/sample_masts.py", last_plt(d), "--site", site, "--case", case, "--out", out, "--z0", str(S["z0"]),
           "--masts", W + "/masts.csv"]
    cmd += ["--terrain", W + "/" + S["xyz"]] if "xyz" in S else ["--erf_terrain", S["erf"]]
    if imm: cmd += ["--flat", str(S["flat"]), "--zlevels", d + "/zlevels.txt"]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, check=True)
    return list(csv.DictReader(open(out)))


class Run:
    """Columns of a run as (ground, height above the ground, speed, k, along-wind velocity), the ground in the terrain's frame."""
    def __init__(self, d, imm, ground, flat, ang):
        self.imm, self.flat = imm, flat
        names = ["u", "v", "tke"] + (["building_beta"] if imm else ["z"])
        h, f = erf_plotfile.read_fields(last_plt(d), names)
        self.f = {k: np.asarray(v) for k, v in f.items()}
        self.lo, self.dx = h["prob_lo"], h["dx"]
        self.n = self.f["u"].shape
        self.ground, self.c, self.s = ground, math.cos(ang), math.sin(ang)
        if imm: self.lev = np.array([float(x) for x in open(d + "/zlevels.txt").read().split()])

    def column(self, x, y):
        i = min(max(int((x - self.lo[0]) / self.dx[0]), 0), self.n[0] - 1)
        j = min(max(int((y - self.lo[1]) / self.dx[1]), 0), self.n[1] - 1)
        xc, yc = self.lo[0] + (i + 0.5) * self.dx[0], self.lo[1] + (j + 0.5) * self.dx[1]
        g = 0.25 * sum(self.ground(xc + sx * self.dx[0] / 2, yc + sy * self.dx[1] / 2) for sx in (-1, 1) for sy in (-1, 1))
        u, v = self.f["u"][i, j, :], self.f["v"][i, j, :]
        s = np.hypot(u, v); k = self.f["tke"][i, j, :]; a = self.c * u + self.s * v
        if not self.imm:
            return g, self.f["z"][i, j, :] - g, s, k, a
        gi = g + self.flat; lo, hi = self.lev[:-1], self.lev[1:]
        op = np.clip(hi - np.maximum(lo, gi), 0, None); live = self.f["building_beta"][i, j, :] < 1
        return g, (hi - 0.5 * op - gi)[live], s[live], k[live], a[live]


def main(case, dbf, dim, out, data, bolund_zip=None):
    global W, ZIP
    W, ZIP = data, bolund_zip
    site, cs = case.split("_"); S = SITES[site]
    if "xyz" in S:
        ground = wn.Idw(W + "/" + S["xyz"])
    else:
        grid = wn.al.read_terrain(S["erf"]); ground = lambda x, y: wn.al.ground(grid, x, y)
    ang = math.radians(270.0 - float(cs))          # the inflow's direction of travel, counter-clockwise from x
    obs = [r for r in csv.DictReader(open(W + "/masts.csv")) if r["site"] == site and r["case"] == cs]
    S = dict(S, flat=float(re.search(r"terrain.flat_height\s*=\s*(\S+)", open(dim + "/deck").read()).group(1)))   # the run's own
    sb, si = sample(dbf, site, cs, S, False), sample(dim, site, cs, S, True)
    names = sorted({r["line"] for r in obs}); sec = S.get("section", names[-1])
    rb, ri = Run(dbf, False, ground, S["flat"], ang), Run(dim, True, ground, S["flat"], ang)
    nl = len(names)
    fig = plt.figure(figsize=(13, 15.5 + 3 * nl))
    gs = fig.add_gridspec(nl + 5, 4, height_ratios=[1.1] * nl + [1, 1, 1, 1.2, 1.2], hspace=0.42, wspace=0.3)
    print("%-8s %8s | %8s %8s %8s   speed-up at %g m" % ("mast", "d", "obs", "bf", "imm", S["agl"]))
    for n, name in enumerate(names):
        ax = fig.add_subplot(gs[n, :]); o = [r for r in obs if r["line"] == name]
        for rows, c, lab, d in ((sb, CB, "body-fitted, same settings", dbf), (si, CI, "immersed", dim)):
            l = [r for r in rows if not r["mast"] and r["line"] == name]
            ax.plot([float(r["d"]) for r in l], [float(r["speedup"]) for r in l], color=c, lw=2, label="%s (%s iterations)" % (lab, iters(d)))
        ax.plot([float(r["d"]) for r in o], [float(r["speedup"]) for r in o], "ko", ms=6, label="masts, %g m" % S["agl"])
        if site == "bolund":
            for r in o: ax.annotate(r["mast"][:2], (float(r["d"]), float(r["speedup"])), textcoords="offset points", xytext=(5, 6), fontsize=9)
        mb = {r["mast"]: float(r["speedup"]) for r in sb if r["mast"]}; mi = {r["mast"]: float(r["speedup"]) for r in si if r["mast"]}
        eb, ei = [], []
        for r in o:
            if r["mast"] in mb and r["mast"] in mi:
                print("%-8s %8.0f | %8.3f %8.3f %8.3f" % (r["mast"], float(r["d"]), float(r["speedup"]), mb[r["mast"]], mi[r["mast"]]))
                eb.append(mb[r["mast"]] - float(r["speedup"])); ei.append(mi[r["mast"]] - float(r["speedup"]))
        print("line %s: rms error of the speed-up, body-fitted %.3f, immersed %.3f (%d masts)" % (name, np.sqrt(np.mean(np.square(eb))), np.sqrt(np.mean(np.square(ei))), len(eb)))
        if name == sec or nl == 1: ax.set_xlim(*S["drange"])
        ax.set_xlabel("distance along the line (m)"); ax.set_ylabel("speed-up at %g m" % S["agl"])
        ax.set_title("%s %s°: speed-up along line %s" % (S["name"], cs, name), loc="left"); ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=9)
    line = [r for r in sb if not r["mast"] and r["line"] == sec]
    d = np.array([float(r["d"]) for r in line]); X = [(float(r["x"]), float(r["y"])) for r in line]
    keep = (d > S["drange"][0]) & (d < S["drange"][1])
    zz = np.linspace(0.25, S["ztop"], 120)
    cols = [(rb.column(*p), ri.column(*p)) for p, k in zip(X, keep) if k]; dk = d[keep]
    g = np.array([c[0][0] for c in cols]); top = g.max() + S["zshow"]

    def field(which, idx):
        return np.array([np.interp(zz, c[which][1], c[which][idx], left=np.nan) for c in cols]).T
    DD, ZZ = np.meshgrid(dk, zz); ZA = ZZ + g[None, :]

    def deco(a, title):
        a.fill_between(dk, 0, g, color="0.75"); a.set_ylim(0, top); a.set_xlim(dk[0], dk[-1]); a.set_title(title, loc="left", fontsize=10)
    for row, (idx, nm, vmax, cmap) in enumerate([(2, "speed (m/s)", S["smax"], "viridis"), (3, "k (m$^2$/s$^2$)", S["kmax"], "magma")]):
        for col, (which, lab) in enumerate([(0, "body-fitted"), (1, "immersed")]):
            a = fig.add_subplot(gs[nl + row, 2 * col:2 * col + 2])
            m = a.pcolormesh(DD, ZA, field(which, idx), vmin=0, vmax=vmax, cmap=cmap, shading="auto")
            deco(a, "line %s, %s: %s" % (sec, lab, nm)); a.set_ylabel("height (m)"); fig.colorbar(m, ax=a, pad=0.01)
    a = fig.add_subplot(gs[nl + 2, 0:2])
    m = a.pcolormesh(DD, ZA, field(1, 2) - field(0, 2), vmin=-S["dmax"], vmax=S["dmax"], cmap="RdBu_r", shading="auto")
    deco(a, "speed, immersed minus body-fitted (m/s)"); a.set_xlabel("distance along the line (m)"); a.set_ylabel("height (m)"); fig.colorbar(m, ax=a, pad=0.01)
    a = fig.add_subplot(gs[nl + 2, 2:4])
    for which, c in ((0, CB), (1, CI)):
        f = field(which, 4)
        if np.nanmin(f) < 0: a.contourf(DD, ZA, np.where(f < 0, 1.0, np.nan), levels=[0.5, 1.5], colors=[c], alpha=0.6)
    deco(a, "wind against the inflow direction (blue body-fitted, orange immersed)"); a.set_xlabel("distance along the line (m)")
    # Mast profiles: up to four masts of the section's line, spread along it.
    o = sorted([r for r in obs if r["line"] == sec], key=lambda r: float(r["d"]))
    pick = [o[int(round(t))] for t in np.linspace(0, len(o) - 1, min(4, len(o)))]
    ms = measured(cs) if site == "bolund" and ZIP else []
    if ms:
        ref = [m for m in ms if m["mast"] == REF[cs] and m["h"] == 5][0]; sobs = ref["S"]
        sref = {id(run): float(np.interp(5.0, *run.column(ref["x"], ref["y"])[1:3])) for run in (rb, ri)}
        print("reference %s at 5 m: measured %.2f, body-fitted %.2f, immersed %.2f m/s" % (REF[cs], sobs, sref[id(rb)], sref[id(ri)]))
        print("%-4s %3s | %8s %8s %8s | %8s %8s %8s   (S/Sref, then 100 k/Sref^2)" % ("mast", "h", "obs", "bf", "imm", "obs", "bf", "imm"))
        for m in ms:
            v = []
            for run in (rb, ri):
                _, z, sp, k, _a = run.column(m["x"], m["y"]); v.append((np.interp(m["h"], z, sp) / sref[id(run)], 100 * np.interp(m["h"], z, k) / sref[id(run)] ** 2))
            print("%-4s %3d | %8.3f %8.3f %8.3f | %8.2f %8.2f %8.2f" % (m["mast"], m["h"], m["S"] / sobs, v[0][0], v[1][0], 100 * m["k"] / sobs ** 2, v[0][1], v[1][1]))
    elif "uref" in S:
        sref = {id(rb): S["uref"], id(ri): S["uref"]}
    else:   # Bolund without the measurement archive: each run's own speed at 5 m at the reference mast
        r0 = ([r for r in obs if r["mast"].startswith(REF[cs])] or pick)[0]
        sref = {id(run): float(np.interp(5.0, *run.column(float(r0["x"]), float(r0["y"]))[1:3])) for run in (rb, ri)}
    for n, r in enumerate(pick):
        p = (float(r["x"]), float(r["y"])); lab = r["mast"][:2] if site == "bolund" else r["mast"]
        a = fig.add_subplot(gs[nl + 3, n]); b = fig.add_subplot(gs[nl + 4, n])
        for run, c in ((rb, CB), (ri, CI)):
            _, z, sp, k, _a = run.column(*p)
            a.plot(sp, z, color=c, lw=2, marker=".", ms=4); b.plot(100 * k / sref[id(run)] ** 2, z, color=c, lw=2, marker=".", ms=4)
        a.plot([float(r["speed"])], [float(r["agl"])], "ko", ms=7); a.set_ylim(0, S["prof"]); a.set_xlim(0, S["smax"])
        a.set_title("%s (d = %.0f m)" % (lab, float(r["d"])), loc="left", fontsize=10); a.set_xlabel("speed (m/s)"); a.grid(alpha=0.3)
        mm = [m for m in ms if m["mast"] == r["mast"][:2]]
        if mm: b.plot([100 * m["k"] / sobs ** 2 for m in mm], [m["h"] for m in mm], "ko", ms=7)
        b.set_ylim(0, S["prof"]); b.set_xlim(0, S["kn"]); b.set_xlabel("100 k / S$_{ref}^2$"); b.grid(alpha=0.3)
        b.set_title("%s: k%s" % (lab, " and sonic masts" if mm else ""), loc="left", fontsize=10)
        if n == 0: a.set_ylabel("height above ground (m)"); b.set_ylabel("height above ground (m)")
    fig.savefig(out, dpi=110, bbox_inches="tight"); print(out)


if __name__ == "__main__":
    import warnings; warnings.simplefilter("ignore")
    main(*sys.argv[1:7])
