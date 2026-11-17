#!/usr/bin/env python3
"""Two-way nesting over the Columbia River Gorge: prepare, run, compare.

    columbia_gorge.py prepare [--work DIR]
    columbia_gorge.py run --gust GUST [--work DIR] [--mpi "mpiexec -n 4"] [--relax F ...] [--relax-only]
    columbia_gorge.py compare [--work DIR]

prepare: fetches the terrain (tools/fetch_terrain.py, SRTM 1 arc-second
         from the AWS Terrain Tiles; nothing is downloaded when the tiles
         are in the cache) to DIR/dem.asc, and writes the outer deck
         DIR/gorge_outer/inputs from case.yaml (tools/make_gust_case.py),
         with turb.nut_relax_end = 0.01 added, plus inputs_ref, the same
         deck on 100 m cells.
run:     runs, one after another in DIR/gorge_outer: the uniform 100 m
         reference (ref_*), the outer domain and the nest one way
         (nest.feedback = none, ow_*), and the joint run (nest.feedback =
         defect, jt_*); with --relax F ..., then one joint run per factor
         with nest.feedback = relax, nest.feedback_factor = F (rx<F>_*);
         --relax-only skips the first three (already run).
compare: the RMS of u, v, w and k of each run against the reference: over
         the nest, over the outer cells outside the nest's footprint, over
         the footprint, and in the outer cells east of the nest (downstream);
         then the cells in the lowest 200 m whose horizontal speed differs
         from the reference by more than 0.25, 0.5 and 1 m/s
         (tools/speed_difference.py), and with --plot PNG a terrain map of
         each column's largest difference. The nest is compared at its own
         cell heights (a joint nest's levels differ from the reference's).
         The relax runs found in DIR (rx<F>_*) are compared as well.

The nest (NEST) covers the river bend and the south face of the Columbia
Hills: its north face cuts the ridge and its east face crosses the river.
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "extern_erf"))
sys.path.insert(0, os.path.join(ROOT, "tools"))

BOX = ["45.40", "45.92", "-121.25", "-120.35"]   # south north west east, degrees
NEST = ["nest.names=rg", "nest.rg.lo=-4000 -2400", "nest.rg.hi=5600 4800"]
OUTER = dict(lo=(-9600.0, -8000.0), dx=200.0)    # outer domain corner and cell size, m
FIELDS = ["u", "v", "w", "tke"]


def sh(cmd, cwd, log=None):
    print("+", " ".join(cmd), flush=True)
    if log is None:
        subprocess.run(cmd, cwd=cwd, check=True)
        return
    with open(os.path.join(cwd, log), "w") as fh:
        rc = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT).returncode
    if rc != 0:
        raise SystemExit(f"columbia_gorge: {' '.join(cmd)} exited {rc}; see {os.path.join(cwd, log)}")


def prepare(work):
    os.makedirs(work, exist_ok=True)
    sh([sys.executable, os.path.join(ROOT, "tools", "fetch_terrain.py"), "--box", *BOX,
        "--out", os.path.join(work, "dem.asc")], work)
    shutil.copy(os.path.join(HERE, "case.yaml"), os.path.join(work, "case.yaml"))
    sh([sys.executable, os.path.join(ROOT, "tools", "make_gust_case.py"), "case.yaml", "--overwrite"], work)
    case = os.path.join(work, "gorge_outer")
    deck = os.path.join(case, "inputs")
    with open(deck) as fh:
        text = fh.read()
    # Nests need the outer run's checkpoint; the k plotfiles are not needed here.
    text = re.sub(r"^io\.turb_plot_file\s*=.*$", "", text, flags=re.M)
    if not re.search(r"^io\.write_checkpoint\s*=", text, re.M):
        text += "io.write_checkpoint = 1\n"
    # Without the end of the nu_t relaxation the outer run stalls in a two-iteration cycle (README, Purpose).
    if not re.search(r"^turb\.nut_relax_end\s*=", text, re.M):
        text += "turb.nut_relax_end = 0.01\n"
    with open(deck, "w") as fh:
        fh.write(text)
    ref, n = re.subn(r"^mesh\.n_cell\s*=.*$", "mesh.n_cell = 192 160 37", text, flags=re.M)
    if n != 1:
        raise SystemExit("columbia_gorge: the deck has no single mesh.n_cell line")
    with open(os.path.join(case, "inputs_ref"), "w") as fh:
        fh.write(ref)
    print(f"columbia_gorge: wrote {deck} and inputs_ref")


def run(work, gust, mpi, relax=(), relax_only=False):
    case = os.path.join(work, "gorge_outer")
    launch = mpi.split() + [os.path.abspath(gust)]
    if not relax_only:
        sh(launch + ["inputs_ref", "io.plot_file=ref_plt", "io.chk_file=ref_chk"], case, "ref.log")
        sh(launch + ["inputs", *NEST, "nest.feedback=none", "io.plot_file=ow_plt", "io.chk_file=ow_chk"], case, "ow.log")
        sh(launch + ["inputs", *NEST, "nest.feedback=defect", "io.plot_file=jt_plt", "io.chk_file=jt_chk"], case, "jt.log")
    for f in relax:
        tag = f"rx{f:g}"
        sh(launch + ["inputs", *NEST, "nest.feedback=relax", f"nest.feedback_factor={f:g}", f"io.plot_file={tag}_plt",
                     f"io.chk_file={tag}_chk"], case, tag + ".log")


def relax_runs(case):
    """The tags of the relax runs in case (rx<F>), by factor."""
    tags = {os.path.basename(p).split("_plt")[0] for p in glob.glob(os.path.join(case, "rx*_plt_rg[0-9]*"))}
    return sorted(tags, key=lambda t: float(t[2:]))


def last(case, prefix):
    found = []
    for p in glob.glob(os.path.join(case, prefix + "[0-9]*")):
        m = re.fullmatch(re.escape(prefix) + r"(\d+)", os.path.basename(p))
        if m and os.path.isdir(p):
            found.append((int(m.group(1)), p))
    if not found:
        raise SystemExit(f"columbia_gorge: no plotfile {prefix}* in {case}")
    return max(found)[1]


def compare(work, plot=None):
    import erf_plotfile
    import speed_difference as sd
    case = os.path.join(work, "gorge_outer")
    rd = lambda p: {k: np.array(v, dtype=float) for k, v in erf_plotfile.read_fields(p, FIELDS + ["z"])[1].items()}
    R = rd(last(case, "ref_plt"))
    OO, OJ = rd(last(case, "ow_plt")), rd(last(case, "jt_plt"))
    NO, NJ = rd(last(case, "ow_plt_rg")), rd(last(case, "jt_plt_rg"))
    RX = {t: (rd(last(case, t + "_plt")), rd(last(case, t + "_plt_rg"))) for t in relax_runs(case)}
    lo = [float(s) for s in NEST[1].split("=")[1].split()]
    hi = [float(s) for s in NEST[2].split("=")[1].split()]
    i0, j0 = (int(round((lo[d] - OUTER["lo"][d]) / OUTER["dx"])) for d in range(2))
    ni, nj = (int(round((hi[d] - lo[d]) / OUTER["dx"])) for d in range(2))
    rms = lambda a: float(np.sqrt(np.mean(a * a)))
    sl = (slice(2 * i0, 2 * (i0 + ni)), slice(2 * j0, 2 * (j0 + nj)))
    zr = R["z"][sl]
    out = np.ones((R["z"].shape[0] // 2, R["z"].shape[1] // 2), bool)
    out[i0:i0 + ni, j0:j0 + nj] = False
    fp = (slice(i0, i0 + ni), slice(j0, j0 + nj))
    dn = (slice(i0 + ni, None), slice(None))
    print(f"footprint: outer cells i {i0}..{i0 + ni - 1}, j {j0}..{j0 + nj - 1}")
    print(f"{'field':5s} {'nest ow':>10s} {'nest jt':>10s} {'out ow':>10s} {'out jt':>10s} "
          f"{'foot ow':>10s} {'foot jt':>10s} {'down ow':>10s} {'down jt':>10s}")
    for f in FIELDS:
        rc = sd.coarsen(R[f])
        # A joint nest's levels follow the parent's aloft: compare it at its own cell heights, never by index.
        ro, rj = sd.at_height(NO["z"], zr, R[f][sl]), sd.at_height(NJ["z"], zr, R[f][sl])
        o3 = np.broadcast_to(out[:, :, None], rc.shape)
        vals = [rms(NO[f] - ro), rms(NJ[f] - rj), rms((OO[f] - rc)[o3]), rms((OJ[f] - rc)[o3]),
                rms((OO[f] - rc)[fp]), rms((OJ[f] - rc)[fp]), rms((OO[f] - rc)[dn]), rms((OJ[f] - rc)[dn])]
        print(f"{f:5s} " + " ".join(f"{v:10.4e}" for v in vals))
    for t, (O, N) in RX.items():
        print(f"\n{t} (nest.feedback = relax, nest.feedback_factor = {t[2:]}):")
        print(f"{'field':5s} {'nest':>10s} {'out':>10s} {'foot':>10s} {'down':>10s}")
        for f in FIELDS:
            rc = sd.coarsen(R[f])
            rn = sd.at_height(N["z"], zr, R[f][sl])
            o3 = np.broadcast_to(out[:, :, None], rc.shape)
            vals = [rms(N[f] - rn), rms((O[f] - rc)[o3]), rms((O[f] - rc)[fp]), rms((O[f] - rc)[dn])]
            print(f"{f:5s} " + " ".join(f"{v:10.4e}" for v in vals))

    # Wind-speed differences in m/s, in the lowest sd.AGL m: counts above the thresholds and a map.
    Rs = sd.speed(R)
    Rc, Zc = sd.coarsen(Rs), sd.coarsen(R["z"])
    mc = sd.near_ground(Zc)
    print(f"\nhorizontal speed |S - S_ref| below {sd.AGL:g} m above the ground "
          f"(reference mean {Rs[sd.near_ground(R['z'])].mean():.2f} m/s)")
    rows, layers = [], {}
    for name, N, O in [("one-way", NO, OO), ("joint", NJ, OJ)] + [(t, N, O) for t, (O, N) in RX.items()]:
        d_n = sd.speed(N) - sd.at_height(N["z"], zr, Rs[sl])
        m_n = sd.near_ground(N["z"])
        d_o = sd.speed(O) - Rc
        rows += [(name, "nest", d_n[m_n]), (name, "foot", d_o[fp][mc[fp]]),
                 (name, "outer", d_o[mc & out[:, :, None]])]
        layers[name] = (sd.column_max(d_n, m_n), sd.column_max(d_o, mc))
    sd.print_table(rows)
    for name in ["ref", "ow", "jt"] + list(RX):
        log = os.path.join(case, name + ".log")
        if os.path.exists(log):
            with open(log) as fh:
                for line in fh:
                    if re.match(r"gust: (stage account|joint run: steady)", line):
                        print(f"{name}: {line.strip()}")
    if plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        dxf = OUTER["dx"] / 2
        h = R["z"][:, :, 0]
        xf = OUTER["lo"][0] + dxf * (np.arange(h.shape[0]) + 0.5)
        yf = OUTER["lo"][1] + dxf * (np.arange(h.shape[1]) + 0.5)
        xc = OUTER["lo"][0] + OUTER["dx"] * (np.arange(out.shape[0]) + 0.5)
        yc = OUTER["lo"][1] + OUTER["dx"] * (np.arange(out.shape[1]) + 0.5)
        Xc, Yc = np.meshgrid(xc, yc, indexing="ij")
        Xn, Yn = np.meshgrid(xf[sl[0]], yf[sl[1]], indexing="ij")
        fig, ax = plt.subplots(1, 2, figsize=(12, 5.2))
        for a, name in zip(ax, ("one-way", "joint")):
            cn, co = layers[name]
            co = co.copy()
            co[fp] = np.nan   # the nest covers its footprint on the map
            sd.plot_map(a, (xf, yf, h), [(Xc, Yc, co), (Xn, Yn, cn)],
                        (lo[0], lo[1], hi[0] - lo[0], hi[1] - lo[1]), f"Columbia Gorge, {name}")
        sd.legend(fig)
        fig.tight_layout(rect=(0, 0.06, 1, 1))
        fig.savefig(plot, dpi=110)
        print(f"columbia_gorge: wrote {plot}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["prepare", "run", "compare"])
    ap.add_argument("--work", default="columbia_gorge_run", help="working directory (default columbia_gorge_run)")
    ap.add_argument("--gust", help="the gust executable (run)")
    ap.add_argument("--plot", help="compare: also write the map of the speed differences to this PNG")
    ap.add_argument("--mpi", default="mpiexec -n 4", help='MPI launcher (default "mpiexec -n 4")')
    ap.add_argument("--relax", type=float, nargs="+", default=[], metavar="F",
                    help="run: also the joint run with nest.feedback = relax at each nest.feedback_factor F")
    ap.add_argument("--relax-only", action="store_true", help="run: only the relax runs (the others already ran)")
    a = ap.parse_args()
    work = os.path.abspath(a.work)
    if a.step == "prepare":
        prepare(work)
    elif a.step == "run":
        if not a.gust:
            ap.error("run needs --gust")
        if a.relax_only and not a.relax:
            ap.error("--relax-only needs --relax")
        run(work, a.gust, a.mpi, a.relax, a.relax_only)
    else:
        compare(work, a.plot)
    return 0


if __name__ == "__main__":
    sys.exit(main())
