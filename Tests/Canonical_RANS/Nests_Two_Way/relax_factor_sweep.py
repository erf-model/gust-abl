#!/usr/bin/env python3
"""The nest.feedback_factor of nest.feedback = relax on four nests (tracker N-114, item 13).

    relax_factor_sweep.py run --gust GUST [--work DIR] [--cases C ...] [--factors F ...] [--keep]
                              [--mpi2 "mpiexec -n 2"] [--mpi4 "mpiexec -n 4"] [--timeout S]
    relax_factor_sweep.py report [--work DIR] [--cases C ...] [--plot PNG]

The cases (CASES), each against its own uniform fine run (the nest's cells
everywhere) and with the recommended setting (docs/turbulence_model.rst):

  hill       Nests_Two_Way, the nest hill (the joint run of rt_nests_two_way), 2 ranks;
  lee        Nests_Two_Way, the nest crest alone with its downstream face on
             the lee slope (x = 900 m; rt_nests_feedback_flux), 2 ranks;
  gorge      Tests/Real_Terrain/Columbia_Gorge (prepare: the terrain is
             fetched and the decks written), 4 ranks;
  askervein  Tests/Real_Terrain/Askervein_Nest, the two_way layout (a 20 m
             nest of half-width 1000 m on the 40 m outer domain), 4 ranks.

run:    in DIR/<case>/: the fine reference (ref), the outer domain and the
        nest one way (nest.feedback = none, ow), the defect run (jt), and one
        relax run per factor (rx<F>). hill and lee share DIR/hill/ref.
        --keep skips a run whose log already has its verdict. A run that
        stops without the steady verdict (a stall, the iteration cap or
        --timeout) is recorded and the sweep goes on.
report: for every run, the exchanges to the joint verdict and the RMS of u,
        v, w and k against the reference over the nest (at the nest's own
        cell heights) and over the outer cells outside the footprint (the
        reference averaged over 2 x 2 cells); for relax, the last logged RMS
        of u of the outer domain minus the restricted nest. Prints the
        tables, writes DIR/relax_factor_sweep.md and .csv, and with --plot
        a figure of the exchanges and the u errors over the one-way ones
        against the factor.
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.join(ROOT, "extern_erf"))
sys.path.insert(0, os.path.join(ROOT, "tools"))

FACTORS = (0.1, 0.25, 0.5, 0.75, 1.0, 2.0)
FIELDS = ("u", "v", "w", "tke")
# The recommended setting (docs/turbulence_model.rst); the blend height is 2.5 terrain heights (100 m hill).
RECOMMENDED_HILL = ["turb.realizability=0.309", "turb.nut_relax=0.5", "turb.max_geom_lscale=100.",
                    "convection.hybrid_sigma_max=0.1", "convection.hybrid_height=250."]
HILL_DECK = os.path.join(HERE, "inputs_outer")
GORGE = os.path.join(ROOT, "Tests", "Real_Terrain", "Columbia_Gorge")
ASKERVEIN = os.path.join(ROOT, "Tests", "Real_Terrain", "Askervein_Nest")


def case_table():
    """name: dict(dir, mpi (2 or 4), nest (the nest keys), foot (outer cells i0, j0, ni, nj), ref_dir)."""
    sys.path.insert(0, GORGE)
    sys.path.insert(0, ASKERVEIN)
    import columbia_gorge as cg
    lo = [float(s) for s in cg.NEST[1].split("=")[1].split()]
    hi = [float(s) for s in cg.NEST[2].split("=")[1].split()]
    gfoot = tuple(int(round((lo[d] - cg.OUTER["lo"][d]) / cg.OUTER["dx"])) for d in range(2)) \
        + tuple(int(round((hi[d] - lo[d]) / cg.OUTER["dx"])) for d in range(2))
    a_half, a_dx = 1000.0, 40.0
    a_i0, a_ni = round((3000.0 - a_half) / a_dx), round(2 * a_half / a_dx)
    return {
        "hill": dict(mpi=2, nest=["nest.names=hill", "nest.hill.lo=400 300", "nest.hill.hi=1200 900"],
                     foot=(8, 6, 16, 12), ref="hill", name="hill"),
        "lee": dict(mpi=2, nest=["nest.names=crest", "nest.crest.lo=300 300", "nest.crest.hi=900 900"],
                    foot=(6, 6, 12, 12), ref="hill", name="crest"),
        "gorge": dict(mpi=4, nest=list(cg.NEST), foot=gfoot, ref="gorge", name="rg"),
        "askervein": dict(mpi=4, nest=["nest.names=ht", f"nest.ht.lo={-a_half:g} {-a_half:g}",
                                       f"nest.ht.hi={a_half:g} {a_half:g}"],
                          foot=(a_i0, a_i0, a_ni, a_ni), ref="askervein", name="ht"),
    }


def tag(f):
    return f"rx{f:g}"


VERDICT = re.compile(r"^gust: (?:joint run: )?(steady|not steady|diverged) after (\d+) (exchanges|iterations)", re.M)


def verdict(log):
    """(state, count) of the last verdict line in the log file, or None."""
    if not os.path.exists(log):
        return None
    with open(log, errors="replace") as fh:
        found = VERDICT.findall(fh.read())
    return (found[-1][0], int(found[-1][1])) if found else None


def launch(where, name, cmd, keep, timeout):
    os.makedirs(where, exist_ok=True)
    log = os.path.join(where, name + ".log")
    if keep and verdict(log) is not None and glob.glob(os.path.join(where, name + "_plt[0-9]*")):
        print(f"relax_factor_sweep: {where}/{name}: kept ({verdict(log)})", flush=True)
        return
    print("+ (" + where + ") " + " ".join(cmd), flush=True)
    with open(log, "w") as fh:
        try:
            rc = subprocess.run(cmd, cwd=where, stdout=fh, stderr=subprocess.STDOUT, timeout=timeout).returncode
        except subprocess.TimeoutExpired:
            fh.write(f"\nrelax_factor_sweep: stopped after the {timeout} s timeout\n")
            rc = "timeout"
    print(f"relax_factor_sweep: {where}/{name}: exit {rc}, verdict {verdict(log)}", flush=True)


def run(a):
    cases = case_table()
    gust = os.path.abspath(a.gust)
    for c in a.cases:
        spec = cases[c]
        where = os.path.join(a.work, c)
        mpi = (a.mpi2 if spec["mpi"] == 2 else a.mpi4).split() + [gust]
        if c in ("hill", "lee"):
            ref_dir = os.path.join(a.work, "hill")
            os.makedirs(ref_dir, exist_ok=True)
            with open(HILL_DECK) as fh:
                text = fh.read()
            text, n = re.subn(r"^mesh\.n_cell\s*=.*$", "mesh.n_cell = 64 48 16", text, flags=re.M)
            if n != 1:
                raise SystemExit("relax_factor_sweep: inputs_outer has no single mesh.n_cell line")
            text = re.sub(r"^io\.(plot|chk)_file\s*=.*$", "", text, flags=re.M)
            with open(os.path.join(ref_dir, "inputs_ref"), "w") as fh:
                fh.write(text)
            with open(HILL_DECK) as fh:   # the deck's output names would be given twice
                text = re.sub(r"^io\.(plot|chk)_file\s*=.*$", "", fh.read(), flags=re.M)
            deck = os.path.join(ref_dir, "inputs_outer")
            with open(deck, "w") as fh:
                fh.write(text)
            out = (lambda t: [f"io.plot_file={t}_plt", f"io.chk_file={t}_chk"])
            launch(ref_dir, "ref", mpi + ["inputs_ref"] + RECOMMENDED_HILL + out("ref"), True, a.timeout)
            base = mpi + [deck] + RECOMMENDED_HILL + spec["nest"]
        elif c == "gorge":
            import columbia_gorge as cg
            if not os.path.exists(os.path.join(where, "gorge_outer", "inputs")):
                cg.prepare(where)
            where = os.path.join(where, "gorge_outer")
            out = (lambda t: [f"io.plot_file={t}_plt", f"io.chk_file={t}_chk"])
            launch(where, "ref", mpi + ["inputs_ref"] + out("ref"), a.keep, a.timeout)
            base = mpi + ["inputs"] + spec["nest"]
        else:
            import askervein_nest as an
            deck = os.path.join(ASKERVEIN, "inputs_askervein_domain")
            common = mpi + [deck, f"terrain.file={an.TERRAIN}", f"bc.inflow_file={an.INFLOW}"] + an.TWO_WAY
            out = (lambda t: [f"io.plot_file={t}_plt", f"io.chk_file={t}_chk"])
            launch(where, "ref", common + ["mesh.n_cell=300 300 18"] + out("ref"), a.keep, a.timeout)
            base = common + ["mesh.n_cell=150 150 18"] + spec["nest"]
        runs = [("ow", ["nest.feedback=none"]), ("jt", ["nest.feedback=defect"])]
        runs += [(tag(f), ["nest.feedback=relax", f"nest.feedback_factor={f:g}"]) for f in a.factors]
        for name, extra in runs:
            launch(where, name, base + extra + out(name), a.keep, a.timeout)
    return 0


def run_dir(work, c):
    return os.path.join(work, c, "gorge_outer") if c == "gorge" else os.path.join(work, c)


def last_plt(where, prefix):
    found = []
    for p in glob.glob(os.path.join(where, prefix + "[0-9]*")):
        m = re.fullmatch(re.escape(prefix) + r"(\d+)", os.path.basename(p))
        if m and os.path.isdir(p):
            found.append((int(m.group(1)), p))
    return max(found)[1] if found else None


def outer_minus_nest(log):
    with open(log, errors="replace") as fh:
        found = re.findall(r"^gust: joint +\d+: outer - nest under the nests: .*?RMS u (\S+)", fh.read(), re.M)
    return float(found[-1].rstrip(",")) if found else float("nan")


def measure(work, cases, factors):
    """Rows: dict(case, run, factor, state, exchanges, nest_<f>, outer_<f>, diff)."""
    import numpy as np
    import erf_plotfile
    import speed_difference as sd
    spec = case_table()
    rd = lambda p: {k: np.array(v, dtype=float) for k, v in erf_plotfile.read_fields(p, list(FIELDS) + ["z"])[1].items()}
    rms = lambda x: float(np.sqrt(np.mean(x * x)))
    rows = []
    for c in cases:
        s = spec[c]
        where = run_dir(work, c)
        ref = last_plt(run_dir(work, s["ref"]), "ref_plt")
        if ref is None:
            print(f"relax_factor_sweep: {c}: no reference plotfile; skipped")
            continue
        R = rd(ref)
        i0, j0, ni, nj = s["foot"]
        sl = (slice(2 * i0, 2 * (i0 + ni)), slice(2 * j0, 2 * (j0 + nj)))
        rc = {f: sd.coarsen(R[f]) for f in FIELDS}
        out = np.ones(rc["u"].shape[:2], bool)
        out[i0:i0 + ni, j0:j0 + nj] = False
        o3 = np.broadcast_to(out[:, :, None], rc["u"].shape)
        for name, fac in [("ow", None), ("jt", None)] + [(tag(f), f) for f in factors]:
            log = os.path.join(where, name + ".log")
            v = verdict(log)
            row = dict(case=c, run=name, factor=fac, state=v[0] if v else "missing", exchanges=v[1] if v else None)
            po, pn = last_plt(where, name + "_plt"), last_plt(where, f"{name}_plt_{s['name']}")
            if v and v[0] == "steady" and po and pn:
                O, N = rd(po), rd(pn)
                for f in FIELDS:
                    row["nest_" + f] = rms(N[f] - sd.at_height(N["z"], R["z"][sl], R[f][sl]))
                    row["outer_" + f] = rms((O[f] - rc[f])[o3])
            if name.startswith("rx") and os.path.exists(log):
                row["diff"] = outer_minus_nest(log)
            if name == "ow":   # the one-way run's own verdict is per domain: count nothing
                row["exchanges"] = None
            rows.append(row)
    return rows


def fmt(x, spec="{:.3e}"):
    return "" if x is None or (isinstance(x, float) and x != x) else spec.format(x)


def report(a):
    rows = measure(a.work, a.cases, a.factors)
    keys = ["case", "run", "factor", "state", "exchanges", "diff"] + [f"{r}_{f}" for r in ("nest", "outer")
                                                                       for f in FIELDS]
    with open(os.path.join(a.work, "relax_factor_sweep.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, keys, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    label = {"ow": "one way", "jt": "defect"}
    lines = []
    for c in a.cases:
        cr = [r for r in rows if r["case"] == c]
        if not cr:
            continue
        ow = next((r for r in cr if r["run"] == "ow"), {})
        lines += [f"### {c}", "",
                  "| run | exchanges | outer - nest u | nest u | nest v | nest w | nest k | outer u | outer v | "
                  "outer w | outer k | nest u / one way | outer u / one way |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in cr:
            name = label.get(r["run"], f"relax, f = {r['factor']:g}" if r["factor"] is not None else r["run"])
            ex = str(r["exchanges"]) if r["state"] == "steady" and r["exchanges"] is not None else (
                "" if r["run"] == "ow" and r["state"] == "steady" else f"{r['state']} {fmt(r['exchanges'], '{}')}")
            ratio = lambda k: fmt(r.get(k) / ow[k], "{:.2f}") if r.get(k) is not None and ow.get(k) else ""
            lines.append(f"| {name} | {ex.strip()} | {fmt(r.get('diff'))} | "
                         + " | ".join(fmt(r.get(f"{g}_{f}")) for g in ("nest", "outer") for f in FIELDS)
                         + f" | {ratio('nest_u')} | {ratio('outer_u')} |")
        lines.append("")
    text = "\n".join(lines)
    print(text)
    with open(os.path.join(a.work, "relax_factor_sweep.md"), "w") as fh:
        fh.write(text + "\n")
    if a.plot:
        plot(rows, a.cases, a.plot)
    return 0


def plot(rows, cases, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    colours = {"hill": "#2a6fdb", "lee": "#d9480f", "gorge": "#2b8a3e", "askervein": "#7048e8"}
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.0))
    for c in cases:
        cr = [r for r in rows if r["case"] == c]
        ow = next((r for r in cr if r["run"] == "ow"), {})
        jt = next((r for r in cr if r["run"] == "jt"), {})
        rx = sorted((r for r in cr if r["factor"] is not None and r["state"] == "steady"), key=lambda r: r["factor"])
        if not rx:
            continue
        f = [r["factor"] for r in rx]
        col = colours[c]
        ax[0].plot(f, [r["exchanges"] for r in rx], "o-", color=col, label=c)
        if jt.get("exchanges"):
            ax[0].axhline(jt["exchanges"], color=col, ls=":", lw=1)
        for k, axis in (("nest_u", ax[1]), ("outer_u", ax[2])):
            if ow.get(k):
                axis.plot(f, [r[k] / ow[k] for r in rx], "o-", color=col, label=c)
                if jt.get(k):
                    axis.axhline(jt[k] / ow[k], color=col, ls=":", lw=1)
    titles = ("Exchanges to steady (dotted: defect)", "Nest u RMS / one way (dotted: defect)",
              "Outer u RMS outside / one way (dotted: defect)")
    for axis, t in zip(ax, titles):
        axis.set_xscale("log")
        axis.set_xticks(FACTORS)
        axis.set_xticklabels([f"{f:g}" for f in FACTORS])
        axis.set_xlabel("nest.feedback_factor f")
        axis.set_title(t, fontsize=10)
        axis.grid(True, color="#dddddd", lw=0.6)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
    for axis in ax[1:]:   # the lee defect run is 3 times one way: a log scale keeps the others readable
        axis.set_yscale("log")
    ax[1].axhline(1.0, color="#888888", lw=0.8)
    ax[2].axhline(1.0, color="#888888", lw=0.8)
    ax[0].legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    print(f"relax_factor_sweep: wrote {path}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=["run", "report"])
    ap.add_argument("--work", default="relax_factor_sweep_run", help="working directory")
    ap.add_argument("--cases", nargs="+", default=["hill", "lee", "gorge", "askervein"],
                    choices=["hill", "lee", "gorge", "askervein"])
    ap.add_argument("--factors", type=float, nargs="+", default=list(FACTORS))
    ap.add_argument("--gust", help="run: the gust executable")
    ap.add_argument("--mpi2", default="mpiexec -n 2", help="run: launcher of the hill and lee cases")
    ap.add_argument("--mpi4", default="mpiexec -n 4", help="run: launcher of the gorge and askervein cases")
    ap.add_argument("--timeout", type=float, default=6 * 3600.0, help="run: seconds per run (default 6 h)")
    ap.add_argument("--keep", action="store_true", help="run: skip runs whose log has a verdict")
    ap.add_argument("--plot", help="report: write the figure to this PNG")
    a = ap.parse_args()
    a.work = os.path.abspath(a.work)
    if a.step == "run":
        if not a.gust:
            ap.error("run needs --gust")
        return run(a)
    return report(a)


if __name__ == "__main__":
    sys.exit(main())
