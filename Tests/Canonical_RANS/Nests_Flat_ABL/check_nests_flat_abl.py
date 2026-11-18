#!/usr/bin/env python3
"""Nests on a flat, horizontally uniform ABL: the nest against the outer domain, one way and joint.

    check_nests_flat_abl.py [--one-way] -- COMMAND [ARGS ...]

COMMAND runs one of the case's decks (inputs_neutral_s, inputs_buoyant_b, inputs_stable_gabls1). The script copies
the deck's bc.inflow_file table (if any) from the deck's folder to the current directory and runs COMMAND there
with nest.feedback = defect: the outer domain, the nest one way (their first plotfiles), then the joint run
(their last plotfiles). With --one-way it runs nest.feedback = none and checks the one-way plotfiles only.

The outer domain is doubly periodic and starts from the steady column (init.first_guess = column), so the
steady 3-D state is the column itself, the same at every (x, y), and a nest inside it must return that column.

Checked, for the one-way plotfiles and (without --one-way) the joint ones:
a. the run exits 0; the log has "steady after" for the outer domain and the nest and, joint, "joint run: steady
   after N exchanges";
b. the outer domain is horizontally uniform: on every level, u, v and k differ from the level mean by at most
   UNIFORM (m/s, m^2/s^2). For a deck without nest.damping = outer (inputs_stable_gabls1: u and v are not damped
   aloft) only on the levels whose mean k is above K_FLOOR times the largest: above them the joint run leaves a
   neutral roll. The decks with nest.damping = outer use UNIFORM_DAMPED (and NEST_UV_DAMPED in d.);
c. the outer level means stay near the column's (the column_guess.dat the run writes): u and v within COLUMN_UV
   (m/s), k within COLUMN_K times the column's largest k, theta within COLUMN_THETA (K). The column and the 3-D
   stage do not discretize the capped, stratified column alike (the 3-D run's first residuals are 0.06 to 0.3
   of the first-iteration scale), so this bound is loose: it catches a 3-D state far from the column, such as
   theta frozen at the wrong profile;
d. the nest's level means are the outer domain's (the nest has the outer cells' heights: ratio 2 2 1): u and v
   within NEST_UV (m/s), k within NEST_K times the largest k, theta within NEST_THETA (K), and the nest is
   horizontally uniform to the bound of b.

The measured values are in README.md. The exit code is the verdict.
"""
import glob
import os
import re
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402

# Measured (README.md), largest over the three decks and both stages:
UNIFORM = 1.0e-4        # m/s, m^2/s^2 (stable, turbulent levels): measured 2.4e-5 (one-way nest), 2.1e-5 (joint)
K_FLOOR = 1.0e-4        # a deck without nest.damping = outer (stable): uniformity is bounded only on levels whose
                        # mean k is above K_FLOOR times the largest; above them k is about 0, u and v are not
                        # damped, and an x-uniform u(y) in geostrophic balance is a neutral mode (measured 4.0e-4
                        # after the joint run; 2.1e-5 below)
COLUMN_UV = 0.3         # m/s: measured 0.149 (neutral, v in the mixed layer), 0.061 (stable), 0.005 (buoyant)
COLUMN_K = 0.1          # relative: measured 0.035 (stable), 0.027 (neutral), 0.014 (buoyant)
COLUMN_THETA = 0.02     # K: measured 0.012 (stable joint, theta solved), 2.0e-3 (stable one way); 2e-12 frozen
NEST_UV = 2.0e-4        # m/s: measured 6.5e-5 (stable, joint); neutral and buoyant <= 1.2e-7
NEST_K = 1.0e-4         # relative: measured 3.5e-6 (stable); neutral and buoyant <= 6e-8
NEST_THETA = 1.0e-5     # K: measured 1.6e-6 (stable); 0 with theta frozen
# The decks with nest.damping = outer (neutral, buoyant) are uniform and nested far more tightly; their bounds:
UNIFORM_DAMPED = 1.0e-6  # measured 2.0e-7 (neutral), 3.7e-8 (buoyant)
NEST_UV_DAMPED = 1.0e-6  # m/s: measured 1.1e-7 (neutral, one way), 2e-8 (buoyant)
FIELDS = ["u", "v", "theta", "tke"]


FAILURES = []


def fail(msg):
    raise SystemExit("check_nests_flat_abl: " + msg)


def note(msg):
    """A failed bound: printed now, the verdict at the end (every measure is printed first)."""
    print("check_nests_flat_abl: FAIL: " + msg)
    FAILURES.append(msg)


def deck_of(cmd):
    for a in cmd:
        if os.path.basename(a).startswith("inputs_") and os.path.isfile(a):
            return a
    fail("COMMAND has no inputs_* deck")


def deck_value(deck, key):
    with open(deck) as fh:
        for line in fh:
            m = re.match(r"\s*" + re.escape(key) + r"\s*=\s*(\S+)", line.split("#", 1)[0])
            if m:
                return m.group(1)
    return None


def plotfiles(prefix):
    found = []
    for p in glob.glob(prefix + "[0-9]*"):
        m = re.fullmatch(re.escape(prefix) + r"(\d+)", p)
        if m and os.path.isdir(p):
            found.append((int(m.group(1)), p))
    return [p for _, p in sorted(found)]


def fields(plt):
    _, f = erf_plotfile.read_fields(plt, FIELDS)
    return {k: np.array(f[k], dtype=float) for k in FIELDS}


def means(f):
    return {k: v.mean(axis=(0, 1)) for k, v in f.items()}


def uniform(f, what, bound, turbulent_only):
    lev = np.max([np.abs(f[k] - f[k].mean(axis=(0, 1), keepdims=True)).max(axis=(0, 1)) for k in ("u", "v", "tke")],
                 axis=0)
    if turbulent_only:
        km = f["tke"].mean(axis=(0, 1))
        sel = km > K_FLOOR * km.max()
        dev = float(lev[sel].max())
        above = float(lev[~sel].max()) if (~sel).any() else 0.0
        print(f"check_nests_flat_abl: {what}: largest departure from the level mean (u, v, k) {dev:.3e} on the "
              f"{int(sel.sum())} levels with k > {K_FLOOR:g} k_max; {above:.3e} above them (not bounded)")
    else:
        dev = float(lev.max())
        print(f"check_nests_flat_abl: {what}: largest departure from the level mean (u, v, k) {dev:.3e}")
    if not dev <= bound:
        note(f"{what} is not horizontally uniform: {dev:.3e} > {bound:g}")


def column():
    rows = np.loadtxt("column_guess.dat", comments="#")
    return {"u": rows[1:, 1], "v": rows[1:, 2], "theta": rows[1:, 3], "tke": rows[1:, 4]}


def compare(a, b, uv, k_rel, th, what):
    du = max(float(np.abs(a[k] - b[k]).max()) for k in ("u", "v"))
    dk = float(np.abs(a["tke"] - b["tke"]).max()) / max(float(np.abs(b["tke"]).max()), 1e-30)
    dt = float(np.abs(a["theta"] - b["theta"]).max())
    print(f"check_nests_flat_abl: {what}: max |du|, |dv| {du:.3e} m/s, k {dk:.3e} (relative), theta {dt:.3e} K")
    if not (du <= uv and dk <= k_rel and dt <= th):
        note(f"{what}: beyond u, v {uv:g} m/s, k {k_rel:g}, theta {th:g} K")


def check_stage(outer_plt, nest_plt, col, stage, bound, nest_uv, turbulent_only):
    o, n = fields(outer_plt), fields(nest_plt)
    uniform(o, f"{stage} outer domain ({outer_plt})", bound, turbulent_only)
    uniform(n, f"{stage} nest ({nest_plt})", bound, turbulent_only)
    mo = means(o)
    compare(mo, col, COLUMN_UV, COLUMN_K, COLUMN_THETA, f"{stage} outer domain against the column")
    compare(means(n), mo, nest_uv, NEST_K, NEST_THETA, f"{stage} nest against the outer domain")


def main():
    argv = sys.argv[1:]
    if "--" not in argv:
        fail("usage: check_nests_flat_abl.py [--one-way] -- COMMAND [ARGS ...]")
    opts, cmd = argv[:argv.index("--")], argv[argv.index("--") + 1:]
    one_way = "--one-way" in opts
    deck = deck_of(cmd)
    table = deck_value(deck, "bc.inflow_file")
    damped = deck_value(deck, "nest.damping") == "outer"
    if table:
        shutil.copy(os.path.join(os.path.dirname(os.path.abspath(deck)), table), table)
    for p in glob.glob("outer_plt*") + glob.glob("outer_chk*") + ["column_guess.dat"]:
        if os.path.isdir(p):
            shutil.rmtree(p)
        elif os.path.exists(p):
            os.remove(p)
    proc = subprocess.run(cmd + ["nest.feedback=" + ("none" if one_way else "defect")], stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, text=True, timeout=3600)
    log = proc.stdout
    print(log[-3000:])
    # a.
    if proc.returncode != 0:
        fail(f"the run exited {proc.returncode}")
    if len(re.findall(r"^gust: steady after \d+ iterations", log, re.M)) < 2:
        fail("the log does not have 'steady after' for the outer domain and the nest")
    if not one_way and not re.search(r"^gust: joint run: steady after \d+ exchanges", log, re.M):
        fail("the joint run is not steady")
    outer, nest = plotfiles("outer_plt"), plotfiles("outer_plt_mid")
    want = 1 if one_way else 2
    if len(outer) != want or len(nest) != want:
        fail(f"expected {want} plotfile(s) each of the outer domain and the nest, found {outer} and {nest}")
    col = column()
    # b.-d.
    uni, nuv = (UNIFORM_DAMPED, NEST_UV_DAMPED) if damped else (UNIFORM, NEST_UV)
    check_stage(outer[0], nest[0], col, "one way", uni, nuv, not damped)
    if not one_way:
        check_stage(outer[1], nest[1], col, "joint", uni, nuv, not damped)
    if FAILURES:
        fail(f"{len(FAILURES)} bound(s) failed (above)")
    print("check_nests_flat_abl: passed")


if __name__ == "__main__":
    main()
