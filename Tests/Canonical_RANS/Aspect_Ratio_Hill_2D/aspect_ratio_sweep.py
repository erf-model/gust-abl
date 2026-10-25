#!/usr/bin/env python3
"""Aspect-ratio sweep on the 2-D hill: iterations to a relative residual of 1e-4 for dx/dz1 = 10 ... 5000.

    python3 aspect_ratio_sweep.py --gust build-rel/gust --out ar-runs [--np 2] [--ratios 10 25 50 100 200 500]
        [--max-iterations 8000] [--check] [--mpiexec "mpirun -np"]

For each ratio the deck inputs_aspect_ratio_hill is copied with mesh.dz0 = dx/ratio and enough layers (stretched by
mesh.stretch_ratio) to reach the first level at or above 1.5 km (1.5-1.65 km). The run's log gives the first iteration with r < 1e-4 and
the steady verdict. Iteration lines are printed every 10 iterations, so a count is exact only where the run stops at
the target; otherwise it is rounded up to a multiple of 10. With --check the script exits non-zero unless every ratio up to 200 reaches r < 1e-4 within
the iteration limit (the CTest).
"""
import argparse
import math
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
DECK = os.path.join(HERE, "inputs_aspect_ratio_hill")
TOP = 1500.0
CHECKED_UP_TO = 200


def deck_get(text, key):
    m = re.search(r"^" + re.escape(key) + r"\s*=\s*(.*)$", text, re.M)
    return m.group(1).split("#")[0].strip()


def write_deck(dst, ratio, max_it):
    text = open(DECK).read()
    nx, ny, _ = deck_get(text, "mesh.n_cell").split()
    dx = (float(deck_get(text, "mesh.prob_hi").split()[0]) - float(deck_get(text, "mesh.prob_lo").split()[0])) / int(nx)
    r = float(deck_get(text, "mesh.stretch_ratio"))
    dz = dx / ratio
    nz = math.ceil(math.log(1 + (r - 1) * TOP / dz) / math.log(r))
    for key, val in (("mesh.n_cell", f"{nx} {ny} {nz}"), ("mesh.dz0", f"{dz:g}"), ("steady.max_iterations", str(max_it))):
        text = re.sub(r"^" + re.escape(key) + r"\s*=.*$", f"{key} = {val}", text, flags=re.M)
    open(dst, "w").write(text)
    return dz, nz


def first_below(log, target):
    for line in open(log):
        m = re.match(r"gust: it\s+(\d+)\s+cfl\s+\S+\s+r\s+(\S+)", line)
        if m and float(m.group(2)) < target:
            return int(m.group(1))
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gust", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--np", type=int, default=2)
    ap.add_argument("--ratios", type=float, nargs="+", default=[10, 25, 50, 100, 200, 500])
    ap.add_argument("--max-iterations", type=int, default=8000)
    ap.add_argument("--target", type=float, default=1e-4)
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--mpiexec", default="mpirun -np", help="launcher and its rank flag; empty for a serial build")
    ap.add_argument("--gust-arg", action="append", default=[], help="extra key=value for gust (repeatable)")
    a = ap.parse_args()
    ok = True
    print("| dx/dz1 | dz1 (m) | layers | iterations to r < 1e-4 | verdict |")
    print("|---|---|---|---|---|")
    for ratio in a.ratios:
        d = os.path.join(os.path.abspath(a.out), f"ar{ratio:g}")
        os.makedirs(d, exist_ok=True)
        dz, nz = write_deck(os.path.join(d, "inputs"), ratio, a.max_iterations)
        with open(os.path.join(d, "log.txt"), "w") as lg:
            launch = a.mpiexec.split() + [str(a.np)] if a.mpiexec.strip() else []
            rc = subprocess.run(launch + [os.path.abspath(a.gust), "inputs"] + a.gust_arg, cwd=d,
                                stdin=subprocess.DEVNULL, stdout=lg, stderr=subprocess.STDOUT).returncode
        log = os.path.join(d, "log.txt")
        it = first_below(log, a.target)
        verdict = next((l.split("gust: ")[1].strip() for l in open(log) if re.search(r"steady after|not steady", l)), "no verdict")
        print(f"| {ratio:g} | {dz:g} | {nz} | {it if it is not None else 'not reached'} | {verdict} |", flush=True)
        if ratio <= CHECKED_UP_TO and (it is None or rc != 0):
            ok = False
    if a.check and not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
