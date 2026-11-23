#!/usr/bin/env python3
"""Time-averaged unsteady RANS against the steady answer.

    urans_averages.py deck STEADY_DECK OUT
    urans_averages.py cost LOG [LOG ...]
    urans_averages.py fields STEADY_PLT UNSTEADY_PLT [--every N] [--max_rel X]

deck    writes the unsteady deck of a steady deck: steady.stage = unsteady,
        steady.sequencing = 0, mesh.hierarchy_levels = 1, and without the keys
        that have no effect in such a run and abort it (each one is printed):
        steady.max_iterations, residual_drop_orders, monitor_rel_tol,
        monitor_window, confidence, confidence_height; steady.sequencing_*;
        terrain.hierarchy_passes; io.iteration_log; and, when the deck has no
        blending, convection.hybrid_height, hybrid_relax, hybrid_freeze.
        The time step, the stop time and the averaging start go on the command
        line of the run (unsteady.dt, unsteady.stop_time, unsteady.average_start).
        A deck with convection.hybrid_sigma_max > 0 is refused: the blending is
        not available in an unsteady run.
cost    prints, for each log, the iterations of a steady run and its verdict
        (with the iterations on the coarser levels of a sequenced run, which
        the verdict line does not count), or the physical steps and the inner
        iterations of an unsteady run (the steps of this log: a restarted run
        is counted from its restart), with the steps that stopped at
        unsteady.max_inner or diverged.
fields  compares the running means u_mean, v_mean of an unsteady plotfile with
        u, v of a steady plotfile on the same mesh (the same index range,
        origin, spacing and cell heights), in every N-th column (default 1; the
        last column of each direction is always included): the largest and the
        RMS difference of the horizontal velocity, as a share of the largest
        steady speed. Also printed: the largest difference of the last step's u,
        v from the steady field, and the largest standard deviation of u over
        the averaging window (the square root of u_var). A NaN fails. With --max_rel the exit code
        is 1 when the largest difference is above it.
Pure Python.
"""
import argparse
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
import askervein_lines as al  # noqa: E402

KEY = re.compile(r"^\s*([A-Za-z_][\w.]*)\s*=\s*(.*?)\s*(#.*)?$")


# Keys that abort a run with steady.stage = unsteady, steady.sequencing = 0 and
# mesh.hierarchy_levels = 1 ("has no effect"): App.cpp, Sequencing.cpp,
# MeshHierarchy.cpp, Convection.cpp.
DROP = {"steady.max_iterations", "steady.residual_drop_orders", "steady.monitor_rel_tol", "steady.monitor_window",
        "steady.confidence", "steady.confidence_height", "terrain.hierarchy_passes", "io.iteration_log",
        "convection.hybrid_height", "convection.hybrid_relax", "convection.hybrid_freeze"}


def cmd_deck(a):
    out, seen = [], set()
    for line in open(a.steady_deck):
        m = KEY.match(line)
        if not m:
            out.append(line)
            continue
        key, val = m.group(1), m.group(2)
        seen.add(key)
        if key == "convection.hybrid_sigma_max":
            try:
                sigma = float(val)
            except ValueError:
                sys.exit(f"{a.steady_deck}: convection.hybrid_sigma_max = '{val}' is not a number")
            if sigma > 0:
                sys.exit(f"{a.steady_deck}: convection.hybrid_sigma_max = {val}; the blending is not available with "
                         "steady.stage = unsteady. Use a deck without it.")
        if key in DROP or key.startswith("steady.sequencing_"):
            print(f"dropped {key} = {val}")
            continue
        if key == "steady.stage":
            line = "steady.stage = unsteady\n"
        elif key == "steady.sequencing":
            line = "steady.sequencing = 0\n"
        elif key == "mesh.hierarchy_levels":
            line = "mesh.hierarchy_levels = 1\n"
        out.append(line)
    if "steady.stage" not in seen:
        out.append("steady.stage = unsteady\n")
    with open(a.out, "w") as f:
        f.writelines(out)
    print(f"wrote {a.out}")


STEADY = re.compile(r"gust: (not steady|steady|diverged) after (\d+) iterations")
COARSE = re.compile(r"gust: level (\d+): loose stop after (\d+) iterations")
STEP = re.compile(r"gust: step (\d+): t = (\S+) s, dt = (\S+) s, (\d+) inner iterations")


def cmd_cost(a):
    for log in a.logs:
        steps, first, last, inner, t, verdict, coarse, at_max, div = 0, None, None, 0, 0.0, None, [], 0, 0
        for line in open(log, errors="replace"):
            m = STEP.match(line)
            if m:
                last, t = int(m.group(1)), float(m.group(2))
                first = last if first is None else first
                steps += 1
                inner += int(m.group(4))
                at_max += "unsteady.max_inner reached" in line
                div += "(diverged)" in line
                continue
            m = COARSE.match(line)
            if m:
                coarse.append((int(m.group(1)), int(m.group(2))))
                continue
            m = STEADY.match(line)
            if m:
                verdict = (m.group(1), int(m.group(2)))
        if steps:
            print(f"{log}: unsteady, {steps} steps ({first} to {last}) to t = {t:g} s, {inner} inner iterations "
                  f"({inner / steps:.2f} per step); {at_max} steps at unsteady.max_inner, {div} diverged")
        elif verdict:
            extra = "".join(f"; level {lev}: {n}" for lev, n in coarse)
            print(f"{log}: {verdict[0]} after {verdict[1]} iterations" + (f" (coarser levels{extra})" if coarse else ""))
        else:
            sys.exit(f"{log}: neither a steady verdict nor an unsteady step")


def cmd_fields(a):
    st, un = al.Columns(a.steady_plt), al.Columns(a.unsteady_plt, mean=True)
    h = st.hdr
    if (h["lo"], h["hi"]) != (un.hdr["lo"], un.hdr["hi"]):
        sys.exit(f"the plotfiles have different index ranges: {h['lo']}..{h['hi']} and {un.hdr['lo']}..{un.hdr['hi']}")
    for q in ("prob_lo", "dx"):
        if any(abs(x - y) > 1e-9 * max(abs(x), abs(y), 1.0) for x, y in zip(h[q], un.hdr[q])):
            sys.exit(f"the plotfiles have different {q}: {h[q]} and {un.hdr[q]}")
    nx, ny = h["hi"][0] - h["lo"][0] + 1, h["hi"][1] - h["lo"][1] + 1
    worst, s2, n, scale, where, zero = 0.0, 0.0, 0, 0.0, None, True
    cols_j = sorted(set(range(0, ny, a.every)) | {ny - 1})
    cols_i = sorted(set(range(0, nx, a.every)) | {nx - 1})
    for j in cols_j:
        for i in cols_i:
            us, vs, zs = st.column(i, j)
            um, vm, zu = un.column(i, j)
            if len(zu) != len(zs) or max(abs(p - q) for p, q in zip(zs, zu)) > 1e-6 * (zs[-1] - zs[0]):
                sys.exit(f"the plotfiles have different cell heights in column ({i}, {j}): not the same mesh")
            for k in range(len(us)):
                d = math.hypot(um[k] - us[k], vm[k] - vs[k])
                if not math.isfinite(d):
                    sys.exit(f"FAIL: a field is not finite at cell ({i}, {j}, {k})")
                if d > worst:
                    worst, where = d, (i, j, k)
                s2 += d * d
                n += 1
                scale = max(scale, math.hypot(us[k], vs[k]))
                zero = zero and um[k] == 0 and vm[k] == 0
    if n == 0 or not scale > 0:
        sys.exit("FAIL: no cells compared, or a steady field at rest")
    if zero:
        sys.exit(f"FAIL: u_mean and v_mean of {a.unsteady_plt} are zero: it was written before unsteady.average_start")
    # The last step against the steady field, and the standard deviation over the window.
    last, var = al.Columns(a.unsteady_plt), al.Columns(a.unsteady_plt, fields=("u_var", "v_var"))
    dl, wl, sd, ws = 0.0, None, 0.0, None
    for j in cols_j:
        for i in cols_i:
            us, vs, _ = st.column(i, j)
            ul, vl, _ = last.column(i, j)
            uv, _, _ = var.column(i, j)
            for k in range(len(us)):
                d = math.hypot(ul[k] - us[k], vl[k] - vs[k])
                if not (math.isfinite(d) and uv[k] >= 0):
                    sys.exit(f"FAIL: the last step or u_var is not finite, or u_var < 0, at cell ({i}, {j}, {k})")
                if d > dl:
                    dl, wl = d, (i, j, k)
                if math.sqrt(uv[k]) > sd:
                    sd, ws = math.sqrt(uv[k]), (i, j, k)
    rel = worst / scale
    print(f"time mean against steady, {n} cells (every {a.every} columns): largest difference {worst:.4e} m/s "
          f"({100 * rel:.4f} % of the largest steady speed {scale:.3f} m/s) at cell {where}; "
          f"RMS {math.sqrt(s2 / n):.4e} m/s")
    print(f"last step against steady: largest difference {dl:.4e} m/s ({100 * dl / scale:.4f} %) at cell {wl}; "
          f"largest standard deviation of u over the averaging window {sd:.4e} m/s at cell {ws}")
    if a.max_rel is not None and not rel <= a.max_rel:
        sys.exit(f"FAIL: above --max_rel {a.max_rel:g}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("deck")
    d.add_argument("steady_deck")
    d.add_argument("out")
    d.set_defaults(fn=cmd_deck)
    c = sub.add_parser("cost")
    c.add_argument("logs", nargs="+")
    c.set_defaults(fn=cmd_cost)
    f = sub.add_parser("fields")
    f.add_argument("steady_plt")
    f.add_argument("unsteady_plt")
    f.add_argument("--every", type=int, default=1)
    f.add_argument("--max_rel", type=float)
    f.set_defaults(fn=cmd_fields)
    a = p.parse_args()
    if a.cmd == "fields" and a.every < 1:
        sys.exit("--every must be >= 1")
    if a.cmd == "fields" and a.max_rel is not None and not (a.max_rel > 0 and math.isfinite(a.max_rel)):
        sys.exit("--max_rel must be finite and > 0")
    a.fn(a)


if __name__ == "__main__":
    main()
