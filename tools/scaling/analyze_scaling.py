#!/usr/bin/env python3
"""Scaling plots and summary from run_scaling.sh's results.csv.

    python3 analyze_scaling.py results.csv --out plots

Writes, into --out:
  cpu_seconds_per_iter.png, cpu_strong_scaling.png, cpu_efficiency.png, cpu_throughput_per_rank.png,
  gpu_seconds_per_iter.png, gpu_strong_scaling.png, gpu_efficiency.png, gpu_throughput_per_gpu.png,
  node_cpu_vs_gpu.png, seconds_per_iter.png    (only when matplotlib is available)
  summary.csv            one row per point with the derived numbers
and, next to results.csv:
  readme_scaling.md      the README section: the CPU and GPU seconds-per-iteration
                         plots and a link to the documentation page
  scaling.rst            the documentation page (docs/scaling.rst): set-up,
                         minimum cells per rank and per GPU, every table and
                         every plot
For the PR, copy the PNGs to docs/figures/scaling/ and scaling.rst to docs/ (the
default paths of both); paste readme_scaling.md into README.md.

Throughput is cell-iterations per second: cells / (seconds per iteration). The
per-rank (per-GPU) throughput is that divided by the rank (GPU) count. The
minimum cells per rank (GPU) at a fraction f is the smallest load such that every
point at that load or larger keeps at least f of the best per-rank (per-GPU)
throughput measured.
Python 3.6 or later; matplotlib is optional.
"""
import argparse
import csv
import math
import os
import sys

# Categorical slots in fixed order, one per mesh (the same mesh has the same
# colour in every plot); markers as a second encoding.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "h"]
INK, MUTED, GRID = "#1f1f1e", "#6b6a63", "#e4e3dd"
FRACTIONS = (0.8, 0.5)


def mesh_name(r):
    return "%dx%dx%d" % (r["nx"], r["ny"], r["nz"])


def load(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            for k in ("nx", "ny", "nz", "cells", "nproc", "cells_per_proc"):
                r[k] = int(r[k])
            try:
                r["s_per_iter"] = float(r["s_per_iter"])
            except ValueError:
                r["s_per_iter"] = float("nan")
            r["ok"] = r["status"] == "ok" and r["s_per_iter"] > 0
            rows.append(r)
    # A point run again keeps its last ok result (its last result if none is ok).
    last = {}
    for r in rows:
        key = (r["device"], r["nx"], r["ny"], r["nz"], r["nproc"])
        if r["ok"] or key not in last or not last[key]["ok"]:
            last[key] = r
    return list(last.values())


def derive(rows):
    for r in rows:
        if r["ok"]:
            r["throughput"] = r["cells"] / r["s_per_iter"]
            r["throughput_per_proc"] = r["throughput"] / r["nproc"]
    for dev in ("cpu", "gpu"):
        for m in sorted({mesh_name(r) for r in rows if r["device"] == dev}):
            pts = sorted((r for r in rows if r["device"] == dev and mesh_name(r) == m and r["ok"]),
                         key=lambda r: r["nproc"])
            if not pts:
                continue
            base = pts[0]
            for r in pts:
                r["speedup"] = base["s_per_iter"] / r["s_per_iter"] * base["nproc"]
                r["efficiency"] = r["speedup"] / r["nproc"]
                r["base_nproc"] = base["nproc"]


def min_cells(points, fraction):
    """Smallest cells per proc such that every point at that load or larger keeps fraction of the peak."""
    pts = [r for r in points if r["ok"]]
    if not pts:
        return None, None
    peak = max(r["throughput_per_proc"] for r in pts)
    loads = sorted({r["cells_per_proc"] for r in pts})
    for c in loads:
        if all(r["throughput_per_proc"] >= fraction * peak for r in pts if r["cells_per_proc"] >= c):
            return c, peak
    return None, peak


def levelled(points):
    """Whether the per-device throughput levels off: the largest load is within
    5 % of the next one measured (best value at each load)."""
    best = {}
    for r in points:
        best[r["cells_per_proc"]] = max(best.get(r["cells_per_proc"], 0), r["throughput_per_proc"])
    loads = sorted(best)
    if len(loads) < 2:
        return "one load only"
    gain = best[loads[-1]] / best[loads[-2]] - 1
    if gain <= -0.05:
        return "no, falls (%+.0f %% at the largest load)" % (100 * gain)
    if gain < 0.05:
        return "yes (%+.0f %% at the largest load)" % (100 * gain)
    return "no (%+.0f %% at the largest load): add a larger mesh" % (100 * gain)


def fmt_cells(c):
    if c is None:
        return "not reached"
    for unit, v in (("M", 1 << 20), ("k", 1 << 10)):
        if c >= v:
            return "%g%s (%d)" % (c / v, unit, c)
    return str(c)


def fmt_rate(x):
    return "%.3g" % (x / 1e6)


PLOTS = [("cpu_seconds_per_iter.png", "CPU seconds per iteration against mesh size, one line per rank count."),
         ("gpu_seconds_per_iter.png", "GPU seconds per iteration against mesh size, one line per GPU count."),
         ("cpu_strong_scaling.png", "CPU strong scaling: speedup against the smallest rank count of each mesh."),
         ("cpu_efficiency.png", "CPU parallel efficiency."),
         ("cpu_throughput_per_rank.png", "CPU throughput per rank against cells per rank."),
         ("gpu_strong_scaling.png", "GPU strong scaling: speedup against the smallest GPU count of each mesh."),
         ("gpu_efficiency.png", "GPU parallel efficiency."),
         ("gpu_throughput_per_gpu.png", "GPU throughput per GPU against cells per GPU."),
         ("node_cpu_vs_gpu.png", "Total throughput: the most CPU ranks measured against each GPU count."),
         ("seconds_per_iter.png", "Seconds per iteration against mesh size: the most CPU ranks measured and each GPU count.")]

# Where the runs ran (set from the command line).
NODES = {"cpu": "a Kestrel CPU node (104 Intel Xeon Sapphire Rapids cores)",
         "gpu": "a Kestrel GPU node (4 × NVIDIA H100)"}


def max_cpu_ranks(rows):
    return max([r["nproc"] for r in rows if r["device"] == "cpu" and r["ok"]] or [0])


def threshold_groups(rows):
    groups = [("CPU, all rank counts", [r for r in rows if r["device"] == "cpu"]),
              ("CPU, %d ranks (most measured)" % max_cpu_ranks(rows),
               [r for r in rows if r["device"] == "cpu" and r["nproc"] == max_cpu_ranks(rows)])]
    for g in sorted({r["nproc"] for r in rows if r["device"] == "gpu"}):
        groups.append(("GPU, %d GPU%s" % (g, "s" if g > 1 else ""),
                       [r for r in rows if r["device"] == "gpu" and r["nproc"] == g]))
    groups.append(("GPU, all GPU counts", [r for r in rows if r["device"] == "gpu"]))
    out = []
    for name, pts in groups:
        ok = [r for r in pts if r["ok"]]
        if ok:
            c80, peak = min_cells(ok, 0.8)
            c50, _ = min_cells(ok, 0.5)
            out.append((name, len(ok), peak, c80, c50, levelled(ok)))
    return out


def setup_text(rows):
    ok = [r for r in rows if r["ok"]] or rows

    def distinct(k):
        return ", ".join(sorted({r[k] for r in ok if r[k]}))
    return distinct("gust_sha"), distinct("node"), distinct("iters_short"), distinct("iters_long")


def counts_text(rows):
    """'CPU runs on 1 to 64 MPI ranks of ..., and GPU runs on 1, 2 and 4 GPUs of ...', from what was measured."""
    parts = []
    for dev, unit in (("cpu", "MPI ranks"), ("gpu", "GPUs")):
        n = sorted({r["nproc"] for r in rows if r["device"] == dev and r["ok"]})
        if not n:
            continue
        if dev == "cpu":
            c = "%d to %d" % (n[0], n[-1]) if len(n) > 1 else str(n[0])
        else:
            c = ", ".join(str(x) for x in n[:-1]) + " and %d" % n[-1] if len(n) > 1 else str(n[0])
        parts.append("%s runs on %s %s of %s" % (dev.upper(), c, unit, NODES[dev]))
    return ", and ".join(parts)


def readme_md(rows, fig_dir, docs_url, written):
    """The README section: the two seconds-per-iteration plots and a pointer to the docs."""
    shas, nodes, _, _ = setup_text(rows)
    lines = ["## Performance", "",
             "Seconds per iteration on %s and on %s, the 3-D hill on " % (NODES["cpu"], NODES["gpu"]) +
             "power-of-2 meshes (GUST %s). The analysis, the minimum cells per rank and per GPU and the full "
             "tables are in the [scaling page](%s) of the documentation." % (shas, docs_url), ""]
    for f, alt in (("cpu_seconds_per_iter.png", "CPU seconds per iteration"),
                   ("gpu_seconds_per_iter.png", "GPU seconds per iteration")):
        if f in written:
            lines += ["![%s](%s/%s)" % (alt, fig_dir, f), ""]
    return "\n".join(lines)


def rst_table(header, body, widths):
    lines = [".. list-table::", "   :header-rows: 1", "   :widths: %s" % " ".join(str(w) for w in widths), ""]
    for row in [header] + body:
        lines.append("   * - %s" % row[0])
        lines += [("     - %s" % c).rstrip() for c in row[1:]]
    return lines + [""]


def scaling_rst(rows, fig_dir, written):
    """The documentation page: set-up, minimum cells, tables and every plot."""
    shas, nodes, ns, nl = setup_text(rows)
    title = "Performance and scaling"
    lines = [title, "=" * len(title), "",
             "Strong scaling of GUST: %s; one rank per GPU. GUST %s on %s." % (counts_text(rows), shas, nodes),
             "", "Set-up", "------", ""]
    lines += rst_table(["Item", "Value"], [
        ["Case", "``Tests/Canonical_RANS/Neutral_Hill_3D/inputs_steady_hill3d``: the 3-D Agnesi hill, periodic, "
                 "Coriolis on, z-line Jacobi sweeps (``steady.sweep = zjacobi``, 16 sweeps, CFL cap 4, as in the "
                 "deck) on both devices; the CPU default ``steady.sweep = lusgs`` was not timed. ``mesh.n_cell``, "
                 "``mesh.max_grid_size`` and ``steady.max_iterations`` are set for each run on the command line"],
        ["Meshes", ", ".join(sorted({mesh_name(r) for r in rows}, key=lambda m: [int(x) for x in m.split("x")]))],
        ["Boxes", "one per rank: x and y split in powers of 2 (``mesh.max_grid_size``), whole columns in z"],
        ["Launch", "CPU ``srun --cpu-bind=cores``; GPU ``srun --gpus-per-task=1 --gpu-bind=none``"],
        ["Timing", "seconds per iteration from the difference of a %s- and a %s-iteration run, using the "
                   "solver's stage time (``gust: stage account``); set-up, the first guess and CUDA start-up "
                   "are left out" % (ns, nl)],
        ["Throughput", "cell-iterations per second (cells / seconds per iteration), in millions (Mcell-it/s); per "
                       "rank or per GPU, divided by the count"],
    ], [20, 80])
    lines += ["Minimum cells per rank and per GPU", "----------------------------------", "",
              "The smallest load such that every measured point at that load or larger keeps 80 % (50 %) of the "
              "best per-rank (per-GPU) throughput. \"Levelled off\" says whether the best throughput at the "
              "largest load is within 5 % of that at the next largest load. If it is still rising, the best "
              "throughput is not yet reached and a larger mesh would show it; if it falls, larger loads are slower "
              "per cell. The two largest loads can come from meshes of different height (``nz``), so the test "
              "also mixes in the mesh shape.", ""]
    groups = threshold_groups(rows)
    if groups:
        lines += rst_table(["Device", "Points", "Best Mcell-it/s per rank or GPU", "80 % from (cells)",
                            "50 % from (cells)", "Levelled off"],
                           [[n, str(k), fmt_rate(p), fmt_cells(a), fmt_cells(b), lev]
                            for n, k, p, a, b, lev in groups], [24, 8, 14, 18, 18, 18])
    else:
        lines += ["No point finished, so there is no threshold to report.", ""]
    lines += ["On the CPU the per-rank throughput depends on two things. It falls as the box per rank grows "
              "(less of it fits in cache), and it falls as more ranks share the node's memory bandwidth. When it "
              "falls with load, the largest loads can stay below 80 % (50 %) of the best, and the table then says "
              "\"not reached\"; a larger mesh does not change that. The all-rank-count row mixes both effects; the row with the most "
              "ranks is the closest to a full node.", "",
              "An efficiency above 100 % means the base point (the smallest count measured for that mesh) ran "
              "slower per rank than the larger counts, for example a box that does not fit in cache; the "
              "per-rank throughput plot shows which.", ""]
    for dev, label, per in (("cpu", "CPU ranks", "rank"), ("gpu", "GPUs", "GPU")):
        pts = sorted((r for r in rows if r["device"] == dev),
                     key=lambda r: (r["cells"], r["nz"], r["ny"], r["nx"], r["nproc"]))
        if not pts:
            continue
        head = "%s runs" % dev.upper()
        lines += [head, "-" * len(head), ""]
        body = []
        for r in pts:
            if r["ok"]:
                body.append([mesh_name(r), str(r["cells"]), str(r["nproc"]), "%sx%s" % (r["boxes_x"], r["boxes_y"]),
                             str(r["cells_per_proc"]), "%.4g" % r["s_per_iter"], fmt_rate(r["throughput"]),
                             fmt_rate(r["throughput_per_proc"]), "%.3g" % r["speedup"],
                             "%.0f %%" % (100 * r["efficiency"])])
            else:
                body.append([mesh_name(r), str(r["cells"]), str(r["nproc"]), "%sx%s" % (r["boxes_x"], r["boxes_y"]),
                             str(r["cells_per_proc"]), r["status"] if r["status"] != "ok" else "invalid time",
                             "", "", "", ""])
        lines += rst_table(["Mesh", "Cells", label, "Boxes", "Cells per %s" % per, "s per iteration", "Mcell-it/s",
                            "per %s" % per, "Speedup", "Efficiency"], body, [12, 10, 7, 7, 10, 10, 10, 10, 8, 8])
        lines += ["Speedup and efficiency against the smallest %s count measured for each mesh." % label[:-1], ""]
    if written:
        lines += ["Plots", "-----", ""]
        for f, cap in PLOTS:
            if f in written:
                lines += [".. figure:: %s/%s" % (fig_dir, f), "   :width: 90%", "", "   %s" % cap, ""]
    lines += ["Reproducing", "-----------", "",
              "``run_scaling.sh`` (in a one-node allocation; ``CPU_TARGET=spr`` on a CPU node) builds GUST for CPU and "
              "CUDA, runs every point and writes ``results.csv``; ``analyze_scaling.py results.csv`` writes the plots "
              "into ``--out`` and this page and the README section (``scaling.rst``, ``readme_scaling.md``) next to "
              "``results.csv``; copy the plots to ``docs/figures/scaling/``, the page to ``docs/`` and the section "
              "into ``README.md``. Finished points are skipped when ``run_scaling.sh`` runs again. The "
              "measured ``results.csv`` and ``analyze_scaling.py`` are in ``tools/scaling/``; ``run_scaling.sh`` is "
              "not in the repository yet.", ""]
    return "\n".join(lines)


def plots(rows, out):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
    except ImportError:
        return set()
    plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
                         "xtick.color": MUTED, "ytick.color": MUTED, "axes.titlecolor": INK,
                         "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False})
    meshes = sorted({(r["cells"], r["nz"], mesh_name(r)) for r in rows})
    style = {m[2]: (SERIES[i % len(SERIES)], MARKERS[i % len(MARKERS)]) for i, m in enumerate(meshes)}
    int_fmt = FuncFormatter(lambda v, _: "%g" % v)
    cells_fmt = FuncFormatter(lambda v, _: fmt_cells(int(round(v))).split()[0] if v >= 1 else "")
    written = set()

    def done(fig, name):
        fig.tight_layout()
        fig.savefig(os.path.join(out, name))
        plt.close(fig)
        written.add(name)

    def plain_y(ax, subs=(2, 5)):
        # Plain numbers on a log axis; minor labels at subs when the
        # axis spans less than two decades.
        fmt = FuncFormatter(lambda v, _: "%g" % float("%.3g" % v))
        ax.yaxis.set_major_formatter(fmt)
        lo, hi = ax.get_ylim()
        if hi / lo < 100:
            ax.yaxis.set_minor_locator(LogLocator(base=10, subs=subs))
            ax.yaxis.set_minor_formatter(fmt)
        else:
            ax.yaxis.set_minor_formatter(NullFormatter())

    def axes(title, xlabel, ylabel, logx=True, logy=True):
        fig, ax = plt.subplots(figsize=(7.6, 4.4), dpi=150)
        ax.set_title(title, loc="left", fontsize=11, fontweight="bold")
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if logx:
            try:
                ax.set_xscale("log", base=2)
            except (TypeError, ValueError):
                ax.set_xscale("log", basex=2)   # matplotlib before 3.3
            ax.xaxis.set_major_formatter(int_fmt)
        if logy:
            ax.set_yscale("log")
        ax.grid(True, which="major", color=GRID, linewidth=0.8)
        ax.set_axisbelow(True)
        return fig, ax

    def save(fig, ax, name, legend=True):
        if legend:
            ax.legend(fontsize=8, title="Mesh", title_fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1))
        done(fig, name)

    def by_mesh(dev):
        for _, _, m in meshes:
            pts = sorted((r for r in rows if r["device"] == dev and mesh_name(r) == m and r["ok"]),
                         key=lambda r: r["nproc"])
            if pts:
                yield m, pts

    for dev, unit in (("cpu", "CPU ranks"), ("gpu", "GPUs")):
        if not any(r["device"] == dev and r["ok"] for r in rows):
            continue
        procs = sorted({r["nproc"] for r in rows if r["device"] == dev and r["ok"]})
        # Seconds per iteration against mesh size, one line per rank (GPU) count.
        fig, ax = axes("%s seconds per iteration" % dev.upper(), "Cells", "Seconds per iteration")
        for i, p in enumerate(procs):
            pts = sorted((r for r in rows if r["device"] == dev and r["nproc"] == p and r["ok"]),
                         key=lambda r: r["cells"])
            ax.plot([r["cells"] for r in pts], [r["s_per_iter"] for r in pts], color=SERIES[i % len(SERIES)],
                    marker=MARKERS[i % len(MARKERS)], markersize=6, linewidth=1.5, label=str(p))
        ax.xaxis.set_major_formatter(cells_fmt)
        plain_y(ax)
        ax.legend(fontsize=8, title=unit, title_fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1))
        done(fig, "%s_seconds_per_iter.png" % dev)
        # Strong scaling: speedup against the smallest count, with the ideal line.
        fig, ax = axes("%s strong scaling" % dev.upper(), unit, "Speedup")
        ax.plot(procs, procs, color=MUTED, linestyle="--", linewidth=1, label="ideal")
        for m, pts in by_mesh(dev):
            c, mk = style[m]
            ax.plot([r["nproc"] for r in pts], [r["speedup"] for r in pts], color=c, marker=mk,
                    markersize=6, linewidth=1.5, label="%s (%s cells)" % (m, fmt_cells(pts[0]["cells"]).split()[0]))
        plain_y(ax, subs=(2, 3, 4))
        save(fig, ax, "%s_strong_scaling.png" % dev)
        # Efficiency.
        fig, ax = axes("%s parallel efficiency" % dev.upper(), unit, "Efficiency (%)", logy=False)
        ax.axhline(80, color=MUTED, linestyle="--", linewidth=1)
        for m, pts in by_mesh(dev):
            c, mk = style[m]
            ax.plot([r["nproc"] for r in pts], [100 * r["efficiency"] for r in pts], color=c, marker=mk,
                    markersize=6, linewidth=1.5, label=m)
        ax.set_ylim(0, max(110, ax.get_ylim()[1]))
        save(fig, ax, "%s_efficiency.png" % dev)
        # Throughput per device against cells per device, with the 80 % line.
        per = "rank" if dev == "cpu" else "GPU"
        fig, ax = axes("%s throughput per %s" % (dev.upper(), per), "Cells per %s" % per,
                       "Mcell-it/s per %s" % per)
        ok = [r for r in rows if r["device"] == dev and r["ok"]]
        c80, peak = min_cells(ok, 0.8)
        ax.xaxis.set_major_formatter(cells_fmt)
        for m, pts in by_mesh(dev):
            c, mk = style[m]
            ax.plot([r["cells_per_proc"] for r in pts], [r["throughput_per_proc"] / 1e6 for r in pts],
                    color=c, marker=mk, markersize=6, linewidth=1.5, label=m)
        plain_y(ax)
        ax.axhline(0.8 * peak / 1e6, color=MUTED, linestyle="--", linewidth=1)
        ax.annotate("80 % of best", (min(r["cells_per_proc"] for r in ok), 0.8 * peak / 1e6),
                    textcoords="offset points", xytext=(2, 4), fontsize=8, color=MUTED)
        if c80 is not None:
            ax.axvline(c80, color=MUTED, linestyle=":", linewidth=1)
            ax.annotate("80 %% from %s cells" % fmt_cells(c80).split()[0], (c80, ax.get_ylim()[0]),
                        textcoords="offset points", xytext=(-4, 6), ha="right", fontsize=8, color=MUTED)
        save(fig, ax, "%s_throughput_per_%s.png" % (dev, per.lower()))

    # Total throughput against mesh size: the most CPU ranks measured and each GPU count.
    nmax = max_cpu_ranks(rows)
    lines = [("cpu", nmax, "%d CPU ranks" % nmax, SERIES[0], "o")]
    for i, g in enumerate(sorted({r["nproc"] for r in rows if r["device"] == "gpu"})):
        lines.append(("gpu", g, "%d GPU%s" % (g, "s" if g > 1 else ""), SERIES[(1 + i) % len(SERIES)],
                      MARKERS[(1 + i) % len(MARKERS)]))
    have = [(d, p, l, c, mk) for d, p, l, c, mk in lines
            if any(r["device"] == d and r["nproc"] == p and r["ok"] for r in rows)]
    if have:
        fig, ax = axes("CPU ranks against GPUs", "Cells", "Mcell-it/s")
        for d, p, label, c, mk in have:
            pts = sorted((r for r in rows if r["device"] == d and r["nproc"] == p and r["ok"]),
                         key=lambda r: r["cells"])
            ax.plot([r["cells"] for r in pts], [r["throughput"] / 1e6 for r in pts], color=c, marker=mk,
                    markersize=6, linewidth=1.5, label=label)
        ax.xaxis.set_major_formatter(cells_fmt)
        plain_y(ax)
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1))
        done(fig, "node_cpu_vs_gpu.png")
        # The same lines as wall time: seconds per iteration against mesh size.
        fig, ax = axes("Seconds per iteration", "Cells", "Seconds per iteration")
        for d, p, label, c, mk in have:
            pts = sorted((r for r in rows if r["device"] == d and r["nproc"] == p and r["ok"]),
                         key=lambda r: r["cells"])
            ax.plot([r["cells"] for r in pts], [r["s_per_iter"] for r in pts], color=c, marker=mk,
                    markersize=6, linewidth=1.5, label=label)
        ax.xaxis.set_major_formatter(cells_fmt)
        plain_y(ax)
        ax.legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.02, 1))
        done(fig, "seconds_per_iter.png")
    return written


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("results", help="the results.csv written by run_scaling.sh")
    ap.add_argument("--out", default="plots", help="directory for the PNGs and summary.csv")
    ap.add_argument("--readme-figures", default="docs/figures/scaling",
                    help="where README.md finds the plots (from the repository root)")
    ap.add_argument("--docs-figures", default="figures/scaling",
                    help="where docs/scaling.rst finds the plots (from docs/)")
    ap.add_argument("--docs-url", default="https://erf-model.github.io/gust-abl/scaling.html",
                    help="the scaling page, linked from the README section")
    ap.add_argument("--cpu-node", default=NODES["cpu"], help="where the CPU runs ran, for the text")
    ap.add_argument("--gpu-node", default=NODES["gpu"], help="where the GPU runs ran, for the text")
    a = ap.parse_args()
    NODES["cpu"], NODES["gpu"] = a.cpu_node, a.gpu_node
    rows = load(a.results)
    if not rows:
        sys.exit("no points in %s" % a.results)
    derive(rows)
    os.makedirs(a.out, exist_ok=True)
    written = plots(rows, a.out)
    keys = ["device", "nx", "ny", "nz", "cells", "nproc", "cells_per_proc", "s_per_iter", "throughput",
            "throughput_per_proc", "speedup", "efficiency", "base_nproc", "status"]
    with open(os.path.join(a.out, "summary.csv"), "w") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for r in sorted(rows, key=lambda r: (r["device"], r["cells"], r["nz"], r["ny"], r["nx"], r["nproc"])):
            w.writerow([r.get(k, "") for k in keys])
    base = os.path.dirname(os.path.abspath(a.results))
    with open(os.path.join(base, "readme_scaling.md"), "w") as f:
        f.write(readme_md(rows, a.readme_figures, a.docs_url, written).rstrip("\n") + "\n")
    with open(os.path.join(base, "scaling.rst"), "w") as f:
        f.write(scaling_rst(rows, a.docs_figures, written).rstrip("\n") + "\n")
    print("wrote %s/readme_scaling.md, %s/scaling.rst, %s/summary.csv%s" % (
        base, base, a.out, " and %d plots" % len(written) if written else " (no plots)"))


if __name__ == "__main__":
    main()
