#!/usr/bin/env python3
"""Run the outer domain with two nests in one process and check it against separate runs.

    check_nests_one_process.py FCOMPARE -- COMMAND [ARGS ...]

COMMAND runs the outer deck (an argument ending in inputs_outer); the script
runs, each in its own directory under the current one:

1. joint/: COMMAND with nest.names and the corners of NESTS appended: the
   outer domain, then each nest from the outer run's final checkpoint;
2. outer/: COMMAND alone (no nest);
3. sep_<name>/: for each nest, COMMAND with inputs_outer replaced by
   inputs_nest (the nest written as a deck of its own) and the nest's mesh,
   bc.nest_parent = the final checkpoint of 2 and its own output names.

Checked:
a. every run exits 0 and reaches the steady verdict (each nest's too);
b. the joint log names each nest, in the order of nest.names, with its
   corners and the outer run's final checkpoint, after the outer run;
c. the joint run's outer plotfile is bit for bit that of run 2 (the nests
   change nothing in the outer domain), and each nest's final plotfile is
   bit for bit that of its separate run (the derived settings are those of
   inputs_nest), by fcompare with zero tolerances;
d. each nest's checkpoint and plotfile carry the "_<name>" suffix;
e. each nest's ground (z_phys_nd at k = 0) equals the outer domain's on the
   outer nodes it shares (within GROUND_TOL): the Agnesi hill stays at the
   outer domain's centre (800 m, 600 m), not at the nest's;
f. each nest's u, averaged over the 2 x 2 nest cells of an outer cell, is
   within DU_TOL of the outer u there (the nest refines the outer solution,
   it does not depart from it).

The exit code is the verdict.
"""
import glob
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402

# name, (x lo, y lo), (x hi, y hi) in m, on the 50 m outer cells. hill: outer
# cells 8 to 24 in x and 6 to 18 in y (32 x 24 nest cells of 25 m); edge:
# outer cells 24 to 30 in x and 6 to 12 in y (12 x 12 nest cells), touching
# hill along x = 1200 m and 2 outer cells from the outflow face.
NESTS = [("hill", (400, 300), (1200, 900)),
         ("edge", (1200, 300), (1500, 600))]
OUTER_DX = 50.0
NZ = 16
ZTOP = 800
GROUND_TOL = 1.0e-6     # m; the same formula at the same points
DU_TOL = {"hill": 0.12, "edge": 0.02}   # m/s; measured 0.0959 and 0.0115 (README)


def fail(msg):
    raise SystemExit("check_nests_one_process: " + msg)


def run(cmd, where):
    os.makedirs(where, exist_ok=True)
    proc = subprocess.run(cmd, cwd=where, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1800)
    print(proc.stdout)
    if proc.returncode != 0:
        fail(f"{' '.join(cmd)} exited {proc.returncode} in {where}")
    return proc.stdout


def last(where, prefix):
    """The directory prefix<step> with the largest step in where."""
    found = []
    for p in glob.glob(os.path.join(where, prefix + "[0-9]*")):
        m = re.fullmatch(re.escape(prefix) + r"(\d+)", os.path.basename(p))
        if m and os.path.isdir(p):
            found.append((int(m.group(1)), p))
    if not found:
        fail(f"no {prefix}<step> in {where}")
    return max(found)[1]


def fcompare(exe, a, b):
    proc = subprocess.run([exe, "--abs_tol", "0", "--rel_tol", "0", a, b],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(proc.stdout)
    if proc.returncode != 0:
        fail(f"{a} and {b} differ")
    print(f"check_nests_one_process: {a} and {b} agree bit for bit")


def read_nodes(prefix):
    """A one-component MultiFab written by VisMF as {(i, j, k): value} over its valid boxes."""
    head = open(prefix + "_H").read()
    boxes = re.findall(r"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)", head)
    fabs = re.findall(r"FabOnDisk: (\S+) (\d+)", head)
    out = {}
    for n, (name, off) in enumerate(fabs):
        lo, hi, data = erf_plotfile._read_fab(os.path.join(os.path.dirname(prefix), name), int(off), 1)
        nx, ny = hi[0] - lo[0] + 1, hi[1] - lo[1] + 1
        vlo = [int(x) for x in boxes[n][:3]]
        vhi = [int(x) for x in boxes[n][3:]]
        for k in range(vlo[2], vhi[2] + 1):
            for j in range(vlo[1], vhi[1] + 1):
                for i in range(vlo[0], vhi[0] + 1):
                    out[(i, j, k)] = data[(i - lo[0]) + nx * ((j - lo[1]) + ny * (k - lo[2]))]
    return out


def check_ground(outer_plt, nest_plt, name, lo):
    """e. The nest's ground nodes on the outer nodes it shares."""
    zo = read_nodes(os.path.join(outer_plt, "Level_0", "z_phys_nd"))
    zn = read_nodes(os.path.join(nest_plt, "Level_0", "z_phys_nd"))
    o = [round(lo[d] / OUTER_DX) for d in range(2)]
    worst, count = 0.0, 0
    for (i, j, k), z in zn.items():
        if k != 0 or i % 2 or j % 2:
            continue
        zz = zo[(o[0] + i // 2, o[1] + j // 2, 0)]
        worst = max(worst, abs(z - zz))
        count += 1
    print(f"check_nests_one_process: nest.{name}: max |ground - outer ground| = {worst:.3e} m on {count} shared nodes")
    if count == 0 or worst > GROUND_TOL:
        fail(f"nest.{name}: the ground differs from the outer domain's by {worst:.3e} m (> {GROUND_TOL:g} m)")


def check_u(outer_plt, nest_plt, name, lo):
    """f. The nest's u averaged to the outer cells against the outer u."""
    _, f_outer = erf_plotfile.read_fields(outer_plt, ["u"])
    _, f_nest = erf_plotfile.read_fields(nest_plt, ["u"])
    uo, un = f_outer["u"], f_nest["u"]
    o = [round(lo[d] / OUTER_DX) for d in range(2)]
    worst = 0.0
    for i in range(len(un) // 2):
        for j in range(len(un[0]) // 2):
            for k in range(len(un[0][0])):
                avg = 0.25 * (un[2 * i][2 * j][k] + un[2 * i + 1][2 * j][k]
                              + un[2 * i][2 * j + 1][k] + un[2 * i + 1][2 * j + 1][k])
                worst = max(worst, abs(avg - uo[o[0] + i][o[1] + j][k]))
    print(f"check_nests_one_process: nest.{name}: max |<u>_nest - u_outer| = {worst:.4f} m/s")
    if worst > DU_TOL[name]:
        fail(f"nest.{name}: the nest u departs from the outer u by {worst:.4f} m/s (> {DU_TOL[name]:g} m/s)")


def main():
    if "--" not in sys.argv or sys.argv.index("--") != 2:
        fail("usage: check_nests_one_process.py FCOMPARE -- COMMAND [ARGS ...]")
    exe = sys.argv[1]
    cmd = sys.argv[3:]
    deck = [i for i, a in enumerate(cmd) if a.endswith("inputs_outer")]
    if len(deck) != 1:
        fail("COMMAND must name inputs_outer once")
    here = os.getcwd()
    for d in ["joint", "outer"] + ["sep_" + n for n, _, _ in NESTS]:
        shutil.rmtree(os.path.join(here, d), ignore_errors=True)

    nest_keys = ["nest.names=" + " ".join(n for n, _, _ in NESTS)]
    for n, lo, hi in NESTS:
        nest_keys += [f"nest.{n}.lo={lo[0]} {lo[1]}", f"nest.{n}.hi={hi[0]} {hi[1]}"]
    joint = run(cmd + nest_keys, os.path.join(here, "joint"))
    outer = run(cmd, os.path.join(here, "outer"))

    # b. The nests in order, after the outer run, each from its final checkpoint.
    outer_chk = os.path.basename(last(os.path.join(here, "joint"), "outer_chk"))
    pos = 0
    for n, lo, hi in NESTS:
        line = f"gust: nest.{n}: "
        i = joint.find(line, pos)
        if i < 0:
            fail(f"the joint log does not name nest.{n} after the outer run (or after the nest before it)")
        text = joint[i:joint.find("\n", i)]
        print("check_nests_one_process: " + text)
        if outer_chk not in text:
            fail(f"nest.{n} does not start from the outer run's final checkpoint {outer_chk}: {text}")
        for v in (*lo, *hi):
            if not re.search(r"\b" + re.escape(f"{v:g}") + r"\b", text):
                fail(f"the line for nest.{n} does not give its corner value {v:g}: {text}")
        pos = i + len(line)
    if joint.find("gust: nest.") < joint.find(outer_chk):
        fail("a nest is named before the outer run wrote its final checkpoint")

    # c. The outer domain is that of the run without nests.
    fcompare(exe, last(os.path.join(here, "joint"), "outer_plt"), last(os.path.join(here, "outer"), "outer_plt"))

    parent = last(os.path.join(here, "outer"), "outer_chk")
    nest_deck = os.path.join(os.path.dirname(cmd[deck[0]]), "inputs_nest")
    for n, lo, hi in NESTS:
        cells = [round(2 * (hi[d] - lo[d]) / OUTER_DX) for d in range(2)]
        sep = list(cmd)
        sep[deck[0]] = nest_deck
        sep += [f"mesh.n_cell={cells[0]} {cells[1]} {NZ}", f"mesh.prob_lo={lo[0]} {lo[1]} 0",
                f"mesh.prob_hi={hi[0]} {hi[1]} {ZTOP}", "bc.nest_parent=" + parent,
                f"io.plot_file=outer_plt_{n}", f"io.chk_file=outer_chk_{n}"]
        run(sep, os.path.join(here, "sep_" + n))
        # d. The suffixed names.
        last(os.path.join(here, "joint"), f"outer_chk_{n}")
        nest_plt = last(os.path.join(here, "joint"), f"outer_plt_{n}")
        fcompare(exe, nest_plt, last(os.path.join(here, "sep_" + n), f"outer_plt_{n}"))
        outer_plt = last(os.path.join(here, "joint"), "outer_plt")
        check_ground(outer_plt, nest_plt, n, lo)
        check_u(outer_plt, nest_plt, n, lo)

    # a. The steady verdict: the outer run and each nest in the joint log.
    verdict = re.compile(r"^gust: steady after \d+ iterations", re.M)
    n_outer, n_joint = len(verdict.findall(outer)), len(verdict.findall(joint))
    if n_outer != 1 or n_joint != 1 + len(NESTS):
        fail(f"expected the steady verdict once in the outer log and {1 + len(NESTS)} times in the joint log; "
             f"found {n_outer} and {n_joint}")
    print("check_nests_one_process: PASS")


if __name__ == "__main__":
    main()
