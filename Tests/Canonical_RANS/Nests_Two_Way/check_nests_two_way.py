#!/usr/bin/env python3
"""Two-way nesting (nest.feedback = defect) on a steep hill against the uniform fine run.

    check_nests_two_way.py -- COMMAND [ARGS ...]

COMMAND runs the outer deck (an argument ending in inputs_outer); the script
runs, each in its own directory under the current one:

1. ref/:   the deck with mesh.n_cell = 64 48 16 (25 m cells, the nest's
           size everywhere), to the steady verdict;
2. joint/: COMMAND with the nest hill (NEST) and nest.feedback = defect:
           the outer domain and the nest one way (their first plotfiles),
           then the joint run (their last plotfiles);
3. joint2/: the same with the two touching nests crest and lee (NESTS2).

Checked:
a. both runs exit 0; the joint log has the joint run's start line, the
   largest |S| per exchange, and "joint run: steady";
b. the one-way plotfiles are written before the joint ones (two plotfiles
   each for the outer domain and the nest);
c. against the reference (RMS of u over the nest, and over the outer cells
   outside the footprint, the reference averaged over 2 x 2 cells):
   the joint nest's error is at most NEST_RATIO of the one-way nest's,
   and the joint outer error outside the footprint at most OUTER_RATIO of
   the one-way one's (measured: README.md).

The exit code is the verdict.
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

NEST = ["nest.names=hill", "nest.hill.lo=400 300", "nest.hill.hi=1200 900"]
FOOT = (8, 6, 16, 12)        # outer cells: i0, j0, ni, nj of the nest's footprint
NEST_RATIO = 0.5             # measured 0.24 (u)
OUTER_RATIO = 0.6            # measured 0.34 (u)
# Two sites under one outer domain that touch at x = 1000 m: the windward side and crest, and the lee. The
# crest nest's downstream face stands 1.33 half-widths behind the crest; with it at 900 m (on the lee slope) the
# joint nest error is 3.3 times the one-way one (README.md).
NESTS2 = ["nest.names=crest lee", "nest.crest.lo=300 300", "nest.crest.hi=1000 900",
          "nest.lee.lo=1000 250", "nest.lee.hi=1500 950"]
FOOTS2 = {"crest": (6, 6, 14, 12), "lee": (20, 5, 10, 14)}
NEST2_RATIO = 0.7            # measured 0.42 (crest) and 0.61 (lee), u
OUTER2_RATIO = 0.9           # measured 0.80 (u)


def fail(msg):
    raise SystemExit("check_nests_two_way: " + msg)


def run(cmd, where):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    proc = subprocess.run(cmd, cwd=where, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1800)
    print(proc.stdout[-4000:])
    if proc.returncode != 0:
        fail(f"{' '.join(cmd)} exited {proc.returncode} in {where}")
    return proc.stdout


def plotfiles(where, prefix):
    found = []
    for p in glob.glob(os.path.join(where, prefix + "[0-9]*")):
        m = re.fullmatch(re.escape(prefix) + r"(\d+)", os.path.basename(p))
        if m and os.path.isdir(p):
            found.append((int(m.group(1)), p))
    return [p for _, p in sorted(found)]


def u_of(plt):
    _, f = erf_plotfile.read_fields(plt, ["u"])
    return np.array(f["u"], dtype=float)


def rms(a):
    return float(np.sqrt(np.mean(a * a)))


def check(log, where, foots, r, nest_ratio, outer_ratio):
    """a.-c. for the joint run in where, with the nests foots (name: outer cells i0, j0, ni, nj)."""
    # a.
    for pat in (r"^gust: joint run \(nest\.feedback = defect\): the outer domain and %d nest" % len(foots),
                r"^gust: joint +1: max \|S\| u ", r"^gust: joint run: steady after \d+ exchanges"):
        if not re.search(pat, log, re.M):
            fail(f"{where}: the joint log has no line matching {pat!r}")
    # b.
    outer = plotfiles(where, "outer_plt")
    if len(outer) != 2:
        fail(f"{where}: expected two outer plotfiles (one way, joint), found {outer}")
    # c.
    rc = 0.25 * (r[0::2, 0::2] + r[1::2, 0::2] + r[0::2, 1::2] + r[1::2, 1::2])
    mask = np.ones(rc.shape[:2], bool)
    for name, (i0, j0, ni, nj) in foots.items():
        nest = plotfiles(where, "outer_plt_" + name)
        if len(nest) != 2:
            fail(f"{where}: expected two plotfiles of nest {name} (one way, joint), found {nest}")
        rn = r[2 * i0:2 * (i0 + ni), 2 * j0:2 * (j0 + nj)]
        ow, jt = rms(u_of(nest[0]) - rn), rms(u_of(nest[1]) - rn)
        print(f"check_nests_two_way: {where}: u RMS against the 25 m run over nest {name}: one way {ow:.4e}, "
              f"joint {jt:.4e} (ratio {jt / ow:.3f})")
        if not jt <= nest_ratio * ow:
            fail(f"{where}: the joint error of nest {name} is not at most {nest_ratio} of the one-way one")
        mask[i0:i0 + ni, j0:j0 + nj] = False
    m3 = np.broadcast_to(mask[:, :, None], rc.shape)
    ow, jt = rms((u_of(outer[0]) - rc)[m3]), rms((u_of(outer[1]) - rc)[m3])
    print(f"check_nests_two_way: {where}: u RMS against the 25 m run over the outer cells outside the footprints: "
          f"one way {ow:.4e}, joint {jt:.4e} (ratio {jt / ow:.3f})")
    if not jt <= outer_ratio * ow:
        fail(f"{where}: the joint outer error outside the footprints is not at most {outer_ratio} of the one-way one")


def main():
    if "--" not in sys.argv:
        print(__doc__)
        return 2
    cmd = sys.argv[sys.argv.index("--") + 1:]
    deck = next((a for a in cmd if a.endswith("inputs_outer")), None)
    if deck is None:
        fail("COMMAND has no inputs_outer")
    cmd = [os.path.abspath(a) if a == deck else a for a in cmd]
    ref_deck = os.path.abspath("inputs_ref")
    with open(deck) as fh:
        text = fh.read()
    text, n = re.subn(r"^mesh\.n_cell\s*=.*$", "mesh.n_cell = 64 48 16", text, flags=re.M)
    if n != 1:
        fail("inputs_outer has no single mesh.n_cell line")
    text = re.sub(r"^io\.(plot|chk)_file\s*=.*$", "", text, flags=re.M)
    with open(ref_deck, "w") as fh:
        fh.write(text + "\nio.plot_file = ref_plt\nio.chk_file = ref_chk\n")
    i = cmd.index(os.path.abspath(deck))
    log_ref = run(cmd[:i] + [ref_deck] + cmd[i + 1:], "ref")
    if not re.search(r"^gust: steady after \d+ iterations", log_ref, re.M):
        fail("the reference did not reach the steady verdict")
    ref = plotfiles("ref", "ref_plt")
    if not ref:
        fail("no reference plotfile")
    r = u_of(ref[-1])
    check(run(cmd + NEST + ["nest.feedback=defect"], "joint"), "joint", {"hill": FOOT}, r, NEST_RATIO, OUTER_RATIO)
    check(run(cmd + NESTS2 + ["nest.feedback=defect"], "joint2"), "joint2", FOOTS2, r, NEST2_RATIO, OUTER2_RATIO)
    print("check_nests_two_way: passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
