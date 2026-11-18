#!/usr/bin/env python3
"""Two-way nesting (nest.feedback = defect, relax, relax_consistent or inject) on a steep hill against the uniform fine run.

    check_nests_two_way.py [--relax | --relax-consistent | --inject | --feedback-flux] -- COMMAND [ARGS ...]

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
   largest |S| per exchange, and "joint run: steady"; on the last exchange's
   line the feedback's change summed over steady.monitor_window exchanges
   (the "(window ...)" value) is finite and at most MONITOR_REL_TOL, the
   default steady.monitor_rel_tol the decks use (JointRun.H);
b. the one-way plotfiles are written before the joint ones (two plotfiles
   each for the outer domain and the nest);
c. against the reference (RMS of u over the nest, and over the outer cells
   outside the footprint, the reference averaged over 2 x 2 cells):
   the joint nest's error is at most NEST_RATIO of the one-way nest's,
   and the joint outer error outside the footprint at most OUTER_RATIO of
   the one-way one's (measured: README.md).

With --relax (PR-71, nest.feedback = relax) the script runs ref/ and joint/
as above, then relax1/, relax4/ and relax16/: COMMAND with the nest hill and
nest.feedback = relax at nest.feedback_factor = 1, 4 and 16, and
relax_default/: relax without nest.feedback_factor. Checked:
a. as above for each relax run (its own start line);
d. the last logged RMS of u of the parent minus the restricted nest under the
   nest (where the relaxation's rate is > 0) falls as the factor grows, and at
   16 is at most RELAX_FALL of its value at 1 (measured: README.md);
e. against the reference, the nest and outer errors (as in c.) of the run at
   factor 1 (the default before tracker N-114, item 13): at most
   RELAX_NEST_RATIO and RELAX_OUTER_RATIO of the one-way ones. At 4 and 16
   they are printed, not checked: they grow with the factor (README.md, with
   the defect run's beside them);
h. old against new default: relax_default/, the same run without
   nest.feedback_factor, names factor RELAX_DEFAULT (0.5) on its start line,
   and its nest error is at most RELAX_DEFAULT_NEST of the factor-1 run's,
   its outer error at most RELAX_DEFAULT_OUTER of it (README.md).

With --feedback-flux (tracker N-114, item 1; nest.feedback_flux, old against new)
the script runs ref/ as above, then, with nest.feedback = defect,
flux_outer/ and flux_nest/: the nest crest alone (CREST: its downstream face
at x = 900 m, on the lee slope) with nest.feedback_flux = outer (the default,
the PR-70 source alone) and nest (the ring rows take the nest's face fluxes),
and flux_hill/: the nest hill with nest.feedback_flux = nest. Checked:
a. as above for each run; the nest runs log the first correction
   ("nest.feedback_flux = nest: ... first max |correction|") and end their
   joint lines with the ring's max |C| ("; ring max |C| u ..."), the outer
   run does neither;
f. CREST with outer: the joint nest error is at least FLUX_OLD_RATIO of the
   one-way one (the defect this option is for); with nest: the joint nest and
   outer errors are at most FLUX_NEST_RATIO and FLUX_OUTER_RATIO of the
   one-way ones, and at most FLUX_GAIN of the joint errors with outer;
g. the hill with nest: the bounds of c. (NEST_RATIO, OUTER_RATIO).
Measured: README.md.

With --inject (tracker row N-114, item 2, nest.feedback = inject) the script
runs ref/ and joint/ as above, then inject/: COMMAND with the nest hill and
nest.feedback = inject at the default nest.feedback_factor (0.5). Checked:
a. as above (its own start line and per-exchange line);
f. against the reference, the nest and outer errors (as in c.): at most
   INJECT_NEST_RATIO and INJECT_OUTER_RATIO of the one-way ones. The defect
   run's are printed beside them (README.md).

With --relax-consistent (tracker row N-114, item 12, nest.feedback =
relax_consistent) the script runs ref/ and joint/ as above, then relax1/
(nest.feedback = relax at nest.feedback_factor = 1) and rc1/ and rc4/
(relax_consistent at factors 1 and 4). Checked:
a. as above for each run (its own start line, which for relax_consistent
   ends with "with the weighted residual of the target");
h. the last logged RMS of u of the parent minus the restricted nest under the
   nest: at factor 1 relax_consistent's is at most RC_DIFF of relax's, and it
   falls from factor 1 to 4;
i. against the reference, the nest and outer errors (as in c.) of rc1/: at
   most RC_NEST_RATIO and RC_OUTER_RATIO of the one-way ones. Those of rc4/
   are printed, not checked: they grow with the factor, as with relax
   (README.md, with the defect and relax runs' beside them).

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
# nest.feedback_flux (--feedback-flux): the crest nest alone, its downstream face at x = 900 m on the lee slope.
CREST = ["nest.names=crest", "nest.crest.lo=300 300", "nest.crest.hi=900 900"]
FOOT_CREST = {"crest": (6, 6, 12, 12)}
FLUX_OLD_RATIO = 2.0         # outer: joint nest error over the one-way one, measured 3.31 (u)
FLUX_NEST_RATIO = 1.6        # nest: measured 1.22 (u)
FLUX_OUTER_RATIO = 1.25      # nest: outer error outside the footprint over the one-way one, measured 0.96 (u)
FLUX_GAIN = 0.6              # nest over outer, joint errors: measured 0.37 (nest) and 0.35 (outer), u
FLUX_LOG = r"^gust: joint run: nest\.feedback_flux = nest: .*first max \|correction\|"
FLUX_RING = r"^gust: joint +\d+: max \|S\| .*; ring max \|C\| u \S+ v \S+ w \S+ theta \S+ k \S+$"
RELAX_FACTORS = (1, 4, 16)
RELAX_FALL = 0.35            # parent - restricted nest, u RMS at 16 over that at 1, measured 0.19
RELAX_NEST_RATIO = 0.5       # factor 1: nest error over the one-way one (u), measured 0.41
RELAX_OUTER_RATIO = 0.6      # factor 1: outer error outside the footprint over the one-way one (u), measured 0.35
RELAX_DEFAULT = 0.5          # the default nest.feedback_factor with relax (1 before tracker N-114, item 13)
RELAX_DEFAULT_NEST = 0.8     # default run's nest error over the factor-1 run's (u), measured 0.69
RELAX_DEFAULT_OUTER = 1.05   # default run's outer error over the factor-1 run's (u), measured 1.00
INJECT_FACTOR = 0.5          # the default nest.feedback_factor with inject
INJECT_NEST_RATIO = 0.55     # nest error over the one-way one (u), measured 0.513
INJECT_OUTER_RATIO = 0.42    # outer error outside the footprint over the one-way one (u), measured 0.377
RC_FACTORS = (1, 4)          # relax_consistent runs (--relax-consistent)
RC_DIFF = 0.25               # factor 1: parent - restricted nest, u RMS, relax_consistent over relax; measured 0.164
RC_NEST_RATIO = 0.6          # factor 1: nest error over the one-way one (u), measured 0.515
RC_OUTER_RATIO = 0.45        # factor 1: outer error outside the footprint over the one-way one (u), measured 0.361


MONITOR_REL_TOL = 1.0e-3     # steady.monitor_rel_tol (default; the decks do not set it)

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


def check(log, where, foots, r, nest_ratio, outer_ratio, factor=None, one_way=False, feedback="relax"):
    """a.-c. for the joint run in where, with the nests foots (name: outer cells i0, j0, ni, nj); factor: the
    nest.feedback_factor of a relax or inject run (feedback; None: defect). Returns the exchanges and the joint
    errors (nest name: u RMS, and "outer"; with one_way also "<name> one way" and "outer one way")."""
    # a.
    if factor is None:
        start, line = r"nest\.feedback = defect", r"max \|S\| u "
    else:
        start = r"nest\.feedback = %s, nest\.feedback_factor = %g" % (feedback, factor)
        line = (r"outer - injected nest under the nests: " if feedback == "inject"
                else r"outer - nest under the nests: ")
    for pat in (r"^gust: joint run \(%s\): the outer domain and %d nest" % (start, len(foots)),
                r"^gust: joint +1: " + line, r"^gust: joint run: steady after \d+ exchanges"):
        if not re.search(pat, log, re.M):
            fail(f"{where}: the joint log has no line matching {pat!r}")
    sums = re.findall(r"^gust: joint +\d+: .*?change \S+ \(window (\S+)\)", log, re.M)
    if not sums or not float(sums[-1]) <= MONITOR_REL_TOL:
        fail(f"{where}: the last exchange's feedback change over the window is {sums[-1:] or 'missing'}, "
             f"not at most {MONITOR_REL_TOL}")
    # b.
    outer = plotfiles(where, "outer_plt")
    if len(outer) != 2:
        fail(f"{where}: expected two outer plotfiles (one way, joint), found {outer}")
    # c.
    rc = 0.25 * (r[0::2, 0::2] + r[1::2, 0::2] + r[0::2, 1::2] + r[1::2, 1::2])
    mask = np.ones(rc.shape[:2], bool)
    errors = {}
    for name, (i0, j0, ni, nj) in foots.items():
        nest = plotfiles(where, "outer_plt_" + name)
        if len(nest) != 2:
            fail(f"{where}: expected two plotfiles of nest {name} (one way, joint), found {nest}")
        rn = r[2 * i0:2 * (i0 + ni), 2 * j0:2 * (j0 + nj)]
        ow, jt = rms(u_of(nest[0]) - rn), rms(u_of(nest[1]) - rn)
        print(f"check_nests_two_way: {where}: u RMS against the 25 m run over nest {name}: one way {ow:.4e}, "
              f"joint {jt:.4e} (ratio {jt / ow:.3f})")
        errors[name] = jt
        if one_way:
            errors[name + " one way"] = ow
        if not jt <= nest_ratio * ow:
            fail(f"{where}: the joint error of nest {name} is not at most {nest_ratio} of the one-way one")
        mask[i0:i0 + ni, j0:j0 + nj] = False
    m3 = np.broadcast_to(mask[:, :, None], rc.shape)
    ow, jt = rms((u_of(outer[0]) - rc)[m3]), rms((u_of(outer[1]) - rc)[m3])
    print(f"check_nests_two_way: {where}: u RMS against the 25 m run over the outer cells outside the footprints: "
          f"one way {ow:.4e}, joint {jt:.4e} (ratio {jt / ow:.3f})")
    if not jt <= outer_ratio * ow:
        fail(f"{where}: the joint outer error outside the footprints is not at most {outer_ratio} of the one-way one")
    errors["outer"] = jt
    if one_way:
        errors["outer one way"] = ow
    exchanges = int(re.search(r"^gust: joint run: steady after (\d+) exchanges", log, re.M).group(1))
    return exchanges, errors


def last_difference(log, where):
    """The last logged RMS of u of the parent minus the restricted nest (relax runs)."""
    found = re.findall(r"^gust: joint +\d+: outer - nest under the nests: .*?RMS u (\S+)", log, re.M)
    if not found:
        fail(f"{where}: no 'outer - nest under the nests' line")
    return float(found[-1].rstrip(","))


def feedback_flux(cmd, r):
    """--feedback-flux: f. and g. of the docstring."""
    inf = float("inf")
    runs = {}
    for where, extra, foots in (("flux_outer", CREST + ["nest.feedback_flux=outer"], FOOT_CREST),
                                ("flux_nest", CREST + ["nest.feedback_flux=nest"], FOOT_CREST),
                                ("flux_hill", NEST + ["nest.feedback_flux=nest"], {"hill": FOOT})):
        log = run(cmd + extra + ["nest.feedback=defect"], where)
        logged = re.search(FLUX_LOG, log, re.M) is not None
        if logged != where.endswith(("nest", "hill")):
            fail(f"{where}: the first-correction line is {'' if logged else 'not '}in the log")
        ring = re.search(FLUX_RING, log, re.M) is not None
        if ring != logged:
            fail(f"{where}: the joint lines {'' if ring else 'do not '}end with the ring's max |C|")
        hill = where == "flux_hill"
        runs[where] = check(log, where, foots, r, NEST_RATIO if hill else inf, OUTER_RATIO if hill else inf,
                            one_way=True)
    (_, old), (_, new) = runs["flux_outer"], runs["flux_nest"]
    ratios = {"old nest": old["crest"] / old["crest one way"], "new nest": new["crest"] / new["crest one way"],
              "new outer": new["outer"] / new["outer one way"]}
    print("check_nests_two_way: feedback_flux, crest to x = 900 m: joint over one-way, u: outer "
          f"{ratios['old nest']:.3f}; nest {ratios['new nest']:.3f} (outer domain {ratios['new outer']:.3f}); nest over "
          f"outer, joint: nest {new['crest'] / old['crest']:.3f}, outer domain {new['outer'] / old['outer']:.3f}")
    for ok, text in ((ratios["old nest"] >= FLUX_OLD_RATIO, f"with outer the joint nest error over the one-way one "
                      f"is {ratios['old nest']:.3f}, not at least {FLUX_OLD_RATIO}"),
                     (ratios["new nest"] <= FLUX_NEST_RATIO, f"with nest the joint nest error over the one-way one is "
                      f"{ratios['new nest']:.3f}, not at most {FLUX_NEST_RATIO}"),
                     (ratios["new outer"] <= FLUX_OUTER_RATIO, f"with nest the joint outer error over the one-way one "
                      f"is {ratios['new outer']:.3f}, not at most {FLUX_OUTER_RATIO}"),
                     (new["crest"] <= FLUX_GAIN * old["crest"] and new["outer"] <= FLUX_GAIN * old["outer"],
                      f"the joint errors with nest are not at most {FLUX_GAIN} of those with outer")):
        if not ok:
            fail(text)
    print("check_nests_two_way: passed")
    return 0


def relax_consistent(cmd, r, n_def, e_def):
    """--relax-consistent: h. and i. of the docstring."""
    print(f"check_nests_two_way: defect: {n_def} exchanges, u RMS against the 25 m run: nest {e_def['hill']:.4e}, "
          f"outer {e_def['outer']:.4e}")
    inf = float("inf")
    log = run(cmd + NEST + ["nest.feedback=relax", "nest.feedback_factor=1"], "relax1")
    n, e = check(log, "relax1", {"hill": FOOT}, r, inf, inf, factor=1)
    d_relax = last_difference(log, "relax1")
    print(f"check_nests_two_way: relax factor 1: {n} exchanges, parent - restricted nest u RMS {d_relax:.4e}; "
          f"u RMS against the 25 m run: nest {e['hill']:.4e}, outer {e['outer']:.4e}")
    diff = []
    for f in RC_FACTORS:
        where = f"rc{f}"
        log = run(cmd + NEST + ["nest.feedback=relax_consistent", f"nest.feedback_factor={f}"], where)
        if not re.search(r"^gust: joint run \(nest\.feedback = relax_consistent, .*with the weighted residual of "
                         r"the target$", log, re.M):
            fail(f"{where}: the start line does not name the weighted residual of the target")
        bound = f == 1
        n, e = check(log, where, {"hill": FOOT}, r, RC_NEST_RATIO if bound else inf,
                     RC_OUTER_RATIO if bound else inf, factor=f, feedback="relax_consistent")
        diff.append(last_difference(log, where))
        print(f"check_nests_two_way: relax_consistent factor {f}: {n} exchanges, parent - restricted nest u RMS "
              f"{diff[-1]:.4e} ({diff[-1] / d_relax:.3f} of relax at 1); u RMS against the 25 m run: "
              f"nest {e['hill']:.4e}, outer {e['outer']:.4e}")
    # h.
    if not diff[0] <= RC_DIFF * d_relax:
        fail(f"at factor 1 the parent - restricted nest difference is {diff[0] / d_relax:.3f} of relax's, "
             f"not at most {RC_DIFF}")
    if not diff[1] < diff[0]:
        fail(f"the parent - restricted nest difference does not fall from factor 1 to 4: {diff}")
    print("check_nests_two_way: passed")
    return 0


def main():
    if "--" not in sys.argv:
        print(__doc__)
        return 2
    opts = sys.argv[1:sys.argv.index("--")]
    relax, inject, flux = "--relax" in opts, "--inject" in opts, "--feedback-flux" in opts
    consistent = "--relax-consistent" in opts
    if relax + inject + flux + consistent > 1:
        fail("--relax, --relax-consistent, --inject and --feedback-flux are separate runs")
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
    if flux:
        return feedback_flux(cmd, r)
    n_def, e_def = check(run(cmd + NEST + ["nest.feedback=defect"], "joint"), "joint", {"hill": FOOT}, r, NEST_RATIO,
                         OUTER_RATIO)
    if inject:
        print(f"check_nests_two_way: defect: {n_def} exchanges, u RMS against the 25 m run: nest {e_def['hill']:.4e}, "
              f"outer {e_def['outer']:.4e}")
        log = run(cmd + NEST + ["nest.feedback=inject"], "inject")
        n, e = check(log, "inject", {"hill": FOOT}, r, INJECT_NEST_RATIO, INJECT_OUTER_RATIO, factor=INJECT_FACTOR,
                     feedback="inject")
        print(f"check_nests_two_way: inject factor {INJECT_FACTOR:g}: {n} exchanges; u RMS against the 25 m run: "
              f"nest {e['hill']:.4e}, outer {e['outer']:.4e}")
        print("check_nests_two_way: passed")
        return 0
    if consistent:
        return relax_consistent(cmd, r, n_def, e_def)
    if relax:
        diff, errs = [], {}
        print(f"check_nests_two_way: defect: {n_def} exchanges, u RMS against the 25 m run: nest {e_def['hill']:.4e}, "
              f"outer {e_def['outer']:.4e}")
        for f in RELAX_FACTORS:
            where = f"relax{f}"
            log = run(cmd + NEST + ["nest.feedback=relax", f"nest.feedback_factor={f}"], where)
            bound = f == 1
            n, e = check(log, where, {"hill": FOOT}, r, RELAX_NEST_RATIO if bound else float("inf"),
                         RELAX_OUTER_RATIO if bound else float("inf"), factor=f)
            diff.append(last_difference(log, where))
            errs[f] = e
            print(f"check_nests_two_way: relax factor {f}: {n} exchanges, parent - restricted nest u RMS "
                  f"{diff[-1]:.4e}; u RMS against the 25 m run: nest {e['hill']:.4e}, outer {e['outer']:.4e}")
        # h.
        log = run(cmd + NEST + ["nest.feedback=relax"], "relax_default")
        n, e = check(log, "relax_default", {"hill": FOOT}, r, float("inf"), float("inf"), factor=RELAX_DEFAULT)
        gain = {k: e[k] / errs[1][k] for k in ("hill", "outer")}
        print(f"check_nests_two_way: relax default (factor {RELAX_DEFAULT:g}): {n} exchanges, parent - restricted "
              f"nest u RMS {last_difference(log, 'relax_default'):.4e}; u RMS against the 25 m "
              f"run: nest {e['hill']:.4e}, outer {e['outer']:.4e}; over factor 1: nest {gain['hill']:.3f}, outer "
              f"{gain['outer']:.3f}")
        if not (gain["hill"] <= RELAX_DEFAULT_NEST and gain["outer"] <= RELAX_DEFAULT_OUTER):
            fail(f"the default factor's errors over the factor-1 run's are nest {gain['hill']:.3f}, outer "
                 f"{gain['outer']:.3f}; not at most {RELAX_DEFAULT_NEST} and {RELAX_DEFAULT_OUTER}")
        # d.
        if not all(b < a for a, b in zip(diff, diff[1:])):
            fail(f"the parent - restricted nest difference does not fall with the factor: {diff}")
        if not diff[-1] <= RELAX_FALL * diff[0]:
            fail(f"the difference at factor {RELAX_FACTORS[-1]} is {diff[-1] / diff[0]:.3f} of that at "
                 f"{RELAX_FACTORS[0]}, not at most {RELAX_FALL}")
        print("check_nests_two_way: passed")
        return 0
    check(run(cmd + NESTS2 + ["nest.feedback=defect"], "joint2"), "joint2", FOOTS2, r, NEST2_RATIO, OUTER2_RATIO)
    print("check_nests_two_way: passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
