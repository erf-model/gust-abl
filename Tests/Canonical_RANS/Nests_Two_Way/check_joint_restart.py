#!/usr/bin/env python3
"""A joint run stopped and restarted continues bit for bit.

    check_joint_restart.py [--relax | --relax-consistent | --inject] --stop M -- COMMAND [ARGS ...]

COMMAND runs the outer deck (an argument ending in inputs_outer). The script
runs, each in its own directory under the current one:

1. straight/: COMMAND with the nests crest and lee and nest.feedback = defect
              (with --relax, --relax-consistent or --inject: the nest hill and
              nest.feedback = relax, relax_consistent or inject), to the steady
              verdict;
2. stopped/:  the same with steady.max_iterations = M: the joint run stops
              after M exchanges, not steady, and names the checkpoint to
              restart from;
3. stopped/:  the restart: the deck without its init.* lines (they are
              refused with io.restart), io.restart = that checkpoint.

Checked:
a. the straight run and the restart exit 0, the stopped run does not; the
   stopped run's message names the checkpoint, and its Joint record has the
   feedback, M exchanges, the stop test's window, the factor (the feedback's
   default; 0 with defect) and one line per nest, in
   nest.names order, each naming a checkpoint that exists; the checkpoint
   holds the last feedback field (joint_feedback_*);
b. the restart logs "restarted after exchange M" and reaches steady after the
   straight run's number of exchanges, each domain after the same number of
   iterations over the run;
c. the last plotfile of every domain, and the last checkpoint's state, node
   heights, closure history, nest faces and Joint record, are byte for byte
   the straight run's (the driver record differs in its wall times only and
   is not compared).

With --relax also, before the restart, each for one exchange with
nest.feedback_factor = 0.75 (the record's is 0.5) from a copy of stopped/:

d. other/: the restart logs both values and that the window fills again, and
   its own record holds M + 1 exchanges and the factor 0.75;
e. older/: the record rewritten as version 2 (no factor line): the restart
   logs that it holds no factor, restores the window, and writes a version 3
   record with the factor 0.75.

The exit code is the verdict.
"""
import filecmp
import glob
import os
import re
import shutil
import subprocess
import sys

NESTS2 = ["nest.names=crest lee", "nest.crest.lo=300 300", "nest.crest.hi=1000 900",
          "nest.lee.lo=1000 250", "nest.lee.hi=1500 950"]
NEST = ["nest.names=hill", "nest.hill.lo=400 300", "nest.hill.hi=1200 900"]


def fail(msg):
    raise SystemExit("check_joint_restart: " + msg)


def run(cmd, where, expect_ok, fresh=True):
    if fresh:
        shutil.rmtree(where, ignore_errors=True)
        os.makedirs(where)
    proc = subprocess.run(cmd, cwd=where, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1800)
    print(proc.stdout[-3000:])
    if (proc.returncode == 0) != expect_ok:
        fail(f"{' '.join(cmd)} in {where}: exit {proc.returncode}, expected {'0' if expect_ok else 'not 0'}")
    return proc.stdout


def last(where, prefix):
    """The last plotfile or checkpoint <prefix>NNNNN (digits right after the prefix)."""
    names = [p for p in glob.glob(os.path.join(where, prefix + "[0-9]*")) if re.fullmatch(r"\d+", p[len(where) + 1 + len(prefix):])]
    if not names:
        fail(f"no {prefix}* in {where}")
    return max(names, key=lambda p: int(p[len(where) + 1 + len(prefix):]))


def same_tree(a, b, skip=()):
    """Every file under a is under b with the same bytes, and b has no others (skip: base names left out)."""
    files = lambda d: sorted(os.path.relpath(os.path.join(r, f), d) for r, _, fs in os.walk(d) for f in fs if f not in skip)
    fa, fb = files(a), files(b)
    if fa != fb:
        fail(f"{a} and {b} hold other files: {sorted(set(fa) ^ set(fb))}")
    diff = [f for f in fa if not filecmp.cmp(os.path.join(a, f), os.path.join(b, f), shallow=False)]
    if diff:
        fail(f"{b} differs from {a} in {diff}")
    print(f"check_joint_restart: {os.path.basename(b)}: {len(fa)} files identical to the straight run's")


def one_exchange(where, rcmd, chk, stop, older=False):
    """A restart for one exchange with another factor, from a copy of stopped/; its log."""
    shutil.rmtree(where, ignore_errors=True)
    shutil.copytree("stopped", where)
    if older:
        path = os.path.join(where, chk, "Joint")
        with open(path) as f:
            rec = f.read().split("\n")
        with open(path, "w") as f:
            f.write("\n".join(["GUST_JOINT 2"] + rec[1:5] + rec[6:]))
    log = run(rcmd + ["nest.feedback_factor=0.75", f"steady.max_iterations={stop + 1}"], where, False, fresh=False)
    if "the joint run did not reach steady" not in log:
        fail(f"the restart in {where} does not end with the joint run's not-steady message")
    with open(os.path.join(last(where, "outer_chk"), "Joint")) as f:
        rec = f.read().split("\n")
    if rec[0] != "GUST_JOINT 3" or rec[2] != f"exchanges {stop + 1}" or float.fromhex(rec[5].split()[1]) != 0.75:
        fail(f"the record the restart in {where} writes does not hold {stop + 1} exchanges and the factor 0.75")
    return log


def counts(log):
    """The joint run's exchanges and each domain's iterations over the run."""
    m = re.search(r"joint run: steady after (\d+) exchanges", log)
    if m is None:
        fail("no 'joint run: steady after N exchanges' line")
    its = dict(re.findall(r"joint run: (\S+(?: domain)?) steady after \d+ joint iterations \((\d+) over the run\)", log))
    return int(m.group(1)), its


def main():
    if "--" not in sys.argv:
        fail("usage: check_joint_restart.py [--relax | --relax-consistent | --inject] --stop M -- COMMAND [ARGS ...]")
    opts = sys.argv[1:sys.argv.index("--")]
    cmd = sys.argv[sys.argv.index("--") + 1:]
    relax = "--relax" in opts
    consistent = "--relax-consistent" in opts
    inject = "--inject" in opts
    if relax + consistent + inject > 1:
        fail("--relax, --relax-consistent and --inject are separate runs")
    if "--stop" not in opts or opts.index("--stop") + 1 >= len(opts):
        fail("--stop M is required")
    stop = int(opts[opts.index("--stop") + 1])
    decks = [i for i, a in enumerate(cmd) if a.endswith("inputs_outer")]
    if len(decks) != 1:
        fail("COMMAND must name the deck inputs_outer once")
    feedback = "relax" if relax else "relax_consistent" if consistent else "inject" if inject else "defect"
    nests = (NEST if relax or consistent or inject else NESTS2) + [f"nest.feedback={feedback}"]
    names = nests[0].split("=", 1)[1].split()

    log = run(cmd + nests, "straight", True)
    n_straight, it_straight = counts(log)
    if n_straight <= stop:
        fail(f"the straight run is steady after {n_straight} exchanges, not after --stop {stop}")

    log = run(cmd + nests + [f"steady.max_iterations={stop}"], "stopped", False)
    m = re.search(r"restart with io\.restart = (\S+) ", log)
    if m is None:
        fail("the stopped run does not name the checkpoint to restart from")
    chk = m.group(1)
    with open(os.path.join("stopped", chk, "Joint")) as f:
        rec = f.read().split("\n")
    print("check_joint_restart: Joint record:\n" + "\n".join(rec))
    if rec[:3] != ["GUST_JOINT 3", f"feedback {feedback}", f"exchanges {stop}"] or not rec[3].startswith("outer "):
        fail(f"the Joint record does not start with the version, feedback {feedback} and {stop} exchanges")
    win = rec[4].split()
    if len(win) < 2 or win[0] != "window" or int(win[1]) != len(win) - 2 or int(win[1]) < 1:
        fail("the Joint record's fifth line is not a window line with its changes")
    fac = rec[5].split()
    # The decks give no nest.feedback_factor: the defaults (Nests.H, NestInject.H), and 0 with defect.
    want = {"defect": 0.0, "relax": 0.5, "relax_consistent": 1.0, "inject": 0.5}[feedback]
    if len(fac) != 2 or fac[0] != "factor" or float.fromhex(fac[1]) != want:
        fail(f"the Joint record's sixth line is not the factor {want} of {feedback}")
    for n in range(5):
        if not os.path.isfile(os.path.join("stopped", chk, "Level_0", f"joint_feedback_{('u', 'v', 'w', 'theta', 'tke')[n]}_H")):
            fail(f"the checkpoint {chk} holds no feedback field {n}")
    lines = [r.split() for r in rec[6:] if r]
    if [ln[1] for ln in lines] != names or any(ln[0] != "nest" or len(ln) != 4 for ln in lines):
        fail(f"the Joint record's nest lines are not {names} in order")
    for ln in lines:
        if not os.path.isfile(os.path.join("stopped", ln[3], "Header")):
            fail(f"the Joint record names {ln[3]}, which is not a checkpoint")

    deck = cmd[decks[0]]
    with open(deck) as f:
        kept = [ln for ln in f if not ln.lstrip().startswith("init.")]
    with open(os.path.join("stopped", "inputs_restart"), "w") as f:
        f.writelines(kept)
    rcmd = cmd[:decks[0]] + ["inputs_restart"] + cmd[decks[0] + 1:] + nests + [f"io.restart={chk}"]
    if relax:
        # One exchange from the same checkpoint with another factor: the window is not restored.
        log = one_exchange("other", rcmd, chk, stop)
        refused = "its joint record holds the feedback of nest.feedback_factor = 0.5, not 0.75, so the joint stop test"
        if refused not in log:
            fail("a restart with another nest.feedback_factor does not log that its window fills again")
        # The same from a record of version 2: it has no factor to compare, and the window is restored.
        log = one_exchange("older", rcmd, chk, stop, older=True)
        if "its joint record holds no nest.feedback_factor" not in log or "goes on with 0.75 and" not in log:
            fail("a restart from a version 2 record does not log that it holds no factor")
        if "fills again" in log:
            fail("a restart from a version 2 record does not restore its window")
    log = run(rcmd, "stopped", True, fresh=False)
    if "fills again" in log:
        fail("the restart with the recorded feedback and factor does not restore the window")
    if f"joint run: restarted after exchange {stop}\n" not in log:
        fail(f"the restart does not log 'restarted after exchange {stop}'")
    n_restart, it_restart = counts(log)
    if (n_restart, it_restart) != (n_straight, it_straight):
        fail(f"the restart is steady after {n_restart} exchanges {it_restart}, the straight run after {n_straight} {it_straight}")

    for p in ["outer_plt"] + [f"outer_plt_{n}" for n in names]:
        a = last("straight", p)
        b = os.path.join("stopped", os.path.basename(a))
        if last("stopped", p) != b:
            fail(f"the restart's last {p} is {last('stopped', p)}, the straight run's {a}")
        same_tree(a, b)
    for p in ["outer_chk"] + [f"outer_chk_{n}" for n in names]:
        a = last("straight", p)
        same_tree(a, os.path.join("stopped", os.path.basename(a)), skip=("Driver",))
    print(f"check_joint_restart: stopped after {stop} exchanges and restarted: steady after {n_straight} exchanges, "
          "bit for bit the straight run")


if __name__ == "__main__":
    main()
