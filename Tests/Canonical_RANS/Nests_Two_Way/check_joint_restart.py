#!/usr/bin/env python3
"""A joint run stopped and restarted continues bit for bit (tracker N-114, item 3).

    check_joint_restart.py [--relax] --stop M -- COMMAND [ARGS ...]

COMMAND runs the outer deck (an argument ending in inputs_outer). The script
runs, each in its own directory under the current one:

1. straight/: COMMAND with the nests crest and lee and nest.feedback = defect
              (with --relax: the nest hill and nest.feedback = relax), to the
              steady verdict;
2. stopped/:  the same with steady.max_iterations = M: the joint run stops
              after M exchanges, not steady, and names the checkpoint to
              restart from;
3. stopped/:  the restart: the deck without its init.* lines (they are
              refused with io.restart), io.restart = that checkpoint.

Checked:
a. the straight run and the restart exit 0, the stopped run does not; the
   stopped run's message names the checkpoint, and its Joint record has the
   feedback, M exchanges and one line per nest, in nest.names order, each
   naming a checkpoint that exists;
b. the restart logs "restarted after exchange M" and reaches steady after the
   straight run's number of exchanges, each domain after the same number of
   iterations over the run;
c. the last plotfile of every domain, and the last checkpoint's state, node
   heights, closure history, nest faces and Joint record, are byte for byte
   the straight run's (the driver record differs in its wall times only and
   is not compared).

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


def counts(log):
    """The joint run's exchanges and each domain's iterations over the run."""
    m = re.search(r"joint run: steady after (\d+) exchanges", log)
    if m is None:
        fail("no 'joint run: steady after N exchanges' line")
    its = dict(re.findall(r"joint run: (\S+(?: domain)?) steady after \d+ joint iterations \((\d+) over the run\)", log))
    return int(m.group(1)), its


def main():
    if "--" not in sys.argv:
        fail("usage: check_joint_restart.py [--relax] --stop M -- COMMAND [ARGS ...]")
    opts = sys.argv[1:sys.argv.index("--")]
    cmd = sys.argv[sys.argv.index("--") + 1:]
    relax = "--relax" in opts
    if "--stop" not in opts or opts.index("--stop") + 1 >= len(opts):
        fail("--stop M is required")
    stop = int(opts[opts.index("--stop") + 1])
    decks = [i for i, a in enumerate(cmd) if a.endswith("inputs_outer")]
    if len(decks) != 1:
        fail("COMMAND must name the deck inputs_outer once")
    nests = (NEST + ["nest.feedback=relax"]) if relax else (NESTS2 + ["nest.feedback=defect"])
    names = nests[0].split("=", 1)[1].split()
    feedback = "relax" if relax else "defect"

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
    if rec[:3] != ["GUST_JOINT 1", f"feedback {feedback}", f"exchanges {stop}"] or not rec[3].startswith("outer "):
        fail(f"the Joint record does not start with the version, feedback {feedback} and {stop} exchanges")
    lines = [r.split() for r in rec[4:] if r]
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
    log = run(rcmd, "stopped", True, fresh=False)
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
