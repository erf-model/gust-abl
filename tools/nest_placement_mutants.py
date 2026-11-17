#!/usr/bin/env python3
"""
The one-line mutants of tools/nest_placement.py, each run against its self-test
(tools/test_nest_placement.py) in a scratch copy of the two files. Every mutant must
make the self-test fail. The table is in docs/testing.rst.

    python3 tools/nest_placement_mutants.py [P1 P4 ...]
"""
import os
import shutil
import subprocess
import sys
import tempfile

F = "nest_placement.py"
M = [
 ("P1", "the angle rule off", "                if ang > a_.angle:", "                if False:"),
 ("P2", "fixed thresholds (the dx 200 m values), not scaled", "setattr(a, name, f * a.dx)", "setattr(a, name, f * 200.0)"),
 ("P3", "no steps", "            found = peaks + steps(prof, w, a_.step_rise, bases)", "            found = peaks"),
 ("P4", "only the first wind direction", "        for dr in a_.direction:", "        for dr in a_.direction[:1]:"),
 ("P5", "a bend flags the face again", "                    notes.append(f\"{what} bends", "                    faults.append(f\"{what} bends"),
 ("P6", "no warning when no position is clear", "        print(msg, file=sys.stderr)", "        pass"),
 ("P7", "the margin of before (4 dx)", '"margin": 10.0}', '"margin": 4.0}'),
 ("P8", "side faces flagged too", 'k != "side" and f["faults"]', 'f["faults"]'),
 ("P9", "the misfit gate off", "                if local <= a_.misfit:", "                if False:"),
 ("P10", "a check always exits 0", 'return 1 if any(ev["flagged"] for ev in results) else 0', "return 0"),
 ("P11", "the face-kind cutoff of before (0.3)", 'c > 0.5 else ("inflow" if c < -0.5', 'c > 0.3 else ("inflow" if c < -0.3'),
 ("P12", "terrain smaller than the domain accepted", "        fail(f\"--domain {x0:g}", "        print(f\"--domain {x0:g}"),
 ("P13", "nest corners off the outer nodes accepted", "        if not all(on_cells(v, o, a.dx) for v, o in zip(n,",
  "        if False and all(on_cells(v, o, a.dx) for v, o in zip(n,"),
]


def run(selected):
    here = os.path.dirname(os.path.abspath(__file__))
    src = open(os.path.join(here, F)).read()
    survived = []
    for tag, what, old, new in M:
        if selected and tag not in selected:
            continue
        if src.count(old) != 1:
            raise SystemExit(f"{tag}: the line to mutate is not found once in {F}: {old!r}")
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, F), "w").write(src.replace(old, new))
            shutil.copy(os.path.join(here, "test_nest_placement.py"), d)
            r = subprocess.run([sys.executable, "-m", "unittest", "test_nest_placement"], cwd=d,
                               capture_output=True, text=True, env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        out = r.stdout + r.stderr
        if r.returncode != 0 and "FAILED" not in out:
            raise SystemExit(f"{tag}: the self-test did not run:\n{out[-2000:]}")
        killed = r.returncode != 0
        print(f"{tag:4s} {'killed' if killed else 'SURVIVED':8s} {what}")
        if not killed:
            survived.append(tag)
    if survived:
        raise SystemExit(f"{len(survived)} mutant(s) survived: {' '.join(survived)}")


if __name__ == "__main__":
    run(set(sys.argv[1:]))
