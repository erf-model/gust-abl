#!/usr/bin/env python3
"""The Rampart_2D case: the ground step of the mesh report against the iterations to the steady verdict.

    check_rampart.py -- LAUNCHER... GUST DECK [FLAGS ...]

DECK is inputs_rampart. The script writes the two terrains with make_terrain.py and runs, each in its own
directory (emptied first):
  gentle2  the 750 m wide cliff, 2 m first layers   to the verdict (at most 600 iterations)
  steep8   the 170 m wide cliff, 8 m first layers   to the verdict (at most 4000 iterations)
  steep2   the 170 m wide cliff, 2 m first layers   one iteration (its 15383 iterations to the verdict are in README.md)
It checks:
  - the report's largest ground step between neighbouring columns, in first layers, against the closed form
    476 m or 119 m over the first layer at the higher node, dz0 (1 - h / z_top): 62.85, 63.86 and 255.45 (0.1 %);
  - the report's largest aspect ratio: 109.9 with 2 m layers, 27.5 with 8 m (0.1 %): gentle2 and steep2 have the same;
  - the warning "the ground steps by up to ..." on steep2 only (the other two are below the tested 64);
  - gentle2 is steady within 400 iterations and steep8 within 4000, and steep8 takes more than four times
    gentle2's iterations at the same ground step (260 and 2570 measured on 2 ranks).
The exit code is the verdict.
"""
import math
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_terrain as mt  # noqa: E402

RUNS = (("gentle2", "rampart_gentle.xyz", 750.0, 2.0, 600),
        ("steep8", "rampart_steep.xyz", 170.0, 8.0, 4000),
        ("steep2", "rampart_steep.xyz", 170.0, 2.0, 1))


def exact_step(w, dz0):
    """The largest of |h(x + dx) - h(x)| / (dz0 (1 - max(h) / z_top)) over the node pairs (BTF first layer)."""
    best = 0.0
    for a, b in zip(mt.XS[:-1], mt.XS[1:]):
        ha, hb = mt.cliff(a, w), mt.cliff(b, w)
        best = max(best, abs(hb - ha) / (dz0 * (1.0 - max(ha, hb) / mt.TOP)))
    return best


def value(log, name):
    m = re.search(re.escape(name) + r"\s+= (\S+)", log)
    if not m:
        raise SystemExit("check_rampart: the log has no line '%s'" % name)
    return float(m.group(1))


def main(argv):
    if len(argv) < 4 or argv[1] != "--":
        print(__doc__)
        return 2
    cmd = argv[2:]
    bad = []

    def check(good, what):
        print("  %s: %s" % ("ok" if good else "FAIL", what))
        if not good:
            bad.append(what)

    mt.main_to(os.getcwd())
    its = {}
    for name, terrain, w, dz0, limit in RUNS:
        shutil.rmtree(name, ignore_errors=True)
        os.makedirs(name)
        full = cmd + ["terrain.file=" + os.path.join("..", terrain), "mesh.dz0=%g" % dz0,
                      "mesh.stretch_ratio=%.6f" % mt.ratio(dz0), "steady.max_iterations=%d" % limit]
        print("check_rampart: in %s/: %s" % (name, " ".join(full)))
        out = subprocess.run(full, cwd=name, check=False, capture_output=True, text=True)
        log = out.stdout
        if out.returncode != 0:
            sys.stdout.write(log[-4000:] + out.stderr[-4000:])
            raise SystemExit("check_rampart: the command failed with code %d in %s/" % (out.returncode, name))
        step, ref = value(log, "max ground step / first layer thickness"), exact_step(w, dz0)
        check(abs(step - ref) <= 1e-3 * ref, "%s: ground step %.2f first layers (closed form %.2f, 0.1 %%)" % (name, step, ref))
        aspect, ref = value(log, "max aspect ratio max(dx,dy) / thickness"), mt.DX / (dz0 * (1.0 - mt.H / mt.TOP))
        check(abs(aspect - ref) <= 1e-3 * ref, "%s: aspect ratio %.2f (closed form %.2f, 0.1 %%)" % (name, aspect, ref))
        warned = "WARNING: the ground steps by up to" in log
        check(warned == (name == "steep2"), "%s: %s" % (name, "the report warns" if name == "steep2" else "no warning"))
        m = re.search(r"gust: steady after (\d+) iterations", log)
        its[name] = int(m.group(1)) if m else None
    check(its["gentle2"] is not None and its["gentle2"] <= 400,
          "gentle2 is steady within 400 iterations (%s)" % its["gentle2"])
    check(its["steep8"] is not None, "steep8 is steady within 4000 iterations (%s)" % its["steep8"])
    if its["gentle2"] and its["steep8"]:
        check(its["steep8"] > 4 * its["gentle2"], "steep8 takes more than four times gentle2's iterations (%d and %d)"
              % (its["steep8"], its["gentle2"]))
    if bad:
        print("check_rampart: %d check(s) failed" % len(bad))
        return 1
    print("check_rampart: all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
