#!/usr/bin/env python3
"""Run a parent or a nested patch to steady and check the patch against its parent.

    check_nest_faces.py parent PREFIX -- COMMAND [ARGS ...]
    check_nest_faces.py patch PARENT_PLOTFILE PREFIX OX OY TOL_VEL TOL_TKE NET_LO NET_HI -- COMMAND [ARGS ...]
    check_nest_faces.py open  PARENT_PLOTFILE PREFIX OX OY TOL_VEL TOL_TKE NET_LO NET_HI -- COMMAND [ARGS ...]

COMMAND is run with io.plot_file=PREFIX_plt io.chk_file=PREFIX_chk appended
(the decks do not name their outputs).

parent: COMMAND must reach the steady verdict; its last plotfile and
checkpoint are copied to
PREFIX_plt_final and PREFIX_chk_final, the names the patch decks use (the
step in the name depends on the iteration count).

patch: COMMAND is the patch at refinement 2 2 1
whose cell (0, 0) is the low corner of parent cell (OX, OY). Checked:

1. the log names the nest faces and the net flux through them, with
   NET_LO <= |net flux / flux entering| <= NET_HI, removed before the pressure
   solve;
2. the first pressure correction converges: the divergence after it is at
   most 1e-8 x the divergence before. With zero-gradient p' on every lateral
   face this needs a boundary flux with no net, so it fails if the net flux
   is left in;
3. the run reaches the steady verdict;
4. the patch keeps the parent's flow: the mean of u, v and w over the 2 x 2
   patch cells of every parent cell is within TOL_VEL (m/s) of the parent
   cell's value (PARENT_PLOTFILE), and that of k within TOL_TKE (m2/s2). On
   flat ground the tolerance is the steady tolerance of the two runs; over a
   hill it bounds what the finer cells change.

open: a patch with an outflow face beside its nest faces. Check 1 requires the
log to report the net flux as not removed (the outflow face carries it); the
other checks are those of patch.

The patch's last plotfile is copied to PREFIX_plt_final (for fcompare).
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
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402


def run(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1800)
    print(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(f"check_nest_faces: the run exited {proc.returncode}")
    return proc.stdout


def named(prefix, cmd):
    """cmd with its outputs named; the outputs of an earlier run removed (the last one is found by name)."""
    for d in glob.glob(prefix + "_plt[0-9]*") + glob.glob(prefix + "_chk[0-9]*"):
        shutil.rmtree(d)
    return cmd + [f"io.plot_file={prefix}_plt", f"io.chk_file={prefix}_chk"]


def keep_last(prefix, kind):
    dirs = sorted(d for d in glob.glob(f"{prefix}_{kind}[0-9]*") if os.path.isdir(d))
    if not dirs:
        raise SystemExit(f"check_nest_faces: no {prefix}_{kind}NNNNN was written")
    final = f"{prefix}_{kind}_final"
    shutil.rmtree(final, ignore_errors=True)
    shutil.copytree(dirs[-1], final)
    return final


def need_steady(log, bad):
    if not re.search(r"^gust: steady after \d+ iterations", log, re.M):
        bad.append("the run did not reach the steady verdict")


def parent(prefix, cmd):
    bad = []
    need_steady(run(named(prefix, cmd)), bad)
    if not bad:
        keep_last(prefix, "plt")
        keep_last(prefix, "chk")
    return bad


def patch(parent_plt, prefix, ox, oy, tol_vel, tol_tke, net_lo, net_hi, cmd, removed=True):
    log = run(named(prefix, cmd))
    bad = []
    if not re.search(r"^gust: nest faces: .* take the parent's u, v, w, theta and k", log, re.M):
        bad.append("the log does not name the nest faces")
    verdict = "removed before the pressure solve" if removed else r"not removed \(an outflow face carries it\)"
    m = re.search(r"^gust: nest faces: net outward flux (\S+) of the flux entering; " + verdict, log, re.M)
    if not m:
        bad.append("the log does not report the net flux as " + verdict.replace("\\", ""))
    else:
        net = abs(float(m.group(1)))
        print(f"check_nest_faces: |net flux| / flux entering = {net:.3e} (allowed {net_lo:g} to {net_hi:g})")
        if not net_lo <= net <= net_hi:
            bad.append(f"|net flux| / flux entering = {net:.3e} is outside [{net_lo:g}, {net_hi:g}]")
    m = re.search(r"^gust: initializer: .* divergence (\S+) before, (\S+) after", log, re.M)
    if not m:
        bad.append("the log has no initializer line")
    else:
        before, after = float(m.group(1)), float(m.group(2))
        print(f"check_nest_faces: divergence {before:.3e} before, {after:.3e} after")
        if not after <= 1.0e-8 * before:
            bad.append(f"the first pressure correction left a divergence of {after:.3e} ({before:.3e} before)")
    need_steady(log, bad)
    if bad:
        return bad

    final = keep_last(prefix, "plt")
    fields = ["u", "v", "w", "tke"]
    hp, P = erf_plotfile.read_fields(parent_plt, fields)
    hq, Q = erf_plotfile.read_fields(final, fields)
    nx, ny, nz = [hq["hi"][d] - hq["lo"][d] + 1 for d in range(3)]
    pn = [hp["hi"][d] - hp["lo"][d] + 1 for d in range(3)]
    if nx % 2 or ny % 2 or nz != pn[2] or ox + nx // 2 > pn[0] or oy + ny // 2 > pn[1]:
        return [f"the patch ({nx} x {ny} x {nz}) does not fit the parent {pn} at ({ox}, {oy})"]
    for f in fields:
        worst, where = 0.0, None
        for i in range(nx // 2):
            for j in range(ny // 2):
                for k in range(nz):
                    mean = 0.25 * (Q[f][2 * i][2 * j][k] + Q[f][2 * i + 1][2 * j][k]
                                   + Q[f][2 * i][2 * j + 1][k] + Q[f][2 * i + 1][2 * j + 1][k])
                    d = abs(mean - P[f][i + ox][j + oy][k])
                    if not d <= worst:   # a NaN is kept
                        worst, where = d, (i + ox, j + oy, k)
        tol = tol_tke if f == "tke" else tol_vel
        print(f"check_nest_faces: {f}: max |2 x 2 mean - parent| = {worst:.3e} at parent cell {where} (allowed {tol:g})")
        if not worst <= tol:
            bad.append(f"{f}: the patch differs from its parent by {worst:.3e} at parent cell {where} (allowed {tol:g})")
    return bad


def main(argv):
    if "--" not in argv or not argv:
        print(__doc__)
        return 2
    sep = argv.index("--")
    args, cmd = argv[:sep], argv[sep + 1:]
    if not cmd:
        print(__doc__)
        return 2
    if args[0] == "parent" and len(args) == 2:
        bad = parent(args[1], cmd)
    elif args[0] in ("patch", "open") and len(args) == 9:
        bad = patch(args[1], args[2], int(args[3]), int(args[4]), *[float(a) for a in args[5:9]], cmd,
                    removed=args[0] == "patch")
    else:
        print(__doc__)
        return 2
    for b in bad:
        print("check_nest_faces: FAIL: " + b)
    if not bad:
        print("check_nest_faces: PASS")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
