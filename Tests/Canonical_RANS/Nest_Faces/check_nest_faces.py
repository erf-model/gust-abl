#!/usr/bin/env python3
"""Run a parent or a nested patch to steady and check the patch against its parent.

    check_nest_faces.py parent PREFIX -- COMMAND [ARGS ...]
    check_nest_faces.py patch PARENT_PLOTFILE PREFIX OX OY TOL_VEL TOL_TKE NET_LO NET_HI -- COMMAND [ARGS ...]
    check_nest_faces.py open  PARENT_PLOTFILE PREFIX OX OY TOL_VEL TOL_TKE NET_LO NET_HI -- COMMAND [ARGS ...]
    check_nest_faces.py same FIELD PLOTFILE_A PLOTFILE_B
    check_nest_faces.py frozen PARENT_PLOTFILE PREFIX OX OY TOL_THETA MIN_RANGE TOL_VEL TOL_TKE -- COMMAND [ARGS ...]

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

frozen: a patch with simplec.theta = frozen, which holds the parent run's
theta. COMMAND (which gives neither steady.max_iterations nor
simplec.theta) is run twice with simplec.theta=frozen appended: for one
iteration (PREFIX_one) and to steady (PREFIX). Checked:

1. the log says theta is frozen at the parent run's theta;
2. the long run reaches the steady verdict, and theta's residual is 0 in
   every row of both logs (theta is data);
3. theta of the steady run equals theta of the one-iteration run bit for bit
   (it does not change over the iterations);
4. it is the parent's theta: every patch cell is within TOL_THETA (K) of
   the parent's theta interpolated bilinearly on its level (weights 3/4 and
   1/4 at ratio 2). A shift by one level or one parent cell fails;
5. the test is not vacuous: on every level theta varies by more than
   MIN_RANGE (K) across the patch. On the hill case the held field's least
   range is 0.022 K (top level); a theta set from the sounding by height
   above the ground has 0.0015 K on the first level and differs from the
   held field by up to 0.28 K, so it fails checks 4 and 5;
6. check 4 of patch (u, v, w, k against the parent, TOL_VEL and TOL_TKE).

same: FIELD of two plotfiles on the same cells is equal bit for bit (no run).

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

    return two_by_two(parent_plt, keep_last(prefix, "plt"), ox, oy, tol_vel, tol_tke)


def two_by_two(parent_plt, final, ox, oy, tol_vel, tol_tke):
    """The mean of u, v, w and k over the 2 x 2 patch cells of every parent cell against the parent cell's value."""
    bad = []
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


def frozen(parent_plt, prefix, ox, oy, tol_theta, min_range, tol_vel, tol_tke, cmd):
    bad = []
    logs = {}
    for name, count in ((prefix + "_one", 1), (prefix, 2000)):
        logs[name] = run(named(name, cmd) + ["simplec.theta=frozen", f"steady.max_iterations={count}"])
        if not re.search(r"^gust: theta frozen at the parent run's theta \(bc\.nest_parent\), interpolated: min \S+, "
                         r"max \S+, mean \S+ K$", logs[name], re.M):
            bad.append(f"{name}: the log does not say theta is frozen at the parent run's theta")
        # The iteration table: Iter CFL solver u v w tke theta p r.
        rows = re.findall(r"^\s*\d+\s+\S+\s+[a-z_]+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s+\S+\s*$", logs[name], re.M)
        if not rows:
            bad.append(f"{name}: the log has no iteration rows")
        elif any(float(r) != 0.0 for r in rows):
            bad.append(f"{name}: theta's residual is not 0 in every iteration (theta was solved)")
    need_steady(logs[prefix], bad)
    if bad:
        return bad
    final = keep_last(prefix, "plt")
    one = keep_last(prefix + "_one", "plt")
    hq, Q = erf_plotfile.read_fields(final, ["theta"])
    _, O = erf_plotfile.read_fields(one, ["theta"])
    hp, P = erf_plotfile.read_fields(parent_plt, ["theta"])
    q, o, pt = Q["theta"], O["theta"], P["theta"]
    nx, ny, nz = [hq["hi"][d] - hq["lo"][d] + 1 for d in range(3)]
    pn = [hp["hi"][d] - hp["lo"][d] + 1 for d in range(3)]
    if nx % 2 or ny % 2 or nz != pn[2] or ox + nx // 2 > pn[0] or oy + ny // 2 > pn[1]:
        return [f"the patch ({nx} x {ny} x {nz}) does not fit the parent {pn} at ({ox}, {oy})"]
    moved = sum(1 for i in range(nx) for j in range(ny) for k in range(nz) if not q[i][j][k] == o[i][j][k])
    print(f"check_nest_faces: theta differs from the one-iteration run's in {moved} of {nx * ny * nz} cells")
    if moved:
        bad.append(f"theta changed over the iterations in {moved} cells")
    # The parent's theta, bilinear on each level at the patch's cell centres (ratio 2: weights 3/4 and 1/4).
    def nb(n, o):
        return (o + n // 2, o + n // 2 - 1) if n % 2 == 0 else (o + n // 2, o + n // 2 + 1)
    worst, where = 0.0, None
    for i in range(nx):
        i0, i1 = nb(i, ox)
        for j in range(ny):
            j0, j1 = nb(j, oy)
            for k in range(nz):
                b = 0.75 * (0.75 * pt[i0][j0][k] + 0.25 * pt[i1][j0][k]) + 0.25 * (0.75 * pt[i0][j1][k] + 0.25 * pt[i1][j1][k])
                d = abs(q[i][j][k] - b)
                if not d <= worst:   # a NaN is kept
                    worst, where = d, (i, j, k)
    print(f"check_nest_faces: theta: max |patch - bilinear(parent)| = {worst:.3e} K at patch cell {where} (allowed {tol_theta:g})")
    if not worst <= tol_theta:
        bad.append(f"theta: the patch differs from the parent's interpolated theta by {worst:.3e} K at patch cell {where} "
                   f"(allowed {tol_theta:g})")
    least = min(max(q[i][j][k] for i in range(nx) for j in range(ny)) - min(q[i][j][k] for i in range(nx) for j in range(ny))
                for k in range(nz))
    print(f"check_nest_faces: theta: least range on a level = {least:.3e} K (needed above {min_range:g})")
    if not least > min_range:
        bad.append(f"theta varies by only {least:.3e} K on a level (needed above {min_range:g}): the test is vacuous")
    return bad + two_by_two(parent_plt, final, ox, oy, tol_vel, tol_tke)


def same(field, a, b):
    ha, A = erf_plotfile.read_fields(a, [field])
    hb, B = erf_plotfile.read_fields(b, [field])
    if ha["lo"] != hb["lo"] or ha["hi"] != hb["hi"]:
        return [f"{a} and {b} do not cover the same cells"]
    n = [ha["hi"][d] - ha["lo"][d] + 1 for d in range(3)]
    diff = sum(1 for i in range(n[0]) for j in range(n[1]) for k in range(n[2]) if not A[field][i][j][k] == B[field][i][j][k])
    print(f"check_nest_faces: {field} differs in {diff} of {n[0] * n[1] * n[2]} cells")
    return [f"{field} of {a} and {b} differs in {diff} cells"] if diff else []


def main(argv):
    if len(argv) == 4 and argv[0] == "same":
        bad = same(argv[1], argv[2], argv[3])
        for b in bad:
            print("check_nest_faces: FAIL:", b)
        print("check_nest_faces:", "FAIL" if bad else "PASS")
        return 1 if bad else 0
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
    elif args[0] == "frozen" and len(args) == 9:
        bad = frozen(args[1], args[2], int(args[3]), int(args[4]), float(args[5]), float(args[6]), float(args[7]), float(args[8]), cmd)
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
