#!/usr/bin/env python3
"""The canopy's drag on the momentum equations (canopy.drag_diagonal) on two
forest patches over flat ground.

    check_canopy_drag.py --tol T --max-ratio Q -- COMMAND [ARGS ...]

COMMAND is the application on Tests/Canonical_RANS/Canopy_Drag/inputs_canopy_drag
with canopy.file (and its launcher and flags). The script runs it, each run in
its own directory under the current one (emptied first):

  derivative/  canopy.drag_diagonal = derivative, to the steady verdict;
  picard/      canopy.drag_diagonal = picard, to the steady verdict;
  sequenced/   derivative with grid sequencing over two levels, to the verdict;
  none/        canopy.drag_diagonal = none (the drag explicit), 20 iterations;
  first/       derivative, 40 iterations; then continued from its checkpoint to 80;
  other/       continued from first/'s checkpoint to 80 with canopy.drag_diagonal = picard;
  straight/    derivative, 80 iterations.

and checks:
  - both implicit runs reach the verdict, and the log names their option;
  - u, v, w, p and k of the two steady states agree to T relative to the
    largest |value| of each field: the diagonal changes the pseudo-time path,
    not the fixed point (T covers the convergence error at the verdict);
  - with grid sequencing (the coarse level takes the restricted drag field)
    the steady state is the same to T, and the log reports the coarse
    level's canopy;
  - the drag acts: at the first cell above 5 m, u at the centre of patch 2
    (Lalic-Mihailovic, LAI 6) is at most Q times u at the open-ground column
    farthest from the patches;
  - with the drag explicit the iteration departs: after 20 iterations its
    residual norm r is above 1 and above the derivative run's (the drag rate
    f |U| exceeds 1 / dtau, the explicit limit along the wind, at the
    starting CFL);
  - the continued and straight plotfiles are identical in every field (the
    canopy and its face drag are rebuilt from the inputs on restart), and a
    restart with another canopy.drag_diagonal is accepted (it is not part of
    the checkpoint's canopy hash).
The exit code is the verdict.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

FIELDS = ["u", "v", "w", "p", "theta", "tke"]


def run(cmd, where, extra):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    full = cmd + extra
    print("check_canopy_drag: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_canopy_drag: the command failed with code %d in %s/" % (out.returncode, where))
    return out.stdout


def verdict(stdout):
    for ln in stdout.splitlines():
        if ln.startswith("gust: steady after "):
            return int(ln.split()[3])
    return None


def residual_at(stdout, it):
    """r (the last column of the iteration table) at iteration it."""
    for ln in stdout.splitlines():
        parts = ln.split()
        if len(parts) > 3 and parts[0] == str(it) and parts[2] == "lusgs":
            return float(parts[-1])
    raise SystemExit("check_canopy_drag: no table row for iteration %d" % it)


def last_plotfile(stdout, where):
    # Level 0's plotfile (with sequencing the coarse levels' are <plot_file>_level<l>_<step>).
    wrote = [ln.split()[-1] for ln in stdout.splitlines()
             if ln.startswith("gust: wrote ") and "_plt" in ln and "_level" not in ln]
    if not wrote:
        raise SystemExit("check_canopy_drag: no plotfile written in %s/" % where)
    return os.path.join(where, wrote[-1])


def values(a, f):
    return [x for plane in a[f] for col in plane for x in col]


def main(argv):
    tol, ratio = None, None
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] in ("--tol", "--max-ratio") and i + 1 < len(argv):
            if argv[i] == "--tol":
                tol = float(argv[i + 1])
            else:
                ratio = float(argv[i + 1])
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if tol is None or ratio is None or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    ok = True

    fresh = ["init.first_guess=profile"]
    der = run(cmd, "derivative", ["canopy.drag_diagonal=derivative"] + fresh)
    pic = run(cmd, "picard", ["canopy.drag_diagonal=picard"] + fresh)
    n_d, n_p = verdict(der), verdict(pic)
    good = n_d is not None and n_p is not None
    print("  steady verdict with derivative and picard: %s, %s iterations: %s" % (n_d, n_p, "pass" if good else "FAIL"))
    if not good:
        raise SystemExit("check_canopy_drag: FAIL")
    good = ("canopy.drag_diagonal = derivative" in der) and ("canopy.drag_diagonal = picard" in pic)
    ok = ok and good
    print("  the log names the option: %s" % ("pass" if good else "FAIL"))

    hdr, a = ep.read_fields(last_plotfile(der, "derivative"), FIELDS + ["z"])
    _, b = ep.read_fields(last_plotfile(pic, "picard"), FIELDS)
    for f in ("u", "v", "w", "p", "tke"):
        va, vb = values(a, f), values(b, f)
        big = max(abs(x) for x in va)
        d = max(abs(x - y) for x, y in zip(va, vb)) / big
        good = d <= tol
        ok = ok and good
        print("  max |%s(derivative) - %s(picard)| / max |%s| = %.3e (<= %g): %s"
              % (f, f, f, d, tol, "pass" if good else "FAIL"))

    seq = run(cmd, "sequenced", ["mesh.hierarchy_levels=2", "steady.sequencing=1"] + fresh)
    good = verdict(seq) is not None and "canopy (restricted from level 0)" in seq
    ok = ok and good
    print("  grid sequencing: steady, the coarse level's canopy logged: %s" % ("pass" if good else "FAIL"))
    _, q = ep.read_fields(last_plotfile(seq, "sequenced"), FIELDS)
    for f in ("u", "v", "w", "p", "tke"):
        va, vq = values(a, f), values(q, f)
        d = max(abs(x - y) for x, y in zip(va, vq)) / max(abs(x) for x in va)
        good = d <= tol
        ok = ok and good
        print("  max |%s(sequenced) - %s(derivative)| / max |%s| = %.3e (<= %g): %s"
              % (f, f, f, d, tol, "pass" if good else "FAIL"))

    # Patch 2 is centred at (960, 640): column (11, 7) or (12, 8) of the 80 m
    # columns; the open column farthest from both patches is (0, 0)'s
    # periodic neighbourhood: x = 640, y = 40 (column (7, 0)).
    k = next(kk for kk in range(len(a["z"][0][0])) if a["z"][11][7][kk] > 5.0)
    u_in, u_open = a["u"][11][7][k], a["u"][7][0][k]
    good = u_in <= ratio * u_open
    ok = ok and good
    print("  u at z = %.1f m: %.3f m/s in patch 2, %.3f m/s on open ground (ratio %.2f <= %g): %s"
          % (a["z"][11][7][k], u_in, u_open, u_in / u_open, ratio, "pass" if good else "FAIL"))

    none = run(cmd, "none", ["canopy.drag_diagonal=none", "steady.max_iterations=20"] + fresh)
    r_none, r_der = residual_at(none, 20), residual_at(der, 20)
    good = r_none > 1 and r_none > r_der
    ok = ok and good
    print("  r after 20 iterations: %.3e explicit, %.3e derivative (explicit above 1 and above): %s"
          % (r_none, r_der, "pass" if good else "FAIL"))

    first = run(cmd, "first", ["steady.max_iterations=40"] + fresh)
    chk = [ln.split()[-1] for ln in first.splitlines() if ln.startswith("gust: wrote ") and "_chk" in ln]
    if not chk:
        raise SystemExit("check_canopy_drag: no checkpoint written in first/")
    cont = subprocess.run(cmd + ["steady.max_iterations=80", "io.restart=" + chk[-1]], cwd="first",
                          check=False, capture_output=True, text=True)
    if cont.returncode != 0:
        sys.stdout.write(cont.stdout[-4000:] + cont.stderr[-4000:])
        raise SystemExit("check_canopy_drag: the restart failed with code %d" % cont.returncode)
    other = run(cmd, "other", ["steady.max_iterations=80", "io.restart=" + os.path.join("..", "first", chk[-1]),
                               "canopy.drag_diagonal=picard"])
    good = "canopy.drag_diagonal = picard" in other and "restarted from" in other
    ok = ok and good
    print("  a restart with canopy.drag_diagonal = picard is accepted: %s" % ("pass" if good else "FAIL"))
    straight = run(cmd, "straight", ["steady.max_iterations=80"] + fresh)
    _, c = ep.read_fields(last_plotfile(cont.stdout, "first"), FIELDS)
    _, s = ep.read_fields(last_plotfile(straight, "straight"), FIELDS)
    for f in FIELDS:
        d = max(abs(x - y) for x, y in zip(values(c, f), values(s, f)))
        good = d == 0
        ok = ok and good
        print("  restart: max |%s(continued) - %s(straight)| = %.3e (bit for bit): %s"
              % (f, f, d, "pass" if good else "FAIL"))
    if not ok:
        raise SystemExit("check_canopy_drag: FAIL")
    print("check_canopy_drag: pass")


if __name__ == "__main__":
    main(sys.argv)
