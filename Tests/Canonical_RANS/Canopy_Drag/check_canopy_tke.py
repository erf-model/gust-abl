#!/usr/bin/env python3
"""The canopy terms of the k-L closure (canopy.tke) on two forest patches over
flat ground.

    check_canopy_tke.py --nut-ratio Q -- COMMAND [ARGS ...]

COMMAND is the application on Tests/Canonical_RANS/Canopy_Drag/inputs_canopy_drag
with canopy.file (and its launcher and flags). The script runs it, each run in
its own directory under the current one (emptied first), with the turbulence
plotfile:

  on/   canopy.tke = on (the default: beta_p 1, beta_d 4, length_alpha 0.04), to the verdict;
  off/  canopy.tke = off (the drag only, the behaviour before the k terms), to the verdict.

and checks:
  - both reach the verdict, and the log says whether the k terms are on;
  - the length cap: in every cell with canopy (f > 0) of on/, L <= alpha / f
    (to rounding), and it holds with equality in some cells (the cap acts);
    off/ exceeds alpha / f in some canopy cell (so the check is not vacuous);
  - the k terms act: in the crown of the dense patch (the cell of patch 2's
    centre column with the largest f) the eddy viscosity of on/ is at most Q
    times that of off/.
The local balance of the canopy terms (k = kappa_c |U|^2) is pinned by the
unit test CanopyTke.CappedBalanceIsKappaU2; at the crown the shear layer
above adds production and transport, so it is not checked here.
The exit code is the verdict.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

ALPHA = 0.04


def run(cmd, where, extra):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    full = cmd + extra
    print("check_canopy_tke: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_canopy_tke: the command failed with code %d in %s/" % (out.returncode, where))
    return out.stdout


def verdict(stdout):
    for ln in stdout.splitlines():
        if ln.startswith("gust: steady after "):
            return int(ln.split()[3])
    return None


def last_written(stdout, where, tag):
    wrote = [ln.split()[-1] for ln in stdout.splitlines()
             if ln.startswith("gust: wrote ") and tag in ln and "_level" not in ln]
    if not wrote:
        raise SystemExit("check_canopy_tke: no %s written in %s/" % (tag, where))
    return os.path.join(where, wrote[-1])


def fields(stdout, where):
    _, t = ep.read_fields(last_written(stdout, where, "turb"),
                          ["x_velocity", "y_velocity", "z_velocity", "KE", "Kmv", "Lturb"])
    _, p = ep.read_fields(last_written(stdout, where, "_plt"), ["canopy_drag"])
    t["f"] = p["canopy_drag"]
    return t


def main(argv):
    ratio = None
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] == "--nut-ratio" and i + 1 < len(argv):
            ratio = float(argv[i + 1])
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if ratio is None or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    common = ["init.first_guess=profile", "io.turb_plot_file=turb"]
    ok = True

    on = run(cmd, "on", common)
    off = run(cmd, "off", common + ["canopy.tke=off"])
    n_on, n_off = verdict(on), verdict(off)
    good = n_on is not None and n_off is not None
    print("  steady verdict with the k terms on and off: %s, %s iterations: %s" % (n_on, n_off, "pass" if good else "FAIL"))
    if not good:
        raise SystemExit("check_canopy_tke: FAIL")
    good = "k terms on (canopy.tke)" in on and "k terms off (canopy.tke = off)" in off
    ok = ok and good
    print("  the log says whether the k terms are on: %s" % ("pass" if good else "FAIL"))

    a, b = fields(on, "on"), fields(off, "off")
    nx, ny, nz = len(a["f"]), len(a["f"][0]), len(a["f"][0][0])
    worst, at_cap, over_off = 0.0, 0, 0
    for ii in range(nx):
        for jj in range(ny):
            for kk in range(nz):
                f = a["f"][ii][jj][kk]
                if f <= 0:
                    continue
                cap = ALPHA / f
                worst = max(worst, (a["Lturb"][ii][jj][kk] - cap) / cap)
                if abs(a["Lturb"][ii][jj][kk] - cap) <= 1e-12 * cap:
                    at_cap += 1
                if b["Lturb"][ii][jj][kk] > cap * (1 + 1e-12):
                    over_off += 1
    good = worst <= 1e-12 and at_cap > 0 and over_off > 0
    ok = ok and good
    print("  length cap: max (L - alpha/f)/(alpha/f) = %.2e (<= 1e-12), %d cells at the cap, %d cells of off/ above it: %s"
          % (worst, at_cap, over_off, "pass" if good else "FAIL"))

    # The crown of patch 2 (centre (960, 640): column (11, 7) of the 80 m columns).
    col = [a["f"][11][7][kk] for kk in range(nz)]
    kc = max(range(nz), key=lambda kk: col[kk])
    nut_on, nut_off = a["Kmv"][11][7][kc], b["Kmv"][11][7][kc]
    good = nut_on <= ratio * nut_off
    ok = ok and good
    print("  crown eddy viscosity: %.4f with the k terms, %.4f without (ratio %.3f <= %g): %s"
          % (nut_on, nut_off, nut_on / nut_off, ratio, "pass" if good else "FAIL"))
    if not ok:
        raise SystemExit("check_canopy_tke: FAIL")
    print("check_canopy_tke: pass")


if __name__ == "__main__":
    main(sys.argv)
