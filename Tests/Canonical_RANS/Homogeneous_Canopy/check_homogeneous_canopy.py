#!/usr/bin/env python3
"""A horizontally uniform canopy against the canopy-flow relations (RT-8.6.1, RT-8.6.5):

    check_homogeneous_canopy.py -- LAUNCHER... GUST [ARGS ...]

LAUNCHER... GUST is the application with its launcher, ending with the
application (the script replaces the number after "-n" for the 1- and 4-rank
runs). ARGS are the deck (inputs_homogeneous_canopy) and the flags. The script
writes the tables and runs C_d = 0.2, h = 20 m at LAI 2, 5 and 10 (lai2/, lai5/,
lai10/), and LAI 10 again on 1 and 4 ranks (lai10_np1/, lai10_np4/), and checks:

  - the column's momentum budget: surface stress + integral of f U^2 = F H;
  - deep in a dense canopy (LAI 5, 10; z = h/4) the drag balances the pressure
    gradient, U = sqrt(F / f), a closed form with nothing fitted;
  - at LAI 2 and 5 the upper canopy (U/U_h >= 0.35) is exponential (Inoue 1963),
    with the rate n = f h / (2 beta^2), beta = u*(h) / U_h, u*(h)^2 = F (H - h).
    LAI 10 is printed only: its e-folding depth is about 2.5 layers;
  - beta in [0.25, 0.40], around the values of about 0.3 in Poggi et al. (2008);
  - LAI 10 gives the same answer on 1, 2 and 4 ranks.
The exit code is the verdict.
"""
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

F, HC, CD, Z0, KAPPA, DZ0, RATIO = 1.25e-3, 20.0, 0.2, 0.01, 0.4, 0.25, 1.02
TOL_BUDGET, TOL_DEEP, TOL_N, TOL_FIT, BETA_RANGE, TOL_RANKS = 0.01, 0.02, 0.15, 0.05, (0.25, 0.40), 1e-8


def with_ranks(cmd, n):
    out = list(cmd)
    for i, a in enumerate(out):
        if a in ("-n", "-np"):
            out[i + 1] = str(n)
            return out
    raise SystemExit("check_homogeneous_canopy: no -n in the launcher, cannot run on %d ranks" % n)


def run(cmd, where, lai):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    xs = [-200.0 + 200.0 * i for i in range(8)]
    for name, v in (("h.xyz", HC), ("lai.xyz", lai)):
        with open(os.path.join(where, name), "w") as f:
            for y in xs:
                for x in xs:
                    f.write("%g %g %g\n" % (x, y, v))
    print("check_homogeneous_canopy: in %s/: %s" % (where, " ".join(cmd)))
    out = subprocess.run(cmd, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-3000:] + out.stderr[-3000:])
        raise SystemExit("check_homogeneous_canopy: exit code %d in %s/" % (out.returncode, where))
    text = out.stdout + out.stderr
    if not any(ln.startswith("gust: steady after ") for ln in text.splitlines()):
        raise SystemExit("check_homogeneous_canopy: not steady in %s/" % where)
    plt = [ln.split()[-1] for ln in text.splitlines() if ln.startswith("gust: wrote ") and "_plt" in ln][-1]
    _, a = ep.read_fields(os.path.join(where, plt), ["u", "z", "detJ"])
    nz = len(a["u"][0][0])
    return ([a["z"][0][0][k] for k in range(nz)], [a["u"][0][0][k] for k in range(nz)],
            [a["detJ"][0][0][k] for k in range(nz)], a["u"])


def at(z, u, zq):
    for k in range(len(z) - 1):
        if z[k] <= zq <= z[k + 1]:
            return u[k] + (zq - z[k]) / (z[k + 1] - z[k]) * (u[k + 1] - u[k])
    raise ValueError(zq)


def check(z, u, lai):
    """The verdicts of one LAI, and its printed lines."""
    f = CD * lai / HC
    edges = [0.0]
    for k in range(len(z)):
        edges.append(edges[-1] + DZ0 * RATIO ** k)
    H = edges[-1]   # the mesh's top, from the layers
    drag = 0.0
    for k, zc in enumerate(z):
        top = min(edges[k + 1], HC)
        if edges[k] < HC:
            drag += f * u[k] ** 2 * (top - edges[k])
    tau = (KAPPA * u[0] / math.log(z[0] / Z0)) ** 2
    bud = (tau + drag) / (F * H) - 1.0
    ok = abs(bud) <= TOL_BUDGET
    print("  LAI %g: budget (tau_s %.3e + drag %.3e) / (F H) - 1 = %+.2e (<= %g): %s"
          % (lai, tau, drag, bud, TOL_BUDGET, "pass" if ok else "FAIL"))
    if lai >= 5:
        deep, exact = at(z, u, HC / 4), math.sqrt(F / f)
        good = abs(deep / exact - 1.0) <= TOL_DEEP
        ok = ok and good
        print("  LAI %g: U(h/4) %.4f, sqrt(F/f) %.4f m/s, ratio - 1 = %+.2e (<= %g): %s"
              % (lai, deep, exact, deep / exact - 1.0, TOL_DEEP, "pass" if good else "FAIL"))
    uh = at(z, u, HC)
    beta = math.sqrt(F * (H - HC)) / uh
    pts = [(zc, math.log(u[k] / uh)) for k, zc in enumerate(z) if zc < HC and u[k] / uh >= 0.35]
    m = len(pts)
    if m < 2:
        print("  LAI %g: %d layers with U/U_h >= 0.35, no fit: %s" % (lai, m, "FAIL" if lai <= 5 else "printed only"))
        return ok and lai > 5
    sx, sy = sum(p[0] for p in pts), sum(p[1] for p in pts)
    sxx, sxy = sum(p[0] ** 2 for p in pts), sum(p[0] * p[1] for p in pts)
    slope = (m * sxy - sx * sy) / (m * sxx - sx * sx)
    rms = math.sqrt(sum((p[1] - (sy + slope * (m * p[0] - sx)) / m) ** 2 for p in pts) / m)
    nfit, npred = slope * HC, f * HC / (2.0 * beta ** 2)
    good = abs(nfit / npred - 1.0) <= TOL_N and rms <= TOL_FIT
    if lai <= 5:
        ok = ok and good
    print("  LAI %g: %d layers, n %.2f, f h/(2 beta^2) %.2f (within %g), fit rms %.3f (<= %g): %s"
          % (lai, m, nfit, npred, TOL_N, rms, TOL_FIT, ("pass" if good else "FAIL") if lai <= 5 else "printed only"))
    good = BETA_RANGE[0] <= beta <= BETA_RANGE[1]
    ok = ok and good
    print("  LAI %g: beta %.3f in [%g, %g]: %s" % (lai, beta, BETA_RANGE[0], BETA_RANGE[1], "pass" if good else "FAIL"))
    return ok


def main(argv):
    if len(argv) < 3 or argv[1] != "--":
        raise SystemExit(__doc__)
    cmd = argv[2:]
    ok = True
    res = {}
    for lai in (2.0, 5.0, 10.0):
        res[lai] = run(cmd, "lai%g" % lai, lai)
        ok = check(res[lai][0], res[lai][1], lai) and ok
    ref = res[10.0][3]
    for n in (1, 4):
        other = run(with_ranks(cmd, n), "lai10_np%d" % n, 10.0)[3]
        worst = max(abs(ref[i][j][k] - other[i][j][k]) for i in range(len(ref))
                    for j in range(len(ref[0])) for k in range(len(ref[0][0])))
        good = worst <= TOL_RANKS
        ok = ok and good
        print("  LAI 10 on %d ranks against 2: largest u difference %.2e m/s (<= %g): %s"
              % (n, worst, TOL_RANKS, "pass" if good else "FAIL"))
    print("check_homogeneous_canopy: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
