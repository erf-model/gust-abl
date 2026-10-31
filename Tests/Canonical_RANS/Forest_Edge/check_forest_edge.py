#!/usr/bin/env python3
"""A forest-clearing-forest pattern against the published adjustment length (RT-8.6.2):

    check_forest_edge.py -- LAUNCHER... GUST [ARGS ...]

LAUNCHER... GUST is the application with its launcher, ending with the
application (the script replaces the number after "-n" for the 1- and 4-rank
runs). ARGS are the deck (inputs_forest_edge) and the flags. The script writes
the tables (an 18 m forest of 306 m on the 6 m columns, centres 3 to 303 m, of the 690 m period) and runs
LAI 2 and 5 to the verdict (lai2/, lai5/), and LAI 2 again on 1 and 4 ranks to the
verdict (np1/, np4/), and checks:

  - the period's momentum budget: surface stress + drag = F H, integrated in x,
    from the cell-centred wind (the solver's drag acts on the faces);
  - the volume flux through every column is the same, again from the centres;
  - at LAI 2 the mean upward motion at the canopy top falls below w/u_ref = 0.01
    within [7.5, 10.5] h of the leading edge (Dupont and Brunet 2008: about 9 h),
    with u_ref = u(x = -8.5 h, z = 2 h);
  - the adjustment is shorter at LAI 5 (Belcher et al. 2003; Dupont and Brunet
    2008, sec. 3.2; their denser case 2, with another leaf profile, gives 8 h);
  - the same answer on 1, 2 and 4 ranks, four uneven boxes, at the verdict.
Printed only: the peak w/u_ref and x0 ~ L_c ln(U_h h / (u* L_c)) of Belcher et al.
(2003), eq. 3.3: the internal boundary layer's origin, a different length, an order of magnitude. The exit code is the verdict.
"""
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

F, HC, CD, Z0, KAPPA, DZ0, RATIO = 1.25e-3, 18.0, 0.2, 0.01, 0.4, 1.0, 1.05
LX, FOREST = 690.0, 308.4
W_END, ADJ_RANGE = 0.01, (7.5, 10.5)
TOL_BUDGET, TOL_FLUX, TOL_RANKS = 0.01, 1e-6, 1e-6


def with_ranks(cmd, n):
    out = list(cmd)
    for i, a in enumerate(out):
        if a in ("-n", "-np"):
            out[i + 1] = str(n)
            return out
    raise SystemExit("check_forest_edge: no -n in the launcher, cannot run on %d ranks" % n)


def run(cmd, where, lai):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    with open(os.path.join(where, "h.xyz"), "w") as h, open(os.path.join(where, "lai.xyz"), "w") as a:
        for y in (-9.0, -3.0, 3.0, 9.0, 15.0, 21.0, 27.0, 33.0):   # the columns' y, so every column is an exact hit
            for i in range(-8, 240):
                x = 3.0 * i
                forest = (x % LX) <= FOREST
                h.write("%g %g %g\n" % (x, y, HC if forest else 0.0))
                a.write("%g %g %g\n" % (x, y, lai if forest else 0.0))
    print("check_forest_edge: in %s/: %s" % (where, " ".join(cmd)))
    out = subprocess.run(cmd, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-3000:] + out.stderr[-3000:])
        raise SystemExit("check_forest_edge: exit code %d in %s/" % (out.returncode, where))
    text = out.stdout + out.stderr
    if not any(ln.startswith("gust: steady after ") for ln in text.splitlines()):
        raise SystemExit("check_forest_edge: not steady in %s/" % where)
    plt = [ln.split()[-1] for ln in text.splitlines() if ln.startswith("gust: wrote ") and "_plt" in ln][-1]
    _, a = ep.read_fields(os.path.join(where, plt), ["u", "w", "z"])
    return a


def at(zc, col, zq):
    for k in range(len(zc) - 1):
        if zc[k] <= zq <= zc[k + 1]:
            return col[k] + (zq - zc[k]) / (zc[k + 1] - zc[k]) * (col[k + 1] - col[k])
    raise ValueError(zq)


def check(a, lai):
    """The verdicts of one LAI and the adjustment length in canopy heights."""
    u, w = a["u"], a["w"]
    nx, nz = len(u), len(u[0][0])
    dx = LX / nx
    zc = [a["z"][0][0][k] for k in range(nz)]
    xs = [(i + 0.5) * dx for i in range(nx)]
    edges = [0.0]
    for k in range(nz):
        edges.append(edges[-1] + DZ0 * RATIO ** k)
    H = edges[-1]   # the mesh's top, from the layers
    f = CD * lai / HC
    force, fluxes = 0.0, []
    for i, x in enumerate(xs):
        col = u[i][0]
        force += (KAPPA * col[0] / math.log(zc[0] / Z0)) ** 2 * dx
        fluxes.append(sum(col[k] * (edges[k + 1] - edges[k]) for k in range(nz)))
        if x <= FOREST:
            for k in range(nz):
                if edges[k] < HC:
                    force += f * col[k] ** 2 * (min(edges[k + 1], HC) - edges[k]) * dx
    bud = force / (F * H * LX) - 1.0
    ok = abs(bud) <= TOL_BUDGET
    print("  LAI %g: budget / (F H L) - 1 = %+.2e (<= %g): %s" % (lai, bud, TOL_BUDGET, "pass" if ok else "FAIL"))
    mean = sum(fluxes) / nx
    spread = max(abs(q / mean - 1.0) for q in fluxes)
    good = spread <= TOL_FLUX
    ok = ok and good
    print("  LAI %g: volume flux spread %.2e (<= %g): %s" % (lai, spread, TOL_FLUX, "pass" if good else "FAIL"))
    iref = min(range(nx), key=lambda i: abs(xs[i] - (LX - 8.5 * HC)))
    uref = at(zc, u[iref][0], 2.0 * HC)
    wtop = [at(zc, w[i][0], HC) / uref for i in range(nx)]
    ipeak = max((i for i in range(nx) if xs[i] <= FOREST), key=lambda i: wtop[i])
    adj = None
    for i in range(ipeak, nx - 1):
        if wtop[i] >= W_END > wtop[i + 1]:
            adj = (xs[i] + dx * (wtop[i] - W_END) / (wtop[i] - wtop[i + 1])) / HC
            break
    if adj is None:
        print("  LAI %g: w/u_ref does not fall below %g: FAIL" % (lai, W_END))
        return False, None
    lc = 1.0 / f
    uh = at(zc, u[int(10.0 * HC / dx)][0], HC)
    ust = math.sqrt(F * (H - HC))
    print("  LAI %g: u_ref %.3f m/s, peak w/u_ref %.3f at %.2f h, adjustment %.2f h; L_c %.1f m, "
          "x0 = L_c ln(U_h h / (u* L_c)) = %.2f h (another length, printed only)"
          % (lai, uref, wtop[ipeak], xs[ipeak] / HC, adj, lc, lc * math.log(uh * HC / (ust * lc)) / HC))
    return ok, adj


def main(argv):
    if len(argv) < 3 or argv[1] != "--":
        raise SystemExit(__doc__)
    cmd = argv[2:]
    ok, adj = True, {}
    run2 = run(cmd, "lai2", 2.0)
    good, adj[2.0] = check(run2, 2.0)
    ok = ok and good
    good, adj[5.0] = check(run(cmd, "lai5", 5.0), 5.0)
    ok = ok and good
    if adj[2.0] is not None:
        good = ADJ_RANGE[0] <= adj[2.0] <= ADJ_RANGE[1]
        ok = ok and good
        print("  LAI 2: adjustment %.2f h in [%g, %g] h (Dupont and Brunet: about 9 h): %s"
              % (adj[2.0], ADJ_RANGE[0], ADJ_RANGE[1], "pass" if good else "FAIL"))
    if adj[2.0] is not None and adj[5.0] is not None:
        good = adj[5.0] < adj[2.0]
        ok = ok and good
        print("  adjustment at LAI 5 (%.2f h) shorter than at LAI 2 (%.2f h): %s"
              % (adj[5.0], adj[2.0], "pass" if good else "FAIL"))
    ref = run2
    scale = max(abs(v) for pl in ref["u"] for col in pl for v in col)
    for n in (1, 4):
        other = run(with_ranks(cmd, n), "np%d" % n, 2.0)
        worst = max(abs(ref[q][i][j][k] - other[q][i][j][k]) for q in ("u", "w") for i in range(len(ref[q]))
                    for j in range(len(ref[q][0])) for k in range(len(ref[q][0][0]))) / scale
        good = worst <= TOL_RANKS
        ok = ok and good
        print("  LAI 2, %d ranks against 2: largest difference / max u %.2e (<= %g): %s"
              % (n, worst, TOL_RANKS, "pass" if good else "FAIL"))
    print("check_forest_edge: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
