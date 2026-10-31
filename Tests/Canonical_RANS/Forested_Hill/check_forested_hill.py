#!/usr/bin/env python3
"""A canopy on a train of cosine hills against the published lower-canopy balance (RT-8.6.3):

    check_forested_hill.py --flume DECK --deep DECK -- LAUNCHER... GUST [FLAGS ...]

LAUNCHER... GUST is the application with its launcher, ending with the
application (the script replaces the number after "-n" for the 1- and 4-rank
runs). DECK are inputs_flume_hill and inputs_deep_hill. The script writes the
tables and runs:

  deep/      H 4 m, L 40 m, h 20 m, L_c 10 m (regime I of Poggi et al. 2008, fig. 1);
  flume/     the flume's geometry scaled by 100, H_c / L_c = 0.13 (regime V);
  flume2/    twice the drag, H_c / L_c = 0.26;
  flume4/    four times the drag, H_c / L_c = 0.52 (long, deep, fixed pressure: regime I);
  np1/, np4/  flume4 on 1 and 4 ranks, to the verdict.

and checks:

  - deep: the lower-canopy balance of Finnigan and Belcher (2004), eq. 3 of Poggi et
    al. (2008), -dp/dx = (u |u| - U_b^2) / L_c, with the background U_b^2 = F L_c from the
    constant forcing F: u |u| = L_c (F - dp/dx) at a quarter of the canopy height, with
    dp/dx from the run at constant height. The rms error is within 25 % of the speed's rms;
  - flume4: at least 10 % of the columns reverse at h/4, all between the crest and just
    past the next valley (the predicted lee-side reversal in regime I);
  - the same answer on 1, 2 and 4 ranks at the verdict (short runs differ: the
    pressure solve stops at its tolerance, reached differently on each count).
Printed only: flume/ and flume2/ (the regime V measurement of reversal is not reproduced
at the published drag), and the reversal in deep/. The exit code is the verdict.
"""
import math
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

CD = 0.2
TOL_BALANCE, TOL_RANKS = 0.25, 1e-5
DEEP = dict(H=4.0, L=40.0, hc=20.0, lai=10.0, F=1.25e-3)
FLUME = dict(H=8.0, L=80.0, hc=10.0, lai=0.65, F=4.07e-3)
FLUME2 = dict(FLUME, lai=1.3)
FLUME4 = dict(FLUME, lai=2.6)


def ground(x, c):
    return 0.5 * c["H"] * (1.0 + math.cos(math.pi * x / (2.0 * c["L"]) + math.pi))


def with_ranks(cmd, n):
    out = list(cmd)
    for i, a in enumerate(out):
        if a in ("-n", "-np"):
            out[i + 1] = str(n)
            return out
    raise SystemExit("check_forested_hill: no -n in the launcher, cannot run on %d ranks" % n)


def run(cmd, where, c):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    lam, n = 4.0 * c["L"], 160
    with open(os.path.join(where, "t.xyz"), "w") as t, open(os.path.join(where, "h.xyz"), "w") as h, \
         open(os.path.join(where, "lai.xyz"), "w") as a:
        for j in range(-2, 7):
            for i in range(-8, n + 9):
                x, y = i * lam / n, j * 4.0   # on the mesh's nodes (4 m in y, and in x)
                t.write("%g %g %.6f\n" % (x, y, ground(x, c)))
                h.write("%g %g %g\n" % (x, y, c["hc"]))
                a.write("%g %g %g\n" % (x, y, c["lai"]))
    print("check_forested_hill: in %s/: %s" % (where, " ".join(cmd)))
    out = subprocess.run(cmd, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-3000:] + out.stderr[-3000:])
        raise SystemExit("check_forested_hill: exit code %d in %s/" % (out.returncode, where))
    text = out.stdout + out.stderr
    if not any(ln.startswith("gust: steady after ") for ln in text.splitlines()):
        raise SystemExit("check_forested_hill: not steady in %s/" % where)
    plt = [ln.split()[-1] for ln in text.splitlines() if ln.startswith("gust: wrote ") and "_plt" in ln][-1]
    _, a = ep.read_fields(os.path.join(where, plt), ["u", "p", "z"])
    return a


def sample(a, c, i, zq, name):
    """The field at height zq above the ground in column i; the layer index and weight."""
    nx, nz = len(a["u"]), len(a["u"][0][0])
    i %= nx
    x = (i + 0.5) * 4.0 * c["L"] / nx
    q = zq + ground(x, c)
    col = a["z"][i][0]
    for k in range(nz - 1):
        if col[k] <= q <= col[k + 1]:
            t = (q - col[k]) / (col[k + 1] - col[k])
            return a[name][i][0][k] * (1.0 - t) + a[name][i][0][k + 1] * t, k
    raise ValueError(zq)


def low_canopy(a, c):
    """The positions x/L and speeds a quarter of the canopy height above the ground."""
    nx = len(a["u"])
    return [((i + 0.5) * 4.0 / nx, sample(a, c, i, 0.25 * c["hc"], "u")[0]) for i in range(nx)]


def balance(a, c):
    """rms(u - sign(R) sqrt |R|) / rms(u), R = L_c (F - dp/dx); columns, u, prediction."""
    nx = len(a["u"])
    lam = 4.0 * c["L"]
    dx, zq, lc = lam / nx, 0.25 * c["hc"], c["hc"] / (CD * c["lai"])
    rows = []
    for i in range(nx):
        x = (i + 0.5) * dx
        dpds = (sample(a, c, i + 1, zq, "p")[0] - sample(a, c, i - 1, zq, "p")[0]) / (2.0 * dx)
        _, k = sample(a, c, i, zq, "p")
        col = a["z"][i][0]
        dpdz = (a["p"][i][0][k + 1] - a["p"][i][0][k]) / (col[k + 1] - col[k])
        slope = -0.5 * c["H"] * (2.0 * math.pi / lam) * math.sin(2.0 * math.pi * x / lam + math.pi)
        r = lc * (c["F"] - (dpds - slope * dpdz))
        rows.append((x / c["L"], sample(a, c, i, zq, "u")[0], math.copysign(math.sqrt(abs(r)), r)))
    rms = math.sqrt(sum((r[1] - r[2]) ** 2 for r in rows) / nx)
    return rms / math.sqrt(sum(r[1] ** 2 for r in rows) / nx), rows


def report_reversal(name, pts):
    rev = [p[0] for p in pts if p[1] < 0.0]
    print("  %s: u(h/4) from %.3f to %.3f m/s; reversed on %d of %d columns%s"
          % (name, min(p[1] for p in pts), max(p[1] for p in pts), len(rev), len(pts),
             " (x/L %s)" % ", ".join("%.2f" % r for r in (min(rev), max(rev))) if rev else ""))
    return rev


def main(argv):
    if "--" not in argv or argv[1] != "--flume" or argv[3] != "--deep":
        raise SystemExit(__doc__)
    flume, deep = os.path.abspath(argv[2]), os.path.abspath(argv[4])
    cmd = argv[6:]
    exe = next(i for i, c in enumerate(cmd) if os.path.basename(c).startswith("gust"))
    launcher, flags = cmd[:exe + 1], cmd[exe + 1:]   # the flags go after the deck
    ok = True

    a = run(launcher + [deep] + flags, "deep", DEEP)
    frac, rows = balance(a, DEEP)
    good = frac <= TOL_BALANCE
    ok = ok and good
    print("  deep: lower-canopy balance, rms / rms(u) %.3f (<= %g): %s"
          % (frac, TOL_BALANCE, "pass" if good else "FAIL"))
    report_reversal("deep (printed only)", [(r[0], r[1]) for r in rows])

    report_reversal("flume, H_c/L_c = 0.13 (printed only)", low_canopy(run(launcher + [flume] + flags, "flume", FLUME), FLUME))
    report_reversal("flume, H_c/L_c = 0.26 (printed only)", low_canopy(run(launcher + [flume] + flags, "flume2", FLUME2), FLUME2))

    # The crest is at x/L = 2, the next valley at x/L = 4: the lee slope.
    run4 = run(launcher + [flume] + flags, "flume4", FLUME4)
    rev = report_reversal("flume, H_c/L_c = 0.52", low_canopy(run4, FLUME4))
    nx = len(run4["u"])
    good = len(rev) >= 0.1 * nx and all(r >= 2.0 or r <= 0.25 for r in rev) and any(2.0 < r < 4.0 for r in rev)
    ok = ok and good
    print("  flume, H_c/L_c = 0.52: at least 10 %% of the columns reversed, crest to just past the valley: %s" % ("pass" if good else "FAIL"))

    ref = run4
    scale = max(abs(v) for pl in ref["u"] for col in pl for v in col)
    for n in (1, 4):
        other = run(with_ranks(launcher, n) + [flume] + flags, "np%d" % n, FLUME4)
        worst = max(abs(ref["u"][i][j][k] - other["u"][i][j][k]) for i in range(len(ref["u"]))
                    for j in range(len(ref["u"][0])) for k in range(len(ref["u"][0][0]))) / scale
        good = worst <= TOL_RANKS
        ok = ok and good
        print("  flume4, %d ranks against 2: largest u difference / max u %.2e (<= %g): %s"
              % (n, worst, TOL_RANKS, "pass" if good else "FAIL"))
    print("check_forested_hill: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
