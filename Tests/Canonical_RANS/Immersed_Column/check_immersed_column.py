#!/usr/bin/env python3
"""Run GUST on the immersed column and check the steady state (rt_immersed_column_*).

    check_immersed_column.py --f F --cut-max PCT [--body-fitted DECK] -- COMMAND [ARGS ...]

COMMAND runs the application on inputs_immersed_column with terrain.flat_height = (2 + F) dz, so
the cut cell is layer k = 2 with solid fraction F. The run must reach the steady verdict; then:
  * the force balance: u*^2 of the immersed wall law at the exact distances,
        u*^2 = (1 - b) [kappa u_c / ln((d1 + z0)/z0)]^2 + b [kappa u_a / ln((d2 + z0)/z0)]^2,
    b = F, d1 = (1 - b) dz / 2, d2 = (1 - b) dz + dz / 2, u_c the cut cell's speed and u_a the
    speed of the cell above it (b = 0: the first term alone), equals F_x D to 1e-5 relative,
    F_x = 1.225e-4 m/s^2 the force and D = top - flat_height the fluid depth. The force acts on
    the fluid only and the top is a slip wall, so at steady the wall stress balances it exactly;
    the residual tolerance leaves about 5e-7;
  * the solid layers (k < 2) have u = 0 exactly;
  * the cell above the cut cell is within 1 % of the log law with u* = sqrt(F_x D) at its centre's
    height above the ground, the next two within 2.5 % (the body-fitted column: 0.9 % and 1.3-1.5 %);
  * the cut cell's speed is at or above the log law at its fluid centroid d1 (to 1e-3 %) by at most PCT
    percent (the form-7 law reads the cut cell's speed as at its centroid; the measured excess
    grows with F, README.md);
  * every column is identical (|u - planar mean| below 1e-10 of u*), every
    field finite;
  * --body-fitted DECK (F = 0 only): the same COMMAND with DECK in place of the immersed deck, run
    in the sub-folder body_fitted/, is the body-fitted column of the fluid's height; its u and tke
    in every layer equal the immersed column's layers k >= 2 to 1e-6 relative.
Prints one row per check; the exit code is the verdict.
"""
import math
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402

DZ = 15.625
TOP = 67 * DZ
FORCE = 1.225e-4
Z0 = 0.1
KAPPA = 0.41
KC = 2   # the cut cell's layer


def run(cmd, cwd):
    os.makedirs(cwd, exist_ok=True)
    out = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    sys.stdout.write(out.stdout[-2000:])
    sys.stderr.write(out.stderr[-2000:])
    return out.returncode, out.stdout + out.stderr


def written_plotfile(cwd, log, prefix):
    """The plotfile this run wrote ("gust: wrote <name>"): an earlier run's may lie in the folder."""
    wrote = [ln.split()[-1] for ln in log.splitlines() if ln.startswith("gust: wrote " + prefix)]
    if not wrote:
        raise SystemExit("check_immersed_column: the run wrote no " + prefix + "*")
    return os.path.join(cwd, wrote[-1])


def main(argv):
    if "--" not in argv:
        raise SystemExit(__doc__)
    sep = argv.index("--")
    opts, cmd = argv[:sep], argv[sep + 1:]
    f = cut_max = None
    body_fitted = None
    i = 0
    while i < len(opts):
        if opts[i] == "--f":
            f = float(opts[i + 1])
        elif opts[i] == "--cut-max":
            cut_max = float(opts[i + 1])
        elif opts[i] == "--body-fitted":
            body_fitted = opts[i + 1]
        else:
            raise SystemExit("check_immersed_column: unknown option " + opts[i])
        i += 2
    if f is None or cut_max is None or not 0 <= f < 1:
        raise SystemExit("check_immersed_column: --f in [0, 1) and --cut-max are required")
    height = (KC + f) * DZ
    depth = TOP - height
    us = math.sqrt(FORCE * depth)
    cmd = cmd + ["terrain.flat_height=" + repr(height)]

    ok = True

    def check(name, good, detail):
        nonlocal ok
        ok = ok and good
        print(("PASS " if good else "FAIL ") + name + ": " + detail)

    code, log = run(cmd, os.getcwd())
    check("steady verdict", code == 0 and "gust: steady after" in log, "exit " + str(code))
    if not ok:
        return 1
    plt = written_plotfile(os.getcwd(), log, "immersed_column_plt")
    z, prof, _ = erf_plotfile.planar_averages(plt, ["u", "tke"])
    u = prof["u"]
    _, full = erf_plotfile.read_fields(plt, ["u"])
    spread = max(abs(col[k] - u[k]) for row in full["u"] for col in row for k in range(len(u)))
    check("layers", len(u) == 67 and abs(z[1] - z[0] - DZ) < 1e-9 * DZ, str(len(u)) + " layers")
    check("finite", all(math.isfinite(x) for x in list(u) + list(prof["tke"])), "u and tke")
    check("identical columns", spread < 1e-10 * us, "max |u - planar mean| %.3e" % spread)

    b = f
    d1 = (1 - b) * DZ / 2
    d2 = (1 - b) * DZ + DZ / 2
    us2 = (1 - b) * (KAPPA * u[KC] / math.log((d1 + Z0) / Z0)) ** 2
    if b > 0:
        us2 += b * (KAPPA * u[KC + 1] / math.log((d2 + Z0) / Z0)) ** 2
    err = us2 / (FORCE * depth) - 1
    check("force balance", abs(err) < 1e-5, "u*^2 / (F D) - 1 = %+.3e" % err)
    check("solid layers", all(u[k] == 0 for k in range(KC)), "u = " + ", ".join("%g" % u[k] for k in range(KC)))

    def log_law(d):
        return us / KAPPA * math.log((d + Z0) / Z0)

    cut = 100 * (u[KC] / log_law(d1) - 1)
    check("cut cell", cut >= -1e-3 and cut <= cut_max, "%+.2f %% of the log law at d1 = %.4f m (max %.2f %%)" % (cut, d1, cut_max))
    for n, k in enumerate(range(KC + 1, KC + 4)):
        d = (k + 0.5) * DZ - height
        e = 100 * (u[k] / log_law(d) - 1)
        lim = 1.0 if n == 0 else 2.5
        check("layer %d" % k, abs(e) <= lim, "%+.2f %% of the log law at %.3f m (max %.1f %%)" % (e, d, lim))

    if body_fitted is not None:
        if f != 0:
            raise SystemExit("check_immersed_column: --body-fitted needs --f 0")
        deck = [i for i, a in enumerate(cmd) if a.endswith("inputs_immersed_column")]
        if len(deck) != 1:
            raise SystemExit("check_immersed_column: COMMAND names inputs_immersed_column once")
        bf_cmd = [a for a in cmd if not a.startswith("terrain.flat_height=")]
        bf_cmd[deck[0]] = os.path.abspath(body_fitted)
        sub = os.path.join(os.getcwd(), "body_fitted")
        code, log = run(bf_cmd, sub)
        check("body-fitted steady verdict", code == 0 and "gust: steady after" in log, "exit " + str(code))
        if code == 0:
            _, bprof, _ = erf_plotfile.planar_averages(written_plotfile(sub, log, "body_fitted_column_plt"), ["u", "tke"])
            for name in ("u", "tke"):
                a, r = prof[name][KC:], bprof[name]
                diff = max(abs(x / y - 1) for x, y in zip(a, r)) if len(a) == len(r) else float("inf")
                check("body-fitted " + name, diff < 1e-6, "max relative difference %.3e over %d layers" % (diff, len(r)))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
