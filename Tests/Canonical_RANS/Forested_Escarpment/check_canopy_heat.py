#!/usr/bin/env python3
"""The forested escarpment with the leaves' heat exchange (canopy.heat = fixed).

    check_canopy_heat.py --max-iterations N -- LAUNCHER... GUST [ARGS ...]

LAUNCHER... GUST is the application with its launcher (the script replaces
the number after "-n" or "-np" for the 4-rank run). ARGS are the flags; the
decks, next to this script, are put right after GUST:
inputs_forested_escarpment (neutral, 300 K) and inputs_frozen_theta (a stable
sounding). The script writes the tables with make_case.py in each run
directory (emptied first) and runs, each to the verdict:

  cool/     2 ranks, theta solved, leaves at 298 K (cooler than the 300 K air);
  cool_np4/ the same on 4 ranks;
  frozen/   2 ranks, simplec.theta = frozen, no leaf heat;
  skip/     2 ranks, frozen, leaves at 304 K, canopy.heat_frozen = skip (the default);
  blend/    2 ranks, frozen, leaves at 304 K, canopy.heat_frozen = blend;
  first/    2 ranks, cool, 40 iterations, then continued from its checkpoint to 80;
  straight/ 2 ranks, cool, 80 iterations.

and checks:
  - every run is steady within N iterations;
  - cool: theta stays between the leaves' 298 K and the air's 300 K (no other
    source; 1e-6 K: without the leaf heat theta drifts by 2e-7 K from 300 K
    in rounding), the canopy cools the air (somewhere in the forest theta is
    below 299.5 K), and theta upstream of the forest (x < 2000 m) stays within
    1e-3 K of 300 K;
  - cool_np4 reaches the same answer: u, v, w within 1e-9 of the largest |u|,
    p and k within 1e-9 of their own largest values, theta within 1e-7 K (the
    runs stop at the verdict, r = 1e-4, where theta, the slowest field, still
    differs by a few 1e-9 K between the layouts' line sweeps);
  - skip logs the warning and gives the frozen run's fields bit for bit;
  - blend: outside the canopy (canopy_lad = 0) theta is the frozen run's bit
    for bit; inside, theta lies between
    the frozen theta and 304 K, and the largest blend weight
    (theta - theta_f) / (304 K - theta_f) is above 0.05;
  - the continued and straight plotfiles are identical in every field.
The exit code is the verdict.
"""
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

FIELDS = ["u", "v", "w", "p", "theta", "tke", "z", "canopy_lad"]
LEAF_COOL, LEAF_WARM, AIR = 298.0, 304.0, 300.0
WARNING = "gust: WARNING: canopy.heat = fixed with simplec.theta = frozen and canopy.heat_frozen = skip"


def with_ranks(cmd, n):
    out = list(cmd)
    for i, a in enumerate(out[:-1]):
        if a in ("-n", "-np"):
            out[i + 1] = str(n)
            return out
    if n != 2:
        raise SystemExit("check_canopy_heat: no -n in the launcher, cannot run on %d ranks" % n)
    return out


def run(cmd, where, extra):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    subprocess.run([sys.executable, os.path.join(HERE, "make_case.py"), where], check=True)
    full = cmd + extra
    print("check_canopy_heat: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_canopy_heat: the command failed with code %d in %s/" % (out.returncode, where))
    return out.stdout


def verdict(stdout):
    for ln in stdout.splitlines():
        if ln.startswith("gust: steady after "):
            return int(ln.split()[3])
    return None


def last_plotfile(stdout, where):
    wrote = [ln.split()[-1] for ln in stdout.splitlines()
             if ln.startswith("gust: wrote ") and "_plt" in ln and "_level" not in ln]
    if not wrote:
        raise SystemExit("check_canopy_heat: no plotfile written in %s/" % where)
    return os.path.join(where, wrote[-1])


def flat(a, f):
    return [x for plane in a[f] for col in plane for x in col]


def main(argv):
    nmax = None
    i = 1
    while i < len(argv) and argv[i] != "--":
        if argv[i] == "--max-iterations" and i + 1 < len(argv):
            nmax = int(argv[i + 1])
            i += 2
        else:
            raise SystemExit("unknown option " + argv[i])
    if nmax is None or i >= len(argv):
        raise SystemExit(__doc__)
    cmd = argv[i + 1:]
    neutral = [os.path.join(HERE, "inputs_forested_escarpment")]
    stable = [os.path.join(HERE, "inputs_frozen_theta"), "simplec.theta=frozen"]
    common = ["init.first_guess=profile", "steady.max_iterations=%d" % nmax]
    cool = ["canopy.heat=fixed", "canopy.leaf_theta=%g" % LEAF_COOL]
    warm = ["canopy.heat=fixed", "canopy.leaf_theta=%g" % LEAF_WARM]
    cwd = os.getcwd()
    runs = {
        "cool": (with_ranks(cmd, 2), neutral + common + cool),
        "cool_np4": (with_ranks(cmd, 4), neutral + common + cool),
        "frozen": (with_ranks(cmd, 2), stable + common),
        "skip": (with_ranks(cmd, 2), stable + common + warm),
        "blend": (with_ranks(cmd, 2), stable + common + warm + ["canopy.heat_frozen=blend"]),
    }
    out, data = {}, {}
    for name, (c, args) in runs.items():
        where = os.path.join(cwd, name)
        exe = next(n for n, x in enumerate(c) if os.path.basename(x) in ("gust", "gust.exe"))
        out[name] = run(c[:exe + 1] + args[:1] + c[exe + 1:] + args[1:], where, [])
        # canopy_lad is written with the leaf heat only (not in frozen/).
        fields = FIELDS if name != "frozen" else [f for f in FIELDS if f != "canopy_lad"]
        _, data[name] = ep.read_fields(last_plotfile(out[name], where), fields)
    ok = True

    def check(cond, msg):
        nonlocal ok
        print(("PASS  " if cond else "FAIL  ") + msg)
        ok = ok and cond

    for name in runs:
        n = verdict(out[name])
        check(n is not None and n <= nmax, "%s: steady after %s iterations (at most %d)" % (name, n, nmax))

    # cool: the maximum principle, the forest cools the air, upstream unchanged.
    a = data["cool"]
    th = flat(a, "theta")
    lo, hi = min(th), max(th)
    check(lo >= LEAF_COOL - 1e-6 and hi <= AIR + 1e-6,
          "cool: theta in [%.9f, %.9f] K, within [%g, %g] K" % (lo, hi, LEAF_COOL, AIR))
    nx = len(a["theta"])
    dx = 6000.0 / nx
    in_forest = [a["theta"][i][j][k] for i in range(nx) for j in range(len(a["theta"][i]))
                 for k in range(len(a["theta"][i][j])) if a["canopy_lad"][i][j][k] > 0]
    check(len(in_forest) > 0 and min(in_forest) < 299.5,
          "cool: coolest canopy cell %.4f K, below 299.5 K" % (min(in_forest) if in_forest else float("nan")))
    upstream = [a["theta"][i][j][k] for i in range(nx) if (i + 0.5) * dx < 2000.0
                for j in range(len(a["theta"][i])) for k in range(len(a["theta"][i][j]))]
    dev = max(abs(t - AIR) for t in upstream)
    check(dev <= 1e-6, "cool: theta upstream of the forest (x < 2000 m) within %.3e K of 300 K (1e-6)" % dev)

    # cool_np4 against cool.
    b = data["cool_np4"]
    umax = max(abs(x) for x in flat(a, "u"))
    worst = 0.0
    for f in ["u", "v", "w", "p", "tke"]:
        x, y = flat(a, f), flat(b, f)
        scale = umax if f in ("u", "v", "w") else max(abs(t) for t in x)
        worst = max(worst, max(abs(p - q) for p, q in zip(x, y)) / scale)
    check(worst <= 1e-9, "cool_np4: largest relative difference of u, v, w, p, k from 2 ranks %.2e (1e-9)" % worst)
    dth = max(abs(p - q) for p, q in zip(flat(a, "theta"), flat(b, "theta")))
    check(dth <= 1e-7, "cool_np4: largest theta difference from 2 ranks %.2e K (1e-7 K)" % dth)

    # skip: the warning, and the frozen run bit for bit.
    check(WARNING in out["skip"], "skip: the start-up warning is logged")
    same = all(flat(data["skip"], f) == flat(data["frozen"], f) for f in ["u", "v", "w", "p", "theta", "tke"])
    check(same, "skip: every field equals the frozen run's bit for bit")

    # blend: theta_f outside the canopy bit for bit; between theta_f and the leaves inside.
    t_b, t_f, drag = flat(data["blend"], "theta"), flat(data["frozen"], "theta"), flat(data["blend"], "canopy_lad")
    outside = [(p, q) for p, q, d in zip(t_b, t_f, drag) if d == 0]
    inside = [(p, q) for p, q, d in zip(t_b, t_f, drag) if d > 0]
    check(all(p == q for p, q in outside), "blend: theta outside the canopy is the frozen theta bit for bit (%d cells)"
          % len(outside))
    between = all(q <= p <= LEAF_WARM for p, q in inside)
    wmax = max((p - q) / (LEAF_WARM - q) for p, q in inside) if inside else 0.0
    check(len(inside) > 0 and between, "blend: theta in the canopy between the frozen theta and %g K (%d cells)"
          % (LEAF_WARM, len(inside)))
    check(wmax > 0.05, "blend: largest blend weight %.4f, above 0.05" % wmax)

    # Restart with the leaf heat: 40 iterations, continued from the checkpoint to 80, against 80 straight.
    c2 = with_ranks(cmd, 2)
    exe = next(n for n, x in enumerate(c2) if os.path.basename(x) in ("gust", "gust.exe"))
    base = c2[:exe + 1] + neutral + c2[exe + 1:] + cool
    first = run(base, os.path.join(cwd, "first"), ["init.first_guess=profile", "steady.max_iterations=40"])
    chk = [ln.split()[-1] for ln in first.splitlines() if ln.startswith("gust: wrote ") and "_chk" in ln]
    if not chk:
        raise SystemExit("check_canopy_heat: no checkpoint written in first/")
    cont = os.path.join(cwd, "first")
    out_c = subprocess.run(base + ["io.restart=" + chk[-1], "steady.max_iterations=80"], cwd=cont, check=False,
                           capture_output=True, text=True)
    if out_c.returncode != 0:
        sys.stdout.write(out_c.stdout[-4000:] + out_c.stderr[-4000:])
        raise SystemExit("check_canopy_heat: the restart failed in first/")
    straight = run(base, os.path.join(cwd, "straight"), ["init.first_guess=profile", "steady.max_iterations=80"])
    _, d_c = ep.read_fields(last_plotfile(out_c.stdout, cont), FIELDS)
    _, d_s = ep.read_fields(last_plotfile(straight, os.path.join(cwd, "straight")), FIELDS)
    same = all(flat(d_c, f) == flat(d_s, f) for f in ["u", "v", "w", "p", "theta", "tke"])
    check(same, "restart: 40 + 40 iterations with the leaf heat equal 80 straight, every field bit for bit")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
