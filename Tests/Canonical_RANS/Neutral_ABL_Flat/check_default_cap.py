#!/usr/bin/env python3
"""The default of turb.max_geom_lscale (100 m) through the application, and the earlier 30 m.

    check_default_cap.py [--iterations N] -- COMMAND [ARGS ...]

COMMAND is the application on Tests/Canonical_RANS/Neutral_ABL_Flat/inputs_steady_neutral_flat
(with its launcher and flags). The script writes a copy of the deck without its
turb.max_geom_lscale and steady.max_iterations lines and runs N iterations (default 40)
three times, each run in its own directory under the current one (emptied first):

  default/   no turb.max_geom_lscale: the default;
  cap100/    turb.max_geom_lscale = 100;
  cap30/     turb.max_geom_lscale = 30, the default before November 2026.

and checks:
  - default/ and cap100/ are identical in every field (the default is 100 m);
  - cap30/ differs from default/ in tke by more than 1 % of the largest tke
    (the earlier answers are reached by the input, and the cap acts).
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
DROPPED = ("turb.max_geom_lscale", "steady.max_iterations")


def run(cmd, deck_text, where, extra):
    shutil.rmtree(where, ignore_errors=True)
    os.makedirs(where)
    with open(os.path.join(where, "inputs"), "w") as fh:
        fh.write(deck_text)
    full = cmd + extra
    print("check_default_cap: in %s/: %s" % (where, " ".join(full)))
    out = subprocess.run(full, cwd=where, check=False, capture_output=True, text=True)
    if out.returncode != 0:
        sys.stdout.write(out.stdout[-4000:] + out.stderr[-4000:])
        raise SystemExit("check_default_cap: the command failed with code %d in %s/" % (out.returncode, where))
    wrote = [ln.split()[-1] for ln in out.stdout.splitlines() if ln.startswith("gust: wrote ") and "_plt" in ln]
    if not wrote:
        raise SystemExit("check_default_cap: no plotfile written in %s/" % where)
    return ep.read_fields(os.path.join(where, wrote[-1]), FIELDS)[1]


def max_diff(a, b, f):
    return max(abs(x - y) for pa, pb in zip(a[f], b[f]) for ca, cb in zip(pa, pb) for x, y in zip(ca, cb))


def main(argv):
    if "--" not in argv:
        print(__doc__)
        return 2
    sep = argv.index("--")
    opts, cmd = argv[1:sep], argv[sep + 1:]
    n = 40
    if opts:
        if len(opts) != 2 or opts[0] != "--iterations":
            print(__doc__)
            return 2
        n = int(opts[1])
    at = [i for i, a in enumerate(cmd) if "=" not in a and os.path.basename(a).startswith("inputs") and os.path.isfile(a)]
    if len(at) != 1:
        raise SystemExit("check_default_cap: the command must name one deck (inputs*)")
    with open(cmd[at[0]]) as fh:
        lines = fh.readlines()
    kept = [ln for ln in lines if ln.split("=")[0].strip() not in DROPPED]
    if len(lines) - len(kept) != len(DROPPED):
        raise SystemExit("check_default_cap: the deck must give each of %s once" % ", ".join(DROPPED))
    cmd = cmd[:at[0]] + ["inputs"] + cmd[at[0] + 1:] + ["steady.max_iterations=%d" % n]
    text = "".join(kept)
    default = run(cmd, text, "default", [])
    cap100 = run(cmd, text, "cap100", ["turb.max_geom_lscale=100."])
    cap30 = run(cmd, text, "cap30", ["turb.max_geom_lscale=30."])
    ok = True
    for f in FIELDS:
        d = max_diff(default, cap100, f)
        print("check_default_cap: %-5s default against 100 m: max difference %.3e" % (f, d))
        ok = ok and d == 0.0
    kmax = max(x for p in default["tke"] for c in p for x in c)
    d30 = max_diff(default, cap30, "tke")
    print("check_default_cap: tke   default against 30 m: max difference %.3e (largest tke %.3e)" % (d30, kmax))
    ok = ok and d30 > 0.01 * kmax
    print("check_default_cap: " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
