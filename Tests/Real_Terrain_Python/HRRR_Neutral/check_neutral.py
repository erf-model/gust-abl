#!/usr/bin/env python3
"""The test deck of the C3 options: the strong-cap hour apr05_13 (2023-04-05 13 UTC) with the recommended option,
run from the downloads of prepare.py and checked against the numbers of November 2026.

    python3 check_neutral.py run --gust EXE [--data DIR] [--work DIR] [--mpi "mpiexec -n 4"]
    python3 check_neutral.py check [SCORE_JSON]

"run" writes the case with options.py into --work (default: check_run/ beside this script), runs it with
run_options.sh and checks the score; "check" checks a score.json written earlier. Without the downloads
(awaken.asc and the crop of apr05_13 in --data, default: $GUST_HRRR_NEUTRAL_DATA, else beside this script) "run"
prints what is missing and exits with 77, which CTest reads as a skip. Exit status 0 when every number is inside
its tolerance (and the run is steady within 1.5 times the November iterations, with w max below 1.5 m/s), 1
otherwise.
"""
import argparse
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from prepare import crop_name  # noqa: E402

HOUR = "apr05_13"
OPTION = "no_coriolis"          # the recommended option for C3 hours (README)
SKIP = 77
# key: (November 2026, tolerance)
# The tolerances are meant to cover a re-fetched crop and another machine or MPI library (round-off in the
# pressure solves); that is not yet shown on another machine. The scheme's numbers for the same hour (rms 1.37,
# vector rms below 500 m 1.73) are far outside them.
PINNED = {"rms": (0.93, 0.05), "bias": (-0.48, 0.08), "low_speed_bias": (-1.06, 0.08),
          "low_vector_rms": (1.17, 0.05), "low_dir_diff": (0.1, 1.0)}
ITERATIONS = 447           # November 2026; at most 1.5 times this passes
W_MAX = 1.5                # m/s; November 2026: 0.69 (frozen theta in the same hour: 6.3)


def check(path):
    if not os.path.exists(path):
        print(f"check_neutral.py: no score at {path}: see log.txt beside it and run_options.sh's output")
        return 1
    s = json.load(open(path))
    bad = []
    for k, (ref, tol) in PINNED.items():
        ok = abs(s[k] - ref) <= tol
        print(f"{k:16s} {s[k]:7.3f}  November {ref:5.2f} +- {tol:.2f}  {'ok' if ok else 'OUT'}")
        if not ok:
            bad.append(k)
    print(f"steady {s['steady']} after {s['iterations']} iterations (November: steady after {ITERATIONS})")
    if s["steady"] is not True:
        bad.append("steady")
    if s["iterations"] > 1.5 * ITERATIONS:
        bad.append("iterations")
    print(f"w_max            {s['w_max']:7.3f}  below {W_MAX}  {'ok' if s['w_max'] < W_MAX else 'OUT'}")
    if not s["w_max"] < W_MAX:
        bad.append("w_max")
    if bad:
        print(f"check_neutral.py: outside the tolerance: {', '.join(bad)}")
        return 1
    print("check_neutral.py: PASS")
    return 0


def run(a):
    data = os.path.abspath(a.data)
    missing = [f for f in ("awaken.asc", crop_name(HOUR)) if not os.path.exists(os.path.join(data, f))]
    if missing:
        print(f"check_neutral.py: skipped: {', '.join(missing)} not in {data} (python3 prepare.py downloads them)")
        return SKIP
    work = os.path.abspath(a.work)
    py = sys.executable
    subprocess.run([py, os.path.join(HERE, "options.py"), "--data", data, "--out", work, "--hours", HOUR,
                    "--options", OPTION, "--overwrite"], check=True)
    env = dict(os.environ, GUST=os.path.abspath(a.gust), MPI=a.mpi, PYTHON=py)
    subprocess.run([os.path.join(work, "run_options.sh")], env=env, check=True)
    from score_options import score_path
    return check(score_path(work, HOUR, OPTION))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--gust", required=True)
    r.add_argument("--data", default=os.environ.get("GUST_HRRR_NEUTRAL_DATA", HERE),
                   help="the folder prepare.py wrote (default: $GUST_HRRR_NEUTRAL_DATA, else beside this script)")
    r.add_argument("--work", default=os.path.join(HERE, "check_run"))
    r.add_argument("--mpi", default="mpiexec -n 4")
    c = sub.add_parser("check")
    c.add_argument("score", nargs="?", default=os.path.join(HERE, "check_run", HOUR, f"nest_{OPTION}", "score.json"))
    a = ap.parse_args()
    return run(a) if a.cmd == "run" else check(a.score)


if __name__ == "__main__":
    sys.exit(main())
