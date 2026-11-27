#!/usr/bin/env python3
"""The score table of the C1 and C2 options: one row per hour and option from each run's score.json (written through
tools/hrrr_outer.py score by run_options.sh), as a Markdown table on stdout.

    python3 score_options.py [DIR]          DIR: the folder options.py wrote (default: beside this script)

Columns: nest minus HRRR at all heights (speed bias and rms), below 500 m (speed bias, vector rms, direction),
w (rms and max), k / k_HRRR below 500 m, iterations and the steady verdict. The all-heights columns cover the
heights up to the nest's top. A run without a score is listed as missing. Exit status 0.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from options import NEST_OPTIONS  # noqa: E402
from prepare import HOURS  # noqa: E402

COLS = ("option", "bias", "rms", "bias<500", "vec rms<500", "dir<500", "w rms", "w max", "k ratio", "its",
        "steady")


def score_path(d, tag, option):
    """Where run_options.sh leaves the score of one option of one hour."""
    if option == "scheme":
        return os.path.join(d, tag, "score.json")
    return os.path.join(d, tag, f"nest_{option}", "score.json")


def row(tag, option, path):
    if not os.path.exists(path):
        return f"| {tag} | {option} | " + " | ".join(["missing"] + [""] * (len(COLS) - 2)) + " |"
    s = json.load(open(path))
    st = {True: "yes", False: "no", None: "?"}[s["steady"]]
    return (f"| {tag} | {option} | {s['bias']:.2f} | {s['rms']:.2f} | {s['low_speed_bias']:.2f} | "
            f"{s['low_vector_rms']:.2f} | {s['low_dir_diff']:.1f} | {s['w_rms']:.3f} | {s['w_max']:.2f} | "
            f"{s['k_ratio']:.2f} | {s['iterations'] if s['iterations'] is not None else '?'} | {st} |")


def main():
    d = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else HERE)
    print("| hour | " + " | ".join(COLS) + " |")
    print("|---" * (len(COLS) + 1) + "|")
    for tag in HOURS:
        if not os.path.isdir(os.path.join(d, tag)):
            if os.path.exists(os.path.join(d, f"{tag}.yaml")):      # options.py wrote the case but the tool refused it
                print(f"| {tag} | all | refused by hrrr_outer.py (see options.py's output) |" + " |" * (len(COLS) - 2))
            continue
        for option in ("scheme",) + tuple(NEST_OPTIONS):
            print(row(tag, option, score_path(d, tag, option)))


if __name__ == "__main__":
    main()
