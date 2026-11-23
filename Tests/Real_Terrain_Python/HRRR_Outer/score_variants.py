#!/usr/bin/env python3
"""The score table of the variant matrix: one row per hour and variant from each run's score.json (written by
run_matrix.sh through tools/hrrr_outer.py score), as a Markdown table on stdout.

    python3 score_variants.py [DIR]          DIR: the folder variants.py wrote (default: beside this script)

Columns: nest minus HRRR at all heights (speed bias and rms), below 500 m (speed bias, vector rms, direction),
w (rms and max), k / k_HRRR below 500 m, iterations and the steady verdict. A run without a score is listed as
missing. Exit status 0.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from prepare import HOURS  # noqa: E402

COLS = ("variant", "bias", "rms", "bias<500", "vec rms<500", "dir<500", "w rms", "w max", "k ratio", "its",
        "steady")


def row(tag, v, path):
    if not os.path.exists(path):
        return f"| {tag} | {v} | " + " | ".join(["missing"] + [""] * (len(COLS) - 2)) + " |"
    s = json.load(open(path))
    st = {True: "yes", False: "no", None: "?"}[s["steady"]]
    return (f"| {tag} | {v} | {s['bias']:.2f} | {s['rms']:.2f} | {s['low_speed_bias']:.2f} | "
            f"{s['low_vector_rms']:.2f} | {s['low_dir_diff']:.1f} | {s['w_rms']:.3f} | {s['w_max']:.2f} | "
            f"{s['k_ratio']:.2f} | {s['iterations'] if s['iterations'] is not None else '?'} | {st} |")


def main():
    d = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else HERE)
    print("| hour | " + " | ".join(COLS) + " |")
    print("|---" * (len(COLS) + 1) + "|")
    for tag in HOURS:
        if not os.path.isdir(os.path.join(d, tag)):
            if os.path.exists(os.path.join(d, f"{tag}.yaml")):     # variants.py wrote the case but the tool refused it
                print(f"| {tag} | all | refused by hrrr_outer.py (see variants.py's output) |" + " |" * (len(COLS) - 2))
            continue
        print(row(tag, "S", os.path.join(d, tag, "score.json")))
        for v in "ABCDEF":
            print(row(tag, v, os.path.join(d, tag, f"nest_{v}", "score.json")))


if __name__ == "__main__":
    main()
