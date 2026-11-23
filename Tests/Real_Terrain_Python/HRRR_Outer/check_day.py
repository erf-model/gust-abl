#!/usr/bin/env python3
"""The day hour (2023-05-10 20 UTC) against the November 2026 numbers: the chain of case.yaml, run by its run.sh,
gives a score.json; this script checks it.

    python3 check_day.py [day20/score.json]

The November 2026 run (theta solved, Coriolis on, cap 100 m, the same meshes, written by the scripts this case
replaces) gave: speed rms 0.47 m/s and bias +0.02 m/s at all heights; below 500 m a speed bias of 0.00 m/s, a vector
rms of 0.57 m/s and a direction difference of 0.3 degrees; steady in 239 iterations. The tolerances below cover the
change of horizontal interpolation (bilinear in HRRR's grid, not triangles: up to 0.11 m/s in the first guess) and
a re-fetched crop. Exit status 0 when every number is inside its tolerance.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
# key: (November 2026, tolerance)
PINNED = {"rms": (0.47, 0.05), "bias": (0.02, 0.08), "low_speed_bias": (0.00, 0.08),
          "low_vector_rms": (0.57, 0.05), "low_dir_diff": (0.3, 1.0)}


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "day20", "score.json")
    s = json.load(open(path))
    bad = []
    for k, (ref, tol) in PINNED.items():
        ok = abs(s[k] - ref) <= tol
        print(f"{k:16s} {s[k]:7.3f}  November {ref:5.2f} +- {tol:.2f}  {'ok' if ok else 'OUT'}")
        if not ok:
            bad.append(k)
    print(f"steady {s['steady']} after {s['iterations']} iterations (November: steady after 239)")
    if s["steady"] is not True:
        bad.append("steady")
    if bad:
        sys.exit(f"check_day.py: outside the tolerance: {', '.join(bad)}")


if __name__ == "__main__":
    main()
