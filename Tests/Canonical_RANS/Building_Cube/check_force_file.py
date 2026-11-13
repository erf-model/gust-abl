#!/usr/bin/env python3
"""The force file's unsteady CTests (inputs_cube_unsteady, README.md).

    check_force_file.py same-plotfile A B   the plotfiles A and B hold the same bytes (Level_0 data)
    check_force_file.py same-file A B       the files A and B hold the same bytes
    check_force_file.py tail-file A B T0    B starts at t = T0 (the checkpoint's time) and its rows are A's
                                            from T0 on (a restart that finds no file)

A run with buildings.force_file must not change the solution (same-plotfile on the final plotfiles
of the runs with and without it), and a restart must reproduce the force file of the uninterrupted
run (same-file), the rows after the checkpoint written again.
"""
import filecmp
import glob
import os
import sys


def same_plotfile(a, b):
    fa = sorted(os.path.basename(p) for p in glob.glob(os.path.join(a, "Level_0", "Cell_D_*")))
    fb = sorted(os.path.basename(p) for p in glob.glob(os.path.join(b, "Level_0", "Cell_D_*")))
    if not fa or fa != fb:
        print(f"check_force_file: data files differ or are missing: {fa} vs {fb}")
        return False
    bad = [f for f in fa if not filecmp.cmp(os.path.join(a, "Level_0", f), os.path.join(b, "Level_0", f), shallow=False)]
    for f in bad:
        print(f"check_force_file: {f} differs between {a} and {b}")
    return not bad


def same_file(a, b):
    if not filecmp.cmp(a, b, shallow=False):
        print(f"check_force_file: {a} and {b} differ")
        with open(a) as fa, open(b) as fb:
            print("--- " + a + "\n" + fa.read() + "--- " + b + "\n" + fb.read())
        return False
    return True


def tail_file(a, b, t0):
    with open(a) as fa, open(b) as fb:
        ra, rb = fa.read().splitlines(), fb.read().splitlines()
    if len(rb) < 2 or ra[0] != rb[0]:
        print(f"check_force_file: {b} has no rows or another header")
        return False
    if float(rb[1].split(",")[0]) != float(t0):
        print(f"check_force_file: {b} starts at t = {rb[1].split(',')[0]}, not at the checkpoint's {t0}")
        return False
    tail = [r for r in ra[1:] if float(r.split(",")[0]) >= float(t0)]
    if tail != rb[1:]:
        print(f"check_force_file: the rows of {b} are not those of {a} from t = {t0}:\n" + "\n".join(rb[1:]))
        return False
    return True


MODES = {"same-plotfile": same_plotfile, "same-file": same_file, "tail-file": tail_file}


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in MODES or len(sys.argv) != (5 if sys.argv[1] == "tail-file" else 4):
        sys.exit(__doc__)
    ok = MODES[sys.argv[1]](*sys.argv[2:])
    print("check_force_file: " + ("same" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
