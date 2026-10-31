#!/usr/bin/env python3
"""A WINSENT case folder without the canopy, the forest as roughness (PT-8.6.1):

    make_z0_only.py CASE OUT [--z0 1.5] [--lai-min 2.5]

Copies CASE (a neutral case folder from prepare_gust.py, holding a copy of the deck as ``deck``,
for example of inputs_winsent_z0map) to OUT, sets the roughness map to z0 where canopy_lai.xyz is at least
--lai-min (0.1 of the 15 m trees; the LAI map has small fringe values), and removes the canopy keys
from the deck. plotfiles and checkpoints are not copied.
"""
import argparse
import os
import shutil


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("case")
    ap.add_argument("out")
    ap.add_argument("--z0", type=float, default=1.5)
    ap.add_argument("--lai-min", type=float, default=2.5)
    args = ap.parse_args()
    if args.z0 <= 0.0:
        raise SystemExit("make_z0_only: --z0 must be above 0")
    shutil.copytree(args.case, args.out, ignore=shutil.ignore_patterns("*_plt*", "*_chk*", "*_turb*", "run.log"))
    with open(os.path.join(args.case, "roughness.xyz")) as r, open(os.path.join(args.case, "canopy_lai.xyz")) as a:
        rough, lai = r.readlines(), a.readlines()
    if len(rough) != len(lai):
        raise SystemExit("make_z0_only: roughness.xyz and canopy_lai.xyz differ in length")
    n = 0
    with open(os.path.join(args.out, "roughness.xyz"), "w") as f:
        for lr, la in zip(rough, lai):
            xr, xa = lr.split(), la.split()
            if xr[:2] != xa[:2]:
                raise SystemExit("make_z0_only: the tables' points differ: %s, %s" % (lr.strip(), la.strip()))
            if float(xa[2]) >= args.lai_min:
                n += 1
                f.write("%s %s %g\n" % (xr[0], xr[1], args.z0))
            else:
                f.write(lr)
    deck = os.path.join(args.out, "deck")
    if not os.path.isfile(deck):
        raise SystemExit("make_z0_only: copy the deck to %s/deck first" % args.case)
    with open(deck) as f:
        lines = f.readlines()
    with open(deck, "w") as f:
        f.writelines(ln for ln in lines if not ln.startswith("canopy."))
    print("make_z0_only: z0 = %g m on %d of %d points" % (args.z0, n, len(rough)))


if __name__ == "__main__":
    main()
