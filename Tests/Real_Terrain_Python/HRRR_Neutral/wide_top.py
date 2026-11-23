#!/usr/bin/env python3
"""tools/hrrr_outer.py with its domain-top limits widened, for the top options of this case only.

The tool refuses a domain top below 3000 m (MIN_TOP) and above 3600 m (HTOP, the top of its height table). The
options top1500 and top5000 test those limits: this wrapper sets MIN_TOP to 1000 m, extends the height table to
6000 m (every 2 m, as below 3600 m), adds theta_v profile points above 3 km, then runs the tool's command line
unchanged. Every other check stays, including the refusal of a top above HRRR's highest fetched level. Below
3600 m the table, and so every value the tool writes, is the tool's own.

    python3 wide_top.py <hrrr_outer.py arguments>
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..", "..", "..", "tools")))
import hrrr_outer as ho  # noqa: E402

WIDE_MIN_TOP, WIDE_HTOP = 1000.0, 6000.0


def widen():
    """Widen the limits of the hrrr_outer module in this process; returns what restore() needs."""
    old = (ho.MIN_TOP, ho.HTOP, ho.HG, ho.PROFILE_HEIGHTS)
    ho.MIN_TOP = WIDE_MIN_TOP
    ho.HTOP = WIDE_HTOP
    ho.HG = np.concatenate([old[2], np.arange(old[2][-1] + ho.DH, WIDE_HTOP + 0.5 * ho.DH, ho.DH)])
    ho.PROFILE_HEIGHTS = old[3] + (3500.0, 4000.0, 4500.0, 5000.0, 5500.0)     # theta_v points above 3 km
    return old


def restore(old):
    ho.MIN_TOP, ho.HTOP, ho.HG, ho.PROFILE_HEIGHTS = old


if __name__ == "__main__":
    widen()
    sys.exit(ho.main())
