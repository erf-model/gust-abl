#!/usr/bin/env python3
"""Classes of 16 HRRR hours at the AWAKEN site (May 10-11 2023) with the default thresholds.

Checks, for each hour: the class and modifiers; L, z_i/L, the bulk Richardson number, the jet fall-off, the geostrophic
shear and the wind tendency against the values pinned below (relative 1e-4, absolute 1e-3); and L against the closed
form worked here from the column's surface fields, independently of hrrr_regime.py. Exit status 0 when all agree.
"""
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
import hrrr_regime as hr  # noqa: E402

# valid time, class, modifiers, L [m], z_i/L, Ri_b (lowest 100 m), jet fall-off [m/s], geostrophic shear
# [m/s per km], tendency [m/s per h]
PINNED = [
    ("2023-05-10T04:00Z", "C5", "B+T", 24646.6, 0.0133949, 0.0162795, 3.27193, 9.37666, 1.95063),
    ("2023-05-10T05:00Z", "C5", "B+T", 770.336, 0.42019, 0.0352612, 3.03097, 7.7462, 1.02299),
    ("2023-05-10T06:00Z", "C5", "T", 361.402, 0.937216, 0.0544044, 3.50518, 0.252354, 2.50473),
    ("2023-05-10T11:00Z", "C6", "B+T", 13.0745, 4.6657, 0.461899, 1.84719, 3.71981, 1.271),
    ("2023-05-10T12:00Z", "C8", "B+T", 121.758, 0.674589, 0.369714, 2.45319, 3.60581, 1.60677),
    ("2023-05-10T13:00Z", "C8", "T", -209.195, -1.32133, 0.130654, 2.52447, 1.9988, 1.97645),
    ("2023-05-10T19:00Z", "C1", "B", -48.8387, -30.7713, -3.67864, -0.158666, 2.02736, 0.902737),
    ("2023-05-10T20:00Z", "C1", "", -66.5844, -24.6164, -2.21376, -0.152228, 1.61723, 0.839004),
    ("2023-05-10T21:00Z", "C1", "B", -87.0638, -20.1118, -1.40347, -0.0687756, 2.33473, 0.985915),
    ("2023-05-10T23:00Z", "C1", "B+T", -156.008, -12.1703, -0.564922, -0.0587169, 2.15763, 2.19118),
    ("2023-05-11T00:00Z", "C7", "B+T", -509.914, -2.58278, -0.093815, 0.567649, 3.01161, 2.18431),
    ("2023-05-11T01:00Z", "C7", "B+T", 420.082, 2.38405, 0.0483591, 0.67045, 2.92845, 1.84979),
    ("2023-05-11T02:00Z", "C4", "B+T", 433.827, 0.841511, 0.0376897, 0.71504, 7.26982, 1.51675),
    ("2023-05-11T08:00Z", "C6", "B+T", 55.6158, 1.8447, 0.221136, 1.20423, 4.54538, 1.09201),
    ("2023-05-11T09:00Z", "C6", "B+T", 63.6391, 1.74997, 0.180434, 1.98827, 3.1925, 1.92288),
    ("2023-05-11T10:00Z", "C5", "B+T", 1.7527, 55.3811, 0.233919, 2.6027, 3.24587, 3.45121),
]


def closed_form_L(s):
    """L = -u*^3 theta_v / (kappa g q_v), q_v = H/(rho cp)(1 + 0.61 q) + 0.61 theta E/(rho Lv), rho = p/(Rd Tv)."""
    q, T, p = s["q2m"], s["t2m"], s["psfc"]
    rho = p / (287.04 * T * (1 + 0.61 * q))
    theta = T * (1e5 / p) ** 0.2857
    qv = s["shf"] / (rho * 1005.0) * (1 + 0.61 * q) + 0.61 * theta * s["lhf"] / (rho * 2.5e6)
    return -s["ustar"] ** 3 * theta * (1 + 0.61 * q) / (0.41 * 9.81 * qv)


def close(a, b):
    return abs(a - b) <= max(1e-4 * abs(b), 1e-3)


def main():
    paths = sorted(glob.glob(os.path.join(HERE, "data", "hrrr_*_f01.json")))
    res = hr.catalogue(paths, hr.default_thresholds(), jobs=1)
    bad = []
    if len(res) != len(PINNED):
        bad.append(f"{len(res)} hours classified, {len(PINNED)} pinned")
    for r, p in zip(res, PINNED):
        if "error" in r:
            bad.append(f"{r['file']}: {r['error']}")
            continue
        d = r["diagnostics"]
        got = (r["valid_time"], r["class"], "+".join(r["modifiers"]), d["L"], d["zi_over_L"], d["rib_low"],
               d["jet"]["drop"], d["geo_shear"], d["tendency"])
        names = ("valid time", "class", "modifiers", "L", "z_i/L", "Ri_b", "jet fall-off", "geostrophic shear",
                 "tendency")
        for n, g, e in zip(names, got, p):
            if (isinstance(e, str) and g != e) or (not isinstance(e, str) and not close(g, e)):
                bad.append(f"{p[0]} {n}: {g} (pinned {e})")
        cl = closed_form_L(hr.read_column(r["file"])["surface"])
        if not math.isclose(d["L"], cl, rel_tol=1e-12):
            bad.append(f"{p[0]} L: {d['L']} (closed form {cl})")
        print(f"{got[0]}  {got[1]}  {got[2] or '-':6s} L {got[3]:10.4g}  z_i/L {got[4]:8.4g}")
    print(hr.count_table(res))
    for b in bad:
        print("FAIL", b)
    print("hrrr_regime_awaken:", "FAIL" if bad else "PASS")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
