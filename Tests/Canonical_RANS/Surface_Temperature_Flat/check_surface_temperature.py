#!/usr/bin/env python3
"""Check the wall model of bc.surface_stability = surface_temperature (rt_surface_temperature_flat_check).

    check_surface_temperature.py RUN_DIR [--rel_tol 1e-9]

RUN_DIR holds one plotfile tsurf_flat_plt<step> and one turbulence plotfile tsurf_flat_turb<step> of
inputs_surface_temperature_flat. For every ground column the script takes the plotfile's first-cell wind
U = |(u, v)| and theta, the column's ground temperature theta_s = 297 + 6 i / 15 K (theta_s_ramp.xyz: one point per
cell centre, so the map holds these values exactly) and the first centre's height d, and solves Monin-Obukhov
similarity again here, by bisection on zeta and with no code of GUST:

    Ri_b = g d (theta - theta_s) / (theta_0 U^2), bounded to [-4, 0.19],
    zeta - Ri_b cm^2 / ch = 0,  cm = max(lz - psi_m, m),  ch = max(lz - psi_h, m),  lz = ln((d + z0)/z0),
    u* = kappa U / cm,

with Businger and Dyer's psi (psi = -5 zeta above zero; Paulson's forms with gamma = 16 below), g = 9.81 m/s^2,
theta_0 = 300 K, kappa = 0.41, z0 = 0.1 m and m = min(lz, 1). The run's u* (the turbulence plotfile's ustar at
k = 0) must match to rel_tol in every column, and so must its theta_s (the map's value) and qflux
(q = -kappa u* (theta - theta_s) / ch, upward > 0). The ramp gives stable columns in the west and unstable ones in
the east; the script also requires both kinds, so the check covers both arms of psi. Last, the flux must have reached
theta: the run starts from theta = 300 K, so every column with |theta_s - 300| > 1 K must have its first-cell theta
moved from 300 K towards theta_s.
"""
import glob
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
# No __pycache__ in the source tree: CI fails a run that leaves files behind.
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402  (ERF's reader, extern_erf/)

KAPPA, G, THETA0, Z0 = 0.41, 9.81, 300.0, 0.1
RIB_MIN, RIB_MAX = -4.0, 0.19


def the_plotfile(directory, prefix):
    found = sorted(p for p in glob.glob(os.path.join(directory, prefix + "*"))
                   if os.path.isdir(p) and os.path.basename(p)[len(prefix):].isdigit())
    if len(found) != 1:
        raise SystemExit("check_surface_temperature: expected one %s* in %s, found %d" % (prefix, directory, len(found)))
    return found[0]


def psi_m(z):
    if z > 0:
        return -5.0 * z
    x = (1 - 16 * z) ** 0.25
    return 2 * math.log(0.5 * (1 + x)) + math.log(0.5 * (1 + x * x)) - 2 * math.atan(x) + 0.5 * math.pi


def psi_h(z):
    if z > 0:
        return -5.0 * z
    return 2 * math.log(0.5 * (1 + math.sqrt(1 - 16 * z)))


def ustar_of(U, d, dtheta):
    lz = math.log((d + Z0) / Z0)
    m = min(lz, 1.0)
    rib = max(RIB_MIN, min(RIB_MAX, G * d * dtheta / (THETA0 * U * U)))
    cm = lambda z: max(lz - psi_m(z), m)
    ch = lambda z: max(lz - psi_h(z), m)
    F = lambda z: z - rib * cm(z) ** 2 / ch(z)
    zeta = 0.0
    if rib != 0:
        a, b = (-100.0, 0.0) if rib < 0 else (0.0, 100.0)
        if rib > 0 and F(b) <= 0:
            zeta = b
        elif rib < 0 and F(a) >= 0:
            zeta = a
        else:
            for _ in range(200):
                mid = 0.5 * (a + b)
                if F(mid) < 0:
                    a = mid
                else:
                    b = mid
            zeta = 0.5 * (a + b)
    ustar = KAPPA * U / cm(zeta)
    return ustar, zeta, -KAPPA * ustar / ch(zeta) * dtheta


def main(argv):
    if len(argv) not in (2, 4) or (len(argv) == 4 and argv[2] != "--rel_tol"):
        raise SystemExit(__doc__)
    rel_tol = float(argv[3]) if len(argv) == 4 else 1e-9
    _, s = erf_plotfile.read_fields(the_plotfile(argv[1], "tsurf_flat_plt"), ["u", "v", "theta", "z"])
    _, t = erf_plotfile.read_fields(the_plotfile(argv[1], "tsurf_flat_turb"), ["ustar", "theta_s", "qflux"])
    nx, ny = len(s["u"]), len(s["u"][0])
    worst, worst_q, worst_ts, stable, unstable, moved, ramp = 0.0, 0.0, 0.0, 0, 0, 0, 0
    for i in range(nx):
        for j in range(ny):
            U = max(math.hypot(s["u"][i][j][0], s["v"][i][j][0]), 0.1)
            ts = 297 + 6 * i / 15
            th = s["theta"][i][j][0]
            want, zeta, q = ustar_of(U, s["z"][i][j][0], th - ts)
            worst = max(worst, abs(t["ustar"][i][j][0] - want) / want)
            worst_q = max(worst_q, abs(t["qflux"][i][j][0] - q) / max(abs(q), 1e-6))
            worst_ts = max(worst_ts, abs(t["theta_s"][i][j][0] - ts) / ts)
            stable += zeta > 0.01
            unstable += zeta < -0.01
            if abs(ts - 300) > 1:
                ramp += 1
                moved += (th - 300) * (ts - 300) > 0 and abs(th - 300) < abs(ts - 300)
    ok = (max(worst, worst_q, worst_ts) <= rel_tol and stable > 0 and unstable > 0 and moved == ramp)
    print("check_surface_temperature: %d columns, %d stable and %d unstable (|zeta| > 0.01); max relative "
          "difference from MOST: u* %.3e, q %.3e, theta_s %.3e (tolerance %.1e); theta moved towards theta_s in "
          "%d of %d columns with |theta_s - 300| > 1 K: %s"
          % (nx * ny, stable, unstable, worst, worst_q, worst_ts, rel_tol, moved, ramp, "ok" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
