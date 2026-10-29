#!/usr/bin/env python3
"""A stable inflow column from amrTerrain's 1-D solver with a surface cooling rate (GABLS-style), held as the inflow.

GUST's own column marches to a true steady state, which a cooled stable column does not have (the surface keeps
removing heat). amrTerrain's 1-D solver (src/backend/amr1DSolver.py, the same Axell-Liungman closure) steps in
physical time and has a cooling-rate surface (heat_flux_model mode 3: theta_s = theta_ref + rate t / 3600); its run
stops when u and v stop changing while theta keeps cooling: the quasi-steady state of GABLS1 (Beare et al. 2006),
whose profiles are steady after about 9 h. This script runs it to that state (or --hours), records the bulk
Richardson number between --za and --zb on the way, and writes the final profile as a GUST inflow file
``z u v theta k``, turned at each height to the lidar's direction like column_inflow.py.
"""
import argparse
import math
import os
import sys

import numpy as np


def rib(z, u, v, th, za, zb, theta0=300.0):
    du2 = (np.interp(zb, z, u) - np.interp(za, z, u)) ** 2 + (np.interp(zb, z, v) - np.interp(za, z, v)) ** 2
    return 9.81 / theta0 * (np.interp(zb, z, th) - np.interp(za, z, th)) * (zb - za) / max(du2, 1e-12)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--amrterrain", required=True, help="amrTerrain's src/backend folder")
    p.add_argument("--out", required=True)
    p.add_argument("--ug", type=float, nargs=2, required=True, help="geostrophic wind in this case's frame (m/s)")
    g = p.add_mutually_exclusive_group(required=True)
    g.add_argument("--rate", type=float, help="surface cooling rate (K/h, < 0 cools): mode 3, the stable column")
    g.add_argument("--flux", type=float, help="surface heat flux (K m/s, > 0 heats): mode 1, the unstable column")
    p.add_argument("--lidar", required=True, help="class lidar mean (lidar_mean.py) for the direction")
    p.add_argument("--hours", type=float, default=9.0, help="longest run (h)")
    p.add_argument("--tol", type=float, default=1e-4, help="the solver's convergence on sum(u), sum(v)")
    p.add_argument("--z0", type=float, default=0.1)
    p.add_argument("--lat", type=float, default=48.7)
    p.add_argument("--ztop", type=float, default=2400.0)
    p.add_argument("--nz", type=int, default=241)
    p.add_argument("--inversion", type=float, nargs=3, default=[1000.0, 100.0, 3e-3],
                   metavar=("HEIGHT", "WIDTH", "LAPSE"), help="initial capping inversion (0.08 K/m across WIDTH, "
                   "fixed in the solver) and the lapse rate above (K/m)")
    p.add_argument("--za", type=float, default=46.0)
    p.add_argument("--zb", type=float, default=73.0)
    a = p.parse_args()

    sys.path.insert(0, a.amrterrain)
    from amr1DSolver import amr1dSolver  # noqa: E402

    os.makedirs(a.out, exist_ok=True)
    s = amr1dSolver(a.nz, a.ztop, a.z0, 0.0, os.path.join(a.out, "amr1d_profile.txt"))
    s.initialize_physics(a.ug[0], a.ug[1], 300.0, 300.0, 0.4, 0.25, 1000.0)
    s.initialize_coriolis(a.lat)
    s.temperature_inversion(a.inversion[0], a.inversion[1], 0.08 * a.inversion[1], a.inversion[2])
    s.heat_flux_model(3, a.rate) if a.rate is not None else s.heat_flux_model(1, a.flux)
    what = f"cooling {a.rate:g} K/h" if a.rate is not None else f"heat flux {a.flux:g} K m/s"
    # run_simulation's loop, with a record of R_B and the surface state every 15 min
    s.err, s.converge, s.counter, s.start_time = a.tol, False, 0, 0.0
    s.initialize_simulation()
    s.dt = 0.8 * s.dz / math.hypot(s.ug, s.vg)
    s.errvelx, s.errvely = np.sum(s.ux), np.sum(s.uy)
    s.errtke, s.errnut, s.errtemp = np.sum(s.tke), np.sum(s.nut), np.sum(s.temperature)
    log = open(os.path.join(a.out, "history.txt"), "w")
    log.write("# t (h)  R_B  theta_s (K)  Qh (K m/s)  L (m)  U300 (m/s)\n")
    next_rec = 0.0
    while s.start_time <= a.hours * 3600.0:
        s.start_time += s.dt
        s.compute_similarity()
        for i in range(s.lower + 1, s.nz - 1):
            s.update_windspeed_x(i, s.dt)
            s.update_windspeed_y(i, s.dt)
            s.update_temperature(i, s.dt)
            s.update_turbulence(i, s.dt)
        s.dt = 0.5 * s.dz / math.hypot(max(s.ux), max(s.uy))
        s.counter += 1
        s.compute_average()
        s.compute_error()
        if s.start_time >= next_rec or s.converge:
            r = rib(s.z, s.ux, s.uy, s.temperature, a.za, a.zb)
            u300 = math.hypot(np.interp(300, s.z, s.ux), np.interp(300, s.z, s.uy))
            log.write(f"{s.start_time / 3600:7.3f} {r:+8.4f} {s.temperature[s.lower]:9.3f} {s.Qh:+10.5f} "
                      f"{s.mo_length:12.1f} {u300:6.2f}\n")
            log.flush()
            next_rec += 900.0
        if s.converge:
            break
    state = "converged (u, v steady)" if s.converge else f"stopped at {a.hours:g} h"
    r = rib(s.z, s.ux, s.uy, s.temperature, a.za, a.zb)

    lid = np.loadtxt(a.lidar)
    m = lid[:, 0] >= 100.0
    th = np.interp(s.z, lid[m, 0], np.unwrap(np.arctan2(lid[m, 2], lid[m, 1])))
    sp = np.hypot(s.ux, s.uy)
    out = np.column_stack([s.z, sp * np.cos(th), sp * np.sin(th), s.temperature, np.maximum(s.tke, 1e-6)])
    out[0, 1:3] = 0.0
    np.savetxt(os.path.join(a.out, "inflow.txt"), out, fmt="%10.3f %9.4f %9.4f %9.4f %9.5f",
               header=f"z u v theta k: amrTerrain 1-D column, {what}, {state} at "
                      f"{s.start_time / 3600:.2f} h; theta_s {s.temperature[s.lower]:.3f} K, Qh {s.Qh:.5f} K m/s, "
                      f"L {s.mo_length:.1f} m, R_B({a.za:g}-{a.zb:g} m) {r:.4f}; the lidar's direction")
    print(f"{what}: {state} at {s.start_time / 3600:.2f} h; R_B {r:+.3f}; Qh {s.Qh:+.5f} K m/s; "
          f"L {s.mo_length:.0f} m; U(300 m) {math.hypot(np.interp(300, s.z, s.ux), np.interp(300, s.z, s.uy)):.2f}")


if __name__ == "__main__":
    main()
