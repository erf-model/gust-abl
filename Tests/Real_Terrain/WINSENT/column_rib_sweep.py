#!/usr/bin/env python3
"""Sweep the Obukhov length of GUST's inflow column and report the bulk Richardson number between two heights.

The benchmark classifies each 10-minute period by the bulk Richardson number at the plateau mast MMNW between
46 and 73 m (winsent-benchmark.readthedocs.io, Phase I). This runs the column deck (inputs_winsent_column: flat
ground, Coriolis at 48.7 N, the geostrophic wind found for --target-speed at 300 m) once per Obukhov length, with a
stratified surface (init.column_mol = L, bc.surface_stability = heat_flux, a Rayleigh layer on theta as the heat
sink the stratified column needs) and computes

    R_B = (g / theta_0) (theta_b - theta_a) (z_b - z_a) / ((u_b - u_a)^2 + (v_b - v_a)^2)

between --za and --zb from the converged table (the benchmark does not state its formula; this is the usual bulk
form between two levels). A first check only: a flat column has no plateau, forest or escarpment.
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def rib(table, za, zb, g=9.81, theta0=300.0):
    z, u, v, th = table[:, 0], table[:, 1], table[:, 2], table[:, 3]
    ua, ub = np.interp(za, z, u), np.interp(zb, z, u)
    va, vb = np.interp(za, z, v), np.interp(zb, z, v)
    ta, tb = np.interp(za, z, th), np.interp(zb, z, th)
    du2 = (ub - ua) ** 2 + (vb - va) ** 2
    return g / theta0 * (tb - ta) * (zb - za) / du2, np.hypot(ua, va), np.hypot(ub, vb), tb - ta


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gust", required=True)
    p.add_argument("--zlevels", required=True)
    p.add_argument("--out", required=True, help="folder for the runs")
    p.add_argument("--target-speed", type=float, required=True, help="column speed at 300 m (m/s)")
    p.add_argument("--mol", type=float, nargs="+", default=[], help="Obukhov lengths L (m); 0 = neutral")
    p.add_argument("--gamma", type=float, nargs="+", default=[],
                   help="prescribed-theta columns: theta = 300 K + gamma z, gamma in K/km (the neutral column holds "
                        "theta at this sounding and solves u, v, k with its buoyancy; van der Laan et al. 2024)")
    p.add_argument("--z0", type=float, default=0.1)
    p.add_argument("--za", type=float, default=46.0)
    p.add_argument("--zb", type=float, default=73.0)
    p.add_argument("--rayleigh-z", type=float, default=1900.0, help="start of the Rayleigh layer on w and theta (m)")
    p.add_argument("--rayleigh-coef", type=float, default=0.2, help="its coefficient (1/s)")
    p.add_argument("--extra", nargs="*", default=[], help="further key=value inputs for the stratified columns")
    a = p.parse_args()

    deck = [ln for ln in open(os.path.join(HERE, "inputs_winsent_column")).read().splitlines()
            if not ln.startswith(("init.column_target_u", "init.column_target_height"))]
    print(f"target {a.target_speed:.2f} m/s at 300 m; R_B between {a.za:g} and {a.zb:g} m")
    print("    case         u* (m/s)   q (K m/s)  U(za)  U(zb)  dtheta (K)    R_B")
    runs = [(f"L{L:g}", L, None) for L in a.mol] + [(f"gamma{g:g}", 0.0, g) for g in a.gamma]
    for label, L, gamma in runs:
        d = os.path.join(a.out, label)
        os.makedirs(d, exist_ok=True)
        open(os.path.join(d, "inputs_winsent_column"), "w").write("\n".join(deck) + "\n")
        shutil.copy(a.zlevels, os.path.join(d, "zlevels.txt"))
        cmd = [a.gust, "inputs_winsent_column", f"bc.surface_z0={a.z0}", f"bc.inflow_z0={a.z0}",
               "steady.max_iterations=1", f"init.column_target_u={a.target_speed}", "init.column_target_height=300"]
        if L != 0:
            cmd += [f"init.column_mol={L}", "bc.surface_stability=heat_flux", f"bc.rayleigh_zdamp={a.rayleigh_z}",
                    f"bc.rayleigh_dampcoef={a.rayleigh_coef}", "bc.rayleigh_fields=w theta"] + a.extra
        if gamma is not None:
            top = 300.0 + gamma * 3.6  # the sounding to 3600 m, above the domain top
            cmd += ["bc.theta_profile_z=0 3600", f"bc.theta_profile=300 {top}", f"bc.zhi_theta_grad={gamma / 1000.0}",
                    "init.column_max_iterations=50000"]
        with open(os.path.join(d, "column.log"), "w") as log:
            rc = subprocess.run(cmd, cwd=d, stdout=log, stderr=subprocess.STDOUT).returncode
        if rc != 0 or not os.path.exists(os.path.join(d, "column_guess.dat")):
            print(f"  {label:10s}   failed (see {d}/column.log)")
            continue
        hdr = open(os.path.join(d, "column_guess.dat")).readline()
        t = np.loadtxt(os.path.join(d, "column_guess.dat"))
        r, ua, ub, dth = rib(t, a.za, a.zb)

        def num(key):
            i = hdr.find(key)
            return float(hdr[i + len(key):].split()[0].rstrip(",;")) if i >= 0 else float("nan")
        print(f"  {label:10s}   {num('u*'):9.3f}  {num('heat flux'):10.4f}  {ua:5.2f}  {ub:5.2f}  {dth:+10.3f}  {r:+7.3f}")


if __name__ == "__main__":
    main()
