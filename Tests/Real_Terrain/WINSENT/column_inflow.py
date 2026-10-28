#!/usr/bin/env python3
"""Make a 3-D inflow file from GUST's neutral inflow column, turned to the lidar's direction.

Runs the column deck (inputs_winsent_column: flat ground, Coriolis at 48.7 N, the geostrophic wind found by GUST
so the column has --target-speed at --target-height) and writes ``z u v theta k`` with the column's speed and k and
the lidar case mean's direction at each height (the lidar's 100 m direction below 100 m). With --scale G the target
is G times the benchmark power law at the target height: kynema-sgf's way of driving its precursor faster than the
lidar so the wind that reaches the lidar after the terrain matches it.
"""
import argparse
import os
import shutil
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gust", required=True, help="the gust executable")
    p.add_argument("--zlevels", required=True, help="the 3-D run's zlevels.txt")
    p.add_argument("--lidar", required=True, help="the lidar case mean: z u v theta k (prepare_gust.py's inflow.txt)")
    p.add_argument("--out", required=True, help="folder to run the column in and write inflow.txt to")
    p.add_argument("--powerlaw", type=float, nargs=2, default=[1.033, 0.323], help="benchmark C a: U = C h^a")
    p.add_argument("--scale", type=float, default=1.0, help="G: the target is G times the power law (default 1)")
    p.add_argument("--target-height", type=float, default=300.0)
    p.add_argument("--z0", type=float, default=0.1, help="the column's surface roughness (m)")
    p.add_argument("--mol", type=float, default=0.0,
                   help="Obukhov length L of a stratified column (m, < 0 unstable); 0 = neutral (default)")
    p.add_argument("--rayleigh-z", type=float, default=1900.0, help="stratified column: start of the theta damping (m)")
    p.add_argument("--rayleigh-coef", type=float, default=0.2, help="stratified column: its coefficient (1/s)")
    p.add_argument("--gamma", type=float, default=0.0,
                   help="prescribed stable theta: the neutral column held at theta = 300 K + gamma z (K/km); "
                        "solves u, v, k with its buoyancy (van der Laan et al. 2024); not with --mol")
    p.add_argument("--theta-z", type=float, nargs="+", help="prescribed theta sounding: heights (m), with --theta")
    p.add_argument("--theta", type=float, nargs="+", help="prescribed theta sounding: values (K), with --theta-z")
    a = p.parse_args()
    if (a.theta_z is None) != (a.theta is None) or (a.theta and len(a.theta) != len(a.theta_z)):
        sys.exit("--theta-z and --theta go together, with as many values")
    if a.theta and (a.gamma or a.mol):
        sys.exit("--theta excludes --gamma and --mol")
    if a.gamma != 0 and a.mol != 0:
        sys.exit("--gamma and --mol exclude each other")
    if a.scale <= 0:
        sys.exit("--scale must be > 0")
    target = a.scale * a.powerlaw[0] * a.target_height ** a.powerlaw[1]

    os.makedirs(a.out, exist_ok=True)
    shutil.copy(os.path.join(HERE, "inputs_winsent_column"), a.out)
    shutil.copy(a.zlevels, os.path.join(a.out, "zlevels.txt"))
    cmd = [a.gust, "inputs_winsent_column", f"bc.surface_z0={a.z0}", f"bc.inflow_z0={a.z0}",
           "steady.max_iterations=1", f"init.column_target_height={a.target_height}"]
    if a.mol != 0:
        cmd += [f"init.column_mol={a.mol}", "bc.surface_stability=heat_flux", f"bc.rayleigh_zdamp={a.rayleigh_z}",
                f"bc.rayleigh_dampcoef={a.rayleigh_coef}", "bc.rayleigh_fields=w theta",
                "init.column_max_iterations=50000"]
    if a.gamma != 0:
        cmd += ["bc.theta_profile_z=0 3600", f"bc.theta_profile=300 {300.0 + 3.6 * a.gamma}",
                f"bc.zhi_theta_grad={a.gamma / 1000.0}", "init.column_max_iterations=50000"]
    if a.theta:
        cmd += ["bc.theta_profile_z=" + " ".join(f"{z:g}" for z in a.theta_z),
                "bc.theta_profile=" + " ".join(f"{t:g}" for t in a.theta), "init.column_max_iterations=50000"]
    deck = open(os.path.join(HERE, "inputs_winsent_column")).read()
    if "init.column_target_u" in deck:  # the deck's 6.50 m/s is replaced through a copy without the key
        lines = [ln for ln in deck.splitlines() if not ln.startswith(("init.column_target_u", "init.column_target_height"))]
        open(os.path.join(a.out, "inputs_winsent_column"), "w").write("\n".join(lines) + "\n")
        cmd.append(f"init.column_target_u={target}")
    with open(os.path.join(a.out, "column.log"), "w") as log:
        subprocess.run(cmd, cwd=a.out, stdout=log, stderr=subprocess.STDOUT, check=True)

    c = np.loadtxt(os.path.join(a.out, "column_guess.dat"))
    lid = np.loadtxt(a.lidar)
    z, s = c[:, 0], np.hypot(c[:, 1], c[:, 2])  # theta (column 3) is the column's: stratified with --mol
    m = lid[:, 0] >= 100.0
    th = np.interp(z, lid[m, 0], np.unwrap(np.arctan2(lid[m, 2], lid[m, 1])))
    out = np.column_stack([z, s * np.cos(th), s * np.sin(th), c[:, 3], c[:, 4]])
    np.savetxt(os.path.join(a.out, "inflow.txt"), out, fmt="%10.3f %9.4f %9.4f %8.2f %8.4f",
               header=f"z u v theta k: GUST {'theta held at z ' + str(a.theta_z) + ' K ' + str(a.theta) if a.theta else f'theta = 300 K + {a.gamma:g} K/km z' if a.gamma else ('neutral' if a.mol == 0 else f'L = {a.mol:g} m')} column "
                      f"(z0 {a.z0}, 48.7 N, {target:.3f} m/s at "
                      f"{a.target_height:g} m = {a.scale:g} x the benchmark power law); the column's speed and k, "
                      "the lidar case mean's direction (its 100 m direction below 100 m)")
    hdr = open(os.path.join(a.out, "column_guess.dat")).readline().strip()
    print(hdr)
    for h in (10, 100, 300, 1000):
        i = int(np.argmin(abs(z - h)))
        print(f"G {a.scale:g}: {z[i]:7.1f} m  {s[i]:.2f} m/s  k {c[i, 4]:.3f}")


if __name__ == "__main__":
    main()
