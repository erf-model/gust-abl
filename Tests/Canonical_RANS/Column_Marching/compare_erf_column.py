#!/usr/bin/env python3
"""GUST's time-marching column (init.column_time) against ERF's near-1-D columns of the WINSENT stratified inflow.

For each case, the ERF run folder (erf_column.py of PR #78, open: inputs, input_sounding, plt00000 ... and history.txt)
gives the forcing, the surface flux, the sounding and the initial state (the plotfile at t = 0, the horizontal mean of
its 2 x 2 columns). This script writes that state as the GUST deck's bc.inflow_file, runs a copy of inputs_column_marching
with the case's values, and compares the column at 9 h with ERF's: the bulk Richardson number between
--za and --zb and the speed at 300 m (erf_column.py's definitions), plus theta at 10 and 100 m.
"""
import argparse
import math
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

CASES = ["high_stable_G2/q-0.004", "low_stable_G3/q-0.0025", "high_unstable_lapse0/q0.018",
         "low_unstable_lapse0/q0.03"]


def read_erf_mean(plt):
    """Horizontal means of ERF's plotfile: z (cell centres), u, v, theta, k (erf_column.py read_columns)."""
    hdr = ep._read_header(plt)
    names = hdr["names"]
    ncomp, _, boxes, fabs = ep._read_cell_h(os.path.join(plt, "Level_0"))
    nz = hdr["hi"][2] - hdr["lo"][2] + 1
    acc = {n: np.zeros(nz) for n in ("x_velocity", "y_velocity", "theta", "KE")}
    cnt = np.zeros(nz)
    for (blo, bhi), (fname, off) in zip(boxes, fabs):
        flo, fhi, data = ep._read_fab(os.path.join(plt, "Level_0", fname), off, ncomp)
        n = [fhi[d] - flo[d] + 1 for d in range(3)]
        arr = np.asarray(data).reshape(ncomp, n[2], n[1], n[0])
        for name in acc:
            acc[name][flo[2]:fhi[2] + 1] += arr[names.index(name)].sum(axis=(1, 2))
        cnt[flo[2]:fhi[2] + 1] += n[0] * n[1]
    z = hdr["prob_lo"][2] + (np.arange(nz) + 0.5) * hdr["dx"][2]
    return z, acc["x_velocity"] / cnt, acc["y_velocity"] / cnt, acc["theta"] / cnt, acc["KE"] / cnt


def rib(z, u, v, th, za, zb, theta0=300.0):
    du2 = (np.interp(zb, z, u) - np.interp(za, z, u)) ** 2 + (np.interp(zb, z, v) - np.interp(za, z, v)) ** 2
    return 9.81 / theta0 * (np.interp(zb, z, th) - np.interp(za, z, th)) * (zb - za) / max(du2, 1e-12)


def erf_inputs(folder):
    """The keys of ERF's deck the GUST run needs."""
    keys = {}
    with open(os.path.join(folder, "inputs")) as f:
        for line in f:
            m = re.match(r"\s*([\w.]+)\s*=\s*(.+?)\s*$", line.split("#")[0])
            if m:
                keys[m.group(1)] = m.group(2).strip('"')
    snd = np.loadtxt(os.path.join(folder, "input_sounding"), skiprows=1)
    need = ["erf.abl_geo_wind", "erf.latitude", "erf.most.surf_temp_flux", "erf.most.z0", "zhi.theta_grad", "stop_time",
            "erf.theta_ref"]
    missing = [k for k in need if k not in keys]
    if missing:
        sys.exit(f"{folder}/inputs: no {', '.join(missing)}")
    return keys, snd


def write_deck(path, values):
    """inputs_column_marching with the case's values in place of the deck's (a key given twice is refused)."""
    out, seen = [], set()
    with open(os.path.join(HERE, "inputs_column_marching")) as f:
        for line in f:
            key = line.split("=")[0].strip() if "=" in line and not line.lstrip().startswith("#") else None
            if key in values:
                line = f"{key} = {values[key]}\n"
                seen.add(key)
            out.append(line)
    missing = set(values) - seen
    if missing:
        sys.exit(f"inputs_column_marching has no line for {sorted(missing)}")
    with open(path, "w") as f:
        f.writelines(out)


def summary(z, u, v, th, za, zb, theta0):
    return (rib(z, u, v, th, za, zb, theta0), math.hypot(np.interp(300, z, u), np.interp(300, z, v)), np.interp(10, z, th),
            np.interp(100, z, th))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gust", required=True, help="the gust executable")
    p.add_argument("--erf-root", required=True, help="the folder of ERF's column runs (the case paths below it)")
    p.add_argument("--out", required=True, help="a run folder per case is made below it")
    p.add_argument("--cases", nargs="+", default=CASES)
    p.add_argument("--dt", type=float, default=1.0, help="init.column_dt (s)")
    p.add_argument("--za", type=float, default=46.0)
    p.add_argument("--zb", type=float, default=73.0)
    p.add_argument("--extra", nargs="*", default=[], help="further GUST inputs (key=value)")
    a = p.parse_args()

    rows = []
    for case in a.cases:
        folder = os.path.join(a.erf_root, case)
        keys, snd = erf_inputs(folder)
        run = os.path.join(a.out, case.replace("/", "_"))
        os.makedirs(run, exist_ok=True)
        # ERF's state at t = 0 as the inflow table (z from 0: the first cell's values at the ground row)
        z0_, u0, v0, t0, k0 = read_erf_mean(os.path.join(folder, "plt00000"))
        # the deck's levels: 240 cells of 10 m to 2400 m, as ERF's
        if len(z0_) != 240 or np.max(np.abs(z0_ - (np.arange(240) + 0.5) * 10.0)) > 1e-6:
            sys.exit(f"{folder}: ERF's cell centres are not the deck's 240 x 10 m levels")
        tab = np.column_stack([np.concatenate([[0.0], z0_]), np.concatenate([[u0[0]], u0]),
                               np.concatenate([[v0[0]], v0]), np.concatenate([[t0[0]], t0]),
                               np.maximum(np.concatenate([[k0[0]], k0]), 1e-6)])
        np.savetxt(os.path.join(run, "column_marching_initial.dat"), tab, fmt="%.10g",
                   header=f"z u v theta k: ERF's column at t = 0 ({case})")
        ug = keys["erf.abl_geo_wind"].split()
        zs, thetas = snd[:, 0], snd[:, 1]
        case_keys = {"forcing.ug": ug[0], "forcing.vg": ug[1], "forcing.latitude": keys["erf.latitude"],
                     "bc.surface_heat_flux": keys["erf.most.surf_temp_flux"], "bc.surface_z0": keys["erf.most.z0"],
                     "bc.zhi_theta_grad": keys["zhi.theta_grad"], "forcing.theta0": keys["erf.theta_ref"],
                     "bc.theta_profile_z": " ".join(f"{x:g}" for x in zs),
                     "bc.theta_profile": " ".join(f"{x:g}" for x in thetas),
                     "init.column_time": keys["stop_time"], "init.column_dt": f"{a.dt:g}"}
        write_deck(os.path.join(run, "inputs"), case_keys)
        cmd = [a.gust, "inputs"] + a.extra
        with open(os.path.join(run, "run.log"), "w") as log:
            rc = subprocess.run(cmd, cwd=run, stdout=log, stderr=subprocess.STDOUT).returncode
        if rc != 0:
            sys.exit(f"{run}: gust failed (rc {rc}); see run.log")
        g = np.loadtxt(os.path.join(run, "column_marching_table.dat"))
        theta0 = float(keys["erf.theta_ref"])
        gs = summary(g[1:, 0], g[1:, 1], g[1:, 2], g[1:, 3], a.za, a.zb, theta0)
        plts = sorted((int(d[3:]), d) for d in os.listdir(folder) if re.fullmatch(r"plt\d+", d))
        last = os.path.join(folder, plts[-1][1])
        t_last = ep._read_header(last)["time"]
        if abs(t_last - float(keys["stop_time"])) > 1e-6 * float(keys["stop_time"]):
            sys.exit(f"{folder}: the last plotfile is at {t_last} s, not stop_time = {keys['stop_time']} s")
        z, u, v, th, k = read_erf_mean(last)
        es = summary(z, u, v, th, a.za, a.zb, theta0)
        rows.append((case, keys["erf.most.surf_temp_flux"], es, gs))

    print(f"R_B between {a.za:g} and {a.zb:g} m, speed at 300 m and theta at 10/100 m after the run "
          f"(init.column_dt = {a.dt:g} s)")
    print(f"{'case':30s} {'q (K m/s)':>9s} | {'R_B ERF':>8s} {'GUST':>8s} | {'U300 ERF':>8s} {'GUST':>6s} | "
          f"{'th10 ERF':>8s} {'GUST':>8s} | {'th100 ERF':>9s} {'GUST':>8s}")
    for case, q, e, g in rows:
        print(f"{case:30s} {float(q):9.4f} | {e[0]:+8.4f} {g[0]:+8.4f} | {e[1]:8.2f} {g[1]:6.2f} | "
              f"{e[2]:8.3f} {g[2]:8.3f} | {e[3]:9.3f} {g[3]:8.3f}")


if __name__ == "__main__":
    main()
