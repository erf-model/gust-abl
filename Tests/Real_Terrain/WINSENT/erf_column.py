#!/usr/bin/env python3
"""A near-1-D ERF column with a fixed surface heat flux, and its profile as a GUST inflow (frozen-theta runs).

ERF (the WINSENT build: the one-equation k RANS closure, as GUST's) on a 2 x 2 periodic column of uniform levels to
--ztop, flat ground, MOST with z0 and a fixed kinematic surface heat flux (erf.most.surf_temp_flux), Coriolis at
--lat and a geostrophic wind; the initial theta is 300 K to the inversion, a jump, then a lapse rate (GABLS1's layout,
Cuxart et al. 2006). Plotfiles every --plot-every s; this script reads each, records the bulk Richardson number
between --za and --zb, and writes the last plotfile's profile as ``z u v theta k`` turned at each height to the
lidar's direction (column_inflow.py's rule): the inflow and the frozen theta of the GUST run (simplec.theta = frozen).
"""
import argparse
import glob
import math
import os
import struct
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

DECK = """# near-1-D ERF column for a GUST frozen-theta inflow (written by erf_column.py)
stop_time = {stop}
amrex.fpe_trap_invalid = 0
fabarray.mfiter_tile_size = 1024 1024 1024
geometry.prob_lo = 0 0 0
geometry.prob_hi = 64 64 {ztop}
amr.n_cell = 2 2 {nz}
geometry.is_periodic = 1 1 0
erf.anelastic = 1
erf.anelastic_type = MidPoint
zlo.type = "surface_layer"
erf.most.z0 = {z0}
erf.most.surf_temp_flux = {flux}
zhi.type = "SlipWall"
zhi.theta_grad = {lapse}
erf.init_type = "input_sounding"
erf.sounding_type = Ideal
erf.input_sounding_file = "input_sounding"
erf.fixed_dt = {dt}
erf.sum_interval = -1
erf.v = 0
amr.max_level = 0
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_per_1 = {plot_every}
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE Kmv
erf.dycore_vert_adv_type = "Upwind_3rd"
erf.dryscal_vert_adv_type = "Upwind_3rd"
erf.molec_diff_type = "None"
erf.use_gravity = true
erf.use_coriolis = true
erf.latitude = {lat}
erf.abl_driver_type = "GeostrophicWind"
erf.abl_geo_wind = {ug} {vg} 0.0
erf.les_type = "None"
erf.rans_type = "kEqn"
erf.init_tke_from_ustar = true
erf.implicit_tke_dissipation = true
erf.theta_ref = 300.0
erf.rans_realizability = 0.309
"""


def read_columns(plt):
    """Horizontal means over the 2 x 2 column: z (cell centres), u, v, theta, k."""
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
    lo, dx = hdr["prob_lo"], hdr["dx"]
    z = lo[2] + (np.arange(nz) + 0.5) * dx[2]
    # ERF plots KE as rho k / rho (erf_derKE): k itself
    return z, acc["x_velocity"] / cnt, acc["y_velocity"] / cnt, acc["theta"] / cnt, acc["KE"] / cnt, hdr


def rib(z, u, v, th, za, zb, theta0=300.0):
    du2 = (np.interp(zb, z, u) - np.interp(za, z, u)) ** 2 + (np.interp(zb, z, v) - np.interp(za, z, v)) ** 2
    return 9.81 / theta0 * (np.interp(zb, z, th) - np.interp(za, z, th)) * (zb - za) / max(du2, 1e-12)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--erf", required=True, help="erf_exec")
    p.add_argument("--out", required=True)
    p.add_argument("--flux", type=float, required=True, help="kinematic surface heat flux (K m/s, > 0 heats)")
    p.add_argument("--ug", type=float, nargs=2, required=True, help="geostrophic wind in the case frame (m/s)")
    p.add_argument("--lidar", required=True, help="class lidar mean (lidar_mean.py), for the direction")
    p.add_argument("--hours", type=float, default=9.0)
    p.add_argument("--plot-every", type=float, default=900.0)
    p.add_argument("--dt", type=float, default=1.0)
    p.add_argument("--z0", type=float, default=0.1)
    p.add_argument("--lat", type=float, default=48.7)
    p.add_argument("--ztop", type=float, default=2400.0)
    p.add_argument("--nz", type=int, default=240)
    p.add_argument("--inversion", type=float, nargs=3, default=[1000.0, 100.0, 8.0], metavar=("BASE", "DEPTH", "JUMP"))
    p.add_argument("--lapse", type=float, default=0.003, help="theta gradient above the inversion (K/m)")
    p.add_argument("--za", type=float, default=46.0)
    p.add_argument("--zb", type=float, default=73.0)
    p.add_argument("--ranks", type=int, default=1)
    a = p.parse_args()

    os.makedirs(a.out, exist_ok=True)
    base, depth, jump = a.inversion
    zs = [0.0, base, base + depth, a.ztop]
    thetas = [300.0, 300.0, 300.0 + jump, 300.0 + jump + a.lapse * (a.ztop - base - depth)]
    with open(os.path.join(a.out, "input_sounding"), "w") as f:
        f.write(f"1000.0 300.0 0.0\n")
        for z, th in zip(zs, thetas):
            f.write(f"{z:.1f} {th:.4f} 0.0 {a.ug[0]:.5f} {a.ug[1]:.5f}\n")
    with open(os.path.join(a.out, "inputs"), "w") as f:
        f.write(DECK.format(stop=a.hours * 3600.0, ztop=a.ztop, nz=a.nz, z0=a.z0, flux=a.flux, lapse=a.lapse, dt=a.dt,
                            plot_every=a.plot_every, lat=a.lat, ug=a.ug[0], vg=a.ug[1]))
    cmd = (["mpiexec", "-n", str(a.ranks)] if a.ranks > 1 else []) + [a.erf, "inputs"]
    with open(os.path.join(a.out, "run.log"), "w") as log:
        rc = subprocess.run(cmd, cwd=a.out, stdout=log, stderr=subprocess.STDOUT).returncode
    plts = sorted(glob.glob(os.path.join(a.out, "plt[0-9]*")))
    if rc != 0 or not plts:
        sys.exit(f"{a.out}: ERF failed (rc {rc}); see run.log")
    with open(os.path.join(a.out, "history.txt"), "w") as h:
        h.write("# t (h)  R_B  U300 (m/s)  theta(10 m)  theta(100 m)\n")
        for plt in plts:
            z, u, v, th, k, hdr = read_columns(plt)
            h.write(f"{hdr['time'] / 3600:7.3f} {rib(z, u, v, th, a.za, a.zb):+8.4f} "
                    f"{math.hypot(np.interp(300, z, u), np.interp(300, z, v)):6.2f} "
                    f"{np.interp(10, z, th):9.3f} {np.interp(100, z, th):9.3f}\n")
    z, u, v, th, k, hdr = read_columns(plts[-1])
    lid = np.loadtxt(a.lidar)
    m = lid[:, 0] >= 100.0
    zz = np.concatenate([[0.0], z])
    sp = np.concatenate([[0.0], np.hypot(u, v)])
    tz = np.interp(zz, lid[m, 0], np.unwrap(np.arctan2(lid[m, 2], lid[m, 1])))
    thz = np.concatenate([[th[0]], th])
    kz = np.maximum(np.concatenate([[k[0]], k]), 1e-6)
    out = np.column_stack([zz, sp * np.cos(tz), sp * np.sin(tz), thz, kz])
    r = rib(z, u, v, th, a.za, a.zb)
    np.savetxt(os.path.join(a.out, "inflow.txt"), out, fmt="%10.3f %9.4f %9.4f %9.4f %9.5f",
               header=f"z u v theta k: ERF near-1-D column, surface heat flux {a.flux:g} K m/s, at "
                      f"{hdr['time'] / 3600:.2f} h; R_B({a.za:g}-{a.zb:g} m) {r:.4f}; the lidar's direction")
    print(f"flux {a.flux:+g} K m/s: {hdr['time'] / 3600:.2f} h, R_B {r:+.3f}, "
          f"U(300 m) {math.hypot(np.interp(300, z, u), np.interp(300, z, v)):.2f}, theta(10/100 m) "
          f"{np.interp(10, z, th):.2f}/{np.interp(100, z, th):.2f}")


if __name__ == "__main__":
    main()
