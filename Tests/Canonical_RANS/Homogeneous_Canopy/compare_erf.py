#!/usr/bin/env python3
"""ERF on the homogeneous canopy deck against GUST (RT-8.6.4); not a CTest, it needs an ERF build:

    compare_erf.py --erf ERF_EXEC --gust-runs DIR [--run] [--stop 12000]

ERF (a build with the canopy terms of erf-model/ERF#4340) on a 2 x 2 x 200 periodic column of 1 m layers
with the deck's canopy, drag, pressure gradient and wall, the k-eqn RANS closure with the same canopy
terms, a top at 200 m (the deck: 199.7 m), marching in time from rest to --stop seconds (--run; 0.1 s steps). The script reads ERF's
last plotfile in erf_lai2/, erf_lai5/, and GUST's runs lai2/, lai5/ of check_homogeneous_canopy.py in
DIR, and prints u and k at seven heights, and ERF's drift over its last 6000 s. The exit code is 0.
"""
import argparse
import glob
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile as ep  # noqa: E402

DECK = """erf.prob_name = "ABL"
stop_time = {stop}
amrex.fpe_trap_invalid = 0
fabarray.mfiter_tile_size = 1024 1024 1024
geometry.prob_lo = 0 0 0
geometry.prob_hi = 64 64 200
amr.n_cell = 2 2 200
geometry.is_periodic = 1 1 0
amr.max_level = 0
zlo.type = "surface_layer"
erf.most.z0 = 0.01
zhi.type = "SlipWall"
erf.use_gravity = false
erf.use_coriolis = false
erf.molec_diff_type = "None"
erf.les_type = "None"
erf.rans_type = "kEqn"
erf.vert_implicit = true
erf.init_type = "uniform"
prob.rho_0 = 1.2
prob.T_0 = 300.0
prob.U_0 = 4.0
prob.V_0 = 0.0
prob.W_0 = 0.0
prob.KE_0 = 0.1
erf.abl_driver_type = "PressureGradient"
erf.abl_pressure_grad = -0.0015 0. 0.
erf.forest_file = forest.txt
erf.forest_canopy_tke = true
erf.implicit_tke_dissipation = true
erf.fixed_dt = 0.1
erf.fixed_mri_dt_ratio = 4
erf.v = 0
erf.sum_interval = -1
erf.check_int = -1
erf.plot_file_1 = plt
erf.plot_per_1 = 2000
erf.plot_vars_1 = density x_velocity y_velocity z_velocity theta KE
"""
HEIGHTS = (0.25, 0.5, 0.75, 1.0, 1.5, 3.0, 5.0)
HC = 20.0


def erf_column(plt):
    hdr = ep._read_header(plt)
    names = hdr["names"]
    ncomp, _, boxes, fabs = ep._read_cell_h(os.path.join(plt, "Level_0"))
    nz = hdr["hi"][2] - hdr["lo"][2] + 1
    acc = {n: np.zeros(nz) for n in ("x_velocity", "KE")}
    cnt = np.zeros(nz)
    for _box, (fname, off) in zip(boxes, fabs):
        flo, fhi, data = ep._read_fab(os.path.join(plt, "Level_0", fname), off, ncomp)
        n = [fhi[d] - flo[d] + 1 for d in range(3)]
        arr = np.asarray(data).reshape(ncomp, n[2], n[1], n[0])
        for name in acc:
            acc[name][flo[2]:fhi[2] + 1] += arr[names.index(name)].sum(axis=(1, 2))
        cnt[flo[2]:fhi[2] + 1] += n[0] * n[1]
    return acc["x_velocity"] / cnt, acc["KE"] / cnt


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--erf", required=True)
    ap.add_argument("--gust-runs", required=True)
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--stop", type=float, default=12000.0)
    args = ap.parse_args()
    erf = os.path.abspath(args.erf)
    zerf = np.arange(200) + 0.5
    for lai in (2, 5):
        where = "erf_lai%d" % lai
        if args.run:
            os.makedirs(where, exist_ok=True)
            with open(os.path.join(where, "inputs"), "w") as f:
                f.write(DECK.format(stop=args.stop))
            with open(os.path.join(where, "forest.txt"), "w") as f:
                f.write("1 32 32 %g 100000 0.2 %d 0.8\n" % (HC, lai))
            subprocess.run([erf, "inputs"], cwd=where, check=True, stdout=subprocess.DEVNULL)
        plts = sorted((p for p in glob.glob(os.path.join(where, "plt*")) if re.fullmatch(r"plt\d+", os.path.basename(p))),
                      key=lambda s: int(os.path.basename(s)[3:]))
        if len(plts) < 4:
            raise SystemExit("compare_erf: %s has %d plotfiles, needs 4 (--stop at least 6000)" % (where, len(plts)))
        uerf, ke = erf_column(plts[-1])
        drift, _ = erf_column(plts[-4])
        gp = sorted(p for p in glob.glob(os.path.join(args.gust_runs, "lai%d" % lai, "*_plt*")))[-1]
        _, a = ep.read_fields(gp, ["u", "tke", "z"])
        zg = list(a["z"][0][0])
        print("LAI %d, ERF at %s, GUST %s" % (lai, os.path.basename(plts[-1]), os.path.basename(gp)))
        print("   z/h   u ERF  u GUST   diff %   k ERF  k GUST   diff %   ERF drift (m/s, last 3 outputs)")
        for q in HEIGHTS:
            z = q * HC
            u1, u2 = np.interp(z, zerf, uerf), np.interp(z, zg, a["u"][0][0])
            k1, k2 = np.interp(z, zerf, ke), np.interp(z, zg, a["tke"][0][0])
            print("  %4.2f  %6.3f  %6.3f  %+7.1f  %6.3f  %6.3f  %+7.1f  %8.4f"
                  % (q, u1, u2, 100.0 * (u2 / u1 - 1.0), k1, k2, 100.0 * (k2 / k1 - 1.0),
                     abs(np.interp(z, zerf, uerf) - np.interp(z, zerf, drift))))
    return 0


if __name__ == "__main__":
    sys.exit(main())
