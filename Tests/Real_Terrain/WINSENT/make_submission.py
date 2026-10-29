#!/usr/bin/env python3
"""Write one WINSENT Phase I steady-case submission file (NetCDF) from a converged GUST run.

The benchmark asks for the wind at the plateau mast MMNW (this case's origin) at the heights of
template_phase1_inflow.nc, one file per case, named WINSENTbench_PhaseI_<name>_<institution>_LAI<h|l>_STRAT<u|n|s>.nc
(https://winsent-benchmark.readthedocs.io/en/latest/submissions.html). Values:

- u, v: east and north components (the "lat, lon system"); w vertical (m/s);
- WS: horizontal speed; WD: the direction the mean wind blows from (degrees from north);
- TKE: GUST's modelled k (J/kg = m^2/s^2); TI = sqrt(2k/3)/WS (isotropic turbulence, RANS has no resolved std);
- T_pot: GUST's potential temperature at the atmospheric heights. T, pressure, humidity and precipitation are not
  modelled and stay at the template's fill value (NaN).

A steady RANS solution is one representative value per case. Each value is interpolated as sample_stations.py
does: bilinear in x and y between the four surrounding columns, linear in ln(z + z0) between cell centres in height
above each column's ground.
"""
import argparse
import math
import os
import shutil
import sys

import netCDF4
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
import askervein_lines as al  # noqa: E402

FIELDS = ("u", "v", "w", "tke", "theta", "z")


def column_values(cols, terrain, x, y, height, z0):
    """u, v, w, k, theta at (x, y, height above the ground)."""
    h = cols.hdr
    lo, dx = h["prob_lo"], h["dx"]
    fx, fy = (x - lo[0]) / dx[0] - 0.5, (y - lo[1]) / dx[1] - 0.5
    i0, j0 = int(math.floor(fx)), int(math.floor(fy))
    a, b = fx - i0, fy - j0
    out = np.zeros(5)
    for di, dj, wgt in ((0, 0, (1 - a) * (1 - b)), (1, 0, a * (1 - b)), (0, 1, (1 - a) * b), (1, 1, a * b)):
        i, j = i0 + di, j0 + dj
        col = cols.column(i, j)
        q, z = np.array(col[:5]), np.array(col[5])
        zg = al.ground(terrain, lo[0] + (i + 0.5) * dx[0], lo[1] + (j + 0.5) * dx[1])
        agl = z - zg
        if height < agl[0]:
            # below the first centre: the wall law for the velocity, k and theta held at the first centre
            f = math.log((height + z0) / z0) / math.log((agl[0] + z0) / z0)
            v = q[:, 0].copy()
            v[:3] *= f
        else:
            k, _ = al.bracket(list(agl), height)
            t = (math.log(height + z0) - math.log(agl[k] + z0)) / (math.log(agl[k + 1] + z0) - math.log(agl[k] + z0))
            v = (1 - t) * q[:, k] + t * q[:, k + 1]
        out += wgt * v
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("plotfile", help="the converged run's plotfile, or its prefix (highest step)")
    p.add_argument("--case", required=True, help="folder with terrain.txt and stations.txt")
    p.add_argument("--template", required=True, help="template_phase1_inflow.nc")
    p.add_argument("--lai", required=True, choices=("h", "l"))
    p.add_argument("--strat", dest="stability", required=True, choices=("u", "n", "s"))
    p.add_argument("--name", default="GUST")
    p.add_argument("--institution", default="NLR")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--wd", type=float, default=293.0, help="the frame's wind direction (deg from north)")
    p.add_argument("--z0", type=float, default=0.1)
    a = p.parse_args()

    cols = al.Columns(a.plotfile)
    names = cols.hdr["names"]
    missing = [f for f in FIELDS if f not in names]
    if missing:
        sys.exit(f"{cols.path}: no field(s) {missing}; give the GUST plotfile (u v w tke theta z)")
    cols.comps = [names.index(f) for f in FIELDS]
    terrain = al.read_terrain(os.path.join(a.case, "terrain.txt"))
    stations = {ln.split()[0]: dict(s.split("=") for s in ln.split()[1:])
                for ln in open(os.path.join(a.case, "stations.txt")) if ln.strip()}
    x, y = float(stations["MMNW"]["x"]), float(stations["MMNW"]["y"])

    az = math.radians(a.wd + 180.0)
    ex = (math.sin(az), math.cos(az))  # the frame's x and y axes in (east, north)
    ey = (-ex[1], ex[0])

    os.makedirs(a.out_dir, exist_ok=True)
    out = os.path.join(a.out_dir, f"WINSENTbench_PhaseI_{a.name}_{a.institution}_LAI{a.lai}_STRAT{a.stability}.nc")
    shutil.copy(a.template, out)
    with netCDF4.Dataset(out, "r+") as nc:
        hw = np.asarray(nc["height_wind"][:], dtype=float)
        ha = np.asarray(nc["height_atm"][:], dtype=float)
        rows = []
        for i, h in enumerate(hw):
            u, v, w, k, _ = column_values(cols, terrain, x, y, h, a.z0)
            east, north = u * ex[0] + v * ey[0], u * ex[1] + v * ey[1]
            ws = math.hypot(east, north)
            wd = math.degrees(math.atan2(-east, -north)) % 360.0
            ti = math.sqrt(2.0 * max(k, 0.0) / 3.0) / ws if ws > 0 else float("nan")
            for var, val in (("u", east), ("v", north), ("w", w), ("WS", ws), ("WD", wd), ("TKE", k), ("TI", ti)):
                nc[var][i] = val
            rows.append((h, east, north, w, ws, wd, k, ti))
        for i, h in enumerate(ha):
            nc["T_pot"][i] = column_values(cols, terrain, x, y, h, a.z0)[4]
        nc.title = f"WINSENT benchmark Phase I, steady case LAI {a.lai}, stratification {a.stability}: GUST steady RANS"
        nc.source = f"GUST (erf-model/gust-abl), plotfile {os.path.basename(cols.path)}"
        nc.location = "MMNW"
        nc.averaging = "steady RANS solution: one representative value per case (no time averaging)"
        nc.TI_definition = "sqrt(2 TKE / 3) / WS (modelled TKE, isotropic)"
    print(out)
    print("  z(m)   u_E    v_N      w     WS     WD    TKE    TI")
    for r in rows:
        print("  %4.0f %6.2f %6.2f %6.3f %6.2f %6.1f %6.3f %5.3f" % r)


if __name__ == "__main__":
    main()
