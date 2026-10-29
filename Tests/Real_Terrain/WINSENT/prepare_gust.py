#!/usr/bin/env python3
"""Turn one WINSENT case folder prepared for ERF into GUST's input tables.

The ERF folder comes from ERF's Exec/CanonicalTests/Real_Terrain/WINSENT/prepare_winsent.py (terrain.txt,
forest_lai.nc, forest_height.nc, inflow.txt, zlevels.txt, stations.txt). This script writes, in --out:

- terrain.txt: ERF's structured terrain file, copied as it is (GUST reads it with terrain.type = erf_grid);
- canopy_lai.xyz, canopy_height.xyz: the forest maps as ``x y value`` points (canopy.type = map);
- inflow.txt: ``z u v theta k`` rows (bc.inflow_file); ERF's columns are ``z u v w theta`` and its inflow
  k is the constant xlo.KE, given here with --inflow-k;
- zlevels.txt: the levels of ERF's erf.terrain_z_levels, one per line (mesh.zlevels_file);
- stations.txt: copied, for the sampling after the run.

With --roughness, it also writes roughness.xyz (bc.surface_z0_file) from the roughness map of the kynema-sgf
WINSENT runs (winsent.roughness: ERF's structured-grid format, 10 m, in kynema-sgf's frame). That frame is this
case's frame turned by --frame-rotation degrees (counter-clockwise) and shifted by --frame-shift; the defaults are
fitted by matching the two terrain files (6.8 m RMS height difference over the masts' area). Each point is a
32 m cell centre of the forest maps, with z0 the mean of ln z0 over the 10 m pixels inside the cell.
"""
import argparse
import math
import os
import shutil
import sys

import netCDF4
import numpy as np


def write_xyz(nc_path, var, out_path, uniform=None):
    """Write the grid as x y value points; with uniform, every point whose value is > 0 takes that value instead."""
    with netCDF4.Dataset(nc_path) as nc:
        x = np.asarray(nc.variables["x"][:], dtype=float)
        y = np.asarray(nc.variables["y"][:], dtype=float)
        v = np.asarray(nc.variables[var][:], dtype=float)
    if v.shape != (len(y), len(x)):
        sys.exit(f"{nc_path}: {var} has shape {v.shape}, expected (y, x) = ({len(y)}, {len(x)})")
    if not np.all(np.isfinite(v)) or v.min() < 0:
        sys.exit(f"{nc_path}: {var} must be finite and >= 0")
    if uniform is not None:
        v = np.where(v > 0, uniform, 0.0)
    X, Y = np.meshgrid(x, y)
    np.savetxt(out_path, np.column_stack([X.ravel(), Y.ravel(), v.ravel()]), fmt="%.3f %.3f %.6g")
    return v


def write_roughness(map_path, nc_path, rotation, shift, out_path):
    vals = np.loadtxt(map_path)
    nx, ny = int(vals[0]), int(vals[1])
    xk, yk = vals[2:2 + nx], vals[2 + nx:2 + nx + ny]
    z0 = vals[2 + nx + ny:].reshape(nx, ny)  # y fastest, as ERF's terrain reader
    if not np.all(np.isfinite(z0)) or z0.min() <= 0:
        sys.exit(f"{map_path}: every z0 must be finite and > 0")
    dk = xk[1] - xk[0]
    with netCDF4.Dataset(nc_path) as nc:
        x = np.asarray(nc.variables["x"][:], dtype=float)
        y = np.asarray(nc.variables["y"][:], dtype=float)
    X, Y = np.meshgrid(x, y)
    c, s = math.cos(math.radians(rotation)), math.sin(math.radians(rotation))
    half = 0.5 * (x[1] - x[0])
    offs = np.arange(-half + 0.5 * dk, half, dk)  # the map's pixels inside one cell
    acc = np.zeros_like(X)
    for ox in offs:
        for oy in offs:
            px, py = X + ox, Y + oy
            xr, yr = c * px - s * py + shift[0], s * px + c * py + shift[1]
            i = np.round((xr - xk[0]) / dk).astype(int)
            j = np.round((yr - yk[0]) / dk).astype(int)
            if i.min() < 0 or i.max() >= nx or j.min() < 0 or j.max() >= ny:
                sys.exit(f"{map_path}: the case's domain is not inside the map with this frame")
            acc += np.log(z0[i, j])
    zc = np.exp(acc / len(offs) ** 2)
    np.savetxt(out_path, np.column_stack([X.ravel(), Y.ravel(), zc.ravel()]), fmt="%.3f %.3f %.5g")
    return zc


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--erf-case", required=True, help="ERF case folder (blended_high_neutral, ...)")
    p.add_argument("--out", required=True, help="GUST run folder to write")
    p.add_argument("--inflow-k", type=float, default=0.2, help="inflow k (m^2/s^2), ERF's xlo.KE (default 0.2)")
    p.add_argument("--roughness", help="kynema-sgf's winsent.roughness, to write roughness.xyz")
    p.add_argument("--frame-rotation", type=float, default=-1.0,
                   help="kynema-sgf's frame against this case's, degrees counter-clockwise (default -1: 292 vs 293)")
    p.add_argument("--frame-shift", type=float, nargs=2, default=[1410.0, -400.0],
                   help="this case's origin (MMNW) in kynema-sgf's frame, m (default 1410 -400)")
    p.add_argument("--canopy-height", type=float,
                   help="a uniform tree height (m) on every forested point, in place of the height map "
                        "(the runs of the README use 15, a fitted height)")
    a = p.parse_args()
    if a.canopy_height is not None and not (math.isfinite(a.canopy_height) and a.canopy_height > 0):
        sys.exit("--canopy-height must be finite and > 0 (m)")
    src = a.erf_case
    os.makedirs(a.out, exist_ok=True)

    for f in ("terrain.txt", "stations.txt", "summary.txt"):
        shutil.copy(os.path.join(src, f), os.path.join(a.out, f))

    lai = write_xyz(os.path.join(src, "forest_lai.nc"), "LAI", os.path.join(a.out, "canopy_lai.xyz"))
    h = write_xyz(os.path.join(src, "forest_height.nc"), "height", os.path.join(a.out, "canopy_height.xyz"),
                  a.canopy_height)

    rows = np.loadtxt(os.path.join(src, "inflow.txt"), comments="#")
    if rows.ndim != 2 or rows.shape[1] != 5:
        sys.exit("inflow.txt: expected the columns z u v w theta")
    if np.abs(rows[:, 3]).max() > 0:
        sys.exit("inflow.txt: w is not 0; GUST's inflow file has no w column")
    k = np.full(len(rows), a.inflow_k)
    np.savetxt(os.path.join(a.out, "inflow.txt"), np.column_stack([rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 4], k]),
               fmt="%10.3f %9.4f %9.4f %8.2f %8.4f",
               header="z (m above ground)  u  v  theta  k   (from ERF's inflow.txt; k = ERF's xlo.KE)")

    levels = None
    for line in open(os.path.join(src, "zlevels.txt")):
        if line.strip().startswith("erf.terrain_z_levels"):
            levels = [float(t) for t in line.split("=", 1)[1].split()]
    if not levels:
        sys.exit("zlevels.txt: no erf.terrain_z_levels line")
    np.savetxt(os.path.join(a.out, "zlevels.txt"), np.array(levels), fmt="%.2f")

    print(f"{a.out}: {len(levels) - 1} levels to {levels[-1]:.0f} m; inflow {len(rows)} rows, k {a.inflow_k}; "
          f"LAI max {lai.max():.2f}, canopy height max {h.max():.1f} m, forested fraction {np.mean(lai > 0):.3f}")
    if a.roughness:
        z0 = write_roughness(a.roughness, os.path.join(src, "forest_lai.nc"), a.frame_rotation, a.frame_shift,
                             os.path.join(a.out, "roughness.xyz"))
        print(f"roughness.xyz: z0 {z0.min():.4g} to {z0.max():.4g} m, median {np.median(z0):.3g}; "
              f"under the forest median {np.median(z0[lai > 0]):.3g}")


if __name__ == "__main__":
    main()
