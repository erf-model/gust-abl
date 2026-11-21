"""compare_fields.py: an immersed run against its terrain-following (TF) partner, column by column (README.md).

    python3 compare_fields.py --body-fitted DIR --immersed DIR --flat H --wind-dir DEG [--plot FILE]

Each folder holds one run's plotfiles (plt*, turb*, with or without a prefix; the last of each is read). In every
column the wind speed S and k are read at 10 and 40 m above that column's ground: TF, the lowest node; immersed,
the solid depth (the sum of building_beta times each layer's thickness, the layers taken from the cell centres)
minus --flat (terrain.flat_height), with the cut cell at the centroid of its fluid and solid cells left out. The
speed-up dS = S / S_up - 1 uses each run's own mean over the upwind sixteenth of the domain (the band of columns
furthest upwind along the wind direction, --wind-dir, the direction the wind comes from).

Printed: the ground's largest slope along the wind and rise per column, then per height the hilltop dS and k of
both runs and, over regions (windward slope > 0.15, upper hill, lee slope < -0.15, flat ground), the mean and
lowest immersed minus TF dS and the mean and largest immersed / TF k ratio, and the correlation of the gap with
the slope along the wind on the hill (ground above 0.15 of the top).
"""
import argparse
import glob
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402

AGL = (10.0, 40.0)


def read(d, imm):
    last = lambda pat: sorted(glob.glob(os.path.join(d, pat)))[-1]  # noqa: E731
    _, f = erf_plotfile.read_fields(last("*plt?????"), ["u", "v", "z"] + (["building_beta"] if imm else []))
    f = {k: np.asarray(v) for k, v in f.items()}
    _, t = erf_plotfile.read_fields(last("*turb?????"), ["KE"])
    f["tke"] = np.asarray(t["KE"])
    return f


def fields(d, imm, flat):
    f = read(d, imm)
    nx, ny, _ = f["u"].shape
    S = np.zeros((len(AGL), nx, ny))
    K = np.zeros_like(S)
    G = np.zeros((nx, ny))
    spd = np.hypot(f["u"], f["v"])
    for i in range(nx):
        for j in range(ny):
            z, s, k = f["z"][i, j], spd[i, j], f["tke"][i, j]
            if imm:
                b = f["building_beta"][i, j]
                faces = [0.0]
                for c in z:
                    faces.append(2 * c - faces[-1])
                faces = np.array(faces)
                g = float(np.sum(b * np.diff(faces)))
                zc = np.where(faces[:-1] < g, 0.5 * (g + faces[1:]), z)
                live = b < 1
                z, s, k = zc[live] - g, s[live], k[live]
                G[i, j] = g - flat
            else:
                g = z[0] - 0.5 * (z[1] - z[0])
                z = z - g
                G[i, j] = g
            S[:, i, j] = np.interp(AGL, z, s)
            K[:, i, j] = np.interp(AGL, z, k)
    return S, K, G


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--body-fitted", required=True)
    p.add_argument("--immersed", required=True)
    p.add_argument("--flat", type=float, required=True)
    p.add_argument("--wind-dir", type=float, required=True)
    p.add_argument("--plot")
    a = p.parse_args()
    for d in (a.body_fitted, a.immersed):
        if not glob.glob(os.path.join(d, "*plt?????")) or not glob.glob(os.path.join(d, "*turb?????")):
            sys.exit(f"compare_fields.py: no plt?????/turb????? plotfiles in {d}")
    St, Kt, Gt = fields(a.body_fitted, False, 0.0)
    Si, Ki, Gi = fields(a.immersed, True, a.flat)
    if Gt.shape != Gi.shape:
        sys.exit(f"compare_fields.py: the runs have {Gt.shape} and {Gi.shape} columns; they must match")
    h = erf_plotfile._read_header(sorted(glob.glob(os.path.join(a.body_fitted, "*plt?????")))[-1])
    lo, dx = h["prob_lo"], h["dx"]
    nx, ny = Gt.shape
    X, Y = np.meshgrid(lo[0] + (np.arange(nx) + 0.5) * dx[0], lo[1] + (np.arange(ny) + 0.5) * dx[1], indexing="ij")
    wd = np.deg2rad(a.wind_dir)
    ew = -np.array([np.sin(wd), np.cos(wd)])  # the direction the wind blows toward
    along = X * ew[0] + Y * ew[1]
    up = along < np.quantile(along, 1 / 16)
    dSt = St / St[:, up].mean(1)[:, None, None] - 1
    dSi = Si / Si[:, up].mean(1)[:, None, None] - 1
    if nx < 2 or ny < 2:
        sys.exit("compare_fields.py: needs at least 2 columns in x and in y")
    gx, gy = np.gradient(Gt, dx[0], dx[1])
    slope = gx * ew[0] + gy * ew[1]
    print(f"ground: TF top {Gt.max():.1f} m, immersed top {Gi.max():.1f} m; largest |slope| along the wind "
          f"{abs(slope).max():.2f}, rise per column {abs(slope).max() * max(dx[0], dx[1]):.1f} m")
    it = np.unravel_index(np.argmax(Gt), Gt.shape)
    top, hill = Gt.max(), Gt > 0.15 * Gt.max()
    regions = (("windward slope > 0.15", hill & (slope > 0.15)), ("upper hill", Gt > 0.65 * top),
               ("lee slope < -0.15", hill & (slope < -0.15)), ("flat ground", Gt < 0.04 * top))
    for n, z in enumerate(AGL):
        print(f"\n{int(z)} m above the ground: hilltop (x {X[it]:.0f}, y {Y[it]:.0f}) dS TF {dSt[n][it]:.3f}, immersed "
              f"{dSi[n][it]:.3f}, gap {dSi[n][it] - dSt[n][it]:+.3f}; k TF {Kt[n][it]:.2f}, immersed {Ki[n][it]:.2f}")
        for lab, m in regions:
            if not m.any():
                continue
            gap, kr = (dSi[n] - dSt[n])[m], (Ki[n] / Kt[n])[m]
            print(f"  {lab:22} {m.sum():6d} columns: gap mean {gap.mean():+.3f}, lowest {gap.min():+.3f}; "
                  f"k ratio mean {kr.mean():.2f}, largest {kr.max():.2f}")
        c = np.corrcoef(slope[hill], (dSi[n] - dSt[n])[hill])[0, 1]
        print(f"  correlation of the slope along the wind with the gap, on the hill: {c:+.2f}")
    if a.plot:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(1, 2, figsize=(12, 5), constrained_layout=True)
        for x_, val, cm, lab, lim in ((ax[0], dSi[0] - dSt[0], "RdBu", "dS immersed minus TF, 10 m", (-0.4, 0.4)),
                                      (ax[1], Ki[0] / Kt[0], "PuOr_r", "k immersed / TF, 10 m", (0.5, 1.5))):
            pc = x_.pcolormesh(X, Y, val, cmap=cm, vmin=lim[0], vmax=lim[1], shading="nearest")
            x_.contour(X, Y, Gt, levels=np.linspace(0.15, 0.95, 5) * top, colors="k", linewidths=0.5)
            x_.set_aspect("equal")
            x_.set_title(lab, fontsize=10, loc="left")
            fig.colorbar(pc, ax=x_, shrink=0.8)
        fig.savefig(a.plot, dpi=140)
        print("wrote", a.plot)


if __name__ == "__main__":
    main()
