#!/usr/bin/env python3
"""The flow around one cube in a GUST plotfile.

    check_cube.py PLOTFILE
    check_cube.py --latest PREFIX

Cell values, 5 m cells: the velocity is 0 in every solid cell; the sum of u over each x-plane is the
first plane's within --flux-tol (cell means, so not exact where a plane cuts the block's faces); the
air reverses behind the block and speeds up over its roof; u in a wall cell of the roof and of a side
and k in the roof's wall cell are pinned to this deck's (--roof, --side, --k-roof: the defaults are those
of the deck as it is; buildings.wall_tke = off has 3.833, 2.387 and 8.618, buildings.wall_form =
roof_stress 3.888, 1.505 and 8.589, bc.surface_wall_k = source 2.843, 1.599 and 0.911). k in a side's wall cell of
the ground's first layer is pinned too (--k-ground: 0.0151; 0.0099, 0.0089 and 0.3124 for those three). With buildings.wall_tke = on k in the solid cells is the closure's floor
(--k-solid, the largest allowed; the two other forms solve k there and pass a large value). The exit code
is the verdict.
"""
import argparse
import glob
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("--latest", action="store_true", help="plotfile is a prefix: the last plotfile written with it")
    ap.add_argument("--flux-tol", type=float, default=0.02)
    ap.add_argument("--roof", type=float, default=2.801, help="u in the roof's wall cell (m/s)")
    ap.add_argument("--side", type=float, default=1.491, help="u in a side's wall cell (m/s)")
    ap.add_argument("--k-roof", type=float, default=0.411, help="k in the roof's wall cell (m2/s2), within 2 %%")
    ap.add_argument("--k-ground", type=float, default=0.0151, help="k in a side's wall cell of the first layer (m2/s2), within 2 %%")
    ap.add_argument("--k-solid", type=float, default=1.001e-6, help="largest k allowed in a solid cell (m2/s2); the deck's turb.tke_floor is 1e-6")
    ap.add_argument("--theta-solid", type=float, default=None,
                    help="theta every solid cell must hold (K, within 1e-9): the walls' temperature with buildings.theta_mode = ambient or fixed")
    ap.add_argument("--theta-solid-free", type=float, default=None,
                    help="theta the solid cells must NOT all hold (K; some cell off by more than 1e-6): with "
                         "buildings.theta_mode = none theta is solved inside the block, as before its building terms")
    a = ap.parse_args()
    name = sorted(glob.glob(a.plotfile + "[0-9]*"))[-1] if a.latest else a.plotfile
    hdr, d = erf_plotfile.read_fields(name, ["u", "v", "w", "tke", "theta", "building_beta"])
    u, v, w, beta = d["u"], d["v"], d["w"], d["building_beta"]
    nx, ny, nz = [hdr["hi"][n] - hdr["lo"][n] + 1 for n in range(3)]
    cells = [(i, j, k) for i in range(nx) for j in range(ny) for k in range(nz)]
    bad = []
    top = max(abs(f[i][j][k]) for f in (u, v, w) for i, j, k in cells if beta[i][j][k] >= 1)
    print(f"largest speed in a solid cell: {top:.3e} m/s")
    if top > 1e-10:
        bad.append("the solid cells move")
    flux = [sum(u[i][j][k] for j in range(ny) for k in range(nz)) for i in range(nx)]
    print(f"plane flux / first plane: min {min(flux) / flux[0]:.4f}, max {max(flux) / flux[0]:.4f}")
    if max(abs(f / flux[0] - 1) for f in flux) > a.flux_tol:
        bad.append("the volume flux changes along x")
    jc = ny // 2
    behind = u[17][jc][0]
    over = u[14][jc][8]
    far = u[1][jc][8]
    print(f"u behind the cube {behind:.3f}, 20 m over its roof {over:.3f}, upstream at that height {far:.3f} m/s")
    # The roof's wall cell against this deck's value: moves with the wall law and the closed face.
    wall = u[14][jc][4]
    print(f"u in the roof's wall cell {wall:.3f} m/s (pinned at {a.roof})")
    if abs(wall - a.roof) > 0.02:
        bad.append("the roof's wall cell moved")
    # A wall cell of the side at y = 45 m, second layer: moves with the sides' wall law and their closed faces.
    side = u[14][8][1]
    print(f"u in a side's wall cell {side:.3f} m/s (pinned at {a.side})")
    if abs(side - a.side) > 0.02:
        bad.append("the side's wall cell moved")
    # k in the roof's wall cell: its wall value with buildings.wall_tke = on, the k equation's without.
    k_roof = d["tke"][14][jc][4]
    print(f"k in the roof's wall cell {k_roof:.3f} m2/s2 (pinned at {a.k_roof}); largest k {max(d['tke'][i][j][k] for i, j, k in cells):.2f}")
    if abs(k_roof - a.k_roof) > 0.02 * a.k_roof:
        bad.append("k in the roof's wall cell moved")
    # k in a side's wall cell of the ground's first layer: the ground's wall value and the side's together.
    k_ground = d["tke"][14][8][0]
    print(f"k in a side's wall cell on the ground {k_ground:.4f} m2/s2 (pinned at {a.k_ground})")
    if abs(k_ground - a.k_ground) > 0.02 * a.k_ground:
        bad.append("k in the side's wall cell on the ground moved")
    k_solid = max(d["tke"][i][j][k] for i, j, k in cells if beta[i][j][k] >= 1)
    print(f"largest k in a solid cell {k_solid:.3e} m2/s2 (at most {a.k_solid:.3e})")
    if k_solid > a.k_solid:
        bad.append("k in the solid cells is above its floor")
    # theta in the solid cells: held at the walls' temperature, or solved as before theta's building terms.
    off = max(abs(d["theta"][i][j][k] - (a.theta_solid if a.theta_solid is not None else (a.theta_solid_free or 0.0)))
              for i, j, k in cells if beta[i][j][k] >= 1)
    if a.theta_solid is not None:
        print(f"largest |theta - {a.theta_solid}| in a solid cell {off:.3e} K (held: at most 1e-9)")
        if off > 1e-9:
            bad.append("theta in the solid cells is not the walls' temperature")
    if a.theta_solid_free is not None:
        print(f"largest |theta - {a.theta_solid_free}| in a solid cell {off:.3e} K (solved: more than 1e-6)")
        if off <= 1e-6:
            bad.append("theta in the solid cells is held, not solved")
    if behind >= 0:
        bad.append("no reversal behind the cube")
    if over <= far:
        bad.append("no speed-up over the roof")
    for b in bad:
        print("FAIL:", b)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
