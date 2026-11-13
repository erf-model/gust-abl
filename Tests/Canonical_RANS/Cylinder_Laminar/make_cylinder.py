#!/usr/bin/env python3
"""Write the inputs of the laminar cylinder case (README.md): the staircase cylinder as a
buildings.format = xyz map and the uniform inflow as a bc.inflow_file.

The cylinder lies along y (periodic, two cells thick): diameter D = 1 m, axis at (x, z) = (0, zc),
zc = 6 - d with d = 1/n the cell size: a cell face, so the staircase is symmetric about the axis
and spans the full diameter (an even number of cells; an axis on a cell centre gives an odd number
and a frontal height of D - d), and the axis lies one cell closer to the ground than to the top
(12 m). That asymmetry, d in 12 m, is the seed of the vortex shedding at Re = 100; the published
values are for an axis midway between far walls.
The flow plane is x-z. A building map gives each column one solid segment from zstart to
height (docs/buildings.rst), which is what a circle cut by vertical columns is: in the column
of cell centre x_i the solid cells are those whose centre (x_i, z_k) lies inside the circle,
and the segment runs from the lower face of the lowest to the upper face of the highest. Every
cell is then solid or fluid, and the wall is a staircase. The map holds one point per cell
centre of [-1, 1] in x (height 0 where the column misses the circle) and per y cell centre of the
decks' two cells in y (0.1 m in all).

  python3 make_cylinder.py --cells-per-d 16 --out cylinder_16.txt
  python3 make_cylinder.py --inflow inflow_uniform.txt
"""
import argparse
import math

X0 = -6.0     # domain's lower x (m); the decks' mesh.prob_lo
YS = (0.025, 0.075)   # the two cell centres in y: the decks' 2 cells, 0.1 m (a map needs 2 y values)
H = 12.0      # the domain height (m), mesh.prob_hi z


def axis_height(n):
    """zc = H/2 - d (m): a cell face one cell below mid-height."""
    return 0.5 * H - 1.0 / n


def cylinder_map(n, x0=X0):
    """Lines 'x y height zstart' and the solid area (m^2) of the staircase."""
    d = 1.0 / n
    zc = axis_height(n)
    i0 = math.floor((-1.0 - x0) / d)
    i1 = math.ceil((1.0 - x0) / d)
    lines = []
    area = 0.0
    for i in range(i0, i1):
        x = x0 + (i + 0.5) * d
        # Cell k spans [k d, (k + 1) d]; its centre is inside when (x^2 + (zk - zc)^2) < 1/4.
        ks = [k for k in range(int((zc - 1) / d), int((zc + 1) / d) + 2)
              if x * x + ((k + 0.5) * d - zc) ** 2 < 0.25]
        if ks:
            lo, hi = min(ks) * d, (max(ks) + 1) * d
            area += (hi - lo) * d
            lines += [f"{x:.12g} {y} {hi:.12g} {lo:.12g}" for y in YS]
        else:
            lines += [f"{x:.12g} {y} 0 0" for y in YS]
    return lines, area


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cells-per-d", type=int, help="cells per diameter (map)")
    ap.add_argument("--out", help="map file")
    ap.add_argument("--inflow", help="inflow file (u = 1 m/s, v = 0, theta = 300 K, k = 0)")
    ap.add_argument("--top", type=float, default=H, help="domain top (m) for the inflow file")
    a = ap.parse_args()
    if a.out:
        if not a.cells_per_d or a.cells_per_d < 2:
            ap.error("--cells-per-d >= 2 is required with --out")
        lines, area = cylinder_map(a.cells_per_d)
        with open(a.out, "w") as f:
            f.write(f"# Staircase cylinder D = 1 m along y at (x, z) = (0, {axis_height(a.cells_per_d):.12g}), {a.cells_per_d} cells per D "
                    f"(make_cylinder.py): x y height zstart; solid area {area:.6f} m^2 (circle {math.pi / 4:.6f})\n")
            f.write("\n".join(lines) + "\n")
        print(f"{a.out}: {len(lines)} points, solid area {area:.6f} m^2 (pi/4 = {math.pi / 4:.6f})")
    if a.inflow:
        with open(a.inflow, "w") as f:
            f.write("# Uniform inflow of the laminar cylinder (make_cylinder.py): z u v theta k\n")
            f.write(f"0 1 0 300 0\n{a.top:.12g} 1 0 300 0\n")
        print(f"{a.inflow}: u = 1 m/s up to z = {a.top:g} m")


if __name__ == "__main__":
    main()
