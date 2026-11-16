#!/usr/bin/env python3
"""The Terrain_Shapes_3D case: analytic terrain shapes and the probe lattice.

    make_terrain.py [OUT_DIR] [--dx DX] [--only NAME[,NAME...]]

Writes <name>.xyz for each shape and probes.txt into OUT_DIR (default: the current directory), and prints, for each
shape, the highest ground and the largest step between two nodes, then the layer ratios for mesh.dz0 = 2 and 8 m and
the log-law speed of the lid.

A body is h = H/2 (1 - tanh(d / W)), H = 760 m: d is the distance outside an ellipse with semi-axes ax and ay (to
first order: (r_e - 1) / |grad r_e| with r_e the normalised radius), turned by rot degrees; W is Ww on the windward
side and Wl on the lee side, blended with the cosine of the angle to the wind (the wind is along +x). Several bodies
combine by their maximum. A crater of depth Dc and radius Rc takes Dc/2 (1 - tanh((r - Rc) / Wc)) off the first body.
The domain is 16 km x 12.8 km (y periodic) on DX = 100 m cells (--dx changes it; the check uses 200 m).
"""
import math
import os
import sys

LX, LY, TOP, NZ = 16000.0, 12800.0, 6000.0, 40
H = 760.0
XC, YC = 6000.0, 6400.0
USTAR, Z0, KAPPA = 0.3031, 0.0002, 0.41
INF = 1.0e7

# name: (bodies, crater); a body is (xc, yc, ax, ay, Ww, Wl, rot), the crater (depth, radius, width) or None
SHAPES = {
    "round_gentle":   ([(XC, YC, 2000.0, 2000.0, 600.0, 600.0, 0.0)], None),
    "round_mid":      ([(XC, YC, 2000.0, 2000.0, 300.0, 300.0, 0.0)], None),
    "round_steep":    ([(XC, YC, 2000.0, 2000.0, 150.0, 150.0, 0.0)], None),
    "round_sharp":    ([(XC, YC, 2000.0, 2000.0, 90.0, 90.0, 0.0)], None),
    "steep_lee":      ([(XC, YC, 2000.0, 2000.0, 600.0, 150.0, 0.0)], None),
    "steep_windward": ([(XC, YC, 2000.0, 2000.0, 150.0, 600.0, 0.0)], None),
    "ridge_cross":    ([(XC, YC, 800.0, 3000.0, 150.0, 150.0, 0.0)], None),
    "ridge_oblique":  ([(XC, YC, 800.0, 3000.0, 150.0, 150.0, 30.0)], None),
    "ridge_stream":   ([(XC, YC, 3000.0, 800.0, 150.0, 150.0, 0.0)], None),
    "two_gap":        ([(XC, YC - 2200.0, 1500.0, 1500.0, 150.0, 150.0, 0.0),
                        (XC, YC + 2200.0, 1500.0, 1500.0, 150.0, 150.0, 0.0)], None),
    "butte_canyon":   ([(XC - 1000.0, YC, 2000.0, 2000.0, 150.0, 150.0, 0.0),
                        (XC + 3200.0, YC, 600.0, 3000.0, 150.0, 150.0, 0.0)], None),
    "narrow":         ([(XC, YC, 800.0, 800.0, 150.0, 150.0, 0.0)], None),
    "rim":            ([(XC, YC, 2000.0, 2000.0, 150.0, 150.0, 0.0)], (400.0, 1200.0, 150.0)),
    "ridge_2d":       ([(XC, YC, 1000.0, INF, 150.0, 150.0, 0.0)], None),
}


def body(x, y, xc, yc, ax, ay, ww, wl, rot):
    c, s = math.cos(math.radians(rot)), math.sin(math.radians(rot))
    dx, dy = x - xc, y - yc
    xr, yr = c * dx + s * dy, -s * dx + c * dy
    re = math.sqrt((xr / ax) ** 2 + (yr / ay) ** 2) + 1e-12
    grad = math.sqrt((xr / ax ** 2) ** 2 + (yr / ay ** 2) ** 2) / re + 1e-12
    d = (re - 1.0) / max(grad, 0.25 / max(ax, ay))
    cw = dx / max(math.hypot(dx, dy), 1e-9)          # +1 on the lee side, -1 on the windward side
    w = ww + (wl - ww) * 0.5 * (1.0 + cw)
    return 0.5 * H * (1.0 - math.tanh(d / w))


def height(name, x, y):
    bodies, crater = SHAPES[name]
    h = max(body(x, y, *b) for b in bodies)
    if crater:
        depth, rc, wc = crater
        h -= 0.5 * depth * (1.0 - math.tanh((math.hypot(x - bodies[0][0], y - bodies[0][1]) - rc) / wc))
    return h


def ratio(dz0):
    """The layer ratio that puts NZ layers from dz0 under TOP."""
    lo, hi = 1.0 + 1e-9, 2.0
    for _ in range(200):
        r = 0.5 * (lo + hi)
        lo, hi = (r, hi) if dz0 * (r ** NZ - 1.0) / (r - 1.0) < TOP else (lo, r)
    return 0.5 * (lo + hi)


def lid_speed():
    return USTAR / KAPPA * math.log(TOP / Z0)


def probes(dx):
    """The lattice: every 400 m, at 80 m above the ground and then at 20 m."""
    every = max(1, int(round(400.0 / dx)))
    nx, ny = int(round(LX / dx)), int(round(LY / dx))
    pts = []
    for agl in (80.0, 20.0):
        for j in range(every // 2, ny, every):
            for i in range(every // 2, nx, every):
                pts.append(((i + 0.5) * dx, (j + 0.5) * dx, agl))
    return pts


def write(out, dx=100.0, only=None):
    nx, ny = int(round(LX / dx)), int(round(LY / dx))
    xs = [i * dx for i in range(nx + 1)]
    ys = [j * dx for j in range(-1, ny + 2)]
    report = {}
    for name in SHAPES:
        if only and name not in only:
            continue
        hs = [[height(name, x, y) for x in xs] for y in ys]
        with open(os.path.join(out, name + ".xyz"), "w") as fh:
            fh.write("# %s\n" % name)
            for y, row in zip(ys, hs):
                for x, h in zip(xs, row):
                    fh.write("%.1f %.1f %.3f\n" % (x, y, h))
        step = max(max(abs(b - a) for a, b in zip(row[:-1], row[1:])) for row in hs)
        step = max(step, max(abs(b - a) for ra, rb in zip(hs[:-1], hs[1:]) for a, b in zip(ra, rb)))
        report[name] = (max(max(row) for row in hs), step)
    pts = probes(dx)
    with open(os.path.join(out, "probes.txt"), "w") as fh:
        fh.write("# lattice every 400 m: probes 1-%d at 80 m above the ground, %d-%d at 20 m\n"
                 % (len(pts) // 2, len(pts) // 2 + 1, len(pts)))
        for p in pts:
            fh.write("%.1f %.1f %g\n" % p)
    return report


def main(argv):
    args, dx, only = [], 100.0, None
    it = iter(argv[1:])
    for a in it:
        if a == "--dx":
            dx = float(next(it))
        elif a == "--only":
            only = next(it).split(",")
        elif a in ("-h", "--help"):
            print(__doc__)
            return 0
        else:
            args.append(a)
    out = args[0] if args else "."
    os.makedirs(out, exist_ok=True)
    for name, (top, step) in write(out, dx, only).items():
        print("%-16s highest ground %.0f m, largest node step %.0f m (%.0f deg)"
              % (name, top, step, math.degrees(math.atan(step / dx))))
    for dz0 in (2.0, 8.0):
        print("mesh.dz0 = %g: mesh.stretch_ratio = %.6f" % (dz0, ratio(dz0)))
    print("lid speed %.3f m/s (bc.zhi_wall_velocity)" % lid_speed())
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
