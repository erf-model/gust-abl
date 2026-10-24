#!/usr/bin/env python3
"""Profile drift in the empty-domain ABL test (Section 5.1 of the paper), from a GUST turbulence plotfile.

    python3 -I hhabl_drift.py <turb plotfile> <reader dir> --law power_law --speed 8.7502 --zref 10 --alpha 0.2113
    python3 -I hhabl_drift.py <turb plotfile> <reader dir> --law log_law --ustar 0.61726

<reader dir> is the directory of ERF's erf_plotfile.py (extern_erf/).

At each station x (default: the first cell, the hill position x = 0 of the Askervein deck, and the last cell) the
script averages the columns over y (the deck is periodic in y), interpolates linearly in height, and prints the
wind speed U, k and nu_t against the inflow law (the analytic profile of Source/InflowProfile.H) and against the
first column. It also prints the ground friction velocity along x against the law's u* at z0. With
--max-u-drift it exits non-zero when |U/U_law - 1| at any requested height and station exceeds the bound.

Standard library only, like the regression check scripts.
"""
import argparse
import math
import sys


def law(args, z):
    """S, k, nu_t of the inflow law at height z above ground (Source/InflowProfile.H)."""
    z0, kap, c0 = args.z0, args.kappa, args.cmu0
    if args.law == "log_law":
        s = args.ustar / kap * math.log((z + z0) / z0)
        us = args.ustar
    else:
        zz = max(z, z0)   # below z0 the law's slope at z0 sets k (its own slope diverges)
        s = args.speed * (z / args.zref) ** args.alpha
        dsdz = args.alpha * args.speed * (zz / args.zref) ** args.alpha / zz
        us = kap * (z + z0) * dsdz
    return s, us * us / (c0 * c0), kap * us * (z + z0)


def interp(zs, vs, z):
    if z <= zs[0]:
        return vs[0]
    for a in range(len(zs) - 1):
        if zs[a] <= z <= zs[a + 1]:
            t = (z - zs[a]) / (zs[a + 1] - zs[a])
            return vs[a] + t * (vs[a + 1] - vs[a])
    return vs[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("plotfile")
    ap.add_argument("readerdir")
    ap.add_argument("--law", choices=["power_law", "log_law"], required=True)
    ap.add_argument("--speed", type=float, default=8.7502)
    ap.add_argument("--zref", type=float, default=10.0)
    ap.add_argument("--alpha", type=float, default=0.2113)
    ap.add_argument("--ustar", type=float, default=0.61726)
    ap.add_argument("--z0", type=float, default=0.03)
    ap.add_argument("--kappa", type=float, default=0.41)
    ap.add_argument("--cmu0", type=float, default=0.5562)
    ap.add_argument("--stations", type=float, nargs="+", default=None,
                    help="x positions (m); default: first cell, 0, last cell")
    ap.add_argument("--heights", type=float, nargs="+", default=[10.0, 20.0, 40.0, 80.0, 160.0, 300.0])
    ap.add_argument("--max-u-drift", type=float, default=None,
                    help="fail when |U/U_law - 1| exceeds this at any height and station")
    args = ap.parse_args()

    sys.path.insert(0, args.readerdir)
    import erf_plotfile as ep  # noqa: E402

    fields = ["x_velocity", "y_velocity", "KE", "Kmv", "z_phys", "ustar"]
    hdr, d = ep.read_fields(args.plotfile, fields)
    lo, hi = hdr["lo"], hdr["hi"]
    nx, ny, nz = [hi[i] - lo[i] + 1 for i in range(3)]
    x0, dx = hdr["prob_lo"][0], hdr["dx"][0]
    xc = [x0 + (i + 0.5) * dx for i in range(nx)]

    def column(i):
        mean = lambda f, k: sum(d[f][i][j][k] for j in range(ny)) / ny  # noqa: E731
        z = [mean("z_phys", k) for k in range(nz)]
        u = [math.hypot(mean("x_velocity", k), mean("y_velocity", k)) for k in range(nz)]
        return z, u, [mean("KE", k) for k in range(nz)], [mean("Kmv", k) for k in range(nz)], mean("ustar", 0)

    stations = args.stations or [xc[0], 0.0, xc[-1]]
    i_of = lambda x: min(range(nx), key=lambda i: abs(xc[i] - x))  # noqa: E731
    z1, u1, k1, n1, _ = column(0)

    worst = 0.0
    print("law %s, z0 %g m; columns averaged over %d cells in y" % (args.law, args.z0, ny))
    print("%9s %7s %8s %8s %8s %8s %8s %8s %8s" % ("x (m)", "z (m)", "U", "U/law-1", "U/in-1",
                                                   "k", "k/law-1", "k/in-1", "nut/law-1"))
    for x in stations:
        i = i_of(x)
        z, u, k, n, _ = column(i)
        for h in args.heights:
            if h < z[0] or h > z[-1]:
                continue
            s_l, k_l, n_l = law(args, h)
            U, K, N = interp(z, u, h), interp(z, k, h), interp(z, n, h)
            du = U / s_l - 1.0
            worst = max(worst, abs(du))
            print("%9.1f %7.1f %8.4f %+8.4f %+8.4f %8.4f %+8.4f %+8.4f %+8.4f" % (
                xc[i], h, U, du, U / interp(z1, u1, h) - 1.0,
                K, K / k_l - 1.0, K / interp(z1, k1, h) - 1.0, N / n_l - 1.0))
    _, us_l, _ = law(args, 0.0)
    us_l = math.sqrt(us_l) * args.cmu0
    us = [column(i)[4] for i in range(nx)]
    print("ground u*: law %.4f m/s; along x min %.4f max %.4f, first %.4f, x=0 %.4f, last %.4f" % (
        us_l, min(us), max(us), us[0], us[i_of(0.0)], us[-1]))
    print("largest |U/law - 1| at the listed heights and stations: %.4f" % worst)
    if args.max_u_drift is not None and worst > args.max_u_drift:
        print("hhabl_drift: FAIL (bound %g)" % args.max_u_drift)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
