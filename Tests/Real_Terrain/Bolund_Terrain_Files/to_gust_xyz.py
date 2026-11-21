#!/usr/bin/env python3
"""Terrain and roughness maps in other formats as GUST point files: "x y z" per line for
terrain.file (terrain.type = xyz) and "x y z0" per line for bc.surface_z0_file.

Input formats, chosen by the file's extension:
  .tif .tiff   single-band GeoTIFF (read by tools/tiff_to_xyz.py: uncompressed or LZW strips, no GDAL)
  .asc         ESRI ASCII grid (ncols, nrows, xllcorner or xllcenter, yllcorner or yllcenter, cellsize,
               NODATA_value; rows from north to south)
  .grd         Surfer ASCII grid (DSAA; nodes from xlo to xhi and ylo to yhi, rows from south to north)
  other        points, "x y value" per line ('#' comments, blanks or commas; other lines skipped)

Commands:
  to_gust_xyz.py terrain IN OUT [options]            heights (m)
  to_gust_xyz.py z0 IN OUT [options]                 IN holds roughness lengths (m)
  to_gust_xyz.py z0 IN OUT --classes 11:0.0003,71:0.015 [options]
                                                     IN holds land-cover codes; each code's z0 is given
  to_gust_xyz.py z0 IN OUT --water-level 0.75 [--z0-water 0.0003] [--z0-land 0.015] [--land-x-min 327]
                                                     IN holds heights: water at or below the level, land above
                                                     it and (with --land-x-min) for x above that value
  to_gust_xyz.py selftest                            the readers and the three z0 modes on small made-up grids

Options for terrain and z0: --stride N (every N-th row and column of a grid from the file's first row and
column, or every N-th point),
--shift DX DY DZ (subtracted from x, y and, for terrain, the height: a local frame),
--crop XLO XHI YLO YHI (kept points, in the output frame).
"""
import argparse
import math
import os
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ME = "to_gust_xyz.py"


def stop(msg):
    sys.exit(f"{ME}: {msg}")


def read_asc(path):
    """ESRI ASCII grid: (x, y, value) of the cell centres, no-data cells left out."""
    head, rows = {}, []
    with open(path) as f:
        for line in f:
            tok = line.split()
            if not tok:
                continue
            if tok[0][0].isalpha():
                head[tok[0].lower()] = float(tok[1])
            else:
                rows.extend(float(t) for t in tok)
    for key in ("ncols", "nrows", "cellsize"):
        if key not in head:
            stop(f"{path}: no {key} line (not an ESRI ASCII grid)")
    nx, ny, d = int(head["ncols"]), int(head["nrows"]), head["cellsize"]
    if len(rows) != nx * ny:
        stop(f"{path}: {len(rows)} values, expected ncols x nrows = {nx * ny}")
    x0 = head["xllcorner"] + 0.5 * d if "xllcorner" in head else head.get("xllcenter")
    y0 = head["yllcorner"] + 0.5 * d if "yllcorner" in head else head.get("yllcenter")
    if x0 is None or y0 is None:
        stop(f"{path}: no xllcorner / xllcenter or yllcorner / yllcenter line")
    a = np.array(rows).reshape(ny, nx)[::-1, :]     # the file's first row is the northern one
    if "nodata_value" in head:
        a = np.where(a == head["nodata_value"], np.nan, a)
    return x0 + d * np.arange(nx), y0 + d * np.arange(ny), a


def read_grd(path):
    """Surfer ASCII grid (DSAA): (x, y, value) of the nodes, blanked nodes left out."""
    with open(path) as f:
        tok = f.read().split()
    if not tok or tok[0] != "DSAA":
        stop(f"{path}: does not start with DSAA (only Surfer ASCII grids are read; convert a binary grid first)")
    nx, ny = int(tok[1]), int(tok[2])
    xlo, xhi, ylo, yhi = (float(t) for t in tok[3:7])
    vals = np.array([float(t) for t in tok[9:]])
    if vals.size != nx * ny:
        stop(f"{path}: {vals.size} values, expected nx x ny = {nx * ny}")
    a = vals.reshape(ny, nx)                        # the file's first row is the southern one
    a = np.where(a >= 1.70141e38, np.nan, a)
    return np.linspace(xlo, xhi, nx), np.linspace(ylo, yhi, ny), a


def read_tif(path):
    sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "tools"))
    import tiff_to_xyz
    a, x0, y0, sx, sy = tiff_to_xyz.read_tiff(path)
    ny, nx = a.shape
    a = np.where(a > -1e30, a, np.nan)[::-1, :]     # the file's first row is the northern one
    return x0 + sx * (np.arange(nx) + 0.5), y0 - sy * (ny - np.arange(ny) - 0.5), a


def read_points(path, stride):
    pts = []
    with open(path) as f:
        for line in f:
            tok = line.split("#")[0].replace(",", " ").split()
            if len(tok) != 3:
                continue
            try:
                pts.append([float(t) for t in tok])
            except ValueError:
                continue
    if not pts:
        stop(f"{path}: holds no 'x y value' line")
    p = np.array(pts)[::stride]
    return p[:, 0], p[:, 1], p[:, 2]


def read_any(path, stride):
    """Flat arrays x, y, value of the points of any supported file, no-data left out."""
    if not os.path.isfile(path):
        stop(f"{path}: no such file")
    ext = os.path.splitext(path)[1].lower()
    reader = {".asc": read_asc, ".grd": read_grd, ".tif": read_tif, ".tiff": read_tif}.get(ext)
    if reader is None:
        x, y, v = read_points(path, stride)
    else:
        xs, ys, a = reader(path)
        rows = np.arange(ys.size)
        if ext != ".grd":                           # the stride starts at the file's first row (the northern one)
            rows = rows[::-1]
        rows = np.sort(rows[::stride])
        X, Y = np.meshgrid(xs[::stride], ys[rows])
        a = a[rows, ::stride]
        x, y, v = X.ravel(), Y.ravel(), a.ravel()
    keep = np.isfinite(v)
    if not keep.any():
        stop(f"{path}: every value is no data")
    return x[keep], y[keep], v[keep]


def parse_classes(text):
    table = {}
    for item in text.split(","):
        if item.count(":") != 1:
            stop(f"--classes: '{item}' is not CODE:Z0")
        code, z0 = item.split(":")
        try:
            table[float(code)] = float(z0)
        except ValueError:
            stop(f"--classes: '{item}' is not CODE:Z0 with two numbers")
    return table


def roughness(x, v, a):
    """z0 of every point from the values of the input, by the mode the options choose."""
    modes = (a.classes is not None) + (a.water_level is not None)
    if modes > 1:
        stop("--classes and --water-level are two ways to give z0; use one")
    if a.water_level is None and (a.land_x_min is not None or a.z0_water is not None or a.z0_land is not None):
        stop("--z0-water, --z0-land and --land-x-min need --water-level")
    if a.classes is not None:
        table = parse_classes(a.classes)
        z0 = np.full(v.shape, np.nan)
        for code, value in table.items():
            z0[v == code] = value
        missing = sorted(set(v[np.isnan(z0)].tolist()))
        if missing:
            stop("--classes: no z0 for the code(s) " + ", ".join(f"{m:g}" for m in missing[:20]))
    elif a.water_level is not None:
        zw = 0.0003 if a.z0_water is None else a.z0_water
        zl = 0.015 if a.z0_land is None else a.z0_land
        water = np.round(v, 3) <= a.water_level      # heights as the terrain file holds them (mm)
        if a.land_x_min is not None:
            water &= ~(x > a.land_x_min)
        z0 = np.where(water, zw, zl)
    else:
        z0 = v
    if not (np.isfinite(z0).all() and (z0 > 0).all()):
        stop(f"z0: every roughness length must be finite and > 0 (minimum found {np.nanmin(z0):g})")
    return z0


def convert(a):
    if a.stride < 1:
        stop(f"--stride = {a.stride}: must be >= 1")
    x, y, v = read_any(a.inp, a.stride)
    x_in = x.copy()
    x, y = x - a.shift[0], y - a.shift[1]
    if a.command == "terrain":
        v = v - a.shift[2]
    else:
        if a.shift[2] != 0:
            stop("--shift: DZ must be 0 for z0 (it shifts heights)")
        v = roughness(x_in - a.shift[0], v, a)
    if a.crop is not None:
        xlo, xhi, ylo, yhi = a.crop
        if not (xlo < xhi and ylo < yhi):
            stop("--crop XLO XHI YLO YHI: needs XLO < XHI and YLO < YHI")
        keep = (x >= xlo) & (x <= xhi) & (y >= ylo) & (y <= yhi)
        x, y, v = x[keep], y[keep], v[keep]
        if x.size == 0:
            stop("--crop leaves no point")
    what = "height (m)" if a.command == "terrain" else "z0 (m)"
    with open(a.out, "w") as out:
        out.write(f"# x y {what} from {os.path.basename(a.inp)} by {ME} {a.command}; {x.size} points\n")
        fmt = "%.3f %.3f %.3f\n" if a.command == "terrain" else "%.3f %.3f %.6g\n"
        for row in zip(x, y, v):
            out.write(fmt % row)
    extra = ""
    if a.command == "z0":
        extra = f", geometric mean {math.exp(float(np.mean(np.log(v)))):.4g}"
    print(f"{a.out}: {x.size} points; x {x.min():.2f} .. {x.max():.2f}, y {y.min():.2f} .. {y.max():.2f}, "
          f"{what} {v.min():.4g} .. {v.max():.4g}{extra}")


def selftest():
    """A 3 x 2 grid h(x, y) = 10 x + y in each format gives the same points; the three z0 modes."""
    def args(command, inp, out, **kw):
        base = dict(command=command, inp=inp, out=out, stride=1, shift=[0.0, 0.0, 0.0], crop=None, classes=None,
                    water_level=None, z0_water=None, z0_land=None, land_x_min=None)
        base.update(kw)
        return argparse.Namespace(**base)

    def load(path):
        p = np.loadtxt(path)
        return p[np.lexsort((p[:, 0], p[:, 1]))]

    xs, ys = [100.0, 102.0, 104.0], [50.0, 52.0]
    want = np.array([[x, y, 10 * x + y] for y in ys for x in xs])
    with tempfile.TemporaryDirectory() as d:
        asc, grd, xyz = (os.path.join(d, n) for n in ("t.asc", "t.grd", "t.txt"))
        with open(asc, "w") as f:
            f.write("ncols 3\nnrows 2\nxllcorner 99\nyllcorner 49\ncellsize 2\nNODATA_value -9999\n")
            f.write("1052 1072 1092\n1050 1070 1090\n")         # north row first
        with open(grd, "w") as f:
            f.write("DSAA\n3 2\n100 104\n50 52\n1050 1092\n1050 1070 1090\n1052 1072 1092\n")   # south row first
        with open(xyz, "w") as f:
            f.write("# points\n" + "".join("%g, %g, %g\n" % tuple(r) for r in want))
        for src in (asc, grd, xyz):
            out = os.path.join(d, "h.xyz")
            convert(args("terrain", src, out))
            got = load(out)
            assert got.shape == want.shape and np.abs(got - want).max() < 1e-9, (src, got)
        out = os.path.join(d, "s.xyz")
        convert(args("terrain", asc, out, shift=[100.0, 50.0, 1050.0], crop=[0.0, 2.0, 0.0, 2.0]))
        got = load(out)
        assert np.abs(got - np.array([[0, 0, 0], [2, 0, 20], [0, 2, 2], [2, 2, 22]])).max() < 1e-9, got
        out = os.path.join(d, "z.xyz")
        convert(args("z0", asc, out, water_level=1070.0, land_x_min=103.0))
        got = load(out)       # water: h <= 1070 and x <= 103
        assert np.allclose(got[:, 2], [0.0003, 0.0003, 0.015, 0.0003, 0.015, 0.015]), got
        cls = os.path.join(d, "c.txt")
        with open(cls, "w") as f:
            f.write("0 0 11\n1 0 71\n")
        convert(args("z0", cls, out, classes="11:0.0003,71:0.015"))
        assert np.allclose(load(out)[:, 2], [0.0003, 0.015])
        for bad in (dict(classes="11:0.0003"), {}):               # a code without z0; raw values 11 and 71 pass,
            try:                                                  # so the second uses a zero below
                if not bad:
                    with open(cls, "w") as f:
                        f.write("0 0 0.03\n1 0 0\n")
                convert(args("z0", cls, out, **bad))
            except SystemExit as e:
                assert ME in str(e), e
            else:
                raise AssertionError("a roughness map with a missing or zero z0 was accepted")
    print(f"{ME} selftest: passed")


def main(argv):
    if argv[:1] == ["selftest"]:
        selftest()
        return
    ap = argparse.ArgumentParser(prog=ME, description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=("terrain", "z0"))
    ap.add_argument("inp", metavar="IN")
    ap.add_argument("out", metavar="OUT")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--shift", type=float, nargs=3, default=[0.0, 0.0, 0.0], metavar=("DX", "DY", "DZ"))
    ap.add_argument("--crop", type=float, nargs=4, metavar=("XLO", "XHI", "YLO", "YHI"))
    ap.add_argument("--classes", help="z0 only: CODE:Z0,CODE:Z0,... for an input of land-cover codes")
    ap.add_argument("--water-level", type=float, help="z0 only: the input holds heights; water at or below this")
    ap.add_argument("--z0-water", type=float, help="with --water-level (default 0.0003 m)")
    ap.add_argument("--z0-land", type=float, help="with --water-level (default 0.015 m)")
    ap.add_argument("--land-x-min", type=float, help="with --water-level: land for x above this (output frame)")
    a = ap.parse_args(argv)
    if a.command == "terrain" and (a.classes is not None or a.water_level is not None or a.z0_water is not None
                                   or a.z0_land is not None or a.land_x_min is not None):
        stop("--classes, --water-level, --z0-water, --z0-land and --land-x-min are options of z0")
    convert(a)


if __name__ == "__main__":
    main(sys.argv[1:])
