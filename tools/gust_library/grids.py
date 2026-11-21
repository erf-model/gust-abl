"""Library fields as WindNinja-format files: ESRI ASCII grids with a .prj, and CSV at points.

    geo = georeference(lib)                       where the case frame sits on the map (from case_resolved.yaml)
    write_grid(path, x, y, field, geo)            one ESRI ASCII grid (and path.prj when geo has one)
    read_grid(path)                               (header dict, field (nx, ny) with NaN for no data, x, y)
    xs, ys, note = grid_of(x, y)                  the square cells of the grids for the columns x, y
    resample(x, y, field, xs, ys)                 a field on those cells
    export_grids(lib, out_dir, prefix, fields, wind, extra, speed_names, units, prj)
                                                  <prefix>_vel.asc, <prefix>_ang.asc and one grid per extra
    export_points(lib, path, points, values, height)                    a CSV at points
    read_points(path)                                                   [(name, x, y)] from 'x y [name]' lines

The grids follow WindNinja's output files (name_vel.asc and name_ang.asc): the speed, and the direction the wind
comes from in degrees clockwise from north. Directions are relative to the grid's north (the y axis of the
case frame, the projection's northing); the meridian convergence of the projection is not removed. The cells
are the library's columns: the corner is half a cell beyond the first column centre and the first row written is
the north row. An ESRI ASCII grid has one cell size: when the columns are not square (the case builder rounds
each cell count on its own), the fields are resampled onto square cells of the smaller side, from the domain's
south-west corner, as many as fit in each direction; each value is bilinear between the column centres, and the
nearest column's value within half a column of the domain's faces. The direction is then taken from the
resampled wind components, not interpolated as an angle. The CSV at points is not resampled.

The case frame sits on the map at the origin of case_resolved.yaml's domain section (the domain centre or its
south-west corner, in the projected frame of domain.center, or in UTM when the centre was given as
domain.center_lat_lon). A .prj (ESRI WKT of WGS 84 / UTM) is written only for a centre given in latitude and
longitude; for a centre in projected metres, give the projection's own .prj file (prj= of export_grids, --prj
of the command line), else none is written.
"""
import math
import os

import numpy as np
import yaml

from .index import LibraryError

NODATA = -9999.0
SQUARE_TOL = 1.0e-6         # relative difference of the cell's two sides
UNIFORM_TOL = 1.0e-6        # relative spread of the column spacing
SPEED_UNITS = {"mps": 1.0, "mph": 3600.0 / 1609.344, "kph": 3.6}
CALM_FRACTION = 1.0e-12     # no direction where the wind is at most this fraction of the speed field


class Georeference:
    """origin: the case frame's (0, 0) in the projected frame; zone, south: the UTM zone, or None; source: text."""

    def __init__(self, origin=(0.0, 0.0), zone=None, south=False, source="the case frame (no georeference)"):
        self.origin = (float(origin[0]), float(origin[1]))
        self.zone = zone
        self.south = bool(south)
        self.source = source

    def prj(self):
        """ESRI WKT of WGS 84 / UTM zone, or None."""
        if self.zone is None:
            return None
        return ('PROJCS["WGS_1984_UTM_Zone_%d%s",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",'
                'SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],'
                'UNIT["Degree",0.0174532925199433]],PROJECTION["Transverse_Mercator"],'
                'PARAMETER["False_Easting",500000.0],PARAMETER["False_Northing",%.1f],'
                'PARAMETER["Central_Meridian",%.1f],PARAMETER["Scale_Factor",0.9996],'
                'PARAMETER["Latitude_Of_Origin",0.0],UNIT["Meter",1.0]]'
                % (self.zone, "S" if self.south else "N", 1.0e7 if self.south else 0.0, 6.0 * self.zone - 183.0))


def georeference(lib):
    """The Georeference of a library from its case folder's case_resolved.yaml (the folder above wind_library/);
    the plain case frame when the file is not there."""
    from gust_case.plan import _resolve_domain          # gust_case sits beside gust_library in tools/
    path = os.path.join(os.path.dirname(os.path.abspath(lib.folder)), "case_resolved.yaml")
    if not os.path.isfile(path):
        return Georeference(source="the case frame (%s not found)" % path)
    try:
        with open(path, encoding="utf-8") as fh:
            case = yaml.safe_load(fh)
        problems = []
        frame = _resolve_domain(case, set(), problems)
    except (OSError, yaml.YAMLError, KeyError, TypeError, ValueError) as err:
        raise LibraryError("%s: cannot read the domain section (%s: %s)" % (path, type(err).__name__, err))
    if problems or frame is None:
        raise LibraryError("%s: domain section refused: %s" % (path, "; ".join(problems) or "no extent"))
    if frame.zone is not None:
        src = "UTM zone %d%s from domain.center_lat_lon" % (frame.zone, "S" if frame.south_hemisphere else "N")
    elif (case.get("domain") or {}).get("center") is not None:
        src = "the projected frame of domain.center (projection not known: give --prj)"
    else:
        src = "the case frame (%s gives no domain centre)" % path
    return Georeference(frame.origin, frame.zone, bool(frame.south_hemisphere), src)


def cell_size(x, y):
    """The cell size of the columns x (nx), y (ny); refuses uneven or non-square cells."""
    if len(x) < 2 or len(y) < 2:
        raise LibraryError("grid: %d x %d columns; an ASCII grid needs 2 or more each way" % (len(x), len(y)))
    dx, dy = np.diff(x), np.diff(y)
    for name, d in (("x", dx), ("y", dy)):
        if not (d.min() > 0.0 and (d.max() - d.min()) <= UNIFORM_TOL * d.mean()):
            raise LibraryError("grid: the column spacing in %s is not uniform (%g to %g m)" % (name, d.min(), d.max()))
    if abs(dx.mean() - dy.mean()) > SQUARE_TOL * dx.mean():
        raise LibraryError("grid: cells of %g m by %g m; an ESRI ASCII grid needs square cells" % (dx.mean(), dy.mean()))
    return float(dx.mean())


def write_grid(path, x, y, field, geo, prj_text=None):
    """field (nx, ny) at the cell centres x, y (case frame; square cells) as an ESRI ASCII grid; NaN as
    NODATA_value. Writes path with .prj for .asc beside it, holding prj_text when given, else geo's WKT when it
    has one."""
    field = np.asarray(field, dtype=float)
    if field.shape != (len(x), len(y)):
        raise LibraryError("grid: a field of shape %s on %d x %d columns" % (field.shape, len(x), len(y)))
    h = cell_size(x, y)
    xll = geo.origin[0] + float(x[0]) - 0.5 * h
    yll = geo.origin[1] + float(y[0]) - 0.5 * h
    with open(path, "w", encoding="ascii") as fh:
        fh.write("ncols %d\nnrows %d\nxllcorner %.6f\nyllcorner %.6f\ncellsize %.9g\nNODATA_value %g\n"
                 % (len(x), len(y), xll, yll, h, NODATA))
        rows = np.where(np.isfinite(field), field, NODATA).T[::-1]          # north row first, west to east
        for row in rows:
            fh.write(" ".join("%.9g" % a for a in row))
            fh.write("\n")
    text = prj_text if prj_text is not None else geo.prj()
    if not text:
        return [path]
    prj_path = os.path.splitext(path)[0] + ".prj"
    with open(prj_path, "w", encoding="utf-8") as fh:
        fh.write(text if text.endswith("\n") else text + "\n")
    return [path, prj_path]


def read_grid(path):
    """(header, field, x, y): the header keys lower case, field (ncols, nrows) west-east by south-north with NaN
    for no data, x and y the cell centres in the grid's own coordinates."""
    with open(path, encoding="ascii") as fh:
        header = {}
        for _ in range(6):
            key, value = fh.readline().split()
            header[key.lower()] = float(value)
        data = np.loadtxt(fh, ndmin=2)
    nx, ny = int(header["ncols"]), int(header["nrows"])
    if data.shape != (ny, nx):
        raise LibraryError("%s: %d x %d values, the header says %d rows of %d" % (path, *data.shape, ny, nx))
    field = data[::-1].T.copy()
    field[field == header["nodata_value"]] = np.nan
    h = header["cellsize"]
    x = header["xllcorner"] + (np.arange(nx) + 0.5) * h
    y = header["yllcorner"] + (np.arange(ny) + 0.5) * h
    return header, field, x, y


def direction_from(u, v, calm=0.0):
    """Degrees from north the wind (u east, v north) comes from; NaN where it is calm (|(u, v)| <= calm)."""
    frm = np.mod(270.0 - np.degrees(np.arctan2(v, u)), 360.0)
    return np.where(np.hypot(u, v) > calm, frm, np.nan)


def grid_of(x, y):
    """(xs, ys, note): the centres of the square cells of the grids for columns centred at x, y; note is empty
    when they are the columns, else says how they were resampled. Refuses uneven columns."""
    if len(x) < 2 or len(y) < 2:
        raise LibraryError("grid: %d x %d columns; an ASCII grid needs 2 or more each way" % (len(x), len(y)))
    d = {}
    for name, c in (("x", x), ("y", y)):
        diff = np.diff(c)
        if not (diff.min() > 0.0 and (diff.max() - diff.min()) <= UNIFORM_TOL * diff.mean()):
            raise LibraryError("grid: the column spacing in %s is not uniform (%g to %g m)"
                               % (name, diff.min(), diff.max()))
        d[name] = float(diff.mean())
    if abs(d["x"] - d["y"]) <= SQUARE_TOL * d["x"]:
        return np.asarray(x, dtype=float), np.asarray(y, dtype=float), ""
    h = min(d["x"], d["y"])
    out = []
    for c, dc in ((x, d["x"]), (y, d["y"])):
        n = int(math.floor(len(c) * dc / h * (1.0 + 1.0e-9)))
        out.append(float(c[0]) - 0.5 * dc + (np.arange(n) + 0.5) * h)
    note = ("columns of %g m by %g m resampled onto %d x %d square cells of %g m (bilinear)"
            % (d["x"], d["y"], len(out[0]), len(out[1]), h))
    return out[0], out[1], note


def resample(x, y, field, xs, ys):
    """field (nx, ny) at the column centres x, y on the cells xs, ys: bilinear (one direction then the other),
    the edge value beyond the outermost centres."""
    field = np.asarray(field, dtype=float)
    if len(xs) == len(x) and len(ys) == len(y) and np.array_equal(xs, x) and np.array_equal(ys, y):
        return field
    along_x = np.array([np.interp(xs, x, field[:, j]) for j in range(len(y))]).T        # (len(xs), ny)
    return np.array([np.interp(ys, y, along_x[i]) for i in range(len(xs))])               # (len(xs), len(ys))


def export_grids(lib, out_dir, prefix, fields, wind, extra=(), speed_names=(), units="mps", prj=None):
    """WindNinja-format grids in out_dir (made when missing): <prefix>_vel.asc (the speed field wind[0], in
    units) and <prefix>_ang.asc (the direction the wind of the components wind[1], wind[2] comes from), and
    <prefix>_<name>.asc for each name of extra (in units when the name is in speed_names). prj: a projection file
    copied beside each grid, else the georeference's own when it has one. Returns (files, georeference, note)."""
    if units not in SPEED_UNITS:
        raise LibraryError("--speed-units: %r; one of %s" % (units, ", ".join(SPEED_UNITS)))
    prj_text = None
    if prj:
        try:
            with open(prj, encoding="utf-8") as fh:                  # read once: it may be one of the outputs
                prj_text = fh.read()
        except OSError as err:
            raise LibraryError("--prj: cannot read %s (%s)" % (prj, err))
    geo = georeference(lib)
    xs, ys, note = grid_of(lib.x, lib.y)
    if os.path.exists(out_dir) and not os.path.isdir(out_dir):
        raise LibraryError("--grids: %s is a file, not a folder" % out_dir)
    os.makedirs(out_dir, exist_ok=True)

    def on_cells(name):
        return resample(lib.x, lib.y, fields[name], xs, ys)
    speed_name, u_name, v_name = wind
    speed = on_cells(speed_name)
    factor = SPEED_UNITS[units]
    files = write_grid(os.path.join(out_dir, prefix + "_vel.asc"), xs, ys, speed * factor, geo, prj_text)
    direction = direction_from(on_cells(u_name), on_cells(v_name), CALM_FRACTION * speed)
    files += write_grid(os.path.join(out_dir, prefix + "_ang.asc"), xs, ys, direction, geo, prj_text)
    for name in extra:
        value = on_cells(name) * (factor if name in speed_names else 1.0)
        files += write_grid(os.path.join(out_dir, "%s_%s.asc" % (prefix, name)), xs, ys, value, geo, prj_text)
    return files, geo, note


def export_points(lib, path, points, values, height):
    """A CSV of values (dict of arrays, one value per point) at points [(name, x, y)]: name, x, y (case frame),
    easting, northing (the case frame's origin added), height, then one column per value."""
    geo = georeference(lib)
    names = list(values)
    columns = ["name", "x", "y", "easting", "northing", "height"] + names
    rows = [[p[0], p[1], p[2], p[1] + geo.origin[0], p[2] + geo.origin[1], height]
            + [float(values[n][i]) for n in names] for i, p in enumerate(points)]
    return write_csv(path, columns, rows)


def write_csv(path, columns, rows):
    """A CSV file: a header line of the column names, then one line per row (numbers with 12 significant digits,
    so a UTM northing keeps 0.1 mm; NaN empty)."""
    def text(a):
        if isinstance(a, str):
            return a
        return "" if not math.isfinite(a) else "%.12g" % a
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(",".join(columns) + "\n")
        for row in rows:
            fh.write(",".join(text(a) for a in row) + "\n")
    return [path]


def read_points(path):
    """[(name, x, y)] from a text file of 'x y [name]' per line (spaces or commas; # starts a comment); a point
    with no name is named by its line number."""
    points = []
    try:
        with open(path, encoding="utf-8") as fh:
            lines = fh.readlines()
    except OSError as err:
        raise LibraryError("--points: cannot read %s (%s)" % (path, err))
    for n, line in enumerate(lines, 1):
        words = line.split("#", 1)[0].replace(",", " ").split()
        if not words:
            continue
        try:
            x, y = float(words[0]), float(words[1])
        except (ValueError, IndexError):
            raise LibraryError("%s:%d: %r; a point is 'x y [name]' (m, case frame)" % (path, n, line.rstrip("\n")))
        if not (math.isfinite(x) and math.isfinite(y)):
            raise LibraryError("%s:%d: (%r, %r) is not a finite point" % (path, n, x, y))
        points.append((" ".join(words[2:]) or "point%d" % n, x, y))
    if not points:
        raise LibraryError("--points: %s has no points" % path)
    return points
