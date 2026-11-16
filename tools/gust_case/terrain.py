"""Terrain helpers: read elevation data, cut it to the domain, flatten the
border and write GUST's ``terrain.type = xyz`` point file.

GUST interpolates the point file to its ground columns by inverse-distance
weighting (``terrain.idw_neighbors``, ``terrain.idw_power``), so the points
need not form a grid: a raster in latitude and longitude is projected point by
point and written as it falls.

The border rule is amrTerrain's (backendInterface.py, the *Slope / *Flat
keys), written for all four sides at once. With d the distance of a point from
a lateral face of the domain, measured inwards,

    w_side = min(max((d - flat_width) / slope_width, 0), 1)
    w      = w_west w_east w_south w_north
    z      = w z_data + (1 - w) z_flat

so the ground is flat in a band ``flat_width`` wide (and outside the domain),
and the data's own from ``flat_width + slope_width`` inwards.
"""
import os

import numpy as np

from . import geo
from .errors import CaseError

NODATA_BELOW = -1.0e30
METERS_PER_DEGREE = 111194.9   # mean Earth radius 6371.0088 km times pi / 180


def geographic_message(file_key, path):
    """The refusal of a file read as longitude and latitude with no domain.center_lat_lon."""
    coords = file_key.rsplit(".", 1)[0] + ".coordinates"
    return ("%s = %s is read as longitude and latitude (%s = auto takes values within 180 and 90 as degrees): give "
            "domain.center_lat_lon, or %s = projected if the file is in metres" % (file_key, path, coords, coords))


class Raster:
    """A regular grid of values: x and y of the cell centres (increasing), values[j, i]."""

    def __init__(self, x, y, values, source):
        self.x = np.asarray(x, dtype=float)
        self.y = np.asarray(y, dtype=float)
        self.values = np.asarray(values, dtype=float)
        self.source = source

    def is_geographic(self):
        """True when the extent fits in longitude ±180, latitude ±90."""
        return bool(abs(self.x[0]) <= 180.0 and abs(self.x[-1]) <= 180.0 and abs(self.y[0]) <= 90.0 and abs(self.y[-1]) <= 90.0)

    def spacing(self):
        dx = (self.x[-1] - self.x[0]) / max(len(self.x) - 1, 1)
        dy = (self.y[-1] - self.y[0]) / max(len(self.y) - 1, 1)
        return dx, dy


def guess_format(path, key):
    """geotiff, esri_ascii or xyz from the file's extension."""
    ext = os.path.splitext(path)[1].lower()
    if ext in (".tif", ".tiff"):
        return "geotiff"
    if ext in (".asc", ".grd"):
        return "esri_ascii"
    if ext in (".xyz", ".txt", ".dat", ".csv"):
        return "xyz"
    raise CaseError("%s = %s: cannot tell the format from the extension '%s'; give the format key beside it"
                    % (key, path, ext))


def _oriented(x, y, values, source):
    """A Raster with x and y increasing (rasters are stored north row first)."""
    if len(x) > 1 and x[1] < x[0]:
        x, values = x[::-1], values[:, ::-1]
    if len(y) > 1 and y[1] < y[0]:
        y, values = y[::-1], values[::-1, :]
    return Raster(x, y, values, source)


def read_esri_ascii(path, key):
    """An ESRI ASCII grid (ncols, nrows, xllcorner or xllcenter, cellsize, NODATA_value)."""
    header = {}
    try:
        with open(path, "r", encoding="utf-8") as handle:
            rows_start = 0
            for line in handle:
                parts = line.split()
                if len(parts) == 2 and parts[0][0].isalpha():
                    header[parts[0].lower()] = float(parts[1])
                    rows_start += 1
                else:
                    break
        values = np.loadtxt(path, skiprows=rows_start, dtype=float, ndmin=2)
    except (OSError, ValueError) as err:
        raise CaseError("%s = %s: cannot read the ESRI ASCII grid (%s)" % (key, path, err)) from None
    try:
        nx, ny, cell = int(header["ncols"]), int(header["nrows"]), header["cellsize"]
        if "xllcenter" in header:
            x0, y0 = header["xllcenter"], header["yllcenter"]
        else:
            x0, y0 = header["xllcorner"] + 0.5 * cell, header["yllcorner"] + 0.5 * cell
    except KeyError as err:
        raise CaseError("%s = %s: the ESRI ASCII header has no %s" % (key, path, err)) from None
    if values.shape != (ny, nx):
        raise CaseError("%s = %s: %d x %d values, but the header says %d rows of %d"
                        % (key, path, values.shape[0], values.shape[1], ny, nx))
    if "nodata_value" in header:
        values = np.where(values == header["nodata_value"], np.nan, values)
    x = x0 + cell * np.arange(nx)
    y = y0 + cell * np.arange(ny)
    return Raster(x, y, values[::-1, :], path)   # the file's first row is the northern one


def read_geotiff(path, key):
    """A single-band GeoTIFF: with rasterio when it is installed, else the reader
    of tools/tiff_to_xyz.py (little-endian integer or float strips, none or LZW; GDAL_NODATA)."""
    try:
        import rasterio   # optional
    except ImportError:
        rasterio = None
    if rasterio is not None:
        try:
            with rasterio.open(path) as data:
                values = data.read(1).astype(float)
                t = data.transform
                nodata = data.nodata
                if t.b != 0.0 or t.d != 0.0:
                    raise CaseError("%s = %s: the raster is rotated; reproject it to a north-up grid" % (key, path))
                x = t.c + t.a * (np.arange(data.width) + 0.5)
                y = t.f + t.e * (np.arange(data.height) + 0.5)
        except CaseError:
            raise
        except Exception as err:   # rasterio raises its own error types
            raise CaseError("%s = %s: cannot read the GeoTIFF (%s)" % (key, path, err)) from None
        if nodata is not None:
            values = np.where(values == nodata, np.nan, values)
    else:
        import sys
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        try:
            from tiff_to_xyz import read_tiff
            values, x0, y0, sx, sy = read_tiff(path)
        except Exception as err:   # the small reader asserts on what it does not cover
            raise CaseError("%s = %s: the built-in reader takes little-endian strips of 8 to 32 bit integers or "
                            "32/64 bit floats, uncompressed or LZW without a predictor (%s: %s); install rasterio "
                            "for other GeoTIFFs, or convert the file"
                            % (key, path, type(err).__name__, err)) from None
        finally:
            sys.path.pop(0)
        values = np.array(values, dtype=float)
        ny, nx = values.shape
        x = x0 + sx * (np.arange(nx) + 0.5)
        y = y0 - sy * (np.arange(ny) + 0.5)
    values = np.where(values <= NODATA_BELOW, np.nan, values)
    return _oriented(x, y, values, path)


def read_points(path, key, columns=3):
    """An x y value point file (# starts a comment; commas allowed)."""
    try:
        with open(path, "r", encoding="utf-8") as handle:
            text = handle.read().replace(",", " ")
        rows = [line.split("#", 1)[0].split() for line in text.splitlines()]
        data = np.array([[float(v) for v in row[:columns]] for row in rows if row], dtype=float)
    except (OSError, ValueError) as err:
        raise CaseError("%s = %s: cannot read the point file (%s)" % (key, path, err)) from None
    if data.ndim != 2 or data.shape[1] != columns or data.shape[0] == 0:
        raise CaseError("%s = %s: expected lines of %d numbers" % (key, path, columns))
    return data[:, 0], data[:, 1], data[:, 2]


def read_raster(path, fmt, key):
    if not os.path.isfile(path):
        raise CaseError("%s = %s: no such file" % (key, path))
    if fmt == "geotiff":
        return read_geotiff(path, key)
    if fmt == "esri_ascii":
        return read_esri_ascii(path, key)
    raise AssertionError(fmt)


class Frame:
    """The case frame: where the domain sits in the data's projected frame.

    ``center`` is the domain centre in the data's frame (UTM when the centre
    was given as latitude and longitude); ``x_lo`` .. ``y_hi`` are the domain's
    faces in the case frame, whose origin is ``origin`` in the data's frame.
    """

    def __init__(self, center, west, east, south, north, origin_at="center", lat_lon=None, zone=None, south_hemisphere=None):
        self.center = (float(center[0]), float(center[1]))
        self.lat_lon = lat_lon
        self.zone = zone
        self.south_hemisphere = south_hemisphere
        if origin_at == "center":
            self.origin = self.center
            self.x_lo, self.x_hi, self.y_lo, self.y_hi = -west, east, -south, north
        else:
            self.origin = (self.center[0] - west, self.center[1] - south)
            self.x_lo, self.x_hi, self.y_lo, self.y_hi = 0.0, west + east, 0.0, south + north

    def to_case(self, x, y):
        return np.asarray(x, dtype=float) - self.origin[0], np.asarray(y, dtype=float) - self.origin[1]

    def from_lat_lon(self, lat, lon):
        """Case-frame x, y of latitude and longitude, in the centre's UTM zone."""
        east, north, _, _ = geo.latlon_to_utm(lat, lon, zone=self.zone, south=self.south_hemisphere)
        return self.to_case(east, north)


def raster_points(raster, frame, geographic, margin, key, stride_spacing=None):
    """The raster's points inside the domain plus ``margin``, in the case frame.

    ``stride_spacing``: thin the raster to about this spacing (m), never finer
    than the data. Returns x, y, values (1-D) and the spacing kept (m).
    """
    x_lo, x_hi = frame.x_lo - margin, frame.x_hi + margin
    y_lo, y_hi = frame.y_lo - margin, frame.y_hi + margin
    dx, dy = raster.spacing()
    if geographic:
        if frame.lat_lon is None:
            raise CaseError(geographic_message(key, raster.source))
        lat0, lon0 = frame.lat_lon
        # the farthest face from the centre, with 10 % for the projection's convergence and scale
        cx0, cy0 = frame.center[0] - frame.origin[0], frame.center[1] - frame.origin[1]
        reach = 1.1 * max(cx0 - x_lo, x_hi - cx0, cy0 - y_lo, y_hi - cy0) * np.sqrt(2.0)
        dlat = reach / METERS_PER_DEGREE
        dlon = reach / (METERS_PER_DEGREE * max(np.cos(np.radians(lat0)), 0.01))
        ii = np.nonzero((raster.x >= lon0 - dlon) & (raster.x <= lon0 + dlon))[0]
        jj = np.nonzero((raster.y >= lat0 - dlat) & (raster.y <= lat0 + dlat))[0]
        native = min(dx * METERS_PER_DEGREE * np.cos(np.radians(lat0)), dy * METERS_PER_DEGREE)
    else:
        ii = np.nonzero((raster.x >= frame.origin[0] + x_lo) & (raster.x <= frame.origin[0] + x_hi))[0]
        jj = np.nonzero((raster.y >= frame.origin[1] + y_lo) & (raster.y <= frame.origin[1] + y_hi))[0]
        native = min(dx, dy)
    if ii.size == 0 or jj.size == 0:
        raise CaseError("%s = %s: no data inside the domain; the file spans x %.6g .. %.6g, y %.6g .. %.6g "
                        "(check domain.center or domain.center_lat_lon)"
                        % (key, raster.source, raster.x[0], raster.x[-1], raster.y[0], raster.y[-1]))
    stride = 1
    if stride_spacing is not None and native > 0.0:
        stride = max(1, int(np.floor(stride_spacing / native + 1.0e-9)))
    ii, jj = ii[::stride], jj[::stride]
    gx, gy = np.meshgrid(raster.x[ii], raster.y[jj])
    values = raster.values[np.ix_(jj, ii)]
    if geographic:
        cx, cy = frame.from_lat_lon(gy, gx)
    else:
        cx, cy = frame.to_case(gx, gy)
    keep = (cx >= x_lo) & (cx <= x_hi) & (cy >= y_lo) & (cy <= y_hi)
    return cx[keep], cy[keep], values[keep], native * stride


def check_cover(x, y, frame, spacing, key, source):
    """Refuse data that does not reach all four faces of the domain."""
    slack = 1.5 * spacing
    gaps = []
    if x.size == 0:
        raise CaseError("%s = %s: no data inside the domain (check domain.center or domain.center_lat_lon)" % (key, source))
    if x.min() > frame.x_lo + slack:
        gaps.append("west (data from x = %.6g, face at %.6g)" % (x.min(), frame.x_lo))
    if x.max() < frame.x_hi - slack:
        gaps.append("east (data to x = %.6g, face at %.6g)" % (x.max(), frame.x_hi))
    if y.min() > frame.y_lo + slack:
        gaps.append("south (data from y = %.6g, face at %.6g)" % (y.min(), frame.y_lo))
    if y.max() < frame.y_hi - slack:
        gaps.append("north (data to y = %.6g, face at %.6g)" % (y.max(), frame.y_hi))
    if gaps:
        raise CaseError("%s = %s does not cover the domain on the %s; make the domain smaller or give a larger file"
                        % (key, source, "; ".join(gaps)))


def border_weight(x, y, frame, flat_width, slope_width):
    """The weight of the data (1: untouched, 0: flat) at points of the case frame."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    w = np.ones_like(x)
    for d in (x - frame.x_lo, frame.x_hi - x, y - frame.y_lo, frame.y_hi - y):
        if slope_width > 0.0:
            w = w * np.clip((d - flat_width) / slope_width, 0.0, 1.0)
        else:
            w = w * (d >= flat_width)
    return w


class Terrain:
    """The prepared ground: points in the case frame, and what was done to them."""

    def __init__(self, x, y, z, spacing, zero_shift, flat_level, notes):
        self.x, self.y, self.z = x, y, z
        self.spacing = spacing
        self.zero_shift = zero_shift      # subtracted from the data's elevations
        self.flat_level = flat_level      # elevation of the flat border in the case frame, or None
        self.notes = notes
        self.lowest = float(np.min(z))
        self.highest = float(np.max(z))

    @property
    def relief(self):
        return self.highest - self.lowest

    def write(self, path):
        with open(path, "w", encoding="utf-8") as out:
            out.write("# x y z (m), case frame; %d points, spacing about %.6g m; z = elevation - %.6g m\n"
                      % (self.x.size, self.spacing, self.zero_shift))
            np.savetxt(out, np.column_stack([self.x, self.y, self.z]), fmt="%.3f")


def prepare(path, fmt, coordinates, frame, margin, point_spacing, zero_level, flat_width, slope_width, border_level,
            key="terrain.file", border_frame=None):
    """Read ``path`` and return the Terrain of the domain.

    zero_level: 'lowest', 'none' or a number. border_level: 'lowest', 'mean' or
    a number in the file's datum; used when flat_width or slope_width is > 0.
    border_frame: the frame whose faces the border is measured from (a nested
    patch takes its parent's, with the parent's levels as numbers).
    """
    notes = []
    if fmt == "auto":
        fmt = guess_format(path, key)
    if fmt == "xyz":
        if not os.path.isfile(path):
            raise CaseError("%s = %s: no such file" % (key, path))
        px, py, pz = read_points(path, key)
        geographic = coordinates == "geographic" or (
            coordinates == "auto" and bool(np.all(np.abs(px) <= 180.0) and np.all(np.abs(py) <= 90.0)))
        if geographic:
            if frame.lat_lon is None:
                raise CaseError(geographic_message(key, path))
            cx, cy = frame.from_lat_lon(py, px)
        else:
            cx, cy = frame.to_case(px, py)
        keep = ((cx >= frame.x_lo - margin) & (cx <= frame.x_hi + margin) &
                (cy >= frame.y_lo - margin) & (cy <= frame.y_hi + margin))
        cx, cy, cz = cx[keep], cy[keep], pz[keep]
        area = (frame.x_hi - frame.x_lo + 2 * margin) * (frame.y_hi - frame.y_lo + 2 * margin)
        spacing = float(np.sqrt(area / max(cx.size, 1)))
        if point_spacing is not None:
            notes.append("terrain.point_spacing does not thin a point file; all %d points are kept" % cx.size)
    else:
        raster = read_raster(path, fmt, key)
        geographic = coordinates == "geographic" or (coordinates == "auto" and raster.is_geographic())
        cx, cy, cz, spacing = raster_points(raster, frame, geographic, margin, key, stride_spacing=point_spacing)
    check_cover(cx, cy, frame, spacing, key, path)
    if not np.all(np.isfinite(cz)):
        raise CaseError("%s = %s: %d of the %d points in the domain have no data; fill the gaps or move the domain"
                        % (key, path, int(np.sum(~np.isfinite(cz))), cz.size))
    notes.append("%s: %s coordinates, %d points kept at about %.4g m" % (os.path.basename(path),
                                                                          "geographic" if geographic else "projected",
                                                                          cx.size, spacing))
    flat = None
    if flat_width > 0.0 or slope_width > 0.0:
        if border_level == "lowest":
            flat = float(np.min(cz))
        elif border_level == "mean":
            flat = float(np.mean(cz))
        else:
            flat = float(border_level)
        w = border_weight(cx, cy, border_frame or frame, flat_width, slope_width)
        cz = w * cz + (1.0 - w) * flat
        notes.append("border: flat at %.6g m over %.6g m, then %.6g m of slope to the data" % (flat, flat_width, slope_width))
    if zero_level == "lowest":
        shift = float(np.min(cz))
    elif zero_level == "none":
        shift = 0.0
    else:
        shift = float(zero_level)
    cz = cz - shift
    return Terrain(cx, cy, cz, spacing, shift, None if flat is None else flat - shift, notes)
