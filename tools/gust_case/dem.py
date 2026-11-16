"""Download terrain for a case: SRTM-format tiles from the AWS Terrain Tiles.

The tiles are the "skadi" set of the Terrain Tiles on the AWS Open Data
Registry (https://registry.opendata.aws/terrain-tiles/, Mapzen / Tilezen
joerd): one gzip-compressed SRTM ``.hgt`` file per 1 x 1 degree, no key or
account needed. A tile is a square of big-endian 16-bit integers (3601 x 3601
at 1 arc-second, 1201 x 1201 at 3), the first row the northern edge, the
first column the western edge; neighbouring tiles share their edge row and
column; -32768 marks a void. Heights are metres above the EGM96 geoid.

The tiles merge several sources (SRTM, GMTED2010, 3DEP, ETOPO1 and national
models; tilezen/joerd docs/attribution.md lists them), so over the sea they
hold ETOPO1's sea floor, below 0 m. ``sea_floor="zero"`` sets those points to
0 m (the sea surface), the GUI's choice; ``"keep"`` leaves them, which keeps
land below sea level (Death Valley, the Dead Sea).

The download writes an ESRI ASCII grid in longitude and latitude, which the
case builder reads without rasterio (``terrain.format = esri_ascii``,
``terrain.coordinates`` auto finds degrees), and a ``.sources.txt`` file with
the attribution the data need. Only the standard library and numpy are used.

The approach follows FastWindTerrain's ``cases/casegen.py`` (the same tiles,
a cache of raw ``.hgt`` files and no partial file left behind) and the tile
layout of massconsistent_amr's ``tools/terrain_processing/terrain_reader_srtm.py``.
"""
import gzip
import http.client
import math
import os
import urllib.error
import urllib.request
import zlib

import numpy as np

from .errors import CaseError

TILE_URL = "https://s3.amazonaws.com/elevation-tiles-prod/skadi/{lat_dir}/{name}.hgt.gz"
VOID = -32768
TILE_SIZES = (3601, 1201)          # 1 and 3 arc-seconds
MAX_TILES = 16                     # a download larger than 4 x 4 degrees is refused (about 400 MB of tiles)
SEA_FLOOR = ("zero", "keep")

ATTRIBUTION = """\
Terrain data: the Terrain Tiles on the AWS Open Data Registry (https://registry.opendata.aws/terrain-tiles/),
skadi (SRTM .hgt) format, made by Mapzen / Tilezen joerd. Accessed {date}.
The tiles merge the sources below; credit the ones that cover your area
(https://github.com/tilezen/joerd/blob/master/docs/attribution.md):
  - United States 3DEP (formerly NED) and global GMTED2010 and SRTM terrain data courtesy of the U.S. Geological Survey.
  - Global ETOPO1 terrain data U.S. National Oceanic and Atmospheric Administration.
  - ArcticDEM terrain data DEM(s) were created from DigitalGlobe, Inc., imagery and funded under National Science
    Foundation awards 1043681, 1559691, and 1542736.
  - Australia terrain data (c) Commonwealth of Australia (Geoscience Australia) 2017.
  - Austria terrain data (c) offene Daten Oesterreichs - Digitales Gelaendemodell (DGM) Oesterreich.
  - Canada terrain data contains information licensed under the Open Government Licence - Canada.
  - Europe terrain data produced using Copernicus data and information funded by the European Union - EU-DEM layers.
  - Mexico terrain data source: INEGI, Continental relief, 2016.
  - New Zealand terrain data Copyright 2011 Crown copyright (c) Land Information New Zealand and the New Zealand
    Government (All rights reserved).
  - Norway terrain data (c) Kartverket.
  - United Kingdom terrain data (c) Environment Agency copyright and/or database right 2015. All rights reserved.
"""


def tile_name(lat, lon):
    """The tile whose south-west corner is at integer degrees (lat, lon): N57W008 for (57, -8)."""
    return "%s%02d%s%03d" % ("N" if lat >= 0 else "S", abs(lat), "E" if lon >= 0 else "W", abs(lon))


def tile_url(lat, lon):
    name = tile_name(lat, lon)
    return TILE_URL.format(lat_dir=name[:3], name=name)


def _tile_range(lo, hi, top):
    """The whole degrees of the tiles that hold a point of [lo, hi]: a tile whose only point in the box would be
    its shared edge (within one arc-second, the finest tile spacing) is left out, its neighbour has that edge."""
    eps = 1.0 / 3600.0 - 1e-9
    first = int(math.floor(lo + eps))
    last = min(int(math.ceil(hi - eps)), top)
    if last <= first:                             # the box lies within one arc-second of a whole degree
        first = min(int(math.floor(lo + 0.5 * (hi - lo))), top - 1)
        last = first + 1
    return range(first, last)


def tiles_for_box(south, north, west, east):
    """(lat, lon) of the south-west corners of the tiles that cover the box, south to north, west to east."""
    if not (-90.0 <= south < north <= 90.0):
        raise CaseError("terrain download: the box's latitudes %g .. %g must increase within -90 .. 90" % (south, north))
    if not (-180.0 <= west < east <= 180.0):
        raise CaseError("terrain download: the box's longitudes %g .. %g must increase within -180 .. 180 (a box "
                        "across the 180th meridian is not supported)" % (west, east))
    lats = _tile_range(south, north, 90)
    lons = _tile_range(west, east, 180)
    tiles = [(a, o) for a in lats for o in lons]
    if len(tiles) > MAX_TILES:
        raise CaseError("terrain download: the box needs %d tiles of 1 x 1 degree; at most %d are fetched at once"
                        % (len(tiles), MAX_TILES))
    return tiles


def box_for_frame(frame, margin):
    """(south, north, west, east) in degrees around the domain of ``frame`` plus ``margin``.

    The same box the builder reads from a geographic raster (terrain.geographic_reach), plus one 1 arc-second
    cell on each side, so the points the builder keeps are all in the file."""
    from .terrain import geographic_reach
    if frame.lat_lon is None:
        raise CaseError("terrain download: needs domain.center_lat_lon (where the domain is on the Earth)")
    dlat, dlon = geographic_reach(frame, margin)
    lat0, lon0 = frame.lat_lon
    pad = 1.0 / 3600.0
    box = (lat0 - dlat - pad, lat0 + dlat + pad, lon0 - dlon - pad, lon0 + dlon + pad)
    if box[0] < -90.0 or box[1] > 90.0:
        raise CaseError("terrain download: the domain reaches latitude %.4f .. %.4f, past a pole"
                        % (box[0], box[1]))
    if box[2] < -180.0 or box[3] > 180.0:
        raise CaseError("terrain download: the domain reaches longitude %.4f .. %.4f, across the 180th meridian "
                        "(not supported)" % (box[2], box[3]))
    return box


def box_for_case(case, given):
    """The download box of a validated case description: its domain plus the builder's margin
    (terrain.margin, or 4 cells of the cell counts the builder will use)."""
    from . import grid
    from .plan import MARGIN_CELLS, _resolve_domain
    problems = []
    frame = _resolve_domain(case, given, problems)
    if problems:
        raise CaseError(problems)
    margin = case["terrain"]["margin"]
    if margin == "auto":
        cell = case["mesh"]["cell_size"]
        if cell is None:
            raise CaseError("terrain download: mesh.cell_size is needed for the margin (terrain.margin = auto)")
        dx = (frame.x_hi - frame.x_lo) / grid.horizontal_cells(frame.x_hi - frame.x_lo, cell)
        dy = (frame.y_hi - frame.y_lo) / grid.horizontal_cells(frame.y_hi - frame.y_lo, cell)
        margin = MARGIN_CELLS * max(dx, dy)
    return box_for_frame(frame, margin)


def read_hgt(data, name="tile"):
    """The heights of one .hgt tile (bytes) as an int16 array, first row north; its size gives the resolution."""
    n = int(round(math.sqrt(len(data) / 2)))
    if n not in TILE_SIZES or 2 * n * n != len(data):
        raise CaseError("terrain download: %s holds %d bytes, not a 3601 x 3601 or 1201 x 1201 SRTM tile"
                        % (name, len(data)))
    return np.frombuffer(data, dtype=">i2").reshape(n, n).astype(np.int16)


def _default_opener(url, timeout):
    return urllib.request.urlopen(url, timeout=timeout)


def fetch_tile(lat, lon, cache_dir, opener=None, timeout=300):
    """The local .hgt file of one tile, downloaded into ``cache_dir`` if it is not there.

    ``opener(url, timeout)`` returns a file-like response (urllib's urlopen by default; the tests pass their
    own). The tile is written to a .part file and renamed only when complete, so a cut download never looks
    cached."""
    name = tile_name(lat, lon)
    path = os.path.join(cache_dir, name + ".hgt")
    if os.path.isfile(path) and os.path.getsize(path) > 0:
        return path                               # a damaged cached tile is refused by read_hgt, by its path
    url = tile_url(lat, lon)
    part = path + ".part"
    try:
        os.makedirs(cache_dir, exist_ok=True)
        response = (opener or _default_opener)(url, timeout)
        try:
            data = gzip.decompress(response.read())
        finally:
            response.close()
        read_hgt(data, name + ".hgt.gz")          # refuse a tile of the wrong size before caching it
        with open(part, "wb") as out:
            out.write(data)
        os.replace(part, path)
    except CaseError:
        raise
    except (urllib.error.URLError, http.client.HTTPException, zlib.error, OSError, EOFError, ValueError) as err:
        raise CaseError("terrain download: %s: %s" % (url, getattr(err, "reason", None) or err)) from None
    finally:
        if os.path.isfile(part):
            os.remove(part)
    return path


class Download:
    """What :func:`download` wrote: the grid file, the attribution file, its extent and notes."""

    def __init__(self):
        self.path = None
        self.sources = None
        self.tiles = []
        self.shape = (0, 0)
        self.lowest = self.highest = None
        self.voids = 0
        self.sea_points = 0
        self.notes = []


def mosaic(box, tile_files):
    """Longitudes and latitudes of the columns and rows (increasing) and the heights (rows south to north, float,
    voids NaN) of the tiles cut to ``box``. ``tile_files`` maps (lat, lon) to a .hgt path; every tile of the box
    is needed, all at one resolution."""
    south, north, west, east = box
    tiles = tiles_for_box(south, north, west, east)
    grids = {}
    for corner in tiles:
        with open(tile_files[corner], "rb") as handle:
            data = handle.read()
        try:
            grids[corner] = read_hgt(data, tile_files[corner])
        except CaseError as err:
            raise CaseError("%s; delete it to fetch the tile again" % err.problems[0]) from None
    sizes = {g.shape[0] for g in grids.values()}
    if len(sizes) != 1:
        raise CaseError("terrain download: the tiles of the box have different resolutions (%s)"
                        % ", ".join("%d" % s for s in sorted(sizes)))
    n = sizes.pop()
    step = n - 1                                  # cells per degree
    lat0, lon0 = tiles[0]
    nlat = max(a for a, _ in tiles) - lat0 + 1
    nlon = max(o for _, o in tiles) - lon0 + 1
    big = np.full((nlat * step + 1, nlon * step + 1), VOID, dtype=np.int16)   # first row south
    for (a, o), g in grids.items():
        r0, c0 = (a - lat0) * step, (o - lon0) * step
        big[r0:r0 + n, c0:c0 + n] = g[::-1, :]    # shared edges hold the same heights; the later tile wins
    # the nodes in the box by index (a box edge on a node keeps that node, whatever its rounding)
    r0 = max(int(math.ceil((south - lat0) * step - 1e-6)), 0)
    r1 = min(int(math.floor((north - lat0) * step + 1e-6)), big.shape[0] - 1)
    c0 = max(int(math.ceil((west - lon0) * step - 1e-6)), 0)
    c1 = min(int(math.floor((east - lon0) * step + 1e-6)), big.shape[1] - 1)
    z = big[r0:r1 + 1, c0:c1 + 1].astype(float)
    z[z == VOID] = np.nan
    return lon0 + np.arange(c0, c1 + 1) / step, lat0 + np.arange(r0, r1 + 1) / step, z


def write_esri_ascii(path, lon, lat, z, nodata=VOID):
    """An ESRI ASCII grid of ``z`` (rows south to north, NaN for no data) at cell centres ``lon``, ``lat``."""
    cell = (lon[-1] - lon[0]) / max(len(lon) - 1, 1)
    with open(path, "w", encoding="utf-8") as out:
        out.write("ncols %d\nnrows %d\nxllcenter %.10f\nyllcenter %.10f\ncellsize %.12g\nNODATA_value %d\n"
                  % (len(lon), len(lat), lon[0], lat[0], cell, nodata))
        for row in z[::-1, :]:                    # the file's first row is the northern one
            out.write(" ".join("%d" % nodata if not np.isfinite(v) else "%.6g" % v for v in row) + "\n")


def download(box, out_path, cache_dir, sea_floor="zero", opener=None, date="(not recorded)", progress=None):
    """Fetch the tiles of ``box`` (south, north, west, east), cut them to it and write ``out_path`` (.asc).

    ``sea_floor``: zero sets heights below 0 m to 0 m, keep leaves them. ``progress(done, total, name)`` is
    called after each tile. Returns a Download."""
    if sea_floor not in SEA_FLOOR:
        raise CaseError("terrain download: sea_floor = '%s': must be zero or keep" % sea_floor)
    if os.path.splitext(out_path)[1].lower() != ".asc":
        raise CaseError("terrain download: %s: the output is an ESRI ASCII grid and needs the extension .asc"
                        % out_path)
    tiles = tiles_for_box(*box)
    files = {}
    for i, (a, o) in enumerate(tiles):
        files[(a, o)] = fetch_tile(a, o, cache_dir, opener=opener)
        if progress is not None:
            progress(i + 1, len(tiles), tile_name(a, o))
    lon, lat, z = mosaic(box, files)
    result = Download()
    result.tiles = [tile_name(a, o) for a, o in tiles]
    result.voids = int(np.sum(~np.isfinite(z)))
    if sea_floor == "zero":
        below = np.isfinite(z) & (z < 0.0)
        result.sea_points = int(np.sum(below))
        z[below] = 0.0
        if result.sea_points:
            result.notes.append("%d of %d points were below 0 m and were set to 0 m (sea floor = zero)"
                                % (result.sea_points, z.size))
    elif np.any(np.isfinite(z) & (z < 0.0)):
        result.notes.append("%d points are below 0 m (sea floor = keep: over the sea these are the sea floor)"
                            % int(np.sum(np.isfinite(z) & (z < 0.0))))
    if result.voids:
        result.notes.append("%d of %d points are voids (written as %d); the builder refuses a domain that holds one"
                            % (result.voids, z.size, VOID))
    result.path = out_path
    result.sources = os.path.splitext(out_path)[0] + ".sources.txt"
    # both files are written beside their names and renamed when complete: a failed write leaves the grid that
    # was there (a terrain.file in use) as it was
    parts = [result.path + ".part", result.sources + ".part"]
    try:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        write_esri_ascii(parts[0], lon, lat, z)
        with open(parts[1], "w", encoding="utf-8") as out:
            out.write(ATTRIBUTION.format(date=date))
            out.write("\nTiles: %s\nBox: latitude %.6f .. %.6f, longitude %.6f .. %.6f\nSea floor: %s\n"
                      % (" ".join(result.tiles), box[0], box[1], box[2], box[3], sea_floor))
        os.replace(parts[0], result.path)
        os.replace(parts[1], result.sources)
    except OSError as err:
        raise CaseError("terrain download: cannot write %s: %s" % (out_path, err.strerror or err)) from None
    finally:
        for part in parts:
            if os.path.isfile(part):
                os.remove(part)
    result.shape = z.shape
    finite = z[np.isfinite(z)]
    if finite.size:
        result.lowest, result.highest = float(finite.min()), float(finite.max())
    return result


def default_cache_dir():
    """Where tiles are kept between downloads: GUST_DEM_CACHE, else ~/.cache/gust/terrain_tiles."""
    return os.environ.get("GUST_DEM_CACHE") or os.path.join(os.path.expanduser("~"), ".cache", "gust", "terrain_tiles")
