#!/usr/bin/env python3
"""Write the terrain and roughness point files of the four Sites_128m cases.

    python3 prepare_sites.py                  all four sites
    python3 prepare_sites.py perdigao         one site (folder name or key)
    python3 prepare_sites.py --cache DIR      tile cache (default GUST_DEM_CACHE or ~/.cache/gust/terrain_tiles)
    python3 prepare_sites.py --faces          only print the inflow-face roughness of the committed files (no network)

For each site the script

1. downloads the terrain (the AWS Terrain Tiles, tools/fetch_terrain.py's download) for the domain and a margin;
2. writes ``terrain.xyz`` (``x y z`` in metres, the UTM zone of the site centre, origin at the centre) through
   tools/make_gust_case.py: points at the spacing of the tiles (about 45 to 49 m) over the domain and a margin,
   heights from the lowest point kept, and for
   the Perdigao case a flat border at the mean level, so its inflow faces are flat like the column of its table;
3. reads the ESA WorldCover 2021 land cover (10 m, a window of the cloud-optimized GeoTIFF on the AWS Open Data
   Registry, through rasterio) and writes ``roughness.xyz`` (``x y z0``): each class becomes a roughness length
   (CLASS_Z0 below) and the 10 m values are averaged in ln z0 over 64 m squares (two per 128 m cell).

The committed point files are the output of this script (README.md, "Data and acknowledgements"). Needs numpy,
PyYAML, rasterio and network access. Exit status 0 when every site was written.
"""
import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "tools"))

from gust_case import geo  # noqa: E402

CELL = 128.0            # the decks' horizontal cell (m)
SPACING = 64.0          # the point spacing of both files (m): two points per cell
MARGIN = 256.0          # points kept beyond each face (m), so every ghost column has points around it

# key: folder, centre (lat, lon), cells per side, border (flat width, slope width) or None
SITES = {
    "butte": ("BigSouthernButte_LogLaw", (43.3960, -113.0220), 80, None),
    "altamont": ("AltamontPass_PowerLaw", (37.7400, -121.6600), 64, None),
    "perdigao": ("Perdigao_TargetWind", (39.7130, -7.7370), 96, (768.0, 1536.0)),
    "flatirons": ("Flatirons_InflowFile", (39.9250, -105.2850), 64, None),
}

# ESA WorldCover classes (Zanaga et al. 2022, product user manual v2.0) and roughness lengths (m), each the class
# of the updated Davenport classification (Wieringa 1992, J. Wind Eng. Ind. Aerodyn. 41, 357-368) that the cover
# describes; trees, built-up, water and bare ground as the ESA CCI table of tools/gust_case/roughness.py.
CLASS_Z0 = {
    10: 1.0,      # tree cover: closed
    20: 0.25,     # shrubland: rough
    30: 0.03,     # grassland: open
    40: 0.1,      # cropland: roughly open
    50: 0.8,      # built-up (ESA CCI 190)
    60: 0.005,    # bare or sparse vegetation (ESA CCI 200)
    70: 0.001,    # snow and ice
    80: 1.0e-4,   # permanent water (ESA CCI 210)
    90: 0.2,      # herbaceous wetland (ESA CCI 180)
    95: 1.0,      # mangroves (ESA CCI 170)
    100: 0.03,    # moss and lichen: open
}
WORLDCOVER = ("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
              "ESA_WorldCover_10m_2021_v200_{tile}_Map.tif")


def worldcover_tile(lat, lon):
    """The 3 x 3 degree tile name that holds (lat, lon): its south-west corner."""
    la = int(math.floor(lat / 3.0) * 3)
    lo = int(math.floor(lon / 3.0) * 3)
    return "%s%02d%s%03d" % ("N" if la >= 0 else "S", abs(la), "E" if lo >= 0 else "W", abs(lo))


def case_yaml(key, centre, cells, border):
    half = 0.5 * cells * CELL
    lines = ["case:", "  name: %s" % key,
             "domain:", "  center_lat_lon: [%.6f, %.6f]" % centre, "  half_width: %.1f" % half, "  top: 3000.",
             "terrain:", "  file: dem.asc", "  point_spacing: %.1f" % SPACING]
    if border is not None:
        lines += ["  border:", "    flat_width: %.1f" % border[0], "    slope_width: %.1f" % border[1],
                  "    level: mean"]
    lines += ["mesh:", "  cell_size: %.1f" % CELL, "  first_layer: 2.", "  stretch_ratio: 1.15",
              "wind:", "  speed: 8.", "  height: 10.", "  direction: 270.",
              "surface:", "  roughness: 0.1"]
    return "\n".join(lines) + "\n"


def write_terrain(key, centre, cells, border, work, out, cache):
    """dem.asc from the tile download, then terrain.xyz through make_gust_case.py."""
    from gust_case import dem
    half = 0.5 * cells * CELL + MARGIN + 2000.0
    dlat = 1.5 * half / 111320.0
    dlon = 1.5 * half / (111320.0 * math.cos(math.radians(centre[0])))
    box = (centre[0] - dlat, centre[0] + dlat, centre[1] - dlon, centre[1] + dlon)
    dem.download(box, os.path.join(work, "dem.asc"), cache, sea_floor="zero", date="2026-11-23")
    with open(os.path.join(work, "case.yaml"), "w") as f:
        f.write(case_yaml(key, centre, cells, border))
    subprocess.run([sys.executable, os.path.join(ROOT, "tools", "make_gust_case.py"), os.path.join(work, "case.yaml"),
                    "--no-check", "--quiet", "--output-dir", os.path.join(work, "case")], check=True)
    found = [os.path.join(d, "terrain.xyz") for d, _, files in os.walk(os.path.join(work, "case"))
             if "terrain.xyz" in files]
    if len(found) != 1:
        raise RuntimeError("%s: make_gust_case.py wrote %d terrain.xyz files" % (key, len(found)))
    shutil.copyfile(found[0], os.path.join(out, "terrain.xyz"))
    shutil.copyfile(os.path.join(work, "dem.sources.txt"), os.path.join(out, "terrain.sources.txt"))
    pts = np.loadtxt(os.path.join(out, "terrain.xyz"))
    return pts


def write_roughness(key, centre, cells, out):
    import rasterio
    from rasterio.windows import from_bounds
    e0, n0, zone, south = geo.latlon_to_utm(np.array([centre[0]]), np.array([centre[1]]))
    e0, n0 = float(np.ravel(e0)[0]), float(np.ravel(n0)[0])
    half = 0.5 * cells * CELL + MARGIN
    dlat = 1.3 * half / 111320.0
    dlon = 1.3 * half / (111320.0 * math.cos(math.radians(centre[0])))
    url = WORLDCOVER.format(tile=worldcover_tile(*centre))
    with rasterio.open(url) as src:
        win = from_bounds(centre[1] - dlon, centre[0] - dlat, centre[1] + dlon, centre[0] + dlat, src.transform)
        codes = src.read(1, window=win)
        tr = src.window_transform(win)
    rows, cols = np.indices(codes.shape)
    lon = tr.c + (cols + 0.5) * tr.a
    lat = tr.f + (rows + 0.5) * tr.e
    unknown = sorted(set(np.unique(codes).tolist()) - set(CLASS_Z0))
    if unknown:
        raise RuntimeError("%s: WorldCover classes without a roughness: %s" % (key, unknown))
    lnz0 = np.vectorize(lambda c: math.log(CLASS_Z0[int(c)]))(codes)
    e, n, _, _ = geo.latlon_to_utm(lat.ravel(), lon.ravel(), zone=zone, south=south)
    x, y = np.asarray(e) - e0, np.asarray(n) - n0
    # 64 m squares centred on the points (i + 1/2) * SPACING - half
    nb = int(round(2 * half / SPACING))
    i = np.floor((x + half) / SPACING).astype(int)
    j = np.floor((y + half) / SPACING).astype(int)
    ok = (i >= 0) & (i < nb) & (j >= 0) & (j < nb)
    s = np.zeros((nb, nb))
    c = np.zeros((nb, nb))
    np.add.at(s, (j[ok], i[ok]), lnz0.ravel()[ok])
    np.add.at(c, (j[ok], i[ok]), 1.0)
    if np.any(c == 0):
        raise RuntimeError("%s: %d squares hold no land-cover pixel" % (key, int(np.sum(c == 0))))
    z0 = np.exp(s / c)
    xc = (np.arange(nb) + 0.5) * SPACING - half
    X, Y = np.meshgrid(xc, xc)
    with open(os.path.join(out, "roughness.xyz"), "w") as f:
        f.write("# x y z0 (m): ESA WorldCover 2021 v200 (tile %s), classes to z0 by prepare_sites.py, ln-mean over "
                "%.0f m squares; frame: UTM zone %d%s, origin at %.6f %.6f\n"
                % (worldcover_tile(*centre), SPACING, zone, "S" if south else "N", centre[0], centre[1]))
        for xv, yv, zv in zip(X.ravel(), Y.ravel(), z0.ravel()):
            f.write("%.1f %.1f %.6g\n" % (xv, yv, zv))
    counts = {k: int(np.sum(codes == k)) for k in CLASS_Z0 if np.any(codes == k)}
    return z0, counts


def face_z0(path, cells):
    """Geometric means of z0 (m) along the xlo and ylo faces: the points of the first 128 m cell inside the face
    ("inner") and the points within 64 m on either side of it ("straddle"); "both" joins the two faces."""
    x, y, z0 = np.loadtxt(path).T
    h = 0.5 * cells * CELL
    inside = lambda v: (v > -h) & (v < h)
    sel = {"inner": lambda v: (v > -h) & (v < -h + CELL),
           "straddle": lambda v: (v > -h - 0.5 * CELL) & (v < -h + 0.5 * CELL)}
    out = {}
    for name, f in sel.items():
        mx, my = f(x) & inside(y), f(y) & inside(x)
        for face, m in (("xlo", mx), ("ylo", my), ("both", mx | my)):
            out[name + " " + face] = math.exp(np.mean(np.log(z0[m])))
    return out


def print_faces(folder, out, cells):
    faces = face_z0(os.path.join(out, "roughness.xyz"), cells)
    print("%s: inflow-face z0 (m) %s" % (folder, ", ".join("%s %.3g" % kv for kv in faces.items())))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("sites", nargs="*", help="site keys or folder names (default: all)")
    p.add_argument("--cache", default=os.environ.get("GUST_DEM_CACHE",
                                                     os.path.expanduser("~/.cache/gust/terrain_tiles")))
    p.add_argument("--faces", action="store_true", help="only print the inflow-face z0 of the committed roughness.xyz")
    a = p.parse_args(argv)
    keys = []
    for s in a.sites or list(SITES):
        k = s if s in SITES else next((kk for kk, v in SITES.items() if v[0] == s), None)
        if k is None:
            p.error("unknown site %s (keys: %s)" % (s, ", ".join(SITES)))
        keys.append(k)
    for key in keys:
        folder, centre, cells, border = SITES[key]
        out = os.path.join(HERE, folder)
        if a.faces:
            print_faces(folder, out, cells)
            continue
        os.makedirs(out, exist_ok=True)
        with tempfile.TemporaryDirectory() as work:
            pts = write_terrain(key, centre, cells, border, work, out, a.cache)
        z0, counts = write_roughness(key, centre, cells, out)
        print("%s: terrain.xyz %d points, x %.0f .. %.0f m, heights %.1f .. %.1f m; roughness.xyz %d points, "
              "z0 %.3g .. %.3g m, geometric mean %.3g m; WorldCover pixels per class %s"
              % (folder, pts.shape[0], pts[:, 0].min(), pts[:, 0].max(), pts[:, 2].min(), pts[:, 2].max(),
                 z0.size, z0.min(), z0.max(), math.exp(np.mean(np.log(z0))), counts))
        print_faces(folder, out, cells)
    return 0


if __name__ == "__main__":
    sys.exit(main())
