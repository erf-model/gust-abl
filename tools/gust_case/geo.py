"""Latitude and longitude to UTM, without pyproj.

The Krueger series in the third flattening n (Karney 2011, J. Geodesy 85,
475-485, eqs. 7-11 and 35), kept to n^3: the truncation is below 1 mm within
a UTM zone. WGS 84.

    zone = floor((lon + 180) / 6) + 1,   lon0 = 6 zone - 183
    E = 500000 + k0 A (eta' + sum_j alpha_j cos(2 j xi') sinh(2 j eta'))
    N = N0     + k0 A (xi'  + sum_j alpha_j sin(2 j xi') cosh(2 j eta'))

with k0 = 0.9996 and N0 = 10^7 m in the southern hemisphere. A domain that
crosses a zone boundary is projected in the zone of its centre (give ``zone``).
The y axis of the projection is grid north; it differs from true north by the
meridian convergence, (lon - lon0) sin(lat), which is not removed.
"""
import math

import numpy as np

WGS84_A = 6378137.0
WGS84_F = 1.0 / 298.257223563
K0 = 0.9996


def utm_zone(lon):
    """The UTM zone (1 to 60) of a longitude in degrees."""
    return int(math.floor((float(lon) + 180.0) / 6.0)) % 60 + 1


def latlon_to_utm(lat, lon, zone=None, south=None):
    """Easting and northing (m) of latitude and longitude (degrees; arrays allowed).

    ``zone`` defaults to that of the mean longitude; ``south`` (a false
    northing of 10^7 m) to the sign of the mean latitude. Returns
    (easting, northing, zone, south).
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    if zone is None:
        zone = utm_zone(float(np.mean(lon)))
    if south is None:
        south = bool(np.mean(lat) < 0.0)
    n = WGS84_F / (2.0 - WGS84_F)
    big_a = WGS84_A / (1.0 + n) * (1.0 + n * n / 4.0 + n ** 4 / 64.0)
    alpha = (n / 2.0 - 2.0 * n * n / 3.0 + 5.0 * n ** 3 / 16.0,
             13.0 * n * n / 48.0 - 3.0 * n ** 3 / 5.0,
             61.0 * n ** 3 / 240.0)
    phi = np.radians(lat)
    dlam = np.radians(lon - (6.0 * zone - 183.0))
    e = 2.0 * math.sqrt(n) / (1.0 + n)
    t = np.sinh(np.arctanh(np.sin(phi)) - e * np.arctanh(e * np.sin(phi)))
    xi = np.arctan2(t, np.cos(dlam))
    eta = np.arctanh(np.sin(dlam) / np.sqrt(1.0 + t * t))
    east = eta.copy()
    north = xi.copy()
    for j, a in enumerate(alpha, start=1):
        east = east + a * np.cos(2 * j * xi) * np.sinh(2 * j * eta)
        north = north + a * np.sin(2 * j * xi) * np.cosh(2 * j * eta)
    east = 500000.0 + K0 * big_a * east
    north = (1.0e7 if south else 0.0) + K0 * big_a * north
    return east, north, zone, south
