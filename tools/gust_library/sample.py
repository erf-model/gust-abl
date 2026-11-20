"""Fields at fixed heights above the ground, and at a point.

    ground(z_nodes)                      the ground under each column's centre: (nx, ny)
    heights_above_ground(z_cc, z_nodes)  each cell centre's height above its column's ground: (nx, ny, nz)
    at_heights(field, agl, heights)      the field at each height above the ground: (len(heights), nx, ny)
    at_point(xc, yc, field2d, x, y)      bilinear between the column centres at (x, y)

The ground under a column's centre is the mean of the four node heights at
the bottom of the column: the bilinear value of the mesh's ground at the
centre. Between two cell centres a field is linear in the height above the
ground. A height below the first cell centre of any column, or above the last,
is refused: below the first centre the profile depends on the wall law, and
the library does not extrapolate.
"""
import math

import numpy as np


class SampleError(ValueError):
    """A height or point the fields cannot be sampled at; the message names the option."""


def ground(z_nodes):
    z0 = z_nodes[:, :, 0]
    return 0.25 * (z0[:-1, :-1] + z0[1:, :-1] + z0[:-1, 1:] + z0[1:, 1:])


def heights_above_ground(z_cc, z_nodes):
    return z_cc - ground(z_nodes)[:, :, None]


def check_heights(agl, heights, option="--heights"):
    """Refuse heights outside every column's cell centres (message names ``option``)."""
    first = float(agl[:, :, 0].max())
    last = float(agl[:, :, -1].min())
    seen = set()
    for h in heights:
        if not (isinstance(h, (int, float)) and math.isfinite(h) and h > 0.0):
            raise SampleError("%s: %r is not a finite height above 0 m" % (option, h))
        if h in seen:
            raise SampleError("%s: %g m is given twice" % (option, h))
        seen.add(h)
        if h < first:
            raise SampleError("%s: %g m is below the first cell centre of some column (the highest first centre is "
                              "%.4g m above its ground; the library does not extrapolate below it)" % (option, h, first))
        if h > last:
            raise SampleError("%s: %g m is above the last cell centre of some column (the lowest is %.6g m above "
                              "its ground)" % (option, h, last))


def at_heights(field, agl, heights):
    """field and agl: (nx, ny, nz); returns (len(heights), nx, ny), linear in the height between centres."""
    check_heights(agl, heights)
    nz = agl.shape[2]
    out = np.empty((len(heights),) + agl.shape[:2])
    for n, h in enumerate(heights):
        k = np.clip((agl <= h).sum(axis=2) - 1, 0, nz - 2)[:, :, None]
        z_lo = np.take_along_axis(agl, k, axis=2)[:, :, 0]
        z_hi = np.take_along_axis(agl, k + 1, axis=2)[:, :, 0]
        f_lo = np.take_along_axis(field, k, axis=2)[:, :, 0]
        f_hi = np.take_along_axis(field, k + 1, axis=2)[:, :, 0]
        t = (h - z_lo) / (z_hi - z_lo)
        out[n] = (1.0 - t) * f_lo + t * f_hi
    return out


def at_point(xc, yc, field2d, x, y, option="--reference-point"):
    """Bilinear between the column centres xc (nx), yc (ny); refused outside the outermost centres."""
    if not (math.isfinite(x) and math.isfinite(y)):
        raise SampleError("%s: (%r, %r) is not a finite point" % (option, x, y))
    if not (xc[0] <= x <= xc[-1] and yc[0] <= y <= yc[-1]):
        raise SampleError("%s: (%g, %g) m is outside the column centres, x %g to %g m and y %g to %g m (case frame)"
                          % (option, x, y, xc[0], xc[-1], yc[0], yc[-1]))
    i = min(int(np.searchsorted(xc, x, side="right")) - 1, len(xc) - 2)
    j = min(int(np.searchsorted(yc, y, side="right")) - 1, len(yc) - 2)
    tx = (x - xc[i]) / (xc[i + 1] - xc[i])
    ty = (y - yc[j]) / (yc[j + 1] - yc[j])
    return ((1 - tx) * (1 - ty) * field2d[..., i, j] + tx * (1 - ty) * field2d[..., i + 1, j]
            + (1 - tx) * ty * field2d[..., i, j + 1] + tx * ty * field2d[..., i + 1, j + 1])
