"""Query a wind library: the wind for any direction and speed at a stored height.

    lib = Library(case_dir or library_dir)
    w = lib.wind(direction, speed, height)     dict: u, v, w, k (nx, ny), x, y, ground, weights

The two sectors either side of ``direction`` (degrees the wind comes from,
clockwise from north) are blended linearly in their stored u, v, w and k, with
weights in proportion to the angle (the vector, not the speed and direction).
The blended u, v, w are then divided by the length of the blended reference
wind and multiplied by ``speed``, so that the wind at the reference point and
first height has the speed ``speed``: two unit vectors an angle apart blend to
a vector cos(half the angle) long, which would otherwise take 0.4 % off a
10-degree library midway and 1.5 % off a 20-degree one. That shortfall comes
from the two directions partly cancelling, not from a weaker flow, so k (a
scalar, stored per unit reference speed squared) is multiplied by speed**2
only. A direction within
1e-6 degrees of a sector takes that sector alone. ``direction`` is the
sweep's inflow direction (the weights use the sectors' nominal directions);
``speed`` is the local wind speed at the library's reference point and first
height, not the inflow speed (over a crest the two differ).

Refused: a height the library did not store; a sector that did not reach
GUST's steady verdict, when the blend uses it; two neighbours 180 degrees or
more apart (a library that does not cover the direction); a library of one
sector away from its direction.

The blend is a linear interpolation: a feature present in one sector and not
the next, such as a lee recirculation, is smeared between them. Scaling by the
speed is exact for a neutral flow over rough ground driven by an inflow
profile (the solution is then independent of the speed, molecular viscosity
aside); with a geostrophic drive the turning angle depends on the speed
(Rossby number), so a library of that kind holds at the speed it was run at.
"""
import json
import math
import numbers
import os

import numpy as np

from .index import FOLDER, FORMAT_VERSION, INDEX, LibraryError

SAME_DIRECTION = 1.0e-6   # degrees
HEIGHT_TOL = 1.0e-9       # m


class Library:
    def __init__(self, path):
        """path: the case folder or its wind_library/ folder."""
        folder = path if os.path.isfile(os.path.join(path, INDEX)) else os.path.join(path, FOLDER)
        index = os.path.join(folder, INDEX)
        if not os.path.isfile(index):
            raise LibraryError("%s: not found; build the library first (make_gust_library.py build)" % index)
        try:
            with open(index, encoding="utf-8") as fh:
                self.index = json.load(fh)
            if self.index.get("format_version") != FORMAT_VERSION:
                raise LibraryError("%s: format_version %r; this tool reads %d (build the library again)"
                                   % (index, self.index.get("format_version"), FORMAT_VERSION))
            self.folder = folder
            self.sectors = sorted(self.index["sectors"], key=lambda s: s["direction"])
            with np.load(os.path.join(folder, self.index["grid"])) as grid:
                self.x, self.y, self.ground = grid["x"], grid["y"], grid["ground"]
        except (ValueError, KeyError, TypeError, OSError) as err:
            if isinstance(err, LibraryError):
                raise
            raise LibraryError("%s: not a library index written by make_gust_library.py (%s: %s)"
                               % (index, type(err).__name__, err))
        if not self.sectors:
            raise LibraryError("%s: no sectors" % index)
        self._cache = {}

    @property
    def heights(self):
        return list(self.index["heights"])

    def _fields(self, s):
        if s["file"] not in self._cache:
            path = os.path.join(self.folder, s["file"])
            try:
                with np.load(path) as data:
                    fields = {f: data[f] for f in ("u", "v", "w", "k")}
            except (OSError, KeyError, ValueError) as err:
                raise LibraryError("%s: cannot read the sector (%s: %s)" % (path, type(err).__name__, err))
            shape = (len(self.heights), len(self.x), len(self.y))
            if any(a.shape != shape for a in fields.values()):
                raise LibraryError("%s: arrays of shape %s, the index asks for %s (build the library again)"
                                   % (path, fields["u"].shape, shape))
            self._cache[s["file"]] = fields
        return self._cache[s["file"]]

    def weights(self, direction):
        """[(sector entry, weight)] for a direction (degrees from); refuses a direction the library does not cover."""
        if not (_real(direction) and math.isfinite(direction)):
            raise LibraryError("--direction: %r is not a finite angle" % (direction,))
        d = float(direction) % 360.0
        for s in self.sectors:
            gap = abs((d - s["direction"] + 180.0) % 360.0 - 180.0)
            if gap <= SAME_DIRECTION:
                return [(s, 1.0)]
        if len(self.sectors) == 1:
            raise LibraryError("--direction: %g; the library has one sector, %g degrees" % (d, self.sectors[0]["direction"]))
        dirs = [s["direction"] for s in self.sectors]
        n = int(np.searchsorted(dirs, d)) % len(dirs)       # the first sector clockwise of d (wraps past 360)
        a, b = self.sectors[n - 1], self.sectors[n]
        width = (b["direction"] - a["direction"]) % 360.0
        if width >= 180.0:
            raise LibraryError("--direction: %g is between the sectors %g and %g degrees, %g degrees apart; the library "
                               "blends neighbours less than 180 degrees apart" % (d, a["direction"], b["direction"], width))
        t = ((d - a["direction"]) % 360.0) / width
        return [(a, 1.0 - t), (b, t)]

    def wind(self, direction, speed, height):
        """The wind for ``direction`` (degrees from) with ``speed`` (m/s) at the reference point and the first
        height, at ``height`` (m above the ground, one of the library's): dict u, v, w, k (nx, ny), x, y,
        ground and weights [(direction, weight)]."""
        if not (_real(speed) and math.isfinite(speed) and speed > 0.0):
            raise LibraryError("--speed: %r is not a finite speed above 0 m/s" % (speed,))
        matches = [n for n, h in enumerate(self.heights) if _real(height) and abs(h - height) <= HEIGHT_TOL]
        if not matches:
            raise LibraryError("--height: %r m is not a height of the library (%s m)"
                               % (height, ", ".join("%g" % h for h in self.heights)))
        speed = float(speed)
        level = matches[0]
        used = self.weights(direction)
        bad = [s for s, w in used if not s["converged"]]
        if bad:
            raise LibraryError("--direction %g: the sector %s (%s) did not reach GUST's steady verdict (%s); "
                               "it is not used" % (direction % 360.0, bad[0]["folder"], "%g degrees" % bad[0]["direction"],
                                                   bad[0]["status"]))
        # Each sector's reference wind is a unit vector at the reference point; the blend of two is shorter
        # (cos of half the angle between them), so the blend is divided by its length: the reference speed
        # of the answer is ``speed``. The length depends only on the angle between the two, so the
        # meteorological angles serve as they are.
        ex = sum(w * math.cos(math.radians(s["reference_direction"])) for s, w in used)
        ey = sum(w * math.sin(math.radians(s["reference_direction"])) for s, w in used)
        scale = speed / math.hypot(ex, ey)
        out = {}
        for f in ("u", "v", "w"):
            out[f] = sum(w * self._fields(s)[f][level] for s, w in used) * scale
        out["k"] = sum(w * self._fields(s)["k"][level] for s, w in used) * speed ** 2
        out.update(x=self.x, y=self.y, ground=self.ground,
                   weights=[(s["direction"], w) for s, w in used])
        return out


def _real(x):
    """A real number (numpy scalars included), not a bool."""
    return isinstance(x, numbers.Real) and not isinstance(x, (bool, np.bool_))
