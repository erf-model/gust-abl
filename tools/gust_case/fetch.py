"""Upstream roughness per wind direction from a regional roughness map (``wind.upstream_roughness = auto``).

The inflow profile of a direction needs one roughness length, that of the ground the wind crossed before it reached
the domain. For a wind from a direction, transects start on the faces the wind enters through and run upwind. Along
each one, ln z0 of the regional map is averaged with the weight exp(-x/L), x the distance upwind of the face and L
``wind.fetch_length``:

    ln z0_up = sum_k w_k ln z0(x_k) / sum_k w_k,    w_k = exp(-x_k/L) - exp(-x_(k+1)/L)

over segments [x_k, x_(k+1)] out to REACH L, with ln z0 taken at each segment's middle. A segment the map does not
cover is left out of both sums; the covered share of the weight is reported, and below MIN_COVER the case is refused.
The transects of a direction are spaced along its inflow faces and weighted by the width of flow they stand for
(face length times the wind's component normal to the face), so the direction's value is the flux-weighted mean of
ln z0 over its inflow faces.

The faces are the deck's inflow faces: those of the wind asked for, or for a geostrophic drive those of the direction
turned by ``turning`` (the middle of the drag law's turn, Builder._deck_faces), whose component normal to each face
also gives the flux weights. The transects run against the wind asked for (wind.direction, the wind near the
ground), whose ground sets z0.

The map is binned on a square grid of its point spacing (exact for a regular raster along x and y), ln z0
averaged in each bin, and two passes fill an empty bin from its filled neighbours (a raster turned into the case frame leaves a few empty bins). A sample takes
the bin it falls in; outside the binned area it is not covered.
"""
import math

import numpy as np

from . import wind
from .errors import CaseError

FETCH_PER_HEIGHT = 100.0       # default wind.fetch_length: 100 times wind.height
REACH = 5.0                    # transects run to 5 L (the weight beyond is exp(-5), 0.7 %)
MIN_COVER = 1.0 - math.exp(-1.0)   # refused below: the map covers less than one fetch length on average
NOTE_COVER = 0.95              # noted below
FILL_PASSES = 2
MAX_POINTS_PER_FACE = 400
MIN_POINTS_PER_FACE = 8
MAX_SEGMENTS = 2000
SEGMENTS_PER_FETCH = 20


class BinnedMap:
    """ln z0 of a point map, binned on a square grid of about its point spacing."""

    def __init__(self, x, y, z0, key):
        x, y, z0 = (np.asarray(v, dtype=float).ravel() for v in (x, y, z0))
        if x.size < 4:
            raise CaseError("%s.file: %d points; the map needs at least 4 to be binned" % (key, x.size))
        self.x_lo, self.y_lo = float(x.min()), float(y.min())
        width, height = float(x.max()) - self.x_lo, float(y.max()) - self.y_lo
        if not (width > 0.0 and height > 0.0):
            raise CaseError("%s.file: the points lie on a line (x span %g m, y span %g m); the map needs an area"
                            % (key, width, height))
        # (W + h)(H + h) = N h^2: exactly the spacing of a regular nx x ny raster (W = (nx - 1) h, N = nx ny)
        n = x.size
        self.h = ((width + height) + math.sqrt((width + height) ** 2 + 4.0 * (n - 1) * width * height)) / (2.0 * (n - 1))
        i = _nearest(x - self.x_lo, self.h)
        j = _nearest(y - self.y_lo, self.h)
        self.nx, self.ny = int(i.max()) + 1, int(j.max()) + 1
        total = np.zeros((self.ny, self.nx))
        count = np.zeros((self.ny, self.nx))
        np.add.at(total, (j, i), np.log(z0))
        np.add.at(count, (j, i), 1.0)
        grid = np.full((self.ny, self.nx), np.nan)
        filled = count > 0.0
        grid[filled] = total[filled] / count[filled]
        for _ in range(FILL_PASSES):
            grid = _fill_once(grid)
        self.ln_z0 = grid

    def sample(self, x, y):
        """ln z0 at the points (x, y); NaN where the map does not reach."""
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        fi, fj = _nearest(x - self.x_lo, self.h), _nearest(y - self.y_lo, self.h)
        inside = (fi >= 0) & (fi < self.nx) & (fj >= 0) & (fj < self.ny)
        out = np.full(x.shape, np.nan)
        out[inside] = self.ln_z0[fj[inside], fi[inside]]
        return out

    def sample_nearest_edge(self, x, y):
        """ln z0 at the points, a point beyond the map taking the nearest edge bin (NaN where that bin is empty).
        For the ground along a face: the inner map may start just inside the face."""
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        fi = np.clip(_nearest(x - self.x_lo, self.h), 0, self.nx - 1)
        fj = np.clip(_nearest(y - self.y_lo, self.h), 0, self.ny - 1)
        return self.ln_z0[fj, fi]

    def within(self, x, y):
        """True where (x, y) lies within the map's extent (its bins), filled or not."""
        x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
        fi, fj = _nearest(x - self.x_lo, self.h), _nearest(y - self.y_lo, self.h)
        return (fi >= 0) & (fi < self.nx) & (fj >= 0) & (fj < self.ny)


def _nearest(offset, h):
    """Index of the nearest bin centre; a point half-way between two goes up (the same at both edges)."""
    return np.floor(np.asarray(offset, dtype=float) / h + 0.5).astype(np.int64)


def _fill_once(grid):
    """Each empty bin with a filled 8-neighbour takes the mean of its filled 8-neighbours."""
    padded = np.pad(grid, 1, constant_values=np.nan)
    total = np.zeros(grid.shape)
    count = np.zeros(grid.shape)
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            block = padded[1 + dj:1 + dj + grid.shape[0], 1 + di:1 + di + grid.shape[1]]
            have = np.isfinite(block)
            total[have] += block[have]
            count[have] += 1.0
    out = grid.copy()
    empty = ~np.isfinite(grid) & (count > 0.0)
    out[empty] = total[empty] / count[empty]
    return out


def segment_weights(reach, fetch_length, segments):
    """Edges, middles and exp(-x/L) weights of ``segments`` equal segments over [0, reach]."""
    edges = np.linspace(0.0, reach, segments + 1)
    weights = np.exp(-edges[:-1] / fetch_length) - np.exp(-edges[1:] / fetch_length)
    return edges, 0.5 * (edges[:-1] + edges[1:]), weights


class FaceResult:
    def __init__(self, face, z0_upwind, z0_ground, cover):
        self.face = face
        self.z0_upwind = z0_upwind      # fetch-weighted z0 of this face's transects (m), or None if uncovered
        self.z0_ground = z0_ground      # geometric mean of the ground's z0 along the face (m), or None if unknown
        self.cover = cover              # covered share of this face's weight

    def mismatch(self):
        """ln(z0 upwind / z0 of the ground at the face), or None."""
        if self.z0_upwind is None or self.z0_ground is None:
            return None
        return math.log(self.z0_upwind / self.z0_ground)


class DirectionResult:
    def __init__(self, direction, z0, cover, faces, fetch_length, turning=0.0):
        self.direction = direction
        self.turning = turning
        self.z0 = z0
        self.cover = cover
        self.faces = faces
        self.fetch_length = fetch_length


def _face_geometry(frame, face):
    """(start point, end point, outward normal) of a lateral face of the domain."""
    x0, x1, y0, y1 = frame.x_lo, frame.x_hi, frame.y_lo, frame.y_hi
    return {"xlo": ((x0, y0), (x0, y1), (-1.0, 0.0)), "xhi": ((x1, y0), (x1, y1), (1.0, 0.0)),
            "ylo": ((x0, y0), (x1, y0), (0.0, -1.0)), "yhi": ((x0, y1), (x1, y1), (0.0, 1.0))}[face]


def upstream_roughness(regional, frame, direction, fetch_length, ground, key="surface.regional_roughness_map",
                       turning=0.0, side_faces="slip_wall"):
    """The DirectionResult of a wind from ``direction`` (degrees, meteorological).

    ``regional`` is the BinnedMap of the regional map; ``ground`` gives ln z0 of the ground inside the domain at
    points (a BinnedMap's sample, or a function of a uniform value), for the mismatch at each face. ``turning`` and
    ``side_faces`` choose the inflow faces as wind.lateral_faces does for the deck."""
    if not (math.isfinite(fetch_length) and fetch_length > 0.0):
        raise ValueError("fetch_length must be finite and > 0")
    towards = wind.towards_from_meteorological(direction)
    downwind = np.array(wind.unit_vector(towards))
    flux = np.array(wind.unit_vector(towards + turning))
    faces = wind.lateral_faces(towards, side_faces, turning)
    reach = REACH * fetch_length
    segments = int(min(MAX_SEGMENTS, max(SEGMENTS_PER_FETCH * REACH, math.ceil(2.0 * reach / regional.h))))
    _, middles, weights = segment_weights(reach, fetch_length, segments)
    total_weight = float(weights.sum())
    face_results = []
    num = den = width_total = uncovered_within = 0.0
    for face in ("xlo", "xhi", "ylo", "yhi"):
        if faces[face] != "inflow":
            continue
        (ax, ay), (bx, by), normal = _face_geometry(frame, face)
        length = math.hypot(bx - ax, by - ay)
        normal_speed = -float(np.dot(flux, normal))             # > 0 on an inflow face
        n = int(min(MAX_POINTS_PER_FACE, max(MIN_POINTS_PER_FACE, math.ceil(length / regional.h))))
        t = (np.arange(n) + 0.5) / n
        px, py = ax + t * (bx - ax), ay + t * (by - ay)
        sx = px[:, None] - downwind[0] * middles[None, :]
        sy = py[:, None] - downwind[1] * middles[None, :]
        ln_z0 = regional.sample(sx, sy)
        covered = np.isfinite(ln_z0)
        w = np.where(covered, weights[None, :], 0.0)
        face_num = float(np.sum(w * np.where(covered, ln_z0, 0.0)))
        face_den = float(np.sum(w))
        width = length * normal_speed
        missing = ~covered & regional.within(sx, sy)
        uncovered_within += width * float(np.sum(np.where(missing, weights[None, :], 0.0))) / n
        ground_ln = np.asarray(ground(px, py), dtype=float)
        z0_ground = math.exp(float(np.mean(ground_ln))) if np.all(np.isfinite(ground_ln)) else None
        face_results.append(FaceResult(face, math.exp(face_num / face_den) if face_den > 0.0 else None,
                                       z0_ground, face_den / (n * total_weight)))
        num += width * face_num / n
        den += width * face_den / n
        width_total += width
    cover = den / (width_total * total_weight)
    if cover < MIN_COVER and uncovered_within > 0.5 * (width_total * total_weight - den):
        raise CaseError("%s.file covers %.0f %% of the fetch weight upwind of the wind from %g degrees (at least "
                        "%.0f %% is needed), and most of the rest lies within the map's extent, where its points are "
                        "further apart than its mean spacing (%.4g m, the bin size) and leave empty bins; give the "
                        "map at one spacing (for example resampled to the coarsest)"
                        % (key, 100.0 * cover, direction, 100.0 * MIN_COVER, regional.h))
    if cover < MIN_COVER:
        raise CaseError("%s.file covers %.0f %% of the fetch weight upwind of the wind from %g degrees (at least "
                        "%.0f %%, one wind.fetch_length = %g m on average, is needed); give a larger map or a "
                        "shorter wind.fetch_length" % (key, 100.0 * cover, direction, 100.0 * MIN_COVER, fetch_length))
    return DirectionResult(direction, math.exp(num / den), cover, face_results, fetch_length, turning)


def uniform_ground(z0):
    """``ground`` for upstream_roughness: a uniform roughness length."""
    value = math.log(z0)
    return lambda x, y: np.full(np.shape(x), value)


def report_lines(results, ustars=None):
    """The Upstream roughness lines of case_summary.txt; ``ustars`` maps a direction to its log-law u* (or None)."""
    if not results:
        return []
    first = results[0]
    lines = ["Upstream wind.upstream_roughness = auto: fetch-weighted ln z0 of %s, weight exp(-x/L), "
             "L = %.6g m" % ("surface.regional_roughness_map", first.fetch_length),
             "        ln mismatch = ln(z0 upwind of the face / z0 of the ground along the face)"]
    for r in results:
        ustar = (ustars or {}).get(r.direction)
        lines.append("        wind from %g: z0 = %.6g m, map covers %.1f %% of the weight%s%s"
                     % (r.direction, r.z0, 100.0 * r.cover, "" if ustar is None else ", log law u* = %.6g m/s" % ustar,
                        "" if r.turning == 0.0 else ", faces turned %+.2f degrees" % r.turning))
        for f in r.faces:
            ground = "not known" if f.z0_ground is None else "%.6g m" % f.z0_ground
            if f.z0_upwind is None:
                lines.append("          %s: not covered upwind; ground %s" % (f.face, ground))
            elif f.z0_ground is None:
                lines.append("          %s: upwind %.6g m, ground not known along the face" % (f.face, f.z0_upwind))
            else:
                lines.append("          %s: upwind %.6g m, ground %.6g m, ln mismatch %+.3f"
                             % (f.face, f.z0_upwind, f.z0_ground, f.mismatch()))
            if f.cover < MIN_COVER:
                lines.append("            %s: the map covers %.0f %% of this face's weight" % (f.face, 100.0 * f.cover))
    return lines
