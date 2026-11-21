"""Statistics of the wind over a distribution of directions and speeds, from a wind library (no new runs).

    members = from_file(path)  or  parse(lines)  or  product(directions, direction_weights, speeds, speed_weights)
    stats = ensemble(lib, members, height, percentiles=(10, 50, 90), thresholds=(5,))
    stats = at_points(lib, members, height, points, percentiles, thresholds)
    values = wind_at_points(lib, direction, speed, height, points)     one wind (Library.wind) at points

A member is (direction, speed, weight): the direction the wind comes from (degrees, clockwise from north; the
sweep's inflow direction, as in Library.wind) and the speed at the library's reference point and first height
(m/s), with a weight. The weights are divided by their sum. Each member's wind is Library.wind(direction, speed,
height); the wind of a member is its direction's wind for a speed of 1 m/s times its speed, so each direction
is blended once.

Per column (or per point) at the height, with S_i the member's local speed, (u_i, v_i) its wind and w_i its
normalised weight:

    mean_speed      sum w_i S_i
    std_speed       sqrt(sum w_i (S_i - mean_speed)^2)            (the spread; a population standard deviation)
    mean_u, mean_v  sum w_i u_i, sum w_i v_i                       (the vector mean)
    mean_direction  the direction the vector mean comes from (degrees from north); NaN where its length is
                    at most 1e-12 of mean_speed (directions that cancel)
    p<P>            the P-th weighted percentile of S_i: the smallest S_i of a member with a weight above 0 whose
                    cumulative weight, in order of speed, reaches P / 100 (no interpolation between members)
    pexceed_<T>     sum of w_i over the members with S_i > T (m/s)
    preversal       sum of w_i over the members whose local wind opposes their own inflow direction:
                    (u_i, v_i) . e_i < 0, with e_i the unit vector towards which the member's inflow blows

At points, the statistics are taken over the members' winds at the point (each interpolated bilinearly between
the column centres), not interpolated from the column statistics. The columns are taken a few rows at a time, so
the memory held is about CHUNK_VALUES values per member array, whatever the number of members.
"""
import math
import os
import re

import numpy as np

from .grids import direction_from
from .index import LibraryError
from .query import _real
from .sample import at_point

CUMULATIVE_TOL = 1.0e-12   # the cumulative weight reaches P / 100 within this, or within 4 n eps for n members
                           # (a sum of n weights rounds by up to about n eps)
CALM_FRACTION = 1.0e-12    # the vector mean is calm (no direction) at or below this fraction of the mean speed
CHUNK_VALUES = 4_000_000   # members x columns held at once by ensemble()


def _number(value, what):
    if not (_real(value) and math.isfinite(value)):
        raise LibraryError("%s: %r is not a finite number" % (what, value))
    return float(value)


def normalise(members, source="members"):
    """[(direction, speed, weight)] with the weights divided by their sum; refuses a negative or non-finite weight,
    a speed not above 0 and weights that sum to 0."""
    members = list(members)
    if not members:
        raise LibraryError("%s: no members" % source)
    out = []
    for n, m in enumerate(members):
        if len(m) != 3:
            raise LibraryError("%s: member %d is %r, not (direction, speed, weight)" % (source, n + 1, m))
        d = _number(m[0], "%s: direction of member %d" % (source, n + 1))
        s = _number(m[1], "%s: speed of member %d" % (source, n + 1))
        w = _number(m[2], "%s: weight of member %d" % (source, n + 1))
        if s <= 0.0:
            raise LibraryError("%s: speed of member %d is %g m/s; it must be above 0" % (source, n + 1, s))
        if w < 0.0:
            raise LibraryError("%s: weight of member %d is %g; weights must be 0 or more" % (source, n + 1, w))
        d = d % 360.0
        out.append((0.0 if d >= 360.0 else d, s, w))          # -1e-20 % 360 is 360.0
    total = sum(m[2] for m in out)
    if not total > 0.0:
        raise LibraryError("%s: the weights sum to 0" % source)
    return [(d, s, w / total) for d, s, w in out]


def product(directions, direction_weights, speeds, speed_weights):
    """Members for every direction with every speed, weight = direction weight x speed weight (independent
    distributions); a list of weights left out (None) is uniform."""
    def weights(values, given, what):
        if given is None:
            return [1.0] * len(values)
        if len(given) != len(values):
            raise LibraryError("--%s-weights: %d weights for %d %ss" % (what, len(given), len(values), what))
        return list(given)
    if not directions:
        raise LibraryError("--directions: no directions")
    if not speeds:
        raise LibraryError("--speeds: no speeds")
    wd = weights(directions, direction_weights, "direction")
    ws = weights(speeds, speed_weights, "speed")
    for name, ws_ in (("--direction-weights", wd), ("--speed-weights", ws)):
        for w in ws_:
            if not (_real(w) and math.isfinite(w) and w >= 0.0):
                raise LibraryError("%s: %r is not a finite weight of 0 or more" % (name, w))
        if not sum(ws_) > 0.0:
            raise LibraryError("%s: the weights sum to 0" % name)
    return normalise([(d, s, a * b) for d, a in zip(directions, wd) for s, b in zip(speeds, ws)],
                     "--directions and --speeds")


def from_file(path):
    """Members from a text file of the form parse() reads."""
    try:
        with open(path, encoding="utf-8") as fh:
            lines = fh.readlines()
    except OSError as err:
        raise LibraryError("--members: cannot read %s (%s)" % (path, err))
    return parse(lines, os.path.basename(path))


def parse(lines, source="members"):
    """Members from lines of text: one member per line, ``direction speed [weight]`` (weight 1 when left out),
    separated by spaces or commas; text after # is a comment. The weights are normalised."""
    members = []
    for n, line in enumerate(lines, 1):
        words = line.split("#", 1)[0].replace(",", " ").split()
        if not words:
            continue
        if len(words) not in (2, 3):
            raise LibraryError("%s:%d: %r; a member is 'direction speed [weight]'" % (source, n, line.rstrip("\n")))
        try:
            values = [float(w) for w in words]
        except ValueError:
            raise LibraryError("%s:%d: %r is not numbers" % (source, n, line.rstrip("\n")))
        members.append(tuple(values) if len(values) == 3 else (values[0], values[1], 1.0))
    return normalise(members, source)


def _towards_unit(direction):
    """(east, north) components of the unit vector towards which a wind from ``direction`` blows."""
    a = math.radians(direction)
    return -math.sin(a), -math.cos(a)


def statistics(u, v, members, percentiles=(), thresholds=()):
    """The statistics of the docstring from the members' winds u, v (shape (members, ...)) and their normalised
    weights; dict of arrays of shape u.shape[1:]."""
    w = np.array([m[2] for m in members]).reshape((-1,) + (1,) * (u.ndim - 1))
    speed = np.hypot(u, v)
    out = {}
    out["mean_speed"] = (w * speed).sum(axis=0)
    out["std_speed"] = np.sqrt(np.maximum((w * (speed - out["mean_speed"]) ** 2).sum(axis=0), 0.0))
    out["mean_u"] = (w * u).sum(axis=0)
    out["mean_v"] = (w * v).sum(axis=0)
    out["mean_direction"] = direction_from(out["mean_u"], out["mean_v"], CALM_FRACTION * out["mean_speed"])
    if percentiles:
        order = np.argsort(speed, axis=0, kind="stable")
        ranked = np.take_along_axis(speed, order, axis=0)
        wr = np.take_along_axis(np.broadcast_to(w, speed.shape), order, axis=0)
        cum = np.cumsum(wr, axis=0)
        tol = max(CUMULATIVE_TOL, 4.0 * len(members) * np.finfo(float).eps)
        for p in percentiles:
            # p / 100 of the sum as summed (not of 1): P = 100 is the last member with a weight above 0
            first = np.argmax((cum >= p / 100.0 * cum[-1] - tol) & (wr > 0.0), axis=0)
            out[percentile_name(p)] = np.take_along_axis(ranked, first[None], axis=0)[0]
    for t in thresholds:
        out[exceed_name(t)] = (w * (speed > t)).sum(axis=0)
    e = np.array([_towards_unit(m[0]) for m in members])
    shape = (-1,) + (1,) * (u.ndim - 1)
    against = u * e[:, 0].reshape(shape) + v * e[:, 1].reshape(shape) < 0.0
    out["preversal"] = (w * against).sum(axis=0)
    return out


def percentile_name(p):
    return "p%s" % ("%.12g" % p).replace(".", "p")


def exceed_name(t):
    return "pexceed_%s" % ("%.12g" % t).replace(".", "p")


def speed_names(percentiles):
    """The statistics that are speeds (m/s), besides mean_speed: converted with --speed-units in the grids."""
    return ["std_speed"] + [percentile_name(p) for p in percentiles]


def check_options(percentiles, thresholds):
    for p in percentiles:
        if not (_real(p) and math.isfinite(p) and 0.0 < p <= 100.0):
            raise LibraryError("--percentiles: %r is not in (0, 100]" % (p,))
    for t in thresholds:
        if not (_real(t) and math.isfinite(t) and t >= 0.0):
            raise LibraryError("--thresholds: %r is not a finite speed of 0 m/s or more" % (t,))
    names = [percentile_name(p) for p in percentiles] + [exceed_name(t) for t in thresholds]
    if len(set(names)) != len(names):
        raise LibraryError("--percentiles / --thresholds: a value is given twice (names are kept to 12 digits)")


def _unit_winds(lib, members, height, source):
    """{direction: (u, v)} of each direction of the members for 1 m/s at the reference point; a refusal of the
    library names the member."""
    unit = {}
    for n, (d, _, _) in enumerate(members):
        if d in unit:
            continue
        try:
            w = lib.wind(d, 1.0, height)
        except LibraryError as err:
            text = re.sub(r"^--direction(?: \S+)?: ", "", str(err))      # the member names the direction
            raise LibraryError("%s: member %d (from %g degrees): %s" % (source, n + 1, d, text))
        unit[d] = (w["u"], w["v"])
    return unit


def _member_winds(unit, members, rows=slice(None)):
    """u, v of every member in the rows of x, shape (members, rows, ny)."""
    u = np.array([s * unit[d][0][rows] for d, s, _ in members])
    v = np.array([s * unit[d][1][rows] for d, s, _ in members])
    return u, v


def ensemble(lib, members, height, percentiles=(10.0, 50.0, 90.0), thresholds=(), source="members"):
    """Column statistics at ``height`` (one of the library's): dict of (nx, ny) arrays, plus x, y."""
    check_options(percentiles, thresholds)
    members = normalise(members, source)
    unit = _unit_winds(lib, members, height, source)
    nx, ny = len(lib.x), len(lib.y)
    step = max(1, CHUNK_VALUES // (len(members) * ny))
    out = {}
    for i0 in range(0, nx, step):
        rows = slice(i0, min(nx, i0 + step))
        part = statistics(*_member_winds(unit, members, rows), members, percentiles, thresholds)
        for name, a in part.items():
            out.setdefault(name, np.empty((nx, ny)))[rows] = a
    out.update(x=lib.x, y=lib.y)
    return out


def at_points(lib, members, height, points, percentiles=(10.0, 50.0, 90.0), thresholds=(), source="members"):
    """The statistics at points [(x, y)] (case frame, m): dict of arrays of length len(points)."""
    check_options(percentiles, thresholds)
    members = normalise(members, source)
    unit = _unit_winds(lib, members, height, source)
    # each direction's wind at the points, then scaled per member (interpolation is linear in the speed)
    at = {d: tuple(np.array([at_point(lib.x, lib.y, f, x, y, "--points") for x, y in points]) for f in uv)
          for d, uv in unit.items()}
    pu = np.array([s * at[d][0] for d, s, _ in members])
    pv = np.array([s * at[d][1] for d, s, _ in members])
    return statistics(pu, pv, members, percentiles, thresholds)


def wind_at_points(lib, direction, speed, height, points):
    """Library.wind at points [(x, y)] (case frame, m): dict u, v, w, k, speed, direction (degrees from)."""
    wind = lib.wind(direction, speed, height)
    out = {f: np.array([at_point(lib.x, lib.y, wind[f], x, y, "--points") for x, y in points])
           for f in ("u", "v", "w", "k")}
    out["speed"] = np.hypot(out["u"], out["v"])
    out["direction"] = direction_from(out["u"], out["v"])
    return out
