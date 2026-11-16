"""Mesh helpers: horizontal cell counts that suit the grid hierarchy, and the
layers of a stretched vertical grid.

GUST's stretched grid (docs/terrain_coordinates.rst) has layers
dz_k = dz0 r^k, k = 0 .. n-1, so the top is

    H(r) = dz0 (r^n - 1) / (r - 1)        (n dz0 for r = 1).

With ``mesh.dz0`` alone GUST finds r so that H = ``mesh.prob_hi`` z; with
``mesh.stretch_ratio`` too, the top is H(r).
"""
import math

COUNT_MULTIPLE = 8   # cell counts are multiples of 8: divisible by 2^(L-1) up to L = 4


def horizontal_cells(length, cell_size, multiple=COUNT_MULTIPLE):
    """The cell count nearest ``length / cell_size`` that is a multiple of ``multiple``."""
    n = int(math.floor(length / cell_size / multiple + 0.5)) * multiple   # a half goes up
    return max(n, multiple)


def hierarchy_levels(nx, ny, max_grid_size, wanted=3, coarsest_min=16):
    """The deepest hierarchy, up to ``wanted`` levels, that the counts allow.

    GUST needs ``mesh.n_cell`` and ``mesh.max_grid_size`` in x and y divisible
    by 2^(L-1) (Tests/inputs_reference, mesh.hierarchy_levels); the helper
    also keeps at least ``coarsest_min`` cells per direction on the coarsest
    level.
    """
    levels = 1
    for candidate in range(2, wanted + 1):
        f = 2 ** (candidate - 1)
        if nx % f or ny % f or max_grid_size % f:
            break
        if nx // f < coarsest_min or ny // f < coarsest_min:
            break
        levels = candidate
    return levels


def stretched_top(dz0, ratio, layers):
    """The height where ``layers`` layers end."""
    if abs(ratio - 1.0) < 1.0e-12:
        return dz0 * layers
    return dz0 * (ratio ** layers - 1.0) / (ratio - 1.0)


def stretch_ratio_for_top(dz0, layers, top):
    """The ratio r >= 1 with stretched_top(dz0, r, layers) = top (bisection)."""
    if layers < 1 or dz0 <= 0.0 or top <= 0.0:
        raise ValueError("dz0, layers and top must be positive")
    if dz0 * layers > top * (1.0 + 1.0e-12):
        raise ValueError("%d layers of %g m pass the top %g m" % (layers, dz0, top))
    if layers == 1 or dz0 * layers >= top * (1.0 - 1.0e-12):
        return 1.0
    lo, hi = 1.0, 2.0
    while stretched_top(dz0, hi, layers) < top:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if stretched_top(dz0, mid, layers) < top:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1.0e-14 * hi:
            break
    return 0.5 * (lo + hi)


def layers_for_ratio(dz0, ratio, top):
    """The fewest layers that reach ``top`` with a ratio of at most ``ratio``."""
    if ratio < 1.0 or dz0 <= 0.0 or top <= 0.0:
        raise ValueError("dz0 and top must be positive and ratio >= 1")
    if abs(ratio - 1.0) < 1.0e-12:
        return int(math.ceil(top / dz0 - 1.0e-9))
    n = math.log(1.0 + top * (ratio - 1.0) / dz0) / math.log(ratio)
    return max(1, int(math.ceil(n - 1.0e-9)))


def layer_edges(dz0, ratio, layers):
    """The layers + 1 level heights, starting at 0."""
    z, edges = 0.0, [0.0]
    for k in range(layers):
        z += dz0 * ratio ** k
        edges.append(z)
    return edges
