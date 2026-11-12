#!/usr/bin/env python3
"""Digitize the wind-tunnel measurements around the cube from Figs. A3 and A4 of Wise et al. (2026).

    digitize_wise2026.py FIG_A3.png FIG_A4.png OUT.csv

The figures are the raster images of pages 37 and 38 of the preprint (egusphere-2026-4030, CC BY 4.0), extracted
with `pdfimages -png -f 37 -l 38 egusphere-2026-4030.pdf fig` (fig-000.png and fig-003.png, 2400 x 2400 px). The
black dots are the measurements of Snyder and Lawson (1994); README.md has the sources.

The figures' scales, read from their tick marks and arrows:
- x: x = -0.3 + (col - 474) / 1400 m, in both panels (a tick every 140 px per 0.1 m), x from the cube's upwind face.
- Panel (a), the horizontal plane z = 0.02 m: plotted Y = 0.8 - (row - 288) / 700. The panel has equal aspect
  (1400 px per m of x) while its Y labels step by 0.2 per 140 px, and the cube (H = 0.2 m) is drawn from Y = -0.2
  to 0.2 and the refined region of Fig. A4, stated as y in (-0.22, 0.22) m, from -0.44 to 0.44: the true lateral
  position is y = Y / 2.
- Panel (b), the vertical plane y = 0: z = 0.4 - (row - 1576) / 1400 m.
- Each profile is drawn from its station's grey line at x_s as x = x_s + 0.025 u (the 2 m/s arrow is 70 px long).

Dots are found by eroding the dark pixels with a disk of radius 5 px: lines, text and ticks vanish, a dot leaves a
smaller full disk, whose centroid is the dot's centre. A dot cut by the frame (the top and bottom rows) takes the
frame's coordinate and its x from the centroid. Dots of one level go to the stations in order. The values are
those of Fig. A3, with three from Fig. A4 that A3's legend hides (y = 0.25 m at the three upwind stations), and
the dot at x = -0.3 m, z = 0.2 m, which the 3 m/s arrow's head covers, from its right edge and the dots' radius.
Fig. A4's refined-region box crosses the dots at x = 0.5 m, so A3 is taken there. Away from that station the two
figures agree to 0.003 m/s rms over the dots they share; a pixel is 0.029 m/s.
"""
import csv
import sys

import numpy as np
from PIL import Image

PX = 1400.0
STATIONS = {'a': [-0.3, -0.2, -0.1, 0.04, 0.16, 0.3, 0.4, 0.5, 0.7],
            'b': [-0.3, -0.2, -0.1, 0.06, 0.14, 0.3, 0.4, 0.5, 0.7]}
SCALE = 0.025
FRAME = {'a': (288.0, 1408.0), 'b': (1576.0, 2136.0)}


def xdata(c):
    return -0.3 + (c - 474.0) / PX


def label(mask):
    """4-connected components of a boolean mask: a list of (rows, cols)."""
    H, W = mask.shape
    lab = np.zeros(mask.shape, dtype=np.int32)
    comps = []
    rr, cc = np.nonzero(mask)
    for r0, c0 in zip(rr, cc):
        if lab[r0, c0]:
            continue
        n = len(comps) + 1
        stack, pts = [(r0, c0)], []
        lab[r0, c0] = n
        while stack:
            r, c = stack.pop()
            pts.append((r, c))
            for dr, dc in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                a, b = r + dr, c + dc
                if 0 <= a < H and 0 <= b < W and mask[a, b] and not lab[a, b]:
                    lab[a, b] = n
                    stack.append((a, b))
        p = np.array(pts)
        comps.append((p[:, 0], p[:, 1]))
    return comps


def erode(mask, rad):
    """Binary erosion by a disk of radius rad (outside the image counts as clear)."""
    out = mask.copy()
    H, W = mask.shape
    for dr in range(-rad, rad + 1):
        for dc in range(-rad, rad + 1):
            if dr * dr + dc * dc > rad * rad:
                continue
            sh = np.zeros_like(mask)
            sh[max(0, -dr):H - max(0, dr), max(0, -dc):W - max(0, dc)] = mask[max(0, dr):H - max(0, -dr), max(0, dc):W - max(0, -dc)]
            out &= sh
    return out


def dots(img):
    dark = img.max(axis=2) < 80
    found = []
    for rows, cols in label(erode(dark, 5)):
        n = len(rows)
        h, w = rows.max() - rows.min() + 1, cols.max() - cols.min() + 1
        # The cube, the arrow heads (triangles) and text remnants are not full disks of a dot's size.
        if n < 4 or n > 400 or h > 22 or w > 22 or n / (h * w) < 0.65:
            continue
        found.append((rows.mean(), cols.mean(), n))
    return found


def profiles(img):
    """{(panel, station, coordinate): u}; panel a's coordinate is y, panel b's z."""
    pts = []
    full = {}
    ds = dots(img)
    for panel in 'ab':
        full[panel] = max(n for r, c, n in ds if (r < 1450) == (panel == 'a'))
    for r, c, n in ds:
        panel = 'a' if r < 1450 else 'b'
        top, bot = FRAME[panel]
        clipped = n < 0.5 * full[panel] and min(abs(r - top), abs(r - bot)) < 10
        rr = (top if abs(r - top) < abs(r - bot) else bot) if clipped else r
        coord = (0.8 - (rr - 288.0) / 700.0) / 2 if panel == 'a' else 0.4 - (rr - 1576.0) / PX
        pts.append([panel, None, round(coord, 3), xdata(c)])
    # Level by level: an order-preserving match of the dots (sorted by x) to the stations, each offset in the
    # window, each station one dot at most; the smallest sum of |offset| among the matches that place every dot.
    lo, hi = -0.035, 0.105
    levels = {}
    for p in pts:
        levels.setdefault((p[0], p[2]), []).append(p)
    for (panel, _), ps in levels.items():
        ps.sort(key=lambda p: p[3])
        st = STATIONS[panel]
        inf = float('inf')
        best = [[inf] * (len(st) + 1) for _ in range(len(ps) + 1)]
        take = [[False] * (len(st) + 1) for _ in range(len(ps) + 1)]
        for j in range(len(st) + 1):
            best[0][j] = 0.0
        for i in range(1, len(ps) + 1):
            for j in range(1, len(st) + 1):
                off = ps[i - 1][3] - st[j - 1]
                t = best[i - 1][j - 1] + abs(off) if lo <= off <= hi else inf
                if t <= best[i][j - 1]:
                    best[i][j], take[i][j] = t, True
                else:
                    best[i][j] = best[i][j - 1]
        if best[len(ps)][len(st)] == inf:
            raise SystemExit(f'no station match at level {panel} {ps[0][2]}')
        i, j = len(ps), len(st)
        while i > 0:
            if take[i][j]:
                ps[i - 1][1] = st[j - 1]
                i -= 1
            j -= 1
    return {(p[0], p[1], p[2]): (p[3] - p[1]) / SCALE for p in pts}


def edge_dot(img, row, col_from, radius):
    """x of a dot whose left part is covered (the arrow's head): its right edge on its centre row, less its radius."""
    lum = img.mean(axis=2)
    c = col_from
    while lum[row, c] >= 128:
        c -= 1
    right = c + (128 - lum[row, c]) / (lum[row, c + 1] - lum[row, c])
    return xdata(right - radius)


def dot_radius(img, found):
    """The median half-width of the dots of panel b on their centre rows (sub-pixel at the 128 crossing)."""
    lum = img.mean(axis=2)
    halves = []
    for r, c, n in found:
        if r < 1500 or n < 0.9 * max(m for _, _, m in found):
            continue
        r, a = int(round(r)), int(round(c))
        b = a
        while lum[r, a - 1] < 128:
            a -= 1
        while lum[r, b + 1] < 128:
            b += 1
        left = a - (128 - lum[r, a]) / (lum[r, a - 1] - lum[r, a])
        right = b + (128 - lum[r, b]) / (lum[r, b + 1] - lum[r, b])
        halves.append((right - left) / 2)
    return float(np.median(halves))


def main():
    a3 = np.asarray(Image.open(sys.argv[1]).convert('RGB')).astype(float)
    a4 = np.asarray(Image.open(sys.argv[2]).convert('RGB')).astype(float)
    p3, p4 = profiles(a3), profiles(a4)
    common = sorted(set(p3) & set(p4))
    clear = [k for k in common if k[1] != 0.5]
    rms = (sum((p3[k] - p4[k]) ** 2 for k in clear) / len(clear)) ** 0.5
    print(f'A3 {len(p3)} dots, A4 {len(p4)}, shared {len(common)}; rms difference away from x = 0.5 m: {rms:.4f} m/s')
    data = {k: (u, 'A3') for k, u in p3.items()}
    for k in sorted(set(p4) - set(p3)):
        data[k] = (p4[k], 'A4')
    # The dot under the 3 m/s arrow's head (x = -0.3 m, z = 0.2 m: row 1856).
    x = edge_dot(a3, 1856, 600, dot_radius(a3, dots(a3)))
    data[('b', -0.3, 0.2)] = ((x + 0.3) / SCALE, 'A3 edge')
    with open(sys.argv[3], 'w', newline='') as fh:
        fh.write('# Snyder and Lawson (1994) cube, H = 0.2 m: mean streamwise velocity digitized from Wise et al. (2026),\n'
                 '# egusphere-2026-4030, Figs. A3 and A4 (CC BY 4.0), by digitize_wise2026.py; see README.md.\n'
                 '# x from the cube upwind face, y from the centreline, z above the floor (m); u (m/s); source: the figure.\n')
        w = csv.writer(fh, lineterminator='\n')
        w.writerow(['plane', 'x', 'y', 'z', 'u', 'source'])
        for (panel, xs, coord), (u, src) in sorted(data.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[0][2])):
            y, z = (coord, 0.02) if panel == 'a' else (0.0, coord)
            w.writerow(['horizontal' if panel == 'a' else 'vertical', f'{xs:g}', f'{y:g}', f'{z:g}', f'{u:.3f}', src])
    print(f'wrote {len(data)} values to {sys.argv[3]}')


if __name__ == '__main__':
    main()
