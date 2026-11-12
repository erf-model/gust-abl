#!/usr/bin/env python3
"""Digitize the measurements of Castro et al. (2006) in the cube array from Fig. 3 of Munoz-Esparza et al. (2020).

    digitize_munoz2020.py PAGE7.svg OUT.csv

PAGE7.svg is page 7 of the paper (https://doi.org/10.1029/2020MS002141) converted with
`pdftocairo -svg -f 7 -l 7 paper.pdf PAGE7.svg`. The figure is vector graphics, so every measurement symbol is a
filled path: black (P0 and P1), grey (P3) and blue (P2). README.md has the sources.

The figure's scales, read from the axis frames in the SVG's coordinates (points):
- z / h = 3 (286.045 - y) / (286.045 - 103.865) in all three panels (z / h from 0 to 3).
- Panel (a), x from 94.02 to 215.47 for u / u_r from -0.2 to 0.8; (b), 258.41 to 379.86, the same range;
  (c), 422.75 to 544.21 for -<u'w'> / u_r^2 from -2 to 10 (x 1e-3).

Symbols inside each panel's legend box are left out, and so are paths with 100 or more points (lines and text).
A symbol partly covered by its neighbour is smaller than the full ones: its centre is moved half a full symbol's
width (or height) from its free edge, the edge away from the nearest symbol. A symbol cut by a legend box's left
edge (one, at P3 in panel c near z = 2.16 h) is placed from its left edge. Other fragments under 60 % of a full
symbol's size are dropped (the script lists them). The output's -<u'w'> values are in units of 1e-3 u_r^2, as
plotted. Not a test.
"""
import collections
import re
import sys

Z0, Z3 = 286.045, 103.865
PANELS = {"a": (94.02, 215.47, -0.2, 0.8), "b": (258.41, 379.86, -0.2, 0.8), "c": (422.75, 544.21, -2.0, 10.0)}
LEGEND = {"a": (94, 160, 103, 160), "b": (263.72, 321.57, 108.46, 158.73), "c": (484.84, 538.89, 109.0, 159.27)}
COLOURS = {(1, 1, 1): "black", (60, 60, 60): "grey", (1, 1, 100): "blue"}
SERIES = {("a", "black"): ("P0", "u"), ("a", "grey"): ("P3", "u"), ("b", "black"): ("P1", "u"),
          ("b", "blue"): ("P2", "u"), ("c", "black"): ("P1", "uw"), ("c", "grey"): ("P3", "uw")}


def load(path):
    """Bounding boxes of the SVG's filled paths, with their fill colour (percent, rounded)."""
    out = []
    for m in re.finditer(r"<path ([^>]*)/>", open(path).read()):
        a = m.group(1)
        f = re.search(r'fill="rgb\(([^)]*)\)"', a)
        if not f:
            continue
        col = tuple(round(float(v.strip(" %"))) for v in f.group(1).split(","))
        nums = [float(v) for v in re.findall(r"-?\d+\.?\d*", re.search(r'd="([^"]*)"', a).group(1))]
        xs, ys = nums[0::2], nums[1::2]
        if xs:
            out.append(dict(col=col, x0=min(xs), x1=max(xs), y0=min(ys), y1=max(ys), n=len(xs)))
    return out


def main():
    svg, out_csv = sys.argv[1], sys.argv[2]
    groups = collections.defaultdict(list)
    for p in load(svg):
        c = COLOURS.get(p["col"])
        if not c or p["n"] >= 100:
            continue
        xc, yc = (p["x0"] + p["x1"]) / 2, (p["y0"] + p["y1"]) / 2
        for k, (lft, rgt, _, _) in PANELS.items():
            if not (lft - 3 < xc < rgt + 3 and Z3 - 3 < yc < Z0 + 3):
                continue
            lg = LEGEND[k]
            if lg[0] < xc < lg[1] and lg[2] < yc < lg[3]:
                continue
            groups[(k, c)].append(p)
    rows = []
    for key, ps in sorted(groups.items()):
        if key not in SERIES:
            continue
        wf = max(p["x1"] - p["x0"] for p in ps)
        hf = max(p["y1"] - p["y0"] for p in ps)
        cs = [((p["x0"] + p["x1"]) / 2, (p["y0"] + p["y1"]) / 2) for p in ps]
        lft, rgt, v0, v1 = PANELS[key[0]]
        for i, p in enumerate(ps):
            w, h = p["x1"] - p["x0"], p["y1"] - p["y0"]
            xc, yc = cs[i]
            lg = LEGEND[key[0]]
            if abs(p["x1"] - lg[0]) < 0.05 and h >= 0.6 * hf:
                # Cut by the legend box's left edge: the symbol's left edge is free.
                rows.append((*SERIES[key], round((Z0 - yc) / (Z0 - Z3) * 3, 4),
                             round(v0 + (p["x0"] + wf / 2 - lft) / (rgt - lft) * (v1 - v0), 5)))
                continue
            if w < 0.6 * wf or h < 0.6 * hf:
                print("dropped fragment", key, round(xc, 1), round(yc, 1), file=sys.stderr)
                continue
            j = min((m for m in range(len(cs)) if m != i), key=lambda m: (cs[m][0] - xc) ** 2 + (cs[m][1] - yc) ** 2)
            nx, ny = cs[j]
            if wf - w > 0.05:
                xc = p["x1"] - wf / 2 if nx < xc else p["x0"] + wf / 2
            if hf - h > 0.05:
                yc = p["y1"] - hf / 2 if ny < yc else p["y0"] + hf / 2
            val = v0 + (xc - lft) / (rgt - lft) * (v1 - v0)
            rows.append((*SERIES[key], round((Z0 - yc) / (Z0 - Z3) * 3, 4), round(val, 5)))
    rows.sort()
    with open(out_csv, "w") as fh:
        fh.write("# Castro et al. (2006) staggered cube array (h = 20 mm), as plotted in Fig. 3 of Munoz-Esparza et al.\n"
                 "# (2020, https://doi.org/10.1029/2020MS002141, CC BY-NC-ND 4.0), digitized by digitize_munoz2020.py; the\n"
                 "# measured values only, no part of the figure; see README.md (Data and acknowledgements).\n"
                 "# u is u / u_r (u_r = 10 m/s); uw is -<u'w'> / u_r^2 in units of 1e-3.\n")
        fh.write("station,quantity,z_over_h,value\n")
        for r in rows:
            fh.write(",".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
