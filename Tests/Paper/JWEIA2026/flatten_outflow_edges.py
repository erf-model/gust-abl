#!/usr/bin/env python3
"""Give an xyz terrain no slope normal to the outflow faces (x high and y high) over the last WIDTH metres: the
terrain is blended with a smooth step over the first half of that width towards its value at the start of the
strip, and is constant in the normal direction over the second half. Used for the Big Southern Butte runs with
level outflow edges (Appendix of the paper; Tests/Canonical_RANS/Outflow_Steep_Terrain).

    python3 flatten_outflow_edges.py IN.xyz OUT.xyz [--width 680]
"""
import argparse


def step(t):
    return 0.0 if t <= 0 else 1.0 if t >= 1 else t * t * (3 - 2 * t)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("--width", type=float, default=680.0)
    a = ap.parse_args()
    hdr, pts = [], {}
    for line in open(a.src):
        if line.startswith("#"):
            hdr.append(line)
            continue
        x, y, z = map(float, line.split()[:3])
        pts[(x, y)] = z
    xs = sorted({p[0] for p in pts})
    ys = sorted({p[1] for p in pts})
    xe = min(xs, key=lambda v: abs(v - (xs[-1] - a.width)))
    ye = min(ys, key=lambda v: abs(v - (ys[-1] - a.width)))
    out = {}
    for (x, y), z in pts.items():
        if x > xe:
            w = step((x - xe) / (xs[-1] - xe) * 2)
            z = (1 - w) * z + w * pts[(xe, y)]
        out[(x, y)] = z
    first = dict(out)
    for (x, y), z in first.items():
        if y > ye:
            w = step((y - ye) / (ys[-1] - ye) * 2)
            out[(x, y)] = (1 - w) * z + w * first[(x, ye)]
    with open(a.dst, "w") as f:
        f.writelines(hdr)
        f.write(f"# outflow edges x > {xe:g} and y > {ye:g}: no normal slope (smooth step over the first half)\n")
        for y in sorted(ys, reverse=True):
            for x in xs:
                f.write(f"{x} {y} {out[(x, y)]:.2f}\n")
    print(a.dst, len(out), "points; xe", xe, "ye", ye)


if __name__ == "__main__":
    main()
