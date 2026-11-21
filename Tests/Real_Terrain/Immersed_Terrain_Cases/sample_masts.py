#!/usr/bin/env python3
"""Masts and lines of a case from a body-fitted or an immersed plotfile (README.md; windninja2019.py sample).

    sample_masts.py PLT --site S --case C --masts MASTS.csv (--terrain XYZ | --erf_terrain FILE) --out CSV
                    [--z0 Z0] [--flat H --zlevels FILE]          (--flat: the plotfile is of an immersed run)

The ground under a column is the mean of its four corner nodes' heights for both terrain kinds. An immersed
column gives its fluid cells only, the cut cell at the centre of its open part."""
import sys, os, argparse
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, REPO + "/Tests/Real_Terrain/WindNinja2019"); sys.path.insert(0, REPO + "/tools"); sys.path.insert(0, REPO + "/extern_erf")
sys.dont_write_bytecode = True
import windninja2019 as wn  # noqa: E402
al = wn.al

class ImmCols:
    """Columns of an immersed run: fluid cells only, the cut cell at its fluid centroid."""
    def __init__(self, base, ground, lev):
        self.b, self.g, self.lev, self.hdr = base, ground, lev, base.hdr
    def column(self, i, j):
        u, v, z = self.b.column(i, j); h = self.hdr
        zg = self.g(h["prob_lo"][0] + (i + 0.5) * h["dx"][0], h["prob_lo"][1] + (j + 0.5) * h["dx"][1])
        uo, vo, zo = [], [], []
        for k in range(len(u)):
            lo, hi = self.lev[k], self.lev[k + 1]
            if hi <= zg + 1e-9 * (hi - lo): continue
            uo.append(u[k]); vo.append(v[k]); zo.append(0.5 * (hi + max(lo, zg)))
        return uo, vo, zo

def main():
    p = argparse.ArgumentParser(); p.add_argument("plotfile")
    for k in ("site", "case", "out", "masts"): p.add_argument("--" + k, required=True)
    p.add_argument("--terrain"); p.add_argument("--erf_terrain")
    p.add_argument("--flat", type=float); p.add_argument("--zlevels"); p.add_argument("--z0", type=float, default=0.01)
    a = p.parse_args()
    flat = a.flat or 0.0
    Cols, Idw, erf_ground = al.Columns, wn.Idw, al.ground
    H = {}
    class ShiftIdw(Idw):                       # xyz terrain: windninja2019.sample takes the corner mean itself
        def __call__(self, x, y): return Idw.__call__(self, x, y) + flat
    def erf_shift(grid, x, y):                 # ERF grid: the corner mean here (sample would take the centre)
        dx = H["dx"]
        return flat + 0.25 * sum(erf_ground(grid, x + sx * dx[0] / 2, y + sy * dx[1] / 2) for sx in (-1, 1) for sy in (-1, 1))
    wn.Idw = ShiftIdw; al.ground = erf_shift
    def cols(path):
        b = Cols(path); H["dx"] = b.hdr["dx"]
        if a.flat is None: return b
        lev = [float(x) for x in open(a.zlevels).read().split()]
        if a.terrain:
            idw = ShiftIdw(a.terrain); dx = b.hdr["dx"]
            g = lambda x, y: 0.25 * sum(idw(x + sx * dx[0] / 2, y + sy * dx[1] / 2) for sx in (-1, 1) for sy in (-1, 1))
        else:
            grid = al.read_terrain(a.erf_terrain); g = lambda x, y: erf_shift(grid, x, y)
        return ImmCols(b, g, lev)
    al.Columns = cols
    wn.sample(a)
main()
