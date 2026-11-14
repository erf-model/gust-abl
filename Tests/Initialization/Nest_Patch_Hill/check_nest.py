#!/usr/bin/env python3
"""Run the nested patch twice and check it against its parent (rt_nest_patch_hill).

    check_nest.py PARENT_PLOTFILE OFFSET_X OFFSET_Y RELTOL -- COMMAND [ARGS ...]

COMMAND is the patch run (init.first_guess = nest). It runs twice:

1. with init.projection = 0 appended (the parent interpolated at 2 2 1, as it
   is); its plotfile is renamed nest_noproj_plt00000. theta and k of every patch cell must equal
   the bilinear interpolation, per zeta-level, of the parent's cell values
   (PARENT_PLOTFILE) at the patch cell's centre, to a relative 1e-12 (the
   CTest is double precision only). The parent's cell (I, J) holds patch cells
   2 (I - OFFSET) and 2 (I - OFFSET) + 1, so the patch centre i sits at
   parent index (i + 1/2) / 2 - 1/2 + OFFSET. This is computed here, from the
   plotfiles, independently of the code's Transfer_K.H. The parent's theta
   varies along a level over the hill (a stable sounding), and
   its k with the height above the ground, so neither is a constant;
2. as given: the log names the parent and the parent cells covered, the
   divergence after the correction is at most 10 x RELTOL x the divergence
   before, and no mass crosses the ground: on every ground face of the
   patch's checkpoint (nest_patch_chk00000, ghost rows included)
   Omega = w - z_x u - z_y v, with the code's averages (u over its two x-faces
   and the cell below and above the face, the slopes from the face's four
   nodes), is at most 1e-12 x max |w| there. The correction keeps the ground
   flux fixed, so a first guess with flux there keeps it.

The exit code is the verdict.
"""
import math
import os
import re
import struct
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "extern_erf"))
sys.dont_write_bytecode = True
import erf_plotfile  # noqa: E402


def run(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=1800)
    print(proc.stdout)
    if proc.returncode != 0:
        raise SystemExit(f"check_nest: the run exited {proc.returncode}")
    return proc.stdout


def bilinear(a, X, Y, k):
    I, J = math.floor(X), math.floor(Y)
    tx, ty = X - I, Y - J
    return ((1 - tx) * (1 - ty) * a[I][J][k] + tx * (1 - ty) * a[I + 1][J][k]
            + (1 - tx) * ty * a[I][J + 1][k] + tx * ty * a[I + 1][J + 1][k])


def read_fabs(prefix):
    """A checkpoint MultiFab (double precision) as {(i, j, k): value}, ghosts
    included; the valid boxes are written last, so they win."""
    head = open(prefix + "_H").read()
    boxes = re.findall(r"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)", head)
    fabs = re.findall(r"FabOnDisk: (\S+) (\d+)", head)
    nvalid = len(fabs)
    out, valid = {}, []
    for n, (name, off) in enumerate(fabs):
        with open(os.path.join(os.path.dirname(prefix), name), "rb") as f:
            f.seek(int(off))
            hdr = f.readline().decode()
            m = re.search(r"\(\((-?\d+),(-?\d+),(-?\d+)\) \((-?\d+),(-?\d+),(-?\d+)\)", hdr)
            lo = [int(x) for x in m.groups()[:3]]
            hi = [int(x) for x in m.groups()[3:]]
            nx, ny, nz = [hi[d] - lo[d] + 1 for d in range(3)]
            data = struct.unpack("<%dd" % (nx * ny * nz), f.read(8 * nx * ny * nz))
        vlo = [int(x) for x in boxes[n][:3]] if n < nvalid else lo
        vhi = [int(x) for x in boxes[n][3:]] if n < nvalid else hi
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    p = (lo[0] + i, lo[1] + j, lo[2] + k)
                    v = data[i + nx * (j + ny * k)]
                    if all(vlo[d] <= p[d] <= vhi[d] for d in range(3)):
                        valid.append((p, v))
                    else:
                        out.setdefault(p, v)
    out.update(valid)
    return out


def ground_flux(chk):
    """max |Omega| and max |w| on the ground faces of checkpoint chk."""
    head = open(os.path.join(chk, "Header")).read()
    n_cell = [int(x) for x in re.search(r"^n_cell (\S+) (\S+) (\S+)", head, re.M).groups()]
    plo = [float.fromhex(x) for x in re.search(r"^prob_lo (\S+) (\S+) (\S+)", head, re.M).groups()]
    phi = [float.fromhex(x) for x in re.search(r"^prob_hi (\S+) (\S+) (\S+)", head, re.M).groups()]
    dx = (phi[0] - plo[0]) / n_cell[0]
    dy = (phi[1] - plo[1]) / n_cell[1]
    u, v, w, z = (read_fabs(os.path.join(chk, "Level_0", name)) for name in ("u", "v", "w", "z_phys_nd"))
    worst, wmax = 0.0, 0.0
    for i in range(n_cell[0]):
        for j in range(n_cell[1]):
            zx = 0.5 * ((z[i + 1, j, 0] - z[i, j, 0]) + (z[i + 1, j + 1, 0] - z[i, j + 1, 0])) / dx
            zy = 0.5 * ((z[i, j + 1, 0] - z[i, j, 0]) + (z[i + 1, j + 1, 0] - z[i + 1, j, 0])) / dy
            uf = 0.25 * sum(u[ii, j, kk] for ii in (i, i + 1) for kk in (-1, 0))
            vf = 0.25 * sum(v[i, jj, kk] for jj in (j, j + 1) for kk in (-1, 0))
            worst = max(worst, abs(w[i, j, 0] - zx * uf - zy * vf))
            wmax = max(wmax, abs(w[i, j, 0]))
    return worst, wmax


def main(argv):
    if "--" not in argv or argv.index("--") != 4:
        print(__doc__)
        return 2
    parent, ox, oy, reltol = argv[0], int(argv[1]), int(argv[2]), float(argv[3])
    cmd = argv[5:]
    fails = []

    # 1. The interpolated parent, without the correction.
    out = run(cmd + ["init.projection=0", "io.write_checkpoint=0"])
    w = re.search(r"gust: wrote (\S+)", out)
    if not w:
        raise SystemExit("check_nest: the first run wrote no plotfile")
    shutil.rmtree("nest_noproj_plt00000", ignore_errors=True)
    os.rename(w.group(1), "nest_noproj_plt00000")
    _, pf = erf_plotfile.read_fields(parent, ["theta", "tke"])
    hdr, cf = erf_plotfile.read_fields("nest_noproj_plt00000", ["theta", "tke"])
    nx, ny, nz = [hdr["hi"][d] - hdr["lo"][d] + 1 for d in range(3)]
    tol = 1.0e-12
    for name in ("theta", "tke"):
        worst, spread = 0.0, []
        for i in range(nx):
            X = 0.5 * (i + 0.5) - 0.5 + ox
            for j in range(ny):
                Y = 0.5 * (j + 0.5) - 0.5 + oy
                for k in range(nz):
                    want = bilinear(pf[name], X, Y, k)
                    got = cf[name][i][j][k]
                    worst = max(worst, abs(got - want) / max(abs(want), 1.0e-300))
                    if k == 2:
                        spread.append(got)
        rng = max(spread) - min(spread)
        print(f"check_nest: {name}: max relative difference {worst:.3e} (tolerance {tol:.0e}); "
              f"range along level 2 {rng:.4e}")
        if not worst <= tol:
            fails.append(f"{name} differs from the parent's bilinear interpolation by {worst:.3e} (relative)")
        if not rng > 1.0e-6:
            fails.append(f"{name} is constant along level 2: the check would not see a wrong position")

    # 2. With the correction.
    out = run(cmd)
    if not re.search(rf"gust: nest: parent \S+ \(step \d+, .*cover parent cells {ox} to \d+ in x and {oy} to \d+ in y, "
                     r"refinement 2 2 1", out):
        fails.append("the log has no nest line naming the parent cells covered")
    m = re.search(r"gust: initializer: (\d+) Poisson iterations \((?:\w+, )?reltol (\S+)\), divergence (\S+) before, "
                  r"(\S+) after", out)
    if not m:
        fails.append("the log has no initializer line")
    else:
        before, after = float(m.group(3)), float(m.group(4))
        print(f"check_nest: divergence {before:.3e} before, {after:.3e} after (limit {10 * reltol * before:.3e})")
        if not (before > 0 and after <= 10 * reltol * before):
            fails.append(f"divergence after the correction {after:.3e} above 10 x {reltol:g} x {before:.3e}")

    if not os.path.isdir("nest_patch_chk00000"):
        fails.append("the second run wrote no nest_patch_chk00000")
    else:
        om, wg = ground_flux("nest_patch_chk00000")
        print(f"check_nest: ground faces: max |Omega| {om:.3e}, max |w| {wg:.3e} m/s")
        if not (wg > 0 and om <= 1.0e-12 * wg):
            fails.append(f"flux through the ground: max |Omega| {om:.3e} above 1e-12 x max |w| {wg:.3e}")

    for f in fails:
        print("check_nest: FAIL: " + f)
    print("check_nest: " + ("FAILED" if fails else "passed"))
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
