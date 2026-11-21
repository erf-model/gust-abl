"""GUST plotfiles as numpy arrays: cell fields and the mesh's node heights.

    read_cells(plotfile, fields)   {name: array (nx, ny, nz)} of cell-centred fields, and the header
    read_nodes(plotfile)           the node heights z_phys_nd, array (nx + 1, ny + 1, nz + 1)

GUST writes the cell fields u, v, w, p, theta, tke and z (the cell centre's
height) with AMReX's single-level plotfile writer, and the mesh's node heights
beside them as Level_0/z_phys_nd (Source/StateIO.cpp, write_plotfile). The
FAB data are read with extern_erf/erf_plotfile.py's reader, which takes the
size and byte order from each FAB's own header (float or double).
"""
import os
import re
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "extern_erf"))
import erf_plotfile  # noqa: E402


def _multifab(prefix, comps, ncomp_expected=None):
    """The valid boxes of a VisMF MultiFab written at prefix (prefix + '_H' is its header), components comps only.

    Returns (ncomp, [(lo, hi, {comp: array (nz, ny, nx) of the valid box})]). A header that is not VisMF
    version 1 with one FAB per box (as GUST writes it), or is cut short, is refused naming the file.
    """
    header = prefix + "_H"
    if not os.path.isfile(header):
        raise ValueError("%s: no such file (a GUST plotfile holds it)" % header)
    try:
        return _read_multifab(prefix, header, comps, ncomp_expected)
    except (IndexError, AttributeError, KeyError, struct.error) as err:
        raise ValueError("%s: not a complete VisMF MultiFab (%s: %s)" % (header, type(err).__name__, err))


def _read_multifab(prefix, header, comps, ncomp_expected):
    with open(header) as fh:
        lines = [ln.rstrip("\n") for ln in fh]
    if lines[0].strip() != "1":
        raise ValueError("%s: VisMF header version %s; version 1 expected" % (header, lines[0].strip()))
    ncomp = int(lines[2])
    if ncomp_expected is not None and ncomp != ncomp_expected:
        raise ValueError("%s: %d components, %d expected" % (header, ncomp, ncomp_expected))
    nbox = int(re.match(r"\((\d+)\s+(\d+)", lines[4]).group(1))
    boxes = []
    for n in range(nbox):
        lo, hi = re.findall(r"\(([-\d,]+)\)", lines[5 + n])[:2]
        boxes.append(([int(v) for v in lo.split(",")], [int(v) for v in hi.split(",")]))
    idx = 5 + nbox + 1
    nfab = int(lines[idx])
    fabs = []
    for n in range(nfab):
        m = re.match(r"FabOnDisk:\s+(\S+)\s+(\d+)", lines[idx + 1 + n])
        fabs.append((m.group(1), int(m.group(2))))
    if nfab != nbox:
        raise ValueError("%s: %d boxes and %d FABs" % (header, nbox, nfab))
    out = []
    folder = os.path.dirname(prefix)
    for (blo, bhi), (fname, off) in zip(boxes, fabs):
        flo, fhi, data = erf_plotfile._read_fab(os.path.join(folder, fname), off, ncomp)
        shape = [fhi[d] - flo[d] + 1 for d in range(3)]
        a = np.frombuffer(data, dtype=np.float64 if data.typecode == "d" else np.float32)
        a = a.reshape(ncomp, shape[2], shape[1], shape[0])
        sl = tuple(slice(blo[d] - flo[d], bhi[d] - flo[d] + 1) for d in (2, 1, 0))
        out.append((blo, bhi, {c: a[(c,) + sl].astype(np.float64) for c in comps}))
        del a, data
    return ncomp, out


def _assemble(boxes, shape, comps):
    """{comp: array (nx, ny, nz)} from the valid boxes; every point must be covered."""
    arrs = {c: np.full(shape, np.nan) for c in comps}
    seen = np.zeros(shape, dtype=bool)
    for blo, bhi, data in boxes:
        sl = tuple(slice(blo[d], bhi[d] + 1) for d in range(3))
        for c in comps:
            arrs[c][sl] = data[c].transpose(2, 1, 0)
        seen[sl] = True
    if not seen.all():
        raise ValueError("the boxes do not cover the domain (%d of %d points missing)"
                         % (int((~seen).sum()), seen.size))
    return arrs


def _header(plotfile):
    if not os.path.isfile(os.path.join(plotfile, "Header")):
        raise ValueError("%s: no Header (not a plotfile)" % plotfile)
    try:
        return erf_plotfile._read_header(plotfile)
    except (StopIteration, IndexError, ValueError) as err:
        raise ValueError("%s/Header: not a complete plotfile header (%s)" % (plotfile, type(err).__name__))


def read_cells(plotfile, fields):
    """Cell-centred fields of a single-level plotfile: ({name: array (nx, ny, nz)}, header)."""
    hdr = _header(plotfile)
    if hdr["finest"] != 0:
        raise ValueError("%s: %d refined levels; one level expected" % (plotfile, hdr["finest"]))
    if any(v != 0 for v in hdr["lo"]):
        raise ValueError("%s: the domain does not start at index 0" % plotfile)
    names = hdr["names"]
    missing = [f for f in fields if f not in names]
    if missing:
        raise ValueError("%s: no field %s (has %s)" % (plotfile, ", ".join(missing), ", ".join(names)))
    shape = tuple(h + 1 for h in hdr["hi"])
    comps = [names.index(f) for f in fields]
    _, boxes = _multifab(os.path.join(plotfile, "Level_0", "Cell"), comps, len(names))
    arrs = _assemble(boxes, shape, comps)
    return {f: arrs[names.index(f)] for f in fields}, hdr


def read_nodes(plotfile):
    """The mesh's node heights Level_0/z_phys_nd: array (nx + 1, ny + 1, nz + 1)."""
    hdr = _header(plotfile)
    shape = tuple(h + 2 for h in hdr["hi"])
    _, boxes = _multifab(os.path.join(plotfile, "Level_0", "z_phys_nd"), [0], 1)
    return _assemble(boxes, shape, [0])[0]
