#!/usr/bin/env python3
"""A single-band GeoTIFF (integer or float samples of 8 to 64 bits, by the
SampleFormat tag; uncompressed or LZW, strips, no predictor) as a GUST
terrain.type = xyz file, without GDAL: one "x y z" line per pixel centre
(every STRIDE-th pixel), in the TIFF's own frame (ModelTiepoint +
PixelScale; PixelIsPoint tie points shifted to the pixel corner),
z shifted by -ZSHIFT. Pixels equal to the GDAL_NODATA tag's
value, or at or below -1e30, are no data and skipped.

Usage: tiff_to_xyz.py IN.tif OUT.xyz [STRIDE] [ZSHIFT]
"""
import struct
import sys

import numpy as np


def lzw_decode(data):
    out = bytearray()
    table = [bytes([i]) for i in range(256)] + [b"", b""]
    width, bitpos, prev = 9, 0, None
    nbits = len(data) * 8
    while bitpos + width <= nbits:
        byte, bit = divmod(bitpos, 8)
        chunk = int.from_bytes(data[byte:byte + 3].ljust(3, b"\0"), "big")
        code = (chunk >> (24 - bit - width)) & ((1 << width) - 1)
        bitpos += width
        if code == 256:  # clear
            table = table[:258]
            width, prev = 9, None
            continue
        if code == 257:  # end of information
            break
        if prev is None:
            entry = table[code]
        elif code < len(table):
            entry = table[code]
            table.append(prev + entry[:1])
        else:
            entry = prev + prev[:1]
            table.append(entry)
        out += entry
        prev = entry
        if len(table) + 1 >= (1 << width) and width < 12:  # early change
            width += 1
    return bytes(out)


def read_tiff(path):
    with open(path, "rb") as handle:
        f = handle.read()
    assert f[:2] == b"II", "little-endian TIFF only"
    off = struct.unpack("<I", f[4:8])[0]
    n = struct.unpack("<H", f[off:off + 2])[0]
    size = {1: 1, 2: 1, 3: 2, 4: 4, 8: 2, 9: 4, 11: 4, 12: 8}
    fmt = {1: "B", 2: "s", 3: "H", 4: "I", 8: "h", 9: "i", 11: "f", 12: "d"}
    tags = {}
    for i in range(n):
        tag, typ, cnt, val = struct.unpack("<HHII", f[off + 2 + 12 * i:off + 14 + 12 * i])
        if typ not in fmt:
            continue
        nb = size[typ] * cnt
        raw = f[val:val + nb] if nb > 4 else f[off + 10 + 12 * i:off + 10 + 12 * i + nb]
        if typ == 2:      # ASCII, e.g. GDAL_NODATA
            tags[tag] = raw.split(b"\0")[0].decode("ascii", "replace").strip()
        else:
            tags[tag] = struct.unpack("<" + fmt[typ] * cnt, raw)
    nx, ny, bps = tags[256][0], tags[257][0], tags[258][0]
    comp, pred = tags[259][0], tags.get(317, (1,))[0]
    sample = tags.get(339, (1,))[0]      # SampleFormat: 1 unsigned, 2 signed, 3 float (TIFF 6.0 default 1)
    dtype = {(1, 8): "u1", (1, 16): "<u2", (1, 32): "<u4", (2, 8): "i1", (2, 16): "<i2", (2, 32): "<i4",
             (3, 32): "<f4", (3, 64): "<f8"}.get((sample, bps))
    assert comp in (1, 5) and pred == 1 and dtype is not None, (comp, pred, sample, bps)
    rows = []
    for o, c in zip(tags[273], tags[279]):
        chunk = f[o:o + c]
        rows.append(lzw_decode(chunk) if comp == 5 else chunk)
    a = np.frombuffer(b"".join(rows), dtype=dtype)[:nx * ny].reshape(ny, nx).astype(float)
    if 42113 in tags:     # GDAL_NODATA
        a[a == float(tags[42113])] = np.nan
    sx, sy = tags[33550][0], tags[33550][1]
    tie = tags[33922]
    x0, y0 = tie[3] - tie[0] * sx, tie[4] + tie[1] * sy  # the corner of pixel (0, 0)
    keys = tags.get(34735, ())          # GeoKeyDirectory: GTRasterTypeGeoKey (1025) = 2 is PixelIsPoint, where
    for k in range(4, len(keys) - 3, 4):  # the tie point is the pixel's centre (shifted to its corner, as GDAL does)
        if keys[k] == 1025 and keys[k + 1] == 0 and keys[k + 3] == 2:
            x0, y0 = x0 - 0.5 * sx, y0 + 0.5 * sy
    return a, x0, y0, sx, sy


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    stride = int(argv[2]) if len(argv) > 2 else 1
    zshift = float(argv[3]) if len(argv) > 3 else 0.0
    a, x0, y0, sx, sy = read_tiff(argv[0])
    ny, nx = a.shape
    n = 0
    with open(argv[1], "w") as out:
        out.write(f"# {argv[0]}: {nx} x {ny} pixels of {sx} x {sy} m, stride {stride}, z - {zshift}\n")
        for j in range(0, ny, stride):
            y = y0 - (j + 0.5) * sy
            for i in range(0, nx, stride):
                z = a[j, i]
                if z > -1e30:     # False for NaN (no data) too
                    out.write(f"{x0 + (i + 0.5) * sx:.3f} {y:.3f} {z - zshift:.3f}\n")
                    n += 1
    good = a[a > -1e30]
    print(f"{argv[1]}: {n} points; x {x0 + 0.5 * sx:.2f} .. {x0 + (nx - 0.5) * sx:.2f}, "
          f"y {y0 - (ny - 0.5) * sy:.2f} .. {y0 - 0.5 * sy:.2f}, z {good.min():.2f} .. {good.max():.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
