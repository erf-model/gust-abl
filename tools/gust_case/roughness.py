"""Roughness helpers: a land-cover raster or a point file as GUST's
``bc.surface_z0_file`` (x y z0 points).

The ESA CCI table maps the ``lccs_class`` codes of the ESA Climate Change
Initiative Land Cover maps (ESA, Land Cover CCI Product User Guide, version
2.0, 2017) to roughness lengths. It starts from amrTerrain's table
(src/backend/backendInterface.py) with five changes, each to the class of the
updated Davenport classification (Wieringa, J. Wind Eng. Ind. Aerodyn. 41,
357-368, 1992) that the land cover describes:

- 0 (no data) is not in the table: such points take surface.roughness_map.default
  or are refused (amrTerrain: 0.1 m);
- 12 (tree or shrub cover in cropland) is 0.25 m, rough, above 11 (herbaceous
  cover in cropland, 0.1 m) (amrTerrain: 0.07 m);
- 160 and 170 (tree cover, flooded) are 1.0 m, closed, as the other tree
  classes (amrTerrain: 0.2 m);
- 200, 201, 202 (bare areas) are 0.005 m, smooth (amrTerrain: 0.05 m).

It is a working table, not a standard; give your own with
``surface.roughness_map.classes`` when you have one.
"""
import os

import numpy as np

from . import terrain as terrain_mod
from .errors import CaseError

ESA_CCI_ROUGHNESS = {
    10: 0.1, 11: 0.1, 12: 0.25, 20: 0.2, 30: 0.3, 40: 0.3,
    50: 1.0, 60: 1.0, 61: 1.0, 62: 1.0, 70: 1.0, 71: 1.0, 72: 1.0, 80: 1.0, 81: 1.0, 82: 1.0, 90: 1.0,
    100: 0.15, 110: 0.3, 120: 0.07, 121: 0.07, 122: 0.07, 130: 0.07,
    140: 0.05, 150: 0.05, 151: 0.05, 152: 0.05, 153: 0.05, 160: 1.0, 170: 1.0, 180: 0.2,
    190: 0.8, 200: 0.005, 201: 0.005, 202: 0.005, 210: 1.0e-4, 220: 0.0012,
}


def class_table(classes, key="surface.roughness_map.classes"):
    """The {class code: z0} table of ``classes``: 'esa_cci' or a mapping."""
    if classes == "esa_cci":
        return dict(ESA_CCI_ROUGHNESS)
    if not isinstance(classes, dict) or not classes:
        raise CaseError("%s: must be esa_cci or a mapping of class code to roughness length (m)" % key)
    table = {}
    for code, z0 in classes.items():
        try:
            code_i, z0_f = int(code), float(z0)
        except (TypeError, ValueError):
            raise CaseError("%s: '%s: %s' is not a class code and a roughness length" % (key, code, z0)) from None
        if not (np.isfinite(z0_f) and z0_f > 0.0):
            raise CaseError("%s: class %d has z0 = %s; every roughness length must be finite and > 0" % (key, code_i, z0))
        table[code_i] = z0_f
    return table


def classes_to_roughness(codes, table, default=None, key="surface.roughness_map"):
    """z0 of an array of class codes. A code the table lacks takes ``default``, or is refused."""
    codes = np.asarray(codes)
    if not np.all(np.isfinite(codes)):
        raise CaseError("%s.file: %d points in the domain have no class" % (key, int(np.sum(~np.isfinite(codes)))))
    icodes = np.rint(codes).astype(int)
    z0 = np.full(icodes.shape, np.nan)
    for code, value in table.items():
        z0[icodes == code] = value
    missing = sorted(set(icodes[np.isnan(z0)].tolist()))
    if missing:
        if default is None:
            raise CaseError("%s.file holds the classes %s, which %s.classes does not; add them or give %s.default"
                            % (key, ", ".join(str(m) for m in missing), key, key))
        z0[np.isnan(z0)] = default
    return z0


class RoughnessMap:
    def __init__(self, x, y, z0, notes):
        self.x, self.y, self.z0 = x, y, z0
        self.notes = notes

    def write(self, path):
        with open(path, "w", encoding="utf-8") as out:
            out.write("# x y z0 (m), case frame; %d points\n" % self.x.size)
            np.savetxt(out, np.column_stack([self.x, self.y, self.z0]), fmt=["%.3f", "%.3f", "%.6g"])


def prepare(path, fmt, coordinates, classes, default, frame, margin, point_spacing, key="surface.roughness_map"):
    """The roughness map of the domain, in the case frame."""
    file_key = key + ".file"
    if not os.path.isfile(path):
        raise CaseError("%s = %s: no such file" % (file_key, path))
    if fmt == "auto":
        fmt = "xyz" if terrain_mod.guess_format(path, file_key) == "xyz" else "land_cover"
    if fmt == "xyz":
        px, py, z0 = terrain_mod.read_points(path, file_key)
        geographic = coordinates == "geographic" or (
            coordinates == "auto" and bool(np.all(np.abs(px) <= 180.0) and np.all(np.abs(py) <= 90.0)))
        if geographic:
            if frame.lat_lon is None:
                raise CaseError(terrain_mod.geographic_message(file_key, path))
            cx, cy = frame.from_lat_lon(py, px)
        else:
            cx, cy = frame.to_case(px, py)
        keep = ((cx >= frame.x_lo - margin) & (cx <= frame.x_hi + margin) &
                (cy >= frame.y_lo - margin) & (cy <= frame.y_hi + margin))
        cx, cy, z0 = cx[keep], cy[keep], z0[keep]
        notes = ["%s: %d roughness points kept" % (os.path.basename(path), cx.size)]
    else:
        raster = terrain_mod.read_raster(path, terrain_mod.guess_format(path, file_key), file_key)
        geographic = coordinates == "geographic" or (coordinates == "auto" and raster.is_geographic())
        cx, cy, codes, spacing = terrain_mod.raster_points(raster, frame, geographic, margin, file_key,
                                                           stride_spacing=point_spacing)
        z0 = classes_to_roughness(codes, class_table(classes, key + ".classes"), default, key)
        notes = ["%s: land cover at about %.4g m, %d points, z0 %.4g .. %.4g m"
                 % (os.path.basename(path), spacing, cx.size, z0.min(), z0.max())]
    if cx.size == 0:
        raise CaseError("%s = %s: no points inside the domain" % (file_key, path))
    if not (np.all(np.isfinite(z0)) and np.all(z0 > 0.0)):
        raise CaseError("%s = %s: every roughness length must be finite and > 0 (lowest %s)" % (file_key, path, np.nanmin(z0)))
    return RoughnessMap(cx, cy, z0, notes)
