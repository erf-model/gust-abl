"""How the form shows the keys (no Qt): plain-language labels, the groups of a tab, and the hints that say what to
type where the key's kind alone does not. The YAML keys do not change; each field's tooltip names its key.
"""

# label of a field; a key not listed shows its own name. Nest keys are listed as nests.<key>.
LABELS = {
    "case.name": "Case name",
    "case.output_dir": "Folder the case is written in",
    "case.description": "Description",
    "domain.center_lat_lon": "Centre",
    "domain.utm_zone": "UTM zone",
    "domain.center": "Centre: easting, northing",
    "domain.half_width": "Half-width (centre to each face)",
    "domain.west": "Centre to west face",
    "domain.east": "Centre to east face",
    "domain.south": "Centre to south face",
    "domain.north": "Centre to north face",
    "domain.top": "Height of the top",
    "domain.origin": "x = 0, y = 0 of the case at",
    "terrain.file": "Elevation file",
    "mesh.cell_size": "Cell size",
    "mesh.max_grid_size": "Largest box (cells)",
    "mesh.hierarchy_levels": "Coarser grid levels",
    "mesh.layers": "Number of layers",
    "mesh.first_layer": "First layer thickness",
    "mesh.stretch_ratio": "Growth ratio between layers",
    "wind.speed": "Wind speed",
    "wind.height": "at height",
    "wind.direction": "Wind from (degrees)",
    "wind.directions": "Several directions (degrees)",
    "surface.roughness": "Roughness length",
    "run.ranks": "MPI ranks",
    "run.executable": "gust program",
    "run.mpi_launcher": "MPI launcher",
    "nests.name": "Name",
    "nests.center": "Centre: x, y",
    "nests.center_lat_lon": "Centre",
    "nests.half_width": "Half-width (centre to each face)",
    "nests.blend_width": "Blend width",
    "nests.max_iterations": "Most iterations",
}

# a tab's groups: (title, note, keys); keys not listed follow in the schema's order after the groups
GROUPS = {
    "domain": [
        ("Where on the Earth", "The centre of the domain (WGS 84). Needed for terrain in latitude and longitude "
         "and for Download terrain.", ["center_lat_lon", "utm_zone"]),
        ("Or: where in the terrain file's own coordinates", "For terrain already in metres (projected). Empty "
         "means 0, 0.", ["center"]),
        ("Size", "Half-width: one number for a square domain, or two (east-west, north-south).",
         ["half_width", "west", "east", "south", "north", "top"]),
        ("Case coordinates", None, ["origin"]),
    ],
    "mesh": [
        ("Horizontal", "The cell counts are rounded to multiples of 8.",
         ["cell_size", "hierarchy_levels", "max_grid_size"]),
        ("Vertical", "With a first layer thickness the layers grow from the ground to the top; without it they are "
         "uniform.", ["layers", "first_layer", "stretch_ratio"]),
    ],
}

# the shape of a value where the kind alone does not say it (key name without its section)
SHAPES = {
    "center": "x, y",
    "half_width": "one number, or east-west, north-south",
    "geostrophic_wind": "u, v",
    "heights": "list of numbers, e.g. 0, 500, 1000",
    "fields": "list of words, e.g. theta",
    "directions": "list, e.g. [210, 240]",
    "initial_blend": "2 numbers > 0",
    "executable": "program name or path",
}
SHAPES_FULL = {
    "domain.center": "easting, northing",
    "atmosphere.sounding.theta": "list of numbers, e.g. 300, 303, 310",
    "domain.utm_zone": "1 to 60; empty: from the longitude",
}

# keys shown as two boxes (the list's two entries) with their own hints
PAIRS = {"center_lat_lon": ("latitude (-90 to 90)", "longitude (-180 to 180)")}

# the basic mode: the keys a case needs; every other key keeps its default (the builder works it out). A key a
# file gives, or one the checker refuses, shows in the basic mode too.
BASIC = {
    "case.name", "case.output_dir",
    "domain.center_lat_lon", "domain.half_width", "domain.top",
    "terrain.file",
    "mesh.cell_size", "mesh.layers", "mesh.first_layer",
    "wind.speed", "wind.height", "wind.direction",
    "surface.roughness",
    "run.ranks",
}

FOLDERS = {"case.output_dir"}       # path keys that name a folder; the other path keys name a file


def label(name, key):
    """The field's label: the plain-language one, with the unit."""
    short = name.rsplit(".", 1)[-1]
    nest = name.startswith("nests")
    text = LABELS.get("nests." + short if nest else name, short)
    if key.unit and short not in ("direction", "directions"):
        text += " (%s)" % key.unit
    return text


def groups(section, keys):
    """[(title, note, [keys])] of a tab: its GROUPS, then a last untitled group of the keys they leave out."""
    listed = GROUPS.get(section, [])
    taken = {k for _, _, ks in listed for k in ks}
    rest = [k for k in keys if k not in taken]
    out = [(t, n, [k for k in ks if k in keys]) for t, n, ks in listed]
    if rest:
        out.append((None, None, rest))
    return out
