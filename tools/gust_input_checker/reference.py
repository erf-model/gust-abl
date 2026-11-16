"""Reader for Tests/inputs_reference, the list of every GUST input.

The file is generated from docs/inputs.rst (tools/input_reference.py --write)
and the CTest input_reference_current keeps it equal to the keys read in
Source/, so it is the authority for which keys exist and what their defaults
are. Its layout:

    # ====...
    # mesh                      <- section banner: the namespace
    # ====...
    #
    # mesh.n_cell =    # required
    #     Three integers >= 1: cells in x, y and z.
    #
    # bc.inflow_ustar = 0.4    # m/s
    #     ...

Several key lines may share one indented text block. A note of ``required``
means the key has no default; ``none`` (or a dash) that it is off or empty
unless given.
"""
import os
import re

_KEY = re.compile(r"^# ([A-Za-z_][A-Za-z0-9_]*\.[A-Za-z0-9_]+) =(.*)$")
_BANNER = re.compile(r"^# =+$")


class RefKey:
    """One key of the reference."""

    def __init__(self, key, default, note, section, line):
        self.key = key
        self.default = default      # the text after '=', before the note ('' when there is none)
        self.note = note            # the text after the '#' of the key line
        self.section = section
        self.line = line
        self.text = ""              # the indented block below the key line(s)

    @property
    def namespace(self):
        return self.key.split(".", 1)[0]

    @property
    def required(self):
        return self.note.startswith("required")

    def __repr__(self):
        return "RefKey(%s = %s)" % (self.key, self.default)


class Reference:
    """Every key of Tests/inputs_reference."""

    def __init__(self, path=""):
        self.path = path
        self.keys = {}
        self.sections = []

    def has(self, key):
        return key in self.keys

    def default(self, key):
        """The default as a list of tokens ([] when the key has none or is unknown)."""
        ref = self.keys.get(key)
        return ref.default.split() if ref is not None else []

    def namespace_keys(self, namespace):
        return [k for k in self.keys if k.startswith(namespace + ".")]

    def text(self, key):
        return self.keys[key].text if key in self.keys else ""


def default_reference_path():
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, "..", "..", "Tests", "inputs_reference"))


def parse_reference_text(text, path=""):
    ref = Reference(path)
    lines = text.splitlines()
    section = ""
    group = []      # the key lines waiting for their text block
    block = []
    i = 0

    def close():
        body = " ".join(s.strip() for s in block).strip()
        for k in group:
            k.text = body
        del group[:]
        del block[:]

    while i < len(lines):
        line = lines[i]
        if _BANNER.match(line) and i + 2 < len(lines) and _BANNER.match(lines[i + 2]):
            close()
            section = lines[i + 1][2:].strip()
            ref.sections.append(section)
            i += 3
            continue
        m = _KEY.match(line)
        if m:
            if block:
                close()
            rest = m.group(2)
            default, _, note = rest.partition("#")
            k = RefKey(m.group(1), default.strip(), note.strip(), section, i + 1)
            ref.keys[k.key] = k
            group.append(k)
        elif line.startswith("#     ") and group:
            block.append(line[1:])
        elif group and block:
            close()
        i += 1
    close()
    return ref


_CACHE = {}


def load_reference(path=None):
    """The reference at `path` (default: Tests/inputs_reference of this repository), read once per path."""
    path = os.path.abspath(path or default_reference_path())
    if path not in _CACHE:
        with open(path, "r", encoding="utf-8") as fh:
            _CACHE[path] = parse_reference_text(fh.read(), path)
    return _CACHE[path]
