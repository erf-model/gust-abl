"""Text of a form field <-> value of a case-description key (no Qt).

A field's text is what the YAML file holds for the key, written as YAML in one line; a list of numbers or words
is written ``1, 2`` (``[1, 2]`` and ``1 2`` are read too). An empty field is a key not given, so the builder's default applies. ``string``, ``path`` and ``choice`` keys take the
text as it is (``name: 012`` stays the text 012). ``gust_inputs`` is a block of YAML lines.
"""
import re

import yaml

from gust_case.schema import NestList, SPEC

# a key named in a refusal: section.key, nests[2].key, gust_inputs.namespace.key
_LIST_SEP = re.compile(r"\s*,\s*|\s+")         # between the entries of a list written without brackets
_KEY_AT_START = re.compile(r"^([A-Za-z_]\w*(?:\[\d+\])?(?:\.[A-Za-z_]\w*)*)")


def _flow(value):
    text = yaml.safe_dump(value, default_flow_style=True, width=float("inf"), sort_keys=False).strip()
    if text.endswith("\n..."):                    # a scalar ends with the end-of-document mark
        text = text[:-4].rstrip()
    return text


def to_text(key, value):
    """The field text of ``value`` (None: not given, empty text)."""
    if value is None:
        return ""
    if key.kind == "any" and not isinstance(value, (list, dict)):
        return _flow(value)
    if key.kind == "any":
        return yaml.safe_dump(value, default_flow_style=False, sort_keys=False).rstrip()
    if key.kind in ("string", "path", "choice") and isinstance(value, str):
        return value
    if key.kind in ("numbers", "strings") and isinstance(value, (list, tuple)) and value \
            and all(not isinstance(v, (list, tuple, dict)) for v in value):
        return ", ".join(_flow(v) for v in value)
    return _flow(value)


def from_text(name, key, text):
    """(given, value) of a field's text; a text that is not YAML raises ValueError naming the key."""
    if key.kind == "any":
        if not text.strip():
            return False, None
        try:
            value = yaml.safe_load(text)
        except yaml.YAMLError as err:
            raise ValueError("%s: not valid YAML: %s" % (name, " ".join(str(err).split()))) from None
        return value is not None, value
    text = text.strip()
    if not text:
        return False, None
    if key.kind in ("string", "path", "choice"):
        return True, text
    if key.kind in ("numbers", "strings") and text[0] not in "[{" and _LIST_SEP.search(text):
        text = "[" + ", ".join(t for t in _LIST_SEP.split(text) if t) + "]"
    try:
        value = yaml.safe_load(text)
    except yaml.YAMLError as err:
        raise ValueError("%s = %s: not valid YAML: %s" % (name, text, " ".join(str(err).split()))) from None
    return value is not None, value


def field_names(spec=None, prefix=""):
    """Names of the form's fields outside the nests, in the schema's order (schema.iter_keys without nests[])."""
    for name, entry in (SPEC if spec is None else spec).items():
        if isinstance(entry, dict):
            yield from field_names(entry, prefix + name + ".")
        elif not isinstance(entry, NestList):
            yield prefix + name, entry


def nest_spec():
    return SPEC["nests"].spec


def get_path(raw, dotted):
    """The value at ``a.b.c`` of nested dicts, or None."""
    node = raw
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            return None
        node = node[part]
    return node


def set_path(raw, dotted, value):
    parts = dotted.split(".")
    node = raw
    for part in parts[:-1]:
        node = node.setdefault(part, {})
    node[parts[-1]] = value


def unknown_parts(raw, spec=None):
    """The entries of ``raw`` the schema does not know, as nested dicts (kept when the form saves the file).

    The nests list is the form's; unknown keys inside a nest stay in that nest's own extra entries."""
    spec = SPEC if spec is None else spec
    out = {}
    if not isinstance(raw, dict):
        return out
    for name, value in raw.items():
        entry = spec.get(name) if isinstance(name, str) else None
        if entry is None:
            out[name] = value
        elif isinstance(entry, dict):
            if isinstance(value, dict):
                sub = unknown_parts(value, entry)
                if sub:
                    out[name] = sub
            elif value is not None:
                out[name] = value                 # a section given as a scalar: kept, the builder refuses it
    return out


def merge(into, extra):
    for name, value in extra.items():
        if isinstance(value, dict) and isinstance(into.get(name), dict):
            merge(into[name], value)
        else:
            into[name] = value
    return into


def problem_field(problem, names):
    """The field a refusal names: the longest dotted prefix of the problem's first word that is a field.

    ``mesh.layers = 0: must be ...`` -> mesh.layers; ``nests[1].half_width: required`` -> nests[1].half_width;
    ``gust_inputs.steady.foo ...`` -> gust_inputs. None when the problem names no field."""
    match = _KEY_AT_START.match(problem)
    if not match:
        return None
    parts = match.group(1).split(".")
    for n in range(len(parts), 0, -1):
        candidate = ".".join(parts[:n])
        if candidate in names:
            return candidate
    return None
