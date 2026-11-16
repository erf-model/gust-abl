"""Case set-up for GUST: a YAML description in, a case folder out.

    from gust_case import load_case, build
    case, given = load_case("case.yaml")
    result = build(case, given, base_dir=".")

``tools/make_gust_case.py`` is the command-line front end; docs/case_setup.rst
describes the keys.
"""
from .builder import Result, build
from .errors import CaseError
from .plan import resolve
from .schema import SPEC, apply_override, key_table_rst, key_table_text, load_yaml, validate


def load_case(path, overrides=()):
    """Read and validate a YAML case file; ``overrides`` are section.key=value strings."""
    raw = load_yaml(path)
    for assignment in overrides:
        apply_override(raw, assignment)
    return validate(raw)


__all__ = ["CaseError", "Result", "SPEC", "apply_override", "build", "key_table_rst", "key_table_text", "load_case",
           "load_yaml", "resolve", "validate"]
