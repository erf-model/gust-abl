"""A static, pre-run checker of GUST input decks.

    from gust_input_checker import check_deck
    report = check_deck("inputs", overrides=["steady.max_iterations=500"])
    if report.has_errors():
        print(report.format())

The checker reads a deck as GUST's ParmParse does, derives the facts of the
run (stage, mesh, layers, faces) once, and runs one check class per area. A
finding is an ERROR (GUST would stop at start-up), a WARNING (GUST starts, but
the deck is very likely not what was meant), an ADVICE (a setting that could
be chosen better) or an INFO (a number derived from the deck). The known keys
and their defaults come from Tests/inputs_reference.

The method follows ERF's input checker (Exec/Tools/check_erf_inputs.py).
"""
import os

from .findings import Severity, Finding, Report
from .deck import Deck, DeckError, parse_deck, parse_text
from .context import RunContext, Options
from .checker import Checker, register, registry, run_all_checks
from . import checks  # noqa: F401  (registers the checkers)

__all__ = ["Severity", "Finding", "Report", "Deck", "parse_deck", "parse_text", "RunContext", "Options",
           "Checker", "register", "registry", "run_all_checks", "check_deck", "check_text"]


def _run(deck, path, options, only):
    report = Report(path)
    try:
        ctx = RunContext(deck, options)
    except Exception as exc:    # a checker bug must not hide the deck's own problems
        report.add(Finding(Severity.ERROR, "context", "internal",
                           "The checker could not derive the run's facts (%s: %s)." % (type(exc).__name__, exc),
                           "Report this with the deck.", ""))
        return report
    run_all_checks(deck, ctx, report, only=only)
    return report


def check_deck(path, overrides=(), options=None, only=None):
    """Parse `path` (plus `key=value` strings in `overrides`) and run the checks.

    Never raises on a bad deck: a parse failure is an ERROR finding. `options` is an Options (ranks, gpu,
    reference, cwd); `only` a list of checker names. Returns a Report.
    """
    path = os.fspath(path)
    try:
        deck = parse_deck(path, overrides=overrides, cwd=options.cwd if options else None)
    except (OSError, UnicodeError, DeckError) as exc:
        report = Report(path)
        report.add(Finding(Severity.ERROR, "deck", "deck", "The deck cannot be read: %s." % exc,
                           "Give the path of an input file.", path))
        return report
    return _run(deck, path, options, only)


def check_text(text, path="<text>", overrides=(), options=None, only=None):
    """As check_deck, for a deck held in a string (`path` only names it and anchors relative files)."""
    try:
        deck = parse_text(text, path=path, overrides=overrides, cwd=options.cwd if options else None)
    except DeckError as exc:
        report = Report(path)
        report.add(Finding(Severity.ERROR, "deck", "deck", "The deck cannot be read: %s." % exc, "", path))
        return report
    return _run(deck, path, options, only)
