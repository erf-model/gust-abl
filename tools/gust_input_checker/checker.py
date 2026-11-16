"""The base class of the checks, their registry and the loop that runs them.

One class per area in ``checks/``; a class registers itself with ``@register``
and reports through ``error``, ``warn``, ``advise`` and ``info``. A check that
fails on a deck (a bug in the check) does not hide the rest of the report: the
failure becomes one ERROR finding with the key ``internal``.
"""
from .deck import DeckError
from .findings import Finding, Severity

_REGISTRY = []


def register(cls):
    """Class decorator: add a Checker to the list that run_all_checks runs."""
    if not getattr(cls, "name", ""):
        raise ValueError("a Checker needs a name")
    if any(c.name == cls.name for c in _REGISTRY):
        raise ValueError("two checkers are named " + cls.name)
    _REGISTRY.append(cls)
    return cls


def registry():
    return list(_REGISTRY)


class Checker:
    """One area of checks. Subclasses set `name` and write `run`."""

    name = ""

    def __init__(self, deck, ctx, report):
        self.deck = deck
        self.ctx = ctx
        self.report = report

    @classmethod
    def applies(cls, ctx):
        return True

    def run(self):
        raise NotImplementedError

    # -- reporting --------------------------------------------------------
    def _add(self, severity, key, message, advice):
        where = self.deck.where(key) if self.deck.has(key) else ""
        self.report.add(Finding(severity, self.name, key, message, advice, where))

    def error(self, key, message, advice=""):
        self._add(Severity.ERROR, key, message, advice)

    def warn(self, key, message, advice=""):
        self._add(Severity.WARNING, key, message, advice)

    def advise(self, key, message, advice=""):
        self._add(Severity.ADVICE, key, message, advice)

    def info(self, key, message, advice=""):
        self._add(Severity.INFO, key, message, advice)

    # -- shorthands on the context ---------------------------------------
    def given(self, key):
        return self.ctx.given(key)

    def no_effect(self, key, why, advice=""):
        """ERROR when `key` is given: GUST refuses a key that has no effect."""
        if self.given(key):
            self.error(key, "GUST stops: %s has no effect %s" % (key, why), advice or "Remove %s." % key)
            return True
        return False


def run_all_checks(deck, ctx, report, only=None):
    """Run every registered check that applies (or only those named in `only`)."""
    for where, message in deck.problems:
        report.add(Finding(Severity.ERROR, "deck", "deck", "AMReX stops while it reads the inputs: " + message,
                           "Write one 'namespace.key = value' definition per line.", where))
    for key, message in ctx.problems:
        report.add(Finding(Severity.ERROR, "deck", key, message, "", deck.where(key) if deck.has(key) else ""))
    for cls in registry():
        if only and cls.name not in only:
            continue
        try:
            if not cls.applies(ctx):
                continue
            cls(deck, ctx, report).run()
        except DeckError as e:
            report.add(Finding(Severity.ERROR, cls.name, e.key, e.message, "",
                               deck.where(e.key) if deck.has(e.key) else ""))
        except Exception as e:   # a bug in a check must not hide the other findings
            report.add(Finding(Severity.ERROR, cls.name, "internal",
                               "the %s checker failed on this deck: %s: %s" % (cls.name, type(e).__name__, e),
                               "Please report this with the deck; the other checks still ran."))
    return report
