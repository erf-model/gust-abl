"""Findings and the report of one deck.

Four severities, as in ERF's checker (Exec/Tools/erf_input_checker):

    ERROR    GUST would stop at start-up with this deck
    WARNING  GUST starts, but the run is very likely not what was meant
    ADVICE   the run is valid; a setting could be chosen better
    INFO     numbers derived from the deck, for the user to read
"""
import enum
import json


class Severity(enum.IntEnum):
    INFO = 0
    ADVICE = 1
    WARNING = 2
    ERROR = 3

    @property
    def label(self):
        return self.name


class Finding:
    """One result of one check."""

    def __init__(self, severity, checker, key, message, advice="", where=""):
        self.severity = Severity(severity)
        self.checker = checker
        self.key = key
        self.message = message
        self.advice = advice
        self.where = where

    def as_dict(self):
        return {"severity": self.severity.label, "checker": self.checker, "key": self.key,
                "message": self.message, "advice": self.advice, "where": self.where}

    def format(self):
        head = "[%s] %s" % (self.severity.label, self.key)
        if self.where:
            head += " (%s)" % self.where
        head += " <%s>" % self.checker
        lines = [head, "    " + self.message]
        if self.advice:
            lines.append("    -> " + self.advice)
        return "\n".join(lines)

    def __repr__(self):
        return "Finding(%s, %s, %s)" % (self.severity.label, self.checker, self.key)


class Report:
    """The findings of one deck."""

    def __init__(self, path=""):
        self.path = path
        self.findings = []

    def add(self, finding):
        self.findings.append(finding)

    def count(self, severity):
        return sum(1 for f in self.findings if f.severity == severity)

    def has_errors(self):
        return self.count(Severity.ERROR) > 0

    @property
    def worst(self):
        """The highest severity in the report, or None when it is empty."""
        return max((f.severity for f in self.findings), default=None)

    def select(self, severity=None, key=None, checker=None):
        """The findings that match every argument given (used by the tests)."""
        return [f for f in self.findings
                if (severity is None or f.severity == severity)
                and (key is None or f.key == key)
                and (checker is None or f.checker == checker)]

    def sorted(self, min_severity=Severity.INFO):
        keep = [f for f in self.findings if f.severity >= min_severity]
        return sorted(keep, key=lambda f: (-int(f.severity), f.checker, f.key))

    def format(self, min_severity=Severity.INFO):
        out = ["== %s" % self.path]
        for f in self.sorted(min_severity):
            out.append(f.format())
        out.append("-- %d error(s), %d warning(s), %d advice, %d info"
                   % (self.count(Severity.ERROR), self.count(Severity.WARNING),
                      self.count(Severity.ADVICE), self.count(Severity.INFO)))
        return "\n".join(out)

    def as_json(self):
        return json.dumps({"path": self.path, "findings": [f.as_dict() for f in self.sorted()]}, indent=2)
