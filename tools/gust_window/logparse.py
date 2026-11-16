"""Reading GUST's log as it is written (no Qt): the iteration table and the column run's geostrophic wind.

The table (Source/Sequencing.cpp, IterationPrinter::print) is

      Iter        CFL  Line solver      u-res      v-res      w-res    tke-res  theta-res      p-res          r
    ------  ---------  -----------  ---------  ...
         1  1.000e+00  lusgs        3.525e-03  ...                                                   1.000e+00

with a ``level N: `` prefix on the coarse levels of grid sequencing, ``-`` for a residual a run does not solve
(and for the CFL of a Newton row without a step, ``inf`` for an infinite one) and notes after r.
"""
import math
import re

RESIDUALS = ("u", "v", "w", "tke", "theta", "p")
_LEVEL = re.compile(r"^level (\d+): ")


def _number(text):
    """A cell: ``-`` (not solved) and every spelling of NaN (MSVC writes ``nan(ind)``, ``-nan(ind)``) are NaN."""
    if text == "-" or text.lstrip("+-").lower().startswith("nan"):
        return math.nan
    return float(text)                            # inf, -inf, 1.000e+100


class Row:
    """One row of the iteration table."""

    def __init__(self, level, iteration, cfl, solver, residuals, r, notes, x):
        self.level = level            # the grid-sequencing level (None: the finest)
        self.iteration = iteration
        self.cfl = cfl
        self.solver = solver
        self.residuals = residuals    # name -> value (NaN when not solved)
        self.r = r
        self.notes = notes
        self.x = x                    # iterations since the run started, across the levels


def parse_row(line):
    """(level, iteration, cfl, solver, residuals, r, notes) of a table row, or None for any other line."""
    text = line.rstrip("\n")
    level = None
    m = _LEVEL.match(text)
    if m:
        level = int(m.group(1))
        text = text[m.end():]
    tokens = text.split()
    if len(tokens) < 10 or not tokens[0].isdigit() or not re.fullmatch(r"[a-z_]+", tokens[2]):
        return None
    try:
        cfl = _number(tokens[1])
        values = [_number(t) for t in tokens[3:10]]
    except ValueError:
        return None
    notes = text.split(None, 10)[10] if len(tokens) > 10 else ""
    return level, int(tokens[0]), cfl, tokens[2], dict(zip(RESIDUALS, values[:6])), values[6], notes


def is_header(line):
    return "Iter" in line and "Line solver" in line and "u-res" in line


class LogReader:
    """Feed it the log's text as it arrives; it returns the table rows of each piece.

    A row counts only after a header (so no other line of the log is read as a row). The x of a row runs on
    across the levels of grid sequencing: a new level starts where the previous one stopped."""

    def __init__(self):
        self._partial = ""
        self._seen_header = False
        self._offset = 0
        self._last = 0
        self._level = "start"
        self.rows = []

    def feed(self, text):
        lines = (self._partial + text).split("\n")
        self._partial = lines.pop()
        new = []
        for line in lines:
            line = line.rstrip("\r")
            if is_header(line):
                self._seen_header = True
                continue
            if not self._seen_header:
                continue
            parsed = parse_row(line)
            if parsed is None:
                continue
            level, iteration = parsed[0], parsed[1]
            if level != self._level:
                if self._level != "start":
                    self._offset += self._last
                self._level = level
            self._last = iteration
            row = Row(*parsed, x=self._offset + iteration)
            self.rows.append(row)
            new.append(row)
        return new


def column_wind(text):
    """(ug, vg) as text from the column run's log (the last ``forcing.ug = ...`` and ``forcing.vg = ...``), as
    run.sh reads them; None when either is missing."""
    found = {}
    for line in text.splitlines():
        tokens = line.split()
        if len(tokens) >= 3 and tokens[0] in ("forcing.ug", "forcing.vg") and tokens[1] == "=":
            found[tokens[0]] = tokens[2]
    if "forcing.ug" in found and "forcing.vg" in found:
        return found["forcing.ug"], found["forcing.vg"]
    return None
