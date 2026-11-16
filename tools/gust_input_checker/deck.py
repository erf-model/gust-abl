"""Reader for GUST input decks (AMReX ParmParse syntax).

A deck is a list of definitions ``namespace.key = value ...``, one per line;
``#`` starts a comment, a trailing ``\\`` continues the line and double quotes
keep blanks inside one value. Every definition is kept with its file and line,
so a finding can say where the key is set.

What the reader follows, each checked against the gust executable:

- definitions given on the command line are added after the file's (AMReX
  appends them to the same table); a key the file already sets is then defined
  twice, and GUST refuses a key of its own namespaces that is given more than
  once (docs/running.rst, "a repeated key is refused");
- a line that holds values but no ``key =`` stops AMReX ("Must use \\ for line
  continuation"), as does ``key =`` without a value ("no values for
  definition");
- numbers go through AMReX's parser: ``2*4`` is 8 and ``8.0`` is a valid
  integer, ``8.5`` is not, and Fortran's ``1.d0`` is refused.
"""
import ast
import math
import operator
import os
import re

COMMAND_LINE = "command line"


class DeckError(Exception):
    """A value of the deck cannot be read as the type a check asked for."""

    def __init__(self, key, message):
        Exception.__init__(self, message)
        self.key = key
        self.message = message


class Unparsed(Exception):
    """A value is an expression this reader does not evaluate (AMReX may)."""


_BIN = {ast.Add: operator.add, ast.Sub: operator.sub, ast.Mult: operator.mul, ast.Div: operator.truediv,
        ast.Pow: operator.pow}
_UNA = {ast.UAdd: operator.pos, ast.USub: operator.neg}
_WORD = re.compile(r"^[A-Za-z_][A-Za-z0-9_.\-/]*$")


def _eval(node):
    if isinstance(node, ast.Expression):
        return _eval(node.body)
    if isinstance(node, ast.Constant) and isinstance(node.value, (int, float)) and not isinstance(node.value, bool):
        return float(node.value)
    if isinstance(node, ast.BinOp) and type(node.op) in _BIN:
        return _BIN[type(node.op)](_eval(node.left), _eval(node.right))
    if isinstance(node, ast.UnaryOp) and type(node.op) in _UNA:
        return _UNA[type(node.op)](_eval(node.operand))
    raise Unparsed()


def to_real(token):
    """The number AMReX reads from `token`.

    Raises ValueError when AMReX would stop on it (a word, Fortran's d
    exponent) and Unparsed for an expression this reader cannot evaluate.
    """
    text = token.strip()
    try:
        if re.match(r"^[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?$", text):
            return float(text)
    except ValueError:
        pass
    if text.lower() in ("nan", "inf", "+inf", "-inf", "infinity"):
        raise Unparsed()
    if _WORD.match(text) or re.match(r"^[+-]?(\d+\.?\d*|\.\d+)[dD][+-]?\d+$", text) or not text:
        raise ValueError(token)
    try:
        return float(_eval(ast.parse(text.replace("^", "**"), mode="eval")))
    except (SyntaxError, ValueError, ZeroDivisionError, OverflowError, TypeError, RecursionError):
        raise Unparsed()


def to_int(token):
    """The integer AMReX reads from `token` (8, 8.0 and 2*4 are integers; 8.5 is not)."""
    x = to_real(token)
    if not math.isfinite(x) or x != math.floor(x):
        raise ValueError(token)
    return int(x)


class Entry:
    """One definition: a key, its values, and where it was given."""

    def __init__(self, key, values, path, line):
        self.key = key
        self.values = list(values)
        self.path = path
        self.line = line

    @property
    def where(self):
        if self.path == COMMAND_LINE:
            return COMMAND_LINE
        return "%s:%d" % (os.path.basename(self.path), self.line)

    def __repr__(self):
        return "Entry(%s = %s, %s)" % (self.key, " ".join(self.values), self.where)


def _tokens(line):
    """Split one logical line into tokens; '=' is a token of its own, quotes keep blanks."""
    out, cur, quoted, i = [], None, False, 0
    while i < len(line):
        c = line[i]
        if quoted:
            if c == '"':
                quoted = False
            else:
                cur += c
        elif c == '"':
            quoted = True
            cur = "" if cur is None else cur
        elif c == "#":
            break
        elif c in " \t\r\n\v\f":
            if cur is not None:
                out.append(cur)
                cur = None
        elif c == "=":
            if cur is not None:
                out.append(cur)
                cur = None
            out.append(("=",))
        else:
            cur = c if cur is None else cur + c
        i += 1
    if cur is not None:
        out.append(cur)
    return out


class Deck:
    """The definitions of one deck, in the order GUST sees them."""

    def __init__(self, path=""):
        self.path = path
        self.entries = []
        self.table = {}
        self.problems = []   # (where, message): what stops AMReX while it reads the deck
        self.cwd = None      # where a relative FILE = path is opened (None: the current directory), as AMReX does
        self.included = []   # the files FILE = read, in order

    # -- building ---------------------------------------------------------
    def add(self, entry):
        self.entries.append(entry)
        self.table[entry.key] = entry   # a query returns the last definition

    def read_text(self, text, path):
        pending, start = "", 0
        for lineno, raw in enumerate(text.splitlines(), 1):
            line = raw
            hash_at = _comment_start(line)
            if hash_at >= 0:
                line = line[:hash_at]
            if not pending:
                start = lineno
            if line.rstrip().endswith("\\"):
                pending += line.rstrip()[:-1] + " "
                continue
            self._read_line(pending + line, path, start)
            pending = ""
        if pending.strip():
            self._read_line(pending, path, start)

    def read_command_line(self, items):
        """Add command-line definitions: AMReX joins the arguments with blanks and reads them as one line."""
        line = " ".join(str(i) for i in items)
        if line.strip():
            self._read_line(line, COMMAND_LINE, 0)

    def _read_line(self, line, path, lineno):
        toks = _tokens(line)
        if not toks:
            return
        where = "%s:%d" % (os.path.basename(path), lineno) if path != COMMAND_LINE else COMMAND_LINE
        eq = [i for i, t in enumerate(toks) if t == ("=",)]
        if not eq or eq[0] != 1:
            shown = " ".join(t for t in toks if t != ("=",))
            self.problems.append((where, "'%s' is not a 'key = value' definition; AMReX stops on a line that does not "
                                         "start with 'key =' (a trailing \\ continues a line)" % shown[:60]))
            return
        for n, i in enumerate(eq):
            key = toks[i - 1]
            end = eq[n + 1] - 1 if n + 1 < len(eq) else len(toks)
            values = [t for t in toks[i + 1:end] if t != ("=",)]
            if key == ("=",) or (n + 1 < len(eq) and eq[n + 1] == i + 1):
                self.problems.append((where, "'=' without a key before it"))
                continue
            if not values:
                self.problems.append((where, "%s has no value; AMReX stops (\"no values for definition\")" % key))
                continue
            self.add(Entry(key, values, path, lineno))
            if key == "FILE":
                self._include(values[0], where)

    def _include(self, name, where):
        """FILE = name: AMReX reads the file's definitions in place, the path taken from the run's directory."""
        target = name if os.path.isabs(name) else os.path.join(self.cwd or os.getcwd(), name)
        target = os.path.normpath(target)
        if target in self.included or len(self.included) >= 20:
            self.problems.append((where, "FILE = %s includes itself (directly or through other files)" % name))
            return
        try:
            with open(target, "r", encoding="utf-8", errors="replace") as fh:
                text = fh.read()
        except OSError:
            self.problems.append((where, "FILE = %s cannot be opened from %s; AMReX stops"
                                         % (name, self.cwd or os.getcwd())))
            return
        self.included.append(target)
        self.read_text(text, target)

    # -- queries ----------------------------------------------------------
    def has(self, key):
        return key in self.table

    def count(self, key):
        return sum(1 for e in self.entries if e.key == key)

    def values(self, key):
        return list(self.table[key].values) if key in self.table else []

    def where(self, key):
        return self.table[key].where if key in self.table else ""

    def keys(self, prefix=""):
        """The distinct keys that start with `prefix`, in deck order."""
        seen, out = set(), []
        for e in self.entries:
            if e.key.startswith(prefix) and e.key not in seen:
                seen.add(e.key)
                out.append(e.key)
        return out

    def duplicates(self):
        """Keys defined more than once: {key: [entries]}."""
        groups = {}
        for e in self.entries:
            groups.setdefault(e.key, []).append(e)
        return {k: v for k, v in groups.items() if len(v) > 1}

    # -- typed access (DeckError when the value is not of the type) -------
    def get_str(self, key, default=None):
        if key not in self.table:
            return default
        return self.table[key].values[0]

    def get_reals(self, key, default=None):
        if key not in self.table:
            return default
        out = []
        for v in self.table[key].values:
            try:
                out.append(to_real(v))
            except ValueError:
                raise DeckError(key, "%s = %s: '%s' is not a number" % (key, " ".join(self.table[key].values), v))
            except Unparsed:
                raise DeckError(key, "%s = %s: '%s' is an expression this checker does not evaluate"
                                % (key, " ".join(self.table[key].values), v))
        return out

    def get_real(self, key, default=None):
        vals = self.get_reals(key)
        return default if vals is None else vals[0]

    def get_ints(self, key, default=None):
        if key not in self.table:
            return default
        out = []
        for v in self.table[key].values:
            try:
                out.append(to_int(v))
            except ValueError:
                raise DeckError(key, "%s = %s: '%s' is not an integer" % (key, " ".join(self.table[key].values), v))
            except Unparsed:
                raise DeckError(key, "%s = %s: '%s' is an expression this checker does not evaluate"
                                % (key, " ".join(self.table[key].values), v))
        return out

    def get_int(self, key, default=None):
        vals = self.get_ints(key)
        return default if vals is None else vals[0]


def _comment_start(line):
    quoted = False
    for i, c in enumerate(line):
        if c == '"':
            quoted = not quoted
        elif c == "#" and not quoted:
            return i
    return -1


def parse_deck(path, overrides=(), cwd=None):
    """Read the deck at `path`, then the `key=value` strings of `overrides` as command-line definitions.

    Raises OSError when the file cannot be read; what stops AMReX inside the
    file is kept in `Deck.problems`.
    """
    deck = Deck(path)
    deck.cwd = cwd
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        deck.read_text(fh.read(), path)
    deck.read_command_line(overrides)
    return deck


def parse_text(text, path="<text>", overrides=(), cwd=None):
    """As parse_deck, from a string (used by the self-test and by callers that build a deck in memory)."""
    deck = Deck(path)
    deck.cwd = cwd
    deck.read_text(text, path)
    deck.read_command_line(overrides)
    return deck
