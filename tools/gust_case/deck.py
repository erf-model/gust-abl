"""A GUST deck: ordered ``namespace.key = value`` lines under section headers.

GUST refuses a key given twice, so :meth:`Deck.set` does too; use
:meth:`Deck.replace` or :meth:`Deck.remove` to change one on purpose.
"""
import math


def format_value(value):
    """The text of a deck value: 0/1 for a flag, %.10g numbers, lists space-separated."""
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("a deck value must be finite")
        text = "%.10g" % value
        if "e" not in text and "." not in text:
            text += "."
        return text
    if isinstance(value, (list, tuple)):
        return " ".join(format_value(v) for v in value)
    return str(value)


class Deck:
    def __init__(self, title_lines=()):
        self.title_lines = list(title_lines)
        self._items = []     # ("section", text) or ("key", name)
        self._values = {}
        self._comments = {}

    def section(self, text):
        self._items.append(("section", text))

    def set(self, key, value, comment=None):
        if key in self._values:
            raise KeyError("the deck already has " + key)
        self._items.append(("key", key))
        self._values[key] = format_value(value)
        if comment:
            self._comments[key] = comment

    def replace(self, key, value, comment=None):
        if key not in self._values:
            self.set(key, value, comment)
            return
        self._values[key] = format_value(value)
        if comment:
            self._comments[key] = comment

    def remove(self, key):
        if key in self._values:
            del self._values[key]
            self._items.remove(("key", key))
            self._comments.pop(key, None)

    def __contains__(self, key):
        return key in self._values

    def get(self, key, default=None):
        return self._values.get(key, default)

    def keys(self):
        return [name for kind, name in self._items if kind == "key"]

    def text(self):
        lines = ["# " + line if line else "#" for line in self.title_lines]
        width = max([len(k) for k in self._values] + [1])
        pending = None
        for kind, name in self._items:
            if kind == "section":
                pending = name
                continue
            if pending is not None:
                lines += ["", "# " + pending]
                pending = None
            line = "%-*s = %s" % (width, name, self._values[name])
            if name in self._comments:
                line += "   # " + self._comments[name]
            lines.append(line)
        return "\n".join(lines) + "\n"

    def write(self, path):
        with open(path, "w", encoding="utf-8") as out:
            out.write(self.text())


def read_deck(path):
    """A deck file as {key: [tokens]} (comments dropped). A key given twice raises ValueError."""
    values = {}
    with open(path, "r", encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.split("#", 1)[0].strip()
            if not line:
                continue
            if "=" not in line:
                raise ValueError("%s:%d: not a key = value line" % (path, number))
            key, value = line.split("=", 1)
            key = key.strip()
            if key in values:
                raise ValueError("%s:%d: %s is given twice" % (path, number, key))
            values[key] = value.split()
    return values
