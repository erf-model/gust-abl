"""The case form: one field per key of gust_case.schema.SPEC, one tab per section.

Fields carry plain-language labels and sit in groups (gust_window.layout); each field's tooltip names its YAML key
and holds the key's text, unit and accepted values (the --list-keys table). An empty field is a key not given: what
to type and the default show in grey. A refused key (CaseError) marks its field red and puts the refusal in the
tooltip. The basic mode shows only layout.BASIC, and any other field that holds a value or a refusal.
"""
import copy
import os
import re

from PySide6 import QtCore, QtWidgets

from gust_case import CaseError
from gust_case.schema import SECTION_DOC, SPEC, NestList

from . import layout, values

MARK_STYLE = "border: 2px solid #c0392b;"
_KINDS = (("a number", "number"), ("an integer", "whole number"), ("a list of numbers", "list of numbers"),
          ("a list of words", "list of words"), ("a file path", "file path"), ("see the text", "YAML (see the tooltip)"))


def format_hint(name, key):
    """The grey text of an empty field: what to type (kind and range; the label holds the unit) and the default,
    if any."""
    short = name.rsplit(".", 1)[-1]
    if key.kind == "choice":
        text = "one of: " + ", ".join(key.choices)
    elif name in layout.SHAPES_FULL:
        text = layout.SHAPES_FULL[name]
    elif short in layout.SHAPES:
        text = layout.SHAPES[short]
    else:
        text = key.range_text()
        for long, kind in _KINDS:
            text = text.replace(long, kind)
        if key.kind == "numbers" and key.length:
            text = text.replace("%d numbers" % key.length, ", ".join(["number"] * key.length))
        if key.lo is not None and 0.0 < key.lo < 1e-10:          # tke_floor: "≥ 2.22e-16" reads as "> 0"
            text = text.replace("≥ %g" % key.lo, "> 0")
        text = re.sub(r"≥ (\S+) and ≤ (\S+)", r"\1 to \2", text)
    if key.default is not None and key.default != "":
        text += "; default " + values.to_text(key, key.default).replace("\n", "; ")
    return text


def _tooltip(name, key):
    unit = " (%s)" % key.unit if key.unit else ""
    return "<b>%s</b>%s<br>%s<br><i>%s</i>" % (name, unit, key.doc, key.range_text())


class Field:
    """One key's input, its text and its mark. ``widget`` is the input (a line, a list, a text box, or for a pair
    the box holding its two lines); ``row`` is what the form lays out (the input, with a Browse button for a path)."""

    def __init__(self, name, key, parent=None, base_dir=None):
        self.name = name
        self.key = key
        self.edits = None
        short = name.rsplit(".", 1)[-1]
        if key.kind in ("flag", "choice"):
            w = QtWidgets.QComboBox(parent)
            w.setEditable(True)                    # a value the list does not hold stays, and is refused by name
            w.addItems([""] + (["true", "false"] if key.kind == "flag" else list(key.choices)))
            w.lineEdit().setPlaceholderText(format_hint(name, key))
            self.changed = w.currentTextChanged
        elif key.kind == "any":
            w = QtWidgets.QPlainTextEdit(parent)
            w.setPlaceholderText("namespace.key: value, one per line" if name == "gust_inputs" else
                                 format_hint(name, key))
            w.setMinimumHeight(120)
            self.changed = w.textChanged
        elif short in layout.PAIRS:
            w = QtWidgets.QWidget(parent)
            box = QtWidgets.QHBoxLayout(w)
            box.setContentsMargins(0, 0, 0, 0)
            self.edits = []
            for hint in layout.PAIRS[short]:
                edit = QtWidgets.QLineEdit(w)
                edit.setPlaceholderText(hint)
                box.addWidget(edit)
                self.edits.append(edit)
            self.changed = _Joined(self.edits)
        else:
            w = QtWidgets.QLineEdit(parent)
            w.setPlaceholderText(format_hint(name, key))
            self.changed = w.textChanged
        w.setObjectName(name)
        w.setProperty("gust_key", name)
        self.widget = w
        self.row = w
        if key.kind == "path":
            self.row = QtWidgets.QWidget(parent)
            box = QtWidgets.QHBoxLayout(self.row)
            box.setContentsMargins(0, 0, 0, 0)
            box.addWidget(w)
            button = QtWidgets.QPushButton("Browse...", self.row)
            button.clicked.connect(lambda: self._browse(base_dir() if base_dir else "."))
            box.addWidget(button)
            self.browse_button = button
        self.tooltip = _tooltip(name, key)
        w.setToolTip(self.tooltip)
        self.marked = None
        self._loaded = None                       # (text, value) set from the file: kept while the text is

    def _browse(self, base):
        """A folder or file dialog; the choice is written relative to the case file's folder when inside it."""
        start = os.path.join(base, self.text()) if self.text() else base
        if self.name in layout.FOLDERS:
            path = QtWidgets.QFileDialog.getExistingDirectory(self.row, "Choose the folder", start)
        else:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self.row, "Choose the file", start)
        if path:
            self.set_text(_relative(path, base))

    def text(self):
        w = self.widget
        if self.edits is not None:
            parts = [e.text().strip() for e in self.edits]
            return "" if not any(parts) else ", ".join(parts)
        if isinstance(w, QtWidgets.QComboBox):
            return w.currentText()
        if isinstance(w, QtWidgets.QPlainTextEdit):
            return w.toPlainText()
        return w.text()

    def set_text(self, text):
        w = self.widget
        if self.edits is not None:
            parts = [t.strip() for t in text.strip("[] ").split(",")] if text.strip() else []
            if len(parts) != len(self.edits):
                parts = [text] + [""] * (len(self.edits) - 1)        # not a pair: kept whole, refused by name
            for edit, part in zip(self.edits, parts):
                edit.setText(part)
        elif isinstance(w, QtWidgets.QComboBox):
            w.setCurrentText(text)
        elif isinstance(w, QtWidgets.QPlainTextEdit):
            w.setPlainText(text)
        else:
            w.setText(text)

    def set_value(self, value):
        text = values.to_text(self.key, value)
        self.set_text(text)
        self._loaded = (self.text(), value)

    def value(self):
        """(given, value). A field not edited since the file was read gives the file's value as it was (so a
        value the builder refuses, such as name: yes, is not turned into one it accepts by the round trip)."""
        if self._loaded is not None and self.text() == self._loaded[0]:
            return self._loaded[1] is not None, copy.deepcopy(self._loaded[1])
        return values.from_text(self.name, self.key, self.text())

    def _style(self, style):
        for w in self.edits if self.edits is not None else [self.widget]:
            w.setStyleSheet(style)

    def mark(self, problem):
        self.marked = problem
        self._style(MARK_STYLE)
        self.widget.setToolTip("<span style='color:#c0392b'>%s</span><br><br>%s" % (problem, self.tooltip))

    def unmark(self):
        self.marked = None
        self._style("")
        self.widget.setToolTip(self.tooltip)

    def focus(self):
        (self.edits[0] if self.edits is not None else self.widget).setFocus()


class _Joined:
    """The ``changed`` signal of a pair: connect() connects every line of it."""

    def __init__(self, edits):
        self.edits = edits

    def connect(self, slot):
        for edit in self.edits:
            edit.textChanged.connect(slot)


def _relative(path, base):
    try:
        rel = os.path.relpath(path, base)
    except ValueError:                            # another drive on Windows
        return path
    return path if rel.startswith("..") else rel.replace(os.sep, "/")


class _Rows:
    """The rows of one tab or patch, for the basic mode: (form layout, row index, field) of each field, and the
    groups with the fields in them."""

    def __init__(self):
        self.rows = []
        self.groups = []


def _section_form(parent, spec, prefix, fields, on_change, rows, base_dir, section=None):
    """The fields of ``spec`` laid out in the tab's groups (layout.GROUPS; subsections as headed groups) in a
    QVBoxLayout; fills ``fields`` and ``rows``."""
    outer = QtWidgets.QVBoxLayout()
    keys = [k for k, e in spec.items() if not isinstance(e, NestList)]
    for title, note, names in layout.groups(section, keys):
        box = QtWidgets.QGroupBox(title, parent) if title else QtWidgets.QWidget(parent)
        form = QtWidgets.QFormLayout(box)
        form.setFieldGrowthPolicy(QtWidgets.QFormLayout.AllNonFixedFieldsGrow)
        if note:
            text = QtWidgets.QLabel(note, box)
            text.setWordWrap(True)
            text.setStyleSheet("color: #555;")
            form.addRow(text)
        members = []
        for name in names:
            entry, full = spec[name], prefix + name
            if isinstance(entry, dict):
                sub = QtWidgets.QGroupBox(full, box)
                sub.setLayout(_section_form(sub, entry, full + ".", fields, on_change, rows, base_dir))
                form.addRow(sub)
                members.append((full + ".", sub))
                continue
            field = Field(full, entry, box, base_dir)
            field.changed.connect(on_change)
            fields[full] = field
            label = QtWidgets.QLabel(layout.label(full, entry), box)
            label.setToolTip(field.tooltip)
            form.addRow(label, field.row)
            rows.rows.append((form, form.rowCount() - 1, field))
            members.append(field)
        rows.groups.append((box, members))
        outer.addWidget(box)
    outer.addStretch(1)
    return outer


class NestBox(QtWidgets.QGroupBox):
    """The fields of one nested patch."""

    def __init__(self, form, index):
        super().__init__(form)
        self.form = form
        self.fields = {}
        self.extra = {}
        self.rows = _Rows()
        layout_ = QtWidgets.QVBoxLayout(self)
        inner = QtWidgets.QWidget(self)
        inner.setLayout(_section_form(inner, values.nest_spec(), "", self.fields, form._changed, self.rows,
                                      form.base_dir))
        layout_.addWidget(inner)
        remove = QtWidgets.QPushButton("Remove this patch", self)
        remove.clicked.connect(lambda: form.remove_nest(self))
        layout_.addWidget(remove, alignment=QtCore.Qt.AlignRight)
        self.set_index(index)

    def set_index(self, index):
        self.index = index
        self.setTitle("Patch %d (nests[%d])" % (index + 1, index))
        for name, field in self.fields.items():
            full = "nests[%d].%s" % (index, name)
            field.name = full
            field.tooltip = _tooltip(full, field.key)
            field.widget.setObjectName(full)
            if field.marked is None:
                field.widget.setToolTip(field.tooltip)


class CaseForm(QtWidgets.QTabWidget):
    """The whole case description. ``set_raw``/``raw`` read and write the nested dicts of the YAML file."""

    changed = QtCore.Signal()

    def __init__(self, parent=None, base_dir=None):
        super().__init__(parent)
        self.fields = {}
        self.nests = []
        self.extra = {}
        self._quiet = False
        self.base_dir = base_dir or os.getcwd     # the case file's folder: relative paths start there
        self.advanced = False
        self.rows = _Rows()
        self.pages = {}
        for section, entry in SPEC.items():
            page = QtWidgets.QWidget()
            outer = QtWidgets.QVBoxLayout(page)
            doc = QtWidgets.QLabel(SECTION_DOC.get(section, ""))
            doc.setWordWrap(True)
            outer.addWidget(doc)
            if isinstance(entry, dict):
                body = QtWidgets.QWidget()
                body.setLayout(_section_form(body, entry, section + ".", self.fields, self._changed, self.rows,
                                             self._base, section))
            elif isinstance(entry, NestList):
                body = self._nest_page()
            else:                                  # gust_inputs: one key, the section itself
                body = QtWidgets.QWidget()
                body.setLayout(_section_form(body, {section: entry}, "", self.fields, self._changed, self.rows,
                                             self._base))
            scroll = QtWidgets.QScrollArea()
            scroll.setWidgetResizable(True)
            scroll.setWidget(body)
            outer.addWidget(scroll)
            self.addTab(page, section)
            self.pages[section] = page
        self.set_advanced(False)

    def _base(self):
        return self.base_dir()

    # ---- basic and advanced -------------------------------------------------------------------------------------
    def _shown(self, field):
        if self.advanced or field.marked is not None or field.text().strip():
            return True
        return field.name in layout.BASIC

    def set_advanced(self, advanced):
        """All fields (True), or the basic ones and every field that holds a value or a refusal (False)."""
        self.advanced = bool(advanced)
        self.refresh_visible()

    def refresh_visible(self):
        boxes = [self.rows] + [b.rows for b in self.nests]
        shown = {}
        for rows in boxes:
            for form, index, field in rows.rows:
                visible = self._shown(field) or rows is not self.rows
                form.setRowVisible(index, visible)
                shown[field.name] = visible
            for box, members in rows.groups:              # subsections come before the group holding them
                if rows is not self.rows:
                    continue
                any_shown = False
                for m in members:
                    if isinstance(m, Field):
                        any_shown |= shown[m.name]
                    else:
                        prefix, sub = m
                        sub_shown = any(v for n, v in shown.items() if n.startswith(prefix))
                        sub.setVisible(sub_shown)
                        any_shown |= sub_shown
                box.setVisible(any_shown)
        for i in range(self.count()):
            section = self.tabText(i)
            if section == "nests":
                visible = self.advanced or bool(self.nests)
            else:
                visible = any(v for n, v in shown.items() if n.split(".")[0] == section)
            self.setTabVisible(i, visible)

    # ---- nests --------------------------------------------------------------------------------------------
    def _nest_page(self):
        body = QtWidgets.QWidget()
        self._nest_layout = QtWidgets.QVBoxLayout(body)
        add = QtWidgets.QPushButton("Add a patch", body)
        add.clicked.connect(self.add_nest)
        self._nest_layout.addWidget(add, alignment=QtCore.Qt.AlignLeft)
        self._nest_layout.addStretch(1)
        return body

    def add_nest(self):
        box = NestBox(self, len(self.nests))
        self._nest_layout.insertWidget(self._nest_layout.count() - 1, box)
        self.nests.append(box)
        self.refresh_visible()
        self._changed()
        return box

    def remove_nest(self, box):
        self.nests.remove(box)
        box.setParent(None)
        box.deleteLater()
        for i, b in enumerate(self.nests):
            b.set_index(i)
        self._changed()

    # ---- all fields -----------------------------------------------------------------------------------------
    def all_fields(self):
        """name -> Field of every field, the patches' as nests[i].key."""
        out = dict(self.fields)
        for box in self.nests:
            out.update({f.name: f for f in box.fields.values()})
        return out

    def _changed(self, *args):
        if not self._quiet:
            self.changed.emit()

    def set_raw(self, raw):
        """Show the YAML file's dicts; keys the schema does not know are kept for ``raw``."""
        self._quiet = True
        try:
            for name, field in self.fields.items():
                field.set_value(values.get_path(raw, name))
                field.unmark()
            for box in list(self.nests):
                self.remove_nest(box)
            nests = raw.get("nests") if isinstance(raw, dict) else None
            self.extra = values.unknown_parts(raw)
            if isinstance(nests, list) and all(isinstance(n, dict) for n in nests):
                for item in nests:
                    box = self.add_nest()
                    for name, field in box.fields.items():
                        field.set_value(item.get(name))
                    box.extra = {k: v for k, v in item.items() if k not in box.fields}
            elif nests is not None:
                self.extra["nests"] = nests        # not a list of patches: kept, and refused by name
        finally:
            self._quiet = False
        self.refresh_visible()
        self.changed.emit()

    def raw(self):
        """The nested dicts the form holds (only the keys given), in the schema's order; a field whose text is not
        YAML raises CaseError naming it."""
        out, problems = {}, []
        for name, _ in values.field_names():
            field = self.fields[name]
            try:
                given, value = field.value()
            except ValueError as err:
                problems.append(str(err))
                continue
            if given:
                values.set_path(out, name, value)
        nests = []
        for box in self.nests:
            item = {}
            for name, field in box.fields.items():
                try:
                    given, value = field.value()
                except ValueError as err:
                    problems.append(str(err))
                    continue
                if given:
                    item[name] = value
            item.update(box.extra)
            nests.append(item)
        if nests:
            out["nests"] = nests
        if problems:
            raise CaseError(problems)
        ordered = {name: out[name] for name in SPEC if name in out}
        return values.merge(ordered, copy.deepcopy(self.extra))

    def mark(self, problems):
        """Mark the fields the refusals name; returns [(problem, field name or None)]."""
        fields = self.all_fields()
        for field in fields.values():
            field.unmark()
        out = []
        for problem in problems:
            name = values.problem_field(problem, fields)
            if name is not None:
                fields[name].mark(problem if fields[name].marked is None else fields[name].marked + "<br>" + problem)
            out.append((problem, name))
        self.refresh_visible()
        return out

    def show_field(self, name):
        """Bring the field's tab to the front and give it the focus."""
        field = self.all_fields().get(name)
        if field is None:
            return
        section = name.split(".")[0].split("[")[0]
        for i in range(self.count()):
            if self.tabText(i) == section:
                self.setCurrentIndex(i)
        if not self._shown(field):
            self.set_advanced(True)
        field.focus()
