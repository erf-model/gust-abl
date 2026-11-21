"""The main window: the case form on the left; the map, the findings, the summary, the run and the wind library on
the right."""
import datetime
import os

import yaml
from PySide6 import QtCore, QtGui, QtWidgets

from gust_case import CaseError, build, load_yaml, validate
from gust_case import dem

from . import values
from .form import CaseForm
from .librarytab import LibraryTab
from .mapview import PlanView
from .runner import RunPanel
from .runplan import find_gust

HEADER = "# A GUST case description (tools/gust_gui.py; tools/make_gust_case.py writes the case from it).\n"
FINDING_COLUMNS = ("Severity", "Where", "Key", "Message", "Advice")


def starter_case():
    """The starter case of make_gust_case.py --example."""
    import make_gust_case
    return yaml.safe_load(make_gust_case.EXAMPLE)


def save_yaml(path, raw):
    with open(path, "w", encoding="utf-8") as out:
        out.write(HEADER)
        yaml.safe_dump(raw, out, sort_keys=False, default_flow_style=None, allow_unicode=True)


class DownloadWorker(QtCore.QObject):
    """dem.download in a thread, so the window stays live while tiles arrive."""

    progress = QtCore.Signal(int, int, str)
    done = QtCore.Signal(object, object)           # Download or None, the CaseError or None

    def __init__(self, box, out, cache, sea_floor):
        super().__init__()
        self.args = (box, out, cache, sea_floor)

    def run(self):
        box, out, cache, sea_floor = self.args
        try:
            result = dem.download(box, out, cache, sea_floor=sea_floor, date=datetime.date.today().isoformat(),
                                  progress=lambda i, n, name: self.progress.emit(i, n, name))
        except CaseError as err:
            self.done.emit(None, err)
            return
        except Exception as err:                  # noqa: BLE001  (the thread must end with done in every case)
            self.done.emit(None, CaseError("terrain download: %s: %s" % (type(err).__name__, err)))
            return
        self.done.emit(result, None)


class DownloadDialog(QtWidgets.QDialog):
    """Where to write the terrain grid and what to do with heights below 0 m."""

    def __init__(self, parent, out, box):
        super().__init__(parent)
        self.setWindowTitle("Download terrain")
        layout = QtWidgets.QFormLayout(self)
        layout.addRow(QtWidgets.QLabel(
            "Tiles of the AWS Terrain Tiles (SRTM and other sources, about 30 m) for latitude %.4f .. %.4f and "
            "longitude %.4f .. %.4f: the domain and the builder's margin. No account is needed. The attribution "
            "the data need is written beside the grid (.sources.txt)." % box))
        layout.itemAt(0).widget().setWordWrap(True)
        self.out = QtWidgets.QLineEdit(out)
        browse = QtWidgets.QPushButton("...")
        browse.clicked.connect(self._browse)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.out)
        row.addWidget(browse)
        layout.addRow("grid file (.asc)", row)
        self.sea = QtWidgets.QComboBox()
        self.sea.addItems(dem.SEA_FLOOR)
        self.sea.setToolTip("zero: heights below 0 m (the sea floor over the sea) become 0 m.\n"
                            "keep: left as they are (land below sea level).")
        layout.addRow("sea floor", self.sea)
        self.cache = QtWidgets.QLineEdit(dem.default_cache_dir())
        layout.addRow("tile cache", self.cache)
        buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.Ok | QtWidgets.QDialogButtonBox.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def _browse(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Terrain grid", self.out.text(), "ESRI ASCII grid (*.asc)")
        if path:
            self.out.setText(path)


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, path=None, frozen_dir=None):
        super().__init__()
        self.frozen_dir = frozen_dir
        self.path = None
        self.dirty = False
        self.last = None                           # the Result of the last check or write
        self.form = CaseForm(base_dir=self.base_dir)
        self.form.changed.connect(self._edited)
        self.plan_view = PlanView()
        self.findings = QtWidgets.QTableWidget(0, len(FINDING_COLUMNS))
        self.findings.setHorizontalHeaderLabels(FINDING_COLUMNS)
        self.findings.horizontalHeader().setStretchLastSection(True)
        self.findings.setEditTriggers(QtWidgets.QAbstractItemView.NoEditTriggers)
        self.findings.setWordWrap(True)
        self.findings.cellDoubleClicked.connect(self._finding_clicked)
        self.summary = QtWidgets.QPlainTextEdit()
        self.summary.setReadOnly(True)
        self.runner = RunPanel()
        self.runner.finished.connect(self._run_finished)
        self._thread = self._worker = None         # the terrain download, while it runs
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.plan_view, "Map")
        self.tabs.addTab(self.findings, "Findings")
        self.tabs.addTab(self.summary, "Summary")
        self.tabs.addTab(self.runner, "Run")
        self.library = LibraryTab()
        self.tabs.addTab(self.library, "Library")
        split = QtWidgets.QSplitter()
        split.addWidget(self.form)
        split.addWidget(self.tabs)
        split.setSizes([520, 760])
        self.setCentralWidget(split)
        self._actions()
        self.resize(1320, 860)
        if path:
            self.open_path(path)
        else:
            self.form.set_raw(starter_case())
            self.dirty = False
            self._title()

    # ---- menus ----------------------------------------------------------------------------------------------
    def _actions(self):
        bar = self.addToolBar("case")
        bar.setMovable(False)
        menu = self.menuBar().addMenu("&File")

        def add(text, slot, shortcut=None, tip=""):
            action = QtGui.QAction(text, self)
            if shortcut:
                action.setShortcut(shortcut)
            action.setToolTip(tip or text)
            action.triggered.connect(slot)
            bar.addAction(action)
            return action

        for action in (add("New", self.new_case, QtGui.QKeySequence.New, "A starter case (make_gust_case.py --example)"),
                       add("Open...", self.open_dialog, QtGui.QKeySequence.Open),
                       add("Save", self.save, QtGui.QKeySequence.Save),
                       add("Save as...", self.save_as, QtGui.QKeySequence.SaveAs)):
            menu.addAction(action)
        bar.addSeparator()
        case = self.menuBar().addMenu("&Case")
        self.case_actions = (add("Check", self.check, "Ctrl+K", "Check the case and draw it; writes nothing"),
                       add("Download terrain...", self.download_terrain, None,
                           "Fetch elevation tiles for domain.center_lat_lon and set terrain.file"),
                       add("Write case", self.write_case, "Ctrl+W", "Write the case folder and check its decks"),
                       add("Run", self.run_case, "Ctrl+R", "Write the case and run gust, step by step as run.sh"))
        for action in self.case_actions:
            case.addAction(action)
        bar.addSeparator()
        self.advanced_action = QtGui.QAction("Show all settings", self)
        self.advanced_action.setCheckable(True)
        self.advanced_action.setToolTip("Off: only the settings a case needs (the rest keep their defaults); "
                                        "on: every setting")
        self.advanced_action.toggled.connect(self.form.set_advanced)
        bar.addAction(self.advanced_action)
        view = self.menuBar().addMenu("&View")
        view.addAction(self.advanced_action)

    def _busy(self, busy):
        """Check, Download, Write and Run wait while a run or a download goes on (they write the case folder)."""
        for action in self.case_actions:
            action.setEnabled(not busy)

    def _title(self):
        name = os.path.basename(self.path) if self.path else "new case"
        self.setWindowTitle("GUST case set-up - %s%s" % (name, " *" if self.dirty else ""))

    def _edited(self):
        self.dirty = True
        self._title()

    def _ask_discard(self):
        if not self.dirty:
            return True
        answer = QtWidgets.QMessageBox.question(self, "Unsaved changes", "Discard the changes to the case?")
        return answer == QtWidgets.QMessageBox.Yes

    def closeEvent(self, event):
        if self.runner.running():
            QtWidgets.QMessageBox.information(self, "Running", "Stop the run first.")
            event.ignore()
            return
        if self._thread is not None:
            QtWidgets.QMessageBox.information(self, "Downloading", "Wait for the terrain download to end.")
            event.ignore()
            return
        if self._ask_discard():
            event.accept()
        else:
            event.ignore()

    # ---- files ----------------------------------------------------------------------------------------------
    def new_case(self):
        if not self._ask_discard():
            return
        self.path = None
        self.form.set_raw(starter_case())
        self.dirty = False
        self._title()

    def open_dialog(self):
        if not self._ask_discard():
            return
        path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Open a case file", "", "YAML (*.yaml *.yml);;All (*)")
        if path:
            self.open_path(path)

    def open_path(self, path):
        try:
            raw = load_yaml(path)
        except CaseError as err:
            self._show_refusals(err.problems)
            return False
        self.path = os.path.abspath(path)
        self.form.set_raw(raw)
        self.dirty = False
        self._title()
        self.plan_view.clear()
        self._fill_findings([])
        return True

    def save(self):
        if self.path is None:
            return self.save_as()
        return self._save_to(self.path)

    def save_as(self):
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save the case file", self.path or "case.yaml",
                                                        "YAML (*.yaml *.yml)")
        if not path:
            return False
        return self._save_to(path)

    def _save_to(self, path):
        try:
            raw = self.form.raw()
        except CaseError as err:
            self._show_refusals(err.problems)
            return False
        save_yaml(path, raw)
        self.path = os.path.abspath(path)
        self.dirty = False
        self._title()
        return True

    def base_dir(self):
        return os.path.dirname(self.path) if self.path else os.getcwd()

    # ---- check and write ------------------------------------------------------------------------------------
    def _show_refusals(self, problems):
        marked = self.form.mark(problems)
        rows = [("refused", "case file", name or "", problem, "") for problem, name in marked]
        self._fill_findings(rows)
        self.tabs.setCurrentWidget(self.findings)
        if marked and marked[0][1]:
            self.form.show_field(marked[0][1])
        self.statusBar().showMessage("The case was refused: %d problem%s" % (len(problems), "" if len(problems) == 1 else "s"))

    def _fill_findings(self, rows):
        self.findings.setRowCount(len(rows))
        for i, row in enumerate(rows):
            for j, text in enumerate(row):
                item = QtWidgets.QTableWidgetItem(str(text))
                if j == 0 and text in ("refused", "ERROR"):
                    item.setForeground(QtGui.QColor(192, 57, 43))
                self.findings.setItem(i, j, item)
        self.findings.resizeRowsToContents()

    def _finding_clicked(self, row, column):
        item = self.findings.item(row, 2)
        if item is not None and self.findings.item(row, 0).text() == "refused":
            self.form.show_field(item.text())

    def _build(self, dry_run, overwrite=False):
        """validate + build; marks the refused fields. Returns the Result or None."""
        try:
            raw = self.form.raw()
            case, given = validate(raw)
            result = build(case, given, self.base_dir(), dry_run=dry_run, overwrite=overwrite, check=not dry_run)
        except CaseError as err:
            self._show_refusals(err.problems)
            return None
        except OSError as err:                    # a file in use (Windows) or a folder that is not writable
            self._show_refusals(["cannot write the case: %s" % err])
            return None
        self.form.mark([])
        self.last = result
        self.summary.setPlainText(result.summary)
        self.plan_view.show_plan(result.plan, result.terrain)
        rows = []
        for rel, report in result.reports:
            for f in report.sorted():
                rows.append((f.severity.name, rel, f.key, f.message, f.advice))
        self._fill_findings(rows)
        return result

    def check(self):
        result = self._build(dry_run=True)
        if result is not None:
            self.tabs.setCurrentWidget(self.plan_view)
            self.statusBar().showMessage("The case is accepted; Write case writes %s" % result.case_dir)
        return result

    def _ready_to_write(self):
        """Save the case file (the case folder's relative paths start from it)."""
        if self.path is None or self.dirty:
            if self.path is None:
                QtWidgets.QMessageBox.information(self, "Save first", "Save the case file first: the case folder "
                                                  "and the files it names are relative to it.")
            return self.save()
        return True

    def write_case(self):
        if not self._ready_to_write():
            return None
        probe = self._build(dry_run=True)
        if probe is None:
            return None
        overwrite = False
        if os.path.exists(probe.case_dir):
            answer = QtWidgets.QMessageBox.question(
                self, "Overwrite", "%s exists. Replace the files make_gust_case.py wrote there? Run outputs "
                "(plotfiles, checkpoints, logs) stay." % probe.case_dir)
            if answer != QtWidgets.QMessageBox.Yes:
                return None
            overwrite = True
        result = self._build(dry_run=False, overwrite=overwrite)
        if result is not None:
            self.library.set_folder(result.case_dir)
            self.tabs.setCurrentWidget(self.findings if result.errors or result.warnings else self.summary)
            self.statusBar().showMessage("Wrote %s: the deck checker found %d error%s and %d warning%s"
                                         % (result.case_dir, result.errors, "" if result.errors == 1 else "s",
                                            result.warnings, "" if result.warnings == 1 else "s"))
        return result

    # ---- terrain download -----------------------------------------------------------------------------------
    def download_terrain(self):
        try:
            raw = self.form.raw()
            case, given = validate(raw)
            box = dem.box_for_case(case, given)
        except CaseError as err:
            self._show_refusals(err.problems)
            return
        current = case["terrain"]["file"]
        out = current if current and current.lower().endswith(".asc") else "dem.asc"
        if not os.path.isabs(out):
            out = os.path.join(self.base_dir(), out)
        dialog = DownloadDialog(self, out, box)
        if dialog.exec() != QtWidgets.QDialog.Accepted:
            return
        self._thread = QtCore.QThread(self)
        self._worker = DownloadWorker(box, dialog.out.text(), dialog.cache.text(), dialog.sea.currentText())
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.progress.connect(self._download_progress)    # a slot of the window: runs in its thread
        self._worker.done.connect(self._downloaded)
        self._worker.done.connect(self._thread.quit)
        self._thread.finished.connect(self._download_ended)
        self._busy(True)
        self.statusBar().showMessage("Downloading terrain tiles...")
        self._thread.start()

    def _download_progress(self, i, n, name):
        self.statusBar().showMessage("Downloading tile %d of %d: %s" % (i, n, name))

    def _download_ended(self):
        self._worker.deleteLater()
        self._thread.deleteLater()
        self._thread = self._worker = None
        self._busy(self.runner.running())

    def _downloaded(self, result, error):
        if error is not None:
            self._fill_findings([("refused", "terrain download", "", p, "") for p in error.problems])
            self.tabs.setCurrentWidget(self.findings)
            self.statusBar().showMessage("The terrain download failed")
            return
        path = result.path
        if self.path is not None:
            try:
                rel = os.path.relpath(path, self.base_dir())
            except ValueError:                    # another drive (Windows)
                rel = path
            path = rel if not rel.startswith("..") else path
        self.form.fields["terrain.file"].set_text(path)
        notes = list(result.notes)
        for key in ("terrain.format", "terrain.coordinates"):
            field = self.form.fields[key]
            if field.text().strip() not in ("", "auto"):
                notes.append("%s was %s; cleared (auto reads the downloaded grid)" % (key, field.text().strip()))
                field.set_text("")
        heights = ("%.6g to %.6g m" % (result.lowest, result.highest) if result.lowest is not None
                   else "no heights (every point a void)")
        lines = ["Wrote %s: %d x %d points, %s" % (result.path, result.shape[1], result.shape[0], heights)]
        lines += ["note: " + n for n in notes]
        with open(result.sources, encoding="utf-8") as handle:
            lines += ["", handle.read()]
        self.summary.setPlainText("\n".join(lines))
        self.tabs.setCurrentWidget(self.summary)
        self.statusBar().showMessage("terrain.file = %s; the attribution is in %s" % (path, result.sources))

    # ---- run ------------------------------------------------------------------------------------------------
    def run_case(self):
        if self.runner.running():
            self.statusBar().showMessage("A run is going on (or stopping)")
            return
        result = self.write_case()
        if result is None:
            return
        r = result.plan.case["run"]
        given = values.get_path(self.form.raw(), "run.executable") is not None
        gust = find_gust(r["executable"], result.case_dir, self.frozen_dir, given=given)
        if gust is None:
            QtWidgets.QMessageBox.warning(self, "gust not found", "No gust executable: run.executable = %s is not a "
                                          "file or on PATH%s." % (r["executable"], " (nor beside this program)"
                                                                  if self.frozen_dir else ""))
            return
        if result.errors:
            answer = QtWidgets.QMessageBox.question(self, "Deck checker", "The deck checker found %d error%s. Run "
                                                    "anyway?" % (result.errors, "" if result.errors == 1 else "s"))
            if answer != QtWidgets.QMessageBox.Yes:
                return
        self.tabs.setCurrentWidget(self.runner)
        self._busy(True)
        self.runner.start(result.case_dir, result.runs, gust, r["ranks"], r["mpi_launcher"], r["resume"])

    def _run_finished(self, ok, message):
        self._busy(self._thread is not None)
        self.statusBar().showMessage(message)
