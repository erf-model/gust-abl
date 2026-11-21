"""The Library tab: the wind of a wind library for one direction and speed, or statistics over a distribution of
directions and speeds, drawn as a map and exported as WindNinja-format grids or a CSV at points.

The tab calls tools/gust_library (Library, ensemble, grids), as tools/make_gust_library.py query and ensemble do;
docs/case_setup.rst describes the options. x is east and y north in the case frame, as on the Map tab.
"""
import os

import numpy as np
import pyqtgraph as pg
from PySide6 import QtWidgets

from gust_library import Library, LibraryError, SampleError, ensemble, grids

ONE, DISTRIBUTION = "One wind", "Distribution"
ONE_WIND = ("speed", "u", "v")              # the fields export_grids writes as _vel and _ang
DISTRIBUTION_WIND = ("mean_speed", "mean_u", "mean_v")
EXAMPLE_MEMBERS = "# direction (degrees from)  speed (m/s at the reference point)  weight\n"


class LibraryTab(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.lib = None
        self.fields = None          # the last fields drawn: dict of (nx, ny) arrays
        self.wind = None            # the names of the speed and wind-component fields among them
        self.used = None            # the options the fields were computed from (exports name and sample these)
        layout = QtWidgets.QVBoxLayout(self)

        row = QtWidgets.QHBoxLayout()
        self.folder = QtWidgets.QLineEdit()
        self.folder.setPlaceholderText("case folder (or its wind_library folder)")
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self._browse)
        self.load_button = QtWidgets.QPushButton("Load")
        self.load_button.clicked.connect(lambda: self.load(self.folder.text()))
        row.addWidget(QtWidgets.QLabel("Library"))
        row.addWidget(self.folder, 1)
        row.addWidget(browse)
        row.addWidget(self.load_button)
        layout.addLayout(row)

        form = QtWidgets.QFormLayout()
        self.mode = QtWidgets.QComboBox()
        self.mode.addItems([ONE, DISTRIBUTION])
        self.mode.currentTextChanged.connect(self._mode_changed)
        form.addRow("Wind", self.mode)
        self.height = QtWidgets.QComboBox()
        form.addRow("Height (m above the ground)", self.height)
        self.direction = QtWidgets.QDoubleSpinBox()
        self.direction.setRange(0.0, 359.999)
        self.direction.setDecimals(3)
        self.direction.setWrapping(True)
        form.addRow("Direction (degrees the wind comes from)", self.direction)
        self.speed = QtWidgets.QDoubleSpinBox()
        self.speed.setRange(0.01, 200.0)
        self.speed.setDecimals(3)
        self.speed.setValue(8.0)
        form.addRow("Speed (m/s at the reference point)", self.speed)
        self.members = QtWidgets.QPlainTextEdit(EXAMPLE_MEMBERS)
        self.members.setMaximumHeight(110)
        form.addRow("Members", self.members)
        self.percentiles = QtWidgets.QLineEdit("10 50 90")
        form.addRow("Percentiles of the speed", self.percentiles)
        self.thresholds = QtWidgets.QLineEdit("")
        self.thresholds.setPlaceholderText("m/s, e.g. 10 15")
        form.addRow("Exceedance thresholds (m/s)", self.thresholds)
        layout.addLayout(form)
        self._form = form

        row = QtWidgets.QHBoxLayout()
        self.show_button = QtWidgets.QPushButton("Show")
        self.show_button.clicked.connect(self.compute)
        self.field = QtWidgets.QComboBox()
        self.field.currentTextChanged.connect(self._draw)
        self.units = QtWidgets.QComboBox()
        self.units.addItems(list(grids.SPEED_UNITS))                # mps first: the default, as on the command line
        self.units.setCurrentText("mps")
        self.grids_button = QtWidgets.QPushButton("Export grids…")
        self.grids_button.clicked.connect(self._export_grids_dialog)
        self.points_button = QtWidgets.QPushButton("Export points…")
        self.points_button.clicked.connect(self._export_points_dialog)
        row.addWidget(self.show_button)
        row.addWidget(QtWidgets.QLabel("Field"))
        row.addWidget(self.field, 1)
        row.addWidget(QtWidgets.QLabel("Grid speed units"))
        row.addWidget(self.units)
        row.addWidget(self.grids_button)
        row.addWidget(self.points_button)
        layout.addLayout(row)

        self.caption = QtWidgets.QLabel("Load a library built by make_gust_library.py build.")
        self.caption.setWordWrap(True)
        layout.addWidget(self.caption)
        self.graphics = pg.GraphicsLayoutWidget()
        self.plot = self.graphics.addPlot(row=0, col=0)
        self.plot.setAspectLocked(True)
        self.plot.setLabel("bottom", "x (m, east)")
        self.plot.setLabel("left", "y (m, north)")
        self.color_bar = None
        layout.addWidget(self.graphics, 1)
        self._mode_changed(ONE)
        self._enable()

    # ---- state ------------------------------------------------------------------------------------------------

    def _enable(self):
        have = self.lib is not None
        for w in (self.show_button, self.height):
            w.setEnabled(have)
        for w in (self.grids_button, self.points_button):
            w.setEnabled(self.fields is not None)

    def _mode_changed(self, mode):
        one = mode == ONE
        for w, show in ((self.direction, one), (self.speed, one), (self.members, not one),
                        (self.percentiles, not one), (self.thresholds, not one)):
            w.setVisible(show)
            self._form.labelForField(w).setVisible(show)
        self.fields = None
        self._enable()

    def set_folder(self, path):
        """Fill the folder box (the window does this for the case it writes); loads it when it holds a library."""
        self.folder.setText(path)
        try:
            self.load(path)
        except LibraryError:
            pass

    def load(self, path):
        """Read the library at path (a case folder or its wind_library folder); raises LibraryError."""
        try:
            lib = Library(path)
        except LibraryError as err:
            self.lib = None                                          # no export from the library loaded before
            self.fields = None
            self.height.clear()
            self.plot.clear()
            self._enable()
            self.caption.setText(str(err))
            raise
        self.lib = lib
        self.fields = None
        self.height.clear()
        self.height.addItems(["%g" % h for h in lib.heights])
        dirs = [s["direction"] for s in lib.sectors]
        self.direction.setValue(dirs[0])
        self.caption.setText("%s: %d sectors (%s degrees), heights %s m; reference point (%g, %g) m"
                             % (lib.folder, len(dirs), ", ".join("%g" % d for d in dirs),
                                ", ".join("%g" % h for h in lib.heights), *lib.index["reference_point"]))
        self._enable()
        return lib

    def _members(self):
        return ensemble.parse(self.members.toPlainText().splitlines(), "Members")

    def _numbers(self, box, name):
        try:
            return [float(w) for w in box.text().replace(",", " ").split()]
        except ValueError:
            raise LibraryError("%s: %r is not a list of numbers" % (name, box.text()))

    def compute(self):
        """The fields for the options shown; draws them. Returns the fields, or None after a refusal (shown in the
        caption)."""
        if self.lib is None:
            return None
        h = self.lib.heights[self.height.currentIndex()]
        try:
            used = {"height": h}
            if self.mode.currentText() == ONE:
                used.update(direction=self.direction.value(), speed=self.speed.value())
                w = self.lib.wind(used["direction"], used["speed"], h)
                fields = {"speed": np.hypot(w["u"], w["v"]), "direction": grids.direction_from(w["u"], w["v"]),
                          "u": w["u"], "v": w["v"], "w": w["w"], "k": w["k"]}
                self.wind = ONE_WIND
                text = "blend " + ", ".join("%g degrees x %.4g" % dw for dw in w["weights"])
            else:
                members = self._members()
                used.update(members=members, percentiles=self._numbers(self.percentiles, "Percentiles"),
                            thresholds=self._numbers(self.thresholds, "Exceedance thresholds"))
                fields = ensemble.ensemble(self.lib, members, h, used["percentiles"], used["thresholds"])
                for f in ("x", "y"):
                    fields.pop(f)
                self.wind = DISTRIBUTION_WIND
                text = "%d members, %d directions" % (len(members), len({m[0] for m in members}))
        except (LibraryError, SampleError) as err:
            self.caption.setText(str(err))
            self.fields = None
            self._enable()
            return None
        self.fields = fields
        self.used = used
        self.caption.setText("%g m above the ground: %s" % (h, text))
        current = self.field.currentText()
        self.field.blockSignals(True)
        self.field.clear()
        self.field.addItems(list(fields))
        self.field.setCurrentText(current if current in fields else list(fields)[0])
        self.field.blockSignals(False)
        self._draw(self.field.currentText())
        self._enable()
        return fields

    # ---- drawing ----------------------------------------------------------------------------------------------

    def _draw(self, name):
        if not self.fields or name not in self.fields:
            return
        self.plot.clear()
        if self.color_bar is not None:
            self.graphics.removeItem(self.color_bar)
        image = self.fields[name]
        x, y = self.lib.x, self.lib.y
        dx, dy = x[1] - x[0], y[1] - y[0]
        finite = image[np.isfinite(image)]
        lo, hi = (float(finite.min()), float(finite.max())) if finite.size else (0.0, 1.0)
        hi = hi if hi > lo else lo + 1.0
        item = pg.ImageItem(image, axisOrder="col-major")            # image[i, j]: i east, j north
        cmap = pg.colormap.get("viridis")
        item.setColorMap(cmap)
        item.setLevels((lo, hi))
        item.setRect(float(x[0] - 0.5 * dx), float(y[0] - 0.5 * dy), float(len(x) * dx), float(len(y) * dy))
        self.plot.addItem(item)
        self.color_bar = pg.ColorBarItem(values=(lo, hi), colorMap=cmap, label=name, interactive=False)
        self.graphics.addItem(self.color_bar, row=0, col=1)
        self.plot.autoRange()

    # ---- export -----------------------------------------------------------------------------------------------

    def prefix(self):
        """The file name prefix of the fields shown (from the options of the last Show, not the boxes now)."""
        case = self.lib.index.get("case") or "library"
        g = lambda v: ("%g" % v).replace(".", "p")                   # noqa: E731
        u = self.used
        if self.wind == ONE_WIND:
            return "%s_%sdeg_%smps_%sm_agl" % (case, g(u["direction"]), g(u["speed"]), g(u["height"]))
        return "%s_ensemble_%sm_agl" % (case, g(u["height"]))

    def export_grids(self, out_dir):
        """The fields shown (computed by the last Show) as WindNinja-format grids in out_dir: <prefix>_vel.asc and _ang.asc (and the other
        statistics of a distribution). Returns the files written."""
        if self.wind == ONE_WIND:
            extra, speeds = [], []
        else:
            extra = [n for n in self.fields if n not in self.wind + ("mean_direction",)]
            speeds = ensemble.speed_names(self.used["percentiles"])
        files, geo, note = grids.export_grids(self.lib, out_dir, self.prefix(), self.fields, self.wind, extra, speeds,
                                              self.units.currentText())
        self.caption.setText("wrote %d files to %s (corner and .prj from %s%s)"
                             % (len(files), out_dir, geo.source, "; " + note if note else ""))
        return files

    def export_points(self, points_file, out_dir):
        """A CSV of the values at the points of points_file ('x y [name]' per line) in out_dir."""
        points = grids.read_points(points_file)
        xy = [(x, y) for _, x, y in points]
        u = self.used                                                # the wind shown, not the boxes now
        h = u["height"]
        if self.wind == ONE_WIND:
            values = ensemble.wind_at_points(self.lib, u["direction"], u["speed"], h, xy)
        else:
            values = ensemble.at_points(self.lib, u["members"], h, xy, u["percentiles"], u["thresholds"])
        path = os.path.join(out_dir, self.prefix() + "_points.csv")
        files = grids.export_points(self.lib, path, points, values, h)
        self.caption.setText("wrote %s" % path)
        return files

    def _browse(self):
        path = QtWidgets.QFileDialog.getExistingDirectory(self, "Case folder", self.folder.text())
        if path:
            self.folder.setText(path)
            try:
                self.load(path)
            except LibraryError:
                pass

    def _export_grids_dialog(self):
        out = QtWidgets.QFileDialog.getExistingDirectory(self, "Folder for the grids", self.folder.text())
        if out:
            try:
                self.export_grids(out)
            except (LibraryError, OSError) as err:
                self.caption.setText(str(err))

    def _export_points_dialog(self):
        points, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Points file (x y [name] per line, m, case frame)",
                                                          self.folder.text(), "Text (*.txt *.csv);;All files (*)")
        if not points:
            return
        out = QtWidgets.QFileDialog.getExistingDirectory(self, "Folder for the CSV", os.path.dirname(points))
        if out:
            try:
                self.export_points(points, out)
            except (LibraryError, SampleError, OSError) as err:
                self.caption.setText(str(err))
