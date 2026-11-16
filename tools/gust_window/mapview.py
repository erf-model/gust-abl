"""The plan view of a case: the terrain, the domain, the border blend and the nests, with the layers beside it.

x is east and y north in the case frame (metres from domain.center, or from the south-west corner with
domain.origin = southwest), as in the decks.
"""
import math

import numpy as np
import pyqtgraph as pg
from PySide6 import QtCore, QtWidgets

IMAGE_CELLS = 300                # the terrain image's cells along its longer side

pg.setConfigOptions(background="w", foreground="k", antialias=True)


def terrain_image(terrain, cells=IMAGE_CELLS):
    """The terrain points (scattered: a lon/lat raster in the case frame) as an image: the mean height of the
    points in each cell of a regular grid, rows south to north; cells with no point are NaN (drawn empty).
    Returns the image and its rectangle (x, y, width, height)."""
    x0, x1, y0, y1 = terrain.x.min(), terrain.x.max(), terrain.y.min(), terrain.y.max()
    # the points' own spacing: a lon/lat raster is coarser in one direction than terrain.spacing says
    own = math.sqrt(max((x1 - x0) * (y1 - y0), 1.0) / terrain.x.size)
    size = max(max(x1 - x0, y1 - y0) / cells, 1.5 * own)
    nx, ny = max(1, int(math.ceil((x1 - x0) / size))), max(1, int(math.ceil((y1 - y0) / size)))
    ye, xe = y0 + size * np.arange(ny + 1), x0 + size * np.arange(nx + 1)
    ye[-1], xe[-1] = max(ye[-1], y1), max(xe[-1], x1)      # rounding must not leave the last row of points out
    total, _, _ = np.histogram2d(terrain.y, terrain.x, bins=(ye, xe), weights=terrain.z)
    count, _, _ = np.histogram2d(terrain.y, terrain.x, bins=(ye, xe))
    with np.errstate(invalid="ignore", divide="ignore"):
        image = np.where(count > 0, total / np.maximum(count, 1), np.nan)
    return image, (x0, y0, xe[-1] - x0, ye[-1] - y0)


def layer_faces(plan):
    """Heights of the layer faces, 0 .. top: uniform, or dz0 r^k scaled to the top (the plan's stretch)."""
    n = plan.layers
    if plan.dz0 is None or plan.stretch_ratio == 1.0:
        return np.linspace(0.0, plan.top, n + 1)
    thick = plan.dz0 * plan.stretch_ratio ** np.arange(n)
    faces = np.concatenate([[0.0], np.cumsum(thick)])
    return faces * (plan.top / faces[-1])


def _rect(x0, x1, y0, y1, pen):
    item = pg.PlotDataItem([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], pen=pen)
    return item


class PlanView(QtWidgets.QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        self.caption = QtWidgets.QLabel("Check the case to draw it.")
        self.caption.setWordWrap(True)
        layout.addWidget(self.caption)
        self.graphics = pg.GraphicsLayoutWidget()
        layout.addWidget(self.graphics)
        self.plan_plot = self.graphics.addPlot(row=0, col=0)
        self.plan_plot.setAspectLocked(True)
        self.plan_plot.setLabel("bottom", "x (m, east)")
        self.plan_plot.setLabel("left", "y (m, north)")
        self.plan_plot.addLegend(offset=(5, 5))
        self.layer_plot = self.graphics.addPlot(row=0, col=1)
        self.layer_plot.setLabel("left", "z (m)")
        self.layer_plot.hideAxis("bottom")
        self.layer_plot.setMaximumWidth(140)
        self.layer_plot.setMouseEnabled(x=False, y=True)
        self.color_bar = pg.ColorBarItem(values=(0.0, 1.0), colorMap=pg.colormap.get("viridis"), label="ground (m)",
                                         interactive=False)
        self.graphics.addItem(self.color_bar, row=0, col=2)
        self.color_bar.hide()

    def clear(self, text="Check the case to draw it."):
        self.plan_plot.clear()
        self.layer_plot.clear()
        if self.color_bar is not None:
            self.graphics.removeItem(self.color_bar)
            self.color_bar = pg.ColorBarItem(values=(0.0, 1.0), colorMap=pg.colormap.get("viridis"), label="ground (m)",
                                         interactive=False)
        self.graphics.addItem(self.color_bar, row=0, col=2)
        self.color_bar.hide()
        self.caption.setText(text)

    def show_plan(self, plan, terrain=None):
        """Draw a resolved Plan and its prepared Terrain (or None)."""
        self.clear("")
        f = plan.frame
        p = self.plan_plot
        if terrain is not None and terrain.x.size:
            image, rect = terrain_image(terrain)
            cmap = pg.colormap.get("viridis")
            lo, hi = float(np.nanmin(image)), float(np.nanmax(image))
            item = pg.ImageItem(image, axisOrder="row-major")
            item.setColorMap(cmap)
            item.setLevels((lo, hi if hi > lo else lo + 1.0))
            item.setRect(*rect)
            p.addItem(item)
            self.color_bar.setLevels((lo, hi if hi > lo else lo + 1.0))
            self.color_bar.show()
        p.addItem(_rect(f.x_lo, f.x_hi, f.y_lo, f.y_hi, pg.mkPen("k", width=2)))
        p.legend.addItem(p.items[-1], "domain")
        flat, slope = plan.border_flat or 0.0, plan.border_slope or 0.0
        if flat > 0.0:
            p.addItem(_rect(f.x_lo + flat, f.x_hi - flat, f.y_lo + flat, f.y_hi - flat,
                            pg.mkPen((200, 120, 0), width=1.5)))
            p.legend.addItem(p.items[-1], "flat border")
        if slope > 0.0:
            w = flat + slope
            p.addItem(_rect(f.x_lo + w, f.x_hi - w, f.y_lo + w, f.y_hi - w,
                            pg.mkPen((200, 120, 0), width=1.5, style=QtCore.Qt.DashLine)))
            p.legend.addItem(p.items[-1], "end of the blend")
        for i, patch in enumerate(plan.patches):
            p.addItem(_rect(patch.x_lo, patch.x_hi, patch.y_lo, patch.y_hi, pg.mkPen((190, 30, 45), width=2)))
            if i == 0:
                p.legend.addItem(p.items[-1], "nests")
            b = patch.blend_width or 0.0
            if b > 0.0:
                p.addItem(_rect(patch.x_lo + b, patch.x_hi - b, patch.y_lo + b, patch.y_hi - b,
                                pg.mkPen((190, 30, 45), width=1, style=QtCore.Qt.DashLine)))
            label = pg.TextItem(patch.name, color=(190, 30, 45), anchor=(0, 1))
            label.setPos(patch.x_lo, patch.y_hi)
            p.addItem(label)
        # the wind: an arrow from the domain's centre, pointing where the first direction blows to
        if plan.directions:
            d = math.radians(plan.directions[0])
            cx, cy = 0.5 * (f.x_lo + f.x_hi), 0.5 * (f.y_lo + f.y_hi)
            length = 0.25 * min(f.x_hi - f.x_lo, f.y_hi - f.y_lo)
            tx, ty = -math.sin(d), -math.cos(d)
            tip = (cx + length * tx, cy + length * ty)
            head = 0.25 * length
            pen = pg.mkPen((30, 90, 200), width=2.5)
            p.addItem(pg.PlotDataItem([cx - length * tx, tip[0]], [cy - length * ty, tip[1]], pen=pen))
            for turn in (150.0, -150.0):                  # two barbs at the tip
                c, s_ = math.cos(math.radians(turn)), math.sin(math.radians(turn))
                p.addItem(pg.PlotDataItem([tip[0], tip[0] + head * (c * tx - s_ * ty)],
                                          [tip[1], tip[1] + head * (s_ * tx + c * ty)], pen=pen))
            p.legend.addItem(p.items[-3], "wind")
        p.autoRange()

        faces = layer_faces(plan)
        lp = self.layer_plot
        for zf in faces:
            lp.addItem(pg.PlotDataItem([0.0, 1.0], [zf, zf], pen=pg.mkPen((90, 90, 90), width=1)))
        lp.setXRange(0.0, 1.0)
        lp.setYRange(0.0, plan.top)
        lp.setTitle("%d layers" % plan.layers)

        text = ("%d x %d x %d cells of %.4g x %.4g m; first layer %.4g m, top %.6g m"
                % (plan.nx, plan.ny, plan.layers, plan.dx, plan.dy, faces[1] - faces[0], plan.top))
        if plan.directions:
            text += "; wind from %s degrees" % ", ".join("%g" % d for d in plan.directions)
        if terrain is not None:
            text += "; ground %.5g to %.5g m" % (terrain.lowest, terrain.highest)
        self.caption.setText(text)
