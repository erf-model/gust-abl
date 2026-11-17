"""Nests run in the same process as the outer domain (nest.*).

ERROR (each is an abort of Source/Nests.cpp read_nest_options and
check_nests, or of Source/App.cpp before the outer domain starts):

- nest.names: a name that is not an identifier (a letter, then letters,
  digits or '_') or is 'names' or 'feedback'; a name given twice;
- nest.feedback other than none or defect, or without nest.names; with
  defect: steady.stage = simplec_newton, io.restart, bc.rayleigh_fields with
  u, v or tke, an odd mesh.max_grid_size that splits a nest into boxes
  with an odd number of cells in x or y, or a nest spanning fewer than
  2 DEFECT_MARGIN + 1 outer cells in a non-periodic x or y (Source/App.cpp);
- nest.<name>.lo or nest.<name>.hi missing for a name of nest.names;
  nest.<name>.* for a name not in nest.names (GUST does not read it);
- a corner that is not finite, off the outer cell faces (1e-4 cells), hi not
  above lo, or less than 2 outer cells (kNestMargin) from a non-periodic
  outer face; a nest over the whole outer domain in x and y (both periodic);
  two nests sharing outer cells (touching is allowed);
- with nests: steady.stage other than simplec or simplec_newton,
  io.write_checkpoint = 0, buildings, simplec.theta = frozen,
  terrain.smoothing = 1 (the default) with terrain other than flat;
- steady.line_overlap above a nest's narrowest box width in x and y (the nest
  mesh cut by mesh.max_grid_size; Source/App.cpp checks each nest's layout
  before the outer domain starts).

The value counts (two values per corner) are in `values`.

WARNING: nest.feedback = defect with bc.sponge_width > 0 (the nests have no
side sponge; GUST prints a note).

INFO: each nest's cells (refinement 2 2 1).
"""
import math
import re

from ..checker import Checker, register

MARGIN = 2      # NestParent.H, kNestMargin
RATIO = 2       # NestParent.H, kNestRatio
DEFECT_MARGIN = 2   # NestDefect.H, kDefectMargin
_NAME = re.compile(r"^[A-Za-z][A-Za-z0-9_]*$")


def _whole_cells(x, x0, dx):
    """The position of x in outer cells from x0, or None when it is not a whole number of cells."""
    c = (x - x0) / dx
    if not math.isfinite(c) or abs(c) > 1e9:
        return None
    n = int(round(c))
    return n if abs(c - n) <= 1e-4 else None


@register
class NestsChecker(Checker):
    """The nests' names and corners, their placement in the outer mesh, and what the outer deck needs with nests."""

    name = "nests"

    def run(self):
        c = self.ctx
        names = c.tokens("nest.names") if self.given("nest.names") else []
        valid = []
        for n in names:
            if not _NAME.match(n) or n in ("names", "feedback"):
                self.error("nest.names", "GUST stops: nest.names: '%s' is not a nest name (a letter, then letters, "
                                         "digits or '_'; not 'names' or 'feedback')." % n, "Rename the nest.")
            elif n in valid:
                self.error("nest.names", "GUST stops: nest.names: '%s' is named more than once." % n,
                           "Name each nest once.")
            else:
                valid.append(n)
        for key in self.deck.keys("nest."):
            parts = key.split(".")
            if len(parts) == 3 and parts[1] not in valid:
                self.error(key, "GUST stops: unknown input '%s': '%s' is not in nest.names." % (key, parts[1]),
                           "Add the name to nest.names, or remove the key.")
        defect = self._feedback(names)
        if not names:
            return
        self._outer()
        spans = []
        for n in valid:
            corners = {}
            for end in ("lo", "hi"):
                key = "nest.%s.%s" % (n, end)
                if not self.given(key):
                    self.error(key, "GUST stops: %s is missing (the nest's %s corner: x y in m)."
                               % (key, "low" if end == "lo" else "high"), "Give two values, x y (m).")
                    continue
                vals = c.rs(key)
                if vals is not None and len(vals) == 2:
                    corners[end] = vals
            if len(corners) == 2:
                span = self._place(n, corners["lo"], corners["hi"])
                if span is not None and span == [0, c.n_cell[0], 0, c.n_cell[1]]:
                    self.error("nest.%s.lo" % n, "GUST stops: nest.%s covers the whole outer domain in x and y (a nest "
                                                 "needs at least one nest face)." % n, "Shrink the nest.")
                    span = None
                if span is not None:
                    for other, s in spans:
                        if span[0] < s[1] and s[0] < span[1] and span[2] < s[3] and s[2] < span[3]:
                            self.error("nest.%s.lo" % n, "GUST stops: nest.%s overlaps nest.%s (nests at one level "
                                                         "may touch, not overlap)." % (n, other),
                                       "Move or shrink one of them.")
                    spans.append((n, span))
                    self._layout(n, span)
                    if defect:
                        self._even_boxes(n, span)
                    self.info("nest.%s.lo" % n, "nest.%s: %d x %d x %s cells of %.6g m x %.6g m (refinement 2 2 1)."
                              % (n, RATIO * (span[1] - span[0]), RATIO * (span[3] - span[2]),
                                 c.n_cell[2], c.dx / RATIO, c.dy / RATIO))

    def _feedback(self, names):
        """nest.feedback (Source/Nests.cpp, Source/App.cpp); True for defect."""
        c = self.ctx
        if not self.given("nest.feedback"):
            return False
        fb = c.s("nest.feedback")
        if fb not in ("none", "defect"):
            self.error("nest.feedback", "GUST stops: nest.feedback = '%s' is not none or defect." % fb,
                       "Set nest.feedback = none or defect.")
            return False
        if not names:
            self.error("nest.feedback", "GUST stops: nest.feedback = %s needs nest.names." % fb,
                       "Name the nests, or remove nest.feedback.")
            return False
        if fb != "defect":
            return False
        if c.stage == "simplec_newton":
            self.error("nest.feedback", "GUST stops: nest.feedback = defect needs steady.stage = simplec (not "
                                        "simplec_newton).", "Set steady.stage = simplec.")
        if c.restart:
            self.error("nest.feedback", "GUST stops: nest.feedback = defect does not restart (io.restart): run from "
                                        "the start.", "Remove io.restart.")
        if c.positive("bc.rayleigh_zdamp") and any(c.rayleigh_field(f) for f in ("u", "v", "tke")):
            self.error("nest.feedback", "GUST stops: nest.feedback = defect needs bc.rayleigh_fields without u, v and "
                                        "tke (the nests damp w and theta only).", "Set bc.rayleigh_fields = w theta.")
        if c.positive("bc.sponge_width"):
            self.warn("nest.feedback", "nest.feedback = defect with bc.sponge_width > 0: the nests have no side sponge, "
                                       "so over a footprint inside the sponge the joint run does not damp.",
                      "Place the nests outside bc.sponge_width.")
        return True

    def _even_boxes(self, n, span):
        """With defect: at least 2 DEFECT_MARGIN + 1 outer cells in each non-periodic direction, and every box of the
        nest an even number of cells in x and y (Source/App.cpp)."""
        c = self.ctx
        for d, ax in enumerate("xy"):
            cells = span[2 * d + 1] - span[2 * d]
            periodic = c.periodic and c.periodic[d] == 1 and cells == c.n_cell[d]
            if not periodic and cells < 2 * DEFECT_MARGIN + 1:
                self.error("nest.feedback", "GUST stops: nest.feedback = defect: nest.%s spans %d outer cells in %s; it "
                                            "needs at least %d (the source keeps %d cells from each face)."
                           % (n, cells, ax, 2 * DEFECT_MARGIN + 1, DEFECT_MARGIN), "Widen the nest.")
        if not c.max_grid:
            return
        for d in range(2):
            cells, m = RATIO * (span[2 * d + 1] - span[2 * d]), c.max_grid[d]
            nb = -(-cells // m)
            sizes = {cells // nb, -(-cells // nb)}
            if any(s % 2 for s in sizes):
                self.error("nest.feedback", "GUST stops: nest.feedback = defect: nest.%s splits into boxes with an odd "
                                            "number of cells in x or y at mesh.max_grid_size: give an even "
                                            "mesh.max_grid_size." % n, "Give an even mesh.max_grid_size.")
                return

    def _outer(self):
        """What the outer deck needs with nests (Source/App.cpp)."""
        c = self.ctx
        if c.stage not in ("simplec", "simplec_newton"):
            self.error("nest.names", "GUST stops: nest.names needs steady.stage = simplec or simplec_newton (each "
                                     "nest starts from the outer domain's steady field).",
                       "Set steady.stage, or remove the nests.")
        if c.i("io.write_checkpoint") == 0:
            self.error("io.write_checkpoint", "GUST stops: nest.names needs io.write_checkpoint = 1 (each nest starts "
                                              "from the outer domain's checkpoint).", "Remove io.write_checkpoint = 0.")
        if c.i("terrain.smoothing") == 1 and c.s("terrain.type") not in ("", "flat"):
            self.error("terrain.smoothing", "GUST stops: nest.names is not supported with terrain.smoothing = 1 over "
                                            "terrain (STF smooths the levels over each mesh's own ground, so a nest's "
                                            "heights above the ground differ from the outer domain's at its faces).",
                       "Set terrain.smoothing = 0 or 2.")
        if c.s("buildings.type") not in ("", "none"):
            self.error("nest.names", "GUST stops: nest.names is not supported with buildings.type yet.",
                       "Remove the buildings or the nests.")
        if c.theta_frozen:
            self.error("nest.names", "GUST stops: nest.names is not supported with simplec.theta = frozen (it needs "
                                     "grid sequencing, which a nest does not run).", "Solve theta, or remove the nests.")

    def _layout(self, n, span):
        """steady.line_overlap against the nest's narrowest box (the outer deck's _overlap_layout on the nest mesh)."""
        c = self.ctx
        overlap = c.i("steady.line_overlap")
        if c.sweep != "lusgs" or not overlap or overlap <= 0 or not c.max_grid:
            return
        widths = []
        for d in range(2):
            cells = RATIO * (span[2 * d + 1] - span[2 * d])
            if cells > 1 and c.max_grid[d] >= 1:
                widths.append(cells // -(-cells // c.max_grid[d]))
        if widths and overlap > min(widths):
            self.error("steady.line_overlap", "GUST stops: nest.%s: steady.line_overlap = %d: must be <= the narrowest "
                                              "box's width in x and y (%d cells; mesh.max_grid_size)."
                       % (n, overlap, min(widths)), "Lower steady.line_overlap, or widen the nest.")

    def _place(self, n, lo, hi):
        """The outer cells the nest spans (x0, x1, y0, y1), or None after reporting why it cannot be placed."""
        c = self.ctx
        if not (c.n_cell and c.extent and c.prob_lo):
            return None                 # the mesh checker reports the mesh
        key = "nest.%s.lo" % n
        span = []
        for d, ax in enumerate("xy"):
            if not (math.isfinite(lo[d]) and math.isfinite(hi[d])):
                self.error(key, "GUST stops: nest.%s.lo and nest.%s.hi %s must be finite." % (n, n, ax),
                           "Give finite corners (m).")
                return None
            dx = c.extent[d] / c.n_cell[d]
            o, e = _whole_cells(lo[d], c.prob_lo[d], dx), _whole_cells(hi[d], c.prob_lo[d], dx)
            if o is None or e is None:
                self.error(key, "GUST stops: nest.%s.lo and nest.%s.hi %s must lie on the outer cell faces (cells "
                                "of %.6g m from %.6g m)." % (n, n, ax, dx, c.prob_lo[d]),
                           "Move the corner to a multiple of %.6g m from %.6g m." % (dx, c.prob_lo[d]))
                return None
            if e <= o:
                self.error("nest.%s.hi" % n, "GUST stops: nest.%s.hi %s must be above nest.%s.lo %s." % (n, ax, n, ax),
                           "Give a high corner above the low corner.")
                return None
            m = 0 if c.periodic and c.periodic[d] == 1 else MARGIN
            N = c.n_cell[d]
            if o < m or e > N - m:
                self.error(key, "GUST stops: nest.%s must lie at least %d outer cells inside the outer %s faces "
                                "(outer cells %d to %d; the nest spans %d to %d)." % (n, m, ax, m, N - m, o, e),
                           "Move the nest inwards or enlarge the outer domain.")
                return None
            span += [o, e]
        return span
