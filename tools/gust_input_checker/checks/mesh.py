"""The mesh: domain, vertical levels, grid hierarchy, terrain and box layout.

ERROR (each is an abort of Source/Inputs.cpp validate, of
Source/MeshHierarchy.cpp check_hierarchy_inputs or of Source/Sequencing.cpp
check_sequencing_options):

- mesh.prob_hi not above mesh.prob_lo; a periodic z; mesh.prob_lo z not 0;
- mesh.stretch_ratio without mesh.dz0 > 0; mesh.zlevels_file with mesh.dz0;
  mesh.dz0 above prob_hi z / n_cell z when no ratio is given;
- mesh.hierarchy_levels = L > 1 with n_cell or max_grid_size in x, y not
  divisible by 2^(L-1), or a coarsest level of 1 x 1 cells;
  terrain.hierarchy_passes with one level;
- steady.sequencing = 1 with one level, with a coarsest level of fewer than
  2 cells in x or y, with sequencing_drop_orders above residual_drop_orders,
  or with sequencing_cfl_start = previous and sequencing_hold = 1; the
  steady.sequencing_* keys without steady.sequencing = 1;
- simplec.theta = frozen without grid sequencing on a fresh start
  (Source/App.cpp);
- terrain.file missing for xyz and erf_grid or given for flat and agnesi;
  terrain.gamma_m without terrain.smoothing = 1; the terrain.hill_* keys
  without agnesi and agnesi without them; the terrain.idw_* keys when nothing
  is interpolated from points.

ADVICE: a stretch ratio above 1.2 and a last layer more than 100 first layers
thick (no sentence of the documentation sets a limit, so these are not
warnings); a large terrain mesh solved on one level (docs/howto.rst, "Larger
terrain meshes: add grid sequencing"); levels without steady.sequencing = 1
(reference, mesh.hierarchy_levels: "The solver uses the coarse levels only
with steady.sequencing = 1"); boxes and ranks (docs/howto.rst, "Many small
boxes"; "Run one MPI rank per GPU"; "A GPU build keeps z-line Jacobi").

INFO: cell sizes, layer heights, the top, cell and box counts, the coarsest
level.
"""
import math

from ..checker import Checker, register
from ..context import geometric_ratio, gust_newton_converges

LARGE_MESH = 500000     # cells; above it a terrain run on one level gets the grid-sequencing advice


def _g(x):
    return "%.6g" % x


@register
class MeshChecker(Checker):
    """Domain corners, periodicity and the vertical levels (uniform, stretched or from a file)."""

    name = "mesh"

    def run(self):
        c = self.ctx
        lo, hi = c.prob_lo, c.prob_hi
        if lo and hi:
            for d, name in enumerate("xyz"):
                if not hi[d] > lo[d]:
                    self.error("mesh.prob_hi", "GUST stops: mesh.prob_hi %s = %s is not above mesh.prob_lo %s = %s."
                               % (name, _g(hi[d]), name, _g(lo[d])), "Give a high corner above the low corner.")
            if lo[2] != 0:
                self.error("mesh.prob_lo", "GUST stops: mesh.prob_lo z = %s: the z value must be 0 (node heights "
                                           "start at the flat ground)." % _g(lo[2]), "Set the z value to 0.")
        for key in ("mesh.prob_lo", "mesh.prob_hi"):
            if any(t.lower().lstrip("+-") in ("inf", "infinity", "nan") for t in c.tokens(key)):
                self.error(key, "GUST stops: mesh.prob_lo and mesh.prob_hi must be finite (%s = %s)."
                           % (key, " ".join(c.tokens(key))), "Give finite corners (m).")
        if c.periodic and c.periodic[2] == 1:
            self.error("mesh.is_periodic", "GUST stops: mesh.is_periodic z = 1: z cannot be periodic.",
                       "Set the third value to 0.")
        dz0, ratio = c.r("mesh.dz0"), c.r("mesh.stretch_ratio")
        if self.given("mesh.stretch_ratio") and not (dz0 is not None and dz0 > 0):
            self.error("mesh.stretch_ratio", "GUST stops: mesh.stretch_ratio needs mesh.dz0 > 0.",
                       "Give the first layer's thickness, mesh.dz0.")
        if self.given("mesh.zlevels_file") and dz0 is not None and dz0 > 0:
            self.error("mesh.zlevels_file", "GUST stops: mesh.zlevels_file and mesh.dz0 cannot both be given.",
                       "Keep the file or mesh.dz0 (with mesh.stretch_ratio).")
        if (dz0 is not None and dz0 > 0 and not self.given("mesh.stretch_ratio") and c.n_cell and hi
                and dz0 > hi[2] / c.n_cell[2]):
            self.error("mesh.dz0", "GUST stops: mesh.dz0 = %s exceeds prob_hi z / n_cell z = %s (a stretched grid "
                                   "starts thinner than uniform)." % (_g(dz0), _g(hi[2] / c.n_cell[2])),
                       "Lower mesh.dz0, or give mesh.stretch_ratio to let the layers set the top.")
        elif (dz0 is not None and dz0 > 0 and not self.given("mesh.stretch_ratio") and c.n_cell and hi
              and not gust_newton_converges(dz0, hi[2], c.n_cell[2])):
            r = geometric_ratio(dz0, hi[2], c.n_cell[2])
            self.error("mesh.dz0", "GUST stops: with mesh.dz0 = %s and %d layers to mesh.prob_hi z = %s, GUST's "
                                   "Newton iteration (ERF's geometric_stretch, from 1.03) does not find the ratio%s."
                       % (_g(dz0), c.n_cell[2], _g(hi[2]), "" if r is None else " (%s)" % _g(r)),
                       "Give mesh.stretch_ratio%s, or use more layers or a thicker first layer."
                       % ("" if r is None else " = %.12g" % r))
        self._report(dz0, ratio)

    def _report(self, dz0, ratio):
        c = self.ctx
        if c.n_cell and c.dx:
            n = c.n_cell
            self.info("mesh.n_cell", "%d x %d x %d = %d cells; dx = %s m, dy = %s m."
                      % (n[0], n[1], n[2], n[0] * n[1] * n[2], _g(c.dx), _g(c.dy)))
        z = c.zlevels
        if not z:
            return
        first, last, top = c.dz_first, c.dz_last, z[-1]
        how = {"uniform": "uniform layers", "fixed_top": "mesh.dz0 stretched to mesh.prob_hi z",
               "ratio": "mesh.dz0 with mesh.stretch_ratio", "file": "mesh.zlevels_file"}[c.zlevels_source]
        text = "%s: first layer %s m, last layer %s m, top at %s m" % (how, _g(first), _g(last), _g(top))
        if c.stretch is not None and c.zlevels_source != "uniform":
            text += ", ratio %s" % _g(c.stretch)
        self.info("mesh.dz0" if c.zlevels_source in ("fixed_top", "ratio") else "mesh.n_cell", text + ".")
        if c.zlevels_source == "ratio" and c.prob_hi and abs(top - c.prob_hi[2]) > 1e-6 * max(top, c.prob_hi[2]):
            self.info("mesh.stretch_ratio", "The layers end at %s m, not at mesh.prob_hi z = %s m: with "
                                            "mesh.stretch_ratio the top is where the layers end."
                      % (_g(top), _g(c.prob_hi[2])),
                      "Heights such as bc.rayleigh_zdamp are measured against %s m." % _g(top))
        if c.dx and first > 0:
            self.info("mesh.prob_hi", "Largest cell aspect ratio over flat ground max(dx, dy) / first layer = %s."
                      % _g(max(c.dx, c.dy) / first))
        if c.zlevels_source in ("fixed_top", "ratio") and c.stretch is not None and c.stretch > 1.2:
            self.advise("mesh.stretch_ratio" if c.zlevels_source == "ratio" else "mesh.dz0",
                        "The layers grow by a ratio of %s; the documented cases use 1.08 to 1.12." % _g(c.stretch),
                        "More layers (mesh.n_cell z) or a thicker first layer give a ratio nearer 1.1.")
        elif first > 0 and last / first > 100 and c.zlevels_source != "uniform":
            self.advise("mesh.dz0", "The last layer (%s m) is %s times the first (%s m)."
                        % (_g(last), _g(last / first), _g(first)), "Check that the upper layers resolve what is there.")


@register
class HierarchyChecker(Checker):
    """mesh.hierarchy_levels against the cell counts and boxes; the grid-sequencing keys."""

    name = "hierarchy"

    def run(self):
        c = self.ctx
        L, n, mg = c.levels, c.n_cell, c.max_grid
        r = 2 ** (L - 1)
        if L > 1 and n:
            for d, name in ((0, "x"), (1, "y")):
                if n[d] % r != 0:
                    self.error("mesh.n_cell", "GUST stops: mesh.n_cell %s = %d is not divisible by %d "
                                              "(mesh.hierarchy_levels = %d)." % (name, n[d], r, L),
                               "Use a multiple of %d cells in x and y, or fewer levels." % r)
            if n[0] % r == 0 and n[1] % r == 0 and n[0] // r == 1 and n[1] // r == 1:
                self.error("mesh.hierarchy_levels", "GUST stops: the coarsest of %d levels has 1 x 1 cells." % L,
                           "Use fewer levels.")
            if mg and (mg[0] % r != 0 or mg[1] % r != 0):
                self.error("mesh.max_grid_size", "GUST stops: mesh.max_grid_size = %d %d is not a multiple of %d in x "
                                                 "and y (mesh.hierarchy_levels = %d)." % (mg[0], mg[1], r, L),
                           "Use a multiple of %d: every box holds whole cells of the coarsest level." % r)
            if n[0] % r == 0 and n[1] % r == 0:
                self.info("mesh.hierarchy_levels", "%d levels; the coarsest has %d x %d x %d cells."
                          % (L, n[0] // r, n[1] // r, n[2]))
        if L == 1:
            self.no_effect("terrain.hierarchy_passes", "with mesh.hierarchy_levels = 1",
                           "Remove it, or set mesh.hierarchy_levels.")
        if c.solve:
            self._sequencing(L, n, r)
        if c.theta_frozen and not (c.sequencing and not c.restart and L > 1):
            self.error("simplec.theta", "GUST stops: simplec.theta = frozen needs grid sequencing (steady.sequencing = 1, "
                                        "mesh.hierarchy_levels > 1) on a fresh start.",
                       "Set mesh.hierarchy_levels >= 2 and steady.sequencing = 1.")
        if c.solve and n and L == 1 and c.terrain_type in ("xyz", "erf_grid") and n[0] * n[1] * n[2] >= LARGE_MESH \
                and not c.buildings and not c.unsteady:
            self.advise("mesh.hierarchy_levels", "A terrain mesh of %d cells is solved on one level." % (n[0] * n[1] * n[2]),
                        "Grid sequencing (mesh.hierarchy_levels = 3, steady.sequencing = 1) gives the fine level a "
                        "better start (docs/howto.rst, Converge faster).")
        if c.solve and L > 1 and not c.sequencing and not c.unsteady:
            self.advise("steady.sequencing", "mesh.hierarchy_levels = %d without steady.sequencing = 1: the solver uses "
                                             "the coarse levels only with grid sequencing." % L,
                        "Set steady.sequencing = 1, or mesh.hierarchy_levels = 1.")

    def _sequencing(self, L, n, r):
        c = self.ctx
        if c.sequencing:
            if L < 2:
                self.error("steady.sequencing", "GUST stops: steady.sequencing = 1 needs mesh.hierarchy_levels >= 2.",
                           "Set mesh.hierarchy_levels (n_cell and max_grid_size in x, y divisible by 2^(L-1)).")
            elif n and n[0] % r == 0 and n[1] % r == 0 and (n[0] // r < 2 or n[1] // r < 2):
                self.error("steady.sequencing", "GUST stops: steady.sequencing = 1 needs a coarsest level of at least "
                                                "2 cells in x and y (it has %d x %d)." % (n[0] // r, n[1] // r),
                           "Use fewer levels.")
            drop, total = c.r("steady.sequencing_drop_orders"), c.r("steady.residual_drop_orders")
            if drop is not None and total is not None and drop > total:
                self.error("steady.sequencing_drop_orders", "GUST stops: steady.sequencing_drop_orders = %s is above "
                                                            "steady.residual_drop_orders = %s." % (_g(drop), _g(total)),
                           "A coarse level stops at or before the drop the fine level needs.")
            if c.s("steady.sequencing_cfl_start") == "previous" and c.flag("steady.sequencing_hold"):
                self.error("steady.sequencing_cfl_start", "GUST stops: steady.sequencing_cfl_start = previous needs "
                                                          "steady.sequencing_hold = 0.", "Remove steady.sequencing_hold.")
        else:
            for key in ("sequencing_drop_orders", "sequencing_monitor_tol", "sequencing_hold", "sequencing_cfl_start"):
                self.no_effect("steady." + key, "without steady.sequencing = 1")


@register
class TerrainChecker(Checker):
    """terrain.type and the keys that go with each type."""

    name = "terrain"

    def run(self):
        c = self.ctx
        t = c.terrain_type
        if t in ("xyz", "erf_grid"):
            if not self.given("terrain.file"):
                self.error("terrain.file", "GUST stops: terrain.file is required with terrain.type = %s." % t,
                           "Name the terrain file.")
        elif t in ("flat", "agnesi"):
            self.no_effect("terrain.file", "with terrain.type = %s" % t)
        if c.i("terrain.smoothing") is not None and c.i("terrain.smoothing") != 1:
            self.no_effect("terrain.gamma_m", "without terrain.smoothing = 1 (the STF decay)")
        if t == "agnesi":
            for key in ("terrain.hill_height", "terrain.hill_half_width"):
                if not self.given(key):
                    self.error(key, "GUST stops: %s is required with terrain.type = agnesi." % key, "Add %s (m)." % key)
        elif t in ("flat", "xyz", "erf_grid"):
            for key in ("terrain.hill_height", "terrain.hill_half_width", "terrain.hill_shape"):
                self.no_effect(key, "with terrain.type = %s" % t)
        points = t == "xyz" or (c.reads_bc and self.given("bc.surface_z0_file")) or c.s("canopy.type") == "map"
        if not points and t in ("flat", "agnesi", "erf_grid"):
            for key in ("terrain.idw_neighbors", "terrain.idw_power"):
                self.no_effect(key, "without points to interpolate (terrain.type = xyz, bc.surface_z0_file or "
                                    "canopy.type = map)")


@register
class LayoutChecker(Checker):
    """Boxes from mesh.max_grid_size, and how they fall on the ranks of --ranks and a GPU build (--gpu)."""

    name = "layout"

    def run(self):
        c = self.ctx
        n, mg = c.n_cell, c.max_grid
        if not n or not mg:
            return
        bx, by = math.ceil(n[0] / mg[0]), math.ceil(n[1] / mg[1])
        boxes = bx * by
        sx, sy = math.ceil(n[0] / bx), math.ceil(n[1] / by)
        ranks = c.options.ranks
        text = "%d box%s (%d x %d, about %d x %d x %d cells each; boxes are never split in z)" \
               % (boxes, "" if boxes == 1 else "es", bx, by, sx, sy, n[2])
        if ranks:
            text += "; %s boxes per rank on %d ranks" % (_g(boxes / ranks), ranks)
        self.info("mesh.max_grid_size", text + ".")
        if ranks and boxes < ranks:
            self.advise("mesh.max_grid_size", "%d box%s for %d ranks: %d rank%s would hold no box."
                        % (boxes, "" if boxes == 1 else "es", ranks, ranks - boxes, "" if ranks - boxes == 1 else "s"),
                        "Lower mesh.max_grid_size (x, y) so that there are at least %d boxes, or use fewer ranks." % ranks)
        elif ranks and boxes % ranks != 0:
            self.advise("mesh.max_grid_size", "%d boxes on %d ranks: the ranks hold unequal numbers of boxes."
                        % (boxes, ranks), "A box count that is a multiple of the ranks balances the load.")
        lusgs = c.solve and c.sweep == "lusgs"
        if lusgs and min(sx, sy) < 16 and boxes > 1:
            self.advise("mesh.max_grid_size", "Boxes of about %d x %d cells with steady.sweep = lusgs: the LU-SGS sweeps "
                                              "lag the coupling across box faces, so very small boxes cost extra "
                                              "iterations." % (sx, sy), "Use larger boxes where the ranks allow it.")
        if c.options.gpu:
            if c.solve and not self.given("steady.sweep"):
                self.advise("steady.sweep", "A GPU build takes steady.sweep = zjacobi by default (a CPU build lusgs), "
                                            "with a CFL cap of 10.",
                            "Give steady.sweep, steady.zline_sweeps and steady.cfl_max to run the same settings on "
                            "both builds.")
            if ranks and boxes > ranks:
                self.advise("mesh.max_grid_size", "%d boxes on %d GPU rank%s: one large box per GPU keeps the kernels "
                                                  "few and large." % (boxes, ranks, "" if ranks == 1 else "s"),
                            "Raise mesh.max_grid_size towards one box per rank (one MPI rank per GPU).")
            if c.solve and c.sweep == "lusgs" and c.s("steady.line_order") != "diagonal":
                self.advise("steady.line_order", "steady.sweep = lusgs on a GPU with steady.line_order = sequential.",
                            "steady.line_order = diagonal lets the LU-SGS sweeps run one column per thread.")
