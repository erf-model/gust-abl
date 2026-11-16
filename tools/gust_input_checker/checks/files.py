"""Files the deck names, and the output keys.

The ``files`` check reports, under its own name so that a caller can leave
them out (a test that writes its files at run time), the files and
directories a key names that are not there. GUST opens a relative name in the
directory the run starts in; a file found only beside the deck is a WARNING
that says so.

ERROR:

- a file that does not exist: terrain.file (Source/TerrainSurface.cpp
  "cannot be read"), mesh.zlevels_file (ZLevels.cpp), bc.inflow_file
  (InflowProfile.cpp), bc.surface_z0_file (SurfaceLayer.cpp), canopy.file,
  canopy.height_file, canopy.lai_file, canopy.cd_file (Canopy.cpp),
  buildings.file (Buildings.cpp), unsteady.probes_file and steady.probes_file
  (UnsteadyOutput.cpp);
  io.restart and bc.nest_parent, checkpoint directories with a ``Header``
  (StateIO.cpp, read_checkpoint_header);
- bc.inflow_file (InflowProfile.cpp, parse_inflow_table): rows of 4 or 5
  numbers, the same count in every row, the first height 0, heights that
  increase, theta > 0, k >= 0, at least two rows;
- mesh.zlevels_file (ZLevels.cpp): mesh.n_cell z + 1 numbers, the first 0,
  increasing;
- unsteady.probes_file and steady.probes_file (UnsteadyOutput.cpp,
  parse_probes): rows of three numbers "x y h" with h > 0, at least one
  probe; steady.probes_file must not be empty (Driver.cpp).

The ``output`` check: ERROR for io.turb_plot_file without a stage,
io.iteration_log with the unsteady stage and io.column_guess_file without
init.first_guess = column (Source/App.cpp), and for an empty io.plot_file or
io.chk_file (Source/Inputs.cpp); ADVICE for a long run without a checkpoint
(docs/running.rst: the checkpoint is what io.restart continues from).
"""
import os
import re

from ..checker import Checker, register

MAX_ROWS = 200000       # rows read from a table before the format check stops
LONG_RUN = 2000         # iterations above which a run without a checkpoint gets the advice


def _number(tok):
    """A token as strtod reads it, or None."""
    if "_" in tok:
        return None
    try:
        return float(tok)
    except ValueError:
        return None


def read_rows(path):
    """(rows, problem): the rows of numbers of a text table ('#' comments), or the first line that is not numbers."""
    rows = []
    try:
        with open(path, "r", encoding="utf-8-sig", errors="replace") as fh:
            for lineno, line in enumerate(fh, 1):
                toks = line.split("#", 1)[0].split()
                if not toks:
                    continue
                vals = [_number(t) for t in toks]
                if None in vals:
                    return rows, "line %d: '%s' is not a number" % (lineno, toks[vals.index(None)])
                rows.append((lineno, vals))
                if len(rows) >= MAX_ROWS:
                    break
    except OSError as e:
        return rows, "cannot be read (%s)" % e.strerror
    return rows, ""


def inflow_table_problem(rows):
    """The first rule of the bc.inflow_file table that `rows` break ('' when none)."""
    ncol, last = 0, None
    for lineno, v in rows:
        if len(v) not in (4, 5):
            return "line %d must hold z u v theta [k] (4 or 5 values)" % lineno
        ncol = ncol or len(v)
        if len(v) != ncol:
            return "line %d has a different number of columns than the first row" % lineno
        if last is not None and not v[0] > last:
            return "heights must increase strictly (line %d)" % lineno
        if last is None and v[0] != 0:
            return "the first height must be 0 (the ground)"
        if v[3] <= 0:
            return "theta must be > 0 (line %d)" % lineno
        if ncol == 5 and v[4] < 0:
            return "k must be >= 0 (line %d)" % lineno
        last = v[0]
    if len(rows) < 2:
        return "needs at least two rows"
    return ""


def inflow_table_theta(ctx):
    """The theta column of bc.inflow_file, or None when the table cannot be read here."""
    path = ctx.find_file(ctx.s("bc.inflow_file"))[0]
    if not path or os.path.isdir(path):
        return None
    rows, problem = read_rows(path)
    if problem or inflow_table_problem(rows):
        return None
    return [v[3] for _, v in rows]


@register
class FilesChecker(Checker):
    """Files and checkpoint directories the deck names: whether they exist, and the format of the small tables."""

    name = "files"

    def run(self):
        c = self.ctx
        if c.terrain_type in ("xyz", "erf_grid"):
            self._file("terrain.file")
        path = self._file("mesh.zlevels_file")
        if path:
            self._zlevels(path)
        if c.reads_bc and c.needs_profile and c.inflow_profile == "file":
            path = self._file("bc.inflow_file")
            if path:
                rows, problem = read_rows(path)
                problem = problem or inflow_table_problem(rows)
                if problem:
                    self._format("bc.inflow_file", problem, "The table holds 'z u v theta [k]' rows from z = 0 up.")
        if c.surface_layer:
            self._file("bc.surface_z0_file")
        kind = c.s("canopy.type")
        if kind == "patches":
            self._file("canopy.file")
        elif kind == "map":
            for key in ("canopy.height_file", "canopy.lai_file", "canopy.cd_file"):
                self._file(key)
        if c.buildings:
            self._file("buildings.file")
        if c.unsteady:
            path = self._file("unsteady.probes_file")
            if path:
                self._probes(path, "unsteady.probes_file")
        elif c.solve:
            if self.given("steady.probes_file") and not c.s("steady.probes_file"):
                self.error("steady.probes_file", "GUST stops: steady.probes_file must name a file.",
                           "Give the file, or remove the key.")
            path = self._file("steady.probes_file")
            if path:
                self._probes(path, "steady.probes_file")
        self._checkpoint("io.restart")
        if c.reads_bc:
            self._checkpoint("bc.nest_parent")

    def _file(self, key):
        """The path of the file `key` names when it is there; reports it when it is not."""
        if not self.deck.has(key):
            return None
        name = self.ctx.s(key)
        if not name:
            return None
        path, where = self.ctx.find_file(name)
        if path is None and self._made_by_run(key, name):
            return None
        if path is None:
            self.error(key, "GUST stops: %s = %s: the file does not exist (looked in %s and beside the deck)."
                       % (key, name, self.ctx.cwd), "Give the path from the directory the run starts in.")
            return None
        if os.path.isdir(path):
            self.error(key, "GUST stops: %s = %s is a directory, not a file." % (key, name), "Name the file.")
            return None
        if where == "deck":
            self.warn(key, "%s = %s is beside the deck, not in the current directory (%s): GUST opens it relative to "
                           "the directory the run starts in." % (key, name, self.ctx.cwd),
                      "Start the run in the deck's directory, or give the path from where the run starts.")
        return path

    def _checkpoint(self, key):
        if not self.deck.has(key):
            return
        name = self.ctx.s(key)
        if not name:
            return
        path, where = self.ctx.find_file(name)
        if path is None and self._made_by_run(key, name):
            return
        if path is None:
            self.error(key, "GUST stops: %s = %s: the checkpoint directory does not exist (looked in %s and beside the "
                            "deck)." % (key, name, self.ctx.cwd), "Give the checkpoint's directory.")
        elif not os.path.isfile(os.path.join(path, "Header")):
            self.error(key, "GUST stops: %s = %s holds no Header: it is not a GUST checkpoint." % (key, name),
                       "Give the directory a run wrote with io.write_checkpoint = 1.")
        elif where == "deck":
            self.warn(key, "%s = %s is beside the deck, not in the current directory (%s)." % (key, name, self.ctx.cwd),
                      "Start the run in the deck's directory, or give the path from where the run starts.")

    def _made_by_run(self, key, name):
        """INFO in place of the missing-file ERROR for a path the caller says the run writes (Options.made_by_run)."""
        if not self.ctx.made_by_run(name):
            return False
        self.info(key, "%s = %s is not checked: written by an earlier step of the run." % (key, name))
        return True

    def _format(self, key, problem, advice):
        self.error(key, "GUST stops: %s = %s: %s." % (key, self.ctx.s(key), problem), advice)

    def _zlevels(self, path):
        rows, problem = read_rows(path)
        advice = "The file holds mesh.n_cell z + 1 increasing node heights, the first 0."
        if problem:
            self._format("mesh.zlevels_file", problem, advice)
            return
        z = [x for _, v in rows for x in v]
        nz = self.ctx.n_cell[2] if self.ctx.n_cell else None
        if nz is not None and len(z) != nz + 1:
            self._format("mesh.zlevels_file", "%d heights for mesh.n_cell z = %d (it needs %d)" % (len(z), nz, nz + 1),
                         advice)
        elif z and z[0] != 0:
            self._format("mesh.zlevels_file", "the first height is %g, not 0" % z[0], advice)
        elif any(b <= a for a, b in zip(z, z[1:])):
            self._format("mesh.zlevels_file", "the heights do not increase strictly", advice)

    def _probes(self, path, key):
        probes = 0
        advice = "The file holds one 'x y h' line per probe, h the height above the ground (m)."
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as fh:
                lines = fh.read().splitlines()[:MAX_ROWS]
        except OSError:
            return
        for lineno, line in enumerate(lines, 1):
            toks = line.split("#", 1)[0].split()
            if not toks:
                continue
            vals = [_number(t) for t in toks]
            if len(toks) != 3:
                self._format(key, "line %d: expected three numbers 'x y h', found %d fields"
                             % (lineno, len(toks)), advice)
                return
            if None in vals:
                self._format(key, "line %d: '%s' is not a number" % (lineno, toks[vals.index(None)]),
                             advice)
                return
            if not vals[2] > 0:
                self._format(key, "line %d: h = %s: the height above the ground must be > 0"
                             % (lineno, toks[2]), advice)
                return
            probes += 1
            lo, hi = self.ctx.prob_lo, self.ctx.prob_hi
            for d, axis in enumerate("xy"):
                if lo and hi and (vals[d] < lo[d] or vals[d] > hi[d]):
                    self._format(key, "probe %d (line %d): %s = %s is outside the domain [%.6g, %.6g]"
                                 % (probes, lineno, axis, toks[d], lo[d], hi[d]),
                                 "Give probes inside mesh.prob_lo and mesh.prob_hi.")
                    return
        if probes == 0:
            self._format(key, "no probe ('x y h' lines)", advice)


@register
class OutputChecker(Checker):
    """The io.* keys: names of the output, and outputs that the run does not write."""

    name = "output"

    def run(self):
        c = self.ctx
        for key in ("io.plot_file", "io.chk_file"):
            if self.deck.has(key) and not c.s(key):
                self.error(key, "GUST stops: %s must not be empty." % key, "Give a name, or remove the key.")
        if not c.solve:
            self.no_effect("io.turb_plot_file", "without steady.stage = simplec")
            self._overwrite_restart()
        if c.unsteady:
            self.no_effect("io.iteration_log", "with steady.stage = unsteady (one line per physical step)")
        if c.first_guess != "column":
            self.no_effect("io.column_guess_file", "without init.first_guess = column")
        if c.solve and not c.unsteady and c.i("io.write_checkpoint") == 0:
            iters = c.i("steady.max_iterations")
            if iters is not None and iters >= LONG_RUN and c.n_cell and c.n_cell[0] * c.n_cell[1] * c.n_cell[2] >= 100000:
                self.advise("io.write_checkpoint", "io.write_checkpoint = 0 on a run of up to %d iterations: a run that "
                                                   "stops early cannot be continued with io.restart." % iters,
                            "Leave io.write_checkpoint at 1 unless the run is short.")

    def _overwrite_restart(self):
        """Without a solve the step stays the restart's, so the final checkpoint <io.chk_file><step> is the
        directory read when io.restart is <io.chk_file> followed by its 5-digit step (App.cpp check_overwrite)."""
        c = self.ctx
        restart = c.s("io.restart")
        if not restart or c.i("io.write_checkpoint") == 0:
            return
        chk = c.s("io.chk_file") or "chk"
        norm = os.path.normpath(restart)
        m = re.match(r"^(.*?)(\d{5,})$", norm)
        if m and os.path.normpath(chk) == os.path.normpath(m.group(1) or "."):
            self.error("io.chk_file", "GUST stops: io.chk_file = %s would overwrite io.restart = %s (without a solve "
                                      "the run's step is the checkpoint's)." % (chk, restart),
                       "Choose another io.chk_file or set io.write_checkpoint = 0.")
