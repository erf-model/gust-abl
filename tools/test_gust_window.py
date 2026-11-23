"""Tests of the case set-up window, tools/gust_gui.py and tools/gust_window/ (CTest gust_window_selftest).

The parts with no Qt (field text, the log reader, the run steps) always run. The window's own tests need PySide6
and pyqtgraph and are skipped without them; they run with QT_QPA_PLATFORM=offscreen (set here when unset).

    - the form has one widget per key of the schema (test_one_widget_per_key)
    - a YAML file read into the form and saved again is unchanged (test_yaml_round_trip_*)
    - a refused key marks its own field (test_refused_key_marks_*)
    - the Library tab: one wind and a distribution drawn and exported; skipped without PySide6
      (LibraryTabTests)

    cd tools && python3 -m unittest -v test_gust_window
"""
import copy
import glob
import os
import stat
import sys
import tempfile
import textwrap
import time
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np  # noqa: E402
import yaml  # noqa: E402

import make_gust_case  # noqa: E402
from gust_case import CaseError, load_yaml, resume, validate  # noqa: E402
from gust_case.schema import SPEC  # noqa: E402
from gust_window import logparse, runplan, values  # noqa: E402

try:
    from PySide6 import QtCore, QtWidgets
    import pyqtgraph  # noqa: F401
    HAVE_QT = True
except ImportError:
    HAVE_QT = False

CASE_FILES = sorted(glob.glob(os.path.join(HERE, "..", "Tests", "Infrastructure", "Case_Setup", "*.yaml")))

# lines 18-24 of the log of Tests/Infrastructure/Case_Setup/flat_nest.yaml (gust at 02058da)
LOG = """\
gust: SIMPLEC stage: 2 lusgs sweeps (ghosts lagged, overlap 0, sequential order), pressure gmres_mlmg reltol 0.001, alpha_p 1, w_vnn_cap 16 (slope margin 4), cfl0 1, cfl_max 16
  Iter        CFL  Line solver      u-res      v-res      w-res    tke-res  theta-res      p-res          r
------  ---------  -----------  ---------  ---------  ---------  ---------  ---------  ---------  ---------
     1  1.000e+00  lusgs        3.525e-03  2.035e-03  8.766e-17  1.435e+00  6.758e-15  8.548e-17  1.000e+00
     2  1.000e+00  lusgs        1.435e-03  8.298e-04  6.336e-05  8.120e-01  1.344e-06  7.237e-08  1.000e+00
    10  4.149e+00  lusgs        6.260e-04  3.599e-04  1.535e-05  1.384e-01  1.776e-07  1.636e-08  2.423e-01
    20  8.775e+00  lusgs        2.060e-05  1.137e-05  6.525e-06  3.655e-03  2.002e-08  2.594e-09  1.030e-01
"""


def every_key_case():
    """A case description that gives every key (nests: two patches), each with a value of its kind."""
    def sample(key):
        if key.kind in ("number", "integer"):
            if key.words:
                return key.words[0]
            lo = key.lo if key.lo is not None else (key.hi - 2 if key.hi is not None else 1)
            hi = key.hi if key.hi is not None else lo + 3
            return int(lo) + 1 if key.kind == "integer" else lo + 0.5 * (hi - lo)
        if key.kind == "flag":
            return True
        if key.kind == "string":
            return "text 1"
        if key.kind == "path":
            return "data/dem.asc"
        if key.kind == "choice":
            return key.choices[-1]
        if key.kind == "numbers":
            lo = key.lo if key.lo is not None else 0.0
            return [lo + 1.0, lo + 2.5][:key.length or 2]
        if key.kind == "strings":
            return list(key.choices[:2]) if key.choices else ["a", "b"]
        if key.kind == "any":
            return {"steady.cfl0": 2.0, "turb.tke_floor": None, "io.plot_file": "out_plt"}
        raise AssertionError(key.kind)
    raw = {}
    for name, key in values.field_names():
        values.set_path(raw, name, sample(key))
    raw["nests"] = [{n: sample(k) for n, k in values.nest_spec().items()} for _ in range(2)]
    raw["nests"][1]["name"] = "second"
    return raw


# ---------------------------------------------------------------------------------------------------------------
# No Qt
# ---------------------------------------------------------------------------------------------------------------
class ValueTextTests(unittest.TestCase):
    def test_text_round_trip_by_kind(self):
        for name, key in list(values.field_names()) + [("nests[0]." + n, k) for n, k in values.nest_spec().items()]:
            raw = every_key_case()
            value = values.get_path(raw, name) if not name.startswith("nests") else raw["nests"][0][name[9:]]
            text = values.to_text(key, value)
            self.assertNotIn("\n...", text, name)
            self.assertEqual(values.from_text(name, key, text), (True, value), name)

    def test_empty_text_is_not_given(self):
        for name, key in values.field_names():
            self.assertEqual(values.from_text(name, key, "  "), (False, None), name)

    def test_texts_kept_as_typed(self):
        key = SPEC["case"]["name"]
        self.assertEqual(values.from_text("case.name", key, "012"), (True, "012"))      # not the number 12
        key = SPEC["model"]["tke_floor"]
        self.assertEqual(values.to_text(key, "1e-6"), "1e-6")                            # PyYAML's text for 1e-6
        self.assertEqual(values.from_text("model.tke_floor", key, "1e-6"), (True, "1e-6"))
        self.assertEqual(values.to_text(SPEC["solver"]["initial_blend"], [1.0, 0.5]), "1.0, 0.5")

    def test_lists_without_brackets(self):
        key = SPEC["domain"]["center_lat_lon"]
        for text in ("57.2, -7.4", "57.2 -7.4", "57.2,-7.4", "[57.2, -7.4]", " 57.2 ,  -7.4 "):
            self.assertEqual(values.from_text("domain.center_lat_lon", key, text), (True, [57.2, -7.4]), text)
        key = SPEC["domain"]["half_width"]
        self.assertEqual(values.from_text("domain.half_width", key, "3000"), (True, 3000))   # one number stays one
        key = SPEC["atmosphere"]["damping"]["fields"]
        self.assertEqual(values.from_text("atmosphere.damping.fields", key, "theta, u"), (True, ["theta", "u"]))

    def test_bad_yaml_names_the_key(self):
        with self.assertRaisesRegex(ValueError, r"^domain\.half_width = \[1,: not valid YAML"):
            values.from_text("domain.half_width", SPEC["domain"]["half_width"], "[1,")

    def test_problem_field(self):
        names = {"mesh.layers", "terrain.border.flat_width", "nests[1].half_width", "gust_inputs", "wind.speed"}
        pf = values.problem_field
        self.assertEqual(pf("mesh.layers = 0: must be an integer ≥ 1", names), "mesh.layers")
        self.assertEqual(pf("terrain.border.flat_width = -1: must be ...", names), "terrain.border.flat_width")
        self.assertEqual(pf("nests[1].half_width: required", names), "nests[1].half_width")
        self.assertEqual(pf("gust_inputs.steady.cfl0 = x: ...", names), "gust_inputs")
        self.assertIsNone(pf("terrain.border: flat_width + slope_width = 900 m reaches the middle", names))
        self.assertIsNone(pf("/tmp/case.yaml: not valid YAML", names))
        self.assertEqual(pf("a.b.c = 1: ...", {"a", "a.b"}), "a.b")                 # the longest field

    def test_unknown_parts(self):
        raw = {"domain": {"top": 1.0, "bogus": 2}, "extra": [1], "nests": [{"name": "a"}], "case": "scalar"}
        self.assertEqual(values.unknown_parts(raw), {"domain": {"bogus": 2}, "extra": [1], "case": "scalar"})


class LogReaderTests(unittest.TestCase):
    def test_rows_of_a_gust_log(self):
        reader = logparse.LogReader()
        rows = []
        for i in range(0, len(LOG), 37):                  # in pieces that split lines, as a pipe delivers them
            rows += reader.feed(LOG[i:i + 37])
        self.assertEqual([r.iteration for r in rows], [1, 2, 10, 20])
        self.assertEqual(rows[2].cfl, 4.149)
        self.assertEqual(rows[2].solver, "lusgs")
        self.assertEqual(rows[2].residuals, {"u": 6.260e-04, "v": 3.599e-04, "w": 1.535e-05, "tke": 1.384e-01,
                                             "theta": 1.776e-07, "p": 1.636e-08})
        self.assertEqual(rows[3].r, 1.030e-01)
        self.assertEqual([r.x for r in rows], [1, 2, 10, 20])

    def test_no_rows_before_a_header(self):
        reader = logparse.LogReader()
        self.assertEqual(reader.feed("     1  1.000e+00  lusgs  " + "  1.0e-3" * 7 + "\n"), [])

    def test_levels_dashes_and_notes(self):
        text = textwrap.dedent("""\
            level 1:   Iter        CFL  Line solver      u-res      v-res      w-res    tke-res  theta-res      p-res          r
            level 1:      1  1.000e+00  lusgs        1.0e-01  1.0e-01  1.0e-01  1.0e-01          -  1.0e-01  1.000e+00
            level 1:     30  1.600e+01  lusgs        1.0e-03  1.0e-03  1.0e-03  1.0e-03          -  1.0e-03  1.000e-02
                 1  1.000e+00  zjacobi      2.0e-02  2.0e-02  2.0e-02  2.0e-02          -  2.0e-02  1.000e+00  confidence 0.95 zmax 300
                 7        inf  newton       3.0e-04  3.0e-04  3.0e-04  3.0e-04          -  3.0e-04  4.000e-03  |R| 1 -> 0.5 lambda 1 attempts 1 gmres 12
            """)
        rows = logparse.LogReader().feed(text)
        self.assertEqual([(r.level, r.iteration, r.x) for r in rows], [(1, 1, 1), (1, 30, 30), (None, 1, 31), (None, 7, 37)])
        self.assertNotEqual(rows[0].residuals["theta"], rows[0].residuals["theta"])     # NaN: theta not solved
        self.assertEqual(rows[2].notes, "confidence 0.95 zmax 300")
        self.assertEqual(rows[3].cfl, float("inf"))
        self.assertEqual(rows[3].solver, "newton")

    def test_msvc_nan_spellings(self):
        head = "  Iter        CFL  Line solver      u-res      v-res      w-res    tke-res  theta-res      p-res          r\n"
        row = "     2  1.000e+00  lusgs      -nan(ind)   nan(ind)        NaN       -nan  1.0e-03  1.0e-03  1.000e+00\n"
        rows = logparse.LogReader().feed(head + row)
        self.assertEqual(len(rows), 1)
        for name in ("u", "v", "w", "tke"):
            self.assertNotEqual(rows[0].residuals[name], rows[0].residuals[name])
        self.assertEqual(rows[0].residuals["theta"], 1.0e-3)

    def test_column_wind(self):
        text = "gust: column\n  forcing.ug = 7.25\n  forcing.vg = -1.5\nlater\n  forcing.ug = 7.5\n"
        self.assertEqual(logparse.column_wind(text), ("7.5", "-1.5"))
        self.assertIsNone(logparse.column_wind("  forcing.ug = 7.5\n"))


class RunPlanTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name
        for sub in ("column", "nest_a", "nest_b"):
            os.makedirs(os.path.join(self.dir, sub))
        for rel in ("inputs", "column/inputs", "nest_a/inputs", "nest_b/inputs", resume.DATA):  # the record's files
            self.write(rel, "# deck\n")
        self.runs = [{"dir": ".", "label": "wind from 240 degrees", "column": True, "wind_from_column": True,
                      "patches": ["a", "b"], "chk": "chk"}]

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, rel, text="x"):
        with open(os.path.join(self.dir, rel), "w", encoding="utf-8") as out:
            out.write(text)

    def test_order_and_what_each_step_takes(self):
        gen = runplan.steps(self.dir, self.runs)
        column = next(gen)
        self.assertEqual((os.path.basename(column.cwd), column.args), ("column", ["inputs"]))
        self.write("column/column_guess.dat")
        self.write("column/log.txt", "  forcing.ug = 7.5\n  forcing.vg = -1.25\n")
        main = next(gen)
        self.assertEqual(os.path.normpath(main.cwd), os.path.normpath(self.dir))
        self.assertEqual(main.args, ["inputs", "forcing.ug=7.5", "forcing.vg=-1.25"])
        for name in ("chk00100", "chk99999", "chk100000", "chkpt"):
            os.makedirs(os.path.join(self.dir, name))
        a = next(gen)
        b = next(gen)
        self.assertEqual(a.args, ["inputs", "bc.nest_parent=../chk100000", "forcing.ug=7.5", "forcing.vg=-1.25"])
        self.assertEqual(os.path.basename(b.cwd), "nest_b")
        self.assertEqual(b.log, os.path.join(b.cwd, "log.txt"))
        with self.assertRaises(StopIteration):
            next(gen)

    def test_finished_direction_is_skipped_unless_resume_is_off(self):
        """A direction whose runs reached the steady verdict with these decks yields no step (run.resume), and runs
        again with skip_finished = False; the record is written when a direction starts."""
        run = dict(self.runs[0], column=False, wind_from_column=False, patches=[])
        gen = runplan.steps(self.dir, [run])
        next(gen)
        self.assertTrue(os.path.isfile(os.path.join(self.dir, resume.RECORD)))
        os.makedirs(os.path.join(self.dir, "plt00012"))
        self.write("plt00012/Header")
        self.write("log.txt", "gust: steady after 12 iterations (12 over the run)\ngust: wrote plt00012\n")
        notes = []
        self.assertEqual(list(runplan.steps(self.dir, [run], note=notes.append)), [])
        self.assertEqual(len(notes), 1)
        self.assertIn("skipped", notes[0])
        self.assertEqual(len(list(runplan.steps(self.dir, [run], skip_finished=False))), 1)

    def test_missing_deck_is_named(self):
        os.remove(os.path.join(self.dir, "nest_b", "inputs"))
        with self.assertRaisesRegex(runplan.RunError, r"wind from 240 degrees: .*nest_b.*write the case again"):
            next(runplan.steps(self.dir, self.runs))

    def test_missing_column_table(self):
        gen = runplan.steps(self.dir, self.runs)
        next(gen)
        with self.assertRaisesRegex(runplan.RunError, "column_guess.dat"):
            next(gen)

    def test_checkpoint_older_than_the_run_is_not_taken(self):
        os.makedirs(os.path.join(self.dir, "chk00500"))
        old = time.time() - 3600
        os.utime(os.path.join(self.dir, "chk00500"), (old, old))
        runs = [dict(self.runs[0], column=False, wind_from_column=False)]
        gen = runplan.steps(self.dir, runs)
        next(gen)
        with self.assertRaisesRegex(runplan.RunError, "wrote no checkpoint chk"):
            next(gen)

    def test_collect_library_after_the_last_step(self):
        self.assertEqual(runplan.collect_library(self.dir), "")          # no case_resolved.yaml: nothing to do
        self.write("case_resolved.yaml", "case: {name: x}\nwind: {drive: inflow_only, direction: 0}\n")
        self.assertEqual(runplan.collect_library(self.dir), "")          # no column run
        self.write("case_resolved.yaml", "case: {name: x}\nwind: {drive: geostrophic_column, direction: 0, "
                   "directions: null, speed: 8, height: 10}\n")
        note = runplan.collect_library(self.dir)                          # a geostrophic drive, its column not run
        self.assertTrue(note.startswith("library.json not written"), note)
        self.assertIn("column_guess.dat", note)

    def test_find_gust(self):
        exe = "gust.exe" if sys.platform.startswith("win") else "gust"
        frozen = os.path.join(self.dir, "frozen")
        os.makedirs(frozen)
        self.write(os.path.join("frozen", exe))
        self.write("my_gust")
        mine = os.path.join(self.dir, "./my_gust")
        # run.executable given: it wins, also in a frozen window; the bundled gust only when it is not found
        self.assertEqual(runplan.find_gust("./my_gust", self.dir, frozen, given=True), mine)
        self.assertEqual(runplan.find_gust("./missing", self.dir, frozen, given=True), os.path.join(frozen, exe))
        # run.executable left at its default (gust): the bundled gust first
        self.assertEqual(runplan.find_gust("gust", self.dir, frozen, given=False), os.path.join(frozen, exe))
        self.assertEqual(runplan.find_gust("./my_gust", self.dir, None), mine)
        self.assertIsNone(runplan.find_gust("./missing", self.dir, None))
        self.assertIsNone(runplan.find_gust("no_such_gust_on_path_9f3", self.dir, None))

    def test_command(self):
        self.assertEqual(runplan.command("/g", ["inputs"], 1, "mpiexec -n {ranks}"), ("/g", ["inputs"]))
        self.assertEqual(runplan.command("/g", ["inputs", "a=1"], 4, "mpiexec -n {ranks}"),
                         ("mpiexec", ["-n", "4", "/g", "inputs", "a=1"]))
        launcher = '"C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe" -n {ranks}'
        with mock.patch.object(runplan.sys, "platform", "win32"):
            self.assertEqual(runplan.command("C:\\g.exe", ["inputs"], 2, launcher),
                             ("C:\\Program Files\\Microsoft MPI\\Bin\\mpiexec.exe", ["-n", "2", "C:\\g.exe", "inputs"]))


class FrozenTests(unittest.TestCase):
    def test_checker_reference_in_a_frozen_build(self):
        """The frozen window's deck checker reads the Tests/inputs_reference that gust_gui.spec carries."""
        from gust_input_checker import reference
        here = reference.default_reference_path()
        self.assertTrue(os.path.isfile(here))
        with mock.patch.object(sys, "frozen", True, create=True), \
                mock.patch.object(sys, "_MEIPASS", "/bundle", create=True):
            self.assertEqual(reference.default_reference_path(), os.path.join("/bundle", "Tests", "inputs_reference"))
        with open(os.path.join(HERE, "gust_gui.spec"), encoding="utf-8") as spec:
            self.assertIn('os.path.join(HERE, "..", "Tests", "inputs_reference"), "Tests"', spec.read())


# ---------------------------------------------------------------------------------------------------------------
# The window (PySide6 and pyqtgraph)
# ---------------------------------------------------------------------------------------------------------------
@unittest.skipUnless(HAVE_QT, "PySide6 or pyqtgraph is not installed")
class FormTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self):
        from gust_window.form import CaseForm
        self.form = CaseForm()

    def tearDown(self):
        self.form.deleteLater()

    def keyed_widgets(self):
        return [w for w in self.form.findChildren(QtWidgets.QWidget) if w.property("gust_key")]

    def test_one_widget_per_key(self):
        """One widget per key, named after it, of the key's kind; each patch has one per nest key."""
        names = [n for n, _ in values.field_names()]
        self.assertEqual(sorted(self.form.fields), sorted(names))
        self.assertEqual(sorted(w.objectName() for w in self.keyed_widgets()), sorted(names))
        kinds = {"flag": QtWidgets.QComboBox, "choice": QtWidgets.QComboBox, "any": QtWidgets.QPlainTextEdit}
        for name, key in values.field_names():
            field = self.form.fields[name]
            if name == "domain.center_lat_lon":           # a pair: two lines in one box
                self.assertEqual([type(e) for e in field.edits], [QtWidgets.QLineEdit] * 2)
                continue
            self.assertIsInstance(field.widget, kinds.get(key.kind, QtWidgets.QLineEdit), name)
            self.assertIn(key.range_text(), self.form.fields[name].widget.toolTip(), name)
        self.form.add_nest()
        self.form.add_nest()
        nest = ["nests[%d].%s" % (i, n) for i in range(2) for n in values.nest_spec()]
        self.assertEqual(sorted(w.objectName() for w in self.keyed_widgets()), sorted(names + nest))
        self.form.remove_nest(self.form.nests[0])
        self.assertEqual(sorted(self.form.all_fields()), sorted(names + nest[:len(values.nest_spec())]))

    def round_trip(self, raw):
        from gust_window.window import save_yaml
        self.form.set_raw(copy.deepcopy(raw))
        self.assertEqual(self.form.raw(), raw)
        known = [k for k in SPEC if k in raw]
        self.assertEqual(list(self.form.raw()), known + [k for k in raw if k not in SPEC])   # the schema's order
        for section in known:                             # and inside each section, whatever the form's groups
            if isinstance(SPEC[section], dict) and isinstance(raw[section], dict):
                order = [k for k in SPEC[section] if k in raw[section]]
                self.assertEqual([k for k in self.form.raw()[section] if k in SPEC[section]], order, section)
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "case.yaml")
            save_yaml(path, self.form.raw())
            self.assertEqual(load_yaml(path), raw)

    def test_yaml_round_trip_case_files(self):
        """The repository's case files, read into the form and saved, read back unchanged."""
        self.assertGreaterEqual(len(CASE_FILES), 4)
        for path in CASE_FILES:
            with self.subTest(path=os.path.basename(path)):
                self.round_trip(load_yaml(path))
        self.round_trip(yaml.safe_load(make_gust_case.EXAMPLE))

    def test_yaml_round_trip_every_key(self):
        """Every key given (two patches, gust_inputs with a null), and keys the schema does not know."""
        raw = every_key_case()
        validate(raw)                                     # every sample value is one the schema accepts
        self.round_trip(raw)
        raw["domain"]["not_a_key"] = 4
        raw["nests"][0]["not_a_key"] = [1, 2]
        raw["another_section"] = {"a": 1}
        self.round_trip(raw)
        raw = yaml.safe_load(make_gust_case.EXAMPLE)
        raw["case"]["name"] = True                        # name: yes, refused: the round trip must not change that
        raw["mesh"]["layers"] = "40"
        self.round_trip(raw)
        self.assertRaisesRegex(CaseError, "case.name", validate, self.form.raw())

    def test_placeholders_show_the_format_and_default(self):
        fields = self.form.fields
        hint = lambda name: fields[name].widget.placeholderText()          # noqa: E731
        self.assertEqual(hint("domain.center"), "easting, northing")
        self.assertEqual([e.placeholderText() for e in fields["domain.center_lat_lon"].edits],
                         ["latitude (-90 to 90)", "longitude (-180 to 180)"])
        self.assertEqual(hint("mesh.cell_size"), "number > 0")
        self.assertEqual(hint("terrain.margin"), "number ≥ 0, or auto; default auto")
        self.assertEqual(hint("mesh.layers"), "whole number ≥ 1")
        self.assertEqual(hint("domain.utm_zone"), "1 to 60; empty: from the longitude")
        self.assertEqual(hint("domain.half_width"), "one number, or east-west, north-south")
        self.assertEqual(hint("atmosphere.reference_temperature"), "number > 0; default 300.0")
        self.assertEqual(hint("solver.initial_blend"), "2 numbers > 0; default 1.0, 0.5")
        self.assertEqual(hint("surface.roughness_map.classes"), "YAML (see the tooltip); default esa_cci")
        self.assertEqual(fields["terrain.format"].widget.lineEdit().placeholderText(),
                         "one of: auto, geotiff, esri_ascii, xyz; default auto")
        for name, field in fields.items():
            if field.edits is not None:
                continue
            text = field.widget.lineEdit().placeholderText() if isinstance(field.widget, QtWidgets.QComboBox) \
                else field.widget.placeholderText()
            self.assertTrue(text and "not given" not in text, name)
        default = hint("solver.initial_blend").split("; default ")[1]
        self.assertEqual(values.from_text("solver.initial_blend", SPEC["solver"]["initial_blend"], default),
                         (True, [1.0, 0.5]))       # the default as shown is one the field takes

    def shown(self, name):
        field = self.form.all_fields()[name]
        for form, index, f in self.form.rows.rows:
            if f is field:
                return form.isRowVisible(index)
        return True

    def test_basic_mode_shows_what_a_case_needs(self):
        from gust_window import layout
        self.form.set_raw({})
        for name in self.form.fields:
            self.assertEqual(self.shown(name), name in layout.BASIC, name)
        visible = {self.form.tabText(i) for i in range(self.form.count()) if self.form.isTabVisible(i)}
        self.assertEqual(visible, {"case", "domain", "terrain", "mesh", "wind", "surface", "run"})
        self.form.set_advanced(True)
        self.assertTrue(all(self.shown(n) for n in self.form.fields))
        self.assertTrue(all(self.form.isTabVisible(i) for i in range(self.form.count())))
        self.form.set_advanced(False)
        self.assertFalse(self.shown("solver.cfl_max"))

    def test_basic_mode_never_hides_a_value_or_a_refusal(self):
        raw = yaml.safe_load(make_gust_case.EXAMPLE)
        raw.setdefault("solver", {})["cfl_max"] = 8.0
        self.form.set_raw(raw)
        self.assertTrue(self.shown("solver.cfl_max"))
        self.assertFalse(self.shown("solver.sweep"))
        self.assertTrue(any(self.form.tabText(i) == "solver" and self.form.isTabVisible(i)
                            for i in range(self.form.count())))
        self.form.mark(["model.tke_floor = -1: must be 0, or >= 2.22e-16"])
        self.assertTrue(self.shown("model.tke_floor"))
        self.form.show_field("output.prefix")             # a finding on a hidden field shows every setting
        self.assertTrue(self.form.advanced)

    def test_latitude_and_longitude_boxes(self):
        field = self.form.fields["domain.center_lat_lon"]
        field.set_value([57.18855, -7.37267])
        self.assertEqual([e.text() for e in field.edits], ["57.18855", "-7.37267"])
        field.edits[0].setText("45.5")
        self.assertEqual(field.value(), (True, [45.5, -7.37267]))
        field.edits[0].setText("")
        field.edits[1].setText("")
        self.assertEqual(field.value(), (False, None))
        field.edits[0].setText("45.5")                     # one box only: a list of one, refused by the builder
        self.assertEqual(field.value(), (True, [45.5]))
        field.set_value("north")                           # not a pair: kept whole
        self.assertEqual(field.value(), (True, "north"))
        field.set_value([1.0, 2.0, 3.0])                   # three entries: kept whole, not cut to two
        self.assertEqual([e.text() for e in field.edits], ["1.0, 2.0, 3.0", ""])
        self.assertEqual(field.value(), (True, [1.0, 2.0, 3.0]))

    def test_paths_have_a_browse_button(self):
        from gust_window.form import _relative
        for name, key in values.field_names():
            self.assertEqual(hasattr(self.form.fields[name], "browse_button"), key.kind == "path", name)
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(_relative(os.path.join(tmp, "terrain", "dem.asc"), tmp), "terrain/dem.asc")
            self.assertEqual(_relative(os.path.dirname(tmp), tmp), os.path.dirname(tmp))   # outside: absolute
        field = self.form.fields["case.output_dir"]
        with mock.patch.object(QtWidgets.QFileDialog, "getExistingDirectory", return_value=os.getcwd()):
            field.browse_button.click()
        self.assertEqual(field.text(), ".")

    def test_refused_key_marks_its_field(self):
        """The field the refusal names is marked, and no other."""
        raw = yaml.safe_load(make_gust_case.EXAMPLE)
        raw["mesh"]["layers"] = 0
        raw["nests"] = [{"name": "a", "half_width": [200.0]}, {"name": "b"}]
        self.form.set_raw(raw)
        with self.assertRaises(CaseError) as err:
            validate(self.form.raw())
        marked = self.form.mark(err.exception.problems)
        self.assertEqual(sorted(n for _, n in marked), ["mesh.layers", "nests[1].half_width"])
        self.assertEqual(sorted(n for n, f in self.form.all_fields().items() if f.marked),
                         ["mesh.layers", "nests[1].half_width"])
        self.assertIn("mesh.layers = 0", self.form.fields["mesh.layers"].widget.toolTip())
        self.form.mark([])
        self.assertFalse(any(f.marked for f in self.form.all_fields().values()))

    def test_refused_key_marks_text_that_is_not_yaml(self):
        """A field whose text is not YAML is refused by name and marked."""
        self.form.set_raw(yaml.safe_load(make_gust_case.EXAMPLE))
        self.form.fields["domain.half_width"].set_text("[2000,")
        with self.assertRaises(CaseError) as err:
            self.form.raw()
        self.form.mark(err.exception.problems)
        self.assertEqual([n for n, f in self.form.all_fields().items() if f.marked], ["domain.half_width"])


@unittest.skipUnless(HAVE_QT, "PySide6 or pyqtgraph is not installed")
class WindowTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def window(self, raw):
        from gust_window.window import MainWindow, save_yaml
        path = os.path.join(self.dir, "case.yaml")
        save_yaml(path, raw)
        w = MainWindow(path)
        self.addCleanup(w.deleteLater)
        return w

    def test_refused_key_marks_the_field_from_check(self):
        """Through the window: a refusal of the builder (two keys related) marks the key it names."""
        raw = yaml.safe_load(make_gust_case.EXAMPLE)
        raw["nests"] = [{"name": "a", "center": [0.0, 0.0], "half_width": 300.0, "blend_width": 100.0}]
        w = self.window(raw)
        self.assertIsNone(w.check())
        self.assertEqual([n for n, f in w.form.all_fields().items() if f.marked], ["nests[0].blend_width"])
        self.assertEqual(w.findings.item(0, 2).text(), "nests[0].blend_width")
        self.assertIs(w.tabs.currentWidget(), w.findings)

    def test_check_draws_and_write_fills_the_findings(self):
        w = self.window(yaml.safe_load(make_gust_case.EXAMPLE))
        result = w.check()
        self.assertIsNotNone(result)
        self.assertEqual((result.plan.nx, result.plan.ny, result.plan.layers), (80, 80, 32))
        self.assertIn("80 x 80 x 32 cells", w.plan_view.caption.text())
        w.form.fields["solver.cfl_max"].set_text("400")   # the deck checker has something to say about this
        written = w.write_case()
        self.assertTrue(os.path.isfile(os.path.join(self.dir, "flat_west", "inputs")))
        self.assertEqual(load_yaml(w.path)["solver"]["cfl_max"], 400)        # Write saved the case file first
        found = [(w.findings.item(i, 0).text(), w.findings.item(i, 1).text()) for i in range(w.findings.rowCount())]
        expect = [(f.severity.name, rel) for rel, report in written.reports for f in report.sorted()]
        self.assertEqual(found, expect)

    def test_download_worker_ends_on_any_error(self):
        from gust_window.window import DownloadWorker
        worker = DownloadWorker((57.0, 57.1, -7.5, -7.4), os.path.join(self.dir, "d.asc"), self.dir, "zero")
        got = []
        worker.done.connect(lambda result, error: got.append((result, error)))
        with mock.patch("gust_window.window.dem.download", side_effect=RuntimeError("socket closed")):
            worker.run()
        self.assertEqual(len(got), 1)
        self.assertIsNone(got[0][0])
        self.assertEqual(got[0][1].problems, ["terrain download: RuntimeError: socket closed"])

    def test_terrain_image_keeps_the_east_column(self):
        # here x0 + 35 * size falls an ulp short of x1: the last bin edge must reach x1 or that column is lost
        import types
        from gust_window.mapview import terrain_image
        x0, x1 = -23990.255, -7191.023
        X, Y = np.meshgrid(np.linspace(x0, x1, 101), np.linspace(-100.0, -100.0 + 0.3 * (x1 - x0), 31))
        z = np.where(X == x1, 1.0, 0.0)
        image, rect = terrain_image(types.SimpleNamespace(x=X.ravel(), y=Y.ravel(), z=z.ravel()), 35)
        self.assertEqual(image.shape[1], 35)
        self.assertGreater(np.nanmax(image[:, -1]), 0.0)
        self.assertGreaterEqual(rect[2], x1 - x0)          # the image is drawn over every point

    def test_layer_faces(self):
        from gust_window.mapview import layer_faces

        class P:
            layers, top, dz0, stretch_ratio = 4, 15.0, 1.0, 2.0
        self.assertEqual(list(layer_faces(P)), [0.0, 1.0, 3.0, 7.0, 15.0])
        P.dz0 = None
        self.assertEqual(list(layer_faces(P)), [0.0, 3.75, 7.5, 11.25, 15.0])


FAKE_GUST = '''#!{python}
"""A stand-in for gust: prints an iteration table, writes what run.sh expects, and records its arguments."""
import os, sys
here = os.path.basename(os.getcwd())
with open("args.txt", "w") as out:
    out.write(" ".join(sys.argv[1:]))
if here == "column":
    open("column_guess.dat", "w").write("table")
    print("  forcing.ug = 7.5")
    print("  forcing.vg = -1.25")
print("  Iter        CFL  Line solver      u-res      v-res      w-res    tke-res  theta-res      p-res          r")
for i in range(1, 4):
    print("%6d  1.000e+00  lusgs        %9.3e  1.0e-03  1.0e-03  1.0e-03  1.0e-03  1.0e-03  1.000e+00" % (i, 10.0 ** -i))
sys.stdout.flush()
if os.environ.get("FAKE_GUST_SLEEP") == here:
    import time
    time.sleep(60)
if here not in ("column",) and not here.startswith("nest_"):
    os.makedirs("chk00003", exist_ok=True)
sys.exit(3 if os.environ.get("FAKE_GUST_FAIL") == here else 0)
'''


@unittest.skipUnless(HAVE_QT, "PySide6 or pyqtgraph is not installed")
@unittest.skipIf(sys.platform.startswith("win"), "the stand-in gust is a script with a #! line")
class RunPanelTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name
        for sub in ("column", "nest_a"):
            os.makedirs(os.path.join(self.dir, sub))
        for rel in ("inputs", "column/inputs", "nest_a/inputs", resume.DATA):   # the resume record's files
            with open(os.path.join(self.dir, rel), "w", encoding="utf-8") as out:
                out.write("# deck\n")
        self.gust = os.path.join(self.dir, "fake_gust")
        with open(self.gust, "w", encoding="utf-8") as out:
            out.write(FAKE_GUST.format(python=sys.executable))
        os.chmod(self.gust, os.stat(self.gust).st_mode | stat.S_IEXEC)
        self.runs = [{"dir": ".", "label": "wind from 240 degrees", "column": True, "wind_from_column": True,
                      "patches": ["a"], "chk": "chk"}]

    def tearDown(self):
        os.environ.pop("FAKE_GUST_FAIL", None)
        os.environ.pop("FAKE_GUST_SLEEP", None)
        self._tmp.cleanup()

    def run_panel(self, stop_after_ms=None):
        from gust_window.runner import RunPanel
        panel = RunPanel()
        self.addCleanup(panel.deleteLater)
        outcome = []
        loop = QtCore.QEventLoop()
        panel.finished.connect(lambda ok, message: (outcome.append((ok, message)), loop.quit()))
        QtCore.QTimer.singleShot(60000, loop.quit)
        if stop_after_ms is not None:
            QtCore.QTimer.singleShot(stop_after_ms, panel.stop)
        panel.start(self.dir, self.runs, self.gust, 1, "mpiexec -n {ranks}")
        if not outcome:
            loop.exec()
        return panel, outcome

    def read(self, rel):
        with open(os.path.join(self.dir, rel), encoding="utf-8") as handle:
            return handle.read()

    def test_steps_run_in_order_with_what_they_need(self):
        panel, outcome = self.run_panel()
        self.assertEqual(outcome, [(True, "All steps ended with status 0.")])
        self.assertEqual(self.read("column/args.txt"), "inputs")
        self.assertEqual(self.read("args.txt"), "inputs forcing.ug=7.5 forcing.vg=-1.25")
        self.assertEqual(self.read("nest_a/args.txt"), "inputs bc.nest_parent=../chk00003 forcing.ug=7.5 forcing.vg=-1.25")
        self.assertIn("lusgs", self.read("nest_a/log.txt"))
        xs, ys = panel._data["u"]                         # the last step's residuals, as plotted
        self.assertEqual(xs, [1, 2, 3])
        self.assertEqual(ys, [0.1, 0.01, 0.001])

    def test_stop_ends_the_step_and_says_so(self):
        os.environ["FAKE_GUST_SLEEP"] = os.path.basename(self.dir)
        start = time.time()
        panel, outcome = self.run_panel(stop_after_ms=1500)
        self.assertLess(time.time() - start, 5.0)       # terminate, not the kill after STOP_GRACE_MS
        self.assertEqual(outcome, [(False, "wind from 240 degrees stopped (Stop)")])
        self.assertFalse(panel.running())
        self.assertFalse(os.path.exists(os.path.join(self.dir, "nest_a", "args.txt")))

    def test_a_step_that_cannot_start_names_the_program(self):
        self.gust = os.path.join(self.dir, "no_such_gust")
        panel, outcome = self.run_panel()
        self.assertEqual(outcome, [(False, "wind from 240 degrees: column did not start: %s" % self.gust)])

    def test_a_failed_step_stops_the_run(self):
        os.environ["FAKE_GUST_FAIL"] = os.path.basename(self.dir)
        panel, outcome = self.run_panel()
        self.assertEqual(len(outcome), 1)
        self.assertFalse(outcome[0][0])
        self.assertIn("ended with status 3", outcome[0][1])
        self.assertFalse(os.path.exists(os.path.join(self.dir, "nest_a", "args.txt")))


@unittest.skipUnless(HAVE_QT, "PySide6 or pyqtgraph is not installed")
class LibraryTabTests(unittest.TestCase):
    """the Library tab on a hand-written two-sector library (tools/test_gust_ensemble.py)."""

    @classmethod
    def setUpClass(cls):
        cls.app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])

    def setUp(self):
        import test_gust_ensemble as ens
        self.ens = ens
        self._tmp = tempfile.TemporaryDirectory()
        self.case = os.path.join(self._tmp.name, "case")
        ens.write_library(self.case, ens.domain(center_lat_lon=[43.40, -113.03]))
        from gust_window.librarytab import LibraryTab
        self.tab = LibraryTab()
        self.addCleanup(self.tab.deleteLater)

    def tearDown(self):
        self._tmp.cleanup()

    def test_tab_in_the_window(self):
        from gust_window.window import MainWindow
        w = MainWindow()
        self.addCleanup(w.deleteLater)
        self.assertEqual([w.tabs.tabText(i) for i in range(w.tabs.count())],
                         ["Map", "Findings", "Summary", "Run", "Library"])
        self.assertIs(w.tabs.widget(4), w.library)

    def test_one_wind_shows_and_exports(self):
        tab = self.tab
        self.assertFalse(tab.show_button.isEnabled())
        tab.set_folder(self.case)
        self.assertEqual([tab.height.itemText(i) for i in range(tab.height.count())], ["10", "40"])
        tab.direction.setValue(90.0)
        tab.speed.setValue(10.0)
        fields = tab.compute()
        self.assertEqual(list(fields), ["speed", "direction", "u", "v", "w", "k"])
        w = tab.lib.wind(90.0, 10.0, 10.0)
        np.testing.assert_allclose(fields["speed"], np.hypot(w["u"], w["v"]), rtol=1e-14)
        self.assertEqual(fields["direction"][0, 0], 90.0)
        self.assertEqual(fields["direction"][self.ens.REVERSED], 270.0)
        self.assertTrue(tab.grids_button.isEnabled())
        out = os.path.join(self._tmp.name, "out")
        tab.units.setCurrentText("kph")
        files = tab.export_grids(out)
        self.assertEqual(sorted(os.path.basename(f) for f in files),
                         ["twosector_90deg_10mps_10m_agl_%s" % s for s in ("ang.asc", "ang.prj", "vel.asc", "vel.prj")])
        from gust_library import grids
        _, vel, _, _ = grids.read_grid(os.path.join(out, "twosector_90deg_10mps_10m_agl_vel.asc"))
        np.testing.assert_allclose(vel, 3.6 * fields["speed"], rtol=1e-8)

    def test_distribution_and_refusal_in_the_caption(self):
        tab = self.tab
        tab.set_folder(self.case)
        tab.mode.setCurrentText("Distribution")
        self.assertTrue(tab.members.isVisibleTo(tab) and not tab.direction.isVisibleTo(tab))
        tab.members.setPlainText("0 6 1\n90 10 3\n")
        tab.thresholds.setText("6.5")
        fields = tab.compute()
        self.assertEqual(list(fields), ["mean_speed", "std_speed", "mean_u", "mean_v", "mean_direction", "p10", "p50", "p90",
                                        "pexceed_6p5", "preversal"])
        expect = self.ens.ensemble.ensemble(tab.lib, [(0.0, 6.0, 1.0), (90.0, 10.0, 3.0)], 10.0, (10, 50, 90), (6.5,))
        for name in fields:
            np.testing.assert_array_equal(fields[name], expect[name])
        points = os.path.join(self._tmp.name, "pts.txt")
        with open(points, "w") as fh:
            fh.write("12.5 12.5 mast\n")
        tab.export_points(points, self._tmp.name)
        self.assertTrue(os.path.isfile(os.path.join(self._tmp.name, "twosector_ensemble_10m_agl_points.csv")))
        tab.members.setPlainText("0 6 -1\n")
        self.assertIsNone(tab.compute())
        self.assertIn("weights must be 0 or more", tab.caption.text())
        self.assertFalse(tab.grids_button.isEnabled())

    def test_folder_without_a_library(self):
        self.tab.set_folder(self._tmp.name)
        self.assertIn("build the library first", self.tab.caption.text())
        self.assertFalse(self.tab.show_button.isEnabled())

    def test_default_units_are_mps(self):
        self.assertEqual(self.tab.units.currentText(), "mps")

    def test_export_is_the_wind_shown(self):
        """Boxes changed after Show: the files are named and sampled from the wind shown, not the boxes."""
        tab = self.tab
        tab.set_folder(self.case)
        tab.direction.setValue(90.0)
        tab.speed.setValue(10.0)
        fields = tab.compute()
        tab.direction.setValue(0.0)
        tab.speed.setValue(5.0)
        out = os.path.join(self._tmp.name, "out")
        tab.export_grids(out)
        from gust_library import grids
        _, vel, _, _ = grids.read_grid(os.path.join(out, "twosector_90deg_10mps_10m_agl_vel.asc"))
        np.testing.assert_allclose(vel, fields["speed"], rtol=1e-8)
        points = os.path.join(self._tmp.name, "pts.txt")
        with open(points, "w") as fh:
            fh.write("12.5 12.5 mast\n")
        tab.export_points(points, out)
        with open(os.path.join(out, "twosector_90deg_10mps_10m_agl_points.csv")) as fh:
            row = dict(zip(*[line.strip().split(",") for line in fh]))
        w = tab.lib.wind(90.0, 10.0, 10.0)
        i, j = list(self.ens.X).index(12.5), list(self.ens.Y).index(12.5)
        self.assertAlmostEqual(float(row["speed"]), float(np.hypot(w["u"], w["v"])[i, j]), places=9)

    def test_failed_load_drops_the_library_loaded_before(self):
        tab = self.tab
        tab.set_folder(self.case)
        self.assertIsNotNone(tab.compute())
        tab.set_folder(self._tmp.name)
        self.assertIsNone(tab.lib)
        self.assertFalse(tab.show_button.isEnabled())
        self.assertFalse(tab.grids_button.isEnabled())


if __name__ == "__main__":
    unittest.main()
