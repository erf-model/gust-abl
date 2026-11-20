"""Resumable sweeps: run.sh and the set-up window skip a direction that holds its answer (gust_case/resume.py).

    python3 -m unittest -v test_gust_resume        (from tools/)

A fake gust stands in for the solver: it writes the log lines GUST writes (the verdict, then ``gust: wrote``) and the
plotfile's Header, and appends its folder to a call list, so a test sees which runs run.sh started.
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from gust_case import build, resume, validate  # noqa: E402
from test_gust_case import base_case, changed, write_hill  # noqa: E402

FAKE_GUST = r"""#!/bin/bash
# The verdict comes from $FAKE_VERDICT (steady by default); a column run writes the files run.sh reads.
echo "$PWD" >> "$FAKE_CALLS"
if [ "$(basename "$PWD")" = column ]; then
  [ -z "${FAKE_COLUMN_FAILS:-}" ] || exit 1
  echo x > column_guess.dat; echo "  forcing.ug = 7.5"; echo "  forcing.vg = -1.25"; exit 0
fi
prefix=$(LC_ALL=C awk '{sub(/\r$/, "")} $1 == "io.plot_file" {v = $3} END {print (v == "" ? "plt" : v)}' inputs)
echo "level 0:      1  1.000e+00  lusgs  1 1 1 1 1 1  1.000e+00"
echo "gust: ${FAKE_VERDICT:-steady} after 12 iterations (12 over the run)"
mkdir -p "${prefix}00012" && touch "${prefix}00012/Header"
echo "gust: wrote ${prefix}00012"
rm -rf chk00012 && mkdir chk00012      # AMReX moves an existing checkpoint aside: a new folder
"""


@unittest.skipIf(os.name == "nt", "run.sh is a bash script")
class RunScriptResumeTests(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self._tmp.cleanup)
        self.dir = self._tmp.name
        self.fake = os.path.join(self.dir, "fake_gust")
        with open(self.fake, "w") as out:
            out.write(FAKE_GUST)
        os.chmod(self.fake, 0o755)
        self.calls = os.path.join(self.dir, "calls.txt")

    def write(self, raw, overwrite=False):
        case, given = validate(raw)
        self.result = build(case, given, self.dir, overwrite=overwrite, check=False)
        return self.result

    def sweep(self, **changes):
        return changed(base_case(), wind__direction=None, wind__directions=[260.0, 270.0, 280.0], **changes)

    def run_sh(self, *args, verdict="steady", resume=None, column_fails=False):
        """Run run.sh with the fake gust; returns the folders it ran gust in (relative to the case folder)."""
        if os.path.exists(self.calls):
            os.remove(self.calls)
        env = dict(os.environ, GUST=self.fake, MPI="", FAKE_CALLS=self.calls, FAKE_VERDICT=verdict)
        env.pop("RESUME", None)
        env.pop("FAKE_COLUMN_FAILS", None)
        if column_fails:
            env["FAKE_COLUMN_FAILS"] = "1"
        if resume is not None:
            env["RESUME"] = resume
        run = subprocess.run(["bash", os.path.join(self.result.case_dir, "run.sh")] + list(args), env=env,
                             capture_output=True, text=True, errors="replace")
        self.assertEqual(run.returncode != 0, column_fails, run.stderr)
        self.stdout = run.stdout
        if not os.path.exists(self.calls):
            return []
        with open(self.calls) as handle:
            real = os.path.realpath(self.result.case_dir)
            return [os.path.relpath(os.path.realpath(line.strip()), real) for line in handle]

    def python_finished(self):
        """resume.finished of every direction, as the window sees them."""
        return [resume.finished(os.path.join(self.result.case_dir, r["dir"]), r) for r in self.result.runs]

    def test_finished_directions_are_skipped_and_the_others_run(self):
        self.write(self.sweep())
        self.assertEqual(self.run_sh(), ["dir_260", "dir_270", "dir_280"])
        self.assertEqual(self.python_finished(), [True, True, True])
        self.assertEqual(self.run_sh(), [])                             # every direction holds its answer
        self.assertEqual(self.stdout.count("finished earlier"), 3)
        # an interrupted sweep: 270's plotfile is gone, 280 stopped before its verdict
        import shutil
        shutil.rmtree(os.path.join(self.result.case_dir, "dir_270", "plt00012"))
        with open(os.path.join(self.result.case_dir, "dir_280", "log.txt"), "w") as out:
            out.write("level 0:      1  1.000e+00  lusgs  1 1 1 1 1 1  1.000e+00\n")
        self.assertEqual(self.python_finished(), [True, False, False])
        self.assertEqual(self.run_sh(), ["dir_270", "dir_280"])

    def test_a_run_that_is_not_steady_runs_again(self):
        self.write(self.sweep())
        self.run_sh(verdict="not steady")
        self.assertEqual(self.python_finished(), [False, False, False])
        self.assertEqual(self.run_sh(verdict="diverged"), ["dir_260", "dir_270", "dir_280"])
        self.assertEqual(self.python_finished(), [False, False, False])
        self.assertEqual(self.run_sh(), ["dir_260", "dir_270", "dir_280"])

    def test_a_later_verdict_decides(self):
        """A steady verdict followed by a later run's not-steady verdict in the same log is not finished, and a
        plotfile written before the steady verdict does not count."""
        self.write(self.sweep())
        self.run_sh()
        folder = os.path.join(self.result.case_dir, "dir_260")
        with open(os.path.join(folder, "log.txt"), "a") as out:
            out.write("gust: not steady after 800 iterations (812 over the run)\n")
        self.assertFalse(resume.finished(folder, self.result.runs[0]))
        with open(os.path.join(folder, "log.txt"), "w") as out:
            out.write("gust: wrote plt00012\ngust: steady after 12 iterations (12 over the run)\n")
        self.assertFalse(resume.finished(folder, self.result.runs[0]))
        self.assertEqual(self.run_sh(), ["dir_260"])
        # two steady verdicts: the plotfile written after the first does not count for the second
        with open(os.path.join(folder, "log.txt"), "a") as out:
            out.write("gust: steady after 3 iterations (15 over the run)\n")
        self.assertFalse(resume.finished(folder, self.result.runs[0]))
        self.assertEqual(self.run_sh(), ["dir_260"])

    def test_other_arguments_or_decks_run_again(self):
        self.write(self.sweep())
        self.run_sh("steady.cfl_max=8")
        self.assertEqual(self.python_finished(), [False, False, False])                # the window gives no arguments
        self.assertEqual([resume.finished(os.path.join(self.result.case_dir, r["dir"]), r, ["steady.cfl_max=8"])
                          for r in self.result.runs], [True, True, True])
        self.assertEqual(self.run_sh("steady.cfl_max=8"), [])
        self.assertEqual(self.run_sh(), ["dir_260", "dir_270", "dir_280"])     # without the argument: another run
        self.write(self.sweep(), overwrite=True)                                 # written again, the same decks
        self.assertEqual(self.run_sh(), [])
        self.assertEqual(self.python_finished(), [True, True, True])
        self.write(self.sweep(solver__max_iterations=900), overwrite=True)       # another deck
        self.assertEqual(self.python_finished(), [False, False, False])
        self.assertEqual(self.run_sh(), ["dir_260", "dir_270", "dir_280"])

    def test_resume_off_runs_every_direction(self):
        self.write(self.sweep())
        self.run_sh()
        self.assertEqual(self.run_sh(resume="0"), ["dir_260", "dir_270", "dir_280"])
        self.write(self.sweep(run__resume=False), overwrite=True)
        self.run_sh()
        self.assertEqual(self.run_sh(), ["dir_260", "dir_270", "dir_280"])     # run.resume = false: the old behaviour
        self.assertEqual(self.run_sh(resume="1"), [])

    def test_patches_and_the_column_run(self):
        """A direction with a patch is finished only when the patch is; the column run is part of the direction."""
        raw = changed(base_case(), wind__direction=None, wind__directions=[0.0], wind__drive="geostrophic_column",
                      atmosphere={"latitude": 45.0},
                      nests=[{"name": "mid", "center": [0.0, 0.0], "half_width": [750.0]}])
        self.write(raw)
        self.assertEqual(self.run_sh(), ["dir_000/column", "dir_000", "dir_000/nest_mid"])
        self.assertEqual(self.run_sh(), [])
        os.remove(os.path.join(self.result.case_dir, "dir_000", "nest_mid", "plt00012", "Header"))
        self.assertEqual(self.python_finished(), [False])
        self.assertEqual(self.run_sh(), ["dir_000/column", "dir_000", "dir_000/nest_mid"])
        with open(os.path.join(self.result.case_dir, "dir_000", "column", "inputs"), "a") as out:
            out.write("# a changed column deck\n")
        self.assertEqual(self.python_finished(), [False])
        # stopped part-way (the column run fails): the old steady logs are gone with the new record
        self.assertEqual(self.run_sh(column_fails=True), ["dir_000/column"])
        self.assertEqual(self.python_finished(), [False])
        self.assertEqual(self.run_sh(), ["dir_000/column", "dir_000", "dir_000/nest_mid"])
        self.assertEqual(self.run_sh(), [])

    def test_another_data_file_runs_again(self):
        """A case written again with another terrain file and the same decks runs every direction again."""
        hill = os.path.join(self.dir, "hill.asc")
        center = write_hill(hill)
        raw = self.sweep(domain={"center": list(center), "half_width": [2000.0], "top": 1000.0},
                         terrain={"file": hill})
        self.write(raw)
        case_dir = self.result.case_dir
        before = [open(os.path.join(case_dir, d, "inputs"), "rb").read() for d in ("dir_260", "dir_270", "dir_280")]
        self.run_sh()
        self.assertEqual(self.run_sh(), [])
        write_hill(hill, x0=500030.0)                      # the hill 30 m further east: the same relief
        self.write(raw, overwrite=True)
        after = [open(os.path.join(case_dir, d, "inputs"), "rb").read() for d in ("dir_260", "dir_270", "dir_280")]
        self.assertEqual(before, after)                                          # the decks alone do not show it
        self.assertEqual(self.python_finished(), [False, False, False])
        self.assertEqual(self.run_sh(), ["dir_260", "dir_270", "dir_280"])
        self.assertEqual(self.python_finished(), [True, True, True])

    def test_plot_file_of_the_deck(self):
        """The test reads the plotfile prefix from the deck: a plotfile of the default name does not count."""
        self.write(self.sweep(output__prefix="a"))
        self.assertEqual(self.run_sh(), ["dir_260", "dir_270", "dir_280"])
        folder = os.path.join(self.result.case_dir, "dir_260")
        self.assertEqual(resume.plot_prefix(os.path.join(folder, "inputs")), "a_plt")
        self.assertEqual(self.run_sh(), [])
        with open(os.path.join(folder, "log.txt"), "a") as out:                 # a later plotfile of another name
            out.write("gust: steady after 3 iterations (15 over the run)\ngust: wrote plt00015\n")
        os.makedirs(os.path.join(folder, "plt00015"))
        open(os.path.join(folder, "plt00015", "Header"), "w").close()
        self.assertEqual(self.python_finished(), [False, True, True])
        self.assertEqual(self.run_sh(), ["dir_260"])

    def test_windows_logs_and_decks(self):
        """CRLF line ends (gust.exe through the window) and a deck comment that is not UTF-8: run.sh and the window
        agree that the direction is finished. The deck names its plotfile prefix, so the prefix is read through
        the CRLF too."""
        self.write(self.sweep(output__prefix="a"))
        self.run_sh()
        folder = os.path.join(self.result.case_dir, "dir_260")
        for name in ("log.txt", "inputs"):
            path = os.path.join(folder, name)
            with open(path, "rb") as handle:
                data = handle.read()
            with open(path, "wb") as out:
                out.write(data.replace(b"\n", b"\r\n") + (b"# 90\xb0 sector\r\n" if name == "inputs" else b""))
        self.assertFalse(resume.finished(folder, self.result.runs[0]))           # the deck changed: runs again
        self.assertEqual(self.run_sh(), ["dir_260"])
        path = os.path.join(folder, "log.txt")
        with open(path, "rb") as handle:
            data = handle.read()
        with open(path, "wb") as out:                                             # a log line that is not UTF-8
            out.write(b"gust: terrain from \xb0 file\r\n" + data.replace(b"\n", b"\r\n"))
        self.assertEqual(self.python_finished(), [True, True, True])
        self.assertEqual(self.run_sh(), [])

    def test_joint_and_unsteady_verdicts(self):
        """A joint run's verdict decides like the run's own; an unsteady stage after a steady verdict does not count."""
        self.write(self.sweep())
        self.run_sh()
        folder = os.path.join(self.result.case_dir, "dir_260")
        log = os.path.join(folder, "log.txt")
        for tail, done in (("gust: joint run: outer not steady after 40 joint iterations (52 over the run)\n"
                            "gust: wrote plt00012\n", False),
                           ("gust: joint run: outer steady after 40 joint iterations (52 over the run)\n"
                            "gust: wrote plt00012\n", True),
                           ("gust: unsteady run ended after 40 iterations (52 over the run)\ngust: wrote plt00012\n",
                            False)):
            with open(log, "w") as out:
                out.write("gust: steady after 12 iterations (12 over the run)\ngust: wrote plt00012\n" + tail)
            self.assertEqual(resume.finished(folder, self.result.runs[0]), done, tail)
            self.assertEqual(self.run_sh(), [] if done else ["dir_260"], tail)

    def test_missing_deck_stops_before_the_logs_go(self):
        self.write(self.sweep())
        self.run_sh()
        os.remove(os.path.join(self.result.case_dir, "dir_270", "inputs"))
        env = dict(os.environ, GUST=self.fake, MPI="", FAKE_CALLS=self.calls, RESUME="0")
        run = subprocess.run(["bash", os.path.join(self.result.case_dir, "run.sh")], env=env,
                             capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("wind from 270 degrees: inputs is missing; write the case again", run.stderr)
        self.assertTrue(os.path.isfile(os.path.join(self.result.case_dir, "dir_270", "log.txt")))

    def test_resume_takes_0_or_1(self):
        self.write(self.sweep())
        env = dict(os.environ, GUST=self.fake, MPI="", FAKE_CALLS=self.calls, RESUME="yes")
        run = subprocess.run(["bash", os.path.join(self.result.case_dir, "run.sh")], env=env,
                             capture_output=True, text=True)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("RESUME must be 0 or 1, not 'yes'", run.stderr)
        self.assertFalse(os.path.exists(self.calls))

    def test_single_direction_case(self):
        self.write(base_case())
        self.assertEqual(self.run_sh(), ["."])
        self.assertEqual(self.python_finished(), [True])
        self.assertEqual(self.run_sh(), [])


class StartTests(unittest.TestCase):
    def test_start_removes_old_logs_then_writes_the_record(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "nest_a"))
            for rel, text in (("inputs", "io.plot_file = p\n"), ("nest_a/inputs", "x = 1\n"), (resume.DATA, "d\n"),
                              ("log.txt", "gust: steady after 1 iterations\ngust: wrote p00001\n"),
                              ("nest_a/log.txt", "old\n")):
                with open(os.path.join(d, rel), "w") as out:
                    out.write(text)
            run = {"dir": ".", "column": False, "patches": ["a"]}
            resume.start(d, run, ["k=v"])
            self.assertFalse(os.path.exists(os.path.join(d, "log.txt")))
            self.assertFalse(os.path.exists(os.path.join(d, "nest_a", "log.txt")))
            with open(os.path.join(d, resume.RECORD), "rb") as handle:
                self.assertEqual(handle.read(), b"io.plot_file = p\nx = 1\nd\n# arguments: k=v\n")
            self.assertEqual(resume.plot_prefix(os.path.join(d, "inputs")), "p")
            self.assertFalse(resume.finished(d, run, ["k=v"]))          # no log: not finished



class PathTests(unittest.TestCase):
    def test_record_files_use_slashes(self):
        """run.sh names the decks with / on every system (a case written on Windows runs under bash)."""
        import ntpath
        from unittest import mock
        run = {"dir": "dir_270", "column": True, "patches": ["mid"]}
        with mock.patch.object(os, "path", ntpath):
            files = resume.record_files(run)
        self.assertEqual(files, ["inputs", "column/inputs", "nest_mid/inputs", "../" + resume.DATA])
        self.assertEqual(resume.record_files(dict(run, dir="."))[-1], resume.DATA)


if __name__ == "__main__":
    unittest.main()
