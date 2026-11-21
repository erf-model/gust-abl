"""The steps of a case's run, in run.sh's order, and where gust is (no Qt).

For each wind direction: the column run (column/, when the wind is geostrophic), then the terrain run, then each
nested patch from the newest checkpoint the terrain run wrote. :func:`steps` is a generator: the window starts
each step it yields and asks for the next when the step has ended with status 0, so a step can use what the
step before it wrote (the column's geostrophic wind, the parent's checkpoint), as run.sh does.
"""
import os
import re
import shlex
import shutil
import sys

from gust_case import library, resume
from gust_case.builder import COLUMN_TABLE

from .logparse import column_wind

PARENT_START = ".parent_start"


class RunError(Exception):
    """A step cannot start: what the step before it should have written is missing."""


class Step:
    def __init__(self, label, cwd, args):
        self.label = label
        self.cwd = cwd
        self.args = args              # gust's arguments: the deck and key=value
        self.log = os.path.join(cwd, "log.txt")


def newest_checkpoint(folder, prefix, since=None):
    """The name of the checkpoint ``<prefix><step>`` in ``folder`` with the highest step (by number, as run.sh:
    chk100000 after chk99999), among those modified at or after ``since`` (seconds), or None."""
    best, best_step = None, -1
    pattern = re.compile(re.escape(prefix) + r"(\d+)$")
    for name in os.listdir(folder):
        m = pattern.match(name)
        path = os.path.join(folder, name)
        if not m or not os.path.isdir(path):
            continue
        if since is not None and os.path.getmtime(path) < since:
            continue
        step = int(m.group(1))
        if step > best_step:
            best, best_step = name, step
    return best


def steps(case_dir, runs, skip_finished=True, note=None):
    """Yield the Steps of the case written in ``case_dir``; ``runs`` is the builder's Result.runs. With
    ``skip_finished`` (run.resume) a direction that holds its answer for the current decks is skipped, as run.sh
    skips it (gust_case/resume.py), and ``note`` (a callable) gets a line saying so."""
    for run in runs:
        base = os.path.normpath(os.path.join(case_dir, run["dir"]))
        if skip_finished and resume.finished(base, run):
            if note is not None:
                note("%s: finished earlier with these decks and data files, skipped (run.resume = false runs it again)" % run["label"])
            continue
        try:
            resume.start(base, run)
        except FileNotFoundError as err:          # a deck of the direction is missing
            raise RunError("%s: %s; write the case again" % (run["label"], err))
        except OSError as err:                    # an old log held open by another program
            raise RunError("%s: %s" % (run["label"], err))
        extra = []
        if run["column"]:
            yield Step("%s: column" % run["label"], os.path.join(base, "column"), ["inputs"])
            guess = os.path.join(base, COLUMN_TABLE)
            if not (os.path.isfile(guess) and os.path.getsize(guess) > 0):
                raise RunError("the column run wrote no %s (see column/log.txt)" % COLUMN_TABLE)
            if run["wind_from_column"]:
                with open(os.path.join(base, "column", "log.txt"), encoding="utf-8", errors="replace") as handle:
                    wind = column_wind(handle.read())
                if wind is None:
                    raise RunError("no geostrophic wind in %s" % os.path.join(base, "column", "log.txt"))
                extra = ["forcing.ug=%s" % wind[0], "forcing.vg=%s" % wind[1]]
        since = None
        if run["patches"]:
            marker = os.path.join(base, PARENT_START)  # the patches take a checkpoint this run writes
            with open(marker, "w", encoding="utf-8"):
                pass
            since = os.path.getmtime(marker)
        yield Step(run["label"], base, ["inputs"] + extra)
        if run["patches"]:
            parent = newest_checkpoint(base, run["chk"], since)
            if parent is None:
                raise RunError("the run in %s wrote no checkpoint %s* (see its log.txt)" % (base, run["chk"]))
            for patch in run["patches"]:
                yield Step("%s: patch %s" % (run["label"], patch), os.path.join(base, "nest_" + patch),
                           ["inputs", "bc.nest_parent=../%s" % parent] + extra)


def collect_library(case_dir):
    """After the last step: library.json and the case summary's Library section for a geostrophic drive, as run.sh
    writes them (gust_case/library.py). Returns a line for the log ('' when the drive has no column run)."""
    if not library.is_geostrophic(case_dir):
        return ""
    try:
        path = library.write(case_dir, library.collect(case_dir))
    except (library.CaseError, OSError, ValueError, KeyError) as err:
        problems = err.problems if isinstance(err, library.CaseError) else [str(err)]
        return "library.json not written: %s" % "; ".join(problems)
    return "wrote %s" % path


def find_gust(executable, case_dir, frozen_dir=None, given=True):
    """The gust to run, or None. ``run.executable`` (a path, relative to the case folder, or a name on PATH) when
    the case gives it (``given``); a frozen window (PyInstaller) then falls back to the gust beside its own exe.
    When the case leaves run.executable at its default, the gust beside the frozen exe comes first."""
    exe = "gust.exe" if sys.platform.startswith("win") else "gust"
    bundled = os.path.join(frozen_dir, exe) if frozen_dir is not None else None

    def named():
        if not executable:
            return None
        if os.path.dirname(executable):
            path = executable if os.path.isabs(executable) else os.path.join(case_dir, executable)
            return path if os.path.isfile(path) else None
        return shutil.which(executable)

    if bundled is not None and not given and os.path.isfile(bundled):
        return bundled
    found = named()
    if found is None and bundled is not None and os.path.isfile(bundled):
        return bundled
    return found


def command(gust, args, ranks, launcher):
    """(program, arguments) of one step: gust alone with one rank, else under the MPI launcher of run.mpi_launcher
    ({ranks} replaced), as run.sh."""
    if ranks <= 1:
        return gust, list(args)
    windows = sys.platform.startswith("win")
    words = shlex.split(launcher.replace("{ranks}", str(ranks)), posix=not windows)
    if windows:                                   # posix=False keeps the quotes of "C:\Program Files\...\mpiexec.exe"
        words = [w[1:-1] if len(w) > 1 and w[0] == w[-1] == '"' else w for w in words]
    if not words:
        return gust, list(args)
    return words[0], words[1:] + [gust] + list(args)
