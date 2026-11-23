"""Running a written case: each step of runplan.steps under QProcess, its log streamed to the window and to the
step's log.txt, and the residual columns of the iteration table plotted as they come."""
import codecs
import math

import pyqtgraph as pg
from PySide6 import QtCore, QtGui, QtWidgets

from . import mapview  # noqa: F401  (pyqtgraph's white background)
from .logparse import RESIDUALS, LogReader
from .runplan import RunError, collect_library, command, steps

STOP_GRACE_MS = 5000             # Stop: terminate (mpiexec passes it to the ranks), then kill after this

COLORS = {"u": (31, 119, 180), "v": (255, 127, 14), "w": (44, 160, 44), "tke": (214, 39, 40),
          "theta": (148, 103, 189), "p": (140, 86, 75), "r": (0, 0, 0)}


class RunPanel(QtWidgets.QWidget):
    """The log and the residual plot of the step that runs."""

    finished = QtCore.Signal(bool, str)          # all steps ended with status 0; the last message

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QtWidgets.QVBoxLayout(self)
        self.status = QtWidgets.QLabel("Nothing runs.")
        layout.addWidget(self.status)
        split = QtWidgets.QSplitter(QtCore.Qt.Vertical)
        self.plot = pg.PlotWidget()
        self.plot.setLogMode(x=False, y=True)
        self.plot.setLabel("bottom", "iteration")
        self.plot.setLabel("left", "residual")
        self.plot.addLegend(offset=(-5, 5))
        split.addWidget(self.plot)
        self.log = QtWidgets.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setMaximumBlockCount(20000)
        self.log.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        font = self.log.font()
        font.setFamily("Menlo" if QtCore.QSysInfo.productType() == "macos" else "Consolas")
        self.log.setFont(font)
        split.addWidget(self.log)
        layout.addWidget(split)
        self.stop_button = QtWidgets.QPushButton("Stop")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop)
        layout.addWidget(self.stop_button, alignment=QtCore.Qt.AlignRight)
        self.process = None
        self._steps = None
        self._step = None
        self._log_file = None
        self._reader = None
        self._curves = {}
        self._data = {}
        self._decoder = None
        self._stopping = False

    def running(self):
        return self.process is not None

    def start(self, case_dir, runs, gust, ranks, launcher, skip_finished=True):
        """Run every step of the case written in ``case_dir`` with the executable ``gust``; with ``skip_finished``
        (run.resume) the directions that hold their answer are skipped, as run.sh skips them."""
        self.log.clear()
        self._gust, self._ranks, self._launcher = gust, ranks, launcher
        self._steps = steps(case_dir, runs, skip_finished, lambda text: self.log.appendPlainText("==== %s\n" % text))
        self._case_dir = case_dir
        self._stopping = False
        self.stop_button.setEnabled(True)
        self._next(None)

    def _next(self, previous):
        try:
            step = next(self._steps) if previous is None else self._steps.send(None)
        except StopIteration:
            note = collect_library(self._case_dir)
            if note:
                self.log.appendPlainText("==== %s\n" % note)
            self._done(True, "All steps ended with status 0." + (" " + note if note else ""))
            return
        except (RunError, OSError, ValueError) as err:
            self._done(False, str(err))
            return
        try:
            self._begin(step)
        except (OSError, ValueError) as err:     # the log cannot be written, the launcher cannot be read
            self._close_log()
            self.process = None
            self._done(False, "%s did not start: %s" % (step.label, err))

    def _begin(self, step):
        self._step = step
        self._reader = LogReader()
        self.plot.clear()
        self._curves, self._data = {}, {}
        for name in RESIDUALS + ("r",):
            self._data[name] = ([], [])
            self._curves[name] = self.plot.plot([], [], pen=pg.mkPen(COLORS[name], width=1.5),
                                                name=name if name == "r" else name + "-res")
        program, args = command(self._gust, step.args, self._ranks, self._launcher)
        self._program = program
        self._decoder = codecs.getincrementaldecoder("utf-8")(errors="replace")   # a character split between reads
        self._log_file = open(step.log, "w", encoding="utf-8", newline="")         # gust's own line ends
        self.status.setText("%s: %s" % (step.label, step.cwd))
        self.log.appendPlainText("==== %s\n$ cd %s\n$ %s %s\n" % (step.label, step.cwd, program, " ".join(args)))
        self.process = QtCore.QProcess(self)
        self.process.setWorkingDirectory(step.cwd)
        self.process.setProcessChannelMode(QtCore.QProcess.MergedChannels)
        self.process.readyReadStandardOutput.connect(self._read)
        self.process.finished.connect(self._ended)
        self.process.errorOccurred.connect(self._failed)
        self.process.start(program, args)

    def _read(self):
        text = self._decoder.decode(bytes(self.process.readAllStandardOutput()))
        if not text:
            return
        self._log_file.write(text)
        self.log.moveCursor(QtGui.QTextCursor.End)
        self.log.insertPlainText(text)
        for row in self._reader.feed(text):
            for name in RESIDUALS:
                value = row.residuals[name]
                if math.isfinite(value) and value > 0.0:
                    self._data[name][0].append(row.x)
                    self._data[name][1].append(value)
            if math.isfinite(row.r) and row.r > 0.0:
                self._data["r"][0].append(row.x)
                self._data["r"][1].append(row.r)
        for name, (xs, ys) in self._data.items():
            self._curves[name].setData(xs, ys)

    def _close_log(self):
        if self._log_file is not None:
            self._log_file.close()
            self._log_file = None

    def _ended(self, code, status):
        self._read()
        self._close_log()
        self.process.deleteLater()
        self.process = None
        if self._stopping:
            self._done(False, "%s stopped (Stop)" % self._step.label)
            return
        if status != QtCore.QProcess.NormalExit or code != 0:
            self._done(False, "%s ended with status %d (see %s)" % (self._step.label, code, self._step.log))
            return
        self._next(self._step)

    def _failed(self, error):
        if error != QtCore.QProcess.FailedToStart:
            return                                # a crash also ends in _ended, with its status
        self._close_log()
        self.process.deleteLater()
        self.process = None
        self._done(False, "%s did not start: %s" % (self._step.label, self._program))

    def stop(self):
        """Ask the step to end (SIGTERM; mpiexec passes it to its ranks), and kill it if it has not after
        STOP_GRACE_MS. On Windows terminate() only closes windows, so the kill ends a console gust."""
        if self.process is None or self._stopping:
            return
        self._stopping = True
        self.status.setText("Stopping %s..." % self._step.label)
        process = self.process
        process.terminate()
        QtCore.QTimer.singleShot(STOP_GRACE_MS, self, lambda: self._kill(process))

    def _kill(self, process):
        if self.process is process:
            process.kill()

    def _done(self, ok, message):
        self._steps = None
        self.stop_button.setEnabled(False)
        self.status.setText(message)
        self.log.appendPlainText("\n==== " + message)
        self.finished.emit(ok, message)
