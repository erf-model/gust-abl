#!/usr/bin/env python3
"""The case set-up window of GUST: a form over tools/make_gust_case.py and tools/check_gust_inputs.py.

    python tools/gust_gui.py                  a starter case
    python tools/gust_gui.py case.yaml        open a case file

It needs PySide6, pyqtgraph, numpy and PyYAML (pip install PySide6 pyqtgraph numpy pyyaml). The form has one
field per key of the case description; Check draws the domain, Write case writes the case folder (as
make_gust_case.py) and checks its decks, Run runs gust step by step and plots the residuals. docs/case_setup.rst
describes the window. A frozen build (gust_gui.spec, PyInstaller) looks for gust beside its own exe first.
"""
import os
import sys


def _paths():
    """The folder of the gust_case, gust_input_checker and gust_window packages, and the frozen exe's folder."""
    if getattr(sys, "frozen", False):
        return getattr(sys, "_MEIPASS", os.path.dirname(sys.executable)), os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__)), None


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    tools, frozen_dir = _paths()
    sys.path.insert(0, tools)
    try:
        from PySide6 import QtWidgets
        import pyqtgraph  # noqa: F401
    except ImportError as err:
        sys.stderr.write("gust_gui.py: %s; install PySide6 and pyqtgraph (pip install PySide6 pyqtgraph numpy "
                         "pyyaml)\n" % err)
        return 2
    from gust_window.window import MainWindow
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv[:1])
    app.setApplicationName("GUST case set-up")
    window = MainWindow(argv[0] if argv else None, frozen_dir=frozen_dir)
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
