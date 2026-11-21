# PyInstaller spec of the case set-up window (tools/gust_gui.py): a one-folder, windowed build.
#
#     cd tools
#     pyinstaller --noconfirm gust_gui.spec          # dist/gust_gui/gust_gui(.exe)
#
# Put the CPU gust.exe (or gust) in dist/gust_gui/ beside gust_gui.exe: the frozen window looks for gust there
# first when run.executable is not given. rasterio is left out: GeoTIFF terrain goes through the builder's own
# reader. Tests/inputs_reference (the deck checker's key table) is carried along. docs/case_setup.rst describes
# the build.
# -*- mode: python ; coding: utf-8 -*-
import os

HERE = os.path.abspath(SPECPATH)  # noqa: F821  (set by PyInstaller)

a = Analysis(  # noqa: F821
    [os.path.join(HERE, "gust_gui.py")],
    pathex=[HERE],
    binaries=[],
    datas=[(os.path.join(HERE, "..", "Tests", "inputs_reference"), "Tests")],
    hiddenimports=["make_gust_case", "gust_case", "gust_input_checker", "gust_window", "gust_window.window",
                   "gust_window.form", "gust_window.layout", "gust_window.mapview", "gust_window.runner",
                   "gust_window.logparse", "gust_window.runplan", "gust_window.values", "gust_window.librarytab",
                   "gust_library", "gust_library.ensemble", "gust_library.grids", "yaml", "numpy"],
    hookspath=[],
    runtime_hooks=[],
    excludes=["rasterio", "tkinter", "matplotlib", "PyQt5", "PyQt6", "PySide2"],
    noarchive=False,
)
pyz = PYZ(a.pure)  # noqa: F821
exe = EXE(  # noqa: F821
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="gust_gui",
    console=False,                    # --windowed
    upx=False,
)
coll = COLLECT(  # noqa: F821
    exe,
    a.binaries,
    a.datas,
    name="gust_gui",                  # --onedir
    upx=False,
)
