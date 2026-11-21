"""A wind library from a GUST sweep: winds at heights above the ground for every direction.

    index.py      build(case_dir, heights, reference_point, directions): <case_dir>/wind_library/
    query.py      Library(path).wind(direction, speed, height): the two neighbouring sectors blended
    sample.py     fields at heights above the ground and at a point
    plotfile.py   GUST plotfiles (cell fields and node heights) as numpy arrays
    ensemble.py   statistics over a distribution of directions and speeds (no new runs)
    grids.py      WindNinja-format output: ESRI ASCII grids with a .prj, and CSV at points

tools/make_gust_library.py is the command line; docs/case_setup.rst describes
the files and the conventions.
"""
from .index import LibraryError, build  # noqa: F401
from .query import Library  # noqa: F401
from .sample import SampleError  # noqa: F401
