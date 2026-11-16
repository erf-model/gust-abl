"""The checks, one module per area. Importing a module registers its checkers (checker.register)."""
from . import keys, values, mesh, boundary, physics, solver, objects, files, nests  # noqa: F401
