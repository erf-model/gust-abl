#!/usr/bin/env python3
"""
The one-line mutants of the force file in unsteady runs and restarts
(Source/App.cpp, Source/BuildingLoads.cpp), each run against the CTests
force_file* (Tests/CMakeLists.txt, the cube of Canonical_RANS/Building_Cube,
inputs_cube_unsteady). The table is in docs/testing.rst.

    python3 tools/force_file_mutants.py BUILD_DIR [FF1 FF3 ...]

Each mutant is applied to the source, gust is rebuilt and the CTests run; the
source is restored and gust rebuilt at the end.
"""
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = "Source/App.cpp"
L = "Source/BuildingLoads.cpp"
FILL = "bc.fill_ghosts(state, res_in.surface, fopt.buoyancy, false, nullptr);\n            ResidualInputs fin"
M = [
    ("FF1", "the force evaluation changes the state (u + 1e-12)", A, FILL,
     "state.u.plus(amrex::Real(1e-12), 0, 1, 0); " + FILL),
    ("FF2", "a restart drops the force row at the checkpoint's time", L,
     "std::strtod(line.c_str(), nullptr) <= keep_until", "std::strtod(line.c_str(), nullptr) < keep_until"),
    ("FF3", "a restart that finds no force file reported as kept (no first row)", L,
     "return keep_until >= 0 && exists;", "return keep_until >= 0;"),
    ("FF4", "a restart's force file of another header accepted", L,
     "if (line != kHeader) {", "if (false) {"),
]


def build(b):
    return subprocess.run(["cmake", "--build", b, "-j", "4", "--target", "gust"], capture_output=True, text=True).returncode == 0


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    b = os.path.abspath(sys.argv[1])
    only = sys.argv[2:]
    rows = []
    try:
        for name, what, f, old, new in M:
            if only and name not in only:
                continue
            p = os.path.join(ROOT, f)
            with open(p) as fh:
                src = fh.read()
            if src.count(old) != 1:
                sys.exit(f"{name}: the line to mutate is not found once in {f}")
            with open(p, "w") as fh:
                fh.write(src.replace(old, new))
            try:
                if not build(b):
                    rows.append((name, what, "build failed", ""))
                    continue
                r = subprocess.run(["ctest", "--test-dir", b, "-R", "force_file"], capture_output=True, text=True)
                failed = [ln.split(" - ", 1)[1].split()[0] for ln in r.stdout.splitlines() if "(Failed)" in ln]
                rows.append((name, what, "killed" if failed else "SURVIVED", ", ".join(failed)))
            finally:
                with open(p, "w") as fh:
                    fh.write(src)
    finally:
        build(b)
    for row in rows:
        print(" | ".join(row))
    return 0 if all(r[2] == "killed" for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
