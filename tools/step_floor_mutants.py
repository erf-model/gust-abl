#!/usr/bin/env python3
"""One-line mutants of the local step floor (Source/PseudoTime.cpp), the canopy in the column first guess
(Source/ColumnGuess.cpp) and the canopy's column profile (Source/Canopy.cpp), each built and run against the unit
tests PseudoTime.StepFloor*, ColumnCanopy.* and Canopy.ColumnProfile* on two and four ranks (docs/testing.rst).

    step_floor_mutants.py BUILD_DIR [--ranks 2,4]

Run from the repository root with a clean tree: each mutant is applied with Python's str.replace (exactly one
occurrence), the unit tests are rebuilt, and the file is restored with git checkout.
"""
import subprocess
import sys

MUTANTS = [
    ("F1", "Source/PseudoTime.cpp", "The floor ignored",
     "min_rate = amrex::max(min_rate, phi * amrex::max(umax[0], umax[1]) / dh);", "min_rate = min_rate;"),
    ("F2", "Source/PseudoTime.cpp", "U_max from u alone",
     "phi * amrex::max(umax[0], umax[1]) / dh", "phi * umax[0] / dh"),
    ("F3", "Source/PseudoTime.cpp", "dx for min(dx, dy)",
     "const amrex::Real dh = amrex::min(geom.CellSize(0), geom.CellSize(1));", "const amrex::Real dh = geom.CellSize(0);"),
    ("F4", "Source/PseudoTime.cpp", "No all-reduce of U_max",
     "        amrex::ParallelAllReduce::Max(umax, 2, amrex::ParallelDescriptor::Communicator());\n", ""),
    ("F5", "Source/PseudoTime.cpp", "A negative floor accepted", "!(phi >= 0)", "!(phi >= -1)"),
    ("C1", "Source/ColumnGuess.cpp", "The column's drag dropped", "dW[k] += 2.0 * fw;", "dW[k] += 0.0 * fw;"),
    ("C2", "Source/ColumnGuess.cpp", "Picard drag (|W| lagged)",
     "dW[k] += 2.0 * fw;\n            rW[k] += fw * W[k];", "dW[k] += fw;"),
    ("C3", "Source/ColumnGuess.cpp", "No wake production",
     "rK[k] += beta_p * f_at(k) * std::pow(std::abs(W[k]), 3);", "rK[k] += 0.0;"),
    ("C4", "Source/ColumnGuess.cpp", "No short-circuit sink",
     "pat[k] += beta_d * f_at(k) * std::abs(W[k]);", "pat[k] += 0.0;"),
    ("C5", "Source/ColumnGuess.cpp", "No length cap",
     "const amrex::Real alpha = canopy_k ? static_cast<amrex::Real>(canopy.tke.length_alpha) : amrex::Real(0);",
     "const amrex::Real alpha = amrex::Real(0);"),
    ("C6", "Source/ColumnGuess.cpp", "The profile's size unchecked",
     "if (has_canopy && static_cast<int>(canopy.f.size()) != nz) {", "if (false) {"),
    ("C7", "Source/ColumnGuess.cpp", "A non-uniform canopy accepted",
     "if (spread > kColumnCanopySpread) {", "if (spread > 1.0e9) {"),
    ("C8", "Source/ColumnGuess.cpp", "Terrain unchecked", "    if (!flat_terrain) {", "    if (false) {"),
    ("P1", "Source/Canopy.cpp", "The layer sum, not the mean",
     "for (double& v : prof) { v /= ncols; fmean_max = std::max(fmean_max, v); }",
     "for (double& v : prof) { fmean_max = std::max(fmean_max, v); }"),
    ("P2", "Source/Canopy.cpp", "The spread absolute, not relative",
     "spread = (fmean_max > 0) ? worst / fmean_max : 0.0;", "spread = worst;"),
    ("P3", "Source/Canopy.cpp", "The spread from the maximum alone",
     "worst = std::max(worst, static_cast<double>(fmax[static_cast<std::size_t>(k)]) - static_cast<double>(fmin[static_cast<std::size_t>(k)]));",
     "worst = std::max(worst, static_cast<double>(fmax[static_cast<std::size_t>(k)]) - static_cast<double>(fmax[static_cast<std::size_t>(k)]));"),
]
FILTER = "PseudoTime.StepFloor*:ColumnCanopy.*:Canopy.ColumnProfile*"


def main(argv):
    if len(argv) < 2:
        raise SystemExit(__doc__)
    build = argv[1]
    ranks = [2, 4]
    if "--ranks" in argv:
        ranks = [int(x) for x in argv[argv.index("--ranks") + 1].split(",")]
    survived = []
    for tag, path, what, old, new in MUTANTS:
        src = open(path).read()
        if src.count(old) != 1:
            raise SystemExit("%s: the text to mutate occurs %d times in %s" % (tag, src.count(old), path))
        open(path, "w").write(src.replace(old, new))
        try:
            build_ok = subprocess.run(["cmake", "--build", build, "-j", "8", "--target", "unit_tests"],
                                      capture_output=True, text=True).returncode == 0
            killed = False
            if build_ok:
                for n in ranks:
                    out = subprocess.run(["mpiexec", "-n", str(n), build + "/unit_tests", "--gtest_filter=" + FILTER],
                                         capture_output=True, text=True)
                    if out.returncode != 0 or "[  FAILED  ]" in out.stdout:
                        killed = True
            print("%s %-36s %s" % (tag, what, "killed" if killed else ("build failed" if not build_ok else "SURVIVED")))
            if not killed:
                survived.append(tag)
        finally:
            subprocess.run(["git", "checkout", "--", path], check=True)
    subprocess.run(["cmake", "--build", build, "-j", "8", "--target", "unit_tests"], capture_output=True)
    print("survived: %s" % (", ".join(survived) if survived else "none"))
    return 1 if survived else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
