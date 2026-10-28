#!/usr/bin/env python3
"""
The one-line mutants of simplec.theta = frozen (Source/Simplec.cpp,
Source/Residual.cpp, Source/Predictor.cpp, Source/Sequencing.cpp,
Source/App.cpp).

F1-F7 and F16-F19 run against FrozenTheta.* by the shared runner
tools/mutants.py (two ranks, the unit-test flags). F8-F15 and F20-F23 change the
wiring of the application, which
the unit tests do not reach: each rebuilds gust and runs the CTests
rt_frozen_theta and abort_simplec_theta_* (about a minute each). The table is
in docs/testing.rst.

    python3 tools/frozen_theta_mutants.py BUILD_DIR [F1 F9 ...]
"""
import os
import subprocess
import sys

import mutants

S, R, P, Q, A = "Source/Simplec.cpp", "Source/Residual.cpp", "Source/Predictor.cpp", "Source/Sequencing.cpp", "Source/App.cpp"
B = "Source/BoundaryConditions.cpp"
UNIT = [
 ("F1", "the buoyancy left in the momentum residual", R,
  "if (in.buoyancy_in_momentum) { subtract_buoyancy", "if (true) { subtract_buoyancy"),
 ("F2", "the predictor's point-implicit buoyancy left in the w rows", P,
  "in.buoyancy.gravity > 0 && in.buoyancy_in_momentum)", "in.buoyancy.gravity > 0)"),
 ("F3", "theta's residual not zeroed", S,
  "if (m_opt.theta_frozen) { m_R.theta.setVal(0); }", ""),
 ("F4", "theta swept although frozen", S,
  "    if (m_opt.theta_frozen) {\n        // theta is held", "    if (false) {\n        // theta is held"),
 ("F5", "simplec.theta = frozen read as solve", S,
  'if (theta == "frozen") { opt.theta_frozen = true; }', 'if (theta == "frozen") { }'),
 ("F6", "any value of simplec.theta accepted", S,
  'else if (theta != "solve") {', 'else if (false) {'),
 ("F7", "the stage accepts residual inputs that do not match the option", S,
  "if (opt.theta_frozen == in.buoyancy_in_momentum) {", "if (false) {"),
 ("F16", "the ground extension left out of the fill", B,
  "if (m_theta_ground_extension && !homogeneous) {", "if (false) {"),
 ("F17", "the ground extension in homogeneous fills too", B,
  "if (m_theta_ground_extension && !homogeneous) {", "if (m_theta_ground_extension) {"),
 ("F18", "only the first ghost layer extended", B,
  "for (int m = 1; m <= 2; ++m) {", "for (int m = 1; m <= 1; ++m) {"),
 ("F19", "the extension's slope halved", B,
  "(a(i, j, klo + 1) - a(i, j, klo)) / (z1 - z0)", "(a(i, j, klo + 1) - a(i, j, klo)) / (z1 - z0) * amrex::Real(0.5)"),
]
APP = [
 ("F8", "the residual inputs keep the momentum buoyancy with frozen theta", Q,
  "    r.buoyancy_in_momentum = !s.simplec.theta_frozen;\n", ""),
 ("F9", "the coarsest level not reset", Q,
  "            print_initializer(r, iopt, setup.popt, label);\n            reset_frozen_theta(h, l, setup, iopt, ls.bc.profile(), *ql);\n",
  "            print_initializer(r, iopt, setup.popt, label);\n"),
 ("F10", "a prolonged coarse level not reset", Q,
  "            s.prolonged = true;\n            reset_frozen_theta(h, l, setup, iopt, ls.bc.profile(), *ql);\n",
  "            s.prolonged = true;\n"),
 ("F11", "level 0 not reset", Q,
  "    reset_frozen_theta(h, 0, setup, iopt, finest.bc.profile(), q);\n", ""),
 ("F12", "theta from the sounding by absolute height, not by height above the ground", Q,
  "&inflow, true);", "&inflow, false);"),
 ("F13", "frozen accepted without grid sequencing", A,
  "if (simplec.theta_frozen && !(sequencing && in.mesh.hierarchy_levels > 1)) {", "if (false) {"),
 ("F14", "frozen accepted with the Newton stage", A,
  "if (newton && simplec.theta_frozen) {", "if (false) {"),
 ("F15", "frozen accepted with a restart", A,
  "if (simplec.theta_frozen && !in.io.restart.empty()) {", "if (false) {"),
 ("F20", "frozen accepted with a stratified column first guess", A,
  "if (simplec.theta_frozen && column_sets_theta(iopt, surf)) {", "if (false) {"),
 ("F21", "the turbulence plotfile's R_theta not zeroed", A,
  "if (simplec.theta_frozen) { R.theta.setVal(0); }", ""),
 ("F22", "the frozen theta read at the floored height (init.z_agl_min)", Q,
  "    topt.z_agl_min = 0;\n", ""),
 ("F23", "the levels do not extend theta below the ground", Q,
  "    bc.set_theta_ground_extension(s.simplec.theta_frozen);\n", ""),
]
CTESTS = "rt_frozen_theta|abort_simplec_theta"


def ctest(build):
    b = subprocess.run(["cmake", "--build", build, "--target", "gust"], capture_output=True, text=True)
    if b.returncode:
        return "not compiled", []
    r = subprocess.run(["ctest", "--test-dir", build, "-R", CTESTS], capture_output=True, text=True)
    failed = [ln.split()[3] for ln in r.stdout.splitlines() if "***Failed" in ln or "Timeout" in ln]
    return ("killed" if failed else "SURVIVED"), failed


def run_app(build, only):
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    base, failed = ctest(build)
    if base != "SURVIVED":
        sys.exit(f"baseline: {base} {', '.join(failed)}")
    survived = 0
    try:
        for mid, what, f, old, new in APP:
            if only and mid not in only:
                continue
            path = os.path.join(root, f)
            src = open(path).read()
            if src.count(old) != 1:
                sys.exit(f"{mid}: the original text occurs {src.count(old)} times in {f}")
            try:
                open(path, "w").write(src.replace(old, new))
                res, failed = ctest(build)
            finally:
                open(path, "w").write(src)
            survived += res != "killed"
            print(f"{mid} | {what} | {res} | {', '.join(failed)}", flush=True)
    finally:
        subprocess.run(["cmake", "--build", build, "--target", "gust"], capture_output=True)
    return survived


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    ids = sys.argv[2:]
    unknown = sorted(set(ids) - {m[0] for m in UNIT + APP})
    if unknown:
        sys.exit(f"unknown mutant ids: {', '.join(unknown)}")
    app_ids = [i for i in ids if i in {m[0] for m in APP}]
    unit_ids = [i for i in ids if i not in app_ids]
    bad = 0
    if not ids or app_ids:
        bad = run_app(os.path.abspath(sys.argv[1]), app_ids)
    if not ids or unit_ids:
        sys.argv = sys.argv[:2] + unit_ids
        try:
            mutants.run(UNIT, "FrozenTheta.*", set(), __doc__)
        except SystemExit as e:
            bad += int(e.code or 0)
    sys.exit(1 if bad else 0)
