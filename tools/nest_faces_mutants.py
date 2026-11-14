#!/usr/bin/env python3
"""
The one-line mutants of the nest faces of a nested patch (Source/NestFaces.cpp
and its hooks in Source/BoundaryConditions.cpp, Source/BoundaryDamping.cpp and
Source/BoundaryInputs.cpp).

Run against NestFaces.* by the shared runner tools/mutants.py (two ranks, the
unit-test flags). The table is in docs/testing.rst.

    python3 tools/nest_faces_mutants.py BUILD_DIR [F1 F4 ...]
"""
import mutants

F = "Source/NestFaces.cpp"
B = "Source/BoundaryConditions.cpp"
D = "Source/BoundaryDamping.cpp"
I = "Source/BoundaryInputs.cpp"
M = [
 # the net flux
 ("F1", "the net flux summed without the outward sign", F,
  "const amrex::Real V = sgn * v(i + ox, j + oy, k);", "const amrex::Real V = v(i + ox, j + oy, k);"),
 ("F2", "the weight without 1 / dx (wrong between x and y faces)", F,
  "const amrex::Real W = a(i + ox, j + oy, k) * dxinv;", "const amrex::Real W = a(i + ox, j + oy, k);"),
 ("F3", "the high face's flux read on the low side of its cell", F,
  "const int off = lo ? 0 : 1;", "const int off = 0;"),
 ("F4", "the correction added with the wrong sign", F,
  "const amrex::Real change = lo ? delta : -delta;", "const amrex::Real change = lo ? -delta : delta;"),
 ("F5", "the correction not applied", F,
  "if (is_nest(d, lo)) { shift_face(*vn[d], d, lo, geom, m_flux.delta); }", "if (false) { shift_face(*vn[d], d, lo, geom, m_flux.delta); }"),
 ("F6", "the correction applied with an outflow face", F,
  "m_flux.applied = !outflow;", "m_flux.applied = true;"),
 ("F7", "the correction on the valid points only (the ghost points keep the old value)", F,
  "        amrex::Box fb = mfi.fabbox();\n        if (fb.smallEnd(d) > face || fb.bigEnd(d) < face) { continue; }\n        if (lo) {",
  "        amrex::Box fb = mfi.validbox();\n        if (fb.smallEnd(d) > face || fb.bigEnd(d) < face) { continue; }\n        if (lo) {"),
 # the relaxation zone
 ("F8", "the source with the wrong sign", F,
  "r(i, j, k) += c(i, j, k) * (a(i, j, k) - p(i, j, k));", "r(i, j, k) -= c(i, j, k) * (a(i, j, k) - p(i, j, k));"),
 ("F9", "the weight reversed (0 at the face, 1 at the inner end)", F,
  "r(i, j, k) += (amrex::Real(1) - s * invN) * inv_tau", "r(i, j, k) += (s * invN) * inv_tau"),
 ("F10", "abs in place of max(0, U_in): the zone acts where the flow leaves", F,
  "inv_tau * amrex::max(amrex::Real(0), inward * U);", "inv_tau * std::abs(inward * U);"),
 ("F11", "the inward direction of a high face taken as the low face's", F,
  "const amrex::Real inward = lo ? amrex::Real(1) : amrex::Real(-1);", "const amrex::Real inward = amrex::Real(1);"),
 ("F12", "cell centres at the distance of the face points (no half cell)", F,
  "const amrex::Real half = stencil::nodal(T, d) ? amrex::Real(0) : amrex::Real(0.5);", "const amrex::Real half = amrex::Real(0);"),
 ("F13", "the rate without 1 / N", F,
  "const amrex::Real inv_tau = factor / (static_cast<amrex::Real>(N) * geom.CellSize(d));", "const amrex::Real inv_tau = factor / geom.CellSize(d);"),
 ("F14", "bc.nest_relax_factor ignored", F,
  "*vn[d], d, lo, geom, m_relax_width, m_relax_factor);", "*vn[d], d, lo, geom, m_relax_width, amrex::Real(1));"),
 ("F15", "the zone one cell short on a high face", F,
  "        z0 = lo ? face : face - N;", "        z0 = lo ? face : face - N + 1;"),
 ("F16", "the staggered point's U_in from one face point (no mean along the face)", F,
  "const int ie = dj * oe, je = di * oe;", "const int ie = 0, je = 0;"),
 ("F17", "w's U_in without the row below", F,
  "+ (v(fi, fj, k - oz) + v(fi - ie, fj - je, k - oz)));", "+ (v(fi, fj, k) + v(fi - ie, fj - je, k)));"),
 ("F18", "k relaxed with theta's parent field", F,
  "const amrex::MultiFab* pf[5] = {&m_par.u, &m_par.v, &m_par.w, &m_par.theta, &m_par.tke};",
  "const amrex::MultiFab* pf[5] = {&m_par.u, &m_par.v, &m_par.w, &m_par.theta, &m_par.theta};"),
 ("F19", "the diagonal twice the rate", F,
  "a(i, j, k, coef_P) += c(i, j, k);", "a(i, j, k, coef_P) += amrex::Real(2) * c(i, j, k);"),
 ("F20", "the diagonal of w from the cell-centred rate", F,
  "return (loc == Loc::xface) ? 0 : (loc == Loc::yface) ? 1 : (loc == Loc::zface) ? 2 : 3;",
  "return (loc == Loc::xface) ? 0 : (loc == Loc::yface) ? 1 : 3;"),
 ("F21", "the source not added (hook removed)", D,
  "if (m_nest != nullptr) { m_nest->add_relaxation(q, R); }", "if (false) { m_nest->add_relaxation(q, R); }"),
 # the fills and the rows
 ("F22", "one ghost layer filled", F,
  "gb.setSmall(d, lo ? edge - ng : edge + inner);", "gb.setSmall(d, lo ? edge - 1 : edge + inner);"),
 ("F23", "the normal velocity's face point not filled", F,
  "const int inner = nodal_d ? 0 : 1;", "const int inner = 1;"),
 ("F24", "the homogeneous fill copies the parent", F,
  "                if (homogeneous) {\n                    mf[mfi].setVal", "                if (false) {\n                    mf[mfi].setVal"),
 ("F25", "the row without the parent's value", F,
  "{ r(i, j, k) = a(i, j, k) - v(i, j, k); });", "{ r(i, j, k) = a(i, j, k); });"),
 ("F26", "the rows of the high faces missing", F,
  "if (is_nest(d, lo)) { face_rows(", "if (is_nest(d, lo) && lo) { face_rows("),
 ("F27", "the nest fill not called (the profile fill's values stay)", B,
  "if (m_nest != nullptr) { m_nest->fill(q, d, homogeneous); }", "if (false) { m_nest->fill(q, d, homogeneous); }"),
 ("F28", "a nest face without the parent's values not refused", B,
  "if (m_nest == nullptr && has_nest_face(m_in)) {", "if (false) {"),
 # the ground blend
 ("F29", "a linear blend (no smoothstep)", F,
  "const amrex::Real s = t * t * (amrex::Real(3) - amrex::Real(2) * t);", "const amrex::Real s = t;"),
 ("F30", "the blend measured from every lateral face, nest or not", F,
  "dist[n]) * invW + other[n]); }", "dist[n]) * invW); }"),
 ("F31", "the blend's two grounds swapped", F,
  "z(i, j, k) = (amrex::Real(1) - s) * h(ic, jc, k) + s * old;", "z(i, j, k) = s * h(ic, jc, k) + (amrex::Real(1) - s) * old;"),
 # the inputs and the restart record
 ("F32", "bc.nest_relax_width not bounded by the patch", I,
  "bc.nest_relax_width = \" + std::to_string(N) + \": must be below half the patch's cells in \"",
  "bc.nest_relax_width = \" + std::to_string(N + 1) + \": must be below half the patch's cells in \""),
 ("F33", "terrain.nest_blend_width of half the patch accepted", F,
  "if (!(W < amrex::Real(0.5) * L)) {", "if (!(W <= amrex::Real(0.5) * L)) {"),
 ("F34", "a nest face without bc.nest_parent accepted", F,
  "if (any && bin.nest_parent.empty()) {", "if (false) {"),
 ("F35", "a nest face without the parent's field as the first guess accepted", F,
  "if (any && !first_guess_nest && in.io.restart.empty()) {", "if (false) {"),
 ("F36", "the restart record written as another version", F,
  "kNestVersion = 1;", "kNestVersion = 2;"),
 ("F37", "a restart with other nest settings not refused", F,
  "if (active && rec != hash) {", "if (false) {"),
 ("F38", "a record of another version read as this one", F,
  "magic == kNestMagic && version == kNestVersion && key", "magic == kNestMagic && key"),
 ("F39", "the zone's diagonal not added (hook removed)", D,
  "if (m_nest != nullptr) { m_nest->add_relaxation_to_implicit(eq, loc, A); }", "if (false) { m_nest->add_relaxation_to_implicit(eq, loc, A); }"),
 # after the three reviews
 ("F41", "the correction on the face plane only (the normal velocity's ghost layers keep the parent's value)", F,
  "if (lo) { fb.setBig(d, face); } else { fb.setSmall(d, face); }", "fb.setBig(d, face); fb.setSmall(d, face);"),
 ("F42", "the restart record without the parent's content", F,
  's += "\\nparent content " + parent.content_digest();', 's += "\\nparent content ";'),
 ("F43", "the restart record with the nest faces only (not the other faces' types)", F,
  'kLateral[n] + "=" + face_type_name(t);', 'kLateral[n] + "=" + (t == FaceType::nest ? "nest" : "other");'),
 ("F44", "the restart record without the patch's place in the parent", F,
  's += "\\noffset " + std::to_string(parent.offset()[0]) + " " + std::to_string(parent.offset()[1]);', 's += "\\noffset ";'),
 # F45 (the restart record without the parent's step) is equivalent: the content digest covers the Header,
 # which holds the step. It is not run.
 ("F46", "a lateral no-slip wall beside a nest face accepted", I,
  "if (t == FaceType::slip_wall || t == FaceType::no_slip_wall) {\n                    err.emplace_back(std::string(\"bc.\") + kFaceKeys[o] + \" = \" + face_type_name(t) + \" is not supported with a nest \"",
  "if (t == FaceType::slip_wall) {\n                    err.emplace_back(std::string(\"bc.\") + kFaceKeys[o] + \" = \" + face_type_name(t) + \" is not supported with a nest \""),
 ("F47", "a Rayleigh layer on k towards the inflow profile accepted on a patch", I,
  "if (in.rayleigh_zdamp > 0 && (rf[0] || rf[1] || rf[4])) {", "if (in.rayleigh_zdamp > 0 && (rf[0] || rf[1])) {"),
 ("F48", "a sponge towards the inflow profile accepted on a patch", I,
  "if (in.sponge_width > 0) {\n                err.emplace_back(\"bc.sponge_width > 0 is not supported", "if (false) {\n                err.emplace_back(\"bc.sponge_width > 0 is not supported"),
 ("F49", "no warning for a large net flux", F,
  "const bool large = std::abs(m_flux.net) > kNestNetWarning * m_flux.inflow;", "const bool large = false;"),
 ("F50", "the warning whatever the net flux", F,
  "const bool large = std::abs(m_flux.net) > kNestNetWarning * m_flux.inflow;", "const bool large = true;"),
 ("F51", "the warning when nothing is removed (an outflow face)", F,
  "if (m_flux.applied && m_flux.inflow > 0 && large) {", "if (m_flux.inflow > 0 && large) {"),
 ("F52", "the content digest from the Header alone (not the fields' files)", "Source/NestParent.cpp",
  'for (const char* f : {"u", "v", "w", "theta", "tke", "z_phys_nd"}) {', "for (const char* f : std::vector<const char*>{}) {"),
]

# Equivalent mutants found while writing the tests, not in the table:
#  the zone's range starting on the face point of the normal velocity
#  (z0 = face): its rate there is replaced by the Dirichlet row
#  (apply_rows after add_damping in Residual.cpp, the identity row in
#  apply_to_implicit).
EQUIVALENT = set()


if __name__ == "__main__":
    mutants.run(M, "NestFaces.*", EQUIVALENT, __doc__)
