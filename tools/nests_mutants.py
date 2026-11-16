#!/usr/bin/env python3
"""
The one-line mutants of the nest inputs (Source/Nests.cpp): the names and
corners, the placement checks and the settings a nest derives from the
outer deck.

Run against the Nests.* unit tests by the shared runner tools/mutants.py (two
ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/nests_mutants.py BUILD_DIR [NS1 NS4 ...]
"""
import mutants

F = "Source/Nests.cpp"
M = [
 # reading
 ("NS1", "'names' accepted as a nest name", F, '|| n == "names"', ""),
 ("NS2", "a corner given twice not refused", F, "pn.countname(key) > 1", "pn.countname(key) > 2"),
 # placement
 ("NS3", "a margin of 1 outer cell", F, "kNestMargin;", "kNestMargin - 1;"),
 ("NS4", "nests that touch taken as overlapping", F, "c[0] < b[1]", "c[0] <= b[1]"),
 ("NS5", "a corner up to a whole outer cell off the faces accepted", F,
  "std::abs(c - n) <= 1.0e-4", "std::abs(c - n) <= 1.0"),
 ("NS6", "a nest of zero width accepted", F, "if (e <= o)", "if (e < o)"),
 # the derived settings
 ("NS7", "the nest at the outer cell size (no refinement)", F, "kNestRatio * (e - o)", "(e - o)"),
 ("NS8", "the nest periodic where the outer domain is, whatever it spans", F,
  "(outer.is_periodic[d] != 0 && o == 0 && e == outer.n_cell[d]) ? 1 : 0", "outer.is_periodic[d]"),
 ("NS9", "the outer grid hierarchy kept", F, "m.hierarchy_levels = 1;", ""),
 ("NS10", "the outer io.restart kept", F, "in.io.restart.clear();", ""),
 ("NS11", "the Rayleigh layer kept on k", F, "bin.rayleigh_fields[4] = false;", ""),
 ("NS12", "the outer relaxation zone kept", F, "bin.nest_relax_width = 0;", ""),
 # beyond the first review
 ("NS13", "a corner near an outer face kept as given (not moved onto it)", F,
  "m.prob_lo[d] = static_cast<amrex::Real>(static_cast<double>(outer.prob_lo[d]) + o * dx);", "m.prob_lo[d] = nest.lo[d];"),
 ("NS14", "a nest over the whole outer domain accepted", F,
  "if (c[0] == 0 && c[1] == outer.n_cell[0]", "if (false && c[0] == 0 && c[1] == outer.n_cell[0]"),
 ("NS15", "nest faces where the nest is periodic", F,
  "(periodic[d] != 0) ? FaceType::periodic : FaceType::nest", "FaceType::nest"),
 ("NS16", "the Agnesi hill centred on the domain whatever hill_centre says", "Source/TerrainSurface.cpp",
  "std::isnan(t.hill_centre[0])", "true"),
 ("NS17", "a nest on the outer high face not given that face exactly", F,
  "if (e == outer.n_cell[d]) { m.prob_hi[d] = outer.prob_hi[d]; }", ""),
 ("NS18", "the outer forcing taken when the outer run kept none", F,
  "if (!out.forcing_set) { return; }", ""),
 ("NS19", "the outer surface heat flux not taken", F,
  "surf.heat_flux = out.surface_heat_flux;", ""),
 ("NS20", "the outer y pressure gradient not kept", F,
  "out.pressure_gradient = {fopt.pressure_gradient[0], fopt.pressure_gradient[1]};",
  "out.pressure_gradient = {fopt.pressure_gradient[0], 0};"),
]
EQUIVALENT = set()


if __name__ == "__main__":
    mutants.run(M, "Nests.*", EQUIVALENT, __doc__)
