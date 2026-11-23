#!/usr/bin/env python3
"""
The one-line mutants of the two-way nests' injection (Source/NestInject.cpp,
nest.feedback = inject) and of its hold in the SIMPLEC iteration and the
Poisson solve (Source/Simplec.cpp, Source/Poisson.cpp): the factor, the blend,
the residual rows, the hold and the Poisson masks.

Run against the NestInject.* and MaskedPoisson.* unit tests by the shared runner tools/mutants.py
(two ranks, the unit-test flags). The table is in docs/testing.rst.

    python3 tools/nest_inject_mutants.py BUILD_DIR [J1 J4 ...]
"""
import mutants

F = "Source/NestInject.cpp"
S = "Source/Simplec.cpp"
P = "Source/Poisson.cpp"
M = [
 # the factor and the blend
 ("J1", "a factor above 1 accepted", F, "factor > 0 && factor <= 1)", "factor > 0 && factor <= 2)"),
 ("J2", "no blend: every exchange injects fully", F, "        if (!full) {\n", "        if (false) {\n"),
 ("J3", "the blend from q* toward q", F, "ta(i, j, k) = qa(i, j, k) + omega * (ta(i, j, k) - qa(i, j, k));",
  "ta(i, j, k) = qa(i, j, k) + omega * omega * (ta(i, j, k) - qa(i, j, k));"),
 # the residual rows
 ("J4", "k's rows kept on the inner box", F, "{&R.u, &R.v, &R.w, &R.p, &R.theta, &R.tke}) { zero_points(*r); }",
  "{&R.u, &R.v, &R.w, &R.p, &R.theta}) { zero_points(*r); }"),
 # the hold
 ("J5", "no hold after the momentum predictor", S, "    if (held) { m_in.inject->hold(q, 3); }\n    // 4.",
  "    // 4."),
 ("J6", "no hold at the end of the iteration", S, "    if (held) { m_in.inject->hold(q, 6); }\n", ""),
 # the Poisson masks
 ("J7", "the held faces open in the Poisson operator", P,
  "        if (held) { amrex::MultiFab::Multiply(*ab[d], *coeffs.held_open[d], 0, 0, 1, 0); }\n", ""),
 ("J8", "the held cells solved as fluid cells", P,
  "        if (held) { amrex::MultiFab::Multiply(m_fluid, *coeffs.held_fluid, 0, 0, 1, 0); }\n", ""),
 ("J9", "the inner box left out of the fluid mask", F, "        zero_region(g.inner, m_fluid);\n", ""),
 ("J10", "the inner box's faces left out of the open masks", F, "        for (auto& o : m_open) { zero_region(g.inner, o); }\n", ""),
 ("J11", "a flux through the held faces", P,
  "                amrex::MultiFab::Multiply((*fluxes)[d], *coeffs.held_open[d], 0, 0, 1, 0);", ""),
 ("J12", "phi shifted by its mean with a mask", P, "    if (res.singular && (!masked || m_level)) {", "    if (res.singular) {"),
]
EQUIVALENT = {}


if __name__ == "__main__":
    mutants.run(M, "NestInject.*:MaskedPoisson.*", EQUIVALENT, __doc__)
