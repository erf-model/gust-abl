# Askervein_Nest

## Purpose
The one-way nesting validation of `Tests/Canonical_RANS/Nest_Validation` on measured terrain (sub-phase 9.3, beyond
plan rev. 6). The whole 6 km Askervein domain of `Askervein_210` runs on 40 m cells (the parent) and on 20 m cells
(the reference). A 20 m patch of half-width D centred on the hilltop HT refines the parent at 2 2 1 with four `nest`
faces. The measure is the wind speed 10 m above the ground at the 30 masts of lines A and AA: the patch against the
reference, as a function of D, beside the parent against the reference.

## Files
| File | Role |
|---|---|
| `inputs_askervein_domain` | The whole domain: `Askervein_210/inputs_steady_askervein_210` without its cells and output names, plus line LU-SGS and two-level grid sequencing (`steady.sweep`, `steady.sequencing`, `mesh.hierarchy_levels`). |
| `inputs_askervein_patch` | A patch: four `nest` faces, first guess from the parent, `terrain.nest_blend_width = 200` (the 20 m ground meets the parent's 40 m ground on the faces), line LU-SGS without sequencing. |
| `askervein_nest.py` | Runs a whole domain (`run`) or a patch (`nest`), samples the masts, compares (`compare`), does the study (`study`), or runs the two-way comparison (`two_way`). Pure Python; uses `tools/askervein_lines.py` and the run helpers of `Nest_Validation/nest_validation.py`. |

The terrain (`extern_erf/askervein.txt`), the inflow table and the mast positions
(`Askervein_210/askervein_210_inflow.txt`, `askervein_lines.csv`) are those of `Askervein_210`.

## Data and acknowledgements
No file is added by this case. It uses, as `Askervein_210` does:
- the Askervein terrain and inflow of ERF's `Exec/CanonicalTests/Real_Terrain/Askervein` (`extern_erf/`, recorded in
  `extern_erf/EXTRACTED.md`; BSD 3-clause with the LBNL clause, `extern_erf/LICENSE.md`). Thanks to the ERF developers;
- the mast positions of lines A and AA (P. A. Taylor and H. W. Teunissen, "The Askervein Hill project: overview and
  background data", Boundary-Layer Meteorology 39, 1987, doi:10.1007/BF00121863), as tabulated in N. Wagenbrenner's
  windninja-tools (github.com/nwagenbrenner/windninja-tools, `data/askervein/askervein_obs.zip`), used with her
  permission (recorded in `Askervein_210/README.md`). Thanks to Natalie Wagenbrenner and to the Askervein field teams.

All three files are already redistributed in this repository; nothing has to be downloaded.

## How to run
```bash
ctest --test-dir <build> -R "askervein_nest"
```
The study, in a scratch directory (1.6 million cells in the reference):
```bash
python3 askervein_nest.py study -- mpiexec -n 4 <build>/gust
```
Single runs:
```bash
python3 askervein_nest.py run par40 150 -- mpiexec -n 4 <build>/gust inputs_askervein_domain
python3 askervein_nest.py run ref20 300 -- mpiexec -n 4 <build>/gust inputs_askervein_domain
python3 askervein_nest.py nest nest_d1400 par40 1400 -- mpiexec -n 4 <build>/gust inputs_askervein_patch
python3 askervein_nest.py compare nest_d1400 ref20 par40 --ratio 0.5
```
Two-way nesting (PR-70, `nest.feedback = defect`): the 20 m reference, and the 40 m outer domain with a 20 m nest of
half-width 1000 m, one way and joint, one after another, with the recommended setting; prints the table of
`docs/boundary_conditions.rst` (Two-way nests, Results):
```bash
python3 askervein_nest.py two_way -- mpiexec -n 4 <build>/gust inputs_askervein_domain
```
`two_way relax inject` adds the joint runs with `nest.feedback = relax` (factor 1) and `nest.feedback = inject`
(factor 0.5), `relax_consistent` the run with `nest.feedback = relax_consistent` (factor 1; tracker N-114 item 12,
`docs/boundary_conditions.rst`, Pressure-consistent relaxation), and `keep` skips runs whose plotfiles are already there:
```bash
python3 askervein_nest.py two_way relax inject -- mpiexec -n 4 <build>/gust inputs_askervein_domain
```

## Targets
Fixed in this file before any run was compared with the reference. When they were set, the 40 m parent and the
D = 1000 m patch had been run to the steady verdict (103 and 95 iterations); the 20 m reference had not been run in
this study and no mast had been sampled.

| | Target |
|---|---|
| A1 | D = 1400 and 2000 m: RMS over the 30 masts of the relative difference of the 10 m speed from the 20 m run <= 1/2 of the 40 m parent's |
| A2 | The RMS does not grow from D = 1000 to 1400 to 2000 m (5 % allowance) |
| A3 | Every patch reaches the steady verdict, and the log prints no net-flux warning (net <= 1e-3 of the flux entering) |

A target that is missed is reported as missed, with the number.

## Expected Results
Measured (macOS arm64, AppleClang, double precision, 4 ranks). Every run reached the steady verdict. The differences
are of the 10 m speed at the 30 masts, relative to the 20 m run's.

| Run | Cells | Iterations (level 0) | RMS (%) | Largest (%) | Net flux / flux entering |
|---|---|---|---|---|---|
| 20 m, whole domain (the reference) | 1 620 000 | 95 | | | |
| 40 m, whole domain (the parent) | 405 000 | 103 | 1.208 | 3.590 | |
| Patch, D = 1000 m | 180 000 | 95 | 0.076 | 0.222 | -4.9e-6 |
| Patch, D = 1400 m | 352 800 | 95 | 0.035 | 0.083 | -1.3e-6 |
| Patch, D = 2000 m | 720 000 | 107 | 0.015 | 0.030 | -1.2e-7 |

| | Result |
|---|---|
| A1 | Met. 0.035 and 0.015 % against the parent's 1.208 % (1/34 and 1/78) |
| A2 | Met. 0.076, 0.035, 0.015 % |
| A3 | Met. No net-flux warning; the ground blend moved the patch's ground on the faces by up to 1.8, 2.4 and 0.6 m (D = 1000, 1400, 2000 m) |

With the nest faces 130 m beyond the outermost mast in x (D = 1000 m) the patch is within 0.08 % RMS of the uniform 20 m
run, on 11 % of its cells. One mast (line A, -841 m) is then in the ground blend, where the patch's ground differs
from the reference's by up to 0.05 m; without it the RMS is 0.077 %. In cell-iterations, parent included and the
coarse sequencing level left out, the three patches cost 0.38, 0.49 and 0.77 of the reference. Wall times are not reported: the runs shared the machine.

The reference here is GUST's own 20 m run. That run against ERF's URANS and its relation to the measured masts are in
`Askervein_210` and `docs/results.rst`, "Askervein on ERF's mesh"; a patch within 0.08 % of it does not change them.
The numbers are therefore the nesting error, not accuracy; the limits listed in `Nest_Validation/README.md` apply
(neutral flow, 2 2 1, one wind direction, centred square patches). Iteration error: the 20 m run against a run to
six orders differs by 3e-7 RMS at the masts.

Two-way (`two_way`; macOS arm64, AppleClang, 4 ranks, the recommended setting). RMS against the 20 m run; the nest
against the reference at its own cell heights, the outer domain against the reference's 2 x 2 means outside the
footprint. "One way" is a separate `nest.feedback = none` run (the nest's own coordinate surfaces). The joint run was
steady after 156 exchanges.

| Field | Nest, one way | Nest, joint levels | Nest, joint | Outer, one way | Outer, joint |
|---|---|---|---|---|---|
| u | 3.372e-3 | 3.92e-3 | 3.372e-3 | 4.62e-3 | **3.21e-3** |
| v | **5.68e-3** | 6.72e-3 | 6.03e-3 | 9.37e-3 | **6.93e-3** |
| w | 2.87e-3 | 2.84e-3 | **1.48e-3** | 2.57e-3 | **2.51e-3** |
| k | **1.56e-3** | 2.80e-3 | 2.76e-3 | 3.73e-3 | **3.12e-3** |

"Joint levels" is the first pass of the joint run: one way, but with the nest levels of the joint run (they follow
the outer domain's surfaces aloft; they move by 0.10 m RMS, at most 1.3 m). The levels alone raise the nest's u, v
and k errors; the feedback then lowers all four fields, so the nest's u ends where the one-way run left it
(tracker N-114 item 9; `docs/boundary_conditions.rst`, Two-way nests, Results). The feedback lowers the errors within
4 nest cells of the faces (u by 31 %, v by 32 %, w by 62 %) and raises u and v by up to 14 % 8 to 40 cells in. Every
nest error here is at most 0.2 % of the field's RMS (u 8.8 m/s, v 15.2 m/s, k 1.66 m²/s²).

With the conservative faces (`nest.feedback_flux = nest`, added to `two_way`'s joint run by hand; 156 exchanges)
the nest's RMS is 2.98e-3 (u), 4.74e-3 (v), 1.21e-3 (w) and 2.65e-3 (k), and the outer domain's outside the
footprint 3.19e-3, 6.73e-3, 2.44e-3 and 3.05e-3.

Horizontal speed in the lowest 200 m: no nest cell differs from the 20 m run by more than 0.25 m/s in either run
(largest 0.23 m/s one way, 0.14 m/s joint); outside the footprint 24 outer cells one way and 19 joint differ by more
than 0.25 m/s, 1 by more than 0.5 m/s. 10 m speed at the masts, relative to the 20 m run: nest 0.04 % one way and
0.11 % joint, outer domain 1.18 % and 1.05 % (RMS).

With `relax inject` (same machine and setting; relax steady after 191 exchanges, inject after 2153):

| Field | Nest, relax | Nest, inject | Outer, relax | Outer, inject |
|---|---|---|---|---|
| u | 3.53e-3 | 9.92e-3 | 3.75e-3 | 1.37e-2 |
| v | 6.33e-3 | 3.45e-2 | 7.56e-3 | 2.46e-2 |
| w | 1.54e-3 | 2.72e-3 | 2.50e-3 | 2.74e-3 |
| k | 2.67e-3 | 4.00e-3 | 3.16e-3 | 4.06e-3 |

Inject is the furthest of the four runs from the 20 m run in every field, nest and outer, except the nest's w (one
way 2.87e-3), and takes 14 times the exchanges of `defect`. No nest cell differs by more than 0.25 m/s (largest 0.12 m/s relax, 0.15 m/s inject); outside
the footprint 19 cells do with both. 10 m speed at the masts: nest 0.09 % (relax) and 0.22 % (inject), outer
1.06 % and 1.06 % (RMS). The RMS of the speed-up against the observations on lines A and AA is 0.281 and 0.217
(relax), 0.283 and 0.220 (inject) for the outer domain, against 0.281 and 0.217 for `defect`.

## CTest
`rt_askervein_nest_study` (fixture `rt_askervein_nest_clean`; `nightly`, four ranks, double precision only): the
study above; targets A1 to A3 are the verdict. The `nightly` label is outside the `unit|smoke` set of the
pull-request CI, so the decks and the script of this case are not run there.
