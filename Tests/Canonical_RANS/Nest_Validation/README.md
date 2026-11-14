# Nest_Validation

## Purpose
How close a one-way nested patch comes to a uniform fine run of the whole domain, as a function of the distance from
the hill to the patch's `nest` faces (sub-phase 9.3, beyond plan rev. 6). The reference is the whole domain on the
patch's cells; the patch is a square of half-width D centred on the hill, at refinement 2 2 1 of its parent
(`init.first_guess = nest`, four `nest` faces; `Tests/Canonical_RANS/Nest_Faces` tests the faces themselves). Three
questions: one nest from a 50 m parent; a chain of three nests from a 100 m run (ratio 8); hard nest faces against the
relaxation zone (`bc.nest_relax_width`). The same study on measured terrain is `Tests/Real_Terrain/Askervein_Nest`.

## Files
| File | Role |
|---|---|
| `inputs_hill_domain` | The whole domain: 3200 m x 2400 m x 800 m, a 150 m round Witch of Agnesi hill of half-width a = 250 m at the centre, first layer 10 m, periodic in y, log-law inflow (u* 0.4 m/s, z0 0.1 m) from 20 degrees, neutral, SIMPLEC to the steady verdict. The cells come from the script (100, 50, 25 or 12.5 m; 16 layers). |
| `inputs_hill_nest` | A patch: four `nest` faces, first guess from the parent. Place, cells, parent and zone come from the script. |
| `nest_validation.py` | Runs a whole domain (`run`) or a patch (`nest`), compares a patch with the reference (`compare`), or does the whole study (`study`). Pure Python. |

## Data
None. The case is analytic.

## How to run
```bash
ctest --test-dir <build> -R "nest_validation"
```
The study, in a scratch directory (786 432 cells in the largest run):
```bash
python3 nest_validation.py study -- mpiexec -n 4 <build>/gust
python3 nest_validation.py study --stages one,zone -- mpiexec -n 4 <build>/gust
python3 nest_validation.py study --stages ranks -- mpiexec -n 4 <build>/gust
```
Single runs:
```bash
python3 nest_validation.py run ref25 25 -- mpiexec -n 4 <build>/gust inputs_hill_domain
python3 nest_validation.py run par50 50 -- mpiexec -n 4 <build>/gust inputs_hill_domain
python3 nest_validation.py nest one_d600 par50 600 -- mpiexec -n 4 <build>/gust inputs_hill_nest
python3 nest_validation.py compare one_d600 ref25 par50 --rms-hill 1e-3
```

## Measures
The difference between two runs on the same cells is the length of the velocity-vector difference as a fraction of
the speed scale U_ref = (u*/0.4) ln(h/z0) = 7.313 m/s with h = 150 m the hill's height (the inflow speed at h is
7.14 m/s: the solver's kappa is 0.41), and the difference of k as a fraction of U_ref^2. Both are given as RMS and largest value over four regions:

| Region | Cells |
|---|---|
| hill | within a of the hill's centre, layers centred below 100 m over flat ground |
| lee | the square of half-width a centred 2 a downwind of the centre, the same layers (inside patches of D >= 3 a) |
| faces | the two rings of patch cells beside the nest faces, all layers |
| patch | every patch cell |

The parent is compared with the reference averaged over each parent cell, over the same regions and at the same
layer. Its ground is the coarser one (cell centres up to 1.4 m from the averaged reference's for the 50 m parent, 5 m
for the 100 m run), so its difference includes its ground. For a parent, faces is its cells centred within the two
rings, or within one of its own cells of the faces when that is wider. Patch and reference have the same ground and
heights on the same cells.

## Targets
Fixed in this file before the runs they judge, except where stated.

A pilot came first: one nest of D = 400, 600 and 1000 m from the 50 m parent against the 25 m reference, on a 60 m
hill (D = 600 m only) and on the 150 m hill. On the 60 m hill the parent was 0.23 % of U_ref from the reference over
the hill and the patch 0.010 %; the hill was made 150 m so that the parent's error is large enough to measure a
patch against (parent 0.70 %; patch 0.082, 0.019 and 0.008 %). Targets T1 and T2 were set with those numbers known.
T3 to T6, D = 800 m, the chain, the zone, the rank counts and all of `Askervein_Nest` had not been run.

| | Target | Set |
|---|---|---|
| T1 | One nest, D >= 600 m (2.4 a): hill velocity RMS <= 0.1 % of U_ref and <= 1/4 of the parent's | after the pilot |
| T2 | One nest: hill velocity RMS does not grow from D = 400 to 600, 800 and 1000 m (5 % allowance) | after the pilot, D = 800 m not run |
| T3 | Chain 100 -> 50 -> 25 -> 12.5 m (D = 1000, 800, 600 m): hill velocity RMS of the last patch against the 12.5 m reference <= 1/2 of that of the 100 m run, and <= 0.5 % of U_ref | before any chain run |
| T4 | Chain: every patch reaches the steady verdict, with a patch as the parent of the next | before any chain run |
| T5 | T1 holds for D = 600 m on 1, 2 and 4 ranks, and the three hill RMS values agree to 10 % | before the runs; run again in four boxes after review (below) |
| T6 | Zone of 3 and 5 cells (factor 1) at D = 600 m: no target. Reported: the four regions and the iterations. Expected from `docs/boundary_conditions.rst`: no gain over the hill | before the runs |

A target that is missed is reported as missed, with the number.

## Expected Results
Measured (macOS arm64, AppleClang, double precision, 4 ranks). Every run reached the steady verdict. All numbers are
percent of U_ref (velocity) or of U_ref^2 (k).

Whole-domain runs: 100 m 12 288 cells, 102 iterations; 50 m 49 152 cells, 107; 25 m 196 608 cells, 91; 12.5 m
786 432 cells, 113. Boxes: at most 64 x 64 cells for a whole domain and 48 x 48 for a patch, so the 100 m and 50 m
runs and the patches of up to 48 x 48 cells (D = 400 and 600 m, the zone runs, the chain's 50 m patch) were one box,
on one of the four ranks; the others had four boxes or more.

One nest from the 50 m parent, against the 25 m reference:

| D (m) | D / a | Cells | Iterations | Hill RMS | Hill max | Lee RMS | Faces RMS | Faces max | Hill k RMS |
|---|---|---|---|---|---|---|---|---|---|
| (the 50 m parent) | | 49 152 | 107 | 0.700 | 3.226 | 0.541 | 0.133 | 1.826 | 0.0338 |
| 400 | 1.6 | 16 384 | 124 | 0.082 | 0.173 | 0.268 | 0.126 | 1.211 | 0.0019 |
| 600 | 2.4 | 36 864 | 131 | 0.019 | 0.059 | 0.247 | 0.125 | 1.391 | 0.0006 |
| 800 | 3.2 | 65 536 | 125 | 0.011 | 0.025 | 0.102 | 0.128 | 1.551 | 0.0004 |
| 1000 | 4.0 | 102 400 | 122 | 0.008 | 0.013 | 0.031 | 0.131 | 1.690 | 0.0003 |

(The parent's lee and faces numbers are over the D = 1000 m patch's regions. The lee region reaches 720 m downwind
of the centre, so the patches of D = 400 and 600 m hold only part of it, ending on their nest faces.)

Chain of three from the 100 m run, against the 12.5 m reference:

| Run | Cells | Iterations | Hill RMS | Hill max | Lee RMS |
|---|---|---|---|---|---|
| The 100 m run | 12 288 | 102 | 3.110 | 7.334 | 1.486 |
| Patches of 50 m (D = 1000 m), 25 m (800 m), 12.5 m (600 m) | 25 600, 65 536, 147 456 | 128, 121, 106 | 0.044 | 0.092 | 0.172 |

Zone at D = 600 m (factor 1):

| `bc.nest_relax_width` | Iterations | Hill RMS | Lee RMS | Faces RMS | Faces max |
|---|---|---|---|---|---|
| 0 (hard faces) | 131 | 0.019 | 0.247 | 0.125 | 1.391 |
| 3 | 122 | 0.020 | 0.247 | 0.125 | 1.391 |
| 5 | 126 | 0.020 | 0.247 | 0.125 | 1.390 |

Targets:

| | Result |
|---|---|
| T1 | Met. D = 600, 800, 1000 m: 0.019, 0.011, 0.008 % (limit 0.1 %), 1/36 to 1/88 of the parent's 0.700 % |
| T2 | Met. 0.082, 0.019, 0.011, 0.008 % |
| T3 | Met. 0.044 % against the 100 m run's 3.110 % (1/71) |
| T4 | Met. 128, 121 and 106 iterations |
| T5 | Met. In four boxes of 24 x 24 cells (stage `ranks`) on 1, 2 and 4 ranks: 127 iterations and 0.0195 % on each; the fields of 1 and 4 ranks agree to 7e-14 of U_ref. The first T5 runs had the patch in one box at every rank count and tested no decomposition |
| T6 | No target. The zone changes the hill RMS by 0.001 % (0.0195, 0.0202, 0.0204 %). The iterations (-9 and -5) differ in the tail set by the theta residual (below) and say nothing about the zone |

What the numbers say:
- Over the hill the patch is within 0.02 % of the uniform fine run once its faces are 2.4 a from the hill's centre.
- In the two rings of cells beside the nest faces the patch keeps the parent's difference from the reference (0.13 %
  RMS, 1.2 to 1.7 % at the largest point, where the hill's wake leaves the patch). A region of interest must not end on
  a nest face: the lee region is within 0.03 % only at D = 1000 m, with the downwind face 1.1 a beyond it.
- To the steady verdict a patch takes more iterations (122 to 131) than the uniform fine run from the inflow profile
  (91). The verdict is relative to each run's own first residuals, which the parent's field lowers (u: 0.015 against
  0.37), and in these neutral runs the patch's last equation to pass is theta at its floor. The D = 600 m patch
  reaches the fine run's final momentum residuals between iterations 70 and 80 and stops with a u residual 25 times
  lower; the two converge at the same rate. The patch is cheaper through its cells. In cell-iterations,
  parent included: D = 600 m 0.56 of the 25 m reference, D = 1000 m 0.99 (this domain is only 1.6 widths of that patch long);
  the chain 0.32 of the 12.5 m reference.
- The chain's 0.044 % is 2.3 times the single nest's 0.019 % at the same D. No single 12.5 m nest from a uniform 25 m
  parent was run, so what the chain adds is not separated from the finer cells.
- Wall times are not reported: the runs shared the machine.

Limits of the study:
- The reference is the same solver with the same closure: the numbers are the nesting error, not accuracy.
- Neutral flow, refinement 2 2 1 with the same layers in every run, one wind direction, one hill with a = 5 parent
  cells. Whether the face distance scales with the hill or with the parent's cells is not separated.
- The flow is attached (smallest streamwise speed in the 25 m run 1.6 m/s): no recirculation reaches a nest face.
- Parent and reference share the outer boundaries; no patch face is near one. The patches are centred squares.
- The zone was tried with factor 1 at D = 600 m only, where the faces' difference does not reach the hill.
- Iteration error: the 25 m reference against a run to seven orders differs by 7e-7 RMS (1e-5 at the largest point)
  of U_ref over the hill, below every number above.

## CTest
- `rt_nest_validation_parent`, `_reference`, `_d400`, `_d600`, `rt_nest_validation_hill` (fixture
  `rt_nest_validation_clean`; `regression;smoke`, double precision only): the study one level coarser. The 100 m
  parent and the 50 m reference in four boxes on two ranks, 50 m patches of D = 400 m (four boxes, two ranks) and
  600 m (six boxes of 12 x 8 cells, three ranks). `_hill` requires the 600 m patch's hill RMS to be at most 0.2 % of U_ref (measured 0.091 %),
  at most 1/4 of its parent's (2.362 %), and not above the 400 m patch's (0.348 %).
- `rt_nest_validation_study` (fixture `rt_nest_validation_study_clean`; `nightly`, four ranks): the study above;
  targets T1 to T5 are the verdict. The `nightly` label is outside the `unit|smoke` set of the pull-request CI.
