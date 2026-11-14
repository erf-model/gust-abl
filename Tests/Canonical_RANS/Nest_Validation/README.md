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
| `inputs_hill_domain` | The whole domain: 3200 m x 2400 m x 800 m, a 150 m round Witch of Agnesi hill of half-width a = 250 m at the centre, 16 layers (first 10 m), periodic in y, log-law inflow (u* 0.4 m/s, z0 0.1 m) from 20 degrees, neutral, SIMPLEC to the steady verdict. The horizontal cells come from the script (100, 50, 25 or 12.5 m). |
| `inputs_hill_nest` | A patch: four `nest` faces, first guess from the parent. Place, cells, parent and zone come from the script. |
| `nest_validation.py` | Runs a whole domain (`run`) or a patch (`nest`), compares a patch with the reference (`compare`), or does the whole study (`study`). Pure Python. |

## Data
None. The case is analytic.

## How to run
```bash
ctest --test-dir <build> -R "nest_validation"
```
The study, in a scratch directory (about 600 000 cells in the largest run):
```bash
python3 nest_validation.py study -- mpiexec -n 4 <build>/gust
python3 nest_validation.py study --stages one,zone -- mpiexec -n 4 <build>/gust
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
U_ref = (u*/kappa) ln(h/z0) = 7.313 m/s, the inflow speed at the hill's height h = 150 m, and the difference of k as a
fraction of U_ref^2. Both are given as RMS and largest value over four regions:

| Region | Cells |
|---|---|
| hill | within a of the hill's centre, layers centred below 100 m over flat ground |
| lee | the square of half-width a centred 2 a downwind of the centre, the same layers (inside patches of D >= 3 a) |
| faces | the two rings of patch cells beside the nest faces, all layers |
| patch | every patch cell |

The parent is compared with the reference averaged over each parent cell, over the same regions.

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
| T5 | T1 holds for D = 600 m on 1, 2 and 4 ranks, and the three hill RMS values agree to 10 % | before the runs |
| T6 | Zone of 3 and 5 cells (factor 1) at D = 600 m: no target. Reported: the four regions and the iterations. Expected from `docs/boundary_conditions.rst`: no gain over the hill | before the runs |

A target that is missed is reported as missed, with the number.

## Expected Results
To be filled from the study.

## CTest
To be filled.
