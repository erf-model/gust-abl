# Immersed_Ridge_2D

## Purpose
A sloping ground immersed in a flat mesh (`terrain.representation = immersed_forcing`,
[docs](../../../docs/immersed_terrain.rst)) against the same ground on a terrain-following mesh. Three pairs of
decks:

| Pair | Ground | Wind | Question |
|---|---|---|---|
| `along` | Witch of Agnesi ridge, 100 m high, half-width 200 m, slopes up to 18° | along the ridge (y), driven by a force of 1.225e-4 m/s² | the wall law on a slope, with no flow over the ridge |
| `ridge` | The ridge with half-width 500 m, slopes up to 7° | across the ridge (x), geostrophic 10 m/s with Coriolis | crest speed-up and lee minimum |
| `hill` | Round Agnesi hill, 100 m high, half-width 300 m, slopes up to 12° | as `ridge` | the same in 3-D |

All are neutral (uniform θ, no gravity), periodic in x and y, with z0 = 0.1 m and a slip top. The body-fitted mesh
has 64 layers of 15.625 m far from the ground (1000 m). The immersed mesh has the same layers and two more below
the terrain's zero level (`terrain.flat_height` = 31.25 m, 66 layers): the cell count is 3 % higher, and on the
crest the immersed layers are 15.6 m against 14.1 m body-fitted. Ridges: 128 x 1 columns of 40 m. Hill: 32 x 32
columns of 80 m.

Every deck uses the recommended turbulence setting (docs/turbulence_model.rst: `turb.realizability = 0.309`,
`turb.nut_relax = 0.5`, `turb.max_geom_lscale = 100`, `convection.hybrid_sigma_max = 0.1`,
`convection.hybrid_height = 250`) and the default face mean, so a pair differs in the ground's representation only.
The immersed decks set `solver.poisson_backend = gmres_fft` (the multigrid backend is refused with an immersed
ground). The stop rule is four orders of residual drop.

A tilted flat channel was considered first. A uniform slope cannot be periodic, and open faces were not
supported with an immersed ground when these decks were written (they are now: `../Immersed_Open_Faces`), so the
wind along a ridge takes its place: the ground slopes across the wind,
and there is no flow over it.

## Files
| File | Role |
|---|---|
| `inputs_body_fitted_along`, `inputs_immersed_along` | Wind along the ridge. |
| `inputs_body_fitted_ridge`, `inputs_immersed_ridge` | Wind across the ridge. |
| `inputs_immersed_ridge_lusgs` | `inputs_immersed_ridge` with LU-SGS (`steady.sweep = lusgs`) in place of z-line Jacobi. |
| `inputs_body_fitted_hill`, `inputs_immersed_hill` | The round hill. |
| `inputs_hybrid_ground`, `inputs_hybrid_restart` | The hybrid mesh on the ridge with the default `terrain.flat_height` and `terrain.hybrid_half_width`: two steps that write a checkpoint, then a restart from it on another box layout. |
| `check_immersed_ridge.py` | Runs an immersed deck, runs (or reads) its body-fitted reference, and compares the two. `--variant NAME` picks `inputs_immersed_<case>_NAME`; `--reference-immersed DECK --sweep PCT` also runs a second immersed deck and bounds the largest difference of the wind between the two. |

## Data and acknowledgements
None: the cases are analytic.

## How to run
```bash
ctest --test-dir <build> -R "^rt_immersed_(ridge|along)"
ctest --test-dir <build> -R "^rt_immersed_hill"
```
By hand, in an empty folder:
```bash
python3 check_immersed_ridge.py --wind along --max 2 --rms 1.5 --body-fitted inputs_body_fitted_along -- mpirun -np 4 <build>/gust inputs_immersed_along
```
The body-fitted run goes to the sub-folder `body_fitted/`. `--reuse` (by hand; no CTest uses it) reads the plotfile already there.

## Expected Results
Both runs reach the steady verdict. The wind is read at d = 20, 40, 80 and 160 m above the ground in every column,
from the fluid cells only (`building_beta` = 0), interpolated linearly in ln((d + z0)/z0); a column's ground is
the mean of its nodes' heights. Differences are in percent of the body-fitted far-field speed at that height
(the mean of the first sixteenth of the columns), except the crest speed, which is relative to the body-fitted
crest speed.

**Wind along the ridge** (limits: largest difference 2 %, rms 1.5 %). Iterations: body-fitted 5298, immersed 3962.

| d (m) | Far field bf / immersed (m/s) | Largest difference | rms |
|---|---|---|---|
| 20 | 4.526 / 4.580 | 1.32 % | 0.79 % |
| 40 | 5.165 / 5.214 | 1.07 % | 0.62 % |
| 80 | 5.861 / 5.904 | 0.79 % | 0.47 % |
| 160 | 6.679 / 6.718 | 0.60 % | 0.36 % |

The immersed far field is 0.6 to 1.2 % above the body-fitted one, where the ground is flat. The 1-D column shows
the same offset with this turbulence setting (see below), so it is not an effect of the slope. On the ridge the
largest speed at 20 m is 4.832 m/s body-fitted and 4.834 m/s immersed.

**Wind across the ridge** (limits: crest 4 %, lee minimum 16 %, rms 8 %). Iterations: body-fitted 6297, immersed
4822 (re-measured on four ranks with the code of `rt_immersed_ridge_lusgs`; the table before gave immersed 4878 and
speeds up to 0.04 m/s lower, from an earlier version of the code).

| d (m) | Far field bf / immersed (m/s) | Crest bf / immersed (m/s) | Crest difference | Lee minimum bf / immersed (m/s) | Lee difference | rms |
|---|---|---|---|---|---|---|
| 20 | 4.863 / 4.664 | 7.411 / 7.298 | −1.5 % | 4.516 / 3.884 | −13.0 % | 6.4 % |
| 40 | 5.569 / 5.365 | 8.004 / 7.868 | −1.7 % | 5.401 / 5.000 | −7.2 % | 4.4 % |
| 80 | 6.399 / 6.203 | 8.529 / 8.357 | −2.0 % | 6.377 / 6.139 | −3.7 % | 2.9 % |
| 160 | 7.486 / 7.327 | 9.142 / 8.992 | −1.6 % | 7.484 / 7.326 | −2.1 % | 1.9 % |

The crest is in the same column in the two runs. The immersed run is 2 to 4 % slower in the far field and 1.5 to
2 % at the crest, and the lee minimum at 20 m is deeper (3.88 against 4.52 m/s) and four columns (160 m) nearer
the crest. The tables of the wind along the ridge and of the round hill were not re-measured.

**The round hill** (limits: crest 5 %, lee minimum 27 %, rms 13 %), read along the row of columns just north of
the hill's centre. Iterations: body-fitted 2441, immersed 2180.

| d (m) | Far field bf / immersed (m/s) | Crest bf / immersed (m/s) | Crest difference | Lee minimum bf / immersed (m/s) | Lee difference | rms |
|---|---|---|---|---|---|---|
| 20 | 4.871 / 4.672 | 7.242 / 6.957 | −3.9 % | 4.546 / 3.358 | −24.4 % | 11.2 % |
| 40 | 5.619 / 5.437 | 7.736 / 7.515 | −2.9 % | 5.520 / 4.899 | −11.1 % | 5.6 % |
| 80 | 6.514 / 6.373 | 8.211 / 8.019 | −2.3 % | 6.505 / 6.369 | −2.1 % | 2.4 % |
| 160 | 7.604 / 7.489 | 8.807 / 8.662 | −1.6 % | 7.603 / 7.489 | −1.5 % | 1.6 % |

**What the offsets are, and are not.** The geometry is exact for the piecewise-linear ground (unit tests
`ImmersedSurface.*`), and the wall stress balances the force on the fluid (`rt_immersed_column_*`).

- *Flat ground and the wind along the ridge (about 1 %).* The closure's central difference of a logarithmic
  profile over-estimates the shear in the first cells above the wall cell, on a body-fitted mesh too; above a cut
  cell the next cell lies closer to the ground and the error is larger (`Immersed_Column/README.md`; a
  log-consistent shear is a separate change, not in these decks). With the setting of these decks the 1-D
  force-driven column (ground 0.6 m above a cell face, as in the far field here) gives, against the log law of
  the exact friction velocity at the first full cell: body-fitted +0.4 %, immersed +1.8 %. With
  `turb.face_mean = logarithmic` the two are +2.1 % and +4.9 %; the decks therefore keep the default face mean in
  both runs. The geostrophic column of the `ridge` decks without the ridge differs by 0.5 % or less above the cut
  cell.
- *Flow over the ridge and the hill (1.5 to 4 % in the far field, more in the lee near the ground).* The immersed
  runs decelerate more in the lee. The slower far field is consistent with the periodic domain returning that
  flow upstream; this was not tested separately. The far field of the `ridge` decks is not flat either: the
  ridge's tail at the side of the domain is 3.7 m high (0.6 m in the `along` decks). The flat column (ground
  0.6 m above a cell face) does not show the difference. Its cause is not isolated here. Not tested: the
  log-consistent shear on a slope, the convection scheme at the cut faces, and finer layers (the lee minimum at
  20 m lies below the first full cell of the immersed mesh in part of the lee). The limits of the tests are the
  measured values with a margin: they detect a change, they do not assert agreement.
- A wall law along the ground's normal (distance along the normal, speed parallel to the ground) was also run on
  these decks. It differed from the vertical form by at most 0.3 % of the far-field speed (the slopes are 7° to
  18°, where cos α is 0.99 to 0.95) and was removed from the code (see the docs).

## CTest
`rt_immersed_along_vertical`, `rt_immersed_ridge_vertical` (label `regression`, four ranks). Each test runs the
immersed and the body-fitted deck in `<build>/Tests/run/rt_immersed_<pair>`. Double precision only.
`rt_immersed_hill_vertical` (label `nightly`, four ranks; about 40 min).
`rt_immersed_ridge_lusgs` (label `regression`, four ranks) runs `inputs_immersed_ridge_lusgs`, the body-fitted ridge
and `inputs_immersed_ridge` (three runs, about 2.5 min). It has the limits of `rt_immersed_ridge_vertical` against
the body-fitted run. At each height, u and w of the two immersed runs may differ by at most 2e-4 % of the body-fitted
far-field speed, and k by at most 2e-4 % of the body-fitted far-field k; the two runs' logs must name different
sweeps. Measured on four ranks: 4670 iterations with LU-SGS, 4822 with z-line Jacobi. The largest differences are
7.1e-7 m/s in u (1.5e-5 %), 2.0e-7 m/s in w (2.7e-6 %) and 1.1e-7 m²/s² in k (1.9e-5 %): the two sweeps reach the
same steady state. Mutants: with the refusal of LU-SGS on an immersed ground restored, the test fails at the first
run; `bc.surface_z0 = 0.101` in the LU-SGS deck fails the sweep check at every height (u differs by 0.12 % at
20 m); `steady.sweep = zjacobi` in the LU-SGS deck fails the check that the logs name different sweeps.
`immersed_advice_run` (two ranks) runs `inputs_immersed_advice` (`immersed_forcing` over the ridge, two
iterations) and checks that the warning recommending `hybrid` is printed; `immersed_advice_restart` (three ranks)
restarts from its checkpoint (`inputs_immersed_advice_restart`) and checks that it is not.
`hybrid_ground_run` (two ranks) checks that the hybrid defaults are applied and logged (half-width 80 m = 2 cells).
`hybrid_ground_restart` (three ranks, other boxes) restarts from its checkpoint. `abort_hybrid_restart_width` restarts with
`terrain.hybrid_half_width = 100` and must be refused because the immersed ground's inputs differ.

Mutant: with the wall cell's distance d1 halved (`Source/BuildingsSurface.cpp`), `rt_immersed_along_vertical`
fails at every height (largest difference 7.4 % at 20 m against the limit of 2 %).
