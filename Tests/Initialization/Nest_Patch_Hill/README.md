# Nest_Patch_Hill

## Purpose
The first guess of a nested patch from its parent run (`init.first_guess = nest`, plan rev. 6, sub-phase 9.1).
The parent is the mass-consistent initial field over a 60 m round Witch of Agnesi hill on 50 m cells; the patch
covers the hill on 25 m cells (refinement 2 2 1) and starts from the parent's checkpoint interpolated per
zeta-level, then the mass-consistent adjustment. The interpolated field is checked against an interpolation of the
parent's plotfile computed independently in Python.

## Files
| File | Role |
|---|---|
| `inputs_parent` | The parent: 32 x 24 x 16 cells, 1600 m x 1200 m x 800 m, periodic in y, inflow from 20 degrees, stable sounding. Writes `nest_parent_chk00000`. |
| `inputs_patch` | The patch: 32 x 24 x 16 cells from (400, 300) m to (1200, 900) m, i.e. parent cells 8 to 23 in x and 6 to 17 in y, with `bc.nest_parent = nest_parent_chk00000`. |
| `inputs_patch_misplaced` | The patch shifted 10 m in x and with 20 cells in y: it must stop, naming `bc.nest_parent`. |
| `check_nest.py` | Runs the patch twice and checks it (below). |

## CTest
Double precision only.
| CTest | What it checks |
|---|---|
| `rt_nest_patch_hill_clean`, `rt_nest_patch_hill_parent` | Fixtures: remove old outputs; run the parent on 2 ranks. |
| `rt_nest_patch_hill` | `check_nest.py nest_parent_plt00000 8 6 1e-10 -- gust inputs_patch` on 3 ranks. |
| `abort_nest_placement` | The misplaced patch stops, naming `bc.nest_parent` and the faults. |
| `abort_nest_inputs`, `abort_nest_init_keys`, `abort_nest_parent_missing`, `abort_nest_parent_unused`, `abort_nest_parent_not_a_checkpoint` | Input faults around `init.first_guess = nest` stop, naming the key. |

## Expected Results
- Without the adjustment (`init.projection=0`), theta and k of every patch cell equal the bilinear interpolation, per
  zeta-level, of the parent's cell values at the patch cell centre, to a relative 1e-12 (measured 3.8e-16 and 4.4e-16).
  Both fields vary along a level over the hill, so a wrong position or level would be seen.
- The log names the parent and the parent cells covered (8 to 23 in x, 6 to 17 in y).
- No mass crosses the ground: on every ground face of `nest_patch_chk00000`, |Omega| = |w - z_x u - z_y v| (the code's
  averages, ghost row included) is at most 1e-12 x max |w| there (measured 4.7e-16 against 1.52 m/s).
- With the adjustment, the divergence after is at most 10 x RELTOL x the divergence before, with RELTOL = 1e-10
  given to `check_nest.py` (the deck's `solver.poisson_reltol`); measured 0.090 before, 9.5e-13 after.

## How to run
```bash
ctest --test-dir <build> -R "nest"
```
