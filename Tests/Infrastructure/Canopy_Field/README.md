# Canopy_Field

## Purpose
The forest canopy's drag field (`canopy.*`, `docs/canopy.rst`) as the application builds and writes it. Four
ERF-format patches on a stretched, periodic-in-x mesh of 12 boxes on two ranks: two that overlap across box edges,
one across the periodic seam, and a dense top-heavy one, with both leaf area density profiles. The check uses a closed
form: each column's integral of the drag field is the drag coefficient times the leaf area index of the patches over
it. A restart with the same canopy continues bit for bit; a restart with another canopy, or without one, stops; a deck
of canopy mistakes stops with every input named.

The drag field does not act on the flow yet, so no stage runs here. The column integral does not see how the drag
is spread over the layers; the profile shapes are checked by the unit tests (`Canopy.*`), which also cover
`canopy.type = map`.

## Files
| File | Role |
|---|---|
| `inputs_canopy` | The 2-rank run; writes `canopy_plt00000` (with `canopy_drag`) and `canopy_chk00000`. `canopy.file` is given on the command line. |
| `canopy_patches.txt` | The patches, ERF's `type x y height diameter cd lai laimax`. |
| `canopy_patches_changed.txt` | The same with one LAI changed (3.0 to 3.5). |
| `check_canopy.py` | Integrates `canopy_drag * detJ * dzeta` over each column of the plotfile and compares with `cd * lai` of the patches over the column centre (nearest periodic image in x); columns outside every patch must hold 0. |
| `inputs_restart` | Restart on another layout and rank count with the same canopy; writes `canopy_restart_plt00000`. |
| `inputs_abort_restart_no_canopy` | Restart without the canopy: must stop. |
| `inputs_abort_canopy` | A map without a drag coefficient, a missing height file, no LAI file, `laimax` with the uniform profile, the patch file with a map, an unknown key: must stop, naming each. |

## Data and acknowledgements
No external data: the patches are made up for the test. The patch format is ERF's `erf.forest_file`
(`Source/SourceTerms/ERF_ForestDrag.cpp`, erf-model/ERF).

## How to run
```bash
ctest --test-dir <build> -R "^canopy_|^abort_canopy_"
```

## Expected Results
- `check_canopy.py`: 22 columns with canopy; the largest column error is at rounding level (3e-16 of the largest
  `cd * lai` in double precision; the check allows 5e-6 for a single-precision plotfile).
- The same check against `canopy_patches_changed.txt` fails (the check is not vacuous).
- The restart log shows `restarted from canopy_chk00000 at step 0`, and `canopy_plt00000` and
  `canopy_restart_plt00000` agree under `amrex_fcompare` at zero tolerance (`--allow_diff_grids`).
- The changed canopy stops with "the canopy inputs differ from the checkpoint's"; the restart without a canopy with
  "the checkpoint was written with a canopy".

## CTest
`canopy_clean`, `canopy_run`, `canopy_column_integrals`, `canopy_check_detects_change` (expected to fail),
`canopy_restart`, `canopy_restart_matches`, `abort_canopy_restart_changed`, `abort_canopy_restart_none`,
`abort_canopy_inputs`. All but the last share the run directory `<build>/Tests/run/canopy`.
