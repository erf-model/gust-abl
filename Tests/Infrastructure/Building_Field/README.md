# Building_Field

## Purpose
The building geometry (`buildings.*`, `docs/buildings.rst`) in a run without a stage: the plotfile holds
`building_beta` and `building_wall`, and each column's solid volume equals the height map's mean over the column. A
restart on another layout gives the same plotfile, and a restart with another map or without buildings is refused.

## Files
| File | Role |
|---|---|
| `inputs_buildings` | 32 x 24 x 20 cells, 8 x 8 boxes, stretched layers, periodic in x. `buildings.file` is given on the command line. |
| `inputs_restart` | Restart on 16 x 12 boxes with the same buildings. |
| `inputs_restart_none` | Restart on the deck's boxes without buildings: refused. |
| `building_map.txt` | ERF format on 20 m points: a 30 m block, and a 17.5 m block across the periodic seam. |
| `building_map_changed.txt` | The same with a 33 m block: the check must fail on it, and a restart must refuse it. |
| `plot_building_kinds.py` | Runs GUST on a map of building kinds, flat and on a hill, and writes `docs/figures/buildings/`. Not a test. |
| `check_buildings.py` | The column volumes against the map; the exit code is the verdict. |

## Data and acknowledgements
No external data: the maps are made up for the test. The file format is ERF's `erf.buildings_file_name`
(`read_custom_terrain` in `Source/ERF_ProbCommon.H`, erf-model/ERF).

## How to run
```bash
ctest --test-dir <build> -R buildings --output-on-failure
```

## Expected Results
`buildings_column_volume` reports a largest column error below 1e-5 of the tallest building and pins the mask of two
columns. `buildings_restart` logs `restarted from buildings_chk00000 at step 0`, and `buildings_restart_matches` finds no
difference. The `abort_*` tests each match the refusal's message.

## CTest
`buildings_clean`, `buildings_run`, `buildings_column_volume`, `buildings_check_detects_change` (fails on purpose), `buildings_restart`,
`buildings_restart_matches`, `abort_buildings_restart_changed`, `abort_buildings_restart_none`, `abort_buildings_inputs`,
`abort_buildings_type`, `abort_buildings_wall_inputs`, `abort_buildings_wall_tke_form`, `abort_buildings_z0_with_drag` (also `buildings.wall_speed` with `drag`),
`buildings_defaults_fraction_stress`, `buildings_defaults_roof_stress`, `buildings_defaults_drag` (each with its `_clean`).
