# Building_Cube

## Purpose
One block (30 x 30 m, 20 m high) on flat ground in a log-law wind: the fixed rows, the masked pressure correction and
the form 7 wall law on the roof and the sides (`docs/buildings.rst`) on one level. It pins that the iteration converges
and that the flow keeps out of the block, and it runs the form without the sides' law (`buildings.wall_form =
roof_stress`) against that form's values. It is not a validation; PR-44 compares with measurements.

## Files
| File | Role |
|---|---|
| `inputs_cube` | 48 x 24 x 24 cells of 5 m, inflow at xlo, outflow at xhi, periodic in y, no gravity. |
| `cube.txt` | The block in ERF's format. |
| `inputs_abort_newton`, `inputs_abort_lusgs` | The deck with a stage or sweep that has no buildings yet: refused. |
| `plot_cube.py` | Writes `docs/figures/buildings/block_flow.png`. Not a test. |
| `check_cube.py` | Speeds in the solid, the volume flux along x, reversal behind and speed-up above, and u in a wall cell of the roof and of a side (`--roof`, `--side`). The exit code is the verdict. |

## Data and acknowledgements
No external data: the map is made up for the test. The wall law follows ERF's `erf.if_wall_form = fraction_stress`
and Kynema-SGF's wall_form 7.

## How to run
```bash
ctest --test-dir <build> -R building_cube --output-on-failure
```

## Expected Results
The run stops at the driver's stopping rule. Speed in solid cells is 0 and the plane flux is constant to four digits.
k has no building terms yet, so these values are not validated.

| `buildings.wall_form` | Iterations | u behind the block, near the ground | u 20 m above the roof (5.91 m/s upstream) | u in the roof's wall cell | u in a side's wall cell |
|---|---|---|---|---|---|
| `fraction_stress` (default) | 111 | -0.77 m/s | 6.30 m/s | 3.833 m/s | 2.387 m/s |
| `roof_stress` | 114 | -0.55 m/s | 6.36 m/s | 3.888 m/s | 1.505 m/s |

Measured on 1, 2 and 4 ranks. The tests pin the two wall-cell values and the signs of the others; the iteration counts are
not pinned. The `roof_stress` plotfile equals the one of `main` before the sides had the law (bba99be) at zero tolerance.

## CTest
`building_cube_clean`, `building_cube_run`, `building_cube_flow`, `building_cube_roof_stress_clean`,
`building_cube_roof_stress_run`, `building_cube_roof_stress_flow`, `abort_buildings_newton`, `abort_buildings_lusgs`,
`abort_buildings_levels`, `abort_buildings_mlmg`.
