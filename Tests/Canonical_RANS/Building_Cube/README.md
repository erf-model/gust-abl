# Building_Cube

## Purpose
One block (30 x 30 m, 20 m high) on flat ground in a log-law wind: the fixed rows, the masked pressure correction and
the form 7 wall law on the roof and the sides and the closure's building terms (`docs/buildings.rst`) on one level. It
pins that the iteration converges and that the flow keeps out of the block. It also runs the form without the sides' law
(`buildings.wall_form = roof_stress`) and the one without the closure's terms (`buildings.wall_tke = off`) against those
forms' values, and the deck with the ground's wall k as a source (`bc.surface_wall_k = source`). It is not a validation; PR-44 compares with measurements.

## Files
| File | Role |
|---|---|
| `inputs_cube` | 48 x 24 x 24 cells of 5 m, inflow at xlo, outflow at xhi, periodic in y, no gravity. |
| `cube.txt` | The block in ERF's format. |
| `inputs_abort_newton`, `inputs_abort_lusgs` | The deck with a stage or sweep that has no buildings yet: refused. |
| `plot_cube.py` | Writes `docs/figures/buildings/block_flow.png`. Not a test. |
| `check_cube.py` | Speeds in the solid, the volume flux along x, reversal behind and speed-up above, and u in a wall cell of the roof and of a side and k in the roof's wall cell (`--roof`, `--side`, `--k-roof`). The exit code is the verdict. |

## Data and acknowledgements
No external data: the map is made up for the test. The wall law follows ERF's `erf.if_wall_form = fraction_stress`
and Kynema-SGF's wall_form 7.

## How to run
```bash
ctest --test-dir <build> -R building_cube --output-on-failure
```

## Expected Results
The run stops at the driver's stopping rule. Speed in solid cells is 0 and the plane flux is constant to four digits.
These values are not validated.

| Form | Iterations | u behind the block, near the ground | u 20 m above the roof (5.91 m/s upstream) | u in the roof's wall cell | u in a side's wall cell | k in the roof's wall cell |
|---|---|---|---|---|---|---|
| the deck (`fraction_stress`, `wall_tke = on`) | 105 | -0.76 m/s | 6.39 m/s | 2.801 m/s | 1.491 m/s | 0.411 m²/s² |
| `bc.surface_wall_k = source` | 105 | -0.69 m/s | 6.39 m/s | 2.843 m/s | 1.599 m/s | 0.911 m²/s² |
| `buildings.wall_tke = off` | 111 | -0.77 m/s | 6.30 m/s | 3.833 m/s | 2.387 m/s | 8.618 m²/s² |
| `buildings.wall_form = roof_stress` | 114 | -0.55 m/s | 6.36 m/s | 3.888 m/s | 1.505 m/s | 8.589 m²/s² |

The deck was measured on 1, 2 and 4 ranks, `source` on 2 and 4, the two earlier forms on 2. The tests pin the three wall-cell values, the signs
of the others and, for the deck, k at its floor (1e-6 m²/s²) in the solid cells;
the iteration counts are not pinned. The block has no cell with a fixed face among those of a wall's parallel
components, so `buildings.wall_speed` does not change it: the deck's plotfile equals `main`'s at c57ad06. The
`wall_tke = off` statement below holds for the same reason. The `wall_tke = off` plotfile equals the one of `main` before the closure saw the
buildings (8f1b605), and the `roof_stress` one that of `main` before the sides had the law (bba99be), at zero tolerance.

## CTest
`building_cube_clean`, `building_cube_run`, `building_cube_flow`, `building_cube_roof_stress_clean`,
`building_cube_roof_stress_run`, `building_cube_roof_stress_flow`, `building_cube_wall_tke_off_clean`,
`building_cube_wall_tke_off_run`, `building_cube_wall_tke_off_flow`, `building_cube_wall_k_source_clean`,
`building_cube_wall_k_source_run`, `building_cube_wall_k_source_flow`, `abort_buildings_newton`, `abort_buildings_lusgs`,
`abort_buildings_levels`, `abort_buildings_mlmg`.
