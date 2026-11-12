# Building_Cube

## Purpose
One block (30 x 30 m, 20 m high) on flat ground in a log-law wind: the fixed rows, the masked pressure correction and
the form 7 wall law on the roof and the sides and the closure's building terms (`docs/buildings.rst`) on one level. It
pins that the iteration converges and that the flow keeps out of the block. It also runs the form without the sides' law
(`buildings.wall_form = roof_stress`) and the one without the closure's terms (`buildings.wall_tke = off`) against those
forms' values, the deck with the ground's wall k as a source (`bc.surface_wall_k = source`), and the fixed rows alone
(`buildings.wall_form = drag`, where `buildings.wall_tke` is off when not given), walls 5 K warmer than the air with
gravity on (`buildings.theta_mode = fixed`), theta without its building terms (`buildings.theta_mode = none`), and walls smoother than the ground
(`buildings.z0 = 0.01`), where the closure's length beside a wall takes the wall's roughness. It is not a validation; PR-44 compares with measurements.

## Files
| File | Role |
|---|---|
| `inputs_cube` | 48 x 24 x 24 cells of 5 m, inflow at xlo, outflow at xhi, periodic in y, no gravity. |
| `inputs_cube_drag` | The same with `buildings.wall_form = drag` and no `buildings.z0` or `buildings.wall_tke`. |
| `inputs_cube_smooth_walls` | The same with `buildings.z0 = 0.01` (the ground has 0.1 m) and a turbulence plotfile. |
| `inputs_cube_heated` | The same with walls at 305 K (`buildings.theta_mode = fixed`) and gravity on. |
| `cube.txt` | The block in ERF's format. |
| `cube_24p6.txt` | The block 24.6 m high: its roof ends 4.6 m into a 5 m layer, so the roof's cells are 0.92 solid. |
| `inputs_abort_newton`, `inputs_abort_lusgs` | The deck with a stage or sweep that has no buildings yet: refused. |
| `plot_cube.py` | Writes `docs/figures/buildings/block_flow.png`. Not a test. |
| `check_cube_length.py` | The turbulence plotfile's `walldist` in three wall cells and far upstream: the wall's distance plus `buildings.z0` less the ground's roughness where a wall is the nearer. The exit code is the verdict. |
| `check_cube_heat.py` | The heated block: theta in the solid cells and in the air, and the heat the walls give (the log's line) against the heat carried through the planes downstream and against the deck's value. The exit code is the verdict. |
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
| the deck (`fraction_stress`, `wall_tke = on`) | 109 | -0.76 m/s | 6.39 m/s | 2.801 m/s | 1.491 m/s | 0.411 m²/s² |
| `bc.surface_wall_k = source` | 106 | -0.69 m/s | 6.39 m/s | 2.843 m/s | 1.599 m/s | 0.911 m²/s² |
| `buildings.wall_tke = off` | 111 | -0.77 m/s | 6.30 m/s | 3.833 m/s | 2.387 m/s | 8.618 m²/s² |
| `buildings.wall_form = roof_stress` | 114 | -0.55 m/s | 6.36 m/s | 3.888 m/s | 1.505 m/s | 8.589 m²/s² |
| `buildings.wall_form = drag` (`inputs_cube_drag`) | 115 | -0.54 m/s | 6.37 m/s | 2.136 m/s | 1.524 m/s | 7.283 m²/s² |
| walls at 305 K, gravity on (`inputs_cube_heated`) | 105 | -0.77 m/s | 6.38 m/s | 2.870 m/s | 1.556 m/s | 0.432 m²/s² |
| `buildings.z0 = 0.01` (`inputs_cube_smooth_walls`) | 116 | -0.76 m/s | 6.39 m/s | 2.855 m/s | 1.452 m/s | 0.145 m²/s² |

The deck was measured on 1, 2 and 4 ranks, `source` on 2 and 4, the two earlier forms, `drag` and the smooth walls on 2. The tests pin the three wall-cell values, the signs
of the others and, for the deck, k at its floor (1e-6 m²/s²) in the solid cells;
the iteration counts are not pinned. The block has no cell with a fixed face among those of a wall's parallel
components, so `buildings.wall_speed` does not change it. Since the solid cells' rows left the residual norm of k the
deck stops 4 iterations later than `main` at 8448d15 (105) and `source` one later; the deck's velocities differ from
that plotfile by at most 3e-5 m/s, and its building distances not at all (the block's walls lie on cell faces). The
`wall_tke = off` statement below holds for the same reason. The `wall_tke = off` plotfile equalled the one of `main` before the closure saw the
buildings (8f1b605) at zero tolerance, and the `roof_stress` one equals that of `main` before the sides had the law (bba99be).
Since theta has building terms (`buildings.theta_mode`, `ambient` by default with `fraction_stress`) the deck and
`wall_tke = off` keep their velocities and k bit for bit and their iteration counts (measured against
`buildings.theta_mode = none`, which is theta as on `main` at 5735108); theta differs by up to 0.002 K, its drift from
300 K inside the block without these terms.

The heated block (`inputs_cube_heated`, 2 ranks): the walls give the air 585.5 K m³/s, the planes downstream of the
block carry 583.4 K m³/s (577.7 to 585.2), and the air is at most 1.07 K above the inflow. With gravity off the walls
give 580.9 K m³/s and the wall-cell values are the deck's. With the roof inside a layer (`cube_24p6.txt`) the run
stops in 133 iterations, the walls give 842.7 K m³/s, the planes carry 839.1, and the air is at most 2.32 K above the
inflow: the exchange of a partly solid roof cell depends on where the roof lies in its layer (`docs/buildings.rst`,
Limits). In the heated deck k in the side's wall cell on the ground is 0.0089 m²/s², against 0.0151 with gravity off:
it moves with the buoyancy of the warmed air, not with the walls' terms.

## CTest
`building_cube_clean`, `building_cube_run`, `building_cube_flow`, `building_cube_roof_stress_clean`,
`building_cube_roof_stress_run`, `building_cube_roof_stress_flow`, `building_cube_wall_tke_off_clean`,
`building_cube_wall_tke_off_run`, `building_cube_wall_tke_off_flow`, `building_cube_wall_k_source_clean`,
`building_cube_wall_k_source_run`, `building_cube_wall_k_source_flow`, `building_cube_drag_clean`,
`building_cube_drag_run`, `building_cube_drag_flow`, `building_cube_smooth_walls_clean`,
`building_cube_smooth_walls_run`, `building_cube_smooth_walls_length`, `building_cube_smooth_walls_flow`,
`building_cube_heated_clean`, `building_cube_heated_run`, `building_cube_heated_heat`, `building_cube_heated_flow`,
`building_cube_theta_none_clean`, `building_cube_theta_none_run`, `building_cube_theta_none_flow`,
`building_cube_heated_roof_clean`, `building_cube_heated_roof_run`, `building_cube_heated_roof_heat`, `abort_buildings_theta_frozen`,
`abort_buildings_newton`, `abort_buildings_lusgs`,
`abort_buildings_levels`, `abort_buildings_mlmg`.
