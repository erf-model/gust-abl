# Building_Cube

## Purpose
One block (30 x 30 m, 20 m high) on flat ground in a log-law wind: the fixed rows, the masked pressure correction and
the form 7 roof law (`docs/buildings.rst`) on one level. It pins that the iteration converges and that the flow
keeps out of the block. It is not a validation; PR-44 compares with measurements.

## Files
| File | Role |
|---|---|
| `inputs_cube` | 48 x 24 x 24 cells of 5 m, inflow at xlo, outflow at xhi, periodic in y, no gravity. |
| `cube.txt` | The block in ERF's format. |
| `inputs_abort_newton`, `inputs_abort_lusgs` | The deck with a stage or sweep that has no buildings yet: refused. |
| `plot_cube.py` | Writes `docs/figures/buildings/block_flow.png`. Not a test. |
| `check_cube.py` | Speeds in the solid, the volume flux along x, reversal behind and speed-up above. The exit code is the verdict. |

## Data and acknowledgements
No external data: the map is made up for the test. The wall law follows ERF's `erf.if_wall_form = fraction_stress`
and Kynema-SGF's wall_form 7.

## How to run
```bash
ctest --test-dir <build> -R building_cube --output-on-failure
```

## Expected Results
The run stops at the driver's stopping rule in 114 iterations. Speed in solid cells is 0, the plane flux is constant to
four digits, u is -0.55 m/s behind the block near the ground, and 20 m above the roof it is 6.36 m/s against 5.91 m/s
upstream.

## CTest
`building_cube_clean`, `building_cube_run`, `building_cube_flow`, `abort_buildings_newton`, `abort_buildings_lusgs`,
`abort_buildings_levels`, `abort_buildings_mlmg`.
