# Building_Hill

## Purpose
One block (30 x 30 m, 20 m high above the local ground) on the windward flank of a 30 m ridge, in a log-law wind: the
fixed rows, the masked pressure correction, the form 7 wall law on the roof and the sides and the closure's building
terms (`docs/buildings.rst`) over terrain. It pins that the
iteration converges with the building's floor and roof following the ground. It is not a validation.

## Files
| File | Role |
|---|---|
| `inputs_block_hill` | 96 x 32 x 32 cells of 5 m (flat layers), inflow at xlo, outflow at xhi, periodic in y, no gravity. |
| `block.txt` | The block in ERF's format, heights above the local ground. |
| `check_block_hill.py` | Speeds in the solid, the volume flux along x, reversal behind the block, the roof's wall cell. The exit code is the verdict. |
| `plot_block_hill.py` | Writes `docs/figures/buildings/block_on_hill.png`. Not a test. |

## Data and acknowledgements
No external data: the map and the ridge are made up for the test. The wall law follows ERF's
`erf.if_wall_form = fraction_stress` and Kynema-SGF's wall_form 7.

## How to run
```bash
ctest --test-dir <build> -R building_hill --output-on-failure
```

## Expected Results
The run stops at the driver's stopping rule in 152 iterations, on 1, 2 and 4 ranks. Horizontal speed in solid cells is
0, the plane flux stays within 0.4 % of the first plane's, and the lowest u behind the block is -2.60 m/s. The air
reverses in a thin layer on the roof, one cell deep: u in the roof's wall cell is -0.15 m/s, the same to four digits on
1, 2 and 4 ranks. The test pins u two cells above it, 7.42 m/s, and k in the wall cell, 1.01 m²/s², and requires k at
its floor in the solid cells; the others are measured. These values are not validated.

Earlier forms, measured once (no test runs them): with `buildings.wall_tke = off`, 174 iterations, -2.67 m/s behind,
3.90 m/s in the roof's wall cell, 6.14 m/s two cells above it and k = 13.1 m²/s² in it, the plotfile of `main` at
8f1b605 at zero tolerance; with `buildings.wall_form = roof_stress`, 212 iterations, -1.91 m/s and 3.99 m/s, the
plotfile of `main` at bba99be.

## CTest
`building_hill_clean`, `building_hill_run`, `building_hill_flow`.
