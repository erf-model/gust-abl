# Building_Courtyard

## Purpose
A courtyard open to the wind: a wall 30 m wide with two wings and a low annex behind it, on flat ground in a log-law
wind. Its inner corners and the annex's roof against the wall are cells with two walls, where one face of a velocity
component lies on the other wall. The case pins the wall law's speed in those cells (`docs/buildings.rst`): from the
faces that are not fixed rows (`buildings.wall_speed = open_faces`, the default), and from every face (`all_faces`,
the form before). It is not a validation; PR-44 compares with measurements.

## Files
| File | Role |
|---|---|
| `inputs_courtyard` | 48 x 24 x 24 cells of 5 m, inflow at xlo, outflow at xhi, periodic in y, no gravity. |
| `courtyard.txt` | The buildings in ERF's format: the wall (x 75..85 m, y 45..75 m), the wings (x 65..75 m; y 45..55 m and 65..75 m), all 20 m high, and the annex (x 85..95 m, y 55..65 m, 10 m high). |
| `check_courtyard.py` | Speeds and k in the solid cells, the volume flux along x, and k in the inner corner cells (`--k-corner`, `--k-top`). The exit code is the verdict. |

## Data and acknowledgements
No external data: the map is made up for the test. The wall law follows ERF's `erf.if_wall_form = fraction_stress`
and Kynema-SGF's wall_form 7.

## How to run
```bash
ctest --test-dir <build> -R building_courtyard --output-on-failure
```

## Expected Results
The run stops at the driver's stopping rule. Speed in solid cells is 0, k there is at its floor, and the plane flux is
constant to within 2 %. These values are not validated.

| `buildings.wall_speed` | Iterations | k in an inner corner, second layer | k in an inner corner, fourth layer |
|---|---|---|---|
| `open_faces` (default) | 135 | 0.0227 m²/s² | 0.633 m²/s² |
| `all_faces` | 141 | 0.0081 m²/s² | 0.432 m²/s² |

Measured on 2 and 4 ranks; the tests run on 2 and pin the k values within 5 %. The iteration counts are not pinned. The
velocities of the two forms differ by about 0.02 m/s at most (w by 0.0204 m/s in the inner corner), and k by up to
0.20 m²/s². The `all_faces` plotfile equalled the one of `main` before the
open faces' mean (c57ad06) at zero tolerance. Since the solid cells' rows left the residual norm of k the runs stop at
other iterations than `main` at 8448d15 (138 and 139); the default's velocities differ from that plotfile by at most
1.4e-4 m/s, and its building distances not at all. Theta's building terms (`buildings.theta_mode`, `ambient` by
default) leave the velocities and k bit for bit, in the same 135 iterations (measured against
`buildings.theta_mode = none`); theta differs by up to 0.007 K.

## CTest
`building_courtyard_clean`, `building_courtyard_run`, `building_courtyard_flow`, `building_courtyard_all_faces_clean`,
`building_courtyard_all_faces_run`, `building_courtyard_all_faces_flow`.
