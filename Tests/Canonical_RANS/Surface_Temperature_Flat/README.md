# Surface_Temperature_Flat

## Purpose
The ground condition `bc.surface_stability = surface_temperature` in a short 3-D run. A flat periodic ABL
(16 x 16 x 32 cells over 2560 m x 2560 m x 1000 m, four boxes, two ranks) with the log law (u* = 0.35 m/s) from the
south and no Coriolis force, z0 = 0.1 m and theta = 300 K up to 700 m. The ground temperature comes from a map that
rises from 297 K in the west columns to 303 K in the east ones, so the run holds stable and unstable columns side by
side. The wind blows along the columns of equal ground temperature. With a wind from the west the warm air of the
east columns would cross the periodic seam onto the coldest columns and warm them. The run stops after 20
iterations; it is not meant to reach a steady state.

The check solves Monin-Obukhov similarity again in each column, in Python and with none of GUST's code, from the
plotfile's first-cell wind and theta, and compares the result with the run's u*, theta_s and heat flux. It also
checks that the flux has reached theta: every column with |theta_s - 300| > 1 K must have its first-cell theta moved
from 300 K towards theta_s.

## Files
| File | Role |
|---|---|
| `inputs_surface_temperature_flat` | The deck. The map is given on the command line (`bc.surface_temperature_file`). |
| `theta_s_ramp.xyz` | The ground-temperature map: one `x y theta_s` point at each cell centre, theta_s = 297 + 6 i / 15 K for column i. |
| `check_surface_temperature.py` | Reads the plotfile and the turbulence plotfile; bisection on zeta for every column; compares u*, theta_s and q; checks that theta moved towards theta_s. |

## Data and acknowledgements
No measured data, and nothing has to be downloaded. The check uses ERF's plotfile reader (`extern_erf/erf_plotfile.py`,
recorded in `extern_erf/EXTRACTED.md`; BSD 3-clause with the LBNL clause, `extern_erf/LICENSE.md`). Thanks to the ERF
developers.

## How to run
```bash
mpirun -n 2 ../../../build-rel/gust inputs_surface_temperature_flat bc.surface_temperature_file=theta_s_ramp.xyz
python3 check_surface_temperature.py .
```

## Expected Results
- The log prints the map summary `min 297, max 303, mean 300 K`.
- In each of the 256 columns, the run's u*, theta_s and q match the Python solve to 1e-9 (measured: 1.3e-14,
  0 and 2.6e-14).
- The columns split into 128 stable and 128 unstable ones (|d/L| > 0.01); the check requires both.
- theta has moved towards theta_s in all 160 columns with |theta_s - 300| > 1 K.
- Flipping the sign of theta_r - theta_s in the wall model makes the check fail (u* off by 176 %, q by 1000 %, theta
  moved the wrong way in all 160 columns). Setting the ground flux of theta to zero makes it fail too (theta moved in
  80 of 160 columns).

## CTest
`rt_surface_temperature_flat_clean`, `rt_surface_temperature_flat` (2 ranks) and `rt_surface_temperature_flat_check`
(labels `regression;smoke`). `abort_surface_temperature_inputs` (a value outside 150-350 K) and
`abort_surface_temperature_both` (a value and a map together) check the input aborts.
