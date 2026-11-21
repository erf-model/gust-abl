# Library_Ensemble

## Purpose
Ensemble statistics and WindNinja-format output of the wind library (`tools/make_gust_library.py ensemble`,
`tools/gust_library/ensemble.py` and `grids.py`) on a sweep GUST ran. A library of four directions gives the
statistics over a distribution of directions and speeds, with no new runs. The case checks them against the
same statistics computed from the runs' own plotfiles, and checks that the grids sit on the map where the
terrain came from. The options and conventions are described in `docs/case_setup.rst`.

## Files
| File | Role |
|---|---|
| `library_ensemble.yaml` | A Gaussian hill (60 m high, 200 m wide) at the centre of a 2 km domain placed by `domain.center_lat_lon` (40 N, 105.5 W: UTM zone 13N), 50 m cells in 9 boxes of uneven size, neutral log-law inflow of 8 m/s at 10 m over z0 = 0.03 m, four directions 0, 90, 180 and 270 degrees. |
| `check_library_ensemble.py` | Writes `hill.asc` in UTM metres, writes the case with `make_gust_case.py --strict`, runs `run.sh`, builds the library at 10 m and 50 m, runs `make_gust_library.py ensemble` and checks its output. |

## Data and acknowledgements
No downloaded data. The hill grid is written by `check_library_ensemble.py` from a formula.

## How to run
```bash
ctest --test-dir <build> -R "^rt_library_ensemble$"
```
By hand, from an empty directory:
```bash
python3 <repo>/Tests/Infrastructure/Library_Ensemble/check_library_ensemble.py --gust <build>/gust --mpi "mpiexec -n 4"
```
The four runs take 155 to 159 iterations each on 4 ranks.

## Expected Results
The members are the directions 0, 45, 90, 180 and 270 degrees with weights 0.1, 0.15, 0.2, 0.25 and 0.3, each at
6 and 10 m/s with weights 0.5 and 0.5 (10 members). 45 degrees lies between two runs: the library blends it, and
the check blends the runs from 0 and 90 degrees by hand (half of each, divided by the length of the half-and-half
mix of their unit winds at the reference point). The statistics are the mean speed, its standard deviation, the 10th, 50th and 90th
percentiles, the probability of a speed above 8 m/s and the probability of a reversal.
- `make_gust_case.py --strict` and `run.sh` exit with 0, and each run reaches GUST's steady verdict.
- In every column at 10 m and 50 m, the tool's statistics equal the same statistics computed by plain loops from
  each run's wind scaled to the member's speed (tolerance 1e-12 of the largest mean speed; probabilities 1e-12;
  the percentiles' cumulative weight within max(1e-12, 4 n eps) for n members, as in the tool).
  Measured: at most 5.3e-15 m/s; probabilities 0.
- The `_vel.asc` and `_ang.asc` grids: 40 columns and rows, cell size 50 m, `NODATA_value` -9999; the corner is
  the domain's south-west corner in UTM (the easting and northing of `domain.center_lat_lon` less 1000 m,
  tolerance 1e-6 m; measured 1.2e-7 m); north row first; the values equal the mean speed (1e-7 relative,
  measured 8.5e-10) and the direction of the vector mean (1e-6 degrees, measured 5e-7). The differences are the
  rounding of the file's nine significant digits.
- The `.prj` names WGS 84 / UTM zone 13N with central meridian -105.
- The CSV at the hilltop: the mean speed equals the statistics of the members' winds there (1e-8 relative), and
  the easting is the hilltop's in UTM (1e-6 m).
- The library's ground written as a grid with the same corner: its centroid (weights: the ground above its lowest
  value) sits at the hill's top in UTM within 1 m. Measured: 1.2e-7 m in easting, 0.68 m in northing. The
  northing offset is in GUST's own mesh: `terrain.xyz` is symmetric, and the inverse-distance interpolation of
  `xyz` terrain (6 neighbours) breaks ties between points at equal distance by their order in the file, south
  rows first. The ground is symmetric to 7e-15 m east-west and to 0.98 m north-south (60 m hill).
- Printed, not checked: mean speed 5.70 to 8.01 m/s at 10 m and 7.70 to 9.24 m/s at 50 m, standard deviation up to
  2.00 and 2.31 m/s. The probability of a reversal is 0 in every column: the flow stays attached over this hill.
  A reversal is checked by the unit tests (`tools/test_gust_ensemble.py`).
- Mutants: weights not divided by their sum (10 checks fail), the corner half a cell east and north (10 checks
  fail), the direction written as the direction the wind blows towards (2 checks fail, 180 degrees) and the
  library's blend not divided by its length (18 checks fail).

Measured in double precision; single precision is not measured.

## CTest
`rt_library_ensemble` (label `regression`, 4 processors; not on Windows, where `run.sh` cannot run). It and
`gust_ensemble_selftest` (the unit tests, `tools/test_gust_ensemble.py`) are registered only when the configure's
Python has PyYAML and numpy.
