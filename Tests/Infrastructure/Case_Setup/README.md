# Case_Setup

## Purpose
The Python case set-up tool, `tools/make_gust_case.py`, end to end: a YAML case description becomes a case folder
(decks, `terrain.xyz`, `run.sh`), the decks pass the deck checker (`tools/check_gust_inputs.py`), `run.sh` runs
GUST, and the result meets the closed forms the YAML file asks for. Each YAML file here is also a small example of
the case description (`docs/case_setup.rst`).

## Files
| File | Role |
|---|---|
| `flat_sweep.yaml` | Flat ground, moving lid (`wind.drive: couette`), two directions (240 and 0 degrees): one run folder per direction. |
| `flat_nest.yaml` | Flat ground from 240 degrees with a nested patch at half the cell size; the patch starts from the parent's checkpoint. |
| `column_target.yaml` | `wind.drive: geostrophic_column` at 45 degrees north: GUST's 1-D column finds the geostrophic wind that gives 8 m/s at 10 m from 240 degrees; `run.sh` passes it to the 3-D run. |
| `column_axis.yaml` | `column_target`'s set-up with the 10 m wind along y (from 0 degrees) and along x (from 270 degrees): the geostrophic wind turns with height, so the faces parallel to the 10 m wind are inflow and outflow faces, not slip walls. A third direction, 240 degrees, is off the axes for `library.json`. |
| `hill.yaml` | A Gaussian hill from an ESRI ASCII grid in a projected frame, off the domain centre, on a plane rising to the east. |
| `upstream.yaml` | `wind.upstream_roughness: auto` from 270 degrees over a half-plane regional map (0.03 m west of x = 1200 m, 1 m east), which also sets the ground; `wind.drive: geostrophic_inflow` at 45 degrees north with the recommended lee setting. |
| `check_case_setup.py` | Copies a YAML file to the run directory (writes `hill.asc` for the hill and `regional.xyz` for `upstream`), writes the case, runs `run.sh` and checks the result. |

## Data and acknowledgements
No downloaded data. The hill grid and the regional roughness map are written by `check_case_setup.py` from formulas.

## How to run
```bash
ctest --test-dir <build> -R "^rt_case_setup_"
```
By hand, from an empty directory:
```bash
python3 <repo>/Tests/Infrastructure/Case_Setup/check_case_setup.py flat_sweep --gust <build>/gust --mpi "mpiexec -n 2"
```

## Expected Results
Every YAML file here has a wind of 8 m/s at 10 m over z0 = 0.03 m. The 10 m wind is the planar mean, linear in
ln(z + z0) between the cell centres around 10 m.
- All cases: `make_gust_case.py --strict` exits with 0 (written, the deck checker finds no error or warning) and
  `run.sh` exits with 0.
- `flat_sweep`, `flat_nest` (parent and patch): 8 m/s at 10 m to 3e-4, along the asked direction to 0.05 degrees,
  and GUST's steady verdict. The patch's cells are half the parent's in x and y.
  Measured on 2 ranks: 7.9994 m/s, 0.003 degrees (240); 8.0002 m/s (0); patch 8.0009 m/s, 0.004 degrees.
- `column_target`: the column's table (`column/column_guess.dat`) holds 8 m/s at 10 m to 2e-4 (the column stops at
  |R - 1| < 1e-4) along 240 degrees to 0.05 degrees. The 3-D run, which takes the column's `forcing.ug` and
  `forcing.vg` from `column/log.txt`, has 8 m/s to 1 % and the direction to 0.3 degrees.
  Measured: table 7.9993 m/s, 0.001 degrees; 3-D 7.981 m/s, 0.09 degrees. The 3-D run is a separate solve on its
  own mesh; the column is the 1-D limit of it, so the two differ by the discretisation. In this small domain the
  10 m wind follows the inflow table, so a separate check covers the hand-over: the 3-D run solves its own column
  first guess from the `forcing.ug` and `forcing.vg` that `run.sh` passed, and its wind at the top (printed in
  `log.txt`) matches the table's top row to 1e-4. Measured: equal to the 4 printed decimals.
  `library.json` (written by `run.sh`'s last step, `make_gust_case.py --collect`) holds 8 m/s at 10 m to 2e-4
  along 240 degrees to 0.01 degrees (the stopping rule allows 1e-4 rad, 0.0057 degrees). Its G equals the last
  `forcing.ug` and `forcing.vg` of `column/log.txt` (what `run.sh` passes to the 3-D run), turned clockwise of the
  10 m wind by 10 to 25 degrees: between half the neutral drag law's angle for the column's u* (23.5 degrees) and
  all of it. Measured: 7.99935 m/s, 0.0014 degrees; |G| 15.9089 m/s, turning -18.1405 degrees.
- `column_axis`: no lateral face of either deck is a slip wall, and every cell's wind in the 3-D run is the column
  table's at the cell's height (linear between the table's levels) to 4 % of the table's top wind. Measured on 2
  ranks: 0.61 % (from 0) and 1.97 % (from 270). The largest difference from 270 is in the first cell at the
  x-high, y-low corner: below 10 m the wind is turned the other way (backed), so it enters through the y-low
  outflow face there; the 99th percentile is 1.15 %. The tool before this case made the faces parallel to the
  10 m wind slip walls, which block the flow that crosses them aloft: 72 % (from 270, 197 iterations against 119
  to 129 with the turned faces).
  `library.json`: the same checks as `column_target` from 0, 270 and 240 degrees, and the same |G| and turning
  angle for the three, as the neutral column is rotation-invariant: 1e-6 relative and 1e-4 degrees in double
  precision. A quarter turn is exact in floating point, so 240 degrees is the direction that tests it. With
  `--precision SINGLE` (CTest passes `GUST_PRECISION`) the tolerances are the stopping rule's bound for two
  directions, 2e-4 and 0.0115 degrees. Measured (double, 2 ranks): |G| apart by 2.0e-10 (270) and 1.1e-9 (240),
  turning by 3.4e-8 and 3.7e-8 degrees.
  Mutants: reading the deck's first guess for G gives a turning of -0.0014 degrees (fails the 10 to 25 degree band
  and the log's forcing lines); interpolating linearly in z gives 7.962 m/s (4.6e-3); no `--collect` step in
  `run.sh` fails "library.json written by run.sh".
- `hill`: the ground in GUST's plotfile, from the first cell centre (`z = zeta + h (1 - zeta / top)`), against the
  formula at the cell centres: the crest (after the plane is removed) in the same cell, the relief to 2 % and the
  RMS difference (each field less its mean) below 0.3 % of the relief. Measured: same cell, 0.7 %, 0.12 %. The
  difference is the interpolation from the 20 m grid to the 50 m mesh.
- `upstream`: the 3-D deck's `bc.inflow_z0` and the column deck's `bc.surface_z0` are 0.03 m (the map west of the
  domain; the map has 1 m east of it), GUST's steady verdict, and the 10 m wind (mean over y) more than 4 cells from
  the x faces varies by less than 1 % along x. Measured on 2 ranks: 0.62 % (126 iterations); the last cells at the
  outflow face drop a further 0.8 %, an effect of the outflow face that a uniform roughness gives too (the decks
  are those of `wind.upstream_roughness: 0.03`, comments aside). Mutant: the meteorological direction read as the
  direction the wind blows towards gives 0.715 m in both decks and 11.8 %.
- Each tolerance fails on a one-line mutant of the tool: u* from ln(z / z0) in place of ln((z + z0) / z0) gives
  8.0035 m/s (4.4e-4); a DEM read half a cell off gives an RMS of 0.66 %.

## CTest
`rt_case_setup_flat_sweep`, `rt_case_setup_flat_nest`, `rt_case_setup_column_target`, `rt_case_setup_column_axis`,
`rt_case_setup_hill`, `rt_case_setup_upstream`
(labels `regression;smoke`; not on Windows, where `run.sh` cannot run). They and `gust_case_selftest` (the tool's
unit tests, `tools/test_gust_case.py`) are registered only when the configure's Python has PyYAML and numpy.
`gust_input_checker_selftest` runs the deck checker's unit tests (standard library only).
