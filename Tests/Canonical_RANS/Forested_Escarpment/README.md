# Forested_Escarpment

## Purpose
A dense forest canopy on steep terrain with an inflow face: the WINSENT benchmark's forested Albtrauf in 2-D, the
set-up the WINSENT runs need before they start (the ERF set-up of that benchmark is in ERF's
`Exec/CanonicalTests/Real_Terrain/WINSENT`). A tanh escarpment of 200 m (steepest slope 22° at x = 3000 m) carries the
WINSENT canopy (uniform profile, 25 m, LAI 5, C_d 0.15) on the slope and the plateau's edge, x in [2200, 3600] m.
WINSENT's layers (8 m first layer stretched by 1.12, about three layers in the canopy), neutral log-law inflow on the
low x face, outflow on the high x face, periodic in y, LU-SGS with grid sequencing over three levels.

The check shows that the case converges at `steady.cfl_max = 2` with the same answer on 1, 2 and 4 ranks, that the
canopy slows the wind, that a restart continues bit for bit, and that the inflow-face warning fires when the forest
reaches the inflow face. At the default CFL cap (8) the case cycles (docs/canopy.rst, "Steep terrain and the
boundaries", has the measured caps).

The same escarpment with a stable sounding tests `simplec.theta = frozen` (θ prescribed from the inflow profile by
height above the ground, buoyancy in the closure only; docs/solution_algorithm.rst, "Prescribed potential
temperature"): `inputs_frozen_theta` and `check_frozen_theta.py`.

## Files
| File | Role |
|---|---|
| `make_case.py` | Writes `terrain.xyz`, `canopy_height.xyz` and `canopy_lai.xyz` on 25 m points into a folder; `--forest-at-inflow` extends the forest to the inflow face. |
| `inputs_forested_escarpment` | The deck (240 x 8 x 30 cells, 6000 m x 200 m), reading the tables from the run folder. `init.first_guess` is given on the command line. |
| `check_forested_escarpment.py` | Runs the six runs below, each in its own folder with its tables, and checks them; the exit code is the verdict. |
| `inputs_frozen_theta` | The deck with a stable sounding (300 K to 1000 m, an 8 K inversion to 1100 m, 3 K/km above); `simplec.theta = frozen` on the command line. |
| `inputs_frozen_theta_file` | The same with the inflow from a table (`frozen_inflow.txt`, written by the check) whose θ differs from the sounding. |
| `check_frozen_theta.py` | The frozen-θ runs `frozen/` (2 ranks, to the verdict), `np4/`, `file/` (the table deck, `init.z_agl_min = 5`, 3 iterations) and `solve/` (θ solved, 5 iterations), and their checks (its docstring). |

The runs: `np2/`, `np1/`, `np4/` to the verdict on 2, 1 and 4 ranks; `first/` for 40 iterations, then continued from
its checkpoint to 80; `straight/` for 80 iterations; `inflow/` with the forest at the inflow face, 1 iteration.

## Data and acknowledgements
No external data: the terrain and the forest are made up. The canopy's height, LAI and drag coefficient follow the
WINSENT complex terrain benchmark, Phase I (A. Rettenmeier, S. Lott, H. Knaus, J. Hohaus, Zenodo,
https://doi.org/10.5281/zenodo.23097333, CC BY 4.0; https://winsent-benchmark.readthedocs.io): its prescribed leafy
LAI of 5.0 and the ERF set-up's assumed 25 m beech canopy and C_d of 0.15. We thank the ZSW team and the IEA Wind Task 57
JAM benchmark for the benchmark description.

## How to run
```bash
ctest --test-dir <build> -R rt_forested_escarpment --output-on-failure
```
or by hand, from an empty folder:
```bash
python3 <repo>/Tests/Canonical_RANS/Forested_Escarpment/check_forested_escarpment.py --max-iterations 300 -- \
    mpiexec -n 2 <build>/gust <repo>/Tests/Canonical_RANS/Forested_Escarpment/inputs_forested_escarpment
```

## Expected Results
Measured on macOS arm64 (AppleClang, double precision):

| Check | Measured | Tolerance |
|---|---|---|
| Iterations to the verdict on 2 ranks | 185 | at most 300 |
| Boundary-zone warnings (forest away from the inflow) | none | none |
| 1 and 4 ranks against 2, max difference (velocities over max \|u\|, p and k over their own max) | about 4e-11 (185 iterations each) | 1e-9 |
| u 12 m above the ground, forest (column at x = 2912.5 m) / upstream (x = 1512.5 m) | 0.72 / 3.52 m/s (0.20) | below 0.5 |
| Restart: continued vs straight, every field | 0 | bit for bit |
| Forest at the inflow face: canopy cells in the warning | 24 | more than 0 |

v is rounding noise in this 2-D case (its largest value is about 5e-8), so the velocities are compared against the
largest |u|.

Frozen θ (`check_frozen_theta.py --max-iterations 400`):

| Check | Measured | Tolerance |
|---|---|---|
| Iterations to the verdict on 2 ranks, θ residual | 209, 0 in every row | at most 400, 0 |
| Level 0 θ against the sounding at each cell's height above the ground | 5.7e-14 K | 1e-9 K |
| Logged min, max, mean of θ on levels 0, 1, 2 against the sounding over each level's cells | equal to the printed digits | 1e-9 K |
| Turbulence plotfile `R_theta` | 0 | 0 |
| 4 ranks against 2 | 1.7e-10 | 1e-9 |
| Table profile with `init.z_agl_min = 5`: θ against the table at the unfloored height | 5.7e-14 K | 1e-9 K |
| Table profile: first-cell Rt (the ground ghosts extend the first two cells) | min 0.023 | above 0 |
| θ solved: θ residual | above 0 in every row | above 0 |

## CTest
`rt_forested_escarpment` (label `regression`, double precision only; about 1 minute on two ranks) and
`rt_frozen_theta` (label `regression`; about 30 s on two ranks).
