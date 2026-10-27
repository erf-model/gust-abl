# Forested_Escarpment

## Purpose
A dense forest canopy on steep terrain with an inflow face: the benchmark benchmark's forested escarpment in 2-D, the
set-up the benchmark runs need before they start (the ERF set-up of that benchmark is in ERF's
`Exec/CanonicalTests/Real_Terrain/benchmark`). A tanh escarpment of 200 m (steepest slope 22° at x = 3000 m) carries the
benchmark canopy (uniform profile, 25 m, LAI 5, C_d 0.15) on the slope and the plateau's edge, x in [2200, 3600] m.
benchmark's layers (8 m first layer stretched by 1.12, about three layers in the canopy), neutral log-law inflow on the
low x face, outflow on the high x face, periodic in y, LU-SGS with grid sequencing over three levels.

The check shows that the case converges at `steady.cfl_max = 2` with the same answer on 1, 2 and 4 ranks, that the
canopy slows the wind, that a restart continues bit for bit, and that the inflow-face warning fires when the forest
reaches the inflow face. At the default CFL cap (8) the case cycles (docs/canopy.rst, "Steep terrain and the
boundaries", has the measured caps).

## Files
| File | Role |
|---|---|
| `make_case.py` | Writes `terrain.xyz`, `canopy_height.xyz` and `canopy_lai.xyz` on 25 m points into a folder; `--forest-at-inflow` extends the forest to the inflow face. |
| `inputs_forested_escarpment` | The deck (240 x 8 x 30 cells, 6000 m x 200 m), reading the tables from the run folder. `init.first_guess` is given on the command line. |
| `check_forested_escarpment.py` | Runs the six runs below, each in its own folder with its tables, and checks them; the exit code is the verdict. |

The runs: `np2/`, `np1/`, `np4/` to the verdict on 2, 1 and 4 ranks; `first/` for 40 iterations, then continued from
its checkpoint to 80; `straight/` for 80 iterations; `inflow/` with the forest at the inflow face, 1 iteration.

## Data and acknowledgements
No external data: the terrain and the forest are made up. The canopy's height, LAI and drag coefficient follow the
benchmark complex terrain benchmark, Phase I (the benchmark authors, Zenodo,
(reference removed), CC BY 4.0; https://benchmark-benchmark.readthedocs.io): its prescribed leafy
LAI of 5.0 and the ERF set-up's assumed 25 m beech canopy and C_d of 0.15. We thank the the authors team and the IEA Wind Task 57
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

## CTest
`rt_forested_escarpment` (label `regression`, double precision only; about 1 minute on two ranks).
