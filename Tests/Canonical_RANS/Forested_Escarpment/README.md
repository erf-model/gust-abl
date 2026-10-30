# Forested_Escarpment

## Purpose
A dense forest canopy on steep terrain with an inflow face: the benchmark benchmark's forested escarpment in 2-D, the
set-up the benchmark runs need before they start (the ERF set-up of that benchmark is in ERF's
`Exec/CanonicalTests/Real_Terrain/benchmark`). A tanh escarpment of 200 m (steepest slope 22° at x = 3000 m) carries the
benchmark canopy (uniform profile, 25 m, LAI 5, C_d 0.15) on the slope and the plateau's edge, x in [2200, 3600] m.
benchmark's layers (8 m first layer stretched by 1.12, about three layers in the canopy), neutral log-law inflow on the
low x face, outflow on the high x face, periodic in y, LU-SGS with grid sequencing over three levels.

The check shows that the case converges at the default CFL cap (8) with the default step floor
(`steady.dtau_floor = 0.1`) with the same answer on 1, 2 and 4 ranks, that the canopy slows the wind, that a restart
continues bit for bit, and that the inflow-face warning fires when the forest reaches the inflow face. Without the
floor (`steady.dtau_floor = 0`, the step before PR-62b) the case cycles at caps of 4 and above and needs
`steady.cfl_max` 2 to 2.5 (docs/canopy.rst, "Steep terrain and the boundaries", has the cause and the measurements).

The same escarpment with a stable sounding tests `simplec.theta = frozen` (θ prescribed from the inflow profile by
height above the ground, buoyancy in the closure only; docs/solution_algorithm.rst, "Prescribed potential
temperature"): `inputs_frozen_theta` and `check_frozen_theta.py`.

Both decks also test the leaves' heat exchange (`canopy.heat = fixed`, docs/canopy.rst, "Canopy heat"):
`check_canopy_heat.py` runs θ solved with leaves cooler than the air, and θ frozen with `canopy.heat_frozen = skip`
and `blend`.

## Files
| File | Role |
|---|---|
| `make_case.py` | Writes `terrain.xyz`, `canopy_height.xyz` and `canopy_lai.xyz` on 25 m points into a folder; `--forest-at-inflow` extends the forest to the inflow face. |
| `inputs_forested_escarpment` | The deck (240 x 8 x 30 cells, 6000 m x 200 m, the default CFL cap and step floor), reading the tables from the run folder. `init.first_guess` is given on the command line. |
| `check_forested_escarpment.py` | Runs the six runs below, each in its own folder with its tables, and checks them; the exit code is the verdict. |
| `inputs_frozen_theta` | The deck with a stable sounding (300 K to 1000 m, an 8 K inversion to 1100 m, 3 K/km above); `simplec.theta = frozen` on the command line. |
| `inputs_frozen_theta_file` | The same with the inflow from a table (`frozen_inflow.txt`, written by the check) whose θ differs from the sounding. |
| `check_canopy_heat.py` | The leaf-heat runs `cool/` (θ solved, leaves at 298 K, 2 ranks), `cool_np4/`, `frozen/` (frozen θ, no heat), `skip/` and `blend/` (frozen θ, leaves at 304 K), each to the verdict, and their checks (its docstring). |
| `check_frozen_theta.py` | The frozen-θ runs `frozen/` (2 ranks, to the verdict), `np4/`, `file/` (the table deck, `init.z_agl_min = 5`, 3 iterations) and `solve/` (θ solved, 5 iterations), and their checks (its docstring). |

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
| Iterations to the verdict on 2 ranks | 256 | at most 300 |
| Boundary-zone warnings (forest away from the inflow) | none | none |
| 1 and 4 ranks against 2, max difference (velocities over max \|u\|, p and k over their own max) | 1.6e-11 and 4.4e-11 (256 iterations each) | 1e-9 |
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

Canopy heat (`check_canopy_heat.py --max-iterations 400`):

| Check | Measured | Tolerance |
|---|---|---|
| Iterations to the verdict: cool, cool on 4 ranks, frozen, skip, blend | 201, 201, 209, 209, 212 | at most 400 |
| cool: θ range (leaves 298 K, air 300 K) | 298.00028 to 300.00000003 K | within [298, 300] K ± 1e-6 K |
| cool: θ upstream of the forest (x < 2000 m) against 300 K | 1.5e-8 K | 1e-3 K |
| cool: 4 ranks against 2 (u, v, w, p, k; θ) | 4.2e-11; 5.9e-9 K | 1e-9; 1e-7 K |
| skip: warning, and every field against the frozen run | logged, 0 | logged, bit for bit |
| blend: θ outside the canopy against the frozen run | 0 (56072 cells) | bit for bit |
| blend: θ in the canopy between θ_f and 304 K; largest weight w | yes (1528 cells); 0.69 | yes; above 0.05 |

With the leaves warmer than the air and θ solved (buoyancy in the momentum equations) the run converges up to
301.5 K over the neutral 300 K air (204 iterations) and does not settle at 302 K: r cycles between 0.3 and 0.4, also
at `steady.cfl_max = 2` for 1500 iterations (docs/canopy.rst has the sweep). The heated canopy air is
unstably stratified; the check does not run that case. θ frozen with `blend` is the route for warm leaves.

## CTest
`rt_forested_escarpment` (label `regression`, double precision only; about 1 minute on two ranks) and
`rt_frozen_theta` (label `regression`; about 30 s on two ranks) and `rt_canopy_heat` (label `regression`; about 2 min).
