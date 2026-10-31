# Column_Marching

## Purpose
The time-marching column (`init.column_time` > 0, `init.column_dt`) against ERF's near-1-D columns with a
stratified inflow. The ERF columns are local runs, not in the repository:
2 x 2 x 240 cells of 10 m, ERF's one-equation k closure, MOST with z0 = 0.1 m and a fixed kinematic surface heat
flux, Coriolis at 48.7 N with a geostrophic wind, 300 K to an 8 K inversion at 1000-1100 m, then 3 K/km (stable) or
0 (unstable) above, 9 h in 1 s steps. GUST marches the same column (240 levels to 2400 m, the same forcing, surface
flux and sounding, starting from ERF's state at t = 0) with backward Euler steps to 9 h, and the two are compared at
9 h.

## Files
| File | Role |
|---|---|
| `inputs_column_marching` | The GUST deck (the high_stable_G2 case, q = -0.004 K m/s). |
| `compare_erf_column.py` | Per case: reads ERF's run folder (inputs, input_sounding, plotfiles), writes ERF's state at t = 0 as `bc.inflow_file`, runs a copy of the deck with the case's values, and prints R_B between 46 and 73 m, the speed at 300 m and theta at 10 and 100 m for both. |

## Data and acknowledgements
No downloaded dataset is used. The ERF columns are local runs of the ERF code
(https://github.com/erf-model/ERF); they are not in the repository.

## How to run
```bash
python3 Tests/Canonical_RANS/Column_Marching/compare_erf_column.py --gust <build>/gust \
    --erf-root <folder of ERF's columns> --out <run folder> [--dt 1]
```
The default cases are `high_stable_G2/q-0.004`, `low_stable_G3/q-0.0025`, `high_unstable_lapse0/q0.018` and
`low_unstable_lapse0/q0.03` below `--erf-root`.

## Expected Results
Measured at 9 h with `init.column_dt = 1` s (2 s of wall time per column on a laptop):

| Case | q (K m/s) | R_B ERF | R_B GUST | U(300 m) ERF | U(300 m) GUST |
|---|---|---|---|---|---|
| high_stable_G2 | -0.004 | +0.122 | +0.111 | 7.76 | 7.80 |
| low_stable_G3 | -0.0025 | +0.097 | +0.089 | 6.88 | 6.88 |
| high_unstable_lapse0 | 0.018 | -0.266 | -0.242 | 7.21 | 7.17 |
| low_unstable_lapse0 | 0.03 | -0.360 | -0.328 | 7.55 | 7.49 |

GUST's R_B is 8-9 % smaller in magnitude than ERF's in all four cases; the speed at 300 m is within 1 % and theta at
10 and 100 m within 0.07 K. Steps of 10 s and 60 s give the same numbers to the digits shown (R_B within 0.0004),
so the difference is the two closures and surface treatments, not the time step. Between 46 and 73 m the stable
columns have the same shear as ERF's and a 7-9 % smaller theta difference (GUST's k is 5-8 % higher at 50-200 m);
the unstable columns have the same theta difference and about 5 % more shear. The common factor is the first cell's
k: by default GUST holds the wall value (u*^3 + kappa B z_1)^(2/3)/Cmu0^2 there (B = max(g q/theta0, 0)), about 1.9
times ERF's evolved first-cell k (0.23 vs 0.12 m^2/s^2 in high_stable_G2). The heat budget closes: with
zhi_theta_grad = 0 the column's heat content grows by q t to 4e-13.

With `--extra init.column_wall_k=solved` the first cell's k is solved (no k flux through the ground) instead:

| Case | R_B ERF | R_B GUST, fixed | R_B GUST, solved | k_1 ERF | k_1 fixed | k_1 solved (m^2/s^2) |
|---|---|---|---|---|---|---|
| high_stable_G2 | +0.1220 | +0.1113 | +0.1135 | 0.122 | 0.231 | 0.160 |
| low_stable_G3 | +0.0967 | +0.0886 | +0.0904 | 0.106 | 0.201 | 0.138 |
| high_unstable_lapse0 | -0.2660 | -0.2424 | -0.2453 | 0.207 | 0.393 | 0.254 |
| low_unstable_lapse0 | -0.3600 | -0.3284 | -0.3320 | 0.239 | 0.457 | 0.292 |

The first-cell k comes to 1.2-1.3 times ERF's, but R_B moves only from 8-9 % to 6.5-8 % below ERF's; the speed at 300
m and theta change by at most 0.01 m/s and 0.01 K, and steps of 10 and 60 s give the same numbers (R_B within
0.0004). So the first-cell k is a small part of the R_B difference; the rest is in the closure above it.

## CTest
None: the ERF columns are not in the repository. The marching mode is pinned by the unit tests
`ColumnMarching.*` (the surface-flux budget, the neutral steady limit, the flux surface's closed form,
the solved first-cell k balance of `init.column_wall_k = solved`, the inputs).
