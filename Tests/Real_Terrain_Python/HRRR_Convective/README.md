# HRRR_Convective

## Purpose
Options for convective and weakly unstable HRRR hours (classes C1 and C2 of `tools/hrrr_regime.py`,
`docs/hrrr_regimes.rst`) as the outer level of a GUST nest at the AWAKEN site north of Enid, Oklahoma. The case uses
the chain of `../HRRR_Outer` (`tools/hrrr_outer.py`, `docs/hrrr_outer.rst`): `hrrr_raw` (56 km, 1000 m cells),
`parent` (52 km, 500 m cells, mass-consistent first guess) and the nest (25 km, 125 m cells, four nest faces). It
starts from the scheme that `hrrr_outer.py case` writes (forecast hour 1, top 3 km, the adjusted parent, θv,
Coriolis on, θ solved with a neutral ground, a Rayleigh layer on w and θ, a 100 m length cap) and tries the options
below on three C1 and three C2 hours. The score is the nest against HRRR, as in `../HRRR_Outer`. It is not a
comparison with the AWAKEN observations, and none are used.

## Files
| File | Role |
|---|---|
| `case.yaml` | The meshes and settings of `../HRRR_Outer/case.yaml` (nest half-width 12500 m, 125 m cells, 42 layers from 10 m, ratio 1.08, top 3042 m, z0 = 0.09 m, 4 ranks); `options.py` sets the hour's crop and terrain. |
| `prepare.py` | Downloads `awaken.asc` (`tools/fetch_terrain.py`, 36.05–36.68 N, 97.80–97.01 W) and the crops of the six hours below (`tools/hrrr_outer.py fetch`, 45 km half-width); `--neighbours` adds the hours before and after, for the classifier. |
| `options.py` | Writes the chain of each hour and its nest options, and `run_options.sh`, which runs them one at a time and scores each nest. |
| `score_options.py` | The score table of the options (Markdown). |
| `check_convective.py` | The test deck: the hour may10_19 with the recommended option, checked against the numbers of November 2026 (CTest `rt_hrrr_convective`). |

The hours (valid time, UTC; each is forecast hour 1 of the cycle an hour before). q is HRRR's kinematic virtual heat
flux at the ground and z_i its boundary-layer height, both means over the nest.

| Tag | Valid | Class | q (K m/s) | z_i (m) |
|---|---|---|---|---|
| may10_19 | 2023-05-10 19:00 | C1 | 0.280 | 1503 |
| may10_20 | 2023-05-10 20:00 | C1 | 0.262 | 1635 |
| may10_21 | 2023-05-10 21:00 | C1 | 0.228 | 1743 |
| apr02_18 | 2023-04-02 18:00 | C2 B | 0.311 | 797 |
| apr04_18 | 2023-04-04 18:00 | C2 B | 0.294 | 1539 |
| apr19_21 | 2023-04-19 21:00 | C2 | 0.233 | 1952 |

The options (each changes only the nest, so all share the chain's parent):

| Option | θ | Ground |
|---|---|---|
| scheme | solved | neutral |
| frozen | `simplec.theta = frozen` (the parent's θ; buoyancy only in k) | neutral |
| frozen_flux | frozen | HRRR's heat flux (`bc.surface_stability = heat_flux`) |
| frozen_buoyant | `simplec.theta = frozen_buoyant` (the parent's θ; buoyancy also in w) | neutral |
| frozen_buoyant_flux | frozen_buoyant | HRRR's heat flux |
| solved_flux | solved | HRRR's heat flux (at most 2000 iterations) |
| frozen_tsurf | frozen | HRRR's skin temperature as a θv map (`bc.surface_stability = surface_temperature`); GUST's own MOST flux |
| frozen_buoyant_tsurf | frozen_buoyant | HRRR's skin temperature |
| solved_tsurf | solved | HRRR's skin temperature (at most 2000 iterations) |

The three `_tsurf` options read `surface_theta.xyz`, which `options.py` writes into each of their folders: HRRR's
skin temperature as θv, θv_s = T_s (10^5 / p_s)^0.2857 (1 + 0.61 q_2m), at the columns of the nest
and 6 km around it (`docs/hrrr_outer.rst`). The flux is then GUST's, from the first cell's θ and wind; it is not
HRRR's. Crops fetched before November 2026 have no skin temperature: `hrrr_outer.py` refuses them for these
options; delete them and run `prepare.py` again (it fetches only the crops that are missing).

## Data and acknowledgements
The meteorological fields are from the High-Resolution Rapid Refresh (HRRR) model of NOAA's Global Systems
Laboratory, run operationally by the National Centers for Environmental Prediction: the one-hour forecasts (f01) of
the native-level (`wrfnat`) and surface (`wrfsfc`) files, read from the NOAA Open Data Dissemination archive on AWS
(https://registry.opendata.aws/noaa-hrrr-pds/) with Herbie (Blaylock, https://github.com/blaylockbk/Herbie), which
downloads only the byte ranges of the fields used. We thank NOAA and the HRRR developers for making the forecasts
free to use, Amazon Web Services for hosting the archive, and Brian Blaylock for Herbie. NOAA open data carry no
restriction on use or redistribution; the crops (about 0.8 MB an hour) are not redistributed here, and `prepare.py`
downloads them.

The elevation data are the SRTM-format tiles of the Terrain Tiles data set on the Registry of Open Data on AWS
(https://registry.opendata.aws/terrain-tiles/), made by Mapzen and maintained by the Tilezen project; see
`../README.md` for the providers and our thanks. `fetch_terrain.py` writes the attribution the data need to
`awaken.sources.txt`. No elevation data are redistributed here.

## How to run
In this folder (the scripts find `tools/` from it; `.gitignore` keeps the downloads and runs out of git), with
Herbie installed (`pip install herbie-data`) for the download:
```bash
python3 prepare.py
python3 options.py
GUST=<build>/gust MPI="mpiexec -n 4" options/run_options.sh
python3 score_options.py options
```
The options take about four hours on 4 ranks of a laptop (6 hours × 6 nests, and the 12 coarser runs); solved_flux
is most of it in the C1 hours, where it does not reach a steady state.

## Expected Results
2026-11-25, 4 ranks. Vector rms (m/s) of the nest minus HRRR below 500 m above the ground (`score_options.py` prints
the full table, with bias, rms, direction, the mean speed bias from 5 to 50 m (`bias 5-50`), w, k and iterations):

| Option | may10_19 | may10_20 | may10_21 | apr02_18 | apr04_18 | apr19_21 |
|---|---|---|---|---|---|---|
| scheme | 0.94 | 0.58 | 1.12 | 1.39 | 1.61 | 1.25 |
| frozen | 0.78 | 0.48 | 0.82 | 1.50 | 1.58 | 1.15 |
| frozen_flux | 0.88 | 0.58 | 0.92 | 1.57 | 1.64 | 1.19 |
| **frozen_buoyant** | **0.54** | **0.36** | **0.60** | 1.67 | 1.41 | 1.29 |
| frozen_buoyant_flux | 0.63 | 0.38 | 0.63 | 1.72 | 1.46 | 1.33 |
| solved_flux | 1.93 (not steady) | 1.48 (not steady) | 1.00 (not steady) | **1.00** | **1.30** | **0.90** |

- **frozen_buoyant is the recommended option for C1 hours** (0.12 to 0.24 m/s below frozen, 0.22 to 0.52 m/s below
  the scheme, in each hour). The gain is mostly in speed: the speed bias below 500 m goes from −0.57, −0.36, −0.72
  m/s (frozen) to −0.36, +0.13, −0.06 m/s.
- Its w is about twice HRRR's (rms 0.065 to 0.077 against 0.026 to 0.038 m/s), most of it in the upper half of
  the mixed layer and above it. The held θ was laid by height above GUST's ground (`docs/solution_algorithm.rst`,
  frozen θ with its buoyancy force, limits), which can give a horizontal force in the capping layer that HRRR does
  not have; part of the gain may come from it. The case does not separate the two.
- **solved_flux is the lowest in the three C2 hours** (0.31 to 0.39 m/s below the scheme), and steady in fewer than
  200 iterations. In the three C1 hours it is not steady after 2000 iterations: w reaches 6.6, 5.4 and 4.9 m/s, and
  the residual stays near 5e-3.
- frozen_buoyant does not help in the C2 hours: it is the highest in apr02_18 (the lowest boundary layer, z_i
  797 m), where w reaches 4.2 m/s.
- HRRR's heat flux at the ground adds 0.02 to 0.10 m/s with θ held, in every hour.
- Near the ground (5 to 50 m) frozen_buoyant is slower than HRRR, by 0.5 to 0.6 m/s in may10_19.
  Below HRRR's lowest level the score interpolates HRRR.
- Three hours per class are a small sample.
- may10_19 with frozen_buoyant (the test deck): rms 0.39 m/s, bias −0.18 m/s; below 500 m a speed bias of
  −0.36 m/s, a vector rms of 0.54 m/s and 0.2° of direction; steady after 194 iterations.

## CTest
`rt_hrrr_convective` (label `nightly`, MPI builds) runs `check_convective.py run`: the chain of may10_19 with the
recommended option, from the downloads in `$GUST_HRRR_CONVECTIVE_DATA` (default: this folder). Without them it is
skipped (exit status 77). `hrrr_outer_selftest` tests the tool itself on synthetic crops.
