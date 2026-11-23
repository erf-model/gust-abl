# HRRR_Neutral

## Purpose
Options for neutral and conventionally neutral HRRR hours (class C3 of `tools/hrrr_regime.py`,
`docs/hrrr_regimes.rst`) as the outer level of a GUST nest at the AWAKEN site north of Enid, Oklahoma. The case uses
the chain of `../HRRR_Outer` (`tools/hrrr_outer.py`, `docs/hrrr_outer.rst`): `hrrr_raw` (56 km, 1000 m cells),
`parent` (52 km, 500 m cells, mass-consistent first guess) and the nest (25 km, 125 m cells, four nest faces). It
starts from the scheme that `hrrr_outer.py case` writes (forecast hour 1, top 3 km, the adjusted parent, θv,
Coriolis on, θ solved with a neutral ground, a Rayleigh layer on w and θ, a 100 m length cap) and tries the options
below on four C3 hours, two of them under a strong capping inversion. The score is the nest against HRRR, as in
`../HRRR_Outer`. It is not a comparison with the AWAKEN observations, and none are used. This case is beyond the
development plan.

## Files
| File | Role |
|---|---|
| `case.yaml` | The meshes and settings of `../HRRR_Outer/case.yaml` (nest half-width 12500 m, 125 m cells, 42 layers from 10 m, ratio 1.08, top 3042 m, z0 = 0.09 m, 4 ranks); `options.py` sets the hour's crop and terrain. |
| `prepare.py` | Downloads `awaken.asc` (`tools/fetch_terrain.py`, 36.05–36.68 N, 97.80–97.01 W) and the crops of the five hours below (`tools/hrrr_outer.py fetch`, 45 km half-width, which also covers the 40 km gradient box of `hrrr_regime.py column`); `--neighbours` adds the hours before and after, for the classifier. |
| `options.py` | Writes the chain of each hour, its nest options and its two top chains, and `run_options.sh`, which runs them one at a time and scores each nest. |
| `wide_top.py` | `tools/hrrr_outer.py` with its top limits widened (1000 to 6000 m), for the two top options only. |
| `score_options.py` | The score table of the options (Markdown). |
| `check_neutral.py` | The test deck: the hour apr05_13 with the recommended option, checked against the numbers of November 2026 (CTest `rt_hrrr_neutral`). |

The hours (valid time, UTC; each is forecast hour 1 of the cycle an hour before). The class and modifiers are from
`hrrr_regime.py catalogue` on the columns of the crops with the hours before and after (B: baroclinic, T: wind
tendency); z_i is HRRR's boundary-layer height and N² the stability of the 500 m above it.

| Tag | Valid | Class | L (m) | z_i (m) | N² above z_i (s⁻²) | Note |
|---|---|---|---|---|---|---|
| apr02_23 | 2023-04-02 23:00 | C3 B+T | −629 | 1553 | 8.9e-5 | weak cap |
| apr05_13 | 2023-04-05 13:00 | C3 B | −575 | 536 | 6.7e-4 | strong cap; the night's jet above it (16.9 m/s at 792 m); the test deck |
| apr18_23 | 2023-04-18 23:00 | C3 B | −837 | 1249 | 2.2e-4 | windy (u* 1.09 m/s) |
| apr20_06 | 2023-04-20 06:00 | C3 T | 1198 | 444 | 1.8e-4 | refused by `hrrr_outer.py` as cloudy (cloud water 0.02 g/kg at 862 m) |
| apr20_13 | 2023-04-20 13:00 | C3 B+T | −603 | 530 | 5.0e-4 | strong cap; the night's jet above it (17.7 m/s at 619 m) |

They were found by a screen of the HRRR surface fields at the site for April 2023 and then classified with their
neighbours; hours whose neighbours make them a transition (C7, C8) were left out. `prepare.py` still downloads
apr20_06, and `options.py` reports its refusal, so the options run on four hours.

The options (the nest options share the chain's parent; the top options are chains of their own):

| Option | θ | Coriolis | Top | Rayleigh layer on w and θ |
|---|---|---|---|---|
| scheme | solved, neutral ground | on | 3042 m | yes |
| frozen | frozen (`simplec.theta = frozen`, the parent's) | on | 3042 m | w only (θ held) |
| no_coriolis | solved | off | 3042 m | yes |
| frozen_no_coriolis | frozen | off | 3042 m | w only (θ held) |
| no_rayleigh | solved | on | 3042 m | no |
| top1500 | solved | on | 1500 m (33 layers) | yes |
| top5000 | solved | on | 5000 m (48 layers) | yes |

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
The class table above (the crops of the hours before and after, and their columns):
```bash
python3 prepare.py --neighbours --columns columns
python3 ../../../tools/hrrr_regime.py catalogue columns/*.json
```
The options take about three hours on 4 ranks of a laptop (4 hours × 7 nests, and the 24 coarser runs).

## Expected Results
2026-11-23, 4 ranks: every run is steady except frozen and frozen_no_coriolis in apr20_13 (4000 iterations).
Vector rms (m/s) of the nest minus HRRR below 500 m above the ground (the full table, with bias, rms, direction, w,
k and iterations, is in `docs/hrrr_outer.rst`, Neutral hours; `score_options.py` prints it):

| Option | apr02_23 | apr05_13 | apr18_23 | apr20_13 | Mean |
|---|---|---|---|---|---|
| scheme | 1.28 | 1.73 | 1.42 | 2.24 | 1.67 |
| frozen | 0.84 | 2.48 | 1.59 | 5.64 | 2.64 |
| **no_coriolis** | **1.08** | **1.17** | **1.97** | **1.41** | **1.41** |
| frozen_no_coriolis | 1.06 | 1.53 | 1.85 | 5.06 | 2.38 |
| no_rayleigh | 1.44 | 1.62 | 1.36 | 2.16 | 1.65 |
| top1500 | 1.49 | 1.40 | 1.78 | 1.61 | 1.57 |
| top5000 | 1.52 | 1.34 | 1.33 | 2.33 | 1.63 |

- **no_coriolis is the recommended option for C3 hours under a strong cap** (apr05_13, apr20_13: 0.56 and
  0.83 m/s lower than the scheme). The evening hours are mixed: 0.20 m/s lower in apr02_23, 0.55 m/s higher in
  apr18_23. The gain is in speed (direction moves by at most 1.4°); the reason is not tested. Four hours are a
  small sample.
- frozen is the lowest in the weak-cap hour and the highest under a strong cap, where w reaches 6 and 14 m/s.
- The hours are weakly unstable in HRRR (L −575 to −837 m); every option uses a neutral ground.
- The Rayleigh layer has a weak effect (−0.11 to +0.16 m/s). The tops change the score by up to 0.64 m/s with no
  order common to the hours.
- apr05_13 with no_coriolis (the test deck): rms 0.93 m/s, bias −0.48 m/s; below 500 m a speed bias of −1.06 m/s,
  a vector rms of 1.17 m/s and 0.1° of direction; steady after 447 iterations.

## CTest
`rt_hrrr_neutral` (label `nightly`) runs `check_neutral.py run`: the chain of apr05_13 with the recommended option,
from the downloads in `$GUST_HRRR_NEUTRAL_DATA` (default: this folder). Without them it is skipped (exit status
77). `hrrr_outer_selftest` tests the tool itself on synthetic crops.
