# HRRR_Outer

## Purpose
HRRR as the outer level of a GUST run at the AWAKEN site (American WAKE experimeNt) north of Enid, Oklahoma. One
hour of the High-Resolution Rapid Refresh forecast drives a nest of 25 km with 125 m cells through two coarser GUST
levels: `hrrr_raw` (56 km, 1000 m cells) takes HRRR's u, v, w, θv and k as they are; `parent` (52 km, 500 m cells)
adjusts that field to be mass-consistent on its own terrain; the nest solves the steady flow with its four lateral
faces driven by the parent. `tools/hrrr_outer.py` (`docs/hrrr_outer.rst`) writes the chain from `case.yaml` and
scores the nest against HRRR. The case checks the day hour against the numbers of the study this tool replaces,
and runs a matrix of seven variants over the hours as the baseline for later work. It is not a comparison with
the AWAKEN observations, and none are used.

## Files
| File | Role |
|---|---|
| `case.yaml` | The day hour (cycle 2023-05-10 19 UTC, forecast hour 1): nest half-width 12500 m, 125 m cells, 42 layers from 10 m (ratio 1.08, top 3042 m), z0 = 0.09 m, BTF smoothing, 4 ranks; `hrrr.file` is the crop. |
| `prepare.py` | Downloads `awaken.asc` (`tools/fetch_terrain.py`, 36.05–36.68 N, 97.80–97.01 W) and the crops of the five hours (`tools/hrrr_outer.py fetch`, 32 km half-width). |
| `check_day.py` | Checks the day chain's `score.json` against the November 2026 numbers (rms 0.47 m/s, vector rms below 500 m 0.57 m/s, steady). |
| `variants.py` | Writes the chain of each hour and six variants of its nest (A–F below) beside it, and `run_matrix.sh`. |
| `score_variants.py` | The score table of the matrix (Markdown). |

The hours (valid time, UTC): day20 (2023-05-10 20:00, convective), dawn12 (2023-05-10 12:00), jet05 (2023-05-10
05:00, low-level jet), eve01 (2023-05-11 01:00, evening transition), night09 (2023-05-11 09:00, stable). Each is
forecast hour 1 of the cycle an hour before; the forecast hour 0 analysis is refused by the tool (it carries
almost no vertical motion). jet05 is refused as cloudy (cloud water inside the raw domain below its top: GUST is
dry), so the matrix has four hours; `prepare.py` still downloads it, and `variants.py` reports the refusal.

The variants share each hour's parent:

| Variant | θ | Coriolis | Geostrophic forcing | Other |
|---|---|---|---|---|
| **S** | solved | on | no | the scheme as `hrrr_outer.py case` writes it (the chain's own nest) |
| A | frozen | off | no | |
| B | frozen | on | no | |
| C | frozen | on | HRRR's G (200–1000 m mean) | |
| D | solved | on | HRRR's G | |
| E | frozen | on | HRRR's G | relaxation zone of 8 cells inside the nest faces |
| F | frozen | on | HRRR's G | HRRR's surface heat flux |

## Data and acknowledgements
The meteorological fields are from the High-Resolution Rapid Refresh (HRRR) model of NOAA's Global Systems
Laboratory, run operationally by the National Centers for Environmental Prediction: the one-hour forecasts (f01) of
the native-level (`wrfnat`) and surface (`wrfsfc`) files, read from the NOAA Open Data Dissemination archive on AWS
(https://registry.opendata.aws/noaa-hrrr-pds/) with Herbie (Blaylock, https://github.com/blaylockbk/Herbie), which
downloads only the byte ranges of the fields used. We thank NOAA and the HRRR developers for making the forecasts
free to use, Amazon Web Services for hosting the archive, and Brian Blaylock for Herbie. NOAA open data carry no
restriction on use or redistribution; the crops (about 4 MB an hour) are not redistributed here, and `prepare.py`
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
python3 ../../../tools/hrrr_outer.py case case.yaml
cd day20 && GUST=<build>/gust ./run.sh && cd ..
python3 check_day.py
```
The matrix (4 hours × 7 nests = 28; a few hours on 4 ranks of a laptop; `--out` can be any folder):
```bash
python3 variants.py --out matrix
GUST=<build>/gust matrix/run_matrix.sh
python3 score_variants.py matrix
```

## Expected Results
- The day chain: `hrrr_raw` and `parent` write their first guess only; the nest is steady in about 240 iterations.
  `check_day.py` exits with 0 (2026-11-22: rms 0.48 m/s, bias −0.03 m/s, below 500 m a speed bias of −0.04 m/s, a
  vector rms of 0.58 m/s and 0.3° of direction; w rms 0.085 m/s, largest |w| inside the margin 0.79 m/s; k / k_HRRR
  0.17; steady after 247 iterations).
- The matrix: see `docs/hrrr_outer.rst` for the table of 2026-11-22.

## CTest
None: the case needs the downloads. `hrrr_outer_selftest` tests the tool on synthetic crops, checkpoints and
plotfiles.
