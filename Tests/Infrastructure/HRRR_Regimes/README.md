# HRRR_Regimes

## Purpose
The stratification class of HRRR hours at the AWAKEN site, from `tools/hrrr_regime.py` (`docs/hrrr_regimes.rst`).
Sixteen hours of 10–11 May 2023 at 36.362 N, 97.405 W (an approximate centre of the AWAKEN area north of Enid,
Oklahoma, the centre of `Tests/Real_Terrain_Python/AWAKEN`). They cover a night with a low-level jet, a morning
transition, a convective afternoon, an evening transition and a second, very stable night. The check is a regression pin:
it pins the class of each hour and its diagnostics with the default thresholds. The thresholds are starting values; this case is not
evidence that they are right.

## Files
| File | Role |
|---|---|
| `data/hrrr_<valid>_f01.json` | The columns of the 16 hours: the means over a 25 km square of HRRR's first 20 hybrid levels (to about 6.4 km), the gradient of θv over an 80 km square, and the surface fields. Written by `hrrr_regime.py fetch`. |
| `check_classes.py` | Classifies the 16 hours and checks them against the values it pins (class, modifiers, L, z_i/L, Ri_b, jet fall-off, geostrophic shear, wind tendency) and L against its closed form worked from each column's surface fields. |

## Data and acknowledgements
The columns are made from the High-Resolution Rapid Refresh (HRRR) model of NOAA's Global Systems Laboratory, run
operationally by the National Centers for Environmental Prediction: the one-hour forecasts (f01) of the native-level
(`wrfnat`) and surface (`wrfsfc`) files, read from the NOAA Open Data Dissemination archive on AWS
(https://registry.opendata.aws/noaa-hrrr-pds/) with Herbie (Blaylock, https://github.com/blaylockbk/Herbie). We thank
NOAA and the HRRR developers for making the forecasts free to use, Amazon Web Services for hosting the archive, and
Brian Blaylock for Herbie. NOAA open data carry no restriction on use or redistribution. The 16 small columns
derived from them (64 kB) are redistributed here; the GRIB files are not (about 190 MB an hour).

## How to run
```bash
cd Tests/Infrastructure/HRRR_Regimes && python3 check_classes.py
```
To make the columns again (needs herbie, xarray and cfgrib, and an internet connection; about 3 GB of GRIB subsets):
```bash
python3 tools/hrrr_regime.py fetch 2023-05-10T04:00Z 2023-05-10T05:00Z ... --lat 36.362 --lon -97.405 \
    --save <grib folder> --out-dir data --jobs 4
```
with the 16 valid times of the table below.

## Expected Results
| Valid (UTC) | Class | Modifiers | L (m) | z_i/L | Note |
|---|---|---|---|---|---|
| 10 May 04 | C5 | B, T | 24 650 | 0.013 | Jet nose 14.0 m/s at 866 m, 3.3 m/s fall-off; weak cooling (q_v = −4e-4 K m/s, sensible −2.6 W/m² against latent +30 W/m²) |
| 05 | C5 | B, T | 770 | 0.42 | Jet nose 16.6 m/s at 455 m, 3.0 m/s fall-off |
| 06 | C5 | T | 361 | 0.94 | |
| 11 | C6 | B, T | 13.1 | 4.7 | |
| 12 | C8 | B, T | 122 | 0.67 | Surface flux turns upward by 13 UTC |
| 13 | C8 | T | −209 | −1.3 | |
| 19 | C1 | B | −48.8 | −30.8 | |
| 20 | C1 | | −66.6 | −24.6 | |
| 21 | C1 | B | −87.1 | −20.1 | |
| 23 | C1 | B, T | −156 | −12.2 | |
| 11 May 00 | C7 | B, T | −510 | −2.6 | Surface flux turns downward by 01 UTC |
| 01 | C7 | B, T | 420 | 2.4 | |
| 02 | C4 | B, T | 434 | 0.84 | |
| 08 | C6 | B, T | 55.6 | 1.8 | |
| 09 | C6 | B, T | 63.6 | 1.8 | |
| 10 | C5 | B, T | 1.75 | 55 | Jet with 2.6 m/s fall-off; q_v = −7e-4 K m/s, u* = 0.026 m/s |

The modifiers B and T are each set on 13 of the 16 hours with the starting thresholds (2 m/s per km and 1 m/s
per hour); they are to be set again from a longer catalogue. Some classes sit close to a threshold: the jet
fall-off at 11 May 09 UTC is 1.99 m/s against 2 m/s (with `jet_drop` = 1, 08 and 09 UTC become C5); z_i/L at 08 and
09 UTC is 1.8 against 1 (`very_stable_ratio` = 2 gives C4); −z_i/L at 23 UTC is 12.2 against 10. The two
transitions keep their class for `flux_eps` from 0.001 to 0.02, and the class counts do not change for
`neutral_length` from 500 to 2000 m. The L of the two calmest hours (u* of 0.026 and 0.065 m/s) is uncertain:
L varies as u*³, so 10 % in u* is 30 % in L. The hours at 05, 12 (as C6 without its neighbours),
20, 01 and 09 UTC are the hours run with an HRRR outer level in November 2026; the labels given to them then were
by eye.

The column of 05 UTC made by `fetch` agrees with the one made by `hrrr_regime.py column` from the crop saved for the
earlier runs (50 levels): every profile and surface field equal, the θv gradient within 3e-8 K/m.

## CTest
`hrrr_regime_awaken` (label `unit`; registered when Python has numpy). The tool's own tests are
`hrrr_regime_selftest`.
