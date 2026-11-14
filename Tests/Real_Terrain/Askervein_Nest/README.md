# Askervein_Nest

## Purpose
The one-way nesting validation of `Tests/Canonical_RANS/Nest_Validation` on measured terrain (sub-phase 9.3, beyond
plan rev. 6). The whole 6 km Askervein domain of `Askervein_210` runs on 40 m cells (the parent) and on 20 m cells
(the reference). A 20 m patch of half-width D centred on the hilltop HT refines the parent at 2 2 1 with four `nest`
faces. The measure is the wind speed 10 m above the ground at the 30 masts of lines A and AA: the patch against the
reference, as a function of D, beside the parent against the reference.

## Files
| File | Role |
|---|---|
| `inputs_askervein_domain` | The whole domain: `Askervein_210/inputs_steady_askervein_210` without its cells and output names, line LU-SGS, two-level grid sequencing. |
| `inputs_askervein_patch` | A patch: four `nest` faces, first guess from the parent, `terrain.nest_blend_width = 200` (the 20 m ground meets the parent's 40 m ground on the faces), line LU-SGS without sequencing. |
| `askervein_nest.py` | Runs a whole domain (`run`) or a patch (`nest`), samples the masts, compares (`compare`), or does the study (`study`). Pure Python; uses `tools/askervein_lines.py` and the run helpers of `Nest_Validation/nest_validation.py`. |

The terrain (`extern_erf/askervein.txt`), the inflow table and the mast positions
(`Askervein_210/askervein_210_inflow.txt`, `askervein_lines.csv`) are those of `Askervein_210`.

## Data and acknowledgements
No file is added by this case. It uses, as `Askervein_210` does:
- the Askervein terrain and inflow of ERF's `Exec/CanonicalTests/Real_Terrain/Askervein` (`extern_erf/`, recorded in
  `extern_erf/EXTRACTED.md`). Thanks to the ERF developers;
- the mast positions of lines A and AA (P. A. Taylor and H. W. Teunissen, "The Askervein Hill project: overview and
  background data", Boundary-Layer Meteorology 39, 1987, doi:10.1007/BF00121863), as tabulated in N. Wagenbrenner's
  windninja-tools (github.com/nwagenbrenner/windninja-tools, `data/askervein/askervein_obs.zip`), used with her
  permission. Thanks to Natalie Wagenbrenner and to the Askervein field teams.

All three files are already redistributed in this repository; nothing has to be downloaded.

## How to run
```bash
ctest --test-dir <build> -R "askervein_nest"
```
The study, in a scratch directory (1.6 million cells in the reference):
```bash
python3 askervein_nest.py study -- mpiexec -n 4 <build>/gust
```
Single runs:
```bash
python3 askervein_nest.py run par40 150 -- mpiexec -n 4 <build>/gust inputs_askervein_domain
python3 askervein_nest.py run ref20 300 -- mpiexec -n 4 <build>/gust inputs_askervein_domain
python3 askervein_nest.py nest nest_d1400 par40 1400 -- mpiexec -n 4 <build>/gust inputs_askervein_patch
python3 askervein_nest.py compare nest_d1400 ref20 par40 --ratio 0.5
```

## Targets
Fixed in this file before any run was compared with the reference. When they were set, the 40 m parent and the
D = 1000 m patch had been run to the steady verdict (103 and 95 iterations); the 20 m reference had not been run in
this study and no mast had been sampled.

| | Target |
|---|---|
| A1 | D = 1400 and 2000 m: RMS over the 30 masts of the relative difference of the 10 m speed from the 20 m run <= 1/2 of the 40 m parent's |
| A2 | The RMS does not grow from D = 1000 to 1400 to 2000 m (5 % allowance) |
| A3 | Every patch reaches the steady verdict, and the log prints no net-flux warning (net <= 1e-3 of the flux entering) |

A target that is missed is reported as missed, with the number.

## Expected Results
To be filled from the study.

## CTest
To be filled.
