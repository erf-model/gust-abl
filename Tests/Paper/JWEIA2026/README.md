# JWEIA2026

## Purpose
Every run, table and figure of the GUST paper on the realizability sensor (Gopalan, Quon, Zhang, Lattanzi and
Almgren, submitted to the *Journal of Wind Engineering and Industrial Aerodynamics*, 2026). The paper compares
GUST - CFD without and with the realizability sensor against wind-tunnel data (the Ishihara–Zhou ridge, the RUSHIL
hills) and field data (Askervein, Bolund cases 1, 3 and 4, Big Southern Butte). The cases themselves, their decks
and their prepare / compare / sample scripts are in `Tests/Real_Terrain/IshiharaZhou`, `Tests/Real_Terrain/RUSHIL`
and `Tests/Real_Terrain/WindNinja2019`. This folder holds the driver that runs all the variants of the paper one
after another, the empty-domain decks, and the scripts that draw the paper's figures.

## Files
| File | Role |
|---|---|
| `paper_runs.py` | `prepare` (all terrain, inflow, mast and roughness files from the downloaded data), `run` (every variant, one at a time, each with its own deck copy), `stats` (windward / ridge / lee errors of the field runs). Groups: M main comparison, A limit and blending, H empty domain, S sensitivity, C C_r = 0.245, E grids. |
| `inputs_hhabl_base`, `inputs_hhabl_consistent` | Empty-domain test (Section 5.1): flat ground with the Askervein vertical grid; log-law inflow with a slip top and damping, and the consistent set-up (uncapped length, lid moving at the inflow's top speed). |
| `hhabl_drift.py` | Drift of U, k and ν_t between the inflow and the hill position in the empty-domain runs. |
| `plot_paper.py` | Figures of the comparison: observations, GUST - CFD without and with the sensor (field sites, ridge, hills). |
| `plot_hhabl.py` | Figure of the empty-domain profiles. |
| `plot_sitemaps.py` | Site maps of Askervein, Big Southern Butte and Bolund with the masts and wind directions. |
| `bolund_masts.py` | Bolund mast by mast (Appendix C): speed and k at 2 m and 5 m at every sonic mast. |
| `flatten_outflow_edges.py` | Gives the Big Southern Butte terrain no slope normal to the outflow faces over the last 680 m (`prepare` writes `bsb_flatedge.xyz`). |

The recommended setting is the one in the case decks: `turb.realizability = 0.309`, `turb.nut_relax = 0.5`,
`convection.hybrid_sigma_max = 0.1`, `convection.hybrid_height` = 2.5 times the terrain height (100 m ridge,
292 m hills, 290 m Askervein, 30 m Bolund, 1965 m Big Southern Butte). The runs without the sensor set
`turb.realizability = 0`, `turb.nut_relax = 1` and `convection.hybrid_sigma_max = 0`. Bolund uses the water/land
roughness map (`bc.surface_z0_file`, PR #59) in place of the decks' single roughness.

The Kato–Launder runs of the paper used an exploratory option (`turb.production`) of a development version that was
removed before release because it did not change the separation; they are not repeated by the driver.

## Data and acknowledgements
Downloaded by the user (not in the repository):
- **windninja-tools** (N. Wagenbrenner, <https://github.com/nwagenbrenner/windninja-tools>): the Askervein masts,
  the Bolund masts and terrain (`Bolund_Measurements.zip`, `bolund_hill_extended.tif`) and the Big Southern Butte
  terrain (`big_butte.asc`). See `Tests/Real_Terrain/WindNinja2019/README.md`.
- **Ishihara and Zhou (2025)**, the supplementary spreadsheets mmc1.xls–mmc4.xls of the ridge wind-tunnel data. See
  `Tests/Real_Terrain/IshiharaZhou/README.md`.
- **ERCOFTAC Classic Collection, case AC5-05 (RUSHIL)**, Khurshudyan, Snyder and Nekrasov (1981), US EPA Fluid
  Modeling Facility. See `Tests/Real_Terrain/RUSHIL/README.md`.
- The tunnel points of Hills 3 and 5 as given by Manning, Hancock and Whiting (2010), *Wind Engineering* 34, 477,
  digitized by the user (`plot_paper.py --manning`); not distributed.

In the repository (CC BY 4.0): the Big Southern Butte observations of Wagenbrenner et al. (2019), *Atmosphere* 10,
672, doi:10.3390/atmos10110672, and the Big Southern Butte station data, Zenodo, doi:10.5281/zenodo.8226638 (in
`Tests/Real_Terrain/WindNinja2019`).

We thank Natalie Wagenbrenner for windninja-tools, the Big Southern Butte database and the open-access paper; the
Missoula Fire Sciences Laboratory for the Big Southern Butte campaign; DTU Wind Energy and the Bolund team for the
Bolund data; the Askervein Hill project for the Askervein data; Takeshi Ishihara and Tong Zhou for publishing the
ridge data; and ERCOFTAC, the QNET-CFD Knowledge Base and the US EPA Fluid Modeling Facility for the RUSHIL data.

## How to run
Prerequisites: Python 3.8 or later with numpy, scipy, matplotlib and xlrd (`ishihara_zhou.py prepare`), and an MPI
launcher (`--mpiexec`, default `mpirun -np`; empty for a serial build).

```bash
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release && cmake --build build-rel -j 8
python3 Tests/Paper/JWEIA2026/paper_runs.py prepare --out paper-runs --wntools ~/windninja-tools \
    --iz-data ~/data/ishihara_zhou --rushil-data ~/data/rushil
python3 Tests/Paper/JWEIA2026/paper_runs.py run --gust build-rel/gust --out paper-runs --np 8          # all groups
python3 Tests/Paper/JWEIA2026/paper_runs.py run --gust build-rel/gust --out paper-runs --np 8 M        # main only
python3 Tests/Paper/JWEIA2026/paper_runs.py stats --out paper-runs --wntools ~/windninja-tools
python3 Tests/Paper/JWEIA2026/plot_paper.py --runs paper-runs --outdir figures
python3 Tests/Paper/JWEIA2026/plot_sitemaps.py --runs paper-runs --out figures/sitemaps.png
python3 Tests/Paper/JWEIA2026/bolund_masts.py --runs paper-runs --wntools ~/windninja-tools
python3 Tests/Paper/JWEIA2026/plot_hhabl.py extern_erf figures/hhabl.png \
    "GUST - CFD without sensor=$(ls -d paper-runs/runs/H_consistent/turb[0-9]* | tail -1)" \
    "GUST - CFD with sensor=$(ls -d paper-runs/runs/H_consistent_sensor/turb[0-9]* | tail -1)"
python3 -I Tests/Paper/JWEIA2026/hhabl_drift.py $(ls -d paper-runs/runs/H_consistent_sensor/turb[0-9]* | tail -1) \
    extern_erf --law log_law --ustar 0.6278 --stations 0 2988 --heights 2 10
```
Each run writes `paper-runs/runs/<group>_<variant>_<case>/` (deck, log, plotfiles) and, for the field sites,
`paper-runs/samples/<run>.csv`. A run with a `done` file is skipped, so the driver can be restarted. Stalled runs
stop at 4000 iterations (`--max-iterations`). The whole set is about 130 runs; most take minutes on 8 ranks, the
Big Southern Butte runs about half an hour.

## Expected Results
All runs of the paper reproduced on this code (`paper_runs.py run`, 8 ranks, 2026-10-25). Field sites: lee / all
MAPE; wind tunnel: extent of reversed flow on the near-surface line and rms error in % of the reference velocity.

| Case | Without sensor | With sensor | Limit only | Upwind, no limit | Upwind, limit |
|---|---|---|---|---|---|
| Smooth ridge | none; 20 % | 42-170 m; 7 % | 41-179 m; 6 % | none | 54-119 m |
| Rough ridge | 68-102 m; 23 % | 47-173 m; 15 % | 47-178 m | 75-95 m | 51-138 m |
| RUSHIL Hill 3 | 0.87a-0.98a; 20 % | 0.58a-1.45a; 17 % | 0.58a-1.48a | none | 0.62a-1.23a |
| RUSHIL Hills 5 / 8 | none / none | none / none | none / none | none / none | none / none |
| Askervein | 25 / 8 % | 14 / 8 % | 13 / 8 % | 24 / 9 % | 14 / 8 % |
| Bolund case 1 | 54 / 18 % | 36 / 14 % | 33 / 13 % | 54 / 18 % | 44 / 15 % |
| Bolund case 3 | 20 / 8 % | 18 / 8 % | 25 / 9 % | 28 / 10 % | 10 / 6 % |
| Bolund case 4 | 54 / 15 % | 81 / 21 % | 82 / 21 % | 35 / 10 % | 61 / 16 % |
| Big Southern Butte | not converged (r 2.2e-3) | 32 / 28 % | not converged (r 3.5e-3) | 38 / 30 % | 42 / 32 % |

- Big Southern Butte with no normal slope at the outflow faces (`A_*_bsb_flatedge`): limit only 32 / 28 % (793
  iterations), with sensor 32 / 28 % (852; 853 on the original terrain).
- C_r (Appendices A, B), smooth ridge / Hill 3 / Bolund case 1 lee: 0.245: 31-279 m / 0.42a-2.13a / 2 %;
  0.33: 47-145 m / 0.63a-1.32a / 40 %; 0.346 (Durbin): 51-130 m / 0.67a-1.23a / 43 %; 0.35: 53-127 m / 0.67a-1.21a /
  43 %; 0.40: 68-100 m / 0.77a-1.06a / 47 %. Big Southern Butte at 0.245: not converged; at 0.346: 36 / 29 %.
- Empty domain (consistent set-up, 3 km): U(2 m) -1.2 % and u* +1.0 % with the sensor; -1.9 % / +2.1 % at C_r = 0.346;
  nu_t -43 % and u* -18.7 % at C_r = 0.245.
- Grids, coarse / base / fine with nested levels and the same top (sensor; limit only): smooth ridge 49-152, 42-170,
  45-156 m; 48-164, 41-179, 45-160 m. Hill 3 0.54a-1.52a, 0.58a-1.45a, 0.59a-1.41a; 0.54a-1.59a, 0.58a-1.48a,
  0.60a-1.42a (the fine limit-only run stops at a fixed state, r 2.6e-2; the sensor run converges). Bolund case 1 with
  the horizontal spacing doubled / halved: lee 46 / 36 / 35 %, all 19 / 14 / 13 %.

Differences from the paper's first runs (made before GUST PR #57, the ghost fill over terrain as a fixed point): the
smooth ridge without the sensor converges (it stopped at r 1.6e-2 before), the extents of reversed flow change by at
most 1 m or 0.01a, the RMSE at the field sites by at most 0.01 m/s and the MAPE by at most one percentage point.

## CTest
None: the runs take hours in total. The options they use are covered by the unit and regression tests of
`Tests/Real_Terrain`.
