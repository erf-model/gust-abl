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

The recommended setting is the one in the case decks: `turb.realizability = 0.309`, `turb.nut_relax = 0.5`,
`convection.hybrid_sigma_max = 0.1`, `convection.hybrid_height` = 2.5 times the terrain height (100 m ridge,
292 m hills, 290 m Askervein, 30 m Bolund, 1965 m Big Southern Butte). The runs without the sensor set
`turb.realizability = 0`, `turb.nut_relax = 1` and `convection.hybrid_sigma_max = 0`. Bolund uses the water/land
roughness map (`bc.surface_z0_file`, PR #59) in place of the decks' single roughness.

The Kato–Launder runs of Section 5.7 used an exploratory option (`turb.production`, commit dfb0443) that was
removed again in cf867cf because it did not change the separation; they are not repeated by the driver.

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
```bash
cmake -S . -B build-rel -DCMAKE_BUILD_TYPE=Release && cmake --build build-rel -j 8
python3 Tests/Paper/JWEIA2026/paper_runs.py prepare --out paper-runs --wntools ~/windninja-tools \
    --iz-data ~/data/ishihara_zhou --rushil-data ~/data/rushil
python3 Tests/Paper/JWEIA2026/paper_runs.py run --gust build-rel/gust --out paper-runs --np 8          # all groups
python3 Tests/Paper/JWEIA2026/paper_runs.py run --gust build-rel/gust --out paper-runs --np 8 M        # main only
python3 Tests/Paper/JWEIA2026/paper_runs.py stats --out paper-runs
python3 Tests/Paper/JWEIA2026/plot_paper.py --runs paper-runs --outdir figures
python3 Tests/Paper/JWEIA2026/plot_sitemaps.py --runs paper-runs --out figures/sitemaps.png
python3 Tests/Paper/JWEIA2026/bolund_masts.py --runs paper-runs --wntools ~/windninja-tools
python3 Tests/Paper/JWEIA2026/plot_hhabl.py extern_erf figures/hhabl.png \
    "without sensor=paper-runs/runs/H_default/turb<last>" "with sensor=paper-runs/runs/H_sensor/turb<last>"
```
Each run writes `paper-runs/runs/<group>_<variant>_<case>/` (deck, log, plotfiles) and, for the field sites,
`paper-runs/samples/<run>.csv`. A run with a `done` file is skipped, so the driver can be restarted. Stalled runs
stop at 4000 iterations (`--max-iterations`). The whole set is about 130 runs; most take minutes on 8 ranks, the
Big Southern Butte runs about half an hour.

## Expected Results
EXPECTED_RESULTS

## CTest
None: the runs take hours in total. The options they use are covered by the unit and regression tests of
`Tests/Real_Terrain`.
