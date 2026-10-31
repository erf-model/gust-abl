# WINSENT

## Purpose
The six steady cases of WINSENT Phase I: high LAI (5.0) and low LAI (0.7), each neutral, stable and unstable. The
flow runs from the valley lidar L140 up the forested Albtrauf escarpment to the plateau mast MMNW. GUST runs the same
domain, terrain, forest and levels as the ERF runs of the same cases, so the two codes can be compared at the
stations. The results are the GUST submission to the benchmark (one NetCDF file per case at MMNW).

The neutral cases solve θ (it stays at the sounding). The stable and unstable cases prescribe θ from a near-1-D ERF
column with a fixed surface heat flux (`simplec.theta = frozen`, `docs/solution_algorithm.rst`, "Prescribed
potential temperature"), with buoyancy in the turbulence closure only, the steady form of van der Laan et al. (2024).
With θ solved, or with buoyancy in the momentum equations, these four cases did not converge.

## Files
| File | Role |
|---|---|
| `prepare_gust.py` | Turns one ERF case folder (from ERF's `Exec/CanonicalTests/Real_Terrain/WINSENT/prepare_winsent.py`) into GUST's tables; `--roughness` adds the roughness map, `--canopy-height 15` the uniform 15 m forest of the runs below. |
| `inputs_winsent_z0map` | The deck of all six runs: the roughness map of the kynema-sgf WINSENT runs (`roughness.xyz`). `steady.cfl_max`, `steady.max_iterations` and, for the stratified cases, `simplec.theta` go on the command line. |
| `make_z0_only.py` | A copy of a neutral case folder with the forest as roughness and no canopy keys. |
| `inputs_winsent` | The same deck with ERF's uniform z0 = 0.1 m (the first neutral runs). |
| `inputs_winsent_column` | The neutral inflow column: flat ground, Coriolis at 48.7° N, the geostrophic wind found for a target wind at 300 m. |
| `column_inflow.py` | Runs the neutral column (`--scale G`: G times the benchmark power law at 300 m) and writes the inflow `z u v theta k`, turned to the lidar's direction. |
| `lidar_mean.py` | The vector-mean lidar profile of one class, for the direction of the stratified inflows. |
| `erf_column.py` | A near-1-D ERF column with a fixed surface heat flux, 9 h; writes the inflow and frozen θ of a stratified case, and logs the bulk Richardson number R_B between 46 and 73 m. |
| `sample_stations.py` | Speed and direction at MMNW and L140 at ERF's station heights, from a GUST or an ERF plotfile. |
| `plot_lidar.py` | A run against the lidar at L140 (bias and RMS over 100–400 m), and its MMNW profile. |
| `make_submission.py` | The benchmark's NetCDF file for one case, from a converged run. |

Tried and not used for the submission (kept so the record is complete):

| File | What it did, and why it was not used |
|---|---|
| `correct_inflow.py` | Scaled the inflow at each height until the run matched the lidar at L140. The corrected inflow was flat from 50 to 400 m with uniform k, which is not a realistic upstream profile. |
| `column_rib_sweep.py` | Swept the Obukhov length of GUST's own stratified column. A steady column needs a Rayleigh layer on θ as its heat sink, and its R_B did not reach the class medians. |
| `stable_1d.py` | amrTerrain's 1-D solver with a cooling rate or a heat flux. Its surface heat flux was 3–4 % of the requested one (the bottom face's diffusivity is u* κ z0/φ_m), so its columns stayed near neutral. |

The roughness map is kynema-sgf's `winsent.roughness` (10 m, z0 0.001–0.5 m, 0.1 m under the forest). Its frame
is this case's frame turned by −1° (292° against 293°) with MMNW at (1410, −400) m, found by matching
`winsent.amrwind` to ERF's `terrain.txt` (6.8 m RMS difference). `prepare_gust.py` averages ln z0 over each 32 m
cell.

## Set-up
- **Mesh:** ERF's blended domain, 416 × 320 × 32 cells (dx = 32 m), 32 stretched levels to 2400 m, frame rotated so
  +x points downwind (wind from 293°), origin at MMNW. Grid sequencing over three levels (128, 64, 32 m).
- **Forest:** the benchmark's LAI map (5.0 or 0.7) with a uniform height of 15 m, Cd 0.15, the canopy k terms. 15 m
  is a fitted height, not a measured one: the neutral sweep found it gives the lowest L140 error (the map's 25 m
  gave a larger one).
- **Ground:** the roughness map; the deck's neutral wall law in all six cases (no surface heat flux is set; with θ
  frozen the stratification comes from the inflow column, not from the wall).
- **Inflow (xlo):** a column profile `z u v theta k`, turned at each height to the lidar's direction. Outflow on xhi,
  periodic in y, slip top, Rayleigh layer on w above 1900 m.
- **Solver:** SIMPLEC with LU-SGS lines, `steady.cfl_max = 2` (the dense canopy on the steep slope,
  `docs/canopy.rst`), realizability 0.309, stop at the default tolerances. These are the settings of the submitted
  runs. On the current step `steady.dtau_floor=1` replaces `steady.cfl_max=2` (cap 8): neutral 221 to 226 iterations,
  stable 211, and the neutral high LAI case keeps its L140 score (100–400 m, bias 0.00, RMS 0.49 m/s).

| Case | Inflow column | Geostrophic wind (m/s, case frame) | Surface flux (K m/s) | Column R_B (46–73 m) | `simplec.theta` |
|---|---|---|---|---|---|
| High LAI, neutral | GUST, `column_inflow.py --scale 1.2` (power law 1.033 h^0.323) | 7.717, −1.644 | 0 | 0 | `solve` |
| Low LAI, neutral | GUST, `column_inflow.py --scale 1.2 --powerlaw 0.423 0.479` | 7.691, −1.635 | 0 | 0 | `solve` |
| High LAI, stable | ERF, `erf_column.py` | 6.697, −1.427 | −0.004 | +0.122 | `frozen` |
| Low LAI, stable | ERF, `erf_column.py` | 5.963, −1.268 | −0.0025 | +0.097 | `frozen` |
| High LAI, unstable | ERF, `erf_column.py --lapse 0` | 7.256, −1.594 | +0.018 | −0.266 | `frozen` |
| Low LAI, unstable | ERF, `erf_column.py --lapse 0` | 7.646, −1.748 | +0.030 | −0.360 | `frozen` |

The ERF columns: 2 × 2 × 240 cells (10 m), the k-equation closure, MOST with z0 = 0.1 m and the fixed flux,
Coriolis at 48.7° N, 300 K to an 8 K inversion at 1000–1100 m, then 3 K/km (stable) or 0 (unstable), 9 h. The
flux and the geostrophic wind of each were chosen from a sweep so that R_B approaches the class median and the
column's 300 m wind stays near the neutral one's.

## Data and acknowledgements
The terrain, forest and lidar tables are made by ERF's WINSENT preparation script from these data. We thank the
benchmark's organisers and the providers of each dataset.
- WINSENT benchmark: A. Rettenmeier, S. Lott, H. Knaus, J. Hohaus, *WINSENT complex terrain benchmark - Phase I*,
  Zenodo, https://doi.org/10.5281/zenodo.23097333. The record is licensed CC BY 4.0; its files are released on
  registration with the benchmark (restricted access), so they are not copied here. The submission template
  (`template_phase1_inflow.nc`, in `submission_phase1_templates.zip`) is part of the record.
- Terrain: NASADEM Merged DEM Global 1 arc second V001, NASA LP DAAC,
  https://doi.org/10.5067/MEASURES/NASADEM/NASADEM_HGT.001 (NASA Earth science data, no restrictions on reuse).
- Forest mask: Copernicus CLCplus Backbone 2023, https://doi.org/10.2909/b0bd43c6-1fa1-4d88-9c45-98b13a95d0b2,
  Copernicus open-data licence.
- Inflow: the valley lidar (L140) 10-minute profiles of the benchmark record, and its class definitions (the bulk
  Richardson number at MMNW between 46 and 73 m).
- LAI: the benchmark's prescribed means, 5.0 (leafy) and 0.7 (leafless).
- Ground roughness: `winsent.roughness` and `winsent.amrwind` from the auxiliary files of the kynema-sgf WINSENT
  set-up (with `winsent_processing.ipynb`). The files name no authors or licence; we thank the kynema-sgf WINSENT
  team for them. They are not in this repository.

None of these data are redistributed here, as the benchmark record is access-restricted. The user registers with
the benchmark, downloads the record, and runs ERF's preparation script. Only the scripts and the decks are in this
repository. The stratified form follows M. P. van der Laan, M. Kelly, M. Baungaard, A. Dicholkar and E. L. Hodgson,
*Wind Energy Science* 9, 1985–2000 (2024), https://doi.org/10.5194/wes-9-1985-2024.

## How to run
From a working folder, with `W` this folder, `GUST` the executable and `ERF` the WINSENT ERF build:
```bash
# Tables (once per LAI): terrain, forest at 15 m, roughness map, levels, stations.
python3 $W/prepare_gust.py --erf-case <erf>/blended_high_neutral --out high --canopy-height 15 \
    --roughness <kynema-sgf>/auxiliary_files/winsent.roughness
# Neutral inflow: GUST's column at 1.2 times the power law (low LAI: add --powerlaw 0.423 0.479).
cd high && cp inflow.txt lidar.txt   # the lidar case mean, kept for plot_lidar.py
python3 $W/column_inflow.py --gust $GUST --zlevels zlevels.txt --lidar lidar.txt --out column_G12 --scale 1.2
cp column_G12/inflow.txt inflow.txt
mpiexec -n 8 $GUST $W/inputs_winsent_z0map steady.cfl_max=2 steady.max_iterations=3000
```
A stratified case takes its inflow from an ERF column (the class lidar mean gives the direction), with the values of
the table above, and runs with θ frozen:
```bash
python3 $W/lidar_mean.py <record>/lidar_data/high_LAI/stable --out lidar_mean.txt
python3 $W/erf_column.py --erf $ERF --out col --flux -0.004 --ug 6.69712 -1.42683 --lidar lidar_mean.txt
cp inflow.txt lidar.txt && cp col/inflow.txt inflow.txt   # unstable: add --lapse 0 to erf_column.py
mpiexec -n 8 $GUST $W/inputs_winsent_z0map steady.cfl_max=2 steady.max_iterations=3000 simplec.theta=frozen
```
Then, from the case folder, on the last plotfile (`winsent_plt` is the prefix; the highest step is read):
```bash
python3 $W/sample_stations.py winsent_plt --case . > stations.csv
python3 $W/plot_lidar.py stations.csv --lidar lidar.txt --out profiles.png   # lidar.txt: the class mean
python3 $W/make_submission.py winsent_plt --case . --template <record>/template_phase1_inflow.nc \
    --lai h --strat s --out-dir submission
```

## Expected Results
Iterations of the fine level to the default stop; the bias and RMS error of the wind speed at L140 over 100–400 m
against the lidar class mean; R_B at MMNW (46–73 m) against the class median. Runs on a Mac (Apple M-series): the
neutral cases two at a time with 4 ranks each (the low-LAI one alone on 8), the stratified ones one at a time on 8
ranks, so the wall times are indicative only.

| Case | Iterations | Wall time | L140 bias (m/s) | L140 RMS (m/s) | R_B at MMNW | Class median |
|---|---|---|---|---|---|---|
| High LAI, neutral | 360 | 19 min | 0.00 | 0.49 | – | – |
| Low LAI, neutral | 361 | 10 min | +0.19 | 0.32 | – | – |
| **High LAI, stable** | 398 | 9 min | +0.75 | 0.77 | +0.085 | +0.132 |
| Low LAI, stable | 443 | 10 min | −0.08 | 0.21 | +0.113 | +0.161 |
| High LAI, unstable | 338 | 9 min | +0.08 | 0.53 | −0.059 | −0.27 |
| Low LAI, unstable | 341 | 12 min | +0.13 | 0.39 | −0.134 | −0.35 |

High LAI, stable has a second option, not submitted: the column at 0.9 times the geostrophic wind with −0.005 K m/s
(6.027, −1.284 m/s; column R_B higher). Its residual oscillates between 1.8e-4 and 2.7e-4 up to 3000 iterations
(oscillatory convergence); L140 bias +0.54, RMS 0.94 m/s; R_B at MMNW +0.338 (run before the first-cell change
below, not rerun). The steady option (bold) is the
submitted one.

What the numbers show:
- In 3-D the MMNW R_B of the unstable cases is weaker than the columns': the θ difference (prescribed by height
  above the ground) stays the column's, so the ratio of the two R_B is the square of the shear ratio, 2.1 (high LAI)
  and 1.6 (low LAI) times the column's shear at MMNW. The columns have flat ground and a uniform surface; MMNW sits
  on the escarpment's edge, behind the forest.
- The first cell's stratification: these runs extend θ below the ground from the first two cells, so the first cell's
  closure sees the column's gradient. Earlier runs had the neutral wall's mirror there (N² = 0 in the first cell).
  The change moved each submitted value by at most 0.0013 m/s (speed), 0.08° (direction) and 0.0011 m²/s² (k): in
  the first cell the momentum eddy viscosity stayed the same (at MMNW c_μ fell from 0.556 to 0.442, but ν_t did not
  move, so a limit holds it there, not c_μ), and the heat diffusivity, which fell by 30 %, has no effect with θ
  frozen.
- A 0.9 G column for high LAI, stable, cycles instead of converging; lowering the eddy-viscosity relaxation
  (`turb.nut_relax` 0.6) did not remove the cycle.

**Canopy against roughness (PT-8.6.1).** High LAI, neutral, on 8 ranks with `steady.dtau_floor=1`: the chosen canopy, and the same
case with the forest as z0 = 1.5 m (0.1 of the 15 m trees) where LAI ≥ 2.5 and no canopy keys (`make_z0_only.py`).

| Run | Iterations | L140 bias (m/s) | L140 RMS (m/s) | L140 speed at 100 m | MMNW speed at 11 / 25 / 46 / 98 m |
|---|---|---|---|---|---|
| Canopy | 221, verdict | 0.00 | 0.49 | 3.58 m/s | 0.54 / 2.63 / 4.86 / 6.21 m/s |
| z0 = 1.5 m | 3000, **no verdict** | +0.01 | 0.41 | 3.75 m/s | 2.85 / 4.34 / 5.31 / 6.26 m/s |

The z0 run stays between r = 1.1e-4 and 2.4e-4 from about iteration 380 to the limit; its last field is the one scored.
At L140 the two differ by at most 0.17 m/s, so the lidar does not separate them, and the 15 m height and the inflow are fitted
to this lidar. They differ in the lowest 50 m at MMNW, inside and just above the forest, where this comparison has no measurement.

## CTest
None for the six cases (the data are not in the repository). The option they use, `simplec.theta = frozen`, is
tested on the forested escarpment (`rt_frozen_theta`, `Tests/Canonical_RANS/Forested_Escarpment`) and by the unit
tests `FrozenTheta.*`.
