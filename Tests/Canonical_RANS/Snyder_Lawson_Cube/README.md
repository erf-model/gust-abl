# Snyder_Lawson_Cube

## Purpose
The flow around an isolated cube (H = 0.2 m) in a wind-tunnel boundary layer, against the measurements of Snyder
and Lawson (1994) as Wise et al. (2026, Appendix A) use them to validate ERF's immersed forcing: the mean
streamwise velocity on the horizontal plane z = 0.02 m (0.1 H) and on the vertical centre plane, upwind of, over,
beside and behind the cube. The set-up follows Wise et al.: a 3.2 m x 1.2 m x 1.2 m domain (16 H x 6 H x 6 H), the
cube's upwind face 1.0 m from the inlet and the cube centred across it, inflow and outflow along x, slip walls on
the sides and the top, and the log law with z0 = 1 mm on the floor. GUST solves the steady RANS equations; the paper
time-averages LES. This is the cube validation of tracker row PR-44 (`docs/buildings.rst`, Validation).

The case also compares the closure's lee options (`turb.realizability`, the upwind blend) with and without their
zone around buildings (`buildings.limiter_zone`): a building's sharp edges fix where its flow separates, so the
options that let separation happen behind smooth terrain only shorten the mixing there.

## Files
| File | Role |
|---|---|
| `inputs_snyder_cube` | 160 x 60 x 90 cells: 20 mm across (10 across the cube), 5 mm layers to 0.3 m and layers growing by 1.1 above (`zlevels_5mm.txt`). The approach flow is the tunnel's log law, z0 = 1 mm and 4 m/s at z = 1 m (u* = 0.2374 m/s). The plain closure. |
| `inputs_snyder_cube_lee` | The same with the closure's lee options of the terrain decks: `turb.realizability` 0.309, `turb.nut_relax` 0.5 and the upwind blend (`convection.hybrid_sigma_max` 0.1, `convection.hybrid_height` 0.5 m, 2.5 H). `buildings.limiter_zone` (default 3) fades them out around the cube; `buildings.limiter_zone=0` on the command line keeps them there. |
| `zlevels_5mm.txt` | The 91 level heights of the decks. |
| `cube_snyder.txt` | The cube in ERF's format: x from 1.0 to 1.2 m, y from 0.5 to 0.7 m, 0.2 m high. |
| `snyder_lawson_u.csv` | The 296 measured values, digitized (Data). x from the cube's upwind face, y from the centreline. |
| `digitize_wise2026.py` | Reads the two figures and writes `snyder_lawson_u.csv` (Data). Not a test. |
| `compare_snyder.py` | Interpolates the run's u to the measured points; the bias, rms difference and hit rate per plane and region, u in the recirculation and the reattachment length; `--pin` for the tests. The exit code is the verdict. |
| `plot_snyder.py` | The profiles in the layout of Wise et al.'s Fig. A3 (`docs/figures/buildings/snyder_lawson_profiles.png`). Not a test. |

## Data and acknowledgements
The measurements are those of **Snyder, W. H. and Lawson, R. E. (1994)**, Wind-tunnel measurements of flow fields in
the vicinity of buildings, *8th Joint Conference on Applications of Air Pollution Meteorology*, American
Meteorological Society, made in the U.S. EPA Fluid Modeling Facility wind tunnel. We thank William Snyder and Robert
Lawson, and the EPA Fluid Modeling Facility, for the measurements. The conference paper has no DOI; we know of no
public release or licence of the measured values themselves, made by the U.S. Environmental Protection Agency.

The values in `snyder_lawson_u.csv` are digitized from Figs. A3 and A4 of **Wise, A. S., Gopalan, H., Lattanzi,
A. M., Almgren, A. S., Wiersema, D., Gowardhan, A., Kang, S., Mirocha, J. D., Gardner, D. J., and Lundquist, K. A.
(2026)**, Implementation of a resolved building capability in the Energy Research and Forecasting model for
accelerated high-fidelity urban-scale simulation, *EGUsphere preprint*, https://doi.org/10.5194/egusphere-2026-4030,
licensed CC BY 4.0. We thank the authors; the figures' black dots are the measurements above. The paper's Zenodo
record (https://doi.org/10.5281/zenodo.21385466, CC BY 4.0) holds its input files and time-averaged fields for the
Oklahoma City cases, not these measurements.

`digitize_wise2026.py` reads the raster images of the preprint's pages 37 and 38 (extracted with
`pdfimages -png -f 37 -l 38 egusphere-2026-4030.pdf fig`). The figures' scales come from their tick marks and arrows:
1400 px per m along x in both panels, each profile drawn from its station as x_s + 0.025 u. Panel (a)'s lateral
labels read twice y: the panel has equal aspect, and the cube (0.2 m wide) and the refined region stated as
y in (-0.22, 0.22) m are drawn over twice those spans, so y is the label halved. The values are Fig. A3's, with
three from Fig. A4 that A3's legend hides and one under the 3 m/s arrow's head recovered from the dot's right edge.
The two figures agree to 0.003 m/s rms; a pixel is 0.029 m/s. The measurements are mirror-symmetric about the
centre plane to the pixel, so the lateral profiles likely hold one measured side mirrored. Redistributed here:
the digitized values (from a CC BY 4.0 figure, with attribution). Not redistributed: the preprint's figures.

## How to run
From this folder (the z levels and the building map go on the command line):
```bash
mpiexec -n 4 <build>/gust inputs_snyder_cube buildings.file=cube_snyder.txt mesh.zlevels_file=zlevels_5mm.txt
```
```bash
python3 compare_snyder.py --latest snyder_plt
```
```bash
python3 plot_snyder.py profiles.png "snyder_plt:GUST"
```
The lee options with their zone, and without it (the form before it):
```bash
mpiexec -n 4 <build>/gust inputs_snyder_cube_lee buildings.file=cube_snyder.txt mesh.zlevels_file=zlevels_5mm.txt
```
```bash
mpiexec -n 4 <build>/gust inputs_snyder_cube_lee buildings.file=cube_snyder.txt mesh.zlevels_file=zlevels_5mm.txt buildings.limiter_zone=0
```
The 10 mm run: `mesh.n_cell = 320 120 90` and `mesh.max_grid_size = 40 40` in a copy of `inputs_snyder_cube`.

## Expected Results
The run stops at the driver's stopping rule. Measured against the 296 values (`compare_snyder.py`): the rms
difference of u, its bias, and the hit rate q (the share within 25 % or 0.1 m/s of the measurement, VDI 3783 Part 9's
form; q moves by about 0.01 for D between 0.03 and 0.15 m/s), in all and by region; and the reattachment length on
the centreline at z = 0.02 m. The measured one is 1.36 H behind the cube on the horizontal plane (linear between the
stations at x = 0.4 and 0.5 m) and 1.33 H on the vertical one: about +-0.03 H. The lateral profiles mirror each
other to 0.002 m/s and 7 points are in both planes, so the 296 values hold about 180 independent ones. These are a
validation against measurements, not pinned to a target beyond the tests.

20 mm cells across (10 across the cube), Mac, 4 or 8 ranks:

| Layers | Closure | Iterations | rms (m/s) | Hit rate | Upwind rms | Over/beside rms | Wake rms | u at x = 0.3 m, z = 0.02 m | Reattachment |
|---|---|---|---|---|---|---|---|---|---|
| uniform 20 mm (60) | plain | 534 | 0.256 | 0.882 | 0.069 | 0.612 | 0.163 | -1.15 m/s | 1.32 H |
| **5 mm to 0.3 m (90)** | **plain** | **472** | **0.244** | **0.882** | 0.075 | 0.566 | 0.167 | -1.15 m/s | **1.34 H** |
| 5 mm to 0.3 m | + `turb.realizability` 0.309, `turb.nut_relax` 0.5 (`buildings.limiter_zone=0`) | 466 | 0.319 | 0.716 | 0.170 | 0.474 | 0.348 | -1.12 m/s | 1.88 H |
| 5 mm to 0.3 m | + the upwind blend (`inputs_snyder_cube_lee`, `buildings.limiter_zone=0`) | 466 | 0.317 | 0.713 | 0.159 | 0.484 | 0.345 | -1.11 m/s | 1.89 H |
| **5 mm to 0.3 m** | **the same with the zone, `buildings.limiter_zone` 3 (`inputs_snyder_cube_lee`)** | **471** | **0.244** | **0.882** | 0.075 | 0.567 | 0.167 | -1.15 m/s | **1.35 H** |
| measured | | | | | | | | -0.80 m/s | 1.36 H |

- The thin layers mainly help over and beside the cube (0.566 against 0.612 m/s), for 1.5 times the cells and fewer
  iterations. Near the floor the two read alike once the comparison interpolates the first cells in log height.
- Over the roof the plain closure has no bubble at 0.1 H above it: at x = 0.06 m, z = 0.22 m u is 1.59 m/s against
  0.03 m/s measured (1.51 against 0.74 at x = 0.14 m), and it is slow in the speed-up above (3.65 against 3.88 m/s at
  z = 0.28 m). The roof separates at the front edge in both, but the plain closure's mixing reattaches it at once.
- The eddy-viscosity cap of the terrain decks acts differently by region. Over the roof it gives the bubble the
  measurements show (u -0.15 m/s at x = 0.06 m, z = 0.22 m) and the over/beside rms falls to 0.474 m/s; at
  x = 0.3 m it brings the near-wake bubble's height closer (u at z = 0.16 m 0.07 against 0.14 measured; 0.62 plain).
  Elsewhere it is worse: the shear layer over the near wake is too slow (0.81 against 1.77 m/s at x = 0.3 m,
  z = 0.2 m), the bubble on the floor reaches 1.88 H, the far wake lags (by 0.05 to 0.23 m/s at x = 1.0 to 2.0 m),
  and the horseshoe vortex ahead of the cube grows (u lower by up to 0.35 m/s at x = -0.1 m, z = 0.02 m). On the
  approach flow away from the cube it changes u by at most 0.01 m/s: at C_r = 0.309, just below Cmu0^2 = 0.30936,
  it lowers the eddy viscosity of an equilibrium layer by 0.12 %. In all, rms 0.319 against 0.244 m/s; the blend
  adds little. The cube's edges fix where its flow separates, so the case for the cap that a smooth hill makes does
  not hold here, and the overall rms is what sets the zone.
- `buildings.limiter_zone` lifts the cap within 1.5 H of the cube and fades it back by 3 H (the log: off in 96180
  cells, partly in 203368 of 864000). Every measured point lies within 3 H, 212 of them where the cap is lifted;
  the answer there is the plain closure's within 0.001 m/s and 0.001 m. Beyond 3 H the cap acts again: u 1.5 m
  behind the cube at z = 0.2 m is 2.802 m/s, between the plain closure's 2.838 and the cap's 2.729. c = 3 is about
  the smallest size that keeps the cap off over the whole measured bubble (2.7). A wider building's wake is longer,
  in proportion to its width more than its height, and can reach beyond 3 H.
- What none of these settings gets: the speed-up past the side edges at x = 0.04 and 0.16 m (|y| = 0.12 to 0.2 m),
  with 10 cells across the cube, and the strength of the reversed flow near the floor behind the cube (-1.15 m/s
  against -0.80 m/s). The reattachment point is right at 20 mm (1.34 H against 1.33 to 1.36 H).
- Part of the bias near the floor upwind is the inflow. The tunnel's approach flow is described both as a log law
  with z0 = 1 mm and as a power law with exponent 0.16; the measured far-lateral values at z = 0.02 m (1.86 to 1.91
  m/s) lie between the log law's 1.76 and the power law's 2.14, and GUST, on the log law, is 0.06 to 0.15 m/s low
  there.
- One measured point disagrees between the two panels of the figures: x = -0.1 m, y = 0, z = 0.02 m reads
  0.34 m/s on the horizontal plane and 0.11 m/s on the vertical one, in both Fig. A3 and Fig. A4 (an 8 px
  difference); the other shared points agree to 0.03 m/s. Both values are kept.

10 mm cells across (20 across the cube), 5 mm layers to 0.3 m (320 x 120 x 90 cells, 3.5 M), plain closure, 12 ranks:

| Cells across | Iterations | rms (m/s) | Hit rate | Upwind rms | Over/beside rms | Wake rms | u at x = 0.3 m, z = 0.02 m | Reattachment |
|---|---|---|---|---|---|---|---|---|
| 20 mm | 472 | 0.244 | 0.882 | 0.075 | 0.566 | 0.167 | -1.15 m/s | 1.34 H |
| **10 mm** | **1728** | **0.216** | **0.905** | 0.073 | **0.475** | 0.169 | -1.14 m/s | 1.29 H |
| measured | | | | | | | -0.80 m/s | 1.33 to 1.36 H |

- Halving the cells across the cube lowers the error over and beside it (0.475 against 0.566 m/s) and raises the hit
  rate on the horizontal plane to 0.952; the vertical plane's rms is 0.259 m/s, most of it over the roof.
- The reattachment moves from 1.34 to 1.29 H, so it is not grid-converged; the measured value lies between the two.
- The run took 1728 iterations, 4.0 h on 12 cores, in 150-iteration segments restarted from their checkpoints. The k
  residual was the last to pass `steady.residual_tol`; u at the measured points changed by less than 1e-4 m/s over
  its last 800 iterations.

Wise et al. (2026) report, for ERF's LES with the TKE-1.5 closure and 5 mm cells, the best agreement in their
figures; they give no rms or hit rate, so the two are compared in the figure only.

## CTest
`snyder_cube_plain_clean`, `snyder_cube_plain_run`, `snyder_cube_plain_compare`, `snyder_cube_lee_clean`,
`snyder_cube_lee_run`, `snyder_cube_lee_compare`, `snyder_cube_lee_nozone_clean`, `snyder_cube_lee_nozone_run`,
`snyder_cube_lee_nozone_compare` (minutes each on 4 ranks; label `regression`, not `smoke`); `abort_buildings_limiter_zone`,
`abort_buildings_limiter_zone_wall_tke`, `abort_buildings_limiter_zone_no_cap`.
