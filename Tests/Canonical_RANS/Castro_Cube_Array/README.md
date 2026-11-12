# Castro_Cube_Array

## Purpose
The flow over a staggered array of cubes (h = 20 mm, 25 % plan area) in a wind-tunnel boundary layer, against the
measurements of Castro et al. (2006) as Muñoz-Esparza et al. (2020, section 3) use them to validate FastEddy's
immersed body force: profiles of the mean streamwise velocity u above a cube (P0), one h behind it (P1), one h
beside it (P2) and one h behind and beside it (P3), and of the momentum flux -<u'w'> at P1 and P3. The set-up
follows Muñoz-Esparza et al.: one 4h x 4h tile of the array, periodic in x and y, 10h deep, driven by a constant
pressure gradient. GUST solves the steady RANS equations; the paper time-averages LES. This is the array validation
of tracker row PR-44 (`docs/buildings.rst`, Validation).

## Files
| File | Role |
|---|---|
| `inputs_cube_array` | 40 x 40 x 59 cells: 2 mm across (10 across a cube), 1 mm layers to 1.5 h and layers growing by 1.1 above (`zlevels_1mm.txt`). Pressure driven with a slip lid, F = u_τ²/L_z = 0.7²/0.2 = 2.45 m/s² (the tunnel's u_τ, no tuning). The first guess is the log law with z0 = 0.1 h (Expected Results). 3000 iterations. |
| `zlevels_1mm.txt` | The 60 level heights of the deck: 1 mm to 0.03 m, then growing by 1.1; the top layer is what is left to 0.2 m (22.4 mm, 1.55 times the one below), under the slip lid. |
| `array_tile.xyz` | The tile's four cubes as 2 mm columns (`buildings.format = xyz`): centred at (0.02, 0.02), (0.02, 0.06), (0.06, 0.04) and (0.06, 0) m, the last cut by the periodic seam. |
| `castro2006_fig3.csv` | The 185 measured values, digitized (Data). u is u/u_r with u_r = 10 m/s; uw is -<u'w'>/u_r² in units of 1e-3. |
| `digitize_munoz2020.py` | Reads the figure's vector paths and writes `castro2006_fig3.csv` (Data). Not a test. |
| `compare_array.py` | Interpolates the run's u and the closure's -<u'w'> = ν_t (∂u/∂z + ∂w/∂x) to the measured points; the bias, rms difference and hit rate per station, the rms after scaling u to the measured u at P0's top point, and the stress peaks near the roof; `--pin` for the tests. The exit code is the verdict. |
| `plot_array.py` | The profiles in the layout of Muñoz-Esparza et al.'s Fig. 3. Not a test. |

## Data and acknowledgements
The measurements are those of **Castro, I. P., Cheng, H., and Reynolds, R. (2006)**, Turbulence over urban-type
roughness: deductions from wind-tunnel measurements, *Boundary-Layer Meteorology*, 118, 109–131,
https://doi.org/10.1007/s10546-005-5747-7, made in the EnFlo wind tunnel of the University of Surrey. We thank Ian
Castro, Hong Cheng and Robert Reynolds for the measurements. We know of no public release or licence of the
measured values themselves.

The values in `castro2006_fig3.csv` are digitized from Fig. 3 of **Muñoz-Esparza, D., Sauer, J. A., Shin, H. H.,
Sharman, R., Kosović, B., Meech, S., García-Sánchez, C., Steiner, M., Knievel, J., Pinto, J., and Swerdlin, S.
(2020)**, Inclusion of building-resolving capabilities into the FastEddy® GPU-LES model using an immersed body force
method, *Journal of Advances in Modeling Earth Systems*, 12, e2020MS002141, https://doi.org/10.1029/2020MS002141,
licensed CC BY-NC-ND 4.0. We thank the authors; the figure's symbols are the measurements above. The paper's
Zenodo record (https://doi.org/10.5281/zenodo.4124164, CC BY 4.0) holds its processed model output for the Oklahoma
City cases (`U_WD_TKE_comparison.nc`, `SF6_comparison.nc` and a video), not these measurements.

`digitize_munoz2020.py` reads page 7 of the paper converted to SVG (`pdftocairo -svg -f 7 -l 7 paper.pdf p7.svg`).
The figure is vector graphics: each symbol is a filled path, found by its colour, and the scales come from the
axis frames. A symbol partly covered by its neighbour is placed from its free edge; two symbols cut by panel (c)'s
legend box (P1 and P3 near z = 2.15 h) from their left edges; one fragment too small to place is dropped. One point of the
SVG is 0.008 u/u_r, 0.1 x 1e-3 u_r² and 0.016 h; a symbol's centre is read to about half a point.

Redistributed here: the 185 digitized values, with attribution. They are the measured values of Castro et al.
(2006) as the figure plots them, read from it; the file holds no part of the figure itself. The paper is licensed
CC BY-NC-ND 4.0, which covers its text and figures; its Zenodo record is CC BY 4.0 but does not hold these values.
Not redistributed: the paper's figures and text. If the authors or the measurements' owners publish the values
under a licence, the file should be replaced by theirs.

## How to run
From this folder (the z levels and the building map go on the command line):
```bash
mpiexec -n 4 <build>/gust inputs_cube_array buildings.file=array_tile.xyz mesh.zlevels_file=zlevels_1mm.txt io.turb_plot_file=turb
```
```bash
python3 compare_array.py --latest array_plt turb castro2006_fig3.csv
```
```bash
python3 plot_array.py castro2006_fig3.csv profiles.png "2mm:array_plt03000:turb03000"
```

## Expected Results
The run stops after 3000 iterations without meeting the steady criteria, and more iterations do not meet them:
from 3000 to 8000 the u residual stays at 3.0 to 3.2e-4, while the k residual still falls (about 13 times in 5000
iterations). Below 2.2 h, the highest measured point, the answer has settled: between 3000 and 8000 iterations u
at the measured points changes by at most 0.001 u_r and the stress peaks by 0.07 x 1e-3 u_r². The outer layer is
still evolving: the plane-mean k falls by 13 % at 4.2 h and by 43 % at 9.4 h, and the scale to P0's top point
moves from 1.054 to 1.049. Measured against the 185 values (`compare_array.py`): the rms
difference and bias of u/u_r and the hit rate q (the share within 25 % or 0.01 u_r of the measurement: VDI 3783 Part
9's form, with D = 0.01 u_r our choice) per station and in all; the rms of -<u'w'>/u_r²; and its peak near the roof (0.8 h < z < 1.3 h) at P1
and P3. These are a validation against measurements, not pinned to a target beyond the test.

2 mm cells across (10 across a cube), Mac, 4 ranks:

| First guess | Iterations | Time | u rms | u bias | q | Scale to P0 top | Scaled rms | -<u'w'> rms | Peak P1 | Peak P3 |
|---|---|---|---|---|---|---|---|---|---|---|
| log law, z0 = 1e-4 m (the floor's) | 6000 | 27 min | 0.154 | +0.115 | 0.220 | 0.699 | 0.035 | 7.11 | 23.3 | 13.6 |
| **log law, z0 = 0.1 h** | **3000** | **12 min** | **0.042** | **-0.032** | **0.724** | 1.054 | 0.033 | 1.71 | 10.45 | 5.97 |
| log law, z0 = 0.1 h | 8000 | 45 min | 0.041 | -0.031 | 0.724 | 1.049 | 0.033 | 1.72 | 10.52 | 6.01 |
| measured | | | | | | | | | 8.49 | 5.80 |

u in u_r, -<u'w'> in 1e-3 u_r².

The grid, rough first guess, 4 ranks (4000 iterations at 1 mm, 80 x 80 x 59 cells, 20 across a cube; the same
layers):

| Cells across | Iterations | r at the end | u rms | u bias | q | Scale to P0 top | -<u'w'> rms | Peak P1 | Peak P3 |
|---|---|---|---|---|---|---|---|---|---|
| **2 mm** | **3000** | 5.7e-3 | **0.042** | -0.032 | **0.724** | 1.054 | 1.71 | 10.45 | 5.97 |
| 1 mm | 4000 | 9.1e-3 | 0.045 | -0.036 | 0.732 | 1.072 | 1.68 | 10.30 | 5.73 |
| measured | | | | | | | | 8.49 | 5.80 |

- Halving the cells changes u at the measured points by 0.003 u_r rms and the stress peaks by 1 to 4 %: at these
  points the 2 mm answer is grid-independent within the scatter of the comparison. The 1 mm run is less converged
  (its u residual is 1.6e-3 and still falling, against 3.2e-4 at 2 mm), and its free stream is 2 % slower.
- 2 ranks give the 4-rank pins to the digits printed.

- The first guess matters. From the smooth floor's log law the flow starts too fast over the array and slows
  over thousands of iterations: after 6000 the plane-mean u at z = 1.3 h is still 50 % above that of the rough
  start, and the stresses are more than twice the measured ones. From z0 = 0.1 h (an array-like roughness) it settles
  near the measurements: the u rms is 0.048 after 1000 iterations and 0.042 after 3000.
- By station (3000 iterations): P0 over the cube q 1.00 (rms 0.029); P1 behind the cube q 0.55 (rms 0.033), where
  u is small and the 25 % band narrow; P2 beside the cube rms 0.056, u too slow by 0.054 on average; P3 rms 0.031.
- The momentum flux has the measured shape. Its peak just above the roof is 23 % high at P1 and 3 % high at P3.
- With F set from the tunnel's u_τ, u at P0's top point (2.2 h) comes out 5 % low (the scale is 1.054).
  Muñoz-Esparza et al. adjusted their pressure gradient to match that point. Scaled the same way, GUST's rms is
  0.033 and q 0.764.
- Their force is 36 % larger (3.875 Pa/m, 3.34 m/s²), but their no-slip lid carries part of it: a smooth-wall log
  law at 10 m/s near the lid takes 0.18 to 0.19 m²/s², which leaves 0.47 to 0.49 m²/s² on the floor and the array,
  u_τ 0.69 to 0.70 m/s, as in GUST (0.49 m²/s²). The two floor stresses agree, so the 5 % comes from GUST, not
  from the forcing. With the slip lid the total stress at the roof level is F (L_z - h), u_* = 0.664 m/s there; if
  the tunnel's 0.7 m/s is a roof-level value, F = 0.49 / (L_z - h) = 2.72 m/s² would match it. With F = 2.72
  (3000 iterations) the u rms is 0.039, q 0.732, the scale 1.042 and the stress peaks 10.70 and 6.10: the
  free-stream deficit shrinks from 5 % to 4 % and the peaks rise. The deck keeps FastEddy's F = u_τ²/L_z.
- FastEddy's set-up has no-slip walls at the floor and the top, a 1 mm grid and the viscous terms. GUST's
  pressure-driven set-up needs a slip lid; the force F = u_τ²/L_z is that of a half channel of depth L_z.
- With buildings in a doubly periodic domain the pressure correction is singular and the unmasked multigrid
  preconditioner of `gmres_mlmg` needs about 300 GMRES iterations per solve; the deck uses `gmres_fft`, which gives
  the same answer about 7 times faster.

## CTest
`castro_cube_array_clean`, `castro_cube_array_run`, `castro_cube_array_compare` (about 12 minutes on 4 ranks; label
`regression`, not `smoke`). The run test expects "gust: not steady after 3000 iterations"; the compare test pins
the u rms (0.042), q (0.724), the -<u'w'> rms (1.707) and the peaks at P1 (10.45) and P3 (5.97), each to 2 %.
