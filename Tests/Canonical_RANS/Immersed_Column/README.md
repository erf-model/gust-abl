# Immersed_Column

## Purpose
The 1-D neutral column of [Neutral_Column](../Neutral_Column/README.md) with the ground immersed
(`terrain.representation = immersed_forcing`, [docs](../../../docs/immersed_terrain.rst)). The mesh is flat:
4 x 4 x 67 layers of dz = 15.625 m (1046.875 m), doubly periodic, two boxes in x. The ground lies at
`terrain.flat_height` = (2 + f) dz, inside layer k = 2, so that f is the cut cell's solid fraction.
z0 = 0.1 m, slip top, uniform theta. The force F = 1.225e-4 m/s² acts on the fluid only.

At steady the immersed wall stress balances the force on the fluid exactly:

    u*² = F D,   D = 1046.875 m - flat_height.

u*² is recomputed from the plotfile with the wall law at the exact distances (ERF form 7):

    u*² = (1 - f) [κ u_c / ln((d1 + z0)/z0)]² + f [κ u_a / ln((d2 + z0)/z0)]²,
    d1 = (1 - f) dz / 2,   d2 = (1 - f) dz + dz / 2,

where u_c is the cut cell's speed, u_a the speed of the cell above it, and κ = 0.41. With f = 0 the ground lies on a
cell face, and the column must repeat the body-fitted column of the fluid's height (65 layers, 1015.625 m).

## Files
| File | Role |
|---|---|
| `inputs_immersed_column` | The immersed column. `terrain.flat_height` is given on the command line. |
| `inputs_body_fitted_column` | The body-fitted column of 65 layers: the reference for f = 0. |
| `check_immersed_column.py` | Runs the application at a given f and checks the items below. |
| `inputs_immersed_ground` | One iteration of the column: writes the checkpoint the restart tests read. |
| `inputs_immersed_restart` | Restart from that checkpoint on another box layout. |
| `inputs_restart_body_fitted` | Restart from that checkpoint with a body-fitted ground (refused). |
| `inputs_abort_log_shear_wall` | The body-fitted column with a no-slip ground and `turb.shear_stencil = log_law` (refused). |

Both decks set `solver.poisson_backend = gmres_fft`. With the default `gmres_mlmg`, f = 0.3, 0.7 and 0.9 do not
reach the steady verdict. In those runs the divergence is round-off but not zero, and the multigrid-preconditioned
Poisson solve returns a pressure that differs between columns by round-off. v and w then grow to about 1e-9, while
the stop rule needs them below 1e-14 here (they start at zero, so their reference is the floor `steady.residual_tol`).
The mesh of an immersed ground is always flat, so the FFT-preconditioned solve applies. With it v stays exactly zero.

## Data
None: the case is analytic.

## How to run
```bash
ctest --test-dir <build> -R "^rt_immersed_column"
```
By hand, for f = 0.7:
```bash
python3 check_immersed_column.py --f 0.7 --cut-max 2.5 -- <build>/gust inputs_immersed_column
```
With the log-consistent shear (`turb.shear_stencil = log_law`, docs/turbulence_model.rst):
```bash
ctest --test-dir <build> -R "^rt_log_shear_column"
```
```bash
python3 check_immersed_column.py --f 0.7 --cut-max 0.2 --log-shear -- <build>/gust inputs_immersed_column turb.shear_stencil=log_law
```

## Expected Results
- The steady verdict, every field finite, every column identical (|u − planar mean| < 1e-10 u*).
- u*² of the exact law equals F D to 1e-5 (measured within 5.3e-7: the residual tolerance).
- The solid layers k = 0, 1 have u = 0 exactly.
- The cell above the cut cell within 1 % of the log law with u* = sqrt(F D), the next two within 2.5 %.
- The cut cell at or above the log law at its fluid centroid d1, by at most the margin in the table. The form-7 law
  takes the cut cell's speed as the speed at d1. The excess grows with f.
- f = 0: the body-fitted column's u and k in every layer equal the immersed column's layers k ≥ 2 to 1e-6.

Measured (two ranks, double precision):

| f | Iterations | u*²/(F D) − 1 | Cut cell vs log law at d1 (margin) | Layer 3 | Layer 4 | Layer 5 |
|---|---|---|---|---|---|---|
| 0 | 2266 | −5.3e-7 | −0.00 % (0.1 %) | −0.89 % | −1.33 % | −1.51 % |
| 0.3 | 2284 | −4.1e-7 | +0.36 % (1 %) | −0.85 % | −1.44 % | −1.66 % |
| 0.7 | 2363 | −2.0e-7 | +1.77 % (2.5 %) | −0.77 % | −1.75 % | −2.03 % |
| 0.9 | 2509 | −6.9e-8 | +5.64 % (6.5 %) | −0.65 % | −2.01 % | −2.35 % |

f = 0 against the body-fitted column: u within 4.2e-9 and k within 1.0e-8 in every layer.

The cut cell's excess is the cell above's deficit times f / (1 − f): u*² is fixed at F D, so
(1 − f)(1 + e1)² + f (1 + e2)² = 1.

### With `turb.shear_stencil = log_law`
The closure's central difference of the log profile over-estimates du/dz in the cells next to the wall cell, and k
there is 13 to 22 % above the local balance u*² (1 − d/D) / Cmu0² (shear production = dissipation at the stress
u*² (1 − d/D)) in the first cell above the cut cell, 7 to 14 % in the second and 4 to 7 % in the third. With the
option the checks are (`--log-shear`):
- the steady verdict, the force balance to 1e-5, the solid layers at rest and identical columns, as above;
- the cut cell within 0.2 % of the log law at d1, on either side;
- the three cells above the cut cell within 0.5 % of the log law;
- their k within 1 %, 1.5 % and 2 % of the local balance;
- f = 0: the body-fitted column with the option equals the immersed column to 1e-6.

Measured (two ranks, double precision):

| f | Iterations | u*²/(F D) − 1 | Cut cell vs log law at d1 | Layers 3, 4, 5 vs log law | k of layers 3, 4, 5 vs local balance |
|---|---|---|---|---|---|
| 0 | 2313 | −5.5e-7 | −0.00 % | −0.14, −0.22, −0.30 % | −0.54, −1.11, −1.56 % |
| 0.3 | 2339 | −4.3e-7 | +0.03 % | −0.07, −0.16, −0.23 % | −0.45, −1.00, −1.46 % |
| 0.7 | 2440 | −2.1e-7 | +0.01 % | −0.00, −0.10, −0.19 % | −0.17, −0.80, −1.29 % |
| 0.9 | 2606 | −7.2e-8 | −0.12 % | +0.01, −0.10, −0.18 % | +0.07, −0.69, −1.20 % |

These decks have no length cap (`turb.max_geom_lscale = 1e12`). With a capped length the cut cell's speed falls
below the log law with the option (docs/turbulence_model.rst, Limits).

Not covered: f = 0.98 and 0.99 do not reach the steady verdict in 6000 iterations with these decks, with either
Poisson backend (f = 0.95 does, in 2612). The w residual rises to about 2e-8 at iteration 40 and the step-size rule
then keeps the CFL number low. f = 0.99 reaches the verdict with `simplec.pressure_reltol = 1e-8` (2863 iterations) or
`steady.cfl_max = 10` (3175); the cause is not identified. Its layers 3-5 are within −0.5 to −2.6 % of the log law
and its cut cell 43 % above it at d1 (the factor 99). `terrain.small_volfrac` (default 0.005) treats only
fractions above 0.995 as solid, so such cells occur.

## CTest
`rt_immersed_column_f0` (with the body-fitted comparison), `rt_immersed_column_f03`, `rt_immersed_column_f07`,
`rt_immersed_column_f09`. Double precision only (six orders and 1e-8 are below single's floor). About 30 s each on
two ranks; f = 0 runs twice.

With `turb.shear_stencil = log_law`: `rt_log_shear_column_f0` (with the body-fitted comparison),
`rt_log_shear_column_f03`, `rt_log_shear_column_f07`, `rt_log_shear_column_f09`, the same cost.
`rt_log_shear_column_f09` solves the pressure to `simplec.pressure_reltol = 1e-6` (2703 iterations): from the
first guess, with the cut layer's centroid height, the default 1e-3 leaves a uniform w near 1e-10 after each
solve, and the u and w residuals level off near 2e-8, above `steady.residual_tol = 1e-8`. The option's
refusals (smoke): `abort_shear_stencil_value`, `abort_shear_stencil_buildings`, `abort_shear_stencil_tucker` and
`abort_shear_stencil_no_surface_layer` stop naming the input.

The ground's wiring in the application (smoke, a few seconds each): `immersed_ground_clean`, `immersed_ground_run`
(2 ranks) and `immersed_ground_restart` (3 ranks, another box layout: the checkpoint's ground is accepted);
`abort_immersed_restart_changed` (another `terrain.flat_height`) and `abort_immersed_restart_body_fitted` stop
naming `io.restart`; `abort_immersed_first_layer`, `abort_immersed_buildings` and `abort_immersed_surface_zref` stop
naming the input.
