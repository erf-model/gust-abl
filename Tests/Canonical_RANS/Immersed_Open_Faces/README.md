# Immersed_Open_Faces

## Purpose
Flat ground immersed in a flat mesh (`terrain.representation = immersed_forcing`,
[docs](../../../docs/immersed_terrain.rst)) between an inflow face and an outflow face, against the same flow on a
body-fitted mesh. On flat ground the log-law inflow is the solution, so the run should carry it to the outflow
face unchanged, whatever fraction F of its layer the ground fills.

The channel is 2560 m long and 1000 m high above the ground, neutral (uniform θ, no gravity), with z0 = 0.1 m,
`bc.xlo = inflow` (`bc.inflow_profile = log_law`, u* = 0.35 m/s), `bc.xhi = outflow`, periodic in y and a slip
top. The body-fitted mesh has 64 x 4 x 64 cells with layers of 15.625 m. The immersed mesh has the same layers and
three more below the ground: the ground lies in the third layer, at solid fraction F of it.

Three parts of the code are exercised that the periodic cases (`../Immersed_Column`, `../Immersed_Ridge_2D`) do
not reach:

- the inflow face: solid faces carry no flow, and the profile's height is measured from the immersed ground (a
  point in the cut layer takes the centroid of its open part);
- the outflow face and the mass balance between the two faces;
- the pressure solve with solid cells at open faces.

## Files
| File | Role |
|---|---|
| `inputs_body_fitted_flat` | The body-fitted channel. |
| `inputs_immersed_flat` | The immersed channel. `terrain.flat_height` is set by the scripts. |
| `sweep_flat.py` | Runs the body-fitted deck and the immersed deck at each F, and compares the outflow columns. With `--limits` it is the test. |
| `sweep_band.py` | The same on layers with a fine band, the ground placed in a coarse layer, in the transition, in the band and at its top. Not a CTest. |

## Data
None: the case is analytic.

## How to run
```bash
ctest --test-dir <build> -R "^rt_immersed_open_faces"
```
By hand:
```bash
python3 sweep_flat.py <build>/gust <folder> [--limits DU DK] [--fresh] [key=value ...] [F ...]
python3 sweep_band.py <build>/gust <folder> [key=value ...]
```
`key=value` arguments replace or add deck keys in every run (each run folder holds its deck as `deck`; `key=` with no value removes a key);
`terrain.*` keys go to the immersed runs only. With no F the sweep is 0 to 0.9 in steps of 0.1, 0.99 and 1. A run already in the folder is read, not repeated,
unless `--fresh` is given.

## Expected Results
Every run reaches the steady verdict, and the mass flux through the outflow face equals that through the inflow
face to the printed precision (the test's limit is 1e-6 of the flux). `du` and `dk` are the largest differences of
u and of k from the body-fitted outflow column at equal height above the ground, in percent of the body-fitted
value. The body-fitted column is interpolated in ln(z + z0), u as its ratio to the log law, which is exact between
two layers for a logarithmic profile. `cut` is the cut cell's speed over the log law at its centroid.

Default settings (`turb.shear_stencil = central`, `terrain.small_volfrac = 0.005`); the body-fitted run takes 122
iterations:

| F | iterations | du (%) | dk (%) |
|---|---|---|---|
| 0 | 122 | 0.10 | 0.67 |
| 0.1 | 122 | 0.12 | 1.37 |
| 0.3 | 122 | 0.19 | 3.35 |
| 0.5 | 123 | 0.33 | 4.75 |
| 0.7 | 124 | 0.55 | 3.97 |
| 0.9 | 128 | 1.25 | 6.37 |
| 0.99 | 139 | 3.45 | 25.4 |
| 1 | 122 | 0.00 | 0.02 |

At F = 0.99 the cut cell is 0.16 m thick and its speed is 3.8 times the log law's: the wall law's weight on the
cut cell grows as F/(1 - F) (docs, "Wall law"). `terrain.small_volfrac` merges such a cell into the ground.

The setting of the real-terrain cases (`turb.shear_stencil=log_law terrain.small_volfrac=0.1`), which is the
CTest:

| F | iterations | du (%) | dk (%) |
|---|---|---|---|
| 0 | 122 | 0.10 | 0.67 |
| 0.5 | 127 | 0.81 | 5.52 |
| 0.95 (merged into the ground) | 122 | 0.13 | 0.55 |
| 1 | 122 | 0.00 | 0.03 |

`sweep_band.py` with `turb.shear_stencil=log_law` (and the default `terrain.small_volfrac`), F = 0 to 0.9 in
steps of 0.1 and 0.99; the body-fitted runs take 120 to 124 iterations:

| Ground in | Layer (m) | Iterations, F ≤ 0.9 | du (%), F ≤ 0.9 | dk (%), F ≤ 0.9 | At F = 0.99: iterations, du, dk |
|---|---|---|---|---|---|
| uniform layers | 15.6 | 122 to 127 | up to 2.3 | up to 23 | 153, 4.7, 41 |
| a coarse layer below the band | 16 | 123 to 136 | up to 2.5 | up to 24 | 211, 5.0, 41 |
| the transition to the band | 12 | 123 to 141 | up to 1.1 | up to 14 | 212, 1.2, 8.0 |
| the band | 4 | 120 to 138 | up to 0.7 | up to 23 | 127, 4.7, 30 |
| the band's last layer | 4 | 119 to 124 | up to 1.1 | up to 23 | 135, 5.3, 30 |

The flux through the two faces is equal in all 60 runs. dk grows with F in every zone: the thin cut cell moves
faster than the log law, and the k of the cell above it falls (docs, "Limits"). The real-terrain cases set
`terrain.small_volfrac = 0.1` for this reason.

Open-face options, second setting, each against the body-fitted channel with the same option (the option CTests):

| Option | Body-fitted iterations | F = 0.5: iterations, du, dk (%) | F = 0.95: iterations, du, dk (%) |
|---|---|---|---|
| sponge, 320 m, 0.05 | 121 | 118, 0.28, 2.68 | 121, 0.07, 2.05 |
| power-law inflow | 125 | 131, 0.83, 5.20 | 125, 0.13, 0.54 |
| side walls moving at 8 m/s | 153 | 153, 1.41, 8.73 | 153, 0.45, 4.62 |
| Rayleigh damping above 700 m, absolute height | 121 | 126, 0.78, 5.69 | 121, 0.13, 0.56 |

With the side walls the columns differ across the channel; the script takes the fluxes as means over y, and they
agree to round-off. Under the sponge the solid cells stay at rest (|u| = 0). `bc.outflow_transport =
from_velocity` gave du = 6.9 % at F = 0.5 (the cut cell at 0.71 of the log law and w up to 0.19 m/s in the
outflow column) and is refused with an immersed ground.

## CTest
| Name | Label | What it checks |
|---|---|---|
| `rt_immersed_open_faces` | regression | F = 0, 0.5, 0.95 and 1 with the second setting, on one rank: the steady verdict, flux out = flux in within 1e-6, du < 1.5 %, dk < 8 %. |
| `rt_immersed_open_faces_sponge` | regression | F = 0.5 and 0.95, second setting with `bc.sponge_width = 320`, `bc.sponge_strength = 0.05`; same checks. |
| `rt_immersed_open_faces_power_law` | regression | F = 0.5 and 0.95, second setting with `bc.inflow_profile = power_law` (8 m/s at 10 m, alpha 0.143); same checks. |
| `rt_immersed_open_faces_wall_velocity` | regression | F = 0.5 and 0.95, second setting with no-slip side walls moving at 8 m/s; same checks with du < 2 %, dk < 10 %. |
| `rt_immersed_open_faces_rayleigh_absolute` | regression | F = 0.5 and 0.95, second setting with `bc.rayleigh_zdamp = 700` on u, v, w and k at the absolute height; same checks. |

The option tests stop at 1000 iterations (`steady.max_iterations`); every run takes 118 to 153.

Mutant: with the call that sets solid inflow faces to zero removed (`Source/Sequencing.cpp`, `set_solid_faces`),
`rt_immersed_open_faces` fails: at F = 0.95 the run has no steady verdict in 20000 iterations and the outflow
flux differs from the inflow flux by 7.6e-5 of it. The same change fails the four option tests at F = 0.95: no steady verdict in
1000 iterations, and a flux difference of 7.6e-5 (1.5e-4 with the power-law inflow).
