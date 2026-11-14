# Canonical_RANS

## Purpose
The SIMPLEC stage run to the steady verdict on the canonical RANS cases. The folder names follow ERF's
`Exec/CanonicalTests/Canonical_RANS/`, and each case runs ERF's own check script (`extern_erf/`) on the
turbulence plotfile where ERF has one.

## Cases
| Folder | CTest | Description |
|---|---|---|
| `Neutral_Column` | `rt_steady_column` | The 1-D neutral column driven by a pressure-gradient force; a restart chain that must continue bit for bit. |
| `Neutral_ABL_Flat` | `rt_steady_neutral_flat` | ERF's Neutral_ABL_Flat to steady, with ERF's `check_neutral.py --physics`. |
| `Neutral_Hill_2D` | `rt_steady_hill2d`, `rt_steady_hill2d_lusgs`, `rt_steady_hill2d_column` | ERF's Neutral_Hill_2D with Jacobi sweeps, line LU-SGS sweeps and the column first guess. |
| `Neutral_Hill_3D` | `rt_steady_hill3d` | ERF's Neutral_Hill_3D (64 x 64 x 20); label `nightly`. |
| `Truly_Neutral_Hill_3D` | `rt_steady_hill3d_neutral_zjacobi`, `rt_steady_hill3d_neutral_lusgs` | A truly neutral Hill_3D to round-off with z-line Jacobi and box-local line LU-SGS sweeps; the two states agree to 1e-8. Label `nightly`. |
| `Thin_Cell_Flat` | `rt_thin_cell_w_vnn_cap` | A flat inflow/outflow strip with a 0.7 m first cell, with and without `simplec.w_vnn_cap`: fewer iterations, the same steady state, a restart bit for bit. |
| `Unsteady_Hill_2D` | `rt_unsteady_hill2d` | The 2-D hill as unsteady RANS (dual time stepping): ten steps of 60 s on 1, 2 and 4 ranks, the same steps and the same answer (1e-9). |
| `Steep_Ridge_Thin_Cell` | `rt_steep_ridge_w_vnn_cap`, `rt_steep_ridge_w_vnn_cap_margin_off` | A steep 2-D ridge (32.5°) with a 0.4 m first cell: `simplec.w_vnn_cap = 16` with the default slope margin is steady within 3500 iterations, without the margin it is not. |
| `Snyder_Lawson_Cube` | `snyder_cube_plain_*`, `snyder_cube_lee_*`, `snyder_cube_lee_nozone_*` | A cube in a wind-tunnel boundary layer against the measurements of Snyder and Lawson (1994); the closure's lee options with and without `buildings.limiter_zone`. |
| `Castro_Cube_Array` | `castro_cube_array_*` | A staggered cube array (one periodic tile) against the measurements of Castro et al. (2006), as Muñoz-Esparza et al. (2020) set it up. |
| `Nest_Validation` | `rt_nest_validation_*` | One-way nesting against a uniform fine run of the whole domain over a 150 m hill: one nest at four patch sizes, a chain of three nests, hard faces against the relaxation zone. |
| `Nest_Faces` | `rt_nest_faces_*`, `abort_nest_faces_*` | A patch at 2 2 1 with its four lateral faces from the parent run (`bc.* = nest`), on flat ground and over a hill, to steady on 1, 2 and 4 ranks; the relaxation zone, the ground blend and the restart record. |
| `Smooth_Plateau_2D` | (comparison case) | A 2-D plateau run in both codes; not a CTest. |

The hill cases use `Supporting_Files/check_steady_hill.py` (shared by the 2-D and 3-D folders).
