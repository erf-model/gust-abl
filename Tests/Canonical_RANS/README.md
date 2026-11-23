# Canonical_RANS

## Purpose
The canonical cases: the SIMPLEC stage run to the steady verdict, the unsteady stage, forest canopies, buildings
and nests. The folders that ERF also has keep the names of ERF's `Exec/CanonicalTests/Canonical_RANS/`, and those
cases run ERF's own check script (`extern_erf/`) on the turbulence plotfile. Each folder's `README.md` gives its
set-up, expected results and the full list of its CTests.

## Cases
| Folder | CTest | Description |
|---|---|---|
| `Neutral_Column` | `rt_steady_column` | The 1-D neutral column driven by a pressure-gradient force; a restart chain that must continue bit for bit. |
| `Immersed_Column` | `rt_immersed_column_*` | The neutral column with the ground immersed in a flat mesh (`terrain.representation = immersed_forcing`) at four cut-cell solid fractions: the exact wall stress, the log law, and the body-fitted column at f = 0. |
| `Immersed_Ridge_2D` | `rt_immersed_along_*`, `rt_immersed_ridge_*`, `rt_immersed_hill_*` (nightly) | A sloping immersed ground (ridge with the wind along and across it, round hill) against the same ground on a terrain-following mesh. |
| `Immersed_Ridge_Resolution` | none | The crest speed-up of a steep 2-D immersed ridge against its terrain-following partner as the column width and layer thickness change: the record behind the docs' grid-resolution guide. |
| `Immersed_Open_Faces` | `rt_immersed_open_faces` | Flat immersed ground between an inflow and an outflow face at four cut-cell solid fractions, against the body-fitted channel: the steady verdict, the mass balance of the two faces, u and k in the outflow column. |
| `Immersed_Roughness_Step` | `rt_immersed_roughness_step`, `*_z0_map_*` | A roughness map under an immersed ground and under the ground of a hybrid mesh, against the same flow on a body-fitted mesh with the same map; the map's restart and refusals. |
| `Immersed_Canopy` | `rt_immersed_canopy_*` | A forest canopy over an immersed ground and on a hybrid mesh, against the same flow on a body-fitted mesh, on the decks of three canopy cases. |
| `Neutral_ABL_Flat` | `rt_steady_neutral_flat` | ERF's Neutral_ABL_Flat to steady, with ERF's `check_neutral.py --physics`. |
| `Neutral_Hill_2D` | `rt_steady_hill2d`, `rt_steady_hill2d_lusgs`, `rt_steady_hill2d_column` | ERF's Neutral_Hill_2D with Jacobi sweeps, line LU-SGS sweeps and the column first guess. |
| `Neutral_Hill_3D` | `rt_steady_hill3d` | ERF's Neutral_Hill_3D (64 x 64 x 20); label `nightly`. |
| `Truly_Neutral_Hill_3D` | `rt_steady_hill3d_neutral_zjacobi`, `rt_steady_hill3d_neutral_lusgs` | A truly neutral Hill_3D to round-off with z-line Jacobi and box-local line LU-SGS sweeps; the two states agree to 1e-8. Label `nightly`. |
| `Thin_Cell_Flat` | `rt_thin_cell_w_vnn_cap` | A flat inflow/outflow strip with a 0.7 m first cell, with and without `simplec.w_vnn_cap`: fewer iterations, the same steady state, a restart bit for bit. |
| `Unsteady_Hill_2D` | `rt_unsteady_hill2d` | The 2-D hill as unsteady RANS (dual time stepping): ten steps of 60 s on 1, 2 and 4 ranks, the same steps and the same answer (1e-9). |
| `Rampart_2D` | `rt_rampart_ground_step` | A 900 m analytic cliff on 200 m cells: the mesh report's ground step between neighbouring columns against closed forms, its warning above 64, and the iterations of a gentle and a steep cliff at the same ground step (260 and 2570). |
| `Terrain_Shapes_3D` | `rt_terrain_shapes_settle` | Fourteen analytic 760 m shapes (buttes of 32° to 75°, ridges, a gap, a canyon, a crater rim) with 2560 probes each: which shapes the steady stage settles on, by the residual and by the wind at the points (0.5 and 0.25 m/s). The CTest runs a gentle and a steep butte on 200 m cells. |
| `Steep_Ridge_Thin_Cell` | `rt_steep_ridge_w_vnn_cap`, `rt_steep_ridge_w_vnn_cap_margin_off` | A steep 2-D ridge (32.5°) with a 0.4 m first cell: `simplec.w_vnn_cap = 16` with the default slope margin is steady within 3500 iterations, without the margin it is not. |
| `Snyder_Lawson_Cube` | `snyder_cube_plain_*`, `snyder_cube_lee_*`, `snyder_cube_lee_nozone_*` | A cube in a wind-tunnel boundary layer against the measurements of Snyder and Lawson (1994); the closure's lee options with and without `buildings.limiter_zone`. |
| `Castro_Cube_Array` | `castro_cube_array_*` | A staggered cube array (one periodic tile) against the measurements of Castro et al. (2006), as Muñoz-Esparza et al. (2020) set it up. |
| `Nest_Validation` | `rt_nest_validation_*` | One-way nesting against a uniform fine run of the whole domain over a 150 m hill: one nest at four patch sizes, a chain of three nests, hard faces against the relaxation zone. |
| `Nest_Faces` | `rt_nest_faces_*`, `abort_nest_faces_*` | A patch at 2 2 1 with its four lateral faces from the parent run (`bc.* = nest`), on flat ground and over a hill, to steady on 1, 2 and 4 ranks; the relaxation zone, the ground blend and the restart record; a patch that holds its parent's theta (`simplec.theta = frozen`). |
| `Nest_Shapes` | `rt_nest_shapes_*` | One-way nested patches at 2 2 1 around six canonical shapes against a uniform fine run of the whole domain; the cube runs are CTests. |
| `Nests_One_Process` | `rt_nests_one_process`, `abort_nests_*` | Nests named in the outer deck (`nest.names`) and run in the same process; the refusals of misplaced nests. |
| `Nests_Two_Way` | `rt_nests_two_way`, `rt_nests_joint_restart*` | Two-way nesting as a joint run (`nest.feedback = defect`, `relax`, `relax_consistent`, `inject`) and its restart. |
| `Nests_Flat_ABL` | `rt_nests_flat_abl_*` | Nests on a flat, horizontally uniform boundary layer, neutral, buoyant and stable, one way and joint. |
| `Nest_Placement_Runs` | (study) | Fifteen runs (three uniform references and twelve joint runs) behind the step rule of `tools/nest_placement.py`; not a CTest. |
| `Canopy_Drag` | `rt_canopy_drag`, `rt_canopy_tke`, `rt_canopy_column` | The canopy's drag with each implicit diagonal, the canopy's k terms and the column first guess through two forest patches. |
| `Homogeneous_Canopy` | `rt_homogeneous_canopy` | A horizontally uniform forest under a pressure gradient at LAI 2, 5 and 10 against published canopy relations. |
| `Forest_Edge` | `rt_forest_edge` | A forest, clearing, forest pattern (the geometry of Dupont and Brunet 2008): the adjustment behind the leading edge. |
| `Forested_Hill` | `rt_forested_hill` | A uniform canopy on a periodic train of cosine hills, in two geometries classified on the regime map of Poggi et al. (2008). |
| `Forested_Escarpment` | `rt_forested_escarpment`, `rt_frozen_theta`, `rt_canopy_heat` | A dense canopy on a 200 m escarpment with an inflow face and grid sequencing; frozen theta; the leaves' heat exchange. |
| `Building_Cube` | `building_cube_*`, `force_file_cube_*` | One block on flat ground in a log-law wind with each wall form; `buildings.force_file` in an unsteady run. |
| `Building_Courtyard` | `building_courtyard_*` | A courtyard open to the wind: the wall law's speed in cells with two walls. |
| `Building_Hill` | `building_hill_*` | One block on the windward flank of a 30 m ridge: the building's floor and roof follow the ground. |
| `Cylinder_Laminar` | `rt_cylinder_re40`, `abort_laminar_*` | Laminar flow past a staircase cylinder (`turb.closure = none`): steady at Re 40 (CTest); the unsteady runs at Re 40 and 100 are comparison runs. |
| `Column_Target` | `rt_column_target`, `rt_column_target_mol` | The column first guess's target wind carried into the 3-D run, neutral and stable. |
| `Column_Marching` | (comparison case) | The time-marching column against ERF's near-1-D columns; not a CTest (unit tests `ColumnMarching.*`). |
| `Aspect_Ratio_Hill_2D` | `rt_aspect_ratio_hill` | A 2-D ridge with first cells of dx/dz1 = 10 to 5000: how thin first cells change the convergence. |
| `Outflow_Steep_Terrain` | `rt_outflow_steep_terrain_*` | A steep slope cut by an outflow face: steady with the default outflow transport, not with the old one. |
| `Smooth_Plateau_2D` | (comparison case) | A 2-D plateau run in both codes; not a CTest. |

The hill cases use `Supporting_Files/check_steady_hill.py` (shared by the 2-D and 3-D folders).
