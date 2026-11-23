# Perdigão, target wind

## Purpose
Perdigão, Portugal (39.713 N, 7.737 W), the double ridge of the Perdigão campaign, 434 m of relief. The wind is given as a target, 10 m/s at 100 m from the south-west (the velocity points 45° from x, towards the north-east, across the NW–SE ridges), not as a geostrophic wind. GUST's flat column run (`inputs_perdigao_column`, `init.column_target_u`, `_v`, `_height`) finds the geostrophic wind that gives the target, (15.3874, 6.43277) m/s at 39.713 N over z0 = 0.398 m (the geometric mean of the map along `xlo` and `ylo`), and writes the converged column as a table. The wind turns from 47° at the ground to 22° aloft, so the table enters through `xlo` and `ylo` at every height; outflow through `xhi` and `yhi`. In the 3-D run Coriolis and that geostrophic wind keep the table in balance (`abl.drive = geostrophic_inflow`). The 3-D run starts from the same column (`init.first_guess = column`), which needs one roughness under the map: `init.column_z0 = 0.398` m, the roughness of the table. The 3-D run keeps the map. The terrain file has a flat border at the mean level (768 m flat, 1536 m slope).

The mesh, the roughness map and the solver settings are those of the group (`../README.md`): 96 x 96 x 56 cells, 12.29 km x 12.29 km,
128 m cells, the 56 layers of `../zlevels_128m.txt` (2 m first layers, top 3063 m), STF terrain-following
coordinates, grid sequencing over 512, 256 and 128 m, and the recommended setting (`turb.realizability = 0.309`,
`turb.nut_relax = 0.5`, the upwind blend below 2.5 terrain heights, `turb.max_geom_lscale = 100`). Neutral: 300 K,
no buoyancy.

The same site, doubly periodic with a geostrophic wind of 10 m/s along x, did not reach a steady state: on the
256 m level the relative residual stayed between 0.07 and 0.2. The ridges' wake came back through the periodic faces
after one pass (12.3 km at 10 m/s, about 20 minutes). With open faces the 256 m level reaches its loose stop after
200 iterations.

## Files
| File | Role |
|---|---|
| `inputs_perdigao_target_wind` | The deck. It names no files; they are given on the command line. |
| `terrain.xyz` | `x y z` terrain heights (m) at the spacing of the tiles (about 45 to 49 m), origin at the site centre. |
| `terrain.sources.txt` | The sources of the terrain tiles. |
| `roughness.xyz` | `x y z0` roughness lengths (m) every 64 m from ESA WorldCover (`bc.surface_z0_file`). |
| `inputs_perdigao_column` | The flat, doubly periodic 4 x 4 x 56 column run with the target wind that writes the table (`io.column_guess_file`). It stops after one SIMPLEC iteration: the table is written when the column guess is made. Its log gives the geostrophic wind ("the run takes geostrophic wind"), which `inputs_perdigao_target_wind` repeats in `forcing.ug` and `forcing.vg`. |
| `perdigao_inflow.txt` | The committed table (`z u v theta k`, u* = 0.7134 m/s). |

## Data and acknowledgements
The terrain is from the AWS Terrain Tiles (Mapzen and Tilezen); `terrain.sources.txt` lists the sources they merge.
The roughness is from ESA WorldCover 10 m 2021 v200 (Zanaga et al. 2022, doi:10.5281/zenodo.7254221, CC BY 4.0)
with the classes of Wieringa (1992). `../README.md` gives the full credits; thanks to the USGS, the Copernicus
programme, ESA and the WorldCover consortium, and Mapzen, Tilezen and AWS. The repository holds the derived
point files; `../prepare_sites.py` downloads the sources again to remake them.

## How to run
```bash
mpirun -np 4 <build>/gust inputs_perdigao_target_wind terrain.file=terrain.xyz bc.surface_z0_file=roughness.xyz bc.inflow_file=perdigao_inflow.txt mesh.zlevels_file=../zlevels_128m.txt
```
Remake the table in an empty folder (the deck writes `perdigao_inflow.txt` where it runs) and compare it with the
committed one. A new target needs the new geostrophic wind of the log in `forcing.ug` and `forcing.vg` of the 3-D deck.
```bash
mkdir new_table && cd new_table
```
```bash
mpirun -np 1 <build>/gust ../inputs_perdigao_column mesh.zlevels_file=../../zlevels_128m.txt
```
```bash
python3 ../../check_column_table.py perdigao_inflow.txt --committed ../perdigao_inflow.txt
```

## Expected Results
Steady after 4253 iterations on 2 ranks: 69 on the 512 m level, 200 on the 256 m level and 3984 on the
128 m level; final r 5.5e-5. `../README.md` compares the four cases.

## CTest
`rt_sites128_perdigao_target_wind` (label `nightly`, 4 ranks) runs the deck to the steady verdict
(`Supporting_Files/expect_abort.py --success "gust: steady after"`).

`rt_sites128_perdigao_column` (1 rank, labels `regression` and `smoke`) runs `inputs_perdigao_column`, and
`rt_sites128_perdigao_table` compares its table with the committed one using `../check_column_table.py` (heights to 1e-9 m; u, v, k and u* to a relative 1e-4; theta to 1e-3 K).
`rt_sites128_perdigao_clean` deletes the old table first.

`rt_sites128_perdigao_column_z0` runs the same column deck with the roughness map in place of `bc.surface_z0` and
`init.column_z0 = 0.398`; `rt_sites128_perdigao_column_z0_table` checks that its table is the committed one
(`init.column_z0 = 0.40` fails it: u differs by 6.5e-4). `rt_sites128_column_z0_needed` checks that the deck without
`init.column_z0` aborts. CMake writes both decks into the build folder from `inputs_perdigao_column`.
