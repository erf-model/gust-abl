# The Flatirons, inflow table

## Purpose
The Flatirons, Colorado (39.925 N, 105.285 W), where the plains meet the Front Range, 1037 m of relief, with an inflow table (`bc.inflow_profile = file`) through `xlo` and `ylo` and outflow through `xhi` and `yhi`. The table is GUST's own column for a geostrophic wind of 10 m/s from the west at 39.925 N over z0 = 0.404 m (the geometric mean of the map along `xlo`); it turns about 23° towards low pressure at the ground. In the 3-D run Coriolis keeps it in balance (`abl.drive = geostrophic_inflow`). Two limits: the z0 is that of `xlo` alone, although `ylo` is an inflow face too (its map gives 0.206 m, both faces 0.286 m); and between 1719 m and 2615 m the table's v is slightly negative (down to -0.147 m/s), so there `ylo` imposes a small outflow.

The mesh, the roughness map and the solver settings are those of the group (`../README.md`): 64 x 64 x 56 cells, 8.19 km x 8.19 km,
128 m cells, the 56 layers of `../zlevels_128m.txt` (2 m first layers, top 3063 m), STF terrain-following
coordinates, grid sequencing over 512, 256 and 128 m, and the recommended setting (`turb.realizability = 0.309`,
`turb.nut_relax = 0.5`, the upwind blend below 2.5 terrain heights, `turb.max_geom_lscale = 100`). Neutral: 300 K,
no buoyancy.

## Files
| File | Role |
|---|---|
| `inputs_flatirons_inflow_file` | The deck. It names no files; they are given on the command line. |
| `terrain.xyz` | `x y z` terrain heights (m) at the spacing of the tiles (about 45 to 49 m), origin at the site centre. |
| `terrain.sources.txt` | The sources of the terrain tiles. |
| `roughness.xyz` | `x y z0` roughness lengths (m) every 64 m from ESA WorldCover (`bc.surface_z0_file`). |
| `inputs_flatirons_column` | The flat, doubly periodic 4 x 4 x 56 column run that writes the table (`io.column_guess_file`). It stops after one SIMPLEC iteration: the table is written when the column guess is made. |
| `flatirons_inflow.txt` | The committed table (`z u v theta k`, u* = 0.4711 m/s). |

## Data and acknowledgements
The terrain is from the AWS Terrain Tiles (Mapzen and Tilezen); `terrain.sources.txt` lists the sources they merge.
The roughness is from ESA WorldCover 10 m 2021 v200 (Zanaga et al. 2022, doi:10.5281/zenodo.7254221, CC BY 4.0)
with the classes of Wieringa (1992). `../README.md` gives the full credits; thanks to the USGS, the Copernicus
programme, ESA and the WorldCover consortium, and Mapzen, Tilezen and AWS. The repository holds the derived
point files; `../prepare_sites.py` downloads the sources again to remake them.

## How to run
```bash
mpirun -np 4 <build>/gust inputs_flatirons_inflow_file terrain.file=terrain.xyz bc.surface_z0_file=roughness.xyz bc.inflow_file=flatirons_inflow.txt mesh.zlevels_file=../zlevels_128m.txt
```
Remake the table in an empty folder (the deck writes `flatirons_inflow.txt` where it runs) and compare it with the
committed one:
```bash
mkdir new_table && cd new_table
```
```bash
mpirun -np 1 <build>/gust ../inputs_flatirons_column mesh.zlevels_file=../../zlevels_128m.txt
```
```bash
python3 ../../check_column_table.py flatirons_inflow.txt --committed ../flatirons_inflow.txt
```

## Expected Results
Steady after 7275 iterations on 2 ranks: 258 on the 512 m level, 388 on the 256 m level and 6629 on the
128 m level; final r 8.8e-5. `../README.md` compares the four cases.

## CTest
`rt_sites128_flatirons_inflow_file` (label `nightly`, 4 ranks) runs the deck to the steady verdict
(`Supporting_Files/expect_abort.py --success "gust: steady after"`).

`rt_sites128_flatirons_column` (1 rank, labels `regression` and `smoke`) runs `inputs_flatirons_column`, and `rt_sites128_flatirons_table` compares its table with the committed one using `../check_column_table.py` (heights to 1e-9 m; u, v, k and u* to a relative 1e-4; theta to 1e-3 K). `rt_sites128_flatirons_clean` deletes the old table first.
