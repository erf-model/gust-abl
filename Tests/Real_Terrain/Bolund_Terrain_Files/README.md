# Bolund_Terrain_Files

## Purpose
A real-terrain case whose deck names both a terrain file and a roughness-map file: Bolund, wind from 270°, with
`terrain.file = bolund.xyz` and `bc.surface_z0_file = bolund_z0.xyz`. The folder also holds the script that
writes the two files from other formats. It is a set-up example; the comparison with the Bolund masts is in
`../WindNinja2019`.

## Files
| File | Role |
|---|---|
| `inputs_bolund_270` | The deck. 10.9 m/s from 270° at 5 m, the log law over water (z0 = 0.0003 m) at the inflow; 256 × 128 × 48 cells over 800 m × 420 m, a 0.5 m first layer; roughness 0.0003 m over water and 0.015 m over land from the map; mesh, solver and lee settings of `../WindNinja2019/inputs_bolund_270`. |
| `to_gust_xyz.py` | Terrain or roughness in another format as GUST point files (`x y z` or `x y z0` per line). Reads a single-band GeoTIFF (`.tif`, through `tools/tiff_to_xyz.py`, no GDAL), an ESRI ASCII grid (`.asc`), a Surfer ASCII grid (`.grd`) and point files (`x y value`, blanks or commas). `terrain` writes heights. `z0` writes roughness from a map of z0 values, from land-cover codes (`--classes CODE:Z0,...`) or from heights and a water level (`--water-level`). `--stride`, `--shift` and `--crop` thin, move and cut the points. `selftest` checks the readers and the three roughness modes. |

`bolund.xyz` and `bolund_z0.xyz` are not in the repository (see below); the script writes them.

## Data and acknowledgements
Downloaded by the user (not in the repository):
- **windninja-tools** (N. Wagenbrenner, <https://github.com/nwagenbrenner/windninja-tools>):
  `data/bolund/bolund_hill_extended.tif`, the Bolund terrain at 0.25 m in Bolund's own frame. The repository
  states no licence, so the file and the two point files made from it are not copied here.

The roughness lengths (0.0003 m water, 0.015 m land), the 0.75 m water level and the land strip east of
x = 327 m are those of the Bolund blind comparison's case description (Bechmann et al. 2011, *Boundary-Layer
Meteorol.* 141, 245); only the numbers are used.

Thanks to Natalie Wagenbrenner for making the terrain file available, and to the Bolund campaign and blind
comparison teams at DTU Wind Energy (Berg et al. 2011, *Boundary-Layer Meteorol.* 141, 219; Bechmann et al. 2011).

## How to run
```bash
T=/path/to/windninja-tools/data/bolund/bolund_hill_extended.tif
python3 to_gust_xyz.py terrain $T bolund.xyz --stride 4
python3 to_gust_xyz.py z0 $T bolund_z0.xyz --stride 4 --water-level 0.75 --land-x-min 327
mpirun -n 4 /path/to/gust inputs_bolund_270
```
Run `gust` from the folder that holds the two point files (the deck names them without a path). The two commands
give the same files as `../WindNinja2019/windninja2019.py prepare` followed by `bolund_z0_map.py`.

Other sources, as examples:
```bash
python3 to_gust_xyz.py terrain dem.asc terrain.xyz --shift 500000 4100000 0 --crop -2000 2000 -2000 2000
python3 to_gust_xyz.py z0 roughness.grd z0.xyz
python3 to_gust_xyz.py z0 landcover.asc z0.xyz --classes 11:0.0003,42:1.0,71:0.03
```
The points need not lie on the mesh: GUST interpolates them to its columns (`terrain.idw_neighbors`,
`terrain.idw_power`), the roughness in ln z0. The terrain and roughness files need not hold the same points.

## Expected Results
`to_gust_xyz.py` reports 336000 points for each file, x −399.75 to 399.25 m, y −209 to 210 m, heights 0.75 to
11.79 m, z0 0.0003 to 0.015 m. At start-up GUST prints
`gust: surface z0 map (bc.surface_z0_file): min 0.0003, max 0.015, geometric mean 0.0005428 m`.
The deck was started for two iterations on one rank (no unused input, plotfiles written); it has not been run to
the steady verdict, and no result is pinned.

## CTest
`rt_bolund_files_converter` (label `unit`): `to_gust_xyz.py selftest`. A 3 × 2 grid in the `.asc`, `.grd` and
point formats gives the same points; `--shift` and `--crop`; the water-level and class roughness modes; a code
without a z0 and a zero z0 are refused. The GeoTIFF reader is `tools/tiff_to_xyz.py`'s. The deck has no CTest,
as its two files are not in the repository.
