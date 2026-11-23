# Sites_128m

## Purpose
Four real terrains on one mesh recipe, each with a different inflow option of GUST:

| Folder | Site | Inflow | Domain | Relief on the 128 m mesh | Steepest slope (128 m) |
|---|---|---|---|---|---|
| `BigSouthernButte_LogLaw` | Big Southern Butte, Idaho (43.396 N, 113.022 W) | log law (`bc.inflow_profile = log_law`) through xlo | 80 x 80 cells, 10.24 km | 771 m | 46° |
| `AltamontPass_PowerLaw` | Altamont Pass, California (37.740 N, 121.660 W) | power law (`power_law`) through xlo | 64 x 64 cells, 8.19 km | 423 m | 28° |
| `Perdigao_TargetWind` | Perdigão, Portugal (39.713 N, 7.737 W) | a target wind of 10 m/s at 100 m from the south-west (the velocity points 45° from x): GUST's column finds the geostrophic wind that gives it (`init.column_target_u`, `_v`, `_height`) and writes a table, which enters through xlo and ylo (`abl.drive = geostrophic_inflow`); the 3-D first guess is the column on one roughness under the map (`init.column_z0`) | 96 x 96 cells, 12.29 km | 434 m | 38° |
| `Flatirons_InflowFile` | The Flatirons, Colorado (39.925 N, 105.285 W) | a table (`file`) through xlo and ylo: GUST's column for a geostrophic wind, kept in balance by Coriolis (`abl.drive = geostrophic_inflow`) | 64 x 64 cells, 8.19 km | 1037 m | 55° |

The mesh is the same for all four:
- 128 m x 128 m cells, terrain-following with the smoothed coordinate surfaces aloft (`terrain.smoothing = 1`, STF,
  the default). Without smoothing (BTF, as in the `WindNinja2019` Butte deck) the Butte stalls on the 128 m level
  at a relative residual near 8e-2;
- 56 layers (`zlevels_128m.txt`): 10 layers of 2 m, then layers that grow by 1.15 until they reach 128 m, then
  128 m layers to 3063 m;
- grid sequencing (`mesh.hierarchy_levels = 3`, `steady.sequencing = 1`): the run starts on 512 m and 256 m levels
  and finishes on the 128 m level.

The solver settings are those of the converged `WindNinja2019` Butte deck: SIMPLEC with LU-SGS
(`steady.zline_sweeps = 2`), `steady.cfl_max = 16`, `simplec.w_vnn_cap = 16` for the thin first cells (aspect ratio 64)
and `turb.max_geom_lscale = 100`.

Each case uses the roughness map of its site (`bc.surface_z0_file`) and the recommended setting of
`docs/turbulence_model.rst`: `turb.realizability = 0.309`, `turb.nut_relax = 0.5` and the blend towards
first-order upwind below 2.5 terrain heights. The flow is truly neutral: 300 K, no buoyancy. Only the two
geostrophic cases have Coriolis.

The decks do not name their files. The terrain, the roughness map, the levels and, for the Flatirons and Perdigão, the
inflow table are given on the command line, as in the other terrain CTests.

## Files
| File | Role |
|---|---|
| `zlevels_128m.txt` | The 57 level heights of the four decks (`mesh.zlevels_file`). |
| `prepare_sites.py` | Writes `terrain.xyz`, `terrain.sources.txt` and `roughness.xyz` of each site (needs numpy, PyYAML, rasterio and network access). |
| `<case>/inputs_*` | The deck of each case. |
| `<case>/terrain.xyz` | `x y z` points (m) at the spacing of the tiles (about 45 to 49 m) over the domain and 256 m beyond it. UTM zone of the site centre, origin at the centre, heights above the lowest point kept. Perdigão has a flat border at the mean level (768 m flat, 1536 m slope), so its inflow faces are flat like the column of its table. |
| `<case>/terrain.sources.txt` | The data sources of the terrain tiles, as the download recorded them. |
| `<case>/roughness.xyz` | `x y z0` points (m) every 64 m, two per 128 m cell (below). |
| `Flatirons_InflowFile/inputs_flatirons_column` | The flat, doubly periodic column run that writes `flatirons_inflow.txt`. |
| `Flatirons_InflowFile/flatirons_inflow.txt` | The inflow table (`z u v theta k`). |
| `Perdigao_TargetWind/inputs_perdigao_column` | The flat, doubly periodic column run with the target wind that writes `perdigao_inflow.txt`. |
| `Perdigao_TargetWind/perdigao_inflow.txt` | The inflow table (`z u v theta k`). |
| `check_column_table.py` | Compares a table that a column run just wrote with the committed one (`--committed`). |

### Roughness maps
`prepare_sites.py` reads the ESA WorldCover 2021 land cover (10 m pixels). It gives each class a roughness length,
using the class of the updated Davenport classification (Wieringa 1992) that the cover describes. Trees, built-up
land, water and bare ground take the values of the ESA CCI table in `tools/gust_case/roughness.py`. The 10 m values
are then averaged in ln z0 over 64 m squares.

| WorldCover class | z0 (m) |
|---|---|
| 10 tree cover | 1.0 |
| 20 shrubland | 0.25 |
| 30 grassland | 0.03 |
| 40 cropland | 0.1 |
| 50 built-up | 0.8 |
| 60 bare or sparse vegetation | 0.005 |
| 70 snow and ice | 0.001 |
| 80 permanent water | 1e-4 |
| 90 herbaceous wetland | 0.2 |
| 95 mangroves | 1.0 |
| 100 moss and lichen | 0.03 |

The inflow roughness of each case is a geometric mean of the map near its inflow faces. `python3 prepare_sites.py --faces` prints the means from the committed files (no network). The Butte, Altamont and the Flatirons use the points of the first 128 m cell inside xlo. Perdigão uses the points within 64 m on either side of xlo and ylo (0.398 m; the first cell inside gives 0.417 m). The Flatirons table enters through xlo and ylo, but its z0 is that of xlo alone (ylo gives 0.206 m, both faces 0.286 m).

| Case | z0 range of the map (m) | Geometric mean (m) | Inflow z0 (m) |
|---|---|---|---|
| Butte | 0.005–1 | 0.032 | 0.0301 (xlo) |
| Altamont | 1e-4–1 | 0.033 | 0.0591 (xlo) |
| Perdigão | 1e-4–1 | 0.519 | 0.398 (xlo and ylo) |
| Flatirons | 1e-4–1 | 0.243 | 0.404 (xlo) |

## Data and acknowledgements
- **Terrain:** the Terrain Tiles on the AWS Open Data Registry (https://registry.opendata.aws/terrain-tiles/),
  made by Mapzen and Tilezen (joerd), accessed 2026-11-23. The three US sites are from the United States 3DEP
  (formerly NED) and SRTM data, courtesy of the U.S. Geological Survey. Perdigão is from the EU-DEM layers,
  produced using Copernicus data and information funded by the European Union. `terrain.sources.txt` in each
  folder lists every source the tiles merge (https://github.com/tilezen/joerd/blob/master/docs/attribution.md).
  Thanks to the USGS, the Copernicus programme, and Mapzen, Tilezen and AWS for the tiles.
- **Land cover:** ESA WorldCover 10 m 2021 v200. D. Zanaga, R. Van De Kerchove, D. Daems, W. De Keersmaecker,
  C. Brockmann, G. Kirches, J. Wevers, O. Cartus, M. Santoro, S. Fritz, M. Lesiv, M. Herold, N.-E. Tsendbazar,
  P. Xu, F. Ramoino and O. Arino (2022), doi:10.5281/zenodo.7254221, licence CC BY 4.0. Read from the
  cloud-optimized GeoTIFFs on the AWS Open Data Registry (https://registry.opendata.aws/esa-worldcover-vito/).
  Thanks to ESA and the WorldCover consortium.
- **Roughness classes:** J. Wieringa (1992), "Updating the Davenport roughness classification", J. Wind Eng. Ind.
  Aerodyn. 41, 357–368, doi:10.1016/0167-6105(92)90434-C.

The repository holds the derived point files: terrain heights at the spacing of the tiles and 64 m roughness lengths, about 6 MB. Both
licences allow this with attribution. Nothing has to be downloaded to run the cases. `prepare_sites.py` downloads
the sources again to remake the files.

## How to run
From a case folder, for example the Butte:
```bash
mpirun -np 4 <build>/gust inputs_butte_log_law terrain.file=terrain.xyz bc.surface_z0_file=roughness.xyz mesh.zlevels_file=../zlevels_128m.txt
```
The Flatirons and Perdigão cases read their committed tables (each case README shows how to remake the table in
another folder, so the committed file is not overwritten):
```bash
mpirun -np 4 <build>/gust inputs_flatirons_inflow_file terrain.file=terrain.xyz bc.surface_z0_file=roughness.xyz bc.inflow_file=flatirons_inflow.txt mesh.zlevels_file=../zlevels_128m.txt
```
```bash
mpirun -np 4 <build>/gust inputs_perdigao_target_wind terrain.file=terrain.xyz bc.surface_z0_file=roughness.xyz bc.inflow_file=perdigao_inflow.txt mesh.zlevels_file=../zlevels_128m.txt
```
Remake the point files (network access):
```bash
python3 prepare_sites.py
```

## Expected Results
Runs on 2 MPI ranks (Apple M-series, Release build), 2026-11-23; the nightly tests use 4 ranks, where the counts
can differ slightly. Each run stops when the relative residual r is
below 1e-4 and every scaled residual is below 5e-5 (`steady.residual_tol`) and 4 orders below its reference. On all
four the last residual to reach its floor is w.

| Case | Inflow | 512 m level | 256 m level | 128 m level | Total | Final r | Final w-res |
|---|---|---|---|---|---|---|---|
| Altamont | power law | 46 | 56 | 515 | **617** | 4.2e-5 | 4.98e-5 |
| Butte | log law | 209 | 130 | 5772 | **6111** | 4.6e-5 | 5.00e-5 |
| Perdigão | target wind | 69 | 200 | 3984 | **4253** | 5.5e-5 | 5.00e-5 |
| Flatirons | inflow table | 258 | 388 | 6629 | **7275** | 8.8e-5 | 5.00e-5 |

The coarse levels stop at 2 orders (the loose stop of grid sequencing). Two set-ups did not reach a steady state
and are not committed:
- the Butte with BTF coordinates (`terrain.smoothing = 0`): r near 8e-2 on the 128 m level after 2150
  iterations. Nearly all of each residual sat in the first 20 m on slopes of 15° or more; with STF those residuals
  are 140 to 60 000 times smaller;
- Perdigão doubly periodic with a geostrophic wind of 10 m/s along x: on the 256 m level r stayed between 0.07
  and 0.2 for 1000 iterations, and on the 128 m level near 0.15. The ridges' wake came back through the periodic
  faces. Most of the residual sat in the first 20 m on steep upper flanks, not in the valley air. With open faces
  (the committed deck) the 256 m level stops after 200 iterations.

Choosing the column first guess (`init.first_guess = column`, `init.column_z0`) for Perdigão instead of the
profile first guess (`init.first_guess = profile`) gives the same iterations (4253), since both start from the
same table.

## CTest
| Test | Ranks | Labels | What it checks |
|---|---|---|---|
| `rt_sites128_butte_log_law`, `rt_sites128_altamont_power_law`, `rt_sites128_perdigao_target_wind`, `rt_sites128_flatirons_inflow_file` | 4 | `nightly` | The deck reaches "gust: steady after" (`Supporting_Files/expect_abort.py --success`, timeout 3600 s). |
| `rt_sites128_flatirons_clean`, `rt_sites128_perdigao_clean`, `rt_sites128_perdigao_column_z0_clean` | – | `regression`, `smoke` | Delete the table of an earlier run, so it cannot pass the table test. |
| `rt_sites128_flatirons_column`, `rt_sites128_perdigao_column` | 1 | `regression`, `smoke` | The column deck runs and writes its table. |
| `rt_sites128_flatirons_table`, `rt_sites128_perdigao_table` | – | `regression`, `smoke` | `check_column_table.py` compares the new table with the committed one (heights to 1e-9 m; u, v, k and u* to a relative 1e-4; theta to 1e-3 K). A target height of 101 m instead of 100 m fails it (u differs by 2.8e-3). |
| `rt_sites128_perdigao_column_z0`, `rt_sites128_perdigao_column_z0_table` | 1 | `regression`, `smoke` | The Perdigão column deck with the roughness map and `init.column_z0 = 0.398` in place of `bc.surface_z0` writes the committed table (`init.column_z0 = 0.40` fails it). |
| `rt_sites128_column_z0_needed` | 1 | `regression`, `smoke` | The same deck without `init.column_z0` aborts. |
