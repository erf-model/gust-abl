# Big Southern Butte, log-law inflow

## Purpose
Big Southern Butte, Idaho (43.396 N, 113.022 W), an isolated volcanic dome 771 m above the plain, with a log-law inflow (`bc.inflow_profile = log_law`) from the west through `xlo`. The inflow has 8 m/s at 10 m over the geometric mean roughness of the inflow face, z0 = 0.0301 m: u* = 0.41 x 8 / ln(10 / 0.0301) = 0.5650 m/s. GUST's log law, u*/κ ln((z + z0)/z0), then gives 8.004 m/s at 10 m. `xhi` is an outflow, the sides and the top are slip walls, and `abl.drive = inflow_only` keeps the inflow profile without a pressure gradient.

The mesh, the roughness map and the solver settings are those of the group (`../README.md`): 80 x 80 x 56 cells, 10.24 km x 10.24 km, the summit at the centre,
128 m cells, the 56 layers of `../zlevels_128m.txt` (2 m first layers, top 3063 m), STF terrain-following
coordinates, grid sequencing over 512, 256 and 128 m, and the recommended setting (`turb.realizability = 0.309`,
`turb.nut_relax = 0.5`, the upwind blend below 2.5 terrain heights, `turb.max_geom_lscale = 100`). Neutral: 300 K,
no buoyancy.

## Files
| File | Role |
|---|---|
| `inputs_butte_log_law` | The deck. It names no files; they are given on the command line. |
| `terrain.xyz` | `x y z` terrain heights (m) at the spacing of the tiles (about 45 to 49 m), origin at the site centre. |
| `terrain.sources.txt` | The sources of the terrain tiles. |
| `roughness.xyz` | `x y z0` roughness lengths (m) every 64 m from ESA WorldCover (`bc.surface_z0_file`). |

## Data and acknowledgements
The terrain is from the AWS Terrain Tiles (Mapzen and Tilezen); `terrain.sources.txt` lists the sources they merge.
The roughness is from ESA WorldCover 10 m 2021 v200 (Zanaga et al. 2022, doi:10.5281/zenodo.7254221, CC BY 4.0)
with the classes of Wieringa (1992). `../README.md` gives the full credits; thanks to the USGS, the Copernicus
programme, ESA and the WorldCover consortium, and Mapzen, Tilezen and AWS. The repository holds the derived
point files; `../prepare_sites.py` downloads the sources again to remake them.

## How to run
```bash
mpirun -np 4 <build>/gust inputs_butte_log_law terrain.file=terrain.xyz bc.surface_z0_file=roughness.xyz mesh.zlevels_file=../zlevels_128m.txt
```

## Expected Results
Steady after 6111 iterations on 2 ranks: 209 on the 512 m level, 130 on the 256 m level and 5772 on the
128 m level; final r 4.6e-5. With BTF coordinates (`terrain.smoothing = 0`) the 128 m level stalls near r = 8e-2.
`../README.md` compares the four cases.

## CTest
`rt_sites128_butte_log_law` (label `nightly`, 4 ranks) runs the deck to the steady verdict
(`Supporting_Files/expect_abort.py --success "gust: steady after"`).
