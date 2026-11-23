# Altamont Pass, power-law inflow

## Purpose
Altamont Pass, California (37.740 N, 121.660 W), rolling hills with 423 m of relief, with a power-law inflow (`bc.inflow_profile = power_law`) from the west through `xlo`: 8 m/s at 10 m, exponent 0.143, and z0 = 0.0591 m (the geometric mean of the map along `xlo`) for the wall law and the inflow k. `xhi` is an outflow, the sides and the top are slip walls, and `abl.drive = inflow_only`.

The mesh, the roughness map and the solver settings are those of the group (`../README.md`): 64 x 64 x 56 cells, 8.19 km x 8.19 km,
128 m cells, the 56 layers of `../zlevels_128m.txt` (2 m first layers, top 3063 m), STF terrain-following
coordinates, grid sequencing over 512, 256 and 128 m, and the recommended setting (`turb.realizability = 0.309`,
`turb.nut_relax = 0.5`, the upwind blend below 2.5 terrain heights, `turb.max_geom_lscale = 100`). Neutral: 300 K,
no buoyancy.

## Files
| File | Role |
|---|---|
| `inputs_altamont_power_law` | The deck. It names no files; they are given on the command line. |
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
mpirun -np 4 <build>/gust inputs_altamont_power_law terrain.file=terrain.xyz bc.surface_z0_file=roughness.xyz mesh.zlevels_file=../zlevels_128m.txt
```

## Expected Results
Steady after 617 iterations on 2 ranks: 46 on the 512 m level, 56 on the 256 m level and 515 on the
128 m level; final r 4.2e-5. `../README.md` compares the four cases.

## CTest
`rt_sites128_altamont_power_law` (label `nightly`, 4 ranks) runs the deck to the steady verdict
(`Supporting_Files/expect_abort.py --success "gust: steady after"`).
