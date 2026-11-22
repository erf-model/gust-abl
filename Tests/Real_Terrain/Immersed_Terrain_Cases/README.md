# Immersed_Terrain_Cases

## Purpose
The measured cases of `../WindNinja2019` and `../IshiharaZhou` run twice in the same settings: on a
terrain-following mesh (`body_fitted`) and with the ground immersed in a flat mesh
(`terrain.representation = immersed_forcing`, [docs](../../../docs/immersed_terrain.rst)). The pair differs in the
ground's representation only, so it shows what the immersed ground costs in cells and iterations and how far its
winds are from the terrain-following run and from the measurements.

| Case | Wind | Body-fitted cells | Immersed cells | Layers of the immersed mesh |
|---|---|---|---|---|
| `bolund_270`, `bolund_239`, `bolund_90` | Bolund cases 1, 3, 4 | 256 x 128 x 48 | 256 x 128 x 74 | 0.5 m up to 14 m, then stretched |
| `askervein_210` | 8.9 m/s from 210° at 10 m | 256 x 256 x 64 | 256 x 256 x 114 | 2 m up to 128 m, then stretched |
| `ridge_smooth`, `ridge_rough` | Ishihara and Zhou (2025), H = 40 m | 600 x 8 x 64 | 600 x 8 x 147 | 0.5 m up to 42 m, then stretched |
| `butte_222` | Big Southern Butte, 8.3 m/s from 222° at 3 m | 136 x 148 x 40 | 136 x 148 x 131 | 8 m up to 808 m, then stretched |

Every deck has one mesh level, `steady.sweep = zjacobi`, `solver.poisson_backend = gmres_fft`,
`turb.shear_stencil = log_law` and the recommended turbulence setting (`docs/turbulence_model.rst`:
`turb.realizability = 0.309`, `turb.nut_relax = 0.5`, the blended convection scheme). The immersed decks add
`terrain.small_volfrac = 0.1` and their layers (`mesh.zlevels_file`): uniform layers over the height range of the
ground, since the ground can cut any of them. The Butte decks have half the horizontal resolution of
`../WindNinja2019/inputs_bsb_srtm` (136 m).

## Files
| File | Role |
|---|---|
| `inputs_<case>_body_fitted`, `inputs_<case>_immersed` | The pairs. |
| `zlevels_bolund.txt`, `zlevels_askervein.txt`, `zlevels_ridge.txt`, `zlevels_butte.txt` | The layers of the immersed meshes. |
| `sample_masts.py` | Speed-up at the masts and along the lines from either kind of plotfile. |
| `compare_masts.py` | Bolund and Askervein: both runs against the masts; sections of speed and k. |
| `compare_ridge.py` | The ridge: both runs against the wind-tunnel stations. |
| `compare_butte.py` | The Butte: runs against the 3 m sensors and against each other. |
| `compare_fields.py` | Any pair: speed-up and k at 10 and 40 m above the ground, column by column, immersed against body-fitted, by region (slopes, upper hill, flat ground). |
| `check_coarse_pair.py` | The CTest: a coarse pair of one case, run and compared. |

## Data and acknowledgements
Terrain and measurements are downloaded by the user; none of it is in this folder. The decks take `terrain.file` (and the ridge's `bc.inflow_file`) on the command line, and the scripts take the
data folder as an argument: these are the folders that `../WindNinja2019/windninja2019.py prepare` and
`../IshiharaZhou/ishihara_zhou.py prepare` write (their READMEs give the steps). The Askervein terrain of the
CTest is `extern_erf/askervein.txt`, which is in the repository.

- **windninja-tools** (N. Wagenbrenner, <https://github.com/nwagenbrenner/windninja-tools>): the Askervein masts,
  the Bolund masts and terrain, the Butte terrain. The repository states no licence, so none of its files is
  copied here.
- **Wagenbrenner, Forthofer, Page and Butler (2019)**, *Atmosphere* 10, 672, doi:10.3390/atmos10110672
  (CC BY 4.0): the cases and their settings.
- **Big Southern Butte**: Butler et al. (2015) and Forthofer et al. (2014), Missoula Fire Sciences Laboratory;
  data at Zenodo, doi:10.5281/zenodo.8226638.
- **Bolund**: Berg et al. (2011) and Bechmann et al. (2011), DTU Wind Energy.
- **Askervein**: Taylor and Teunissen (1987).
- **Ishihara and Zhou (2025)**, *Data in Brief* 63, 112260, doi:10.1016/j.dib.2025.112260 (CC BY 4.0): the ridge
  and its wind-tunnel stations.

We thank Natalie Wagenbrenner and the Missoula Fire Sciences Laboratory, DTU Wind Energy and the Bolund team, the
Askervein Hill project, and Takeshi Ishihara and Tong Zhou for making these data available.

## How to run
Each run has its own folder. The scripts expect in it the deck as `deck`, the solver's output as `log`, the
plotfiles `plt?????` and, for an immersed run, its levels as `zlevels.txt`:
```bash
mkdir bolund_270_immersed && cd bolund_270_immersed
cp ../inputs_bolund_270_immersed deck && cp ../zlevels_bolund.txt zlevels.txt
mpirun -np 4 <build>/gust deck terrain.file=<data>/bolund.xyz mesh.zlevels_file=zlevels.txt > log
```
Then, with `<data>` the folder the `prepare` step wrote:
```bash
python3 compare_masts.py bolund_270 <body-fitted folder> <immersed folder> bolund_270.png <data>
python3 compare_ridge.py smooth <body-fitted folder> <immersed folder> ridge_smooth.png <data>
python3 compare_butte.py butte.png <data> body_fitted=<folder> immersed=<folder>
```
For any pair, the speed-up and k at 10 and 40 m above the ground by region (the flat height is the immersed deck's
`terrain.flat_height`, the wind direction is the case's, in degrees):
```bash
python3 compare_fields.py --body-fitted <body-fitted folder> --immersed <immersed folder> --flat 3 --wind-dir 210
```
The full-size runs take minutes (Bolund, the ridge) to hours (Askervein, the Butte) on four ranks.

## Expected Results
Full-size runs, both representations with the same settings (one level, `zjacobi`,
`turb.shear_stencil = log_law`, the recommended limiter set, `terrain.small_volfrac = 0.1`, no side walls on the
ground, the fluid's mean height above the ground as the wall distance). Speed-up errors are against the measurements; the ridge is against
the wind-tunnel data.

| Case | Measure | Body-fitted | Immersed | Iterations (body-fitted / immersed) |
|---|---|---|---|---|
| Bolund 270 | rms speed-up error at 5 m, line B (4 masts) | 0.113 | 0.032 | 189 / 213 |
| Bolund 239 | rms speed-up error at 5 m, line A (5 masts) | 0.075 | 0.133 | 209 / 210 |
| Bolund 90 | rms speed-up error at 5 m, line B (4 masts) | 0.164 | 0.191 | 222 / 553 |
| Askervein 210 | rms speed-up error at 10 m, line B (16 masts) | 0.120 | 0.259 | 202 / 363 |
| Ridge, smooth | rms U error (198 points) | 0.272 m/s | 0.147 m/s | 537 / 1142 |
| Ridge, rough | rms U error (187 points) | 0.565 m/s | 0.528 m/s | 655 / 1126 |
| Butte 222 | rms speed error at 3 m (40 sensors) | 3.20 m/s | 6.14 m/s | 410 / 999 |

The immersed ground is closer on Bolund 270 and the smooth ridge, about the same on the rough ridge and Bolund 90,
and further on Bolund 239, Askervein and the Butte. On Askervein and the Butte the hill tops are a few cut columns
wide, and the immersed run is slower on the crests. The immersed runs take 1.0 to 2.5 times the iterations.

With the Bolund roughness map (`../WindNinja2019/bolund_z0_map.py`: 0.0003 m over water, 0.015 m over land) in
place of the uniform 0.01 m, and the inflow of the water (`bc.surface_z0_file=<data>/bolund_z0.xyz`; in the deck, `bc.surface_z0` removed, `bc.inflow_z0 = 0.0003` and `bc.inflow_ustar` =
0.4597 for 270, 0.3669 for 239), the same measure is 0.088 (body-fitted, 157 iterations) and 0.064 (immersed, 360)
on Bolund 270, and 0.079 (174) and 0.129 (188) on Bolund 239. The hybrid mesh at its defaults
(`terrain.representation = hybrid`) gives 0.061 (156) and 0.126 (172).

On the Butte, side walls on the ground's steps (a form since removed from the code) changed the speed at 3 m by
0.14 m/s rms (1082 iterations, rms error 6.09 m/s at the sensors). A wall distance from the column centre, also
removed, was tested on the 2-D ridge only (docs, Immersed terrain); its runs on the measured cases were stopped
before they finished.

The rows above were run before wall cells lost their tolerance (docs, Immersed terrain, Geometry). With the current
rule, Askervein (2 m layers) and the two ridges have the same wall cells, and a rerun of Askervein changes `u` by at
most 0.002 m/s (364 iterations). Bolund gains 29 wall cells (34199 to 34228) and the Butte 1058 (25759 to 26817);
those rows have not been rerun.

## CTest
| Name | Label | What it checks |
|---|---|---|
| `rt_immersed_askervein_coarse` | nightly | Askervein at 32 x 32 columns (188 m), both runs on 4 ranks with each deck's layers: the steady verdict, and the speed at 10 m above the ground of the two runs: rms difference over the columns more than 1000 m from the side faces below 22 %, difference on the highest ground below 52 % (of the mean body-fitted speed). |

Measured: the body-fitted run is steady after 128 iterations, the immersed run after 401. The mean body-fitted speed is 8.40 m/s over 400 columns; the rms difference is 20.1 %. On the highest ground (111 m) the body-fitted speed is 12.80 m/s and the immersed speed 8.71 m/s (-48.7 %): at 188 m columns the hill top is a few cut columns wide. The limits detect a change; they are not an accuracy claim. The pair takes 4.5 minutes on four ranks of an idle laptop and up to 50 minutes on a loaded one.

Mutants. With the cut layer's height taken at the layer's centre (`Source/Heights_K.H`, the centroid branch
never taken) the test fails: the immersed run stops in its first pressure solve (no convergence in 200
iterations). Two other changes stay inside the limits, and other tests cover them: with the solid inflow faces
not set to 0 the run passes (the open-faces CTests fail), and with the cut cell's wall-law shear not applied
(`LogShear::apply`) the rms difference is 20.3 % and the highest ground -49.1 %, 0.2 and 0.4 points from the
measured values (the unit tests `LogShear.ImmersedGroundIsExactFromTheCutCellUp` and
`ImmersedTerrain.StrainCorrectionFollowsTheFlatGround` fail).
