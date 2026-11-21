# Real_Terrain

## Purpose
Askervein hill run to steady on ERF's grid (`extern_erf/askervein.txt`). The folder name follows ERF's
`Exec/CanonicalTests/Real_Terrain/Askervein`.

## Cases
| Folder | CTest | Description |
|---|---|---|
| `Askervein` | `rt_steady_askervein` | Log-law inflow, outflow and slip sides, no Coriolis, to the steady verdict. |
| `Askervein_210` | `rt_askervein_210_erf_lines`, `rt_askervein_210_level_history` | ERF's own setup (210° power-law inflow, ERF's full mesh) against ERF's URANS on lines A and AA. |
| `Askervein_Nest` | `rt_askervein_nest_study` | One-way nesting: 20 m patches around the hilltop inside the 40 m run, against the uniform 20 m run on lines A and AA. |
| `Askervein_CNBL` | `rt_precursor_askervein_cnbl`, `rt_steady_askervein_cnbl` | A conventionally neutral boundary layer: a precursor column writes the inflow table, then the 3-D case runs on it. |
| `Columbia_Gorge` | (comparison case) | Two-way nesting over the Columbia River Gorge at the WFIP2 site: the 200 m outer domain and a 100 m nest across the river and the Columbia Hills, one way and joint, against the uniform 100 m run. |
| `WindNinja2019` | (comparison case) | Askervein, Bolund cases 1, 3, 4 and Big Southern Butte against the observations and WindNinja curves of Wagenbrenner et al. (2019); GUST - Fast and GUST - CFD. |
| `Immersed_Terrain_Cases` | `rt_immersed_askervein_coarse` | Bolund, Askervein, the Ishihara-Zhou ridge and Big Southern Butte as pairs in the same settings: terrain-following mesh against the ground immersed in a flat mesh (comparison cases); the CTest runs a coarse Askervein pair. |
| `IshiharaZhou` | (comparison case) | The steep 2-D ridge of Ishihara and Zhou (2025), smooth and rough, against the wind tunnel. |
| `RUSHIL` | (comparison case) | The RUSHIL 2-D hills (a/H = 3, 5, 8) against the wind tunnel (ERCOFTAC AC5-05). |
| `URANS_Averages` | `rt_urans_mean_option` | The time mean of the unsteady run against the steady run of the same deck: RUSHIL Hill 5 and Askervein with ERF's setup (comparison cases); the CTest checks the `--mean` option of the readers. |

The first three use `Supporting_Files/expect_abort.py --success "gust: steady after"` to require the verdict.
