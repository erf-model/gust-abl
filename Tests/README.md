# Tests

## Overview
The GUST tests. Each case has its own folder with a `README.md`, its input decks, its check script and its
gold files (`Gold/`). The layout follows ERF's `Exec/CanonicalTests/<group>/<case>/`. `CMakeLists.txt` in this
folder registers every case with CTest; each test runs in a directory under `<build>/Tests/run/` (its own by default; the
tests of a chain share one), so nothing is written into the source tree.

## Groups
| Folder | Cases | Character |
|---|---|---|
| `Infrastructure` | Smoke run and restart, input validation, the named set-ups (`Abl_Setup`), the canopy drag field, the building field, reference archives, YAML case set-up (`Case_Setup`), the wind library (`Library_Hill`, `Library_Ensemble`), the Kestrel scripts. | unit / regression / smoke |
| `Discretization` | The terrain mesh and its grid hierarchy, boundary fluxes, the k-L closure, the Poisson layer and the wall distance, each on a fixed state. | regression |
| `Initialization` | The mass-consistent initial field: Askervein, the hilltop speed-up trend against FastWindTerrain, and the first guess of a nested patch (`Nest_Patch_Hill`). | regression |
| `Canonical_RANS` | The SIMPLEC stage to steady on ERF's canonical RANS cases: the neutral column (also with the ground immersed, `Immersed_Column`; a ridge and a round hill immersed, `Immersed_Ridge_2D`; flat immersed ground between an inflow and an outflow face, `Immersed_Open_Faces`), Neutral_ABL_Flat, Neutral_Hill_2D, Neutral_Hill_3D and its truly neutral variant, and the smooth plateau comparison case. Also the thin-cell aspect-ratio sweep (`Aspect_Ratio_Hill_2D`) steep terrain at an outflow face (`Outflow_Steep_Terrain`), the forest canopy's drag and k terms (`Canopy_Drag`) a dense canopy on steep terrain with an inflow face (`Forested_Escarpment`), the canopy validation (`Homogeneous_Canopy`, `Forest_Edge`, `Forested_Hill`) and the same decks over an immersed ground and on a hybrid mesh (`Immersed_Canopy`), one block with the building terms on flat ground and on a ridge (`Building_Cube`, `Building_Hill`), a courtyard (`Building_Courtyard`), and the wind-tunnel cube of Snyder and Lawson (1994) against its measurements (`Snyder_Lawson_Cube`). | regression / nightly |
| `Real_Terrain` | Askervein to steady, with a log-law inflow and in a conventionally neutral boundary layer. | regression |
| `Real_Terrain_Python` | Real-terrain cases set up from a YAML file with `tools/fetch_terrain.py` and `tools/make_gust_case.py`: the AWAKEN area and the Columbia River Gorge. | examples (no CTest) |
| `Paper` | `JWEIA2026`: the driver and figure scripts of the paper on the realizability sensor, over the `Real_Terrain` comparison cases. | scripts (no CTest) |
| `Unit_Tests` | The GoogleTest suite (`unit_tests`) and its ERF golds. | unit |
| `Supporting_Files` | Check scripts shared by more than one case, and the test-only drivers (`Tools/`). | assets |

## Running
```bash
ctest --test-dir <build> -L "unit|smoke"
```
`-L regression` adds the longer runs; `-L nightly` runs the long cases (the 3-D hills, Askervein on ERF's mesh, the
nesting studies, the long immersed pairs). `-R <name>` runs one test; each case
README lists its CTest names.

## Notes
- Gold files change only in a reviewed commit that says why. The ERF golds (`Unit_Tests/Gold/erf_hill2d_*`) come
  from ERF; the others are change detectors written by GUST itself.
- `Gold/` folders are exempt from the whitespace and spelling checks (`.gitattributes`, `.codespellrc`).
- The full list of tests and what each checks is in `docs/testing.rst`.
- `inputs_reference` is the master input reference: every input with its default, range and meaning, as comments
  to copy into a deck. It is generated from `docs/inputs.rst` (`tools/input_reference.py --write`) and checked by
  the CTest `input_reference_current`.
- `inputs_yaml_reference` is the same for the YAML case description of `tools/make_gust_case.py`: every key with
  its default, unit, range and meaning, as comments to copy into a case file. It is generated from the tool's key
  table (`tools/make_gust_case.py --write-yaml-reference Tests/inputs_yaml_reference`) and checked by the CTest
  `gust_case_selftest`.
