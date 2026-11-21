# Columbia_Gorge

## Purpose
An example of a real-terrain case set up from Python on steep ground: the Columbia River Gorge near Wasco, Oregon.
The river lies at about 50 m above sea level, the Columbia Hills to the north rise to about 950 m, and the plateau
to the south is at about 300 m. The domain and the wind are those of `Tests/Real_Terrain/Columbia_Gorge`, which
adds a nest and compares one-way and two-way nesting; this folder is the single-domain case alone. It is not a
comparison with observations, and none are used.

## Files
| File | Role |
|---|---|
| `case.yaml` | The case for `tools/make_gust_case.py`: 19.2 km x 16 km centred on 45.70 N, 120.80 W, top 3000 m, 200 m cells, first layer 10 m; wind 8 m/s at 10 m from the west (270 degrees), z0 = 0.05 m, neutral; a 1000 m flat border and a 2000 m slope to the terrain on every side. |

`make_gust_case.py` writes the folder `columbia_gorge`: the deck `inputs` (96 x 80 x 37 = 0.28 million cells, layers
from 10 m to 272 m, three hierarchy levels with grid sequencing, the recommended lee setting), `terrain.xyz`,
`run.sh`, `case_resolved.yaml` and `case_summary.txt`.

## Data and acknowledgements
Elevation: two SRTM-format tiles (N45W122, N45W121) of the AWS Terrain Tiles, downloaded by
`tools/fetch_terrain.py`; see `../README.md` for the providers and our thanks. Nothing is redistributed here.

## How to run
In a copy of this folder:
```bash
python3 <repo>/tools/fetch_terrain.py case.yaml
python3 <repo>/tools/make_gust_case.py case.yaml
cd columbia_gorge && GUST=<build>/gust ./run.sh
```

## Expected Results
- `fetch_terrain.py` writes `dem.asc` (1501 x 1049 points, 38 to 956 m) and `dem.sources.txt`.
- `make_gust_case.py` exits with 0 and the deck checker reports no errors or warnings. The summary gives the log
  law's u* = 0.6185 m/s and terrain from 0 to 906 m after 49 m is subtracted.
- GUST starts from the deck and writes plotfiles. Checked with two iterations on one rank. On this terrain the
  run may need `turb.nut_relax_end = 0.01` to reach the steady verdict
  (`Tests/Real_Terrain/Columbia_Gorge/README.md`); add it under `gust_inputs` in `case.yaml` if the residual
  stalls.

## CTest
None: the case needs a download. `Tests/Infrastructure/Case_Setup` tests the tool itself.
