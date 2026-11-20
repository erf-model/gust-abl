# AWAKEN

## Purpose
An example of a real-terrain case set up from Python on nearly flat ground: the area of the AWAKEN field campaign
(American WAKE experimeNt) north of Enid, Oklahoma. The domain is 12.8 km square, centred on 36.362 N, 97.405 W
(an approximate centre of the area, not a surveyed instrument site). The case shows the workflow. It is not a
comparison with the AWAKEN observations, and none are used. There are no wind turbines in the run.

## Files
| File | Role |
|---|---|
| `case.yaml` | The case for `tools/make_gust_case.py`: half-width 6400 m, top 2000 m, 100 m cells, first layer 5 m; wind 8 m/s at 80 m from the south (180 degrees), z0 = 0.05 m, neutral; a 500 m flat border and a 1000 m slope to the terrain on every side. |

`make_gust_case.py` writes the folder `awaken`: the deck `inputs` (128 x 128 x 39 = 0.64 million cells, layers from
5 m to 186 m, three hierarchy levels with grid sequencing, the recommended lee setting), `terrain.xyz`, `run.sh`,
`case_resolved.yaml` and `case_summary.txt`.

## Data and acknowledgements
Elevation: one SRTM-format tile (N36W098) of the AWS Terrain Tiles, downloaded by `tools/fetch_terrain.py`; see
`../README.md` for the providers and our thanks. Nothing is redistributed here.

## How to run
In a copy of this folder:
```bash
python3 <repo>/tools/fetch_terrain.py case.yaml
python3 <repo>/tools/make_gust_case.py case.yaml
cd awaken && GUST=<build>/gust ./run.sh
```

## Expected Results
- `fetch_terrain.py` writes `dem.asc` (853 x 687 points, 274 to 372 m) and `dem.sources.txt`.
- `make_gust_case.py` exits with 0 and the deck checker reports no errors or warnings. The summary gives the log
  law's u* = 0.4445 m/s and terrain from 0 to 59 m after 290 m is subtracted.
- GUST starts from the deck and writes plotfiles. Checked with two iterations on one rank; the run to the steady
  verdict has not been made.

## CTest
None: the case needs a download. `Tests/Infrastructure/Case_Setup` tests the tool itself.
