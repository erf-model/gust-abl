# Real_Terrain_Python

## Purpose
Real-terrain cases set up from Python: each folder holds one YAML case description, and the tools in `tools/`
write everything GUST reads from it. `tools/fetch_terrain.py` downloads the elevation data for the domain, and
`tools/make_gust_case.py` writes the deck, `terrain.xyz` and a run script (`docs/case_setup.rst`). No deck is kept
in the repository: the YAML file is the case.

The cases in `Tests/Real_Terrain` are the other way round: their decks are written by hand.

## Cases
| Folder | CTest | Description |
|---|---|---|
| `AWAKEN` | (example case) | The AWAKEN area north of Enid, Oklahoma: 12.8 km square, 100 m cells, gently rolling ground (59 m of relief), wind from the south. |
| `Columbia_Gorge` | (example case) | The Columbia River Gorge near Wasco, Oregon: 19.2 km x 16 km, 200 m cells, 906 m of relief, wind from the west. |

## How to run
In a copy of a case folder (the tools write beside the YAML file):
```bash
python3 <repo>/tools/fetch_terrain.py case.yaml
python3 <repo>/tools/make_gust_case.py case.yaml
cd <case name> && GUST=<build>/gust ./run.sh
```
The first command needs an internet connection the first time; the tiles are then kept in a cache
(`~/.cache/gust/terrain_tiles`).

## Data and acknowledgements
The elevation data are the SRTM-format tiles of the Terrain Tiles data set on the Registry of Open Data on AWS
(https://registry.opendata.aws/terrain-tiles/), made by Mapzen and maintained by the Tilezen project. We thank
Mapzen, Tilezen and Amazon Web Services for making them free to use, and NASA and the USGS for the SRTM and 3D
Elevation Program data the tiles hold over these two sites. `fetch_terrain.py` writes the attribution the data
need to `dem.sources.txt` beside `dem.asc`; a product made from the data must carry it (`docs/case_setup.rst`,
the terrain download). No elevation data are redistributed in this repository: `dem.asc` is downloaded by the user.
