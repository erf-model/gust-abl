# Gold files

| Gold | What | How it was made |
|---|---|---|
| `erf_hill2d_btf/`, `erf_hill2d_stf/`, `erf_hill2d_sullivan/` | ERF's `z_phys`, `detJ` and nodal heights (`Level_0/Nu_nd`, component 2 = z - k dz) for `Mesh.MetricsEqualErf` | `tools/make_erf_mesh_gold.sh` with ERF at the commit in `external/VERSIONS` |

These golds are a cross-code check: they come from ERF, not from GUST. The regression golds written by
GUST live in each case folder (`Tests/<group>/<case>/Gold/`), with a "Gold files" table in that case's
README. A gold changes only in a reviewed commit that says why.

## Data and acknowledgements
The three golds are output of the ERF code (https://github.com/erf-model/ERF; BSD 3-clause with the LBNL clause,
`extern_erf/LICENSE.md`) at the commit in `external/VERSIONS`. Thanks to the ERF developers. They are redistributed
in this repository; nothing has to be downloaded.
