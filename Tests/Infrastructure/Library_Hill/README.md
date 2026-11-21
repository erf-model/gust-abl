# Library_Hill

## Purpose
The wind library (`tools/make_gust_library.py`, `tools/gust_library/`) on a sweep GUST ran: a library built from
two directions is blended at the direction midway and compared with the run from that direction. The case shows
how large the blend's error is for sectors 20 and 10 degrees wide over a hill, also away from the midway
direction and for k, and that a stored sector comes back unchanged. The library and its conventions are described in `docs/case_setup.rst`.

## Files
| File | Role |
|---|---|
| `library_hill.yaml` | A Gaussian hill (60 m high, 200 m wide) at the centre of a 2.4 km domain, 50 m cells, neutral log-law inflow of 8 m/s at 10 m over z0 = 0.03 m, five directions 260 to 280 degrees, 5 apart. |
| `check_library_hill.py` | Writes `hill.asc`, writes the case with `make_gust_case.py --strict`, runs `run.sh`, builds the libraries and checks them. |

## Data and acknowledgements
No downloaded data. The hill grid is written by `check_library_hill.py` from a formula.

## How to run
```bash
ctest --test-dir <build> -R "^rt_library_hill$"
```
By hand, from an empty directory:
```bash
python3 <repo>/Tests/Infrastructure/Library_Hill/check_library_hill.py --gust <build>/gust --mpi "mpiexec -n 2"
```
The five runs take about 65 s on 2 ranks (137 to 158 iterations each).

## Expected Results
The reference point is the library's default, the domain centre (the hilltop), and the reference height is 10 m.
- `make_gust_case.py --strict` and `run.sh` exit with 0, and each run reaches GUST's steady verdict.
- A library of 260 and 280 queried at 260 gives the 260 run's wind at 10 m, scaled to the asked 8 m/s, to 1e-12
  of the speed. Measured: 0.
- The blend at 270 has 8 m/s at the reference point and 10 m (1e-12 relative). Measured: 2e-16.
- The blend at 270 against the 270 run, at 10 m and 50 m: the largest difference of the horizontal wind over the
  reference speed, with the 4 cells next to each lateral face left out.

  | Library | Height | Interior max (tolerance) | Interior RMS | All cells max |
  |---|---|---|---|---|
  | 260/280 (20 degrees) | 10 m | 0.67 % (1 %) | 0.20 % | 1.55 % |
  | 260/280 (20 degrees) | 50 m | 0.44 % (1 %) | 0.11 % | 1.26 % |
  | 265/275 (10 degrees) | 10 m | 0.23 % (0.4 %) | 0.10 % | 1.21 % |
  | 265/275 (10 degrees) | 50 m | 0.19 % (0.4 %) | 0.09 % | 0.99 % |

  Measured in double precision on 2 ranks; single precision is not measured. The interior maximum falls by 2.9
  from 20 to 10 degree sectors. The 10 degree values are close to the runs' own convergence level (the 260 and 280
  runs mirror each other to 0.3 % at r = 1e-4), so the ratio is not a measured order of accuracy. The faces are
  left out because from 260 and 280 the y faces are inflow and outflow faces and from 270 they are slip walls; the
  largest differences sit in the outflow corners (1.2 to 1.6 %).
- The 260/280 library queried at 265 and 275 against those runs (tolerance 1 %): 0.51 % at both heights.
- k: the interior mean of k / k_270 - 1 for the blend at 270 (tolerance 0.5 %): +0.11 % (10 m) and +0.33 % (50 m)
  for 260/280, +0.07 % and +0.28 % for 265/275.
- The vertical wind is printed, not checked: blend differences up to 0.02 m/s where |w| reaches 1 m/s.
- Mutants: without the renormalisation to the asked speed (`scale = speed`) the wind checks fail (1.54 % and
  1.73 % for 260/280, 0.51 % and 0.56 % for 265/275, reference speed 7.879 and 7.969 m/s); with k scaled by the
  renormalised speed squared the k checks fail (+3.2 % and +3.4 % for 260/280, +0.8 % and +1.1 % for 265/275).

## CTest
`rt_library_hill` (label `regression`, 2 processors; not on Windows, where `run.sh` cannot run). It and
`gust_library_selftest` (the library's unit tests, `tools/test_gust_library.py`) are registered only when the
configure's Python has PyYAML and numpy.
