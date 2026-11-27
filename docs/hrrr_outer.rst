HRRR as the outer level
=======================

``tools/hrrr_outer.py`` drives a GUST nest with one hour of the High-Resolution Rapid Refresh (HRRR) forecast. It
downloads a crop of the hour, writes the decks of a chain of three GUST levels from one case file, puts HRRR's
fields into the coarsest level's checkpoint, and scores the nest against HRRR. GUST itself reads no HRRR file and
has no HRRR input: the chain uses the nest faces (:doc:`boundary_conditions`) and ``init.first_guess = nest``.

The script needs numpy, and the ``case`` command PyYAML. Downloading an hour also needs
`Herbie <https://github.com/blaylockbk/Herbie>`_, xarray and cfgrib. The case folder
``Tests/Real_Terrain_Python/HRRR_Outer`` runs it at the AWAKEN site; :doc:`hrrr_regimes` classifies the hours.

Commands
--------

.. code-block:: bash

   # the crop of one hour: cycle and forecast hour (1 or later), a square of half-width 32 km
   python3 tools/hrrr_outer.py fetch --cycle 2023-05-10T19:00Z --fxx 1 --lat 36.362 --lon -97.405 \
       --half-km 32 --save grib --out awaken_2023051019_f01.npz
   # the decks of the chain and its run.sh, from a gust_case file with an hrrr section
   python3 tools/hrrr_outer.py case case.yaml
   cd day20 && GUST=<build>/gust ./run.sh
   # the steps run.sh calls
   python3 tools/hrrr_outer.py to-chk CROP.npz hrrr_raw/chk00000 --lat 36.362 --lon -97.405
   python3 tools/hrrr_outer.py balance hrrr_raw/chk00000
   python3 tools/hrrr_outer.py score CROP.npz nest/plt00247 --lat 36.362 --lon -97.405 --log nest/log.txt
   # HRRR's box-mean theta_v profile, as written to the nest's deck
   python3 tools/hrrr_outer.py profile CROP.npz --lat 36.362 --lon -97.405 --half-width 12500

``fetch`` downloads byte-range subsets of the native-level (``wrfnat``) and surface (``wrfsfc``) files: wind,
vertical velocity, height, temperature, humidity, pressure, cloud water and k on the first 30 hybrid levels, and the
ground height, :math:`z_0`, the sensible heat flux, the 2 m temperature and the surface pressure. The crop is
:math:`\sqrt{2}` times ``--half-km`` plus three HRRR cells a side, so a square of that half-width fits at any turn of
HRRR's grid against the UTM grid. It records the cycle and the forecast hour; ``to-chk``, ``profile``, ``score`` and
``case`` read them from it. ``to-chk`` refuses a crop without them. A refusal names the input or the field and
exits with status 1 (status 2 for a malformed command line). ``to-chk`` checks the top, the cloud and every field
before it writes any of them, so a refused hour leaves the checkpoint as it was.

The chain
---------

``case`` writes three GUST cases beside each other and a ``run.sh`` that runs them in turn:

1. ``hrrr_raw``: a coarse mesh larger than the parent that only writes its first guess (no steady stage). Its
   checkpoint is the skeleton: ``to-chk`` overwrites u, v, w, :math:`\theta` and k with HRRR's at every point
   inside the domain (ghost points outside keep the first guess; p and the heights are kept), with w = 0 on the
   ground faces.
2. ``parent``: ``init.first_guess = nest`` from the raw checkpoint, which applies GUST's mass-consistent
   adjustment on the parent's terrain; first guess only.
3. ``nest``: the steady run, with the parent's field on its four lateral faces (``bc.nest_ratio``).

``run.sh`` then scores the nest's last plotfile into ``score.txt`` and ``score.json``. The meshes follow from the
nest's half-width :math:`L` and cell :math:`\Delta`. The nest has :math:`n = 2L/\Delta` cells a side, a multiple of
8. The parent's cells are ``parent_ratio`` :math:`r` times larger. Its half-width is the smallest multiple of four
parent cells at least :math:`2L`. The raw domain's cells are twice the parent's, and its half-width is the smallest
multiple of four raw cells at least the parent's plus ``raw_margin``. Every face of a finer mesh then lies on a face
of the coarser one. The terrain blends into the coarser level's over two of its cells (``terrain.nest_blend_width``).
All three levels share the vertical mesh and the terrain file. For the AWAKEN case (:math:`L` = 12500 m,
:math:`\Delta` = 125 m, :math:`r` = 4):

.. list-table::
   :header-rows: 1

   * - Level
     - Half-width (m)
     - Cell (m)
     - Cells a side
   * - ``hrrr_raw``
     - 28000
     - 1000
     - 56
   * - ``parent``
     - 26000
     - 500
     - 104
   * - ``nest``
     - 12500
     - 125
     - 200

The scheme is fixed, and the keys it writes cannot be given in the case's ``gust_inputs``:

- **A forecast hour of 1 or later.** The f00 analysis carries almost no vertical motion; it is refused by name.
- **A domain top of at least 3000 m.** HRRR's horizontal convergence is balanced near 3 km, not below; a lower
  slip-wall top forces the adjustment to change every face velocity. The top must also be at most 3600 m (the
  height table) and below HRRR's highest fetched level in every column of the crop.
- **The adjusted parent.** The nest is driven by the mass-consistent parent, not by HRRR's raw field.
- **The virtual potential temperature** :math:`\theta_v` as GUST's :math:`\theta`. ``bc.theta_profile`` is HRRR's
  box mean over the nest, ``forcing.theta0`` its ground value and ``bc.zhi_theta_grad`` its gradient over the top
  100 m; the Rayleigh layer acts on w and :math:`\theta`.
- **Coriolis on** at the centre's latitude.
- **The length-scale cap at 100 m** (``turb.max_geom_lscale``).
- **A cloudy hour refused.** Cloud water above ``cloud_mixing_ratio`` in any HRRR column within 3 km of the raw
  domain, below its top, refuses the hour, naming the amount and its height: GUST is dry.

The ``hrrr`` section
--------------------

The case file is a gust_case description of the nest (:doc:`case_setup`) with one more section, ``hrrr``. The
sections ``wind``, ``atmosphere`` and ``nests`` are refused (the chain sets them), as are a domain given by its
edges, a ``domain.utm_zone`` other than the centre's (HRRR is placed in that zone), ``output.prefix`` (``run.sh``
scores ``pltNNNNN``), ``mesh.hierarchy_levels`` other than 1, ``solver.grid_sequencing``, and
``model.length_scale_cap`` other than ``auto`` or 100. The case's ``gust_inputs`` reach the nest's deck only. An unknown key in ``hrrr`` is refused by name; every value is range-checked and never clamped.

.. list-table::
   :header-rows: 1
   :widths: 22 14 64

   * - Input
     - Default
     - Effect
   * - ``hrrr.file``
     - (required)
     - The crop written by ``fetch``, relative to the case file. It carries the cycle and the forecast hour.
   * - ``hrrr.parent_ratio``
     - 4
     - ``bc.nest_ratio`` of the nest: 2, 4, 8 or 16. The nest's faces must lie on parent faces
       (:math:`n/2` a multiple of :math:`r`).
   * - ``hrrr.raw_margin``
     - 2000
     - How far the raw domain reaches beyond the parent, m; > 0.
   * - ``hrrr.theta``
     - ``solved``
     - ``solved``: :math:`\theta` is solved with gravity on. ``frozen``: ``simplec.theta = frozen``.
   * - ``hrrr.surface``
     - ``neutral``
     - ``heat_flux``: ``bc.surface_stability = heat_flux`` with HRRR's box-mean virtual kinematic heat flux
       (below).
   * - ``hrrr.geostrophic``
     - ``false``
     - ``true``: ``forcing.geostrophic`` with HRRR's geostrophic wind (below).
   * - ``hrrr.cloud_mixing_ratio``
     - 1e-5
     - Cloud water (kg/kg) above which an hour is refused; > 0. The classifier's cloud threshold. Write
       ``1.0e-5``: YAML reads ``1e-5`` as text, which is refused.

The raw and parent levels get only the keys GUST reads without a steady stage: ``case`` removes the ``forcing.``,
``convection.``, ``simplec.``, ``turb.`` and ``steady.`` keys, ``io.turb_plot_file``, ``abl.drive`` and
``abl.theta`` from their decks. Their inflow profile is the log law through HRRR's box-mean wind at 80 m; it is
overwritten by ``to-chk`` and by the nest faces.

HRRR's fields at GUST's points
------------------------------

Winds are turned from HRRR's grid (Lambert conformal, true latitude 38.5° N, central longitude 97.5° W) to east and
north, and then to the UTM grid of the domain's centre:

.. math::

   \alpha = \sin(38.5^\circ)\,(\lambda + 97.5^\circ), \qquad
   c = \arctan\!\left(\tan(\lambda - \lambda_0)\,\sin\varphi\right),

where :math:`\lambda` and :math:`\varphi` are the column's longitude and latitude, :math:`\lambda_0 = 6z - 183^\circ`
is the central meridian of UTM zone :math:`z`, :math:`\alpha` turns the grid wind to east and north and :math:`c` (the
meridian convergence) turns that to the UTM grid. At the AWAKEN site the two together are about 1°. At 104.5° W and
40° N they are 4.7°, where leaving them out would put 8 % of the speed across the wind.

The vertical velocity is :math:`w = -\omega R_d T_v / (p\, g)`, with :math:`\omega` HRRR's vertical velocity in
pressure, :math:`T_v = T(1 + 0.61 q)`, :math:`R_d` = 287.04 J/(kg K) and :math:`g` = 9.81 m/s². The potential
temperature is :math:`\theta_v = T (10^5/p)^{0.2857} (1 + 0.61 q)`. Each column goes on a table of heights above
HRRR's own ground: 0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1 and 1.5 m, where the log law curves, then every 2 m to
3600 m. Below HRRR's first level (about 11 m at AWAKEN), u and v follow the neutral log law with HRRR's :math:`z_0`
at the first level's direction, w goes linearly to 0, and :math:`\theta_v` and k keep the first level's values. A GUST point takes the table at its height above GUST's ground: the mean of the point's surrounding nodes
less the mean of their ground nodes. Horizontally the table is interpolated bilinearly in HRRR's grid indices, with
the map from indices to UTM inverted by Newton's method.

The geostrophic wind is a least-squares plane :math:`p = a x + b y + c` at each height
:math:`h` = 200, 300, …, 1000 m above the box-mean ground, fitted over the HRRR columns inside the nest, with
:math:`x, y` the UTM coordinates:

.. math::

   u_g = -\frac{b}{\rho f}, \qquad v_g = \frac{a}{\rho f}, \qquad \rho = \frac{p}{R_d T_v},
   \qquad f = 2\Omega\sin\varphi_0,

averaged over the nine heights. The heights are above the box mean, so in a column whose first level lies above
200 m over that mean the fit would use extrapolated pressures; such a box is refused. The fit varies with height:
at AWAKEN its speed ranges over 7.4–12.8 m/s between 200 and 1000 m in the day hour and 3.0–6.7 m/s in the night
hour (fits at single heights, 100 m apart). The mean over that layer is a choice, not a property of HRRR. The surface heat flux is the virtual kinematic flux:

.. math::

   \overline{w'\theta_v'} = \frac{H}{\rho c_p}(1 + 0.61 q_2) + 0.61\,\theta\,\frac{L_E}{\rho L_v},
   \qquad \rho = \frac{p_s}{R_d T_2 (1 + 0.61 q_2)}, \qquad \theta = T_2 (10^5/p_s)^{0.2857},

with :math:`H` and :math:`L_E` HRRR's sensible and latent heat fluxes (W/m², upward positive), :math:`T_2` and
:math:`q_2` the 2 m temperature and specific humidity, :math:`p_s` the surface pressure,
:math:`c_p` = 1004 J/(kg K) and :math:`L_v` = 2.5e6 J/kg, as a box mean. GUST's :math:`\theta` is :math:`\theta_v`,
so the moist part belongs in its flux. The November 2026 study used :math:`H/(\rho c_p)` alone.

The score
---------

``score`` compares the nest with HRRR on the nest's cells more than ``--margin`` (default 2500 m) from the lateral
faces. With :math:`S` the nest's horizontal speed, :math:`S_H` HRRR's at the same point, and :math:`\mathbf{u}`,
:math:`\mathbf{u}_H` the horizontal vectors:

.. math::

   \text{bias} = \overline{S - S_H}, \qquad \text{rms} = \sqrt{\overline{(S - S_H)^2}}, \qquad
   \text{vector rms} = \sqrt{\overline{|\mathbf{u} - \mathbf{u}_H|^2}},

over all cells, and over the cells below 500 m above the ground for the low-level bias and vector rms. The low-level
direction difference is between the directions of the two mean vectors below 500 m. The table also gives w's rms
and largest magnitude (HRRR's rms beside it), the ratio of the summed k below 500 m to HRRR's, a line per level
(mean height, speeds, bias, rms, vector rms, directions, k), and from the nest's log the iterations and whether the
run was steady (the last "steady", "not steady" or "diverged" line). ``--json`` writes the same as JSON.

Three things about the numbers. The means are over cells, not volumes: GUST's layers grow upward, so the cells
near the ground weigh more than their volume (volume weighting would move the day hour's bias from −0.03 to
−0.12 m/s). The levels below HRRR's first level (about 11 m) compare the nest with the log-law fill above, not
with a HRRR value. Directions are from north on the UTM grid, which at AWAKEN is 0.95° from true north (the same
for the nest and HRRR, so the difference is unaffected).

Results
-------

The day hour at the AWAKEN site (2023-05-10 20 UTC, cycle 19 UTC, forecast hour 1), by the case's ``run.sh``,
against the run of November 2026 with the scripts this tool replaces (same meshes, :math:`\theta` solved, Coriolis
on). The differences come from the horizontal interpolation (bilinear in HRRR's grid rather than over triangles: up
to 0.11 m/s in the first guess) and the re-written checkpoint. The largest :math:`|w|` is over the cells inside the margin;
the November number is over the whole nest.

.. list-table::
   :header-rows: 1

   * - Run
     - Bias (m/s)
     - Rms (m/s)
     - Bias < 500 m
     - Vector rms < 500 m
     - Direction < 500 m
     - w rms / max (m/s)
     - k ratio
     - Iterations
   * - November 2026
     - 0.02
     - 0.47
     - 0.00
     - 0.57
     - 0.3°
     - 0.071 / 1.20
     - 0.17
     - 239, steady
   * - **This tool**
     - −0.03
     - **0.48**
     - −0.04
     - **0.58**
     - 0.3°
     - 0.085 / 0.79
     - 0.17
     - 247, steady

The variant matrix
~~~~~~~~~~~~~~~~~~

``Tests/Real_Terrain_Python/HRRR_Outer/variants.py`` runs the chain of each hour and six variants of its nest on
the chain's parent (2026-11-23, 4 ranks): **S** is the scheme as ``case`` writes it (:math:`\theta` solved,
Coriolis on, neutral ground); A to F are the variants of the November 2026 study (A :math:`\theta` frozen and
Coriolis off; B frozen; C B with HRRR's geostrophic wind; D S with it; E C with an 8-cell relaxation zone; F C with
HRRR's surface heat flux). The low-level jet hour (2023-05-10 05 UTC) is refused as cloudy, so there are four
hours. Every run that started was steady. Bias, rms and vector rms are in m/s; "< 500 m" is below 500 m above the
ground.

.. csv-table::
   :header: "Hour", "Variant", "Bias", "Rms", "Bias < 500 m", "Vector rms < 500 m", "Direction < 500 m (°)", "w rms / max", "k ratio", "Iterations"

   "day20", "**S**", "−0.03", "**0.48**", "−0.04", "**0.58**", "0.3", "0.085 / 0.79", "0.17", "247"
   "day20", "A", "−0.24", "0.39", "−0.42", "0.64", "3.1", "0.043 / 0.43", "0.35", "164"
   "day20", "B = C", "−0.21", "0.35", "−0.36", "0.48", "1.0", "0.040 / 0.43", "0.35", "178"
   "day20", "E", "−0.21", "0.35", "−0.36", "0.48", "0.9", "0.041 / 0.43", "0.35", "160"
   "day20", "F", "−0.25", "0.42", "−0.45", "0.58", "0.9", "0.039 / 0.40", "0.37", "164"
   "eve01", "**S**", "−0.58", "**1.22**", "−0.94", "**1.77**", "4.6", "0.144 / 0.74", "0.79", "336"
   "eve01", "A", "−0.23", "0.69", "−0.33", "1.05", "3.0", "0.118 / 0.87", "0.81", "198"
   "eve01", "B = C", "−0.41", "0.90", "−0.65", "1.19", "2.0", "0.112 / 0.75", "0.78", "222"
   "eve01", "E", "−0.41", "0.89", "−0.65", "1.19", "2.0", "0.113 / 0.75", "0.78", "220"
   "eve01", "F", "−0.40", "0.89", "−0.63", "1.18", "2.0", "0.112 / 0.75", "0.77", "221"
   "night09", "**S**", "−0.15", "**0.53**", "−0.29", "**0.92**", "0.3", "0.075 / 0.66", "1.80", "334"
   "night09", "A", "−0.38", "0.82", "−0.74", "1.16", "0.9", "0.168 / 1.20", "0.76", "302"
   "night09", "B = C", "−0.55", "1.04", "−1.02", "1.42", "0.5", "0.189 / 1.48", "0.68", "343"
   "night09", "E", "−0.54", "1.03", "−1.01", "1.41", "0.5", "0.186 / 1.50", "0.69", "340"
   "night09", "F", "stopped", "", "", "", "", "", "", ""
   "dawn12", "**S**", "0.11", "**0.75**", "0.47", "**1.12**", "1.9", "0.097 / 0.62", "2.59", "355"
   "dawn12", "A", "0.14", "1.42", "0.61", "3.04", "8.2", "0.544 / 3.26", "2.82", "1072"
   "dawn12", "B = C", "0.09", "1.49", "0.62", "3.01", "7.3", "0.587 / 3.41", "2.98", "1065"
   "dawn12", "E", "0.09", "1.46", "0.60", "2.95", "7.2", "0.571 / 3.29", "2.89", "1069"
   "dawn12", "F", "stopped", "", "", "", "", "", "", ""

- **C equals B and D equals S** in every hour (the iterations differ by up to 64). With all four lateral faces
  given, a height-constant geostrophic force is the gradient of a linear pressure field, which the pressure solve
  takes up whole; the velocity does not change. ``hrrr.geostrophic`` therefore has no effect on this chain; it is
  kept for a nest with an open face. November's B and C scores were equal in the same way.
- **E** (the relaxation zone) moves no speed score by more than 0.06 m/s from C.
- **F stops at the start in the two stable hours**: with HRRR's downward heat flux (−0.009 and −0.004 K m/s),
  GUST's surface layer finds no finite friction velocity in some low-wind columns and aborts. November's night F
  run stopped the same way.
- **S** has the lowest rms in the night and dawn hours and the highest in the day and evening hours. With :math:`\theta` frozen
  the dawn hour needs about 1070 iterations and puts 3 m/s of vector rms below 500 m; with :math:`\theta` solved
  it needs 355 and 1.12 m/s.

November's nests of the day and night hours, scored by this tool against the same crops (the other hours' November
plotfiles are not kept), differ from this matrix by up to 0.04 m/s in rms (day B 0.32 against 0.35, day F 0.38
against 0.42, which also has the virtual heat flux, night D 0.51 against 0.53).

Not yet supported
-----------------

- The same score against a lidar profile.
- HRRR hours with cloud; levels above 3600 m.
- HRRR's surface heat flux in a stable hour (variant F; see above).
