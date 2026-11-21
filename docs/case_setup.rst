Case set-up from a YAML file
============================

``tools/make_gust_case.py`` writes a GUST case folder from a short YAML case description: the deck (or decks), the
terrain and roughness files GUST reads, a run script and a summary. It then checks every deck it wrote with the deck
checker, ``tools/check_gust_inputs.py``, which can also be run on any deck by hand. The workflow follows the YAML
case set-up of amrTerrain (a private case set-up tool written for AMR-Wind), with GUST's own
inputs, cleaner key names and the capabilities GUST has that amrTerrain does not set up (named drives, a column
first guess, nested patches, canopy and buildings).

.. contents:: On this page
   :local:
   :depth: 1

Requirements
------------

Python 3.8 or later with PyYAML and numpy. rasterio is optional: without it, GeoTIFF files go through a built-in
reader (little-endian, one band, strips with no compression or LZW, integer or float samples, ``GDAL_NODATA`` and
pixel-is-point); with it, any GeoTIFF that rasterio reads. Latitude and longitude go to UTM through a built-in
conversion (Karney's series, below 1 mm within a zone); pyproj is not needed. The deck checker needs only the
standard library.
The terrain download (``tools/fetch_terrain.py``) needs numpy and an internet connection. The case set-up window
(``tools/gust_gui.py``) also needs PySide6 and pyqtgraph (``pip install PySide6 pyqtgraph numpy pyyaml``).

A first case
------------

.. code-block:: bash

   python3 tools/make_gust_case.py --example > case.yaml
   python3 tools/make_gust_case.py case.yaml
   cd flat_west && GUST=../build-rel/gust ./run.sh

``--example`` prints a starter case: wind from the west, 8 m/s at 10 m, over 4 km of flat ground. Add a
``terrain.file`` and a ``domain.center`` to move it onto real terrain:

.. code-block:: yaml

   case:
     name: butte_222
   domain:
     center: [9345., 10080.]   # in the terrain file's frame (m)
     half_width: 4000.         # m from the centre to each face
     top: 2800.                # m above the lowest point kept
   terrain:
     file: dem.tif
   mesh:
     cell_size: 50.
     layers: 40
     first_layer: 4.
   wind:
     speed: 8.
     height: 10.
     direction: 222.           # the direction the wind comes from, clockwise from north
   surface:
     roughness: 0.01

The folder ``Tests/Infrastructure/Case_Setup`` holds five more examples, each run and checked by CTest: a sweep of
directions, a nested patch, a geostrophic wind from the column (along a diagonal, and along x and y) and a hill from
an ESRI ASCII grid.

What the tool does
------------------

- **Domain and mesh.** The domain is a box around ``domain.center``. The cell counts are the nearest multiples of 8
  to the asked cell size (so the hierarchy coarsens), and ``mesh.hierarchy_levels = auto`` takes up to three levels
  as the counts allow. With ``mesh.first_layer`` the layers stretch geometrically to reach ``domain.top``.
- **Terrain.** The elevation data are cut to the domain plus a margin, thinned to no finer than half a cell,
  shifted so that the lowest point kept is z = 0 (``terrain.zero_level``) and written as ``terrain.xyz``. A
  geographic file is projected to UTM. ``terrain.border`` blends the ground to a flat band inside the lateral faces.
- **Wind.** The meteorological direction becomes ``bc.inflow_direction`` (towards, counter-clockwise from x:
  :math:`270^\circ - \text{from}`). The faces the wind enters through are ``inflow``, those it leaves through are
  ``outflow``; a face parallel to the wind is ``wind.side_faces``. With a geostrophic drive the wind turns with
  height, so the faces are chosen from the middle of the turn: the wind asked for, turned by half the drag law's
  angle :math:`\alpha = \operatorname{atan2}(B,\ \ln(u_*/(|f| z_0)) - A)` towards the geostrophic wind (clockwise
  in the northern hemisphere; :math:`A`, :math:`B` and :math:`f` below). With ``atmosphere.geostrophic_wind`` given,
  the middle is that wind's direction turned back by :math:`\alpha/2`, :math:`u_*` from the drag law, and
  ``wind.direction`` does not set the faces. No face is then a slip wall: a face parallel to the wind near the ground
  takes the flow that crosses it aloft. The faces also differ from those of the wind asked for when it is within
  :math:`\alpha/2` of an axis (11.8 degrees for 8 m/s at 10 m over :math:`z_0 = 0.03` m at 45 degrees north), where
  the wind crosses a face one way near the ground and the other way aloft. ``wind.side_faces`` has no effect and is
  refused. :math:`\alpha/2` is a rule of thumb, not the direction of the mass flux: over a 1.5 to 3 km top at 45
  degrees north the height-integrated flux turns 0.8 to 0.9 :math:`\alpha`, and at 10 degrees north, where the
  boundary layer is deeper than the domain, 6.8 degrees against :math:`\alpha/2 = 11`. A face nearly parallel to the
  wind carries little of the flow either way: from 265 degrees the 3-D run is 2.03 % off the column with the turned
  faces and 1.93 % with those of the wind asked for. Below about 15 degrees of latitude, check the faces against the
  column's table. From 270 degrees in that case, slip walls left the 3-D run up to 72 % of the top wind off the column's;
  the turned faces leave 2.0 % (``rt_case_setup_column_axis``). The log law gives
  :math:`u_* = \kappa U_\mathrm{ref} / \ln((z_\mathrm{ref} + z_0)/z_0)`, :math:`\kappa = 0.41`.
- **Drive** (``abl.drive``, :doc:`howto`). ``inflow_only`` and ``couette`` (a moving lid at the log law's speed at
  the top) need no rotation. ``geostrophic_column`` writes a 1-D column run first: with ``atmosphere.geostrophic_wind``
  not given, the column turns and scales the geostrophic wind until the wind at ``wind.height`` is ``wind.speed``
  from ``wind.direction``. It starts from the neutral drag law,

  .. math::

     G = \frac{u_*}{\kappa}\sqrt{\left(\ln\frac{u_*}{|f|\,z_0} - A\right)^2 + B^2},
     \qquad A = 1.8,\quad B = 4.5,

  with :math:`f = 2\Omega\sin\phi` the Coriolis parameter at latitude :math:`\phi`. ``run.sh`` reads the adjusted
  ``forcing.ug`` and ``forcing.vg`` from the column's log and gives them to the 3-D run, which starts from the
  column (``init.first_guess = column``). ``geostrophic_inflow`` uses the same column for the inflow faces.
- **Potential temperature** (``abl.theta``): ``neutral``, ``capped`` (an inversion and lapse rate, solved),
  ``frozen`` (a sounding held fixed) or ``surface_flux``.
- **Roughness.** One ``surface.roughness``, or a map: ``x y z0`` points or a land-cover raster with a class table
  (ESA CCI by default, or your own mapping). The ESA CCI table is amrTerrain's with five classes moved to the
  updated Davenport classes of Wieringa (1992, *J. Wind Eng. Ind. Aerodyn.* 41, 357–368): no data (0) is refused
  or takes ``surface.roughness_map.default``; tree or shrub cover in cropland (12) is 0.25 m; flooded tree cover
  (160, 170) is 1.0 m, as the other trees; bare areas (200–202) are 0.005 m. The class codes are those of the ESA
  Climate Change Initiative Land Cover maps; thanks to ESA and the Land Cover CCI team for making them public.
  The inflow profile takes ``wind.upstream_roughness`` (default ``surface.roughness``), or with ``auto`` a value per
  direction from a regional map (`Upstream roughness of each direction`_).
- **Closure.** ``model.lee_setting = recommended`` writes the recommended lee setting of :doc:`howto`
  (``turb.realizability = 0.309``, ``turb.nut_relax = 0.5`` and the convection blend below
  ``model.blend_height``).
- **Nests.** Each ``nests`` entry becomes a folder ``nest_<name>`` with a patch deck at half the parent's cell size
  (refinement 2 2 1), its faces moved to the parent's nearest cell faces and at least two parent cells inside the
  parent's faces. ``run.sh`` runs the parent, then each patch with ``bc.nest_parent`` set to the parent's last
  checkpoint.
- **Sweeps.** ``wind.directions`` (a list, or ``start``, ``stop``, ``step``) writes one folder ``dir_<deg>`` per
  direction with its own deck; ``run.sh`` runs them in turn and skips those that already hold their answer
  (`Stopping and going on`_).
- **Anything else.** ``gust_inputs`` writes GUST inputs as they are, after everything else: a key the tool also
  writes is replaced, and ``null`` removes it. They go to the 3-D decks (the run folder's and each patch's), not to
  ``column/inputs``, the 1-D column run, which has its own mesh and keys.

Every key is checked before anything is written: an unknown key, a value out of range or a combination GUST would
refuse stops the tool with a message that names the key. Nothing is written then (exit status 2).

The case folder
---------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - File
     - Content
   * - ``inputs``
     - The deck (``dir_<deg>/inputs`` for a sweep). Change the YAML file and write the case again, rather than this
       file.
   * - ``column/inputs``
     - The column run of a geostrophic drive.
   * - ``nest_<name>/inputs``
     - The deck of a nested patch.
   * - ``terrain.xyz``, ``roughness.xyz``
     - The ground and the roughness map, in the case frame.
   * - ``run.sh``
     - Runs everything in order. ``GUST=/path/to/gust`` and ``MPI="srun -n 64"`` change the executable and the
       launcher (``MPI=`` runs without one); ``./run.sh key=value ...`` passes more inputs to the 3-D runs (for
       example ``steady.max_iterations=1``). A key the deck already holds stops GUST (a key given twice); change it
       in the YAML file, with ``--set``, or with ``gust_inputs``. A direction that already holds its answer is
       skipped (below); ``RESUME=0 ./run.sh`` runs every direction.
   * - ``case_summary.txt``
     - Domain, mesh, wind, notes and the deck checker's findings. With a geostrophic drive, ``run.sh`` adds a
       section "Library" after the runs (below).
   * - ``library.json``
     - Written by ``run.sh`` after the runs of a geostrophic drive (``make_gust_case.py --collect``): for each
       direction the geostrophic wind :math:`G` the column run found, the turning angle, the column's wind at
       ``wind.height``, its :math:`u_*`, and, for a neutral column, the neutral drag law's :math:`G` and angle for
       that :math:`u_*` (below).
   * - ``wind_library/``
     - Written by ``tools/make_gust_library.py build`` after the runs of a sweep: each direction's wind at heights
       above the ground, for the blend between directions (below).
   * - ``case_resolved.yaml``
     - The case description with every default filled in.

Stopping and going on
---------------------

A sweep of many directions can be stopped (a laptop closed, a job out of time) and started again: ``run.sh`` and
the window's **Run** skip each direction that already holds its answer (``run.resume = true``, the default), and
print a line for each. A direction holds its answer when

- its record ``.gust_ran`` (in ``dir_<deg>/``, or in the case folder for one direction) equals its decks now (the
  run's, ``column/inputs`` and each ``nest_<name>/inputs``, byte for byte), the case's data list
  ``.gust_case_data`` (the SHA-256 of each data file the tool wrote: ``terrain.xyz``, ``roughness.xyz``, the patch
  terrains and the copied files) and the ``key=value`` arguments given to ``run.sh`` now (the window gives none); and
- the run's ``log.txt`` and each patch's end in GUST's steady verdict (``gust: steady after N iterations``) followed
  by a plotfile ``<io.plot_file><step>`` that is still there with its ``Header``.

Everything else runs again: a run that stopped at ``steady.max_iterations`` (``not steady``), diverged, was stopped
or crashed; a direction whose plotfile was removed; a case written again with other settings or other terrain,
roughness or canopy files; other arguments. A case written again with the same settings and files keeps its decks
and data list byte for byte, so its finished directions stay skipped.
The whole direction runs again, its column run and patches included. When a direction starts, its old logs are
removed and then its record is written, so a direction stopped part-way never counts as finished.

``run.resume = false`` keeps the earlier behaviour, every direction run each time; ``RESUME=0 ./run.sh`` does the
same for one call (and ``RESUME=1`` skips with ``run.resume = false``); another value of ``RESUME`` stops
``run.sh``. The skip line reads ``run.sh: wind from 260 degrees: finished earlier with these decks, data files and
arguments, skipped``.

The test does not see:

- another GUST executable (``GUST=``, ``run.executable`` or a rebuild) or a data file edited by hand after the case
  was written: use ``RESUME=0`` (or ``run.resume = false``) to run again;
- an ``io.plot_file`` given as a ``run.sh`` argument, a run that writes no plotfile (``io.write_plotfile = 0``)
  and an unsteady run (``gust: unsteady run ended``): these always run again.

A joint run's verdict (``gust: joint run: <name> steady after N joint iterations``) counts like the run's own. The
check is per direction: a patch that stopped part-way runs the direction's column and terrain runs again too.

A steady verdict is what GUST itself accepts: every residual below ``steady.residual_tol`` and four orders below
its reference, and the monitors settled, or the confidence stop when ``steady.confidence`` is set.

Starting a direction from the converged neighbouring direction was measured and is not offered. On the hill of
``Tests/Infrastructure/Library_Hill`` (2 ranks, no convection blend), 260 degrees started with ``io.restart`` from
the 265 degree checkpoint (the same inflow and outflow faces) took 128 level-0 iterations; from the usual first
guess it took 45 iterations on the coarse level and 148 on level 0, about 160 level-0 iterations of work. The
restart saved about 20 % of the work and reached the same answer: the two differ by at most
:math:`4\times10^{-5}` m/s, which is the iteration error at :math:`r = 10^{-4}` (both are within
:math:`4\times10^{-5}` m/s of the same run taken to :math:`r < 10^{-6}`). GUST refuses ``io.restart`` with the
convection blend of the recommended lee setting (``convection.hybrid_sigma_max``), and the blend moves the wind by
about 0.8 % of the reference speed at 10 m (1 % in the first cell), so the restart cannot give the recommended
setting's answer. Where the faces differ (for example 265 to 270 degrees, where the y faces change from inflow and
outflow to slip walls), the flow would first have to turn to the new faces; that case was not measured.

The geostrophic wind of each direction
--------------------------------------

A geostrophic drive runs a 1-D column for every direction before its 3-D run; each column finds the geostrophic
wind :math:`G` that gives ``wind.speed`` at ``wind.height`` from that direction. The last step of ``run.sh`` reads
what the columns left and writes ``library.json`` and a table at the end of ``case_summary.txt``
(``tools/make_gust_case.py --collect <case folder>`` does the same by hand; a second collect replaces the table).
The window (``gust_gui.py``) does the same after its last run.

- ``run.sh`` removes ``library.json`` before its first run, and ``--overwrite`` removes it too, so a file left
  from earlier runs never stays beside new ones.
- A failed collect does not fail ``run.sh``: the runs are done, and ``run.sh`` prints why ``library.json`` was not
  written and how to collect again.
- The collect step needs the case tool's Python packages (PyYAML and numpy). ``GUST_CASE_TOOL`` names another
  ``make_gust_case.py`` (the path written into ``run.sh`` is that of the tool that wrote it) and ``PYTHON``
  another interpreter.
- ``--collect`` reads a written case folder only; it is refused with a case file or any other option but
  ``--quiet``. Exit status: 0 written, 2 refused (nothing written).

``library.json`` holds ``format_version`` (1), ``case``, ``drive``, ``reference_wind`` (``speed``, ``height``:
the wind asked for, null when ``atmosphere.geostrophic_wind`` is given), ``conventions`` and one entry per
direction:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Field
     - Meaning
   * - ``direction``, ``folder``
     - The direction asked for (degrees, from; with ``atmosphere.geostrophic_wind``, ``wind.direction``) and its
       run folder (``dir_<deg>``, or ``.`` without a sweep).
   * - ``geostrophic_wind``
     - :math:`G`: ``ug`` and ``vg`` (east and north, m/s), ``speed``, ``direction`` (degrees, from) and
       ``adjusted_by_column``. With a wind asked for, the column log's last lines ``forcing.ug = ...`` and
       ``forcing.vg = ...`` (9 digits, the values ``run.sh`` passes to the 3-D run), read only after the line
       "target wind reached" (a log without it is refused). With ``atmosphere.geostrophic_wind``, that wind, from
       the column deck.
   * - ``wind_at_height``
     - ``height`` (m), ``u``, ``v``, ``speed`` (m/s) and ``direction`` (degrees, from): the column table's wind at
       ``wind.height``, computed as GUST's column computes it to meet its target (``wind_at`` in
       ``Source/ColumnGuess.cpp``): linear in :math:`\ln(z + z_0)` between the two cell centres around the height;
       below the first centre the log law along the first centre's wind (a stratified column is refused there,
       as its profile needs :math:`\psi_m`). It must match the wind asked for to :math:`10^{-3}` and 0.06 degrees
       (GUST stops at :math:`|R - 1| < 10^{-4}`), or the table is refused as one from another run. It is the
       reference wind of the upstream flat column, not a point on the terrain. Null without a wind asked for.
   * - ``turning_deg``
     - The angle from ``wind_at_height`` to :math:`G` (degrees, counter-clockwise positive): negative in the
       northern hemisphere, where :math:`G` is turned clockwise from the wind near the ground. Null without a
       wind asked for.
   * - ``turning_from_stress_deg``
     - The angle from the wind at the first cell centre (the direction of the surface stress) to :math:`G`.
   * - ``wind_at_top``
     - The table's top level: ``height``, ``u``, ``v``. Not :math:`G`: the column's wind aloft can be
       supergeostrophic and turned from it.
   * - ``ustar``, ``inv_L``, ``surface_heat_flux``
     - From the table's header: :math:`u_*` (m/s), :math:`1/L` (1/m) and the kinematic heat flux (K m/s).
   * - ``roughness``, ``latitude``
     - The column deck's ``bc.surface_z0`` (m) and ``forcing.latitude`` (degrees): the values the column ran with,
       also when the latitude came from ``domain.center_lat_lon``.
   * - ``drag_law``
     - For a neutral column (:math:`1/L = 0` and no heat flux), the neutral drag law's ``speed`` and
       ``turning_deg`` for the column's :math:`u_*`, with ``A`` and ``B``:

       .. math::

          G = \frac{u_*}{\kappa}\sqrt{\left(\ln\frac{u_*}{|f| z_0} - A\right)^2 + B^2}, \qquad
          \alpha = \mathrm{atan2}\left(B,\ \ln\frac{u_*}{|f| z_0} - A\right),

       :math:`\kappa = 0.41`, :math:`f = 2\Omega \sin\varphi` (:math:`\Omega = 2\pi/86400` s\ :sup:`-1`,
       :math:`\varphi` the latitude), ``turning_deg`` :math:`= -\alpha` in the north and :math:`+\alpha` in the
       south. :math:`\alpha` is the angle from the surface stress, so compare it with
       ``turning_from_stress_deg``. Null for a stratified column.

A column and the drag law need not agree: the drag law's constants are empirical, and the column's turning depends
on its turbulence model, its length-scale cap (``turb.max_geom_lscale``) and its top. For 8 m/s at 10 m over
:math:`z_0 = 0.03` m at 45 degrees north (``rt_case_setup_column_axis``) the column gives :math:`|G| = 15.909` m/s
and turns 18.36 degrees from the stress; the drag law gives 15.944 m/s and 23.51 degrees.

Directions follow the case description: the direction the wind comes from, clockwise from north. The neutral
column is rotation-invariant (it is solved for :math:`u + iv`), so :math:`|G|` and the turning angle are the same
in every direction to round-off; ``rt_case_setup_column_axis`` checks this from 0, 270 and 240 degrees.

Upstream roughness of each direction
------------------------------------

The inflow profile of a direction (the log law's :math:`u_*`, ``bc.inflow_z0``, and the column's
``bc.surface_z0`` with a geostrophic drive) needs the roughness of the ground the wind crossed before it reached
the domain. ``wind.upstream_roughness = auto`` takes it from a regional roughness map,
``surface.regional_roughness_map`` (the same formats as ``surface.roughness_map``: ``x y z0`` points or a land-cover
raster), which reaches beyond the domain on the sides the wind comes from. ``surface.roughness`` or
``surface.roughness_map`` still sets the ground inside the domain.

The method is empirical: a weighted mean of the upwind roughness, not a model of the flow over the regional map.
For a wind from a direction, transects start on the faces it enters through and run upwind (towards the side the
wind comes from). Along each, :math:`\ln z_0` is averaged with the weight :math:`e^{-x/L}`:

.. math::

   \ln z_{0,\mathrm{up}} = \frac{\sum_k w_k \ln z_0(x_k)}{\sum_k w_k}, \qquad
   w_k = e^{-x_k/L} - e^{-x_{k+1}/L},

with :math:`x` the distance upwind of the face, :math:`L` ``wind.fetch_length``, segments
:math:`[x_k, x_{k+1}]` out to :math:`5L` (the weight beyond is :math:`e^{-5}`, 0.7 %) and :math:`\ln z_0` taken at
each segment's middle. The transects are spread along the inflow faces and weighted by the flow through them (face
length times the wind's component normal to the face), so the direction's value is the flux-weighted mean over its
inflow faces. The faces are those the deck uses: with a geostrophic drive the deck turns the wind and may add
side faces, so the transects then start on those faces, run along the wind near the ground, and are weighted by the
turned wind. The turning depends on :math:`z_0`, so the faces and :math:`z_{0,\mathrm{up}}` are found together
(at most five passes; the summary gives the turning).

- **Default** ``wind.fetch_length`` = 100 ``wind.height`` (1 km for a 10 m wind), an empirical choice. :math:`L`
  is an e-folding length: a wind at the reference height feels a change of surface from about 10-20 heights
  upwind and is fully adjusted only after a few hundred, so 100 heights lies in that range. Give
  ``wind.fetch_length`` to change it; it is required when ``wind.height`` is not given (for example with
  ``atmosphere.geostrophic_wind``).
- **Geometric mean.** The average is of :math:`\ln z_0`. An average matched to the surface drag would lean towards
  the rough value: half 0.3 m and half 0.03 m gives about 0.13 m against the geometric 0.095 m.
- **Coverage.** A segment outside the map is left out of both sums. Below 95 % of the weight covered, the
  summary says so; below :math:`1 - e^{-1}` (63 %, the map reaches less than one :math:`L` on average) the case is
  refused. The map is binned on a square grid of its point spacing (exact for a regular raster), and two passes
  fill an empty bin from its neighbours.
- **Mismatch.** ``case_summary.txt`` lists each direction's :math:`z_0`, the covered share and :math:`u_*`, and for
  each inflow face the upwind :math:`z_0`, the ground's :math:`z_0` along the face and
  :math:`\ln(z_{0,\mathrm{up}}/z_{0,\mathrm{ground}})`. A large mismatch starts an internal boundary layer at the
  face; put the face where the two agree, or far enough upwind of the area of interest.
- **Refused:** ``auto`` without ``surface.regional_roughness_map.file``, with ``wind.drive = couette`` (one lid
  speed for all directions), with ``wind.profile = table`` or ``power_law`` (the profile's wind is the same for
  every direction), and the regional map or ``wind.fetch_length`` without ``auto``.
- A geostrophic drive runs each direction's column with that direction's :math:`z_0`. With
  ``surface.roughness_map`` the column first guess is refused (it takes one roughness), so such a sweep uses
  ``geostrophic_inflow``.

``rt_case_setup_upstream`` checks the set-up on 2 ranks: a half-plane map (0.03 m west of the domain, 1 m east)
also sets the ground; from 270 degrees both decks take 0.03 m and the 10 m wind more than 4 cells from the x faces
varies by 0.62 % along x. Reading the direction as the one the wind blows towards gives 0.715 m and 11.8 %.

**Sensitivity to** :math:`L`. The same domain (2 km wide, ``geostrophic_inflow`` at 45 degrees north, 8 m/s at
10 m, ground 0.1 m) with 0.3 m for 450 m upwind of the inflow face and 0.03 m beyond (a map bin edge at the step),
on 2 ranks. These runs are not in CTest:

.. list-table::
   :header-rows: 1
   :widths: 25 25 25 25

   * - ``wind.fetch_length`` (m)
     - :math:`z_{0,\mathrm{up}}` (m)
     - 10 m wind at the centre (m/s)
     - 10 m wind, max - min along x
   * - 500
     - 0.1187
     - 8.106
     - 1.2 %
   * - **1000 (default)**
     - **0.0695**
     - **7.834**
     - **3.7 %**
   * - 2000
     - 0.0479
     - 7.665
     - 6.2 %

:math:`z_{0,\mathrm{up}}` is the closed form of the step to the four digits shown. Halving or doubling
:math:`L` moves the 10 m wind at the centre by +3.5 % and -2.2 %; the wind then adjusts to the 0.1 m ground
across the domain.

A wind library from a sweep
---------------------------

A sweep (``wind.directions``) runs one direction per folder. ``tools/make_gust_library.py`` turns the finished
runs into a wind library: each direction's wind at fixed heights above the ground, and the wind for any direction
in between, blended from the two directions either side. It needs PyYAML and numpy.

.. code-block:: bash

   tools/make_gust_library.py build CASE_DIR --heights 10 50 [--reference-point X Y] [--directions D ...]
   tools/make_gust_library.py query CASE_DIR --direction 272 --speed 6 --height 10 [--point X Y] [--out wind.npz]

``build`` reads each ``dir_<deg>/`` run of the case folder and writes ``CASE_DIR/wind_library/``. It is not the
``library.json`` of a geostrophic drive (above), which holds the columns' geostrophic winds. Every run is read
before anything is written; the new library is written beside the old one and replaces it whole, so a refused
build leaves the previous library as it was and a build of fewer directions leaves no stale sector files.

Positions (``--reference-point``, ``--point``, and ``x``, ``y`` in the files) are in the case frame, GUST's mesh
coordinates: :math:`x` and :math:`y` measured from ``domain.center``, heights shifted by ``terrain.zero_level``.

.. list-table::
   :header-rows: 1
   :widths: 28 14 58

   * - Option
     - Default
     - Effect
   * - ``--heights H ...``
     - 10
     - Heights above the ground (m). Each must be at or above the first cell centre and at or below the last
       in every column, and given once; otherwise refused.
   * - ``--reference-point X Y``
     - domain centre
     - Where the reference speed is taken, at the first height (m, case frame). Refused outside the outermost
       column centres.
   * - ``--directions D ...``
     - all
     - A subset of the sweep's directions; a direction not in the sweep, or one given twice, is refused.
   * - ``--quiet``
     - off
     - Print nothing when the build succeeds.

For each run the tool reads ``log.txt`` (GUST's verdict, the iterations, the last relative residual :math:`r`)
and the last plotfile ``<io.plot_file><step>`` written after the verdict (cell fields ``u``, ``v``, ``w``,
``tke``, ``z`` and the node heights ``Level_0/z_phys_nd``; the prefix is the deck's, so an ``io.plot_file`` given on
the command line is not seen). Refused, with the run folder named: a run with no ``inputs``, no verdict, or no
plotfile after the verdict (``io.write_plotfile = 1`` writes one); an unsteady run; a two-way joint run
(``nest.names``); a mesh with fewer than 2 cells in a direction; a sampled field that is not finite; no wind at
the reference point. A run that ended "not steady" or "diverged" is kept, marked ``"converged": false``, and the
query refuses it. All runs must share one mesh.

**Heights above the ground.** The ground under a column is the mean of its four bottom node heights, and a cell
centre's height above the ground is :math:`h = z_c - z_g`. Between two cell centres a field is linear in
:math:`h`. Below the first centre the profile depends on the wall law, so the library does not extrapolate there. Over the
log layer, linear sampling between centres reads low: 0.5 to 0.7 % at 10 m on the test hill below.

**What is stored.** Let :math:`S_i` be the horizontal wind speed of run :math:`i` at the reference point and the
first height (bilinear between column centres). The library stores

.. math::

   \hat u_i = u_i / S_i, \quad \hat v_i = v_i / S_i, \quad \hat w_i = w_i / S_i, \quad \hat k_i = k_i / S_i^2

at each height, in the case frame (:math:`u` towards :math:`+x`, east; :math:`v` towards :math:`+y`, north).
``wind_library/`` holds ``index.json``, ``grid.npz`` (the column centres ``x``, ``y`` and the ``ground``) and one
``sector_<deg>.npz`` per direction (``u``, ``v``, ``w``, ``k``, each of shape (heights, nx, ny)).

``index.json`` holds ``format_version`` (1), ``case``, ``heights``, ``reference_point``, ``reference_height`` (the
first height), ``grid``, ``frame``, ``stored`` and one entry per direction, sorted by direction:

.. list-table::
   :header-rows: 1
   :widths: 28 72

   * - Field
     - Meaning
   * - ``direction``, ``folder``, ``file``
     - The sweep's direction (degrees the wind comes from, clockwise from north), its run folder and its ``.npz``.
   * - ``status``, ``converged``, ``iterations``, ``last_r``
     - GUST's verdict (``steady``, ``not steady``, ``diverged``), whether it is ``steady``, the iterations and
       the last :math:`r` of the iteration table before the verdict (``null`` when the log has none).
   * - ``plotfile``
     - The plotfile read, relative to the case folder.
   * - ``reference_speed``, ``reference_direction``
     - :math:`S_i` (m/s) and the direction of the wind at the reference point (degrees, from). Over terrain this
       differs from ``direction`` where the ground turns the flow.
   * - ``surface_z0``, ``inflow_z0``
     - The deck's ``bc.surface_z0`` and ``bc.inflow_z0`` (m); ``null`` when the deck does not set them.

**The blend.** For a direction :math:`\theta` between the stored directions :math:`\theta_a` and
:math:`\theta_b` (clockwise, across north when needed), with :math:`\Delta = (\theta_b - \theta_a) \bmod 360`
and :math:`t = ((\theta - \theta_a) \bmod 360) / \Delta`, and a speed :math:`S` asked for. :math:`\theta` is
the sweep's inflow direction (the weights use the sweep's directions, not ``reference_direction``), and
:math:`S` is the local horizontal wind speed at the reference point and first height, not the inflow speed (over
a crest the two differ):

.. math::

   \mathbf{u}(\theta, S) = \frac{S}{L}\left[(1 - t)\,\hat{\mathbf{u}}_a + t\,\hat{\mathbf{u}}_b\right], \qquad
   k(\theta, S) = S^2\left[(1 - t)\,\hat k_a + t\,\hat k_b\right],

   L = \left|(1 - t)\,\mathbf{e}_a + t\,\mathbf{e}_b\right|,

where :math:`\mathbf{e}_a` and :math:`\mathbf{e}_b` are the unit vectors of the two reference winds
(``reference_direction``). Each sector's reference wind has length 1, so the blend at the reference point has
length :math:`L` before the division and :math:`S` after it. Without the division the blend midway would be
:math:`\cos(\Delta/2)` short: 0.4 % for 10 degree sectors and 1.5 % for 20 degree ones. That shortfall comes
from the two directions partly cancelling, not from a weaker flow, so :math:`k`, a scalar, is scaled by
:math:`S^2` only (dividing it by :math:`L^2` as well puts it 3 % high at 270 degrees on the test hill). The
blended wind at the reference point points within 0.04 degrees of :math:`\theta` for 20 degree sectors (the
vector blend is on the bisector only midway). A direction within :math:`10^{-6}` degrees of a stored one takes
that sector alone.

Refused by the query (the message names ``--height``, ``--speed`` or ``--direction``): a height the library did not
store; a speed, direction or height that is not a finite number (``bool`` included); a speed that is not above 0; a
library of another ``format_version`` or a sector file of the wrong shape; a direction between two
stored directions 180 degrees or more apart; any direction but its own for a library of one direction; a sector
marked ``"converged": false`` when its weight is not zero (the message names its folder).

``query`` prints the weights and ``u``, ``v``, ``w``, the speed, the direction and ``k`` at ``--point`` (default
the reference point); ``--out`` writes the fields at the height (``u``, ``v``, ``w``, ``k``, ``x``, ``y``,
``ground``) to a ``.npz`` file. The same is available from Python:

.. code-block:: python

   from gust_library import Library        # tools/ on sys.path
   w = Library("CASE_DIR").wind(272.0, 6.0, 10.0)   # direction from, speed at the reference point, height
   w["u"], w["v"], w["k"]                            # arrays (nx, ny)

Exit status of ``make_gust_library.py``: 0 done; 2 refused, with the option or run folder named.

**Ensemble statistics.** ``ensemble`` gives statistics of the wind over a distribution of directions and speeds,
from the library alone (no new runs):

.. code-block:: bash

   tools/make_gust_library.py ensemble CASE_DIR --height 10 --members members.txt [--percentiles 10 50 90]
       [--thresholds 10 15] [--grids DIR] [--points points.txt] [--out stats.npz]
   tools/make_gust_library.py ensemble CASE_DIR --height 10 --directions 240 270 300 --direction-weights 1 2 1
       --speeds 5 10 --speed-weights 3 1

A member is a direction (degrees the wind comes from, the sweep's inflow direction), a speed :math:`S` at the
reference point and first height (as in ``query``) and a weight. ``--members`` reads one member per line,
``direction speed [weight]`` (weight 1 when left out; spaces or commas; ``#`` starts a comment). Otherwise every
``--directions`` is taken with every ``--speeds``, weighted by the product of ``--direction-weights`` and
``--speed-weights`` (independent distributions; equal weights when left out). The weights :math:`w_i` are divided
by their sum. Each member's wind is the query's, :math:`\mathbf{u}_i = \mathbf{u}(\theta_i, S_i)`, with local speed
:math:`s_i = |\mathbf{u}_i|` (horizontal). Per column, at the height:

.. math::

   \bar s = \sum_i w_i s_i, \qquad
   \sigma_s = \Big(\sum_i w_i (s_i - \bar s)^2\Big)^{1/2}, \qquad
   \bar{\mathbf{u}} = \sum_i w_i \mathbf{u}_i,

   s_P = \min\Big\{ s_j : w_j > 0,\ \sum_{i:\, s_i \le s_j} w_i \ge P/100 \Big\}, \qquad
   p_T = \sum_{i:\, s_i > T} w_i, \qquad
   p_\mathrm{rev} = \sum_{i:\, \mathbf{u}_i \cdot \mathbf{e}_i < 0} w_i,

where :math:`\mathbf{e}_i` is the unit vector towards which a wind from :math:`\theta_i` blows. :math:`s_P` is a
weighted percentile without interpolation between members; a member of weight 0 is never a percentile, and the
cumulative weight counts as reaching :math:`P/100` within :math:`\max(10^{-12}, 4 n \epsilon)` for :math:`n`
members (the rounding of a sum of :math:`n` weights). :math:`p_T` is the probability of a speed above
:math:`T`, and :math:`p_\mathrm{rev}` the probability of a reversal: the local wind blows against the member's own
inflow direction, as in a recirculation. The direction of the vector mean :math:`\bar{\mathbf{u}}` is the mean
direction (no data where :math:`|\bar{\mathbf{u}}| \le 10^{-12}\,\bar s`, directions that cancel). The statistics hold only as far as the members' winds do:
each member is a blend and scaling of the library's runs, with the limits below.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Option
     - Default
     - Effect
   * - ``--members FILE``
     - none
     - The members, one per line. Refused together with ``--directions``/``--speeds``.
   * - ``--directions D ...``, ``--speeds S ...``
     - none
     - Every direction with every speed. Each speed must be above 0 m/s.
   * - ``--direction-weights W ...``, ``--speed-weights W ...``
     - equal
     - One weight per direction (per speed); each finite and 0 or more, with a sum above 0.
   * - ``--percentiles P ...``
     - 10 50 90
     - Weighted percentiles of the speed, each in (0, 100], each once (kept to 12 digits in the names).
       ``--percentiles`` with no values gives none.
   * - ``--thresholds T ...``
     - none
     - Speeds (m/s, 0 or more, each once) for the probability of a speed above each.
   * - ``--out FILE``
     - none
     - The statistics (``mean_speed``, ``std_speed``, ``mean_u``, ``mean_v``, ``mean_direction``, ``p<P>``,
       ``pexceed_<T>``, ``preversal``, ``x``, ``y``) as a ``.npz`` file; a ``.`` in :math:`P` or :math:`T` is
       written ``p`` (``pexceed_7p5``). ``mean_u`` and ``mean_v`` are in this file and the CSV, not written as
       grids.

A weight, speed or direction that is not a finite number, a negative weight, weights that sum to 0, and a
percentile or threshold out of range are refused, with the option named (exit status 2). A member whose
direction the library cannot give (a neighbouring run without GUST's steady verdict) is refused with the member's
number and direction. The weights of ``--members`` need not sum to 1. The columns are taken a few rows at a time,
so the memory does not grow with the number of members times the columns.

**WindNinja-format output.** ``query`` and ``ensemble`` write the wind at the height as WindNinja writes its
output grids, for fire models and GIS:

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Option
     - Default
     - Effect
   * - ``--grids DIR``
     - none
     - ESRI ASCII grids in ``DIR``: ``<prefix>_vel.asc`` (speed) and ``<prefix>_ang.asc`` (the direction the wind
       comes from, degrees clockwise from north), each with a ``.prj``. For ``ensemble``, ``_vel`` is the mean
       speed :math:`\bar s`, ``_ang`` the direction of the vector mean, and each other statistic gets
       ``<prefix>_<name>.asc`` (``std_speed``, ``p50``, ``pexceed_10``, ``preversal``, ...).
   * - ``--points FILE``
     - none
     - A CSV ``<prefix>_points.csv`` (in ``DIR``, or the current folder) of the values at the points of ``FILE``:
       one point per line, ``x y [name]`` (m, case frame). Columns: ``name``, ``x``, ``y``, ``easting``,
       ``northing`` (the case frame's origin added), ``height``, then ``u``, ``v``, ``w``, ``k``, ``speed``,
       ``direction`` (``query``) or the statistics (``ensemble``). Speeds in m/s. For ``ensemble`` the statistics
       at a point are those of the members' winds at the point, not interpolated from the column statistics.
   * - ``--prefix TEXT``
     - from the case
     - The start of the file names; by default ``<case>_<D>deg_<S>mps_<H>m_agl`` (``query``) or
       ``<case>_ensemble_<H>m_agl``.
   * - ``--speed-units``
     - ``mps``
     - The units of every speed grid: ``_vel``, and for ``ensemble`` also ``std_speed`` and ``p<P>``: ``mps``,
       ``mph`` or ``kph``. Probabilities and directions have no units; the CSV stays in m/s. WindNinja writes m/s
       by default; some fire models read mph or km/h.
   * - ``--prj FILE``
     - none
     - A projection file copied beside each grid, in place of the one written from ``domain.center_lat_lon``. It
       is read once before any grid is written, so it may be one of the outputs.

``--prj`` and ``--speed-units`` other than ``mps`` without ``--grids``, ``--prefix`` without ``--grids`` or
``--points``, a ``--prj`` file that does not exist and a ``--grids`` that is a file are refused before anything
is read. Everything else ``--grids`` and ``--points`` can refuse (the columns, the points file) is checked before
the first file is written.

The grids' cells are the library's columns: ``cellsize`` is the mesh cell, the corner (``xllcorner``,
``yllcorner``) is half a cell beyond the first column centre, the first row written is the north row, and the
no-data value is -9999. An ASCII grid has square cells: for a mesh with :math:`\Delta x \ne \Delta y` the fields
are resampled bilinearly onto square cells of :math:`h = \min(\Delta x, \Delta y)` from the domain's south-west
corner (the value at the outermost column centre beyond it), and the tool prints how. The direction is taken from
the resampled components. Columns of uneven spacing are refused. The corner is
placed on the map from the case folder's ``case_resolved.yaml``: the case frame's origin (the domain centre, or
its south-west corner with ``domain.origin``) in UTM for a centre given as ``domain.center_lat_lon``, and then the
``.prj`` is the ESRI text of WGS 84 / UTM in that zone; or in the projected frame of ``domain.center``, whose
projection the tool does not know, so no ``.prj`` is written unless ``--prj`` gives one. Without
``case_resolved.yaml`` (a ``wind_library/`` folder copied alone) the corner is in the case frame. Directions are
relative to the grid's north (the projection's northing); the meridian convergence of UTM is not removed. It is
about :math:`\Delta\lambda \sin\phi`: 0.75 degree at 100 km from the central meridian at 40 degrees latitude,
up to about 2 degrees at the edge of the zone, and more with a ``domain.utm_zone`` other than the centre's own.
The members' directions are in the same grid frame, while station winds are reported from true north; the
difference is the same convergence. The same is available from Python:

.. code-block:: python

   from gust_library import Library, ensemble, grids
   lib = Library("CASE_DIR")
   st = ensemble.ensemble(lib, ensemble.parse(["240 5 1", "270 10 2"]), 10.0, percentiles=(50,), thresholds=(8,))
   grids.export_grids(lib, "out", "mycase", st, ("mean_speed", "mean_u", "mean_v"), ["p50", "preversal"],
                      speed_names=ensemble.speed_names((50,)))

``Tests/Infrastructure/Library_Ensemble`` (CTest ``rt_library_ensemble``, 4 ranks) checks the statistics against
the runs and the grids' corner against the map.

**Accuracy and limits.** ``Tests/Infrastructure/Library_Hill`` (CTest ``rt_library_hill``, 2 ranks) blends a
library at 270 degrees and compares it with the 270 degree run, over a 60 m Gaussian hill with 50 m cells. The
horizontal-wind differences, as a percentage of :math:`S`, 4 cells or more from the lateral faces:

.. list-table::
   :header-rows: 1
   :widths: 30 14 28 28

   * - Library
     - Height
     - Largest difference, 4 cells or more from the faces
     - RMS, same cells
   * - 260/280 (20 degrees)
     - 10 m
     - 0.67 % of :math:`S`
     - 0.20 %
   * - 260/280 (20 degrees)
     - 50 m
     - 0.44 %
     - 0.11 %
   * - **265/275 (10 degrees)**
     - 10 m
     - **0.23 %**
     - 0.10 %
   * - 265/275 (10 degrees)
     - 50 m
     - 0.19 %
     - 0.09 %

The largest difference falls by about 2.9 from 20 to 10 degree sectors; the 10 degree values are close to the
runs' own convergence level (the 260 and 280 runs mirror each other to 0.3 % at :math:`r = 10^{-4}`), so this is not
a measured order of accuracy. The 260/280 library queried at 265 and 275 is within 0.51 % of those runs. The
interior mean of :math:`k` at 270 is 0.1 % (10 m) and 0.3 % (50 m) above the run's for the 260/280 library. Next to
the lateral faces the differences are larger (1.2 to 1.6 %), where a face is an inflow face for one direction and
a slip wall for the other. These numbers are for double precision; single precision is not measured.

- The blend is linear: a feature present in one direction and not the next, such as a lee recirculation, is
  smeared between them. Narrow sectors over steep terrain reduce this.
- Scaling by the speed holds for a neutral flow over rough ground driven by an inflow profile: the solution then
  does not depend on the speed, apart from molecular viscosity, the damping layers' rates (fixed in 1/s) and the
  wall law's 0.1 m/s minimum wind, together 0.02 % at half the run's speed. With a geostrophic drive the turning with height
  depends on the speed (Rossby number), so such a library holds at the speed its runs used.
- The fields are those of the final plotfile; a library does not read checkpoints or nested patches.

Command-line options
--------------------

A few keys can be changed on the command line without editing the YAML file; each replaces the file's value.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Key
   * - ``--name``, ``--output-dir``
     - ``case.name``, ``case.output_dir``
   * - ``--speed``, ``--height``, ``--direction``, ``--drive``
     - ``wind.speed``, ``wind.height``, ``wind.direction``, ``wind.drive``
   * - ``--directions DEG [DEG ...]``
     - ``wind.directions`` (replaces ``wind.direction``)
   * - ``--cell-size``, ``--layers``, ``--top``
     - ``mesh.cell_size``, ``mesh.layers``, ``domain.top``
   * - ``--terrain FILE``, ``--roughness``
     - ``terrain.file``, ``surface.roughness``
   * - ``--max-iterations``, ``--ranks``, ``--executable``
     - ``solver.max_iterations``, ``run.ranks``, ``run.executable``
   * - ``--set section.key=value``
     - Any key (the value read as YAML); ``--set gust_inputs.steady.cfl_max=8`` sets one GUST input.

``--dry-run`` checks the case and prints the summary and the decks without writing; ``--overwrite`` writes into a
folder that exists; ``--no-check`` skips the deck checker; ``--strict`` makes a checker warning count as a failure;
``--list-keys`` prints the key table below. Exit status: 0 written, 1 written but the deck checker found an error (or
a warning with ``--strict``), 2 refused.

.. code-block:: bash

   python3 tools/make_gust_case.py case.yaml --direction 240 --name from_240
   python3 tools/make_gust_case.py case.yaml --directions 0 90 180 270 --set mesh.layers=48

Checking a deck before a run
----------------------------

``tools/check_gust_inputs.py`` reads a deck the way GUST does and reports, without running GUST, what would stop the
run and what is likely not what was meant. It follows the design of ERF's ``check_erf_inputs.py``: a deck parser,
derived facts of the run, and independent checkers, each reporting findings of four levels.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Level
     - Meaning
   * - ERROR
     - GUST would stop at start-up (unknown or repeated key, value out of range, a key with no effect, a missing
       file, a combination the named set-ups refuse).
   * - WARNING
     - GUST starts, but the setting is very likely not what was meant (for example an inflow face the wind leaves
       through).
   * - ADVICE
     - A setting that could be chosen better (boxes per rank, the sweep for the build).
   * - INFO
     - Derived numbers: cell sizes, layers, the stretch ratio, boxes, the wind direction.

.. code-block:: bash

   python3 tools/check_gust_inputs.py inputs
   python3 tools/check_gust_inputs.py inputs io.plot_file=test_plt       # command-line keys, as for gust
   python3 tools/check_gust_inputs.py --ranks 16 --min-severity WARNING --strict Tests/*/*/inputs*

A command-line key the deck already holds is an ERROR, as in GUST: the decks this tool writes hold
``steady.max_iterations``, so change it in the YAML file instead. Files named in a deck are looked up from the current
directory, as GUST does: run the checker from the folder GUST will run in (``run.sh`` runs each deck from its own
folder). ``FILE = name`` lines are read in place from that folder, as AMReX does; a file that cannot be
opened, or that includes itself, is an ERROR.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Option
     - Effect
   * - ``--ranks N``
     - The MPI ranks of the run; adds the boxes-per-rank advice.
   * - ``--gpu``
     - The run uses a GPU build: the default sweep and the box-size advice for a GPU.
   * - ``--min-severity LEVEL``
     - Hide findings below ``INFO``, ``ADVICE``, ``WARNING`` or ``ERROR``.
   * - ``--strict``
     - Exit 1 on a WARNING too.
   * - ``--json``
     - Print the findings as JSON.
   * - ``--only NAME``
     - Run only the named checker (repeatable); ``--list-checkers`` lists them.
   * - ``--reference PATH``
     - Another ``inputs_reference`` (default ``Tests/inputs_reference``).
   * - ``--made-by-run PATH``
     - A file or checkpoint an earlier step of the run writes (repeatable).
   * - ``--self-test``
     - Run the checker's own table of cases and exit.

The exit status is 0 when no deck has an ERROR (and, with ``--strict``, no WARNING), 1 otherwise, and 2 for a
command-line error.

The ranges and defaults come from ``Tests/inputs_reference``, which is generated from :doc:`inputs`, so the checker
follows the readers. ``--made-by-run PATH`` names a file or checkpoint an earlier step of the run writes (the case
tool passes the column table and the parent checkpoint this way). The checker reads the deck only; a deck it passes
can still stop in GUST on something only the run knows, such as the terrain file's coverage at every node.

Downloading terrain
-------------------

``tools/fetch_terrain.py`` downloads the elevation for a case's domain and its margin and writes it as an ESRI ASCII
grid in longitude and latitude, which the builder reads without rasterio:

.. code-block:: bash

   python3 tools/fetch_terrain.py case.yaml                        # dem.asc beside case.yaml
   python3 tools/fetch_terrain.py case.yaml --out terrain/dem.asc --sea-floor keep
   python3 tools/fetch_terrain.py --box 57.15 57.23 -7.45 -7.30 --out askervein.asc

The case file needs ``domain.center_lat_lon`` (the download box comes from it); then give the grid as
``terrain.file``. The tiles are the SRTM-format ("skadi") set of the AWS Terrain Tiles: one 1° × 1° tile per file,
1 arc-second (about 30 m) where SRTM covers the ground, heights in metres above the EGM96 geoid. No key or account
is needed.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Option
     - Effect
   * - ``--out PATH``
     - The grid to write (``.asc``); default ``dem.asc`` beside the case file, or in the current folder with
       ``--box``.
   * - ``--sea-floor zero|keep``
     - Over the sea the tiles hold the sea floor (below 0 m). ``zero`` (default) sets points below 0 m to 0 m;
       ``keep`` leaves them, for land below sea level.
   * - ``--cache DIR``
     - Where the raw tiles are kept (default ``GUST_DEM_CACHE`` or ``~/.cache/gust/terrain_tiles``); a cached tile
       is not fetched again.
   * - ``--box S N W E``
     - A box in degrees instead of a case file.

A download of more than 16 tiles (4° × 4°) or one that crosses the antimeridian is refused, and so is a tile the
server does not have (its address is named). Voids are written as the grid's no-data value, not filled; the
download reports their number and the builder refuses a domain that holds one. Beside the grid,
``<name>.sources.txt`` records the tiles, the date and the attribution below.

The case set-up window
----------------------

``tools/gust_gui.py`` is a window over the builder and the checker:

.. code-block:: bash

   python3 tools/gust_gui.py                          # the starter case
   python3 tools/gust_gui.py case.yaml                # a case file

- **Basic mode.** The window opens with only the settings a case needs: the case name and folder, the centre in
  latitude and longitude, the half-width and top, the elevation file, the cell size and layers, the wind speed,
  height and direction, the roughness and the MPI ranks. Every other key keeps its default. **Show all settings**
  (toolbar and View menu) shows every key. A key that a file gives, or that a check refuses, shows in the basic mode
  too, so no value is hidden.
- **Form.** One tab per section of the case description and one field per key, with plain labels in groups (the
  domain's position on the Earth, its position in the terrain file's own coordinates, its size; the horizontal and
  vertical mesh). The field's tooltip names its YAML key and gives its unit, range and help. An empty field is a key
  not given; it shows in grey what to type (``number > 0``, ``easting, northing``, ``one of: auto, xyz``) and the
  default, if any. Latitude and longitude have one box each. Lists are typed with commas (``0, 500, 1000``);
  ``[0, 500, 1000]`` is read too. File and folder keys have a **Browse** button; a choice inside the case file's
  folder is written relative to it. Nested patches are added and removed on the ``nests`` tab. Keys the form does
  not know are kept when the file is saved.
- **Open and Save** read and write the same YAML as ``make_gust_case.py``; a file read and saved again holds the same
  values.
- **Check** resolves the case without writing it. A refused key marks its field in red (the refusal is in its
  tooltip) and the refusals fill the Findings table. The Map tab draws the terrain, the domain, the flat border and
  the end of the border blend, the nests with their blend, and the wind direction (the first one, when there are
  several), with the vertical layers beside it.
- **Write case** saves the file, writes the case folder and fills the Findings table with the deck checker's
  findings for every deck written.
- **Download terrain** runs the download above and sets ``terrain.file``.
- **Run** writes the case and runs gust step by step in ``run.sh``'s order (the column, the run, each patch from
  the run's newest checkpoint). The log streams into the Run tab and into each step's ``log.txt``, and the residual
  columns of the iteration table are plotted as they come (from the aligned table of ``io.iteration_log = table``,
  the default; with ``full`` nothing is plotted). gust is ``run.executable``, a relative path taken from the case
  folder; ``run.ranks`` above 1 uses ``run.mpi_launcher``. **Stop** asks the step to end and ends it after 5 s if
  it has not; the steps after it do not run. As ``run.sh``, it skips the directions that already hold their answer
  (`Stopping and going on`_) unless ``run.resume = false``.
- **Library** tab: the wind library of a case folder (filled in after **Write case**; **Browse** picks another).
  Choose the height and either **One wind** (a direction and a speed at the reference point, as ``query``) or
  **Distribution** (members typed as ``direction speed [weight]`` lines, percentiles and thresholds, as
  ``ensemble``), then **Show** draws the chosen field on the case frame. **Export grids** writes the WindNinja-format
  grids of `A wind library from a sweep`_ (in the speed units chosen) and **Export points** the CSV at the points of a
  file. Both export what the last **Show** drew. A refusal shows above the map. The tab has no ``--prj``: for a
  centre given as ``domain.center`` it writes no ``.prj`` (use ``make_gust_library.py --prj``). The tab reads a
  library; ``make_gust_library.py build`` writes it.

On Windows, the window can be frozen into a folder with PyInstaller (``tools/gust_gui.spec``: one folder, no console
window, rasterio left out):

.. code-block:: bash

   cd tools
   pyinstaller --noconfirm gust_gui.spec              # dist/gust_gui/gust_gui.exe

Put the CPU ``gust.exe`` in ``dist/gust_gui/`` beside ``gust_gui.exe``. When the case does not give
``run.executable``, the frozen window runs that gust; when it does, it runs the one named and falls back to the one
beside it. The frozen window reads GeoTIFF terrain with its built-in reader only.

Data and acknowledgements
-------------------------

The terrain download uses the `Terrain Tiles <https://registry.opendata.aws/terrain-tiles/>`_ on the Registry of
Open Data on AWS, made by Mapzen and maintained by the Tilezen project (`tilezen/joerd
<https://github.com/tilezen/joerd>`_). We thank Mapzen, Tilezen and Amazon Web Services for making them free to use
with the attribution below. The tiles merge several sources; their attribution, which a product made from them must
carry, is listed in `joerd's attribution file <https://github.com/tilezen/joerd/blob/master/docs/attribution.md>`_.
It includes SRTM (NASA and USGS), GMTED2010 (USGS), the 3D Elevation Program (USGS), ETOPO1 (NOAA), and national
models of several countries. No tile is redistributed in this repository; the download fetches them to the user's
cache.

The download follows FastWindTerrain's ``cases/casegen.py`` (the same tiles, a cache of raw tiles and no partial
file left behind) and the tile layout used by massconsistent_amr's ``terrain_reader_srtm.py`` and by amrTerrain.
No code is copied from them.

Case keys
---------

``Tests/inputs_yaml_reference`` holds the same keys as a commented YAML file to copy lines from: remove the leading
``# `` of the lines needed and what is left is a case file, with each key's meaning still a comment. It is written by
``tools/make_gust_case.py --write-yaml-reference Tests/inputs_yaml_reference``; the CTest ``gust_case_selftest``
fails when the file differs from what the tool writes, when the uncommented file is not YAML, or when a key or a
default in it differs from the table below. ``Tests/inputs_reference`` is its counterpart for GUST's own deck
inputs (:doc:`inputs`).

.. include:: case_setup_keys.inc

Limits
------

- Lateral faces are inflow, outflow, slip or periodic; one wind direction per run folder.
- The tool writes the column first guess for the geostrophic drives only; the other drives start from the profile.
- Nested patches are one level deep (a patch inside the parent, not inside another patch) and one-way: each
  patch is its own deck with ``nest`` faces. The tool does not write the nests of one run (``nest.names``) or a
  two-way option (``nest.feedback``).
- The tool writes steady decks (``steady.stage = simplec``). For the Newton stage or an unsteady run, change the
  stage through ``gust_inputs`` and set ``model.lee_setting: none``, since both stages refuse the upwind blend of
  the recommended setting; an unsteady deck also needs the ``unsteady.*`` keys and no steady-stop keys
  (:doc:`howto`). ``tools/check_gust_inputs.py`` checks decks of every stage and with nests.
- A terrain file must cover the domain and its margin; the tool refuses one with gaps rather than filling them.
