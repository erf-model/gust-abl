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
- **Closure.** ``model.lee_setting = recommended`` writes the recommended lee setting of :doc:`howto`
  (``turb.realizability = 0.309``, ``turb.nut_relax = 0.5`` and the convection blend below
  ``model.blend_height``).
- **Nests.** Each ``nests`` entry becomes a folder ``nest_<name>`` with a patch deck at half the parent's cell size
  (refinement 2 2 1), its faces moved to the parent's nearest cell faces and at least two parent cells inside the
  parent's faces. ``run.sh`` runs the parent, then each patch with ``bc.nest_parent`` set to the parent's last
  checkpoint.
- **Sweeps.** ``wind.directions`` (a list, or ``start``, ``stop``, ``step``) writes one folder ``dir_<deg>`` per
  direction with its own deck; ``run.sh`` runs them in turn.
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
       in the YAML file, with ``--set``, or with ``gust_inputs``.
   * - ``case_summary.txt``
     - Domain, mesh, wind, notes and the deck checker's findings. With a geostrophic drive, ``run.sh`` adds a
       section "Library" after the runs (below).
   * - ``library.json``
     - Written by ``run.sh`` after the runs of a geostrophic drive (``make_gust_case.py --collect``): for each
       direction the geostrophic wind :math:`G` the column run found, the turning angle, the column's wind at
       ``wind.height``, its :math:`u_*`, and, for a neutral column, the neutral drag law's :math:`G` and angle for
       that :math:`u_*` (below).
   * - ``case_resolved.yaml``
     - The case description with every default filled in.

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
  it has not; the steps after it do not run.

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
