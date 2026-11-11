Buildings
=========

Buildings on the terrain shape the wind where it enters the wildland–urban interface. GUST represents them by
immersed forcing with the form 7 wall law of ERF (``erf.if_wall_form = fraction_stress``) and Kynema-SGF
(``wall_form 7``). The ground keeps its body-fitted wall law on the terrain-following mesh and is not blanked.

**Status.** The geometry acts on momentum with SIMPLEC and z-line sweeps on one level (`Momentum`_). Roofs,
undersides and sides have the wall law, and the closure sees the buildings (`Turbulence`_): :math:`k` has its wall
value beside them and the length scale their distance. :math:`\theta` has no building terms yet, and answers near
buildings are not validated. The run refuses
``steady.stage = simplec_newton``, ``steady.sweep = lusgs``, ``mesh.hierarchy_levels`` above 1 and
``solver.poisson_backend = mlmg``, and fluid sealed off by buildings.

Height map
----------

``buildings.file`` holds heights above the local ground on a rectangular grid of points, in one of two layouts
(``buildings.format``):

- ``erf``: ERF's ``erf.buildings_file_name``: ``nx``, ``ny``, the ``nx`` x values, the ``ny`` y values, then the
  ``nx ny`` heights, x the slow index.
- ``xyz``: ``x y height`` or ``x y height zstart`` lines holding every grid point once, in any order, comma or
  white-space separated, ``#`` comments allowed. Scattered points are refused, as interpolating them would smear the
  walls. The building is solid from ``zstart`` (default 0) up to ``height``, so ``zstart > 0`` leaves an open void deck
  under it. ``height`` is the top of the building, not its thickness.

Over terrain, ``height`` and ``zstart`` are measured from the ground of each column.

Each point owns a tile reaching halfway to its neighbours, the end tiles as far outward as inward. The height is
constant over a tile and 0 outside the tiles, so a wall lies on a tile edge and is vertical. In a periodic direction
the tiles wrap, and a file listing both ``x = 0`` and ``x = L`` is refused.

ERF interpolates the same points to the cell corners, which makes a ramp over a cell on the edge of a cube. With points
at every cell node, that column holds the same volume as here, in other layers; other point spacings can differ in
volume too.

Solid fractions
---------------

With :math:`h` the tile's height, :math:`z_s` its ``zstart`` and :math:`z_b, z_t` a cell's lower and upper heights above the ground of its column
(four-corner means, as in ``detJ``), the solid fraction :math:`\beta` is the exact mean over the cell's footprint
:math:`A`,

.. math::

   \beta = \frac{1}{|A|} \int_A \frac{\max\bigl(0, \min(h, z_t) - \max(z_s, z_b)\bigr)}{z_t - z_b}\, dA.

So a building follows the ground: on a slope its roof is parallel to the terrain. ``buildings.small_volfrac`` then
sets a fraction below it to 0 and one within it of 1 to 1. A cell is fluid at :math:`\beta = 0`, solid at
:math:`\beta = 1` and partial otherwise.

A vertical face takes the same mean along it, with its two-node mean heights; on a tile edge the taller tile counts. A
horizontal face takes the part of its area inside a building or at its roof or underside. A face touching a solid cell then has
fraction 1 and one between two fluid cells 0.

Wall cells
----------

The difference between the blocked fractions of a cell's opposite faces is the wall area in the cell facing that way,
and in a fluid cell any blocked face is a wall. Along z this gives ERF's ``ib_wall::is_wall_cell``: a partial cell with
fluid above it, or a fluid cell on a solid one. The plotfile's ``building_wall`` holds the sum of the bits of the walls of a cell:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Bit
     - Wall
   * - 1
     - Roof below the cell
   * - 2, 4
     - Side at the cell's low, high x
   * - 8, 16
     - Side at the cell's low, high y
   * - 32
     - Underside of a lifted building above the cell

A roof and a side make an edge cell, and two sides a corner. Along a normal with a wall, with :math:`\Delta_n` the
cell's size along it (:math:`\Delta x`, :math:`\Delta y`, or its height) and :math:`\Delta_r` that of the reference cell,
the next one out, the fluid centroid and the reference cell lie at

.. math::

   d_1 = \tfrac{1}{2}(1 - \beta)\,\Delta_n, \qquad d_2 = (1 - \beta)\,\Delta_n + \tfrac{1}{2}\Delta_r

from the wall; on equal layers :math:`d_2 = (\tfrac{3}{2} - \beta)\Delta_n`, as in ERF. :math:`d_2` is 0 when the reference
cell is solid or outside the domain, such as across a one-cell street or under a building in the lowest layer.

Kinds of building
-----------------

The figures show ``building_beta`` from GUST on a slice along x, with 5 m cells, wall cells outlined
(``Tests/Infrastructure/Building_Field/plot_building_kinds.py``). On the hill the same map follows the ground.

.. figure:: figures/buildings/kinds_flat.png
   :width: 100%

   Building kinds on flat ground.

.. figure:: figures/buildings/kinds_hill.png
   :width: 100%

   The same map on a 30 m hill.

Momentum
--------

**Fixed rows.** A face with fraction 1 has :math:`u = 0`. Its residual is the velocity, its implicit row the identity,
with no :math:`1/\Delta\tau`, and the SIMPLEC denominator is zero. This is the steady limit of ERF's
``if_implicit_projection``. No drag coefficient enters, so the answer does not depend on one.

The ground's ``w`` stays the boundary conditions' row under a building too, as they set it from the terrain's slope.

**Mass fluxes.** A partial face's flux is scaled by its open fraction :math:`1 - \beta_f`. The metric itself is unchanged,
as the pressure gradient reads it.

**Pressure correction.** Blocked faces have zero area, a solid cell's row is :math:`\phi = 0`, and the mean of the singular
problem is over the fluid cells. GMRES runs on ERF's operator. Its FFT or multigrid preconditioner keeps every face
and does not know the mask, so the number of iterations grows with buildings (8 for one, 11 for nine, 3 for none, on the
test's mesh with ``gmres_mlmg``). Under multigrid GMRES is flexible, as the cycle's bottom solver is not a fixed
linear map. Fluid with no open path to the rest has its own null vector, so the run refuses it.

**Wall law** (``buildings.wall_form = fraction_stress``). A wall cell has, neutrally, for each of its walls,

.. math::

   u_*^2 = (1 - \beta) u_{*,\mathrm{own}}^2 + \beta u_{*,\mathrm{ref}}^2, \qquad
   u_{*,x} = \frac{\kappa U_x}{\ln(\max(d_x, 10 z_0) / z_0)},

with :math:`\kappa = 0.41`, :math:`z_0` = ``buildings.z0``, :math:`d_x` the distances :math:`d_1` (own) and :math:`d_2`
(reference) along the wall's normal, :math:`U` the speed parallel to the wall in the cell (own) and in the next one out
(reference), and :math:`u_*^2 = u_{*,\mathrm{own}}^2` where :math:`d_2 = 0`. The speed is taken at the cell centre from the
two components parallel to the wall:

.. list-table::
   :header-rows: 1
   :widths: 34 22 22 22

   * - Wall
     - Normal
     - Parallel components
     - :math:`\Delta_n`
   * - Roof, underside
     - z
     - :math:`u, v`
     - the layer's height
   * - Side at low or high x
     - x
     - :math:`v, w`
     - :math:`\Delta x`
   * - Side at low or high y
     - y
     - :math:`u, w`
     - :math:`\Delta y`

The wall is a sink on each parallel component :math:`u_t`,

.. math::

   S = -r u_t, \qquad r = \frac{f\, u_*^2}{\Delta_n \max(|U_\mathrm{own}|, 10^{-3}\ \mathrm{m/s})},

with :math:`f` the blocked part of the wall's face and :math:`\Delta_n` the cell's size along the normal. The force on the
cell is then the stress :math:`u_*^2` on the wall's area, along the flow parallel to it.

**Edges and corners.** A cell with several walls takes the sum of their rates on each component: at an edge between a
roof and a side normal to x, :math:`u` has the roof's rate, :math:`w` the side's, and :math:`v`, which runs along the edge,
both. At a corner between two sides :math:`w` has both. Each wall uses its own parallel speed and its own distances.

A velocity row takes the mean :math:`r` of its two cells. The viscous flux of a parallel component through the blocked
part of a wall's face is closed, in the residual and in the implicit rows, and :math:`r` goes into :math:`a_P`, not into
the local pseudo-time step.

``buildings.wall_form = roof_stress`` applies the law to roofs and undersides alone: the sides keep their fixed rows,
with viscous fluxes to the zero velocity. This is what ``fraction_stress`` did before the sides had the law; on the two
test cases it gives that form's plotfiles bit for bit. A checkpoint written with that earlier ``fraction_stress``
restarts with ``roof_stress`` and is refused with today's ``fraction_stress``, whose answer differs. ``buildings.wall_form = drag`` keeps the fixed rows alone on every wall, as Kynema-SGF's laminar pathway.

.. figure:: figures/buildings/block_flow.png
   :width: 100%

   One block in a log-law wind (``Tests/Canonical_RANS/Building_Cube``), with the fixed rows, the wall law on the roof
   and the sides, and the closure's building terms. This shows the method working, not a validated flow.

.. figure:: figures/buildings/block_on_hill.png
   :width: 100%

   The same block on a ridge's windward flank (``Tests/Canonical_RANS/Building_Hill``), with its floor and roof following
   the ground. The same limits apply.

Turbulence
----------

With ``buildings.wall_tke = on``, the default with ``buildings.wall_form = fraction_stress``, the :math:`k`–:math:`L`
closure (:doc:`turbulence_model`) sees the buildings in four ways.

**Solid cells** hold :math:`k` at the closure's floor (``turb.tke_floor``): their row is :math:`R_k = k - k_\mathrm{floor}`.

**Wall cells** relax :math:`k` to the wall value of the law's friction velocity, as the ground's first cell does:

.. math::

   W = a\,(k - k_w), \qquad k_w = \max\!\left(\frac{u_*^2}{C_{\mu 0}^2},\, k_\mathrm{floor}\right), \qquad
   a = \frac{\max(u_*, 0.01\ \mathrm{m/s})}{d_w},

with :math:`u_*^2` the largest of the cell's walls, :math:`C_{\mu 0}` = ``turb.Cmu0`` and :math:`d_w` the wall cell's
:math:`d_b` below. :math:`k_w` is the value at which the closure's production equals its dissipation in a log layer
with that :math:`u_*`, the one the ground uses. How :math:`W` enters follows the ground's choice ``bc.surface_wall_k``,
with :math:`f` the largest blocked part of the faces the cell's walls lie on:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ``bc.surface_wall_k``
     - Wall cell's row
   * - ``dirichlet`` (default)
     - :math:`R_k \leftarrow (1 - f) R_k + f W`: a fully blocked face replaces the cell's :math:`k` equation, and a
       slightly blocked one leaves it nearly as it was. In the ground's first cell, which has the ground's own
       relaxation row, :math:`f W` is added to that row, so :math:`k` settles between the two wall values.
   * - ``source``
     - :math:`R_k \leftarrow R_k + f W` in every wall cell, as the ground adds its own.
   * - no surface layer
     - As ``dirichlet``, without the ground's row.

**Fluxes.** The diffusive flux of :math:`k` through the blocked part of a face is closed, in the residual and in the
implicit rows: :math:`k` does not diffuse into a building.

**Length scale.** The closure's wall distance becomes :math:`\min(d_\mathrm{ground}, d_b)`, so its geometric length
:math:`l_g = \kappa (d + z_0)` shortens near a building. :math:`z_0` there is the ground's (``bc.surface_z0`` or its map),
as before, and ``turb.max_geom_lscale`` still caps the length. The turbulence plotfile's ``walldist`` shows this
distance. The hybrid convection blend (``convection.hybrid_height``) keeps the ground's distance. :math:`d_b` is

- in a wall cell, the distance of its fluid centroid from its nearest wall, :math:`\tfrac{1}{2}(1 - \beta)\Delta_n`,
  not below ``buildings.z0``;
- elsewhere, the distance from the cell centre to the nearest wall of a building. Each wall is a rectangle on its
  plane: a roof :math:`\beta \Delta_z` above its cell's floor, a side :math:`\beta \Delta_x` in from its face, over the
  blocked part of the face's height.

The distance is found once, at start-up, by passing each cell's nearest wall to its 26 neighbours until no distance
falls (a vector distance transform). The log gives the number of passes and of cells whose distance the buildings set.
For one block it equals the distance to the block's box to round-off, and so it does for an L with a step, the
tests' second shape, against the smallest distance over its boxes.

``buildings.wall_tke = off`` leaves the closure as it was before these terms: :math:`k` is solved in every cell with
the ground's distance alone. A checkpoint written before them restarts with ``off`` and is refused with ``on``. ``on``
needs the wall law on every wall, so it is refused with ``roof_stress`` and ``drag``.

On the block of the first figure, :math:`k` in the roof's wall cell is 0.41 m²/s² with these terms and 8.6 m²/s²
without, against 0.35 m²/s² upstream at that height, and without them :math:`k` in the solid cells was 3 m²/s². On the
ridge of the second figure the air now reverses in a thin layer on the roof (u in the roof's wall cell -0.15 m/s; it was
3.90 m/s), as the eddy viscosity there is lower.

Limits
------

- A domain face is not a wall: outside a non-periodic side the fields take the nearest domain column. The log warns of a
  building on such a side.
- :math:`d_1` uses the cell's :math:`\beta`. In a cell whose :math:`\beta` comes from another direction, such as an edge
  cell filled only partly in plan and in height, it is an estimate.
- A building thinner than a cell has its opposite sides in one cell, so its faces carry no net wall, and whether it blocks
  a face depends on where it sits against the mesh. The log warns; refine the mesh around it.
- On a slope the operator's cross terms and the pressure gradient's read the solid cells' :math:`\phi = 0` and pressure, and
  a roof's speed is horizontal, not along the roof. On flat ground neither arises. A side is vertical on any ground, so
  its two components are parallel to it.
- A side's reference cell is the next one along x or y. Outside a non-periodic domain face it holds the boundary
  conditions' values.
- A wall's parallel speed is taken at the cell centre, each component the mean of its two faces. In a corner or edge
  cell one of those faces can lie on the other wall, a fixed row, so that component is half the open face's value. ERF
  evaluates the law at each velocity face from that face's own value; with the flow normal to the other wall this
  gives twice the cell's contribution here. The component along an edge and plain side cells are not affected.
- A velocity row takes the mean of its two cells, so where a wall ends partway along a row's direction the sink and
  the closed flux are spread over the two rows beside the end (three quarters and one quarter where the wall ends
  mid-cell). The total force is kept.
- The own distance is floored at :math:`10 z_0`. A wall cell less than :math:`20 z_0` across (2 m with the default
  ``buildings.z0``) is below it: the stress is then :math:`(\kappa / \ln 10)^2 U^2` whatever the cell size, lower than the
  log law at the true distance (about half for a 1 m cell with :math:`z_0` = 0.1 m). The log notes it for the sides.
- The distance to buildings takes one pass over the whole mesh per cell of distance it reaches. A wall as far as the
  mesh's tallest column is not looked for, since the ground is nearer, so the passes are at most that height in
  cells of the smallest size (23 on the block's mesh, 31 on the ridge's). It is done once per run, restarts included.
- The distance treats a roof as covering its cell's whole footprint and a side its face's whole width, so where a
  wall ends partway through a cell it can be short by up to half a cell. The transform passes nearest walls between
  neighbours and is not exact for every shape.
- Outside a non-periodic domain face the closure's ghost cells take the nearest domain column's distance.
- The wall value of :math:`k` is neutral and reads the speed parallel to the wall. On a wall that faces the flow that
  speed is small, so :math:`k` and the eddy viscosity in its wall cells are small too (0.008 m²/s² on the block's
  windward face, against 0.76 m²/s² one cell upstream). The ground's row does the same where the flow reattaches.
- The shear production in the second cell from a wall still reads the wall cell's velocity, and :math:`\theta` has no
  building terms yet.
- The length near a building is :math:`\kappa (d_b + z_0)` with the ground's :math:`z_0`, not ``buildings.z0``. The two
  differ where the ground's roughness is much larger than the walls'.
- A side's rectangle in a cell that is partly solid lies :math:`\beta \Delta_x` in from its face even where the cell's
  :math:`\beta` comes from a roof, so beside a step between two roofs the distance can be off by up to the layer's
  height.
- The thin reversed layer on the roof of the block on the ridge is one cell deep and lies in partly solid cells; the
  block on flat ground, whose roof cells are fluid, has none. Its presence follows the mesh as well as the flow.
- The residual norm of :math:`k` includes the solid cells' rows, :math:`k - k_\mathrm{floor}`, and the wall rows. From a
  first guess with :math:`k` well above the floor in the solid cells they set the first iterations' norm (4.0 against
  0.2 on the block); they fall within a few iterations.
- ERF and Kynema-SGF keep the :math:`k` equation in their wall cells and dissipate there with the flat first-cell
  length (``flat_wall_dissipation``). Here the wall cell's row is the relaxation, where the dissipation does not enter,
  so that length is not used.
- A partial cell's volume is not reduced: the sink and the fluxes are per full cell, which leaves the steady balance as it is.
- On strongly compressed terrain layers a face's fraction can differ a little from its cells', which can add a side bit
  at a roof.
