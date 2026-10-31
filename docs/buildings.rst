Buildings
=========

Buildings on the terrain shape the wind where it enters the wildland–urban interface. GUST represents them by
immersed forcing with the form 7 wall law of ERF (``erf.if_wall_form = fraction_stress``) and Kynema-SGF
(``wall_form 7``). The ground keeps its body-fitted wall law on the terrain-following mesh and is not blanked.

**Status.** The geometry acts on momentum with SIMPLEC and z-line sweeps on one level (`Momentum`_). Roofs and
undersides have the wall law. Sides of buildings have fixed rows and no wall law yet, and :math:`k` and
:math:`\theta` have no building terms yet, so answers near buildings are not validated. The run refuses
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

**Mass fluxes.** A partial face's flux is scaled by its open fraction :math:`1 - \beta_f`. The metric itself is unchanged,
as the pressure gradient reads it.

**Pressure correction.** Blocked faces have zero area, a solid cell's row is :math:`\phi = 0`, and the mean of the singular
problem is over the fluid cells. GMRES runs on ERF's operator. Its FFT or multigrid preconditioner keeps every face
and does not know the mask, so the number of iterations grows with buildings (8 for one, 11 for nine, 3 for none, on the
test's mesh with ``gmres_mlmg``). Under multigrid GMRES is flexible, as the cycle's bottom solver is not a fixed
linear map. Fluid with no open path to the rest has its own null vector, so the run refuses it.

**Wall law** (``buildings.wall_form = fraction_stress``). A roof or underside cell has, neutrally,

.. math::

   u_*^2 = (1 - \beta) u_{*,\mathrm{own}}^2 + \beta u_{*,\mathrm{ref}}^2, \qquad
   u_{*,x} = \frac{\kappa U_x}{\ln(\max(d_x, 10 z_0) / z_0)},

with :math:`\kappa = 0.41`, :math:`z_0` = ``buildings.z0``, :math:`U` the horizontal speed of the cell (own) and of the next
one out (reference), and :math:`u_*^2 = u_{*,\mathrm{own}}^2` where :math:`d_2 = 0`. It enters as a sink on the horizontal velocity,

.. math::

   S = -r u, \qquad r = \frac{f\, u_*^2}{\Delta_z \max(|U_\mathrm{own}|, 10^{-3}\ \mathrm{m/s})},

with :math:`f` the blocked part of the wall's z-face and :math:`\Delta_z` the layer's height. A velocity row takes the mean
:math:`r` of its two cells. The viscous flux through the blocked part is closed, and :math:`r` goes into
:math:`a_P`, not into the local pseudo-time step. ``buildings.wall_form = drag`` keeps the fixed rows alone, with
viscous fluxes to the zero velocity, as Kynema-SGF's laminar pathway.

.. figure:: figures/buildings/block_flow.png
   :width: 100%

   One block in a log-law wind (``Tests/Canonical_RANS/Building_Cube``), with the fixed rows and the roof's wall law. The
   sides have no wall law yet and :math:`k` no building terms, so this shows the method working, not a validated flow.

Limits
------

- A domain face is not a wall: outside a non-periodic side the fields take the nearest domain column. The log warns of a
  building on such a side.
- :math:`d_1` uses the cell's :math:`\beta`. In a cell whose :math:`\beta` comes from another direction, such as an edge
  cell filled only partly in plan and in height, it is an estimate.
- A building thinner than a cell has its opposite sides in one cell, so its faces carry no net wall, and whether it blocks
  a face depends on where it sits against the mesh. The log warns; refine the mesh around it.
- On a slope the operator's cross terms and the pressure gradient's read the solid cells' :math:`\phi = 0` and pressure, and
  a roof's speed is horizontal, not along the roof. On flat ground neither arises.
- A partial cell's volume is not reduced: the sink and the fluxes are per full cell, which leaves the steady balance as it is.
- On strongly compressed terrain layers a face's fraction can differ a little from its cells', which can add a side bit
  at a roof.
