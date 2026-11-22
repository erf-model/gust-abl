Immersed terrain
================

By default the mesh follows the ground (:doc:`terrain_coordinates`). With
``terrain.representation = immersed_forcing`` the mesh is flat and the ground is a surface inside it. Cells below the
ground are solid, and the cells the ground cuts carry the form 7 wall law of ERF
(``erf.if_wall_form = fraction_stress``) and Kynema-SGF (``wall_form 7``), as the walls of :doc:`buildings` do. The
choice is made once, at start-up, and cannot change across a restart.

**Status.** A flat or sloping ground (``terrain.type`` = ``flat``, ``agnesi``, ``xyz`` or ``erf_grid``) in a neutral
run with SIMPLEC and z-line sweeps, laterally periodic or with inflow and outflow faces (`Open faces`_).
The roughness length is one value or a map (``bc.surface_z0_file``, `Roughness map`_). Stratification, buildings,
canopies and nests on an immersed ground are not supported yet; `Limits`_ lists what the run refuses. ``body_fitted`` remains the default and its answers are
unchanged.

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Input
     - Default
     - Effect
   * - ``terrain.representation``
     - ``body_fitted``
     - ``immersed_forcing``: a flat mesh with the ground immersed in it. ``hybrid``: the same ground in a mesh that
       follows its lower envelope (`Hybrid mesh`_).
   * - ``terrain.flat_height``
     - required (``hybrid``: the first layer's thickness)
     - The height of the terrain's zero level above the mesh's lower face (m): the ground itself with
       ``terrain.type = flat``.
   * - ``terrain.hybrid_half_width``
     - two cells (2 max(dx, dy))
     - ``hybrid`` only: the half-width (m) of the square over which the envelope takes its minimum and mean.
   * - ``terrain.small_volfrac``
     - ``0.005``
     - A solid fraction below it becomes 0, and one within it of 1 becomes 1 (≥ 1e-4).

The ranges are in :doc:`inputs`.

Geometry
--------

The ground's height is

.. math::

   h(x, y) = h_0 + h_t(x, y),

with :math:`h_0` = ``terrain.flat_height`` and :math:`h_t` the height ``terrain.type`` gives (0 for ``flat``; the
hill of ``agnesi``; the interpolated file of ``xyz`` or ``erf_grid``). It is evaluated at the mesh's nodes, as the
body-fitted ground is, and the mesh itself stays flat. Over a cell's footprint the ground is piecewise linear on the
four triangles that the footprint's centre (at the mean of the four corner heights) makes with each pair of
neighbouring corners. Along a cell's edge it is the straight line between the two nodes, so the two cells that share
a face see the same ground on it. A smooth hill is therefore a continuous surface, not a staircase of flat tiles
with side walls.

Every quantity is a mean over the footprint (or over a face's edge) of

.. math::

   M_n(t) = \overline{\,(t - h)_+^{\,n}\,}, \qquad n = 0, 1, 2,

for a level :math:`t`, which has a closed form on a triangle and on a segment (``Source/ImmersedSurface_K.H``).
:math:`M_0` is the part of the footprint whose ground lies below :math:`t`, :math:`M_1` the depth of fluid below
:math:`t`, and :math:`M_2` twice that fluid's first moment about the ground. For the layer between the levels
:math:`z_b` and :math:`z_t`, of height :math:`\Delta`:

.. math::

   \beta = 1 - \frac{M_1(z_t) - M_1(z_b)}{\Delta}, \qquad
   D = \frac{M_2(z_t) - M_2(z_b)}{2\,\bigl(M_1(z_t) - M_1(z_b)\bigr)} .

:math:`\beta` is the cell's solid fraction and :math:`D` the mean height of the cell's fluid above the ground. In a
cell wholly above the ground :math:`D` is the cell centre's height above the footprint's mean ground. A layer wholly
below the ground is solid (:math:`\beta = 1`). Velocities in solid cells are held at zero, and :math:`k` at its
floor. The blocked part of an x or y face is :math:`1 - (m_1(z_t) - m_1(z_b))/\Delta` with :math:`m_1` the same mean
along the face's edge, and that of a z face at the level :math:`z` is :math:`1 - M_0(z)`. Fractions within
``terrain.small_volfrac`` of 0 or 1 are rounded; a face next to a solid cell is blocked and one between two fluid
cells is open, as for buildings.

**Wall cells.** A cell is a wall cell where the blocked part of its low z face exceeds that of its high z face (or,
in a fluid cell, where its low face is blocked). The difference is the horizontal area of ground in the cell, as a
part of the footprint; over a column these areas add to 1. On a slope a column has several wall cells. The test has
no tolerance, unlike the buildings' wall cells, which need the difference to exceed ``buildings.small_volfrac``: on a
slope that rises through :math:`n` layers across a column the blocked part changes by about :math:`1/n` per layer,
and a tolerance above that drops most of the column's roof (with a tolerance of 0.1, a slope of 0.15 on 80 m
columns of 1 m layers keeps a third of its footprint as roof; unit test
``ImmersedSurface.RoofsAddToTheFootprintOnASlopeThroughManyLayers``). The ground has no side walls: the
differences of the x and y faces are the same surface, whose area the z faces already count. A vertical step in the terrain's data becomes a slope across one cell.

**Closure of the faces.** Before the rounding, the open parts of a cell's faces close: a uniform wind parallel to
a plane ground has no net flux through them, and the same holds for the four-triangle surface. The rounding of
the fractions (``terrain.small_volfrac``, and the rule that a face next to a solid cell is blocked) changes a face
by at most ``terrain.small_volfrac`` of its area in the cells it rounds.

**Distances.** The two distances of the wall law are the mean heights above the ground of the fluid of the wall cell
and of the cell above it (the reference cell),

.. math::

   d_1 = D_k, \qquad d_2 = D_{k+1} .

For a flat ground in a layer of height :math:`\Delta`, with :math:`\Delta_r` the height of the reference cell, they
are

.. math::

   d_1 = (1 - \beta)\,\frac{\Delta}{2}, \qquad d_2 = (1 - \beta)\,\Delta + \frac{\Delta_r}{2} .

When the ground lies on a cell face, the cell above it is the wall cell with :math:`\beta = 0` and
:math:`d_1 = \Delta / 2`. The closure's wall distance is :math:`D` in every cell (`Turbulence`_).

The ground must lie between the top of the mesh's first layer and the bottom of its top layer everywhere. The
plotfile fields ``building_beta`` and ``building_wall`` show the solid fraction and the wall cells
(:doc:`output`). The node heights in the plotfile are the flat mesh's.

**Cost.** A terrain-following mesh keeps its near-ground layer height over the whole relief. A flat mesh needs
that layer height from the lowest to the highest ground: :math:`H_r/\Delta` more layers for a relief :math:`H_r`,
when the layers are uniform. The ridge of `Verification`_ has the same layer height as its body-fitted reference
far from the ridge (66 against 64 layers); on the crest the body-fitted layers are 10 % thinner.

Wall law
--------

The friction velocity of the cut cell is

.. math::

   u_*^2 = (1 - \beta)\,u_{*1}^2 + \beta\,u_{*2}^2, \qquad
   u_{*n} = \frac{\kappa\,U_n}{\ln\bigl((d_n + z_0)/z_0\bigr)},

with :math:`U_1` and :math:`U_2` the speeds parallel to the ground in the cut cell and in the reference cell,
:math:`\kappa` the von Kármán constant and :math:`z_0` = ``bc.surface_z0``, or the column's value of the map
``bc.surface_z0_file`` (`Roughness map`_). This is the neutral form of the
body-fitted surface layer (:doc:`boundary_conditions`) at the two distances. The distances are exact: the
buildings' law floors them at 10 :math:`z_0`, which the ground's law does not, since
:math:`\ln((d + z_0)/z_0) > 0` for every :math:`d > 0`. With :math:`\beta = 0` the law is the
body-fitted neutral one read at the first cell's centre (the body-fitted default; ``bc.surface_zref`` and
``bc.surface_denominator_floor`` have no counterpart here).
The whole stress is a sink in the cut cell's momentum rows, as for a building's roof; the reference cell gives
its speed only.

.. _immersed-slope-forms:

Slope
-----

On a sloping ground the wall law reads the heights above the ground (:math:`d_1, d_2`) and the horizontal speed
:math:`\sqrt{u^2 + v^2}`, acts on the ground's horizontal area in the cell, and is a sink on :math:`u` and :math:`v`.
This is ERF's form. The ground has no side walls (`Geometry`_), and the closure's distance is :math:`D`
(`Turbulence`_).

**Forms tried and removed.** Three other forms were tried on a 2-D ridge (:math:`H` = 200 m, :math:`a` = 200 m)
against the body-fitted run on the same mesh, at 10 m above the ground. With the form above the rms of dS is 14.5 %
on 80 m columns of 15.6 m layers and 9.8 % on 160 m columns of 7.8 m layers. A wall law along the ground's normal
(distance along the normal, speed parallel to the ground, the slope's area) gave 14.1 % and 9.3 %. On the ridges and
the hill of `Verification`_ it differed from the form above by 0.3 % or less of the far-field speed. Side walls with
a wall law on the steps where the ground drops by more than a layer to the next column gave 14.1 % and 9.1 %, and
lowered the lee overspeed of the 160 m mesh from +9.7 to +5.9 %. The crest deficit of the 80 m mesh stayed (−25.0 %
without, −27.2 % with). A variant with no slip on the step faces was not measured on this ridge. With those side
walls, a closure distance from the column centre's ground, on footprints steeper than the layer, brought
:math:`\nu_t` beside the steps close to the body-fitted value and changed the rms by 0.3 points or less (9.1 to
8.8 % and 14.1 to 14.3 %). The three forms lowered the rms by 0.7 points or less and left the coarse mesh's crest
deficit, and were removed. Other changes (the distance along the ground's normal for the closure, the distance to
the nearest wall rectangle, the cell top's height, no shear production in wall cells beyond a distance, and
:math:`\nu_t` from the local speed and strain) changed the rms by 0.3 points or less, or raised it.

**Diffusion across the cut cell's top face.** The cut cell's speed stands for its fluid centroid, not its cell
centre. The viscous flux of :math:`u` and :math:`v` between the cut cell and the reference cell therefore uses the
distance between the two points where the speeds are known,

.. math::

   G \leftarrow s\,G, \qquad s = \frac{(\Delta + \Delta_r)/2}{d_2 - d_1} \ge 1,

in the residual and in the implicit rows; the diagonal gains the same amount, so the row sums are kept. With the
cell-centre spacing the cut cell came out slow and the cell above it fast, by up to 5 % of the log law. Buildings
keep the cell-centre spacing, and so do :math:`k` and :math:`\theta` across that face: in a log layer the flux of
:math:`k` there is zero.

**Body forces.** The pressure-gradient forcing, the Coriolis force and buoyancy act on the open part of each face's
control volume (the face's open fraction), so a solid layer takes no force and the cut cell takes the force on its
fluid. The steady balance of a column driven by a force :math:`F` is then :math:`u_*^2 = F D` with :math:`D` the
fluid's depth. Buildings keep the force on the whole control volume.

**Heat.** The ground takes no heat. The diffusive flux of :math:`\theta` through the blocked part of every face is
closed, in the residual and in the implicit rows (each coupling times the face's open part, the diagonal kept), as
for :math:`k`. A solid cell's :math:`\theta` therefore reaches no fluid cell, also through the side faces on a steep
slope where a column's ground drops by more than a layer to the next column (unit test
``ImmersedTerrain.ThetaIsNotConductedThroughTheGround``). Before this closure the solid cells conducted
:math:`\theta` as fluid cells do. A uniform ``bc.theta_profile`` is still required (`Limits`_), so the answers of the
runs that are accepted do not change; the closure is in place for a stratified immersed ground.

Turbulence
----------

The closure sees the ground as it sees a building with ``buildings.wall_tke = on`` (:doc:`buildings`, Turbulence):
the cut cell relaxes :math:`k` to the wall value :math:`u_*^2 / C_{\mu 0}^2`, the diffusive flux of :math:`k` into
the ground is closed, and the length scale takes the distance to the ground. That distance is :math:`D`, the mean
height of the cell's fluid above the ground: :math:`d_1` in a wall cell, as in the wall law, and the cell centre's
height above the ground under it elsewhere. It is the vertical height, as the body-fitted closure's. The initial
profile, the inflow profile, the Rayleigh layer and the convection blend (``convection.hybrid_height``) read heights
above the immersed ground too, through the same helper as on a body-fitted mesh. In a layer the ground cuts that
height is the centroid of the layer's open part, half the open height, at cell centres and at the x and y faces:
the inflow profile then gives a cut face the wind at the height the wall law reads it.

Where the ground under a footprint spans more than the cell's layer, :math:`D` averages over the footprint's deep
corner and can exceed the cell's height above the ground below its centre several times; a smaller distance there
did not lower the error (`Slope`_).

**Shear in the cut cell.** The closure's central difference in the cut cell reads the solid cell below it, where
the wind is 0. The vertical gradients of :math:`u` and :math:`v` there are taken from the wall law at the cell's
fluid centroid,

.. math::

   \frac{\partial (u, v)}{\partial z} = \frac{u_*}{\kappa\,(d_1 + z_0)}\,\frac{(u_c, v_c)}{|U_c|},

with :math:`u_*` the cut cell's friction velocity (`Wall law`_), :math:`d_1` its fluid centroid's height and
:math:`(u_c, v_c)` its wind. This is the gradient the body-fitted ground's first cell takes. Where the cut cell's
wind is exactly 0 it has no direction and the gradient is set to 0. The other terms of :math:`S^2` are unchanged. It applies with both values of ``turb.shear_stencil``.

**Shear above the cut cell.** In the first full cell above the cut cell the closure's central difference would place
the cut cell's speed at its cell centre. There the vertical gradients of :math:`u` and :math:`v` are taken as

.. math::

   \frac{\partial u}{\partial z} = R\,\frac{u_{k+1} - u_c}{(d_2 + z_0)\,\ln\dfrac{d_2 + h + z_0}{d_1 + z_0}}, \qquad
   R = \frac{(f_3 + z_0)\,\ln\dfrac{f_4 + z_0}{f_2 + z_0}}{2 \Delta_r},

with :math:`u_c` the cut cell's value, :math:`u_{k+1}` that of the cell above the reference cell, :math:`h` the
spacing between the reference cell's centre and the next, and :math:`f_2 = \Delta/2`,
:math:`f_3 = \Delta + \Delta_r/2`, :math:`f_4 = f_3 + h`. The first factor is the gradient of a log profile through
the two points, at the reference cell. :math:`R` is the ratio of the central difference to the true gradient of a log
profile for a ground on the cut cell's bottom face, so the cell makes the error the body-fitted ground's first
interior cell makes. With :math:`\beta = 0` (the ground on a layer's face) the expression is the closure's central difference
:math:`(u_{k+1} - u_c)/(2 \Delta_r)`, on stretched layers too. This follows the
rule of ERF's ``erf.if_flat_wall_dissipation`` and Kynema-SGF #2067 (an immersed wall repeats the flat ground's
discretisation), applied here to the shear.

With ``turb.shear_stencil = log_law`` (:ref:`log-shear`) this formula is not used: every fluid cell above the cut
cell takes the log-law gradient through its two neighbours, at the fluid centroids' distances, which is the
expression above with :math:`R = 1` in the first full cell. The cut cell takes the wall law's gradient (above)
and the solid cells keep the closure's :math:`S^2`. For a log profile the result is exact from the cut cell
upwards on flat ground.

Restart
-------

A checkpoint of an immersed run holds a file ``ImmersedTerrain`` with a hash of ``terrain.type`` and its keys (the
hill's inputs, or the terrain file's contents and its interpolation inputs), ``terrain.flat_height``,
``terrain.small_volfrac`` and the roughness: ``bc.surface_z0``, or with ``bc.surface_z0_file`` the SHA-256 of the
file's contents with ``terrain.idw_neighbors`` and ``terrain.idw_power`` (the file's path is not in the hash; a
map that cannot be read stops the run, with or without a stage). A
restart with another value of one of them, with a map in place of a value or a value in place of a map, or with the other ``terrain.representation``, stops naming ``io.restart``. The ground's heights are not in
the checkpoint: a restart computes them from the inputs, which the hash pins. The record's version is 3; a
checkpoint written before sloping grounds (version 1) or before the wall-cell test lost its tolerance (version 2) is
refused, with a message naming the two versions. A restart reads the ``bc`` inputs in a stage only, so a restart of
an immersed run names a ``steady.stage``.

Verification
------------

``Tests/Canonical_RANS/Immersed_Column`` is the 1-D neutral column of :doc:`test_cases` on 67 flat layers of
15.625 m, with the ground at :math:`(2 + f)` layers, so that :math:`f` is the cut cell's solid fraction. A force
:math:`F` drives it. The friction velocity recomputed from the plotfile with the law above must equal :math:`F D`,
and the speeds are compared with the log law of :math:`u_* = \sqrt{F D}`.

.. list-table::
   :header-rows: 1
   :widths: 8 14 22 22 17 17

   * - f
     - Iterations
     - :math:`u_*^2/(F D) - 1`
     - Cut cell vs log law at :math:`d_1`
     - Cell above
     - Two cells above
   * - 0
     - 2266
     - −5.3e-7
     - 0.00 %
     - −0.89 %
     - −1.33 %
   * - 0.3
     - 2284
     - −4.1e-7
     - +0.36 %
     - −0.85 %
     - −1.44 %
   * - 0.7
     - 2363
     - −2.0e-7
     - +1.77 %
     - −0.77 %
     - −1.75 %
   * - 0.9
     - 2509
     - −6.9e-8
     - +5.64 %
     - −0.65 %
     - −2.01 %

At :math:`f = 0` the column equals the body-fitted column of 65 layers within 4.2e-9 in :math:`u` and 1.0e-8 in
:math:`k`. The cut cell's excess over the log law follows from the law. With :math:`u_*^2` fixed at :math:`F D` and
:math:`e_1, e_2` the relative differences of the cut cell and the reference cell from the log law,
:math:`(1 - \beta)(1 + e_1)^2 + \beta (1 + e_2)^2 = 1`, so :math:`e_1 \approx -e_2\,\beta/(1 - \beta)`: the reference
cell's small deficit, of the size the body-fitted column has in its first cells, is multiplied by
:math:`\beta/(1 - \beta)` in the cut cell (9 at 0.9, 99 at 0.99).

The same columns with ``turb.shear_stencil = log_law`` (``rt_log_shear_column_*``). The deck has no length cap, so
the log law holds, and :math:`k` is compared with the local balance :math:`u_*^2 (1 - d/D)/C_{\mu 0}^2` (shear
production equal to dissipation at the stress :math:`u_*^2 (1 - d/D)`):

.. list-table::
   :header-rows: 1
   :widths: 8 12 18 20 21 21

   * - f
     - Iterations
     - :math:`u_*^2/(F D) - 1`
     - Cut cell vs log law at :math:`d_1`
     - Speed, three cells above
     - :math:`k`, three cells above
   * - 0
     - 2313
     - −5.5e-7
     - 0.00 %
     - −0.14 to −0.30 %
     - −0.5 to −1.6 %
   * - 0.3
     - 2339
     - −4.3e-7
     - +0.03 %
     - −0.07 to −0.23 %
     - −0.5 to −1.5 %
   * - 0.7
     - 2440
     - −2.1e-7
     - +0.01 %
     - 0.00 to −0.19 %
     - −0.2 to −1.3 %
   * - 0.9
     - 2606
     - −7.2e-8
     - −0.12 %
     - +0.01 to −0.18 %
     - +0.1 to −1.2 %

With ``central`` the same three cells have :math:`k` 13 to 22 %, 7 to 14 % and 4 to 7 % above the local balance.
With a capped length the cut cell's speed is lower than the log law with ``log_law`` (:ref:`log-shear`, Limits).

These numbers are for ``turb.face_mean = logarithmic``, which the decks set: with it and the spacing :math:`s`, the
flux of a log profile across the cut cell's top face is :math:`u_*^2` exactly. With the default ``arithmetic`` the
same flux is 1.10 :math:`u_*^2` at :math:`f = 0` and 1.43 :math:`u_*^2` at 0.9, and the column at 0.9 has its cut cell
24 % above the log law and the three cells above it 3 to 4 % below.


``Tests/Canonical_RANS/Immersed_Ridge_2D`` compares a sloping immersed ground with the same ground on a
terrain-following mesh, with the recommended turbulence setting (:doc:`turbulence_model`) in both. The wind is
read at fixed heights above the ground from the fluid cells. Differences are in percent of the body-fitted
far-field speed at that height (the crest speed: of the body-fitted crest speed).

.. list-table::
   :header-rows: 1
   :widths: 30 14 14 14 14 14

   * - Case
     - Iterations bf / immersed
     - Quantity
     - 20 m
     - 40 m
     - 80 m
   * - Wind along a ridge, slopes to 18°
     - 5298 / 3962
     - Largest difference
     - 1.3 %
     - 1.1 %
     - 0.8 %
   * - Wind across a ridge, slopes to 7°
     - 6297 / 4878
     - Crest speed
     - −2.2 %
     - −2.2 %
     - −2.4 %
   * -
     -
     - Lee minimum
     - −13.8 %
     - −7.9 %
     - −4.2 %
   * - Round hill, slopes to 12°
     - 2441 / 2180
     - Crest speed
     - −3.9 %
     - −2.9 %
     - −2.3 %
   * -
     -
     - Lee minimum
     - −24.4 %
     - −11.1 %
     - −2.1 %

With no flow over the slope (the wind along the ridge) the two representations agree within 1.3 %. About 1 % of
that is present on flat ground: the closure's central difference of a logarithmic profile over-estimates the
shear in the first cells above the wall cell, and more so above a cut cell, whose neighbour lies closer to the
ground (the column above; ``Tests/Canonical_RANS/Immersed_Column/README.md``). With flow over the ridge or the
hill the immersed runs decelerate more in the lee near the ground, and are 2 to 5 % slower in the far field of
the periodic domain (where the ridge's tail is still 3.7 m high). The cause of the lee difference is not isolated; the case's README lists what was not
tested. No growth of the GMRES iteration count was seen with the sloping ground, so the
Poisson preconditioner is unchanged.

Open faces
----------

An immersed ground may meet inflow and outflow faces. Three rules apply at the faces.

- An inflow face that lies in the ground (face fraction 1) carries no flow: its normal velocity is set to 0
  after the inflow profile is filled. A cut face carries the profile over its open part.
- The inflow profile's height is measured from the immersed ground. A point in the layer the ground cuts takes
  the centroid of the layer's open part, which is the height the wall law assigns to the cut cell's wind. The
  damping heights and the first guess (``init.first_guess = profile``) do the same. From this first guess the
  periodic column at :math:`f = 0.9` with ``turb.shear_stencil = log_law`` levels off with its :math:`u` and
  :math:`w` residuals near 2e-8 (each pressure solve leaves a uniform :math:`w` near 1e-10); its CTest solves the
  pressure to ``simplec.pressure_reltol = 1e-6``.
  From the layer's centre instead, the coarse Askervein pair stops in its first pressure solve.
- The outflow face needs nothing new: the solid cells next to it are at rest, and the mass balance of the two
  faces counts open areas.

``Tests/Canonical_RANS/Immersed_Open_Faces`` carries a log-law inflow over flat immersed ground to an outflow
face, 2560 m downstream, and compares the outflow column with the body-fitted channel at equal height above the
ground. With the default settings (``turb.shear_stencil = central``, ``terrain.small_volfrac`` = 0.005):

.. list-table::
   :header-rows: 1
   :widths: 16 16 16 16

   * - :math:`f`
     - iterations
     - :math:`u` (%)
     - :math:`k` (%)
   * - 0
     - 122
     - 0.10
     - 0.67
   * - 0.3
     - 122
     - 0.19
     - 3.35
   * - 0.5
     - 123
     - 0.33
     - 4.75
   * - 0.7
     - 124
     - 0.55
     - 3.97
   * - 0.9
     - 128
     - 1.25
     - 6.37
   * - 0.99
     - 139
     - 3.45
     - 25.4
   * - 1
     - 122
     - 0.00
     - 0.02

The body-fitted channel takes 122 iterations. The flux through the outflow face equals that through the inflow
face in every run. The differences are largest where the cut cell is thin; ``terrain.small_volfrac`` merges such
cells into the ground. CTest ``rt_immersed_open_faces`` runs the setting of the real-terrain cases
(``turb.shear_stencil = log_law``, ``terrain.small_volfrac = 0.1``) at :math:`f` = 0, 0.5, 0.95 and 1: the
differences are 0.8 % in :math:`u` and 5.5 % in :math:`k` at 0.5, and 0.1 % and 0.5 % at 0.95, where the cut
cell is merged into the ground. Its limits (1.5 % and 8 %) are these with a margin. The body-fitted column is
interpolated in :math:`\ln(z + z_0)`.

Four further CTests run the same channel and setting with one option of the open faces each, against the
body-fitted channel with the same option:

.. list-table::
   :header-rows: 1
   :widths: 30 34 18 18

   * - CTest
     - Option
     - :math:`f = 0.5`: :math:`u`, :math:`k` (%)
     - :math:`f = 0.95`: :math:`u`, :math:`k` (%)
   * - ``rt_immersed_open_faces_sponge``
     - ``bc.sponge_width = 320``, ``bc.sponge_strength = 0.05``
     - 0.28, 2.7
     - 0.07, 2.0
   * - ``rt_immersed_open_faces_power_law``
     - ``bc.inflow_profile = power_law`` (8 m/s at 10 m, :math:`\alpha` = 0.143)
     - 0.83, 5.2
     - 0.13, 0.5
   * - ``rt_immersed_open_faces_wall_velocity``
     - wall side faces moving at 8 m/s
     - 1.41, 8.7
     - 0.45, 4.6
   * - ``rt_immersed_open_faces_rayleigh_absolute``
     - ``bc.rayleigh_zdamp = 700`` on :math:`u, v, w, k`, absolute height
     - 0.78, 5.7
     - 0.13, 0.6

The limits are 1.5 % and 8 % (2 % and 10 % with the side walls). With the side walls the columns differ across
the channel, so the fluxes are the means over :math:`y`; they agree to round-off. The sponge relaxes towards the
inflow profile in solid cells too; those cells stay at rest, as the solid cells' equations hold them. With
``bc.rayleigh_height = absolute`` (the default) the height is measured from the terrain's zero, as on a mesh that
follows the terrain: the immersed mesh starts :math:`h_f` = ``terrain.flat_height`` below it, and that offset is
taken out (unit test ``ImmersedTerrain.AbsoluteRayleighHeightIsFromTheTerrainZero``).

Grid resolution
---------------

The horizontal velocity lives on the column faces. Where the ground rises by :math:`\Delta h` across a column, every
layer of that column below the downstream node's ground has its downstream face shut (the faces are exact for the
piecewise-linear ground, `Geometry`_), so the air must turn up inside one column. On a terrain-following mesh it
follows the slope. The immersed ground therefore needs columns narrow enough that :math:`\Delta h` is small near the
crest; the layer thickness matters much less.

``Tests/Canonical_RANS/Immersed_Ridge_Resolution`` measures this on a 2-D Witch of Agnesi ridge, 200 m high with
half-width 500 m (largest slope 0.26), periodic, with the measured cases' settings (``terrain.small_volfrac`` = 0.1,
``turb.shear_stencil = log_law``). The crest speed-up at 10 m above the ground, :math:`\Delta S_{10}`, is the speed
over the run's own far-field speed at that height, minus 1. Layers of 2 m:

.. list-table::
   :header-rows: 1
   :widths: 10 18 26 16 16

   * - dx (m)
     - :math:`\Delta h` per column (m)
     - Crest :math:`\Delta S_{10}`, terrain-following / immersed
     - Gap at 10 m
     - Gap at 40 m
   * - 64
     - 17
     - 1.133 / 0.652
     - −0.48
     - −0.10
   * - 32
     - 8
     - 1.148 / 0.857
     - −0.29
     - −0.05
   * - 16
     - 4
     - 1.159 / 1.032
     - −0.12
     - −0.02

The terrain-following value changes by 0.03 over the sweep. The immersed one shrinks its gap by a factor of 1.6 to
2.4 per halving of dx, about first order. Each run is divided by its own far field, and the immersed far field is
slower by 3.0, 2.6 and 1.6 % at dx = 64, 32 and 16 m; the crest speed at 10 m, immersed over terrain-following, is
0.751, 0.842 and 0.928. At dx = 64 m, layers of 2, 8, 16 and 32 m (aspect ratios dx/dz from 32 to 2) give the same
immersed crest value at 40 m, 0.71, and values from 0.51 to 0.65 at 10 m with no monotone trend. The immersed lee minimum is deeper than the
terrain-following one at every dx and converges more slowly than the crest. The extra :math:`k` above the crest
(1.87 against 0.97 m²/s² at dx = 64 m) follows the shear at the top of the slow layer, not the other way round:
capping the shear production there removed most of the extra :math:`k` and left the speed gap.

Askervein (wind from 210°, ``Tests/Real_Terrain/Immersed_Terrain_Cases``, dx = 23.4 m, largest slope 0.41 along the
wind, :math:`\Delta h` up to 9.7 m) shows the same pattern. At the hilltop, 10 m above the ground, the speed-up is
0.71 terrain-following and 0.48 immersed (gap −0.22) and :math:`k` is 2.0 times larger; at 40 m the gap is −0.04.
Averaged over the columns, the gap at 10 m is −0.12 on the windward slopes (slope above 0.15), −0.34 on the lee slopes
and −0.015 on the flat ground, and it correlates with the slope along the wind (+0.50 over the hill). The steepest
columns rise 9.7 m, about twice the guide below and close to the ridge at dx = 32 m (8 m per column, gap −0.29). With
2 m layers, the wall-cell rule without a tolerance (`Geometry`_) leaves this run unchanged (the same 70827 wall cells;
364 against 363 iterations; largest change in :math:`u` 0.002 m/s). The numbers come from
``Tests/Real_Terrain/Immersed_Terrain_Cases/compare_fields.py``.

**Guide.** Choose dx so that the ground rises by about 5 m or less across one column on the slopes that matter
(dx :math:`\lesssim` 5 m divided by the largest slope); on the ridge this leaves a crest gap of about −0.12 at
10 m. The layers can then be a few metres thick. Expect the lee to need finer columns than the crest. The guide
comes from one ridge, on which the rise per column and the columns per half-width change together.

Hybrid mesh
-----------

On the flat mesh the ground cuts many layers per column on a slope, and `Grid resolution`_ shows the crest error
that this leaves. With ``terrain.representation = hybrid`` the mesh follows a smooth surface :math:`h_s` that lies at
or below the terrain everywhere, and the ground is immersed in it as on the flat mesh:

.. math::

   m(x, y) = \min_{|x' - x| \le w,\ |y' - y| \le w} h_t(x', y'), \qquad
   h_s(x, y) = \min\Big(h_t(x, y),\ \overline{m}^{\,w}(x, y)\Big),

with :math:`h_t` the terrain's height at the nodes, :math:`w` = ``terrain.hybrid_half_width`` and
:math:`\overline{m}^{\,w}` the mean of :math:`m` over the same square. The square holds
:math:`\lfloor w/\Delta x \rfloor` nodes on each side of the node in x and :math:`\lfloor w/\Delta y \rfloor` in y.
It wraps across a periodic face and is cut at any other face. Every node of the square of a node lies in its own
square, so the mean of the minima is at most :math:`h_t` (the outer minimum only removes round-off). The mesh's
lower face is :math:`h_s`, and its levels follow it as a body-fitted mesh follows the ground
(:doc:`terrain_coordinates`; ``terrain.smoothing`` and ``terrain.gamma_m`` apply). The ground is at
:math:`h_0 + h_t` with :math:`h_0` = ``terrain.flat_height``, as on the flat mesh, so the gap

.. math::

   \delta(x, y) = h_0 + h_t(x, y) - h_s(x, y) \ \ge\ h_0

between the lower face and the ground is never smaller than :math:`h_0`. :math:`w = 0` gives :math:`h_s = h_t`, a
body-fitted mesh with the ground :math:`h_0` above its lower face. A flat :math:`h_s` would give the flat mesh's
ground; a wide square approaches it.

**Lower face, column by column.** :math:`h_0 \ge 0`. Where a column's first cell holds solid (:math:`\beta > 0`),
the lower face lies inside the terrain: the surface layer gives no stress there and no wall row for k, so the cell
is held at rest with k at ``turb.tke_floor`` (as every solid cell), and the immersed wall law acts at the ground
above it. Where the first cell holds none, which needs :math:`h_0 = 0` and :math:`h_s = h_t`, the lower face is
the ground and the surface layer acts on it as on a body-fitted mesh, along the mesh's slope. With :math:`h_0 > 0`
every column is of the first kind. With :math:`h_0 = 0` the second kind is where the envelope meets the terrain:
flat parts wider than the square, and next to faces that are not periodic, where the square is cut. The run prints
how many columns are of each kind (unit tests ``SurfaceLayer.ImmersedColumnsOfAHybridMesh`` and
``HybridTerrain.ImmersedColumnsWhereTheFirstCellIsCut``).

**Why.** The flat mesh's crest error grows with the rise of the ground across one column. On the hybrid mesh the
ground cuts the levels at the slope of :math:`\delta`, not of :math:`h_t`, so a narrow square leaves the immersed
ground cutting few layers per column. The gap is then small, so the layers through it can keep the first layer's
thickness: with ``mesh.zlevels_file`` uniform up to above the largest gap, every cut cell and the cells above it
have the same thickness, as on the flat mesh.

**What.** The envelope is computed once, at start-up, on every rank from the whole surface (one value per domain
node), so it does not depend on the boxes or the rank count (``Source/HybridTerrain.cpp``). The cut geometry
(`Geometry`_) measures the ground from each node's own height on a level. On a sloping level this integrates the
height above the ground along the level, and on a flat level it is the flat mesh's (unit tests
``ImmersedSurface.SlopedMeshMatchesTheShearedClosedForm`` and ``ImmersedSurface.HybridKernelsOnAFlatMeshMatchImmersed``).
The wall law, the shear of the first full cell and the closure use the mesh's metric, as on the flat mesh. The run
prints the gap's range and the height to which the layers keep the first one's thickness. The refusals are those of
the flat mesh (`Limits`_), except that ``terrain.smoothing`` and ``terrain.gamma_m`` are accepted. The ground must
lie at or below the bottom of each column's top layer. The restart hash
includes ``terrain.hybrid_half_width`` in metres, so a restart with another width is refused even where both give the
same node count. The square is cut at a face that is not periodic, so on ground that still slopes there the gap is
lopsided: on a ramp it is 0 at the low face and 1.5 times the interior gap at the high face. A width below one cell
gives no envelope (the mesh follows the terrain, as with 0), and flat terrain gives the flat immersed mesh; the run
and the deck checker warn of both.

**Pressure in the solid.** The pressure is 0 in the solid cells, as the Poisson solution leaves it there. On a
sloping level the x face's pressure gradient holds the slope term :math:`-z_x\,\partial p/\partial z`, the mean of
four vertical differences of :math:`p`, and at a cut cell the lower ones span the solid boundary. With
:math:`p \approx -17` Pa in the fluid and :math:`z_x \approx 0.3`, that difference gave a force of about 1 Pa/m,
some 30 times the real gradient, which turned the flow round in the cut cells of the windward slope and raised
:math:`k` above the crest. The momentum residual therefore takes the gradient of a copy of :math:`p` in which each
solid cell holds the value of the first cell above it that is not solid, in the valid columns, one ghost column on
each side and the ghost layer below the domain (``hybrid_solid_pressure``, unit test
``HybridTerrain.SolidPressureMakesTheSlopedGradientExact``). The slope term then takes no difference across the
boundary. On the flat mesh :math:`z_x = 0` and the copy is not made. The copy makes the gradient of a pressure
that varies only along the level exact. With a vertical gradient :math:`b = \partial p/\partial z` it leaves
about :math:`z_x b/4` to :math:`z_x b/2` at the first open faces above the solid, an error that does not shrink
with the layers' thickness; linear extrapolation from the two cells above would remove it. The Poisson operator's
cross terms still take the solid cells' 0 at cut cells, so the correction and the momentum gradient differ there
(the converged answer does not depend on it, the path to it may).

**Ridge.** The ridge of `Grid resolution`_ with an inflow face (log law, :math:`u_* = 0.35` m/s, :math:`z_0 = 0.1`
m), an outflow face, no Coriolis force, layers of 2 m and :math:`w = 128` m. The crest :math:`\Delta S_{10}` minus
the terrain-following value, and :math:`k` at 10 m over the terrain-following value:

.. list-table::
   :header-rows: 1
   :widths: 10 18 18 18 18 18

   * - dx (m)
     - Terrain-following :math:`\Delta S_{10}`
     - Immersed, :math:`h_0 = 2` m
     - **Hybrid**, :math:`h_0 = 2` m
     - Hybrid, :math:`h_0 = 0`
     - :math:`k_{10}` ratio, immersed / hybrid :math:`h_0 = 2` m
   * - 64
     - 0.946
     - −0.487
     - **−0.028**
     - −0.039
     - 2.20 / 1.05
   * - 32
     - 0.959
     - −0.294
     - **−0.011**
     - −0.015
     - 1.74 / 1.01
   * - 16
     - 0.945
     - −0.126
     - **0.000**
     - −0.005
     - 1.25 / 1.00

Without the pressure copy the hybrid values at dx = 64 and 32 m were −0.30 and −0.27, with :math:`k_{10}` 2.9
times the terrain-following value. The hybrid runs took 340 to 420 iterations to :math:`r < 10^{-4}`, the
terrain-following ones 440 to 920.

**Choosing w.** The same ridge (gentle: :math:`H` = 200 m, :math:`a` = 500 m, largest slope 0.26) and a steep one
(:math:`H` = 400 m, :math:`a` = 250 m, slope 1.0), :math:`h_0` = 2 m, against the terrain-following run on the same
mesh. The speed difference is counted at 5, 10, 20, 40, 80, 120, 160 and 200 m above the ground in every column;
"> 0.25" is the number of those points that differ by more than 0.25 m/s; :math:`\Delta S_{10}` is the crest speed
at 10 m minus the terrain-following one (7.6 m/s on the gentle ridge, 12.6 m/s on the steep one at dx = 32 m):

.. list-table::
   :header-rows: 1
   :widths: 16 12 16 16 14 14 12

   * - Ridge, dx (m)
     - :math:`w` (cells)
     - Largest gap / relief
     - Iterations (terrain-following)
     - > 0.25 m/s
     - :math:`\Delta S_{10}` (m/s)
     - Largest difference (m/s)
   * - gentle, 64
     - 0 / 1 / **2**
     - 0 / 0.08 / **0.16**
     - 823 / 569 / **381** (917)
     - 0 / 68 / **19** of 1024
     - +0.04 / −0.22 / **−0.11**
     - 0.07 / 0.88 / **0.61**
   * - gentle, 64
     - 4 / 8 / 16
     - 0.38 / 0.70 / 0.90
     - 896 / 1085 / 992
     - 346 / 503 / 523
     - −1.01 / −1.71 / −1.97
     - 2.6 / 3.6 / 3.7
   * - gentle, 32
     - 0 / 1 / 2 / **4**
     - 0 / 0.04 / 0.08 / **0.16**
     - 436 / 432 / 463 / **360** (710)
     - 0 / 0 / 43 / **1** of 2048
     - +0.05 / −0.02 / −0.15 / **−0.05**
     - 0.07 / 0.24 / 0.49 / **0.27**
   * - gentle, 32
     - 8 / 16
     - 0.38 / 0.69
     - 447 / 485
     - 317 / 644
     - −0.47 / −0.95
     - 1.3 / 2.7
   * - steep, 32
     - 0 / 1
     - 0 / 0.08
     - r near 1e-3 after 20000 / Poisson stop (6171)
     -
     -
     -
   * - steep, 32
     - **2** / 4 / 8 / 16
     - **0.16** / 0.38 / 0.70 / 0.90
     - **3459** / 2606 / 2122 / 1983
     - **67** / 1008 / 1017 / 445
     - **−0.66** / −4.93 / −5.51 / −3.62
     - **1.7** / 7.0 / 6.8 / 5.1

On the gentle ridge every width converges, and up to a gap of 0.16 of the relief (128 m at both spacings: the gap
follows :math:`w` in metres, not in cells) the speeds stay within about 0.25 m/s of terrain-following in fewer
iterations. On the steep ridge widths of 0 and 1 cell do not converge; 2 cells converge in about half the
terrain-following iterations with a crest 0.66 m/s low, and wider ones leave the crest 4 to 6 m/s low. At dx = 64 m
the steep ridge has about four columns per half-width; terrain-following there gives a crest 6 m/s below its
dx = 32 m value, and no width came within 0.25 m/s of it.

``tools/hybrid_width.py`` turns these two rules into a range for a given terrain and mesh, before a run. It builds
:math:`h_s` as the code does for every width up to ``--max-cells`` cells and reports the slope of :math:`h_s`, the
largest gap and the gap over the relief (max − min of :math:`h_t` after a least-squares plane is removed). The
lower end is 2 cells where the terrain's slope exceeds 0.5 and 0 otherwise. The upper end is the widest width whose
largest gap is at most 0.2 of the relief (between the 0.16 and 0.38 of the runs above). It suggests the lower end
when it lies within the range, and otherwise says that no width meets both and gives the spacing at which 2 cells
would fit (exit status 2). The upper end measures agreement with terrain-following, not with measurements: where
terrain-following is itself in error (a cliff, a separated lee) a larger difference need not be worse.

.. code-block:: bash

   python3 tools/hybrid_width.py --deck inputs --terrain-file terrain.xyz

The rules come from two analytic ridges only. A terrain file is interpolated linearly to the nodes (GUST uses
inverse-distance weights) and ``terrain.smoothing`` is not applied.

**Defaults.** Without the keys, ``terrain.flat_height`` is the thickness of the mesh's first layer and
``terrain.hybrid_half_width`` is 2 cells (twice the larger of dx and dy; a direction of one cell does not count).
No single width suits every terrain: 2 cells is the narrowest that converged on the steep ridge, and on the gentle
one it stays within 0.61 m/s of terrain-following. The log names each default it sets. On the measured cases, at
these defaults, against terrain-following on the same mesh (rms error of the speed-up against the measurements, or
of :math:`U` against the tunnel):

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Case
     - Terrain-following → hybrid
     - Iterations (terrain-following / hybrid)
   * - Askervein 210°, lines A / AA / B
     - 0.122 / 0.081 / 0.120 → 0.155 / 0.146 / 0.155
     - 202 / 235
   * - Bolund 270°, line B
     - 0.113 → 0.024
     - 189 / 189
   * - Bolund 239°, line A
     - 0.075 → 0.124
     - 209 / 209
   * - Bolund 90°, line B
     - 0.164 → 0.180
     - 222 / 394
   * - Ishihara smooth ridge, :math:`U` (m/s)
     - 0.27 → 0.17
     - 537 / 619
   * - Ishihara rough ridge, :math:`U` (m/s)
     - 0.57 → 0.55
     - 655 / 664
   * - Big Southern Butte 222°, speed at the 40 sensors at 3 m (m/s)
     - 3.20 → 5.94
     - 410 / 1094

The hybrid mesh is closer to the measurements behind the Bolund cliff and on the smooth ridge, and further from
them on Askervein, the other two Bolund directions and Big Southern Butte. Most masts differ from terrain-following by more than
0.25 m/s at these defaults, so the representation changes the answer on these sites.

On Big Southern Butte (dx = 137 m) the default width gives a largest gap of 294 m on 780 m of relief (0.37), and
the hybrid run is close to the flat immersed one (rms 6.14 m/s): the speed at the two sensors on the top is 10.8 m/s
low, against 3.2 m/s for terrain-following. ``tools/hybrid_width.py`` reports that no width meets both rules on this
mesh. The run reached :math:`r < 10^{-4}` after 570 iterations (terrain-following 260) and stopped on
``steady.residual_tol`` after 1094, in 58 times the terrain-following run's time on 3 ranks (114 layers against 40; part of the run shared the
machine with other jobs).

The 2-cell default comes from 2-D ridges with a half-width of at least 4 cells. A round hill has a larger gap for the
same width (about 0.26 of the relief at a half-width of 8 cells, 0.58 at 4), and on narrow features the default
lies above the tool's upper end. The run warns when the largest gap exceeds 0.2 of the ground's height range and
names ``tools/hybrid_width.py``. The gap is about the slope times the width while the square is narrower than the
feature; once it spans a crest the crest gap is up to 1.2 times that.

Roughness map
-------------

``bc.surface_z0_file`` (:ref:`boundary conditions <roughness-map>`) gives each ground column its own :math:`z_0`, for a surface such as
water beside land. An immersed ground and the ground of a hybrid mesh read the same map as the body-fitted ground:
the map is interpolated once, at start-up, to the centre of each column, and a column's value is that of its cut
cell wherever the ground cuts the column. The value enters:

- the wall law of the cut cell (`Wall law`_), and through its :math:`u_*` the wall value of :math:`k`;
- the strain of the cut cell and the corrected strain of the cell above it (`Turbulence`_);
- with ``turb.shear_stencil = log_law``, the log-law gradient of every cell above the wall cell;
- the closure's length :math:`\kappa (d + z_0)`, as in a body-fitted run;
- on a hybrid mesh, the surface layer of the columns whose lower face is the ground, as in a body-fitted run.

The log names the map as the wall law's :math:`z_0`. ``terrain.idw_neighbors`` and ``terrain.idw_power`` apply to
the map as they do in a body-fitted run.

**Roughness step.** ``Tests/Canonical_RANS/Immersed_Roughness_Step`` is the channel of `Open faces`_ with
:math:`z_0` = 0.03 m up to :math:`x` = 1120 m and 0.3 m beyond. The step lies inside a box, so a box holds both
values. The map changes :math:`u` of the body-fitted run by up to 24.3 %. Each ground is compared with the
body-fitted run with the same map, at equal height above the ground, in percent of the body-fitted value
(``turb.shear_stencil = log_law``, ``terrain.small_volfrac = 0.1``; :math:`F` is the solid fraction of the cut
cell, :math:`H` the ground's height in the first layer of the hybrid mesh):

.. list-table::
   :header-rows: 1
   :widths: 25 15 15 15 15 15

   * - Ground
     - Iterations
     - :math:`u`
     - :math:`k`, away from the step
     - :math:`k`, first 400 m after it
     - :math:`k`, wall cell
   * - body-fitted
     - 121
     - —
     - —
     - —
     - —
   * - immersed, :math:`F = 0`
     - 120
     - 0.07
     - 0.62
     - 0.47
     - 0.21
   * - immersed, :math:`F = 0.3`
     - 121
     - 0.96
     - 8.3
     - 27.1
     - 4.6
   * - immersed, :math:`F = 0.7`
     - 126
     - 2.81
     - 14.0
     - 51.6
     - 10.2
   * - hybrid, :math:`H = 0`
     - 121
     - 0.00
     - 0.00
     - 0.00
     - 0.00
   * - hybrid, :math:`H = 0.5`
     - 125
     - 2.04
     - 12.5
     - 40.8
     - 5.9

A ground on a face (:math:`F = 0`, :math:`H = 0`) reproduces the body-fitted run. A ground that cuts a cell
differs as it does on a uniform ground (`Open faces`_), and more just downstream of the step: there the new
internal layer grows from less than one layer to a few layers, :math:`k` changes across it within one or two
cells, and a cut cell's centre samples it at another height than the body-fitted centres. The limits bound the
cut cell's known difference, not the map. The test ``rt_immersed_roughness_step`` bounds the
first, second and fourth columns at 4, 18 and 13 %; the third is printed only.

**Bolund.** The decks of ``Tests/Real_Terrain/Immersed_Terrain_Cases`` for Bolund 270° and 239°, with the map of
``Tests/Real_Terrain/WindNinja2019/bolund_z0_map.py`` (0.0003 m over water, 0.015 m over land; see the
:ref:`boundary conditions <roughness-map>`) in place of the uniform 0.01 m, and the inflow log law of the water
(``bc.inflow_z0 = 0.0003``, :math:`u_*` for the same reference wind). The hybrid mesh is at its defaults. The
number is the rms error of the speed-up against the masts at 5 m (line B, 4 masts, for 270°; line A, 5 masts, for
239°), with the iterations in brackets:

.. list-table::
   :header-rows: 1
   :widths: 22 26 26 26

   * - Case
     - Terrain-following
     - Immersed
     - Hybrid
   * - 270°, uniform :math:`z_0`
     - 0.113 (189)
     - 0.032 (213)
     - 0.024 (189)
   * - 270°, map
     - 0.088 (157)
     - 0.064 (360)
     - 0.061 (156)
   * - 239°, uniform :math:`z_0`
     - 0.075 (209)
     - 0.133 (210)
     - 0.124 (209)
   * - 239°, map
     - 0.079 (174)
     - 0.129 (188)
     - 0.126 (172)

All six runs with the map reach the steady verdict. With the map the three grounds are closer to each other for
270° than with the uniform value: the terrain-following run moves towards the masts and the immersed and hybrid
runs away from them (the mast behind the cliff, M8: measured -0.56, the three runs -0.46, -0.64 and -0.62). For
239° the map changes each number by at most 0.004. The immersed run for 270° takes 360 iterations: its residual
falls to 2e-4, rises to 9e-4 between iterations 240 and 270, and then falls to the verdict. The uniform rows are
those of ``Tests/Real_Terrain/Immersed_Terrain_Cases/README.md`` and of `Hybrid mesh`_. The hybrid runs with the
map were sampled with another script than those rows; it gives 0.032 (188 iterations) for the hybrid run for 270°
with the uniform value, where the table has 0.024. The runs with the map are not a test of the
repository: the Bolund terrain and masts are not redistributed (``Tests/Real_Terrain/WindNinja2019/README.md``
names the data and its providers).

Limits
------

- The mesh must have one hierarchy level. On a periodic mesh a terrain that is not periodic (a file, a hill whose
  tails have not decayed at the sides) makes a step across the seam, which the geometry treats as a slope across
  one cell.
- ``bc.zlo = surface_layer``, ``bc.surface_stability = neutral``, ``bc.surface_wall_k = dirichlet``, a uniform
  ``bc.theta_profile``, ``turb.closure = kl`` and ``steady.sweep = zjacobi`` are
  required. The refusals in this list apply to a run with no stage too (the sweep and the Poisson backend in a
  stage only).
- ``buildings.type``, nests, ``init.first_guess = column``, ``steady.probes_file``,
  ``steady.stage = simplec_newton`` and ``steady.stage = unsteady`` are refused. The column's and
  the probes' heights are measured from the mesh's lower face, not from the immersed ground.
- A forest canopy (``canopy.type``) stands on the immersed ground (:ref:`the drag field of a cut cell <canopy-immersed-ground>`).
  ``canopy.heat = fixed`` is refused: the run is neutral.
- ``bc.surface_zref`` and ``bc.surface_denominator_floor`` are refused when given: the immersed wall law reads the
  cut cell and the cell above it at their exact distances and has no floor.
- The ``ustar`` field of the turbulence plotfile and the ``mean_ustar`` monitor are those of the mesh's lower
  face, which lies in solid cells (0.0094 m/s in the column at :math:`f = 0.5`, where the ground's :math:`u_*` is
  0.351 m/s). The ground's :math:`u_*` follows from the cut cell's speed and the wall law above; an output for it
  is not in this version.
- ``solver.poisson_backend = mlmg`` is refused (its operator has no solid cells). For a 1-D column use
  ``gmres_fft``: with ``gmres_mlmg`` the pressure of a column differs between columns by round-off, :math:`v` and
  :math:`w` grow to about 1e-9, and a column that starts with :math:`v = w = 0` does not meet the stop rule. The
  mesh of an immersed ground is flat, so the FFT-preconditioned solve always applies.
- ``terrain.smoothing``, ``terrain.gamma_m``, ``terrain.hierarchy_passes`` and ``terrain.nest_blend_width`` are
  refused when given: nothing follows or smooths a ground.
- ``terrain.flat_height`` and the whole ground must lie between the top of the mesh's first layer and the bottom of
  its top layer: a cut cell may not be the body-fitted ground's first cell, whose :math:`k` row would share the
  two wall values, and the law needs a cell above the cut cell.
- On a slope the reference cell of a wall cell is the cell above it, which can be a wall cell itself. Where that
  cell is solid (a slope steeper than the layers resolve), :math:`d_2 = 0` and the law uses the wall cell alone.
- ``bc.outflow_transport = from_velocity`` is refused: at the outflow face it slows the cut cell to 0.71 of the
  log law at :math:`f = 0.5` and lifts :math:`w` to 0.19 m/s in the outflow column. ``mirror`` (the default)
  applies. The first-guess field of a run with no stage is written before the solid inflow faces are set to 0.
- A thin cut cell (:math:`f` above about 0.9) moves faster than the log law (`Wall law`_), which lowers the shear
  and :math:`k` of the cell above it: 0.27 against 0.34 m\ :sup:`2`/s\ :sup:`2` at :math:`f = 0.9` and 0.21 at
  0.99 with ``turb.shear_stencil = log_law`` in the open-face channel. ``terrain.small_volfrac = 0.1`` merges
  these cells into the ground and removes the dip (`Open faces`_).
- The cut cell's shear reads the cell's plain face means of :math:`u` and :math:`v` at the vertical distance
  :math:`d_1`. The wall law of the momentum equations weights the faces by their open part. Next to a blocked
  side face the two differ.
- Stretched layers (``mesh.stretch_ratio``) are accepted; the column is verified on uniform layers only.
- ``turb.face_mean``: the column decks set ``logarithmic``, for which the flux of a log profile across the cut
  cell's top face is exact (`Verification`_). The ridge and hill decks use the recommended turbulence setting and
  keep the default ``arithmetic`` in both runs: with that setting the force-driven column's first full cell is
  1.8 % above the log law with ``arithmetic`` and 4.9 % with ``logarithmic`` (body-fitted: 0.4 % and 2.1 %).
- A cut cell that is nearly solid: with :math:`f = 0.98` and 0.99 the column above does not reach the steady verdict
  in 6000 iterations with the deck's settings, with either backend; :math:`f = 0.95` does (2612 iterations). At
  0.99 it reaches it with ``simplec.pressure_reltol = 1e-8`` (2863 iterations) or ``steady.cfl_max = 10`` (3175).
  The cause is not identified. The cut cell is then 43 % above the log law (the factor :math:`\beta/(1 - \beta)`
  above) and the cells above it within 2.6 %. ``terrain.small_volfrac`` (default 0.005) rounds only fractions above
  0.995 to solid, so such cells can occur.
- The wall value of :math:`k` relaxes at the rate :math:`u_*/d_1`, so the :math:`k` residual of a nearly solid cut
  cell starts large (308 at :math:`f = 0.99` against 2.6 at 0.5) and its drop is measured from there.
- The column is verified with :math:`z_0` = 0.1 m on 15.6 m layers. With :math:`z_0` = 1 m the cells above the cut
  cell are 1.3 to 2.9 % below the log law at :math:`f = 0.7`.
- Against the terrain-following mesh the wind over a ridge and a round hill is 2 to 5 % lower, and lower still in
  the lee near the ground (`Verification`_); the cause is not isolated. Not run on an immersed ground: a terrain
  from a file in the application (the unit tests fill and hash one), stretched layers on a slope, and slopes
  above 18°.
