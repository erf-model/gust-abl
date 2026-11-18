Immersed terrain
================

By default the mesh follows the ground (:doc:`terrain_coordinates`). With
``terrain.representation = immersed_forcing`` the mesh is flat and the ground is a surface inside it. Cells below the
ground are solid, and the cells the ground cuts carry the form 7 wall law of ERF
(``erf.if_wall_form = fraction_stress``) and Kynema-SGF (``wall_form 7``), as the walls of :doc:`buildings` do. The
choice is made once, at start-up, and cannot change across a restart.

**Status.** A flat ground at ``terrain.flat_height`` in a neutral, laterally periodic run with SIMPLEC and z-line
sweeps. Hills, open boundaries, stratification, roughness maps, buildings, canopies and nests on an immersed ground
are not supported yet; `Limits`_ lists what the run refuses. ``body_fitted`` remains the default and its answers are
unchanged.

.. list-table::
   :header-rows: 1
   :widths: 28 22 50

   * - Input
     - Default
     - Effect
   * - ``terrain.representation``
     - ``body_fitted``
     - ``immersed_forcing``: a flat mesh with the ground immersed in it.
   * - ``terrain.flat_height``
     - required
     - The ground's height above the mesh's lower face (m).
   * - ``terrain.small_volfrac``
     - ``0.005``
     - A solid fraction below it becomes 0, and one within it of 1 becomes 1 (≥ 1e-4).

The ranges are in :doc:`inputs`.

Geometry
--------

The ground is one tile of height ``terrain.flat_height`` covering the domain, handled as a building's height map
(:doc:`buildings`, Solid fractions and Wall cells). A layer wholly below the ground is solid
(:math:`\beta = 1`), one wholly above is fluid (:math:`\beta = 0`), and the layer the ground lies in is the cut cell,
with solid fraction :math:`\beta`. Velocities in solid cells are held at zero, and :math:`k` at its floor. With
:math:`\Delta` the cut cell's height and :math:`\Delta_r` that of the cell above it (the reference cell), the two
distances of the wall law are

.. math::

   d_1 = (1 - \beta)\,\frac{\Delta}{2}, \qquad d_2 = (1 - \beta)\,\Delta + \frac{\Delta_r}{2}:

the heights above the ground of the cut cell's fluid centroid and of the reference cell's centre. When the ground
lies on a cell face, the cell above it is the wall cell with :math:`\beta = 0` and :math:`d_1 = \Delta / 2`.

The plotfile fields ``building_beta`` and ``building_wall`` show the solid fraction and the wall cells
(:doc:`output`).

Wall law
--------

The friction velocity of the cut cell is

.. math::

   u_*^2 = (1 - \beta)\,u_{*1}^2 + \beta\,u_{*2}^2, \qquad
   u_{*n} = \frac{\kappa\,U_n}{\ln\bigl((d_n + z_0)/z_0\bigr)},

with :math:`U_1` and :math:`U_2` the speeds parallel to the ground in the cut cell and in the reference cell,
:math:`\kappa` the von Kármán constant and :math:`z_0` = ``bc.surface_z0``. This is the neutral form of the
body-fitted surface layer (:doc:`boundary_conditions`) at the two distances. The distances are exact: the
buildings' law floors them at 10 :math:`z_0`, which the ground's law does not, since
:math:`\ln((d + z_0)/z_0) > 0` for every :math:`d > 0`. With :math:`\beta = 0` the law is the
body-fitted neutral one read at the first cell's centre (the body-fitted default; ``bc.surface_zref`` and
``bc.surface_denominator_floor`` have no counterpart here).
The whole stress is a sink in the cut cell's momentum rows, as for a building's roof; the reference cell gives
its speed only.

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

Turbulence
----------

The closure sees the ground as it sees a building with ``buildings.wall_tke = on`` (:doc:`buildings`, Turbulence):
the cut cell relaxes :math:`k` to the wall value :math:`u_*^2 / C_{\mu 0}^2`, the diffusive flux of :math:`k` into
the ground is closed, and the length scale takes the distance to the ground. The wall cell's distance is
:math:`d_1`, exact, as in the wall law.

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
expression above with :math:`R = 1` in the first full cell. The cut cell and the solid cells keep the closure's
:math:`S^2`.

Restart
-------

A checkpoint of an immersed run holds a file ``ImmersedTerrain`` with a hash of ``terrain.flat_height``,
``terrain.small_volfrac`` and ``bc.surface_z0``. A restart with another value of one of them, or with the other
``terrain.representation``, stops naming ``io.restart``. A restart reads the ``bc`` inputs in a stage only, so a
restart of an immersed run names a ``steady.stage``.

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

Limits
------

- ``terrain.type`` must be ``flat``; the mesh must be periodic in x and y (``mesh.is_periodic = 1 1 0``) and have
  one hierarchy level.
- ``bc.zlo = surface_layer``, ``bc.surface_stability = neutral``, ``bc.surface_wall_k = dirichlet``, a uniform
  ``bc.surface_z0``, a uniform ``bc.theta_profile``, ``turb.closure = kl`` and ``steady.sweep = zjacobi`` are
  required. The refusals in this list apply to a run with no stage too (the sweep and the Poisson backend in a
  stage only).
- ``buildings.type``, ``canopy.type``, nests, ``init.first_guess = column``, ``steady.probes_file``,
  ``steady.stage = simplec_newton`` and ``steady.stage = unsteady`` are refused. The canopy's, the column's and
  the probes' heights are measured from the mesh's lower face, not from the immersed ground.
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
- ``terrain.flat_height`` must lie between the top of the mesh's first layer and the bottom of its top layer: the
  cut cell may not be the body-fitted ground's first cell, whose :math:`k` row would share the two wall values,
  and the law needs a cell above the cut cell.
- Stretched layers (``mesh.stretch_ratio``) are accepted; the column is verified on uniform layers only.
- Use ``turb.face_mean = logarithmic`` (`Verification`_).
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
- The answers are verified on the column only. No hill has been run on an immersed ground.
