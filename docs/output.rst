Output files
============

A run writes up to three kinds of output at the end: a plotfile with the solution, a plotfile with the turbulence
fields, and a checkpoint for restarting. While it runs it prints the iteration log.

Iteration log
-------------

The solver prints a row for the first ten iterations, then every tenth, and the last; in the Newton stage every
step. With ``io.iteration_log = table`` (the default) the rows are aligned columns under a header that is repeated
every 20 rows:

.. code-block:: text

     Iter        CFL  Line solver      u-res      v-res      w-res    tke-res  theta-res      p-res          r
   ------  ---------  -----------  ---------  ---------  ---------  ---------  ---------  ---------  ---------
        1  1.000e+00  zjacobi      4.997e-03  0.000e+00  0.000e+00  9.271e-02  0.000e+00  0.000e+00  1.000e+00
        2  1.000e+00  zjacobi      1.453e-03  0.000e+00  0.000e+00  3.119e-02  0.000e+00  0.000e+00  1.000e+00

``Iter`` is the iteration, ``CFL`` the pseudo-time CFL number and ``Line solver`` the z-line solver of the sweeps
(``steady.sweep``). The six residual columns are the scaled residual norms of :math:`u, v, w, k, \theta` and
:math:`p`. ``r`` is the largest residual over its reference (1 during the reference hold); it drives the CFL ramp
and the stage switch, while the steady criterion tests each residual against ``steady.residual_tol`` and
``steady.residual_drop_orders``. On the coarse levels of grid sequencing each row starts with ``level l:``. A
Newton row names ``newton`` as its line solver, gives the continuation step's CFL (``-`` when it took no step) and
adds the step's outcome after ``r`` (the residual before and after, the line-search factor, the attempts and the
GMRES iterations, ``FAILED`` for a failed step). With ``steady.confidence`` a SIMPLEC row adds after ``r`` the
share of the cells that meet the residual criterion and the height of the highest one that does not.

``io.iteration_log = full`` prints the log of earlier versions: one ``gust: it`` line per row with the CFL, ``r``, the
residuals in the order :math:`u, v, w, p, \theta, k`, the monitors, the Poisson iterations, the divergence before
and after the correction, the clipped :math:`k` cells and, with ``steady.confidence``, the confidence share and height.

Plotfile
--------

``<io.plot_file>NNNNN``, where ``NNNNN`` is the step, is an AMReX plotfile with the cell-centred variables below.
It can be read with yt, ParaView, VisIt or AMReX's ``amrex_fcompare``.

.. list-table::
   :header-rows: 1
   :widths: 15 85

   * - Variable
     - Meaning
   * - ``u``, ``v``, ``w``
     - Velocity components: the mean of the two face values around the cell.
   * - ``p``
     - Perturbation pressure p′.
   * - ``theta``
     - Potential temperature.
   * - ``tke``
     - Turbulent kinetic energy k.
   * - ``z``
     - Height of the cell centre: the mean of the cell's eight node heights.
   * - ``detJ``
     - The cell volume factor, (cell volume) / (Δx Δy Δζ).
   * - ``canopy_drag``
     - With ``canopy.type`` other than ``none``: the canopy's drag field :math:`f = C_d a` (1/m), the mean over the
       cell (:doc:`canopy`).
   * - ``canopy_lad``
     - With ``canopy.heat = fixed``: the leaf area density :math:`a` (1/m) (:doc:`canopy`).
   * - ``building_beta``
     - With ``buildings.type`` other than ``none``: the cell's solid fraction (:doc:`buildings`).
   * - ``building_wall``
     - With ``buildings.type`` other than ``none``: the sum of the bits of the cell's walls (:doc:`buildings`).

The node heights themselves, ghost nodes included, are written to ``Level_0/z_phys_nd`` inside the plotfile as a
VisMF MultiFab.

Grid hierarchy plotfiles
------------------------

With ``mesh.hierarchy_levels`` = L > 1 and ``io.write_plotfile = 1``, every coarse level l = 1 … L − 1 is written next
to the plotfile as ``<io.plot_file>_level<l>_NNNNN``: the cell variables ``z`` and ``detJ`` of that level, and its
node heights in ``Level_0/z_phys_nd``. The run log has one mesh-quality report per level and a summary table
(:doc:`grid_hierarchy`).

Turbulence diagnostics plotfile
-------------------------------

With ``io.turb_plot_file`` set and the SIMPLEC stage on, ``<io.turb_plot_file>NNNNN`` holds the closure's fields at
the final state under ERF's variable names, so ERF's Canonical_RANS check scripts read it unchanged
(``Source/TurbulenceOutput.H``): ``x_velocity``, ``y_velocity``, ``z_velocity`` (cell means of the face velocities),
``theta``, ``KE`` (k), ``Kmv`` (:math:`\nu_t`; the density is 1), ``Khv`` (:math:`\gamma_\theta` less the molecular
diffusivity), ``Lturb``, ``walldist`` (the closure's wall distance; with ``buildings.wall_tke = on`` the smaller of the ground's and the nearest building wall's, the latter plus ``buildings.z0`` less the ground's roughness, which can be below zero, :doc:`buildings`), ``diss`` (:math:`\varepsilon`), ``density`` (1), ``Rt``, ``cmu`` and
``cmu_prime`` (the stability functions of :math:`R_t`, AL01 Eqs. 31–32, as ERF writes them), ``ustar`` (the ground
cells; 0 above), ``z_phys`` (the cell-centre height), and the residual at the final state, ``R_u``, ``R_v``, ``R_w`` (cell means of the face
residuals), ``R_p``, ``R_theta``, ``R_tke``: where the iteration has not converged. Inside it ``surf_hist.dat`` holds ``t u_star`` and one row with the iteration count and the
mean :math:`u_*` over the ground, the file ERF's ``--physics`` checks read.

Checkpoint
----------

``<io.chk_file>NNNNN`` holds a text ``Header`` and every field with its ghost cells: two for the solution fields,
one for the node heights. A solution field's ghost cells that stand for another box's cells (or a periodic image)
hold that box's valid values, refreshed as the file is written; the running state is not changed. (Older
checkpoints with one ghost layer restart unchanged, since only valid values are read back.)

The ``Header`` starts with ``gust-checkpoint`` and lists the format version, the bytes per real (8 or 4), the step, the mesh (``n_cell``, and
``prob_lo`` and ``prob_hi`` in hexadecimal floating point so they read back exactly), ``is_periodic`` and the
BoxArray. The fields are ``u v w`` (faces), ``p theta tke`` (cell centres) and ``z_phys_nd`` (nodes).

A restart restores every valid value exactly. The node heights are rebuilt from the mesh and terrain inputs and
compared with the checkpoint's: a difference beyond 64 ulp of the top level's height stops the run, naming
``io.restart`` (a changed mesh or terrain input); a smaller one, as from another compiler that contracts the mesh
build's multiply-adds differently, is reported in the log, and the run continues on the rebuilt heights.
The restart fills the other ghost cells before the first step, whose local :math:`\Delta\tau` reads them. With the SIMPLEC stage the checkpoint also holds the driver history
(``Driver``, :doc:`solution_algorithm`) and the closure's :math:`\nu_t` history (``Level_0/closure_nut``, :doc:`turbulence_model`).

Checkpoints, driver and stage histories and reference archives written before the rename to GUST (October 2026)
start with ``ablsteady-<kind>`` instead of ``gust-<kind>``. They are read the same way and restart unchanged.
