Input reference
===============

This page lists every input of GUST. The first table is a quick reference; the sections after it give the
allowed range and meaning of each input, grouped by namespace.

The same keys, with their defaults and the text of the detailed tables, are in ``Tests/inputs_reference``: a
commented file to copy lines from into a deck. It is generated from this page (``tools/input_reference.py
--write``); the CTest ``input_reference_current`` fails when the two differ, or when this page and the keys read in
``Source/`` disagree.

Input file format
-----------------

An input file has one input per line, ``namespace.key = value``. Lists are separated by spaces, and ``#`` starts a
comment:

.. code-block:: text

   mesh.n_cell  = 64 64 32     # cells in x, y, z
   bc.zlo       = surface_layer

Inputs can also be given on the command line after the input file (:doc:`running`).

The run stops before doing any work, naming the input, when:

- a key is misspelled or not used by this version;
- a key's namespace is unknown (the AMReX namespaces ``amrex``, ``fab``, ``fabarray``, ``vismf``, ``tiny_profiler``,
  ``DistributionMapping``, ``device`` and ``hypre`` are allowed);
- a key is given twice, or a single-value key is given several values;
- a value is outside its range;
- a key has no effect with the other inputs (for example ``terrain.hill_height`` with ``terrain.type = flat``).

All errors are reported together.

Quick reference
---------------

.. list-table::
   :header-rows: 1
   :widths: 34 20 46

   * - Input
     - Default
     - Description
   * - **Mesh** (`mesh`_)
     -
     -
   * - ``mesh.n_cell``
     - required
     - Number of cells in x, y, z
   * - ``mesh.prob_lo``
     - ``0 0 0``
     - Low corner of the domain (m)
   * - ``mesh.prob_hi``
     - required
     - High corner of the domain (m)
   * - ``mesh.is_periodic``
     - ``0 0 0``
     - Periodic directions
   * - ``mesh.max_grid_size``
     - ``64``
     - Largest box size in x and y
   * - ``mesh.dz0``
     - ``0``
     - First layer thickness of a stretched grid (m)
   * - ``mesh.stretch_ratio``
     - none
     - Fixed ratio between successive layers
   * - ``mesh.zlevels_file``
     - none
     - File with the vertical levels
   * - ``mesh.hierarchy_levels``
     - ``1``
     - Number of grid-hierarchy levels
   * - **Terrain** (`terrain`_)
     -
     -
   * - ``terrain.type``
     - ``flat``
     - Terrain source
   * - ``terrain.file``
     - none
     - Terrain file
   * - ``terrain.smoothing``
     - ``1``
     - Terrain-following mapping (0 BTF, 1 STF, 2 Sullivan)
   * - ``terrain.gamma_m``
     - ``0.5``
     - Smallest layer thickness under STF, relative to flat
   * - ``terrain.nest_blend_width``
     - ``0``
     - Width (m) over which a patch's ground blends to the parent's at ``nest`` faces
   * - ``terrain.hill_height``
     - required for ``agnesi``
     - Hill height (m)
   * - ``terrain.hill_half_width``
     - required for ``agnesi``
     - Hill half-width (m)
   * - ``terrain.hill_shape``
     - ``ridge_x``
     - Ridge or round hill
   * - ``terrain.idw_neighbors``
     - ``6``
     - Points used per height for ``xyz`` terrain
   * - ``terrain.idw_power``
     - ``2``
     - Inverse-distance weight exponent
   * - ``terrain.hierarchy_passes``
     - ``1``
     - Terrain smoothing on coarse hierarchy levels
   * - **Initial conditions** (`init`_)
     -
     -
   * - ``init.first_guess``
     - ``none``
     - Initial profile: ``none``, ``profile``, ``column`` or ``nest``
   * - ``init.projection``
     - ``1``
     - Mass-consistent adjustment on or off
   * - ``init.alpha_h``, ``init.alpha_v``
     - ``1``, ``1``
     - Horizontal and vertical transmissivity
   * - ``init.z_agl_min``
     - ``0``
     - Lowest height at which the profile is evaluated (m)
   * - ``init.k_decay_height``
     - ``0``
     - Height above which the initial TKE decays (m)
   * - ``init.obrien``
     - ``0``
     - O'Brien vertical-velocity correction
   * - ``init.column_tol``
     - ``1e-8``
     - Convergence tolerance of the single-column solve
   * - ``init.column_max_iterations``
     - ``5000``
     - Step limit of the single-column solve
   * - ``init.column_stability``, ``column_mol``, ``column_surface_theta``
     - ``mol``, ``0``, ``0``
     - The column's surface: Obukhov length (0: neutral, θ held), heat flux or surface θ
   * - ``init.column_target_u``, ``column_target_v``, ``column_target_height``
     - ``0``, ``0``, ``0``
     - Target wind of the column: the forcing is adjusted to reach it (0: off)
   * - ``init.column_rib_max``, ``column_rib_min``
     - ``4``, ``-4``
     - Bounds of the bulk Richardson number at the stratified column's surface
   * - ``init.column_time``, ``column_dt``
     - ``0``, ``1``
     - The time-marching column: end time and step (s; 0: the steady column)
   * - ``init.column_wall_k``
     - ``fixed``
     - The time-marching column's first-cell k: the wall value or solved
   * - **Boundary conditions** (`bc`_)
     -
     -
   * - ``bc.xlo`` … ``bc.yhi``
     - ``periodic``
     - Side boundary types
   * - ``bc.zlo``
     - ``surface_layer``
     - Ground boundary type
   * - ``bc.zhi``
     - ``slip_wall``
     - Top boundary type
   * - ``bc.<face>_wall_velocity``
     - ``0 0 0``
     - Velocity of a no-slip wall (m/s)
   * - ``bc.inflow_profile``
     - ``log_law``
     - Inflow profile type
   * - ``bc.inflow_ustar``, ``bc.inflow_z0``
     - 0.4 m/s, 0.1 m
     - Log-law friction velocity and roughness
   * - ``bc.inflow_speed``, ``bc.inflow_zref``, ``bc.inflow_alpha``
     - 10 m/s, 10 m, 0.14
     - Power-law speed, reference height and exponent
   * - ``bc.inflow_direction``
     - ``0``
     - Wind direction (degrees from the x axis)
   * - ``bc.inflow_file``
     - none
     - Tabulated inflow profile
   * - ``bc.inflow_tke_length``
     - ``closure``
     - Length of the inflow profile's equilibrium :math:`k` and :math:`\nu_t`
   * - ``bc.theta_profile_z``, ``bc.theta_profile``
     - ``0``, ``300``
     - Potential-temperature sounding (m, K)
   * - ``bc.zhi_theta_grad``
     - ``0``
     - θ gradient through the top (K/m)
   * - ``bc.rayleigh_zdamp``, ``bc.rayleigh_dampcoef``
     - ``0``, ``0``
     - Rayleigh layer start height (m) and strength (1/s)
   * - ``bc.rayleigh_fields``
     - ``u v w theta tke``
     - Fields the Rayleigh layer damps
   * - ``bc.rayleigh_height``
     - ``absolute``
     - Height used by the Rayleigh layer
   * - ``bc.sponge_width``, ``bc.sponge_strength``
     - ``0``, ``0``
     - Side sponge width (m) and strength (1/s)
   * - ``bc.outflow_transport``
     - ``mirror``
     - Transport fluxes of the ghost cells outside outflow faces
   * - ``bc.surface_z0``
     - 0.1 m
     - Ground roughness length
   * - ``bc.surface_z0_file``
     - none
     - Ground roughness map, ``x y z0`` points (instead of ``bc.surface_z0``)
   * - ``bc.surface_zref``
     - ``0``
     - Wall-model reference height (0: first cell centre)
   * - ``bc.surface_stability``
     - ``neutral``
     - ``neutral`` or ``heat_flux``
   * - ``bc.surface_heat_flux``
     - ``0``
     - Kinematic surface heat flux (K m/s)
   * - ``bc.surface_wall_k``
     - ``dirichlet``
     - Wall treatment of the TKE
   * - ``bc.surface_denominator_floor``
     - ``log_law``
     - Floor of the MOST denominator
   * - ``bc.surface_stable_functions``
     - ``businger_dyer``
     - Stable-side MOST functions: ``businger_dyer`` or ``qnse``
   * - ``bc.surface_qnse_am``, ``_bm``, ``_ah``, ``_bh``
     - ``2.5``, ``0.2``, ``2.5``, ``0.2``
     - QNSE coefficients
   * - ``bc.nest_parent``
     - none
     - Parent run's checkpoint for ``init.first_guess = nest`` and ``nest`` faces
   * - ``bc.nest_relax_width``, ``bc.nest_relax_factor``
     - ``0``, ``1``
     - Relaxation zone inside ``nest`` faces: width (cells) and rate factor
   * - **Forest canopy** (`canopy`_)
     -
     -
   * - ``canopy.type``
     - ``none``
     - ``none``, ``patches`` (ERF's patch file) or ``map`` (point maps)
   * - ``canopy.file``
     - required for ``patches``
     - Patch file, ``type x y height diameter cd lai laimax`` per line
   * - ``canopy.height_file``, ``canopy.lai_file``
     - required for ``map``
     - Canopy height (m) and leaf area index, ``x y value`` points
   * - ``canopy.cd``, ``canopy.cd_file``
     - one required for ``map``
     - Drag coefficient: one value or ``x y value`` points
   * - ``canopy.profile``
     - ``uniform``
     - Leaf area density shape: ``uniform`` or ``lalic_mihailovic``
   * - ``canopy.laimax``
     - ``0.8``
     - Height of the largest leaf area density, as a fraction of the canopy height
   * - ``canopy.drag_diagonal``
     - ``derivative``
     - Implicit drag term on the momentum diagonal: ``derivative``, ``picard`` or ``none``
   * - ``canopy.tke``
     - ``on``
     - Canopy terms of the k–L closure: wake production, short-circuit sink, length cap
   * - ``canopy.beta_p``, ``canopy.beta_d``, ``canopy.length_alpha``
     - ``1``, ``4``, ``0.04``
     - Wake production fraction, short-circuit coefficient, length cap :math:`L \le \alpha/f`
   * - **Buildings** (`buildings`_)
     -
     -
   * - ``buildings.type``
     - ``none``
     - ``none`` or ``immersed_forcing``
   * - ``buildings.file``
     - required when on
     - Height map, heights above the local ground (m)
   * - ``buildings.format``
     - ``erf``
     - ``erf`` (ERF's text layout) or ``xyz`` (``x y height [zstart]`` lines)
   * - ``buildings.force_file``
     - none
     - CSV of the force on the buildings per unit density (m⁴/s²)
   * - ``buildings.small_volfrac``
     - ``0.005``
     - Solid fractions within this of 0 or 1 become 0 or 1
   * - ``buildings.z0``
     - ``0.1``
     - Roughness of the wall law on buildings and, with ``buildings.wall_tke = on``, of the closure's length beside their walls (m)
   * - ``buildings.wall_form``
     - ``fraction_stress``
     - ``fraction_stress`` (form 7 wall law on every wall), ``roof_stress`` (on roofs and undersides alone) or ``drag`` (fixed rows alone)
   * - ``buildings.wall_speed``
     - ``open_faces`` with ``fraction_stress``, ``all_faces`` with ``roof_stress``
     - The wall law's speed in a cell from its faces that are not fixed rows (``open_faces``) or from every face (``all_faces``)
   * - ``buildings.wall_tke``
     - ``on`` with ``fraction_stress``, else ``off``
     - The closure sees the buildings: wall value of k, solid cells at the floor, distance to buildings in the length scale
   * - ``buildings.theta_mode``
     - ``ambient`` with ``fraction_stress``, else ``none``
     - How theta meets a building: no terms, or walls at the reference profile's or a given temperature
   * - ``buildings.theta_surface``
     - none
     - The walls' temperature with ``buildings.theta_mode = fixed`` (K)
   * - ``buildings.limiter_zone``
     - ``3``
     - Zone around each building where ``turb.realizability`` and the upwind blend fade out (building heights)
   * - **Named set-ups** (`abl`_)
     -
     -
   * - ``abl.drive``
     - none
     - What pushes the flow: a name that sets defaults and checks the deck
   * - ``abl.theta``
     - none
     - How θ is treated: a name that sets defaults and checks the deck
   * - **Forcing** (`forcing`_)
     -
     -
   * - ``forcing.gravity``
     - ``9.81``
     - Gravitational acceleration (m/s²)
   * - ``forcing.theta0``
     - ``300``
     - Boussinesq reference temperature (K)
   * - ``forcing.coriolis``
     - ``0``
     - Coriolis force on or off
   * - ``forcing.latitude``
     - ``0``
     - Latitude (degrees)
   * - ``forcing.rotation_period``
     - ``86400``
     - Earth's rotation period (s)
   * - ``forcing.coriolis_3d``
     - ``1``
     - Include the cos(latitude) Coriolis terms
   * - ``forcing.geostrophic``
     - ``0``
     - Geostrophic forcing on or off
   * - ``forcing.ug``, ``forcing.vg``
     - ``0``, ``0``
     - Geostrophic wind (m/s)
   * - ``forcing.pressure_gradient``
     - ``0 0 0``
     - Constant pressure-gradient force (m/s²)
   * - **Turbulence model** (`turb`_)
     -
     -
   * - ``turb.closure``
     - ``kl``
     - ``kl`` (the k-L closure) or ``none`` (laminar)
   * - ``turb.Cmu0``
     - ``0.5562``
     - Model constant :math:`C_{\mu 0}`
   * - ``turb.Cb``
     - ``0.35``
     - Model constant :math:`C_b`
   * - ``turb.Rt_crit``, ``turb.Rt_min``
     - ``-1``, ``-3``
     - Richardson-number smoothing limits
   * - ``turb.max_geom_lscale``
     - 30 m
     - Cap on the geometric length scale
   * - ``turb.tke_floor``
     - ``0``
     - Lower bound on k (0: machine epsilon)
   * - ``turb.sigma_k``
     - ``1``
     - TKE diffusion Schmidt number
   * - ``turb.nu_molecular``, ``turb.alpha_molecular``
     - ``0``, ``0``
     - Molecular viscosity and heat diffusivity (m²/s); ``nu_molecular`` > 0 with ``turb.closure = none``
   * - ``turb.kappa``
     - ``0.41``
     - von Kármán constant
   * - ``turb.nut_relax``
     - ``1``
     - Under-relaxation of :math:`\nu_t`
   * - ``turb.frozen_switch``
     - ``1``
     - Freeze the k-floor decisions in the Newton stage
   * - ``turb.ri_gate``
     - ``0``
     - Stratification gate :math:`Ri_0` (0: off)
   * - ``turb.realizability``
     - ``0``
     - Eddy-viscosity cap :math:`\nu_t \le C_r k/S` (0: off; recommended 0.309 for the lee)
   * - ``turb.face_mean``
     - ``arithmetic``
     - Face average of the diffusivities
   * - ``turb.wall_distance``
     - ``terrain_height``
     - Wall-distance method
   * - **Convection** (`convection`_)
     -
     -
   * - ``convection.limiter``
     - ``van_leer``
     - Slope limiter for momentum and θ
   * - ``convection.smooth_eps_velocity``, ``convection.smooth_eps_theta``
     - ``1e-5``, ``1e-5``
     - Smoothing constants of ``van_albada``
   * - ``convection.hybrid_sigma_max``, ``hybrid_height``, ``hybrid_relax``, ``hybrid_freeze``
     - ``0``, —, ``0.1``, ``1e-2``
     - Blend towards upwind near the ground where the :math:`\nu_t` cap is active (0: off)
   * - **Steady solver** (`steady`_)
     -
     -
   * - ``steady.stage``
     - ``none``
     - ``none`` (build and stop), ``simplec``, ``simplec_newton`` (SIMPLEC, then Newton), or ``unsteady`` (physical
       time steps, `unsteady`_)
   * - ``steady.max_iterations``
     - ``10000``
     - Iteration limit
   * - ``steady.residual_tol``
     - ``5e-5``
     - Scaled residual tolerance
   * - ``steady.residual_drop_orders``
     - ``4``
     - Required residual drop (orders of magnitude)
   * - ``steady.residual_ref_iterations``
     - ``5``
     - Iterations over which the reference residual is taken
   * - ``steady.confidence``, ``steady.confidence_height``
     - ``0``, ``0``
     - Confidence-percentage stop: share of converged cells (0: off), height limit of the others (0: anywhere)
   * - ``steady.monitor_rel_tol``, ``steady.monitor_window``
     - ``1e-3``, ``50``
     - Monitor tolerance and window
   * - ``steady.cfl0``, ``steady.cfl_max``
     - ``1``, ``8`` (``lusgs``) or ``10`` (``zjacobi``)
     - Initial and largest CFL number
   * - ``steady.ser_exponent``, ``steady.ser_max_growth``
     - ``1``, ``2``
     - CFL ramp exponent and largest growth per iteration
   * - ``steady.dtau_floor``
     - ``0.1``
     - Local step floor: a cell steps at most as the level's fastest cell does, 1/``dtau_floor`` times (0: off)
   * - ``steady.switch_window``, ``steady.switch_residual``, ``steady.switch_back_growth``
     - ``20``, ``1e-4``, ``10``
     - Stage-switching rules (``simplec_newton`` only)
   * - ``steady.switch_monitor_tol``, ``steady.newton_fallback_failures``
     - ``1e-2``, ``2``
     - Monitors settled before the switch to Newton; failed Newton steps before the fallback
   * - ``steady.newton_cfl_max``
     - ``128``
     - The Newton stage's CFL cap
   * - ``steady.sweep``
     - ``lusgs`` (CPU build), ``zjacobi`` (GPU build)
     - Line solver: ``zjacobi`` or ``lusgs``
   * - ``steady.zline_sweeps``
     - ``2``
     - Sweeps per iteration
   * - ``steady.line_ghosts``, ``steady.line_overlap``
     - ``lagged``, ``0``
     - LU-SGS coupling across box faces
   * - ``steady.line_order``
     - ``sequential``
     - LU-SGS visiting order
   * - ``steady.line_alternate``
     - ``0``
     - Alternate the LU-SGS sweep direction
   * - ``steady.sequencing``
     - ``0``
     - Grid sequencing, coarse to fine
   * - ``steady.sequencing_drop_orders``
     - ``2``
     - Residual drop for the loose stop of a coarse level
   * - ``steady.sequencing_monitor_tol``
     - ``0``
     - Optional monitor test of the loose stop (0: off)
   * - ``steady.sequencing_hold``
     - ``0``
     - Keep the reference hold on the seeded levels
   * - ``steady.sequencing_cfl_start``
     - ``cfl0``
     - CFL of a seeded level's first iteration
   * - ``steady.machine``, ``steady.nodes``
     - ``local``, ``1``
     - Machine and node count for the allocation accounting
   * - ``steady.newton_*``
     - see `steady`_
     - Newton–Krylov inputs (``simplec_newton`` only)
   * - ``steady.newton_preconditioner``
     - ``simplec``
     - Newton preconditioner: ``identity`` or ``simplec``
   * - **Unsteady runs** (`unsteady`_)
     -
     -
   * - ``unsteady.stop_time``
     - —
     - Physical time the run ends at (s); required with ``steady.stage = unsteady``
   * - ``unsteady.dt``, ``unsteady.cfl``
     - —, ``0``
     - Fixed physical step (s), or a target convective CFL for a global step
   * - ``unsteady.inner_drop_orders``, ``unsteady.max_inner``
     - ``3``, ``50``
     - Inner stop of each physical step: residual drop (orders) and iteration limit
   * - ``unsteady.average_start``
     - ``0``
     - Physical time the running averages and variances start at (s)
   * - ``unsteady.plot_interval``, ``unsteady.checkpoint_interval``
     - ``0``, ``0``
     - A plotfile, a checkpoint every so many physical steps (0: at the end only)
   * - ``unsteady.probes_file``
     - —
     - File of probe points ``x y h`` (h above the ground); one CSV row per step
   * - **SIMPLEC** (`simplec`_)
     -
     -
   * - ``simplec.pressure_reltol``
     - ``1e-3``
     - Tolerance of the pressure-correction solve
   * - ``simplec.alpha_p``
     - ``1``
     - Pressure under-relaxation
   * - ``simplec.rotation_relax``
     - ``0.3``
     - Point-implicit rotation step
   * - ``simplec.w_vnn_cap``
     - ``0``
     - Vertical-velocity cap (0: off)
   * - ``simplec.w_vnn_cap_slope_margin``
     - ``4``
     - The cap's safety margin on sloping ground (0: off)
   * - ``simplec.theta``
     - ``solve``
     - ``frozen``: θ prescribed from the inflow profile, buoyancy in the closure only
   * - **Poisson solver** (`solver`_)
     -
     -
   * - ``solver.poisson_backend``
     - ``gmres_mlmg``
     - Poisson solver
   * - ``solver.poisson_reltol``, ``solver.poisson_abstol``
     - ``1e-10`` (``1e-5`` single), ``0``
     - Relative and absolute tolerance
   * - ``solver.poisson_max_iter``
     - ``200``
     - Iteration limit
   * - ``solver.poisson_restart``
     - ``50``
     - GMRES restart length
   * - ``solver.poisson_precond``
     - ``fft``
     - Preconditioner of ``gmres_fft``
   * - ``solver.poisson_verbose``
     - ``0``
     - Output level
   * - **Input/output** (`io`_)
     -
     -
   * - ``io.plot_file``, ``io.chk_file``
     - ``plt``, ``chk``
     - Plotfile and checkpoint prefixes
   * - ``io.write_plotfile``, ``io.write_checkpoint``
     - ``1``, ``1``
     - Write the plotfile and checkpoint
   * - ``io.turb_plot_file``
     - none
     - Prefix of the turbulence plotfile
   * - ``io.iteration_log``
     - ``table``
     - Layout of the per-iteration log
   * - ``io.column_guess_file``
     - ``column_guess.dat``
     - File of the converged column's table
   * - ``io.restart``
     - none
     - Checkpoint to restart from

``mesh``
--------

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Input
     - Default
     - Range and meaning
   * - ``mesh.n_cell``
     - required
     - Three integers ≥ 1: cells in x, y and z.
   * - ``mesh.prob_lo``
     - ``0 0 0``
     - Three finite values: the low corner of the domain (m). The z value must be 0.
   * - ``mesh.prob_hi``
     - required
     - Three finite values, each greater than ``prob_lo``: the high corner of the domain (m).
   * - ``mesh.is_periodic``
     - ``0 0 0``
     - Three values, each 0 or 1. z must be 0.
   * - ``mesh.max_grid_size``
     - ``64``
     - One or two integers ≥ 1: the largest box size in x and y. Boxes are never split in z.
   * - ``mesh.dz0``
     - ``0``
     - Finite, ≥ 0 (m). First layer thickness of a stretched grid; 0 gives uniform layers. Without
       ``mesh.stretch_ratio`` it must not exceed ``prob_hi`` z / ``n_cell`` z (:doc:`terrain_coordinates`).
   * - ``mesh.stretch_ratio``
     - none
     - Finite, ≥ 1; needs ``mesh.dz0`` > 0. A fixed ratio between successive layers; the top is then where the
       layers end, and ``prob_hi`` z only sets the computational spacing.
   * - ``mesh.zlevels_file``
     - none
     - A file of ``n_cell`` z + 1 increasing heights starting at 0. Not together with ``mesh.dz0``.
   * - ``mesh.hierarchy_levels``
     - ``1``
     - Integer from 1 to 10: the levels :math:`L` of the grid hierarchy (:doc:`grid_hierarchy`). With
       :math:`L > 1`, ``n_cell`` and ``max_grid_size`` in x and y must be divisible by :math:`2^{L-1}`, and the
       coarsest level must be larger than 1 × 1 cells. The box layout of level 0 can then differ from :math:`L = 1`,
       which changes results in the last bits. The solver uses the coarse levels only with ``steady.sequencing =
       1``. May change across a restart.

``terrain``
-----------

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``terrain.type``
     - ``flat``
     - ``flat``, ``agnesi`` (analytic hill), ``xyz`` (point cloud) or ``erf_grid`` (ERF's terrain file)
       (:doc:`terrain_coordinates`).
   * - ``terrain.file``
     - required for ``xyz``, ``erf_grid``
     - The point file (``x y z`` per line) or ERF's structured terrain file.
   * - ``terrain.smoothing``
     - ``1``
     - 0 (BTF), 1 (STF) or 2 (Sullivan).
   * - ``terrain.gamma_m``
     - ``0.5``
     - > 0 and < 1; STF only. Smallest layer thickness as a fraction of the flat layer.
   * - ``terrain.hill_height``
     - required for ``agnesi``
     - Finite, > 0 (m).
   * - ``terrain.hill_half_width``
     - required for ``agnesi``
     - Finite, > 0 (m). Distance from the crest at which the hill is half as high.
   * - ``terrain.hill_shape``
     - ``ridge_x``
     - ``ridge_x`` (a ridge along y, varying in x), ``ridge_y`` or ``round``; ``agnesi`` only.
   * - ``terrain.idw_neighbors``
     - ``6``
     - ≥ 1; ``xyz``, ``bc.surface_z0_file`` or ``canopy.type = map`` only. Points used for each height (and each
       column's :math:`z_0` of the roughness map, and the canopy's maps).
   * - ``terrain.idw_power``
     - ``2``
     - Finite, > 0; ``xyz``, ``bc.surface_z0_file`` or ``canopy.type = map`` only. Weight exponent :math:`p` in
       :math:`d^{-p}`.
   * - ``terrain.hierarchy_passes``
     - ``1``
     - Integer from 0 to 64; only with ``mesh.hierarchy_levels`` > 1. The ground of level :math:`l` is smoothed by
       this number times :math:`(2^l - 1)` passes (:doc:`grid_hierarchy`).
   * - ``terrain.nest_blend_width``
     - ``0``
     - Finite, ≥ 0 (m), below half the patch's width across its ``nest`` faces; 0: off. Refused without a ``nest``
       face. The patch's ground becomes the parent's on the ``nest`` faces and its own from this distance inwards,
       with a smoothstep between (:doc:`boundary_conditions`).

``init``
--------

See :doc:`initial_conditions`. The keys other than ``init.first_guess`` are refused with ``init.first_guess = none``.
A restart cannot be combined with ``init.*`` keys.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``init.first_guess``
     - ``none``
     - ``none`` (fluid at rest), ``profile`` (the inflow profile by height above the ground) or ``column`` (the
       steady single-column profile; needs a horizontal forcing and ``steady.stage = simplec`` or ``simplec_newton``, and
       with a canopy flat terrain and a canopy that is the same in every column; see :doc:`initial_conditions` for
       the full list) or ``nest`` (a parent run's checkpoint, ``bc.nest_parent``, interpolated at 2 2 1; not with
       ``steady.sequencing = 1``; see :doc:`initial_conditions`).
   * - ``init.projection``
     - ``1``
     - 0 or 1: the mass-consistent adjustment. With 0, ``init.alpha_h`` and ``init.alpha_v`` are refused.
   * - ``init.alpha_h``
     - ``1``
     - Finite, > 0. Horizontal transmissivity :math:`\alpha_h`.
   * - ``init.alpha_v``
     - ``1``
     - Finite, > 0. Vertical transmissivity :math:`\alpha_v`. Only :math:`\alpha_v/\alpha_h` matters. Ratios below
       0.1 need more than the default ``solver.poisson_max_iter``. With :math:`\alpha_v \ne \alpha_h` the
       adjustment uses the ``gmres_fft`` back end.
   * - ``init.z_agl_min``
     - ``0``
     - Finite, ≥ 0 (m). Lowest height above the ground at which the profile is evaluated. Not with ``nest``.
   * - ``init.k_decay_height``
     - ``0``
     - Finite, ≥ 0 (m). Height :math:`z_i` above which the initial :math:`k` decays; 0: no decay. Not with
       ``column`` or ``nest``.
   * - ``init.obrien``
     - ``0``
     - 0 or 1. The O'Brien correction of the vertical velocity before the adjustment. Not with ``nest``.
   * - ``init.column_tol``
     - ``1e-8``
     - Finite, > 0; ``column`` only. Convergence tolerance of the single-column solve: the largest relative change
       of :math:`u, v, k` over a step. Raised to at least 100 machine epsilon (1.2e-5 single, 2.2e-14 double).
   * - ``init.column_max_iterations``
     - ``5000``
     - ≥ 1; ``column`` only. The run stops if a column has not converged after this many steps. The limit applies
       to each column solved: the surface search and the target wind solve several (the log gives the total).
   * - ``init.column_stability``
     - ``mol``
     - ``mol``, ``heat_flux`` or ``surface_temperature``; ``column`` only. How the column's surface sets θ
       (:ref:`stratified column <column-stability>`). ``heat_flux`` takes ``bc.surface_heat_flux`` and needs
       ``bc.surface_stability = heat_flux``. Any non-neutral surface needs ``forcing.gravity`` > 0 and, in the steady
       column, a Rayleigh layer on θ (``bc.rayleigh_zdamp`` > 0, ``theta`` in ``bc.rayleigh_fields``).
   * - ``init.column_mol``
     - ``0``
     - 0 (neutral: θ held at the sounding, the column before this input) or finite with \|L\| ≥ 1 (m); ``mol`` only.
       The Obukhov length L (> 0 stable, < 0 unstable). With ``bc.surface_stability = heat_flux`` the 3-D run takes
       the column's flux (``bc.surface_heat_flux`` is then refused).
   * - ``init.column_surface_theta``
     - ``0``
     - Finite, > 0 (K), required with ``surface_temperature`` and only there: the surface potential temperature.
   * - ``init.column_target_u``, ``init.column_target_v``
     - ``0``, ``0``
     - Finite (m/s); ``column`` only. A target wind at ``init.column_target_height``: the geostrophic wind and the
       column's Rayleigh reference wind are turned and scaled by the ratio of the target to the column's wind there,
       the pressure gradient by that ratio times its magnitude, until they agree to 1e-4; the run then uses that
       forcing. Both 0: off (giving them as 0 is refused). Refused with what would hold the 3-D wind to the deck's
       inflow profile: a Rayleigh layer on u and v, a sponge or an inflow face; and with a stratified column unless
       ``bc.surface_stability = heat_flux``.
   * - ``init.column_target_height``
     - ``0``
     - Finite, > 0 (m) with a target wind, and only then: its height above the ground, at most the top cell centre.
       Below the first cell centre the wind is MOST's profile; above, linear in :math:`\ln(z + z_0)` between centres.
   * - ``init.column_rib_max``
     - ``4``
     - Finite, > 0; a stratified column only. The upper bound of the bulk Richardson number of the column's surface
       solve (ERF's value). With Businger–Dyer's stable functions the bound used is min(``init.column_rib_max``, 0.19),
       below their 0.2 limit (:ref:`stratified column <column-stability>`).
   * - ``init.column_rib_min``
     - ``-4``
     - Finite, < 0; a stratified column only. The lower bound (ERF's value). Both bounds are refused with
       ``init.column_time`` > 0 and ``heat_flux`` (the flux is given; no bulk Richardson number is solved).
   * - ``init.column_time``
     - ``0``
     - Finite, ≥ 0 (s); ``column`` only. 0: the steady column. > 0: the time-marching column
       (:ref:`time-marching column <column-marching>`), backward Euler steps of ``init.column_dt`` in physical time
       to this time; the table is the column at that time. Refused with ``init.column_tol``,
       ``init.column_max_iterations``, a target wind and a non-zero ``init.column_mol``.
   * - ``init.column_dt``
     - ``1``
     - Finite, > 0 (s); with ``init.column_time`` > 0 only. The step; the last step ends at ``init.column_time``.
       ``init.column_time / init.column_dt`` must be at most 1e7 steps.
   * - ``init.column_wall_k``
     - ``fixed``
     - ``fixed`` or ``solved``; with ``init.column_time`` > 0 only. ``fixed``: the first cell's k is the surface
       layer's wall value, as in the steady column and the 3-D run. ``solved``: the first cell's k is solved with
       the k equation of the other cells and no k flux through the ground (:ref:`time-marching column
       <column-marching>`).

``bc``
------

See :doc:`boundary_conditions`. Read when ``init.first_guess`` is not ``none`` or with the SIMPLEC stage; refused
otherwise. A key that has no effect for the chosen boundary types is refused. The initial field uses the boundary
types, the inflow profile and the sounding; the other keys act in the SIMPLEC stage.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``bc.xlo``, ``bc.xhi``, ``bc.ylo``, ``bc.yhi``
     - ``periodic``
     - ``periodic`` exactly when ``mesh.is_periodic`` is 1 in that direction; otherwise ``inflow``, ``outflow``,
       ``slip_wall``, ``no_slip_wall`` or ``nest`` (the parent run's field on the face of a nested patch,
       :doc:`boundary_conditions`). ``nest`` needs ``bc.nest_parent`` and ``init.first_guess = nest`` or
       ``io.restart``; it is refused with ``steady.sequencing = 1``, with a ``slip_wall`` or ``no_slip_wall`` lateral
       face, with buildings, and with an ``inflow`` face unless there is an ``outflow`` face too.
   * - ``bc.zlo``
     - ``surface_layer``
     - ``surface_layer``, ``slip_wall`` or ``no_slip_wall`` (``nest`` is for the lateral faces only).
   * - ``bc.zhi``
     - ``slip_wall``
     - ``slip_wall``, ``no_slip_wall`` or ``geostrophic`` (needs Coriolis with geostrophic forcing).
   * - ``bc.<face>_wall_velocity``
     - ``0 0 0``
     - Three finite values (m/s): the velocity of a ``no_slip_wall`` face (``<face>`` is ``xlo`` … ``zhi``). The
       normal component must be 0.
   * - ``bc.inflow_profile``
     - ``log_law``
     - ``log_law``, ``power_law`` or ``file``. The inflow keys need an inflow face, a damping layer,
       ``init.first_guess`` or a restart.
   * - ``bc.inflow_ustar``
     - 0.4 m/s
     - Finite, > 0; ``log_law`` only.
   * - ``bc.inflow_z0``
     - 0.1 m
     - Finite, > 0. Roughness of the log law, and of the equilibrium :math:`k` of every profile type. Normally equal
       to ``bc.surface_z0``.
   * - ``bc.inflow_speed``, ``bc.inflow_zref``, ``bc.inflow_alpha``
     - 10 m/s, 10 m, 0.14
     - ``power_law`` only. Speed finite and ≥ 0; height and exponent finite and > 0.
   * - ``bc.inflow_direction``
     - ``0``
     - Finite (degrees). The direction the wind blows towards, counter-clockwise from the x axis. ``log_law`` and
       ``power_law`` only; refused with ``file``, whose profile gives u and v itself.
   * - ``bc.inflow_file``
     - none
     - Required with ``file``: rows of ``z u v theta [k]`` (4 or 5 columns, at least 2 rows), z increasing from
       0, θ > 0, k ≥ 0. ``#`` starts a comment. Every token must be a complete number with nothing attached (``300.0K``,
       ``1.0D+00`` or ``1,2`` stops the run with the file and line named); a UTF-8 byte-order mark is skipped. Read on the I/O rank and broadcast.
   * - ``bc.inflow_tke_length``
     - ``closure``
     - ``closure`` or ``unbounded``: the length :math:`l` of the equilibrium :math:`k` and :math:`\nu_t`
       (:doc:`boundary_conditions`). ``closure``: the turbulence model's :math:`\kappa(z + z_0)` capped by
       ``turb.max_geom_lscale``. ``unbounded``: :math:`\kappa(z + z_0)`, the form before this input (for the log law,
       :math:`k = u_*^2/C_{\mu 0}^2` at every height). Refused with a ``k`` column in ``bc.inflow_file`` (no effect).
       A checkpoint written before this input continues unchanged only with ``unbounded``: the inflow
       :math:`k` is not stored in the checkpoint.
   * - ``bc.theta_profile_z``, ``bc.theta_profile``
     - ``0``, ``300``
     - Lists of equal length: finite, strictly increasing heights (m) and θ > 0 (K). The reference sounding.
   * - ``bc.zhi_theta_grad``
     - ``0``
     - Finite (K/m): θ gradient through a ``slip_wall`` or ``geostrophic`` top.
   * - ``bc.rayleigh_zdamp``
     - ``0`` (off)
     - In [0, top level): start height of the Rayleigh layer (m); > 0 needs ``bc.rayleigh_dampcoef`` > 0. The top
       level is ``mesh.prob_hi`` z, or where the layers end with ``mesh.stretch_ratio`` or ``mesh.zlevels_file``.
   * - ``bc.rayleigh_dampcoef``
     - ``0``
     - Finite, ≥ 0 (1/s); needs ``bc.rayleigh_zdamp`` > 0.
   * - ``bc.rayleigh_fields``
     - ``u v w theta tke``
     - One or more of ``u v w theta tke``; needs ``bc.rayleigh_zdamp`` > 0. With a ``nest`` face and no ``inflow``
       face only ``w`` and ``theta`` (u, v and k relax towards the ``bc.inflow_*`` profile, not the parent).
   * - ``bc.rayleigh_height``
     - ``absolute``
     - ``absolute`` or ``above_ground``: the height used for the layer's coefficient and for θ's reference
       (:math:`u, v, k` always use the height above the ground).
   * - ``bc.sponge_width``
     - ``0`` (off)
     - Finite, ≥ 0 (m); > 0 needs an inflow or outflow face and ``bc.sponge_strength`` > 0. Refused with a ``nest``
       face and no ``inflow`` face.
   * - ``bc.sponge_strength``
     - ``0``
     - Finite, ≥ 0 (1/s); needs ``bc.sponge_width`` > 0.
   * - ``bc.outflow_transport``
     - ``mirror``
     - ``mirror`` or ``from_velocity``; needs an ``outflow`` face. ``mirror``: each ghost cell outside an outflow face
       has the transport fluxes of the cell inside it, mirrored through the face (the normal flux continued linearly
       through its face value), so it has that cell's divergence. ``from_velocity``: the fluxes are computed from the
       ghost cells' copied velocities and their own geometry, the only choice up to v1.0.0. The velocity ghosts are the
       same either way (:doc:`boundary_conditions`).
   * - ``bc.surface_z0``
     - 0.1 m
     - Finite, > 0. Roughness length of the wall model. The ``bc.surface_*`` keys need ``bc.zlo = surface_layer``.
   * - ``bc.surface_z0_file``
     - none (``bc.surface_z0``)
     - A file of ``x y z0`` points (``terrain.file``'s format), interpolated to the ground columns in
       :math:`\ln z_0` with ``terrain.idw_neighbors`` and ``terrain.idw_power``
       (:ref:`roughness map <roughness-map>`). Not with ``bc.surface_z0`` (both given aborts) or
       ``init.first_guess = column``. Every value finite and > 0, at least ``terrain.idw_neighbors`` points, and
       the file readable, or the run aborts naming the key. ``bc.inflow_z0`` stays the inflow's roughness.
   * - ``bc.surface_zref``
     - ``0``
     - Finite. ≤ 0: each column's first cell centre. Otherwise a height between the first and last cell-centre
       heights of every column (m).
   * - ``bc.surface_stability``
     - ``neutral``
     - ``neutral`` or ``heat_flux``.
   * - ``bc.surface_heat_flux``
     - ``0``
     - Finite (K m/s); ``heat_flux`` only.
   * - ``bc.surface_wall_k``
     - ``dirichlet``
     - ``dirichlet`` or ``source``.
   * - ``bc.surface_denominator_floor``
     - ``log_law``
     - ``log_law`` (:math:`m = \min(\ln((d + z_0)/z_0), 1)`) or ``one`` (:math:`m = 1`): the floor of the MOST
       denominator. They differ only where :math:`\ln((d + z_0)/z_0) < 1`.
   * - ``bc.surface_stable_functions``
     - ``businger_dyer``
     - ``businger_dyer`` or ``qnse``: the stable side (:math:`z/L > 0`) of :math:`\phi_m, \phi_h, \psi_m, \psi_h`
       (:ref:`stable functions <stable-functions>`), for the wall model and the column. ``qnse`` needs
       ``bc.surface_stability = heat_flux`` or a stratified column (``init.column_stability``), otherwise the run
       stops.
   * - ``bc.surface_qnse_am``, ``bc.surface_qnse_bm``
     - ``2.5``, ``0.2``
     - :math:`a_m, b_m` of :math:`\phi_m = (1 + a_m\zeta)/(1 + b_m\zeta)`: :math:`b_m > 0`, :math:`a_m \ge b_m`,
       finite. ``qnse`` only (otherwise the run stops).
   * - ``bc.surface_qnse_ah``, ``bc.surface_qnse_bh``
     - ``2.5``, ``0.2``
     - :math:`a_h, b_h` of :math:`\phi_h`, as :math:`a_m, b_m`. ``qnse`` only.
   * - ``bc.nest_parent``
     - none
     - A checkpoint directory of the parent run, written by this version (checkpoint format 2, which holds the
       z-levels); required with ``init.first_guess = nest`` or a ``nest`` face and refused without both (no effect). The patch must refine
       the parent by 2 2 1, lie on the parent's cell faces (to :math:`10^{-4}` of a parent cell), keep 2 parent cells
       from every non-periodic face of the parent (it may reach a periodic face), span the full width of a direction in
       which it is periodic (``mesh.prob_lo`` and ``mesh.prob_hi`` equal to the parent's there), and have the parent's
       z-levels bit for bit (``mesh.n_cell`` z, ``mesh.prob_lo`` z, ``mesh.prob_hi`` z, ``mesh.dz0``,
       ``mesh.stretch_ratio``, ``mesh.zlevels_file``). Any other placement, a directory without a checkpoint ``Header`` or a
       format-1 checkpoint stops the run, naming this key (:doc:`initial_conditions`).
   * - ``bc.nest_relax_width``
     - ``0``
     - Integer ≥ 0 (cells), below half the patch's cells across its ``nest`` faces; 0: no zone. Refused without a
       ``nest`` face. The number of cells inside each ``nest`` face over which u, v, w, θ and k are relaxed towards
       the parent's field where the parent's flow enters. A robustness setting: it holds the patch to the parent's
       coarser solution over the zone (:doc:`boundary_conditions`).
   * - ``bc.nest_relax_factor``
     - ``1``
     - Finite, > 0. Multiplies the zone's rate; with 1 the rate at the face is the parent's inward velocity over the
       zone's width. About :math:`e^{-f/2}` of a difference from the parent passes the zone (61 % with 1, 0.7 % with
       10). Refused without ``bc.nest_relax_width`` > 0.

The wall model and the inflow profile's equilibrium :math:`k` and :math:`\nu_t` use the turbulence model's
``turb.Cmu0`` and ``turb.kappa``, and the inflow profile its ``turb.max_geom_lscale``. Without the SIMPLEC stage
(``turb.*`` is then refused) they use the defaults, 0.5562, 0.41 and 30 m.

``canopy``
----------

See :doc:`canopy`. The canopy's drag field is built in every run (with or without ``steady.stage``) and written to the
plotfile as ``canopy_drag``; its drag acts on the momentum equations and, with ``canopy.tke``, on the k–L closure. A
key that has no effect with the chosen ``canopy.type``, ``canopy.profile`` or ``canopy.tke`` is refused. A restart
refuses a canopy geometry other than that of the run it continues (``canopy.type`` and the keys that build the drag
field); ``canopy.drag_diagonal``, the closure's canopy terms (``canopy.tke``, ``canopy.beta_p``,
``canopy.beta_d``, ``canopy.length_alpha``) and the leaf heat (``canopy.heat``, ``canopy.leaf_theta``; ``canopy.heat_frozen`` needs frozen θ, which does
not restart) may change across a restart, as the ``turb.*`` inputs may (the closure's
inputs then change the answer, not only the path; a checkpoint of fully frozen closure decisions keeps
``canopy.length_alpha``, :doc:`canopy`).

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``canopy.type``
     - ``none``
     - ``none``, ``patches`` (round patches from ERF's forest file) or ``map`` (point maps of height, LAI and
       optionally :math:`C_d`, interpolated with ``terrain.idw_neighbors`` and ``terrain.idw_power``).
   * - ``canopy.file``
     - required for ``patches``
     - A readable file of patches, one per line, ``type x y height diameter cd lai laimax`` (ERF's
       ``erf.forest_file``); ``#`` starts a comment. type 1 or 2 (uniform or Lalic–Mihailovic), height and diameter
       finite and > 0 (m), cd and lai finite and ≥ 0, laimax in [0, 1) on a type 2 line and finite on a type 1 line.
       Any other line stops the run, naming the file and the line. A patch wider than a periodic direction of the
       domain is refused.
   * - ``canopy.height_file``
     - required for ``map``
     - A readable file of ``x y h`` points (``terrain.file``'s xyz format): the canopy height (m), every value finite
       and ≥ 0, at least ``terrain.idw_neighbors`` points.
   * - ``canopy.lai_file``
     - required for ``map``
     - A readable file of ``x y LAI`` points: the leaf area index, every value finite and ≥ 0, at least
       ``terrain.idw_neighbors`` points.
   * - ``canopy.cd``
     - none
     - Finite, > 0; ``map`` only. The drag coefficient of every column. ``map`` needs this or ``canopy.cd_file``, not
       both (no default, as in ERF). 0 is refused: without a canopy, set ``canopy.type = none``. The maps, by
       contrast, may hold 0 where there is no forest.
   * - ``canopy.cd_file``
     - none
     - ``map`` only. A readable file of ``x y cd`` points, every value finite and ≥ 0, at least
       ``terrain.idw_neighbors`` points.
   * - ``canopy.profile``
     - ``uniform``
     - ``uniform`` or ``lalic_mihailovic``; ``map`` only. The shape of the leaf area density over the canopy height.
   * - ``canopy.laimax``
     - ``0.8``
     - Finite, in [0, 1); ``map`` with ``lalic_mihailovic`` only. The height of the largest leaf area density as a
       fraction of the canopy height (ERF's ``erf.forest_laimax``).
   * - ``canopy.drag_diagonal``
     - ``derivative``
     - ``derivative``, ``picard`` or ``none``; ``patches`` or ``map``. The drag's term on the diagonal of the
       momentum rows: the exact derivative :math:`f_c(|U| + u_c^2/|U|)`, :math:`f_c|U|` with :math:`|U|` frozen,
       or nothing (the drag explicit, as in ERF). The converged answer does not depend on it; the iterations do. May
       change across a restart.
   * - ``canopy.tke``
     - ``on``
     - ``on`` or ``off``; ``patches`` or ``map``. The canopy terms of the k–L closure (:doc:`canopy`): the wake
       production :math:`\beta_p f |U|^3`, the short-circuit sink :math:`\beta_d f |U| k` and the length cap
       :math:`L \le \alpha/f`. ``off`` keeps the drag only.
   * - ``canopy.beta_p``
     - ``1``
     - Finite, in [0, 1]; ``canopy.tke = on`` only. The fraction of the drag work :math:`f|U|^3` that becomes
       :math:`k` (ERF's ``erf.forest_canopy_beta_p``).
   * - ``canopy.beta_d``
     - ``4``
     - Finite, ≥ 0; ``canopy.tke = on`` only. The short-circuit coefficient of the sink :math:`\beta_d f |U| k`
       (ERF's ``erf.forest_canopy_beta_d``).
   * - ``canopy.length_alpha``
     - ``0.04``
     - Finite, in [0, 1]; ``canopy.tke = on`` only. The cap :math:`L \le \alpha/f` on the turbulent length; 0 turns
       the cap off and keeps the sources (ERF's ``erf.forest_canopy_length_alpha``). ``canopy.beta_p``,
       ``canopy.beta_d`` and ``canopy.length_alpha`` all 0 with ``canopy.tke = on`` is refused (it changes
       nothing).
   * - ``canopy.heat``
     - ``off``
     - ``off`` or ``fixed``; ``patches`` or ``map``. The leaves' sensible heat exchange with the air at a fixed leaf
       temperature (:doc:`canopy`): :math:`\theta` gains :math:`\lambda(\theta_l - \theta)`,
       :math:`\lambda = 2a/r_b`. ``off`` leaves :math:`\theta` as without it, bit for bit.
   * - ``canopy.leaf_theta``
     - none
     - In [150, 400] (K: a value in °C is refused); required with ``canopy.heat = fixed`` and refused otherwise. The
       leaves' potential temperature :math:`\theta_l` (ERF's ``erf.forest_leaf_theta_fixed``); "fixed" names this
       fixed leaf temperature.
   * - ``canopy.heat_frozen``
     - ``skip``
     - ``skip`` or ``blend``; ``canopy.heat = fixed`` with ``simplec.theta = frozen`` only (refused otherwise).
       ``skip``: the leaf heat is left out, with a warning. ``blend``: in the canopy :math:`\theta` is set each
       iteration to :math:`\theta_f + w(\theta_l - \theta_f)`, :math:`w = \lambda k/(\varepsilon + \lambda k)`, from
       the frozen :math:`\theta_f` (:doc:`canopy`).

``buildings``
-------------

See :doc:`buildings`. The geometry is built in a run without ``steady.stage`` and written to the plotfile as
``building_beta`` and ``building_wall``. With a stage it acts on momentum: SIMPLEC, z-line sweeps, one level, and
``solver.poisson_backend`` other than ``mlmg``; others are refused. It acts on the closure too (``buildings.wall_tke``) and on theta (``buildings.theta_mode``). A key that has no effect with ``buildings.type = none`` is refused. A restart
refuses buildings other than those of the run it continues.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Meaning and range
   * - ``buildings.type``
     - ``none``
     - ``none`` or ``immersed_forcing``. Refused with a ``nest`` face.
   * - ``buildings.file``
     - required when on
     - The height map, readable. Heights above the local ground (m), finite and ≥ 0. At least one building in the domain, none reaching the top layer, and in a periodic direction the map no wider than the domain.
   * - ``buildings.format``
     - ``erf``
     - ``erf``: ``nx``, ``ny``, the x values, the y values, the heights with x the slow index. ``xyz``: ``x y height`` or ``x y height zstart`` lines holding each point of a rectangular grid once, in any order; solid from ``zstart`` (≥ 0, and below ``height`` where ``height`` > 0; default 0) to ``height``.
   * - ``buildings.force_file``
     - none
     - A file name, not empty: the CSV ``t,Fx,Fy,Fz`` of the force on the buildings per unit density (m⁴/s²), from the momentum residual over their blocked faces (:doc:`buildings`, Force on the buildings). One row after a steady solve; with ``steady.stage = unsteady`` one at the start and one per step, and a restart keeps the rows up to the checkpoint's time. Needs ``steady.stage``.
   * - ``buildings.small_volfrac``
     - ``0.005``
     - In [0, 0.5). Below it a fraction is 0, within it of 1 the fraction is 1. At 0, round-off alone is clipped.
   * - ``buildings.z0``
     - ``0.1``
     - Finite and > 0 (m). Refused with ``buildings.wall_form = drag``. The law reads distances below ``10 z0`` as ``10 z0`` (:doc:`buildings`, Limits). With ``buildings.wall_tke = on`` it is also the roughness in the closure's length beside a wall (:doc:`buildings`, Length scale).
   * - ``buildings.wall_form``
     - ``fraction_stress``
     - ``fraction_stress``: the form 7 wall law on roofs, undersides and sides, with the fixed rows. ``roof_stress``: the law on roofs and undersides alone, the sides with their fixed rows. ``drag``: the fixed rows alone.
   * - ``buildings.wall_speed``
     - ``open_faces`` with ``buildings.wall_form = fraction_stress``, ``all_faces`` with ``roof_stress``
     - ``open_faces`` or ``all_faces``. ``open_faces``: each velocity component of a wall cell's speed is the mean of its two faces that are not fixed rows, so an edge or corner cell takes the open face's value. ``all_faces``: the mean of both faces, the form before. Either may be given with ``roof_stress``. Refused with ``buildings.wall_form = drag`` (no wall law). A restart from a checkpoint written before this key needs ``all_faces``.
   * - ``buildings.wall_tke``
     - ``on`` with ``buildings.wall_form = fraction_stress``, else ``off``
     - ``on`` or ``off``. ``on``: solid cells hold k at ``turb.tke_floor``, wall cells relax k to ``u*^2 / Cmu0^2`` of the wall law, k does not diffuse through blocked faces, and the closure's length takes the nearer of the ground and the nearest building wall, each with its own roughness. ``off``: the closure does not see the buildings. ``on`` is refused unless ``buildings.wall_form = fraction_stress``.
   * - ``buildings.theta_mode``
     - ``ambient`` with ``buildings.wall_form = fraction_stress``, else ``none``
     - ``none``, ``ambient`` or ``fixed``. ``ambient`` and ``fixed``: solid cells hold theta at the walls' temperature, each wall exchanges heat with its cell by the form 7 law for a given surface temperature (neutral, with the wall law's u*), and theta does not diffuse through blocked faces. ``ambient``: the walls are at the reference profile (``bc.theta_profile``) at their height. ``fixed``: at ``buildings.theta_surface``. ``none``: theta does not see the buildings. ``ambient`` and ``fixed`` are refused unless ``buildings.wall_form = fraction_stress``, and with ``simplec.theta = frozen`` when given.
   * - ``buildings.theta_surface``
     - none
     - Finite and > 0 (K). Required with ``buildings.theta_mode = fixed``; refused otherwise.
   * - ``buildings.limiter_zone``
     - ``3``
     - 0 (off) or finite, ≥ 0.5 and ≤ 20 (building heights). With ``turb.realizability`` > 0, the eddy-viscosity cap is lifted within half the zone of a building and fades back in to its edge, measured by the distance to the nearest building wall in heights of that wall's building; the upwind blend, which follows the cap, fades with it. When given, refused without ``buildings.wall_tke = on`` or without ``turb.realizability`` > 0; the default acts only with both.

``abl``
-------

See :ref:`named-setups`. Read with the SIMPLEC stage only; refused otherwise. Without these keys nothing is checked
and no default changes. A name that the deck's other inputs contradict stops the run.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``abl.drive``
     - none
     - ``inflow_only``, ``couette``, ``pressure_driven``, ``geostrophic_column`` or ``geostrophic_inflow``.
   * - ``abl.theta``
     - none
     - ``neutral``, ``capped``, ``frozen`` or ``surface_flux``. ``abl.drive = couette`` needs ``neutral``.

``forcing``
-----------

See :doc:`governing_equations`. Read with the SIMPLEC stage only.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``forcing.gravity``
     - ``9.81``
     - Finite, ≥ 0 (m/s²). 0 switches buoyancy off.
   * - ``forcing.theta0``
     - ``300``
     - Finite, > 0 (K). Boussinesq reference temperature. Refused with ``forcing.gravity = 0``.
   * - ``forcing.coriolis``
     - ``0``
     - 0 or 1: the Coriolis force.
   * - ``forcing.latitude``
     - ``0``
     - In [−90, 90] (degrees). With ``coriolis = 1`` only.
   * - ``forcing.rotation_period``
     - ``86400``
     - Finite, > 0 (s); :math:`f = 4\pi\sin\varphi/T`. With ``coriolis = 1`` only. The default is the solar day (ERF's
       value); the sidereal day, 86164 s, gives an f 0.27 % larger.
   * - ``forcing.coriolis_3d``
     - ``1``
     - 0 or 1: the :math:`\cos\varphi` terms on :math:`u` and :math:`w`. With ``coriolis = 1`` only.
   * - ``forcing.geostrophic``
     - ``0``
     - 0 or 1: geostrophic forcing. Needs ``coriolis = 1`` and a non-zero latitude.
   * - ``forcing.ug``, ``forcing.vg``
     - ``0``
     - Finite (m/s): the geostrophic wind. With ``geostrophic = 1`` only.
   * - ``forcing.pressure_gradient``
     - ``0 0 0``
     - Three finite values (m/s²): a constant pressure-gradient force on :math:`u, v, w`.

``turb``
--------

See :doc:`turbulence_model`. Read with the SIMPLEC stage only.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``turb.closure``
     - ``kl``
     - ``kl`` or ``none``. ``kl``: the k-L closure (:doc:`turbulence_model`). ``none``: laminar flow; the viscosity is ``turb.nu_molecular`` (then required > 0), the heat diffusivity ``turb.alpha_molecular``, :math:`\nu_t` and every closure source are 0, and :math:`k` is a passive scalar diffusing with ``turb.nu_molecular``. With ``none`` only ``turb.nu_molecular``, ``turb.alpha_molecular`` and ``turb.tke_floor`` may be given (``turb.face_mean`` has no effect on constant diffusivities), and a canopy, ``bc.zlo = surface_layer``, ``init.first_guess = column``, ``steady.stage = simplec_newton``, ``buildings.wall_form`` other than ``drag`` and ``buildings.wall_tke = on`` are refused. May change across a restart, as the other ``turb.*`` inputs: ``kl`` to ``none`` leaves the checkpoint's :math:`\nu_t` history unused, ``none`` to ``kl`` starts the closure without one, and an unsteady restart's running averages then hold both models.
   * - ``turb.Cmu0``
     - ``0.5562``
     - Finite, > 0: :math:`C_{\mu 0}`; also used by the wall k and the inflow profile's equilibrium :math:`k`.
   * - ``turb.Cb``
     - ``0.35``
     - Finite, > 0: :math:`C_b` of the stable and unstable length scales.
   * - ``turb.Rt_crit``
     - ``-1``
     - Finite, ≤ 0: :math:`R_t` is smoothed below this value.
   * - ``turb.Rt_min``
     - ``-3``
     - Finite, < ``turb.Rt_crit`` and > −3.598: the smoothed :math:`R_t` stays above it. :math:`c_\mu` has its pole at
       −3.5987 and grows without bound as ``turb.Rt_min`` approaches it (:doc:`turbulence_model`).
   * - ``turb.max_geom_lscale``
     - 30 m
     - Finite, > 0: the harmonic cap of :math:`l_g` (about :math:`\kappa \times 0.1\delta`, the :math:`l_g` at the
       top of the surface layer, :math:`\delta` the boundary-layer depth);
       also the cap of the inflow profile's length with ``bc.inflow_tke_length = closure``.
   * - ``turb.tke_floor``
     - ``0``
     - 0 (machine epsilon) or ≥ machine epsilon (m²/s²): lower bound on :math:`k`.
   * - ``turb.sigma_k``
     - ``1``
     - Finite, > 0: :math:`\gamma_k = \nu_\mathrm{mol} + \nu_t/\sigma_k`.
   * - ``turb.nu_molecular``, ``turb.alpha_molecular``
     - ``0``, ``0``
     - Finite, ≥ 0 (m²/s): molecular viscosity and heat diffusivity. ``turb.nu_molecular`` > 0 with ``turb.closure = none``.
   * - ``turb.kappa``
     - ``0.41``
     - Finite, > 0: von Kármán constant; also used by the wall law and the inflow profile's :math:`\nu_t`.
   * - ``turb.nut_relax``
     - ``1``
     - In (0, 1]: under-relaxation of :math:`\nu_t` (1: none).
   * - ``turb.frozen_switch``
     - ``1``
     - 0 or 1: freeze the k-floor decisions of the closure from the switch to Newton to the end of the Newton stage
       (:doc:`solution_algorithm`). With ``steady.stage = simplec`` it aborts as a key with no effect.
   * - ``turb.ri_gate``
     - ``0``
     - 0 (off) or finite, ≥ 1e-6 and ≤ 1: :math:`Ri_0` of the stratification gate
       :math:`N^2 \to N^2\,Ri^2/(Ri^2 + Ri_0^2)` in the closure (:doc:`turbulence_model`).
   * - ``turb.realizability``
     - ``0``
     - 0 (off) or finite, ≥ 0.01 and ≤ 10: :math:`C_r` of the cap :math:`\nu_t \le C_r k/S`, applied to the new
       :math:`\nu_t` before ``turb.nut_relax`` (:ref:`nut-cap`). 0.309 is about :math:`C_{\mu 0}^2` (the log-layer
       value) for the default ``turb.Cmu0``.
   * - ``turb.face_mean``
     - ``arithmetic``
     - ``arithmetic`` or ``logarithmic`` (:doc:`discretization`).
   * - ``turb.wall_distance``
     - ``terrain_height``
     - ``terrain_height`` or ``tucker``.

``convection``
--------------

See :doc:`discretization`. Read with the SIMPLEC stage only. :math:`k` always uses first-order upwinding.

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Input
     - Default
     - Range and meaning
   * - ``convection.limiter``
     - ``van_leer``
     - ``van_leer``, ``van_albada`` or ``upwind``, for momentum and θ.
   * - ``convection.smooth_eps_velocity``
     - ``1e-5``
     - > 0 (m/s), with a finite square: the :math:`\epsilon` of ``van_albada`` for momentum.
   * - ``convection.smooth_eps_theta``
     - ``1e-5``
     - The same for θ (K).
   * - ``convection.hybrid_sigma_max``
     - ``0``
     - 0 (off) or finite, > 0 and ≤ 1: the largest blend :math:`\sigma` of the momentum and θ convection towards
       first-order upwind (:ref:`hybrid-convection`). Needs ``turb.realizability`` > 0; not with
       ``convection.limiter = upwind``, ``steady.stage = simplec_newton`` or ``io.restart``.
   * - ``convection.hybrid_height``
     - none
     - Finite, > 0 (m), required with ``hybrid_sigma_max`` > 0: the blend acts below this height above the ground
       and tapers to 0 at 1.5 times it (recommended: 2.5 terrain heights).
   * - ``convection.hybrid_relax``
     - ``0.1``
     - In (0, 1]: under-relaxation of :math:`\sigma` per iteration. Refused unless ``hybrid_sigma_max`` > 0.
   * - ``convection.hybrid_freeze``
     - ``1e-2``
     - 0 (never) or in (0, 1): :math:`\sigma` is frozen once the relative residual falls below it. Refused unless
       ``hybrid_sigma_max`` > 0.

``steady``
----------

See :doc:`solution_algorithm`. With ``steady.stage = none`` the application builds the mesh and initial field,
writes them and stops; the other ``steady.*`` keys are then refused.

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Input
     - Default
     - Range and meaning
   * - ``steady.stage``
     - ``none``
     - ``none``, ``simplec``, ``simplec_newton`` or ``unsteady``. ``simplec_newton``: SIMPLEC, then the Newton stage
       once the switch rules allow it, with the fallback to SIMPLEC (:doc:`solution_algorithm`). ``unsteady``:
       backward-Euler physical steps, each converged by SIMPLEC pseudo-time iterations (`unsteady`_).
   * - ``steady.max_iterations``
     - ``10000``
     - ≥ 1. Iteration limit over the whole run, restarts included. With grid sequencing it applies to each level.
   * - ``steady.residual_tol``
     - ``5e-5``
     - Finite, > 0. Every scaled residual must be below it. The scales divide by the depth of the layers, not
       ``mesh.prob_hi`` z (:ref:`convergence-criteria`).
   * - ``steady.residual_drop_orders``
     - ``4``
     - ≥ 0 and ≤ the largest decimal exponent of the precision (308 double, 38 single). Every residual must also
       drop this many orders below its reference. Round-off limits the reachable drop to about 15 orders in double
       and 6 to 7 in single; a larger value runs to ``steady.max_iterations``.
   * - ``steady.residual_ref_iterations``
     - ``5``
     - ≥ 1. The reference is the largest residual over this many first iterations. During them the CFL stays at
       ``steady.cfl0`` (for grid sequencing see ``steady.sequencing_hold``).
   * - ``steady.monitor_rel_tol``
     - ``1e-3``
     - Finite, > 0. Relative tolerance of the monitors.
   * - ``steady.monitor_window``
     - ``50``
     - ≥ 2. Monitor window (iterations).
   * - ``steady.confidence``
     - ``0``
     - 0 (off) or finite, > 0 and ≤ 1: the run is also steady when this share of the cells meet the residual
       criterion cell by cell and every monitor criterion is met (:ref:`confidence stop <confidence-stop>`). SIMPLEC stage only.
   * - ``steady.confidence_height``
     - ``0``
     - 0 (anywhere) or finite, > 0 (m), only with ``steady.confidence`` > 0: every cell that does not meet the
       criterion must also lie within this height above the ground.
   * - ``steady.cfl0``
     - ``1``
     - Finite, > 0. CFL number of the first iterations.
   * - ``steady.cfl_max``
     - ``8`` with ``lusgs``, ``10`` with ``zjacobi``
     - Finite, ≥ ``steady.cfl0``. Largest CFL number. The default depends on ``steady.sweep`` (the sweep and
       CFL-cap study on Askervein, :doc:`results`); a given value is used with either sweep. With the step floor off
       (``steady.dtau_floor = 0``) a dense canopy needs a lower cap (:doc:`canopy`, steep terrain and the boundaries).
   * - ``steady.ser_exponent``
     - ``1``
     - Finite, ≥ 0. Exponent :math:`p` of the CFL ramp; 0 holds the CFL at ``cfl0``.
   * - ``steady.ser_max_growth``
     - ``2``
     - Finite, ≥ 1. Largest CFL growth factor per iteration.
   * - ``steady.dtau_floor``
     - ``0.1``
     - Finite, ≥ 0. The local pseudo-time rate is kept at or above this fraction of the level's fastest convective
       rate, :math:`U_{\max} / \min(\Delta x, \Delta y)`, so a cell at rest steps at most ``1/dtau_floor`` times
       as the fastest cell does (:doc:`solution_algorithm`, local time step). 0 is the plain local step, the form
       before the floor (a cell at rest then steps ``cfl`` / 1e-12 s).
   * - ``steady.switch_window``
     - ``20``
     - ≥ 1. Iterations below ``switch_residual`` before a switch to the next stage; no other switch within this many
       iterations of the last one (hysteresis), except the Newton stage's own hand-back.
   * - ``steady.switch_residual``
     - ``1e-4``
     - > 0 and ≤ 1. Relative residual at which the Newton stage starts (it was ``1e-2`` and unused
       before the Newton stage).
   * - ``steady.switch_back_growth``
     - ``10``
     - Finite, > 1. Residual growth over a window that returns the run to the previous stage.
   * - ``steady.switch_monitor_tol``
     - ``1e-2``
     - 0 (off), or > 0 and < 1. The switch to Newton also needs every monitor to have changed by at most this fraction
       of max(\|value\|, scale) over the window.
   * - ``steady.newton_fallback_failures``
     - ``2``
     - Integer, ≥ 1. Consecutive failed Newton steps (every retry used) after which the run goes back to SIMPLEC.
       After three such fallbacks the run stays in SIMPLEC.
   * - ``steady.newton_cfl_max``
     - ``128``
     - Finite, ≥ ``steady.cfl_max`` (with ``simplec_newton``). The Newton stage's CFL cap: from ``steady.cfl_max`` the
       CFL is multiplied by ``steady.ser_max_growth`` per accepted step up to this cap, then (with
       ``steady.newton_dtau_inf_at_cap = 1``) :math:`\Delta\tau` is infinite. The preconditioner of an infinite
       :math:`\Delta\tau` uses this CFL.
   * - ``steady.sweep``
     - ``lusgs`` (CPU build), ``zjacobi`` (GPU build)
     - ``zjacobi`` (z-line Jacobi) or ``lusgs`` (line LU-SGS). The default depends on the build (the sweep and
       CFL-cap study on Askervein, :doc:`results`): a GPU build keeps ``zjacobi`` until the sweeps are measured on H100.
       ``zjacobi`` reproduces the solver used before line LU-SGS became the CPU default. A restart takes the current default, not the sweep of the run
       that wrote the checkpoint: give ``steady.sweep = zjacobi`` to continue a run that used the old default. The
       ``steady.line_*`` keys need ``lusgs``, so a deck that gives them should also give ``steady.sweep``.
   * - ``steady.zline_sweeps``
     - ``2``
     - Sweeps per iteration: 2 to 100 with ``zjacobi``, 1 to 100 with ``lusgs``. With ``zjacobi`` the stable CFL
       grows with the sweeps: 0.9 with 2 (the default), 2.4 with 4, 3.7 with 8, 6.1 with 16 (a Fourier analysis,
       :doc:`results`). The ``zjacobi`` default cap of 10 is above the 2-sweep value; see :doc:`results` for the
       runs that set it.
   * - ``steady.line_ghosts``
     - ``lagged``
     - ``lagged`` or ``zero``; ``lusgs`` only. ``lagged``: neighbours across a box face take the previous pass's
       values. ``zero``: each box sweeps on its own, without communication.
   * - ``steady.line_overlap``
     - ``0``
     - Integer from 0 to 8, and at most the narrowest box width (on every hierarchy level with grid sequencing);
       > 0 only with ``line_ghosts = zero``. Columns by which each box's swept region is grown.
   * - ``steady.line_order``
     - ``sequential``
     - ``sequential`` or ``diagonal``; ``lusgs`` only. ``diagonal`` runs one column per GPU thread with the same
       result.
   * - ``steady.line_alternate``
     - ``0``
     - 0 or 1; ``lusgs`` only, with ``zline_sweeps`` ≥ 2. 1 reverses the y direction on every second sweep.
   * - ``steady.sequencing``
     - ``0``
     - 0 or 1. 1: grid sequencing over the ``mesh.hierarchy_levels`` levels (:doc:`solution_algorithm`). Needs
       ``mesh.hierarchy_levels`` ≥ 2 and at least 2 cells in x and y on the coarsest level. ``init.*`` then acts on
       the coarsest level only, and ``steady.max_iterations`` applies to each level. On a restart the coarse levels
       are skipped and the ``steady.sequencing*`` keys have no effect (the log says so).
   * - ``steady.sequencing_drop_orders``
     - ``2``
     - Finite, > 0 and ≤ ``steady.residual_drop_orders``; with ``sequencing = 1`` only. A coarse level stops once
       every residual is this many orders below its reference, or below ``steady.residual_tol``.
   * - ``steady.sequencing_monitor_tol``
     - ``0`` (off)
     - 0, or finite, > 0 and < 1; with ``sequencing = 1`` only. A coarse level also stops once every monitor changed
       by at most this fraction over ``steady.monitor_window`` iterations.
   * - ``steady.sequencing_hold``
     - ``0``
     - 0 or 1; with ``sequencing = 1`` only. 0: the seeded levels (all but the coarsest) skip the reference hold and
       ramp the CFL from their first iteration. 1: every level holds the CFL at ``steady.cfl0`` for its first
       ``steady.residual_ref_iterations`` iterations. A restart keeps the setting of the checkpoint.
   * - ``steady.sequencing_cfl_start``
     - ``cfl0``
     - ``cfl0`` or ``previous``; with ``sequencing = 1`` only. The CFL of a seeded level's first iteration:
       ``steady.cfl0``, or the CFL the previous level reached (at least ``steady.cfl0``). ``previous`` needs
       ``steady.sequencing_hold = 0``.
   * - ``steady.machine``
     - ``local``
     - ``local``, ``kestrel-cpu`` or ``kestrel-gpu``: allocation-unit rate (0, 10, 100 AU per node-hour).
   * - ``steady.nodes``
     - ``1``
     - ≥ 1. Nodes of the run, for the allocation units.
   * - ``steady.newton_gmres_restart``
     - ``30``
     - Integer from 5 to 200. GMRES restart length of the Newton step.
   * - ``steady.newton_fd_eps_scale``
     - :math:`\sqrt{\epsilon_\mathrm{mach}}`
     - ≥ :math:`\epsilon_\mathrm{mach}^{3/4}` and ≤ 1e-2. Relative step of the Jacobian–vector product.
   * - ``steady.newton_krylov_rtol``
     - ``1e-2``
     - > 0 and < 1. GMRES tolerance in each Newton step.
   * - ``steady.newton_mem_budget_gb``
     - ``70``
     - > 0 and ≤ 80. Krylov memory per rank (GB); the run stops if the estimate exceeds it.
   * - ``steady.newton_gmres_max_iters``
     - ``2000``
     - Integer, ≥ 1. GMRES iteration limit of one Newton step (all restarts). A solve that reaches it is unconverged;
       its update goes through the usual usability rule.
   * - ``steady.newton_preconditioner``
     - ``simplec``
     - ``identity`` or ``simplec``; any other value aborts and names the key. The preconditioner of each GMRES
       solve of the Newton step (:doc:`solution_algorithm`). ``simplec``: one SIMPLEC iteration with
       ``steady.sweep`` and ``steady.zline_sweeps``. ``identity``: none.
   * - ``steady.newton_armijo_c``
     - ``1e-4``
     - > 0 and < 0.5. Armijo constant :math:`c` of the continuation's line search.
   * - ``steady.newton_backtrack``
     - ``0.5``
     - > 0 and < 1. Backtracking factor :math:`\beta`: :math:`\lambda = 1, \beta, \beta^2, \ldots`
   * - ``steady.newton_lambda_min``
     - ``1e-2``
     - > 0 and < 1. The smallest :math:`\lambda` tried; below it the attempt is rejected.
   * - ``steady.newton_cfl_cut``
     - ``0.25``
     - > 0 and < 1. A rejected or skipped attempt multiplies the CFL by this factor (from infinite
       :math:`\Delta\tau`: ``steady.newton_cfl_max`` times it) and retries.
   * - ``steady.newton_max_retries``
     - ``4``
     - Integer, ≥ 1. Retries of one continuation step; after the last the step returns *failed*.
   * - ``steady.newton_dtau_inf_at_cap``
     - ``1``
     - 0 or 1. 1: once the CFL reaches ``steady.newton_cfl_max``, the next steps run without a pseudo-time term
       (:math:`D = 0`, plain Newton) until an attempt fails. 0: the CFL stays at the cap.
   * - ``steady.newton_ptc_growth``
     - ``1``
     - Finite, ≥ 1. At a finite :math:`\Delta\tau`, when no :math:`\lambda` meets the Armijo condition, accept the
       full step if :math:`|R|` grows by at most this factor. 1: off. Not used at infinite :math:`\Delta\tau`.

The ``steady.newton_*`` and ``steady.switch_*`` keys are read with ``steady.stage = simplec_newton``; with
``simplec`` they abort as keys with no effect. Their defaults are provisional. The last seven rows are the pseudo-transient
continuation of :doc:`solution_algorithm`. In the Newton stage its CFL ramps from ``steady.cfl_max`` to
``steady.newton_cfl_max`` with the SER keys; the SER law's residual term takes ``steady.newton_cfl_max`` as its
starting CFL, so ``steady.ser_exponent`` rarely limits the ramp there.

``unsteady``
------------

See :ref:`time-driver`. Read with ``steady.stage = unsteady``; without it every ``unsteady.*`` key aborts as a key
with no effect. In an unsteady run the steady stop does not apply, so ``steady.max_iterations``,
``steady.residual_drop_orders``, ``steady.monitor_rel_tol``, ``steady.monitor_window``, ``steady.confidence``,
``steady.confidence_height`` and ``io.iteration_log`` abort as keys with no effect; ``steady.residual_tol`` is the
floor of the inner stop. ``steady.residual_ref_iterations`` sets how many inner iterations, counted over the run,
hold the CFL at ``steady.cfl0``; the other ``steady.*`` keys of the SIMPLEC stage act as in a steady run.
``steady.sequencing = 1`` and ``convection.hybrid_sigma_max`` > 0 are not supported yet and abort.
A run may start from a steady checkpoint (``io.restart``): its physical time starts at 0 and the driver's history
(CFL, residual reference, hold) continues from the checkpoint (:ref:`time-driver`). From a checkpoint of an unsteady
run it continues that run: its time, step count and averages. The outputs (averages, probes, plotfiles and checkpoints
during the run) and the restart are described in :ref:`unsteady-outputs`.

.. list-table::
   :header-rows: 1
   :widths: 34 12 54

   * - Input
     - Default
     - Range and meaning
   * - ``unsteady.stop_time``
     - —
     - Required; finite and > 0 (s). The physical time the run ends at; the last step is shortened to end there.
   * - ``unsteady.dt``
     - —
     - Finite and > 0 with a finite 1/dt (s). The physical step. Required unless ``unsteady.cfl`` > 0; refused with it.
       The last step takes the time left when that is within 0.1 % of dt.
   * - ``unsteady.cfl``
     - ``0``
     - 0 (use ``unsteady.dt``) or finite and > 0. The target convective CFL: each step's :math:`\Delta t` is this
       times the inverse of the largest convective rate over the cells. A field at rest has no such rate and stops
       the run at the first step.
   * - ``unsteady.inner_drop_orders``
     - ``3``
     - Finite, > 0 and ≤ 15. A step ends when every scaled residual is at or below :math:`10^{-d}` times its value at
       the step's start, or below ``steady.residual_tol``. 3 keeps the algebraic error below the time error when the
       inner iteration converges slowly; 1 does not (:ref:`time-driver`).
   * - ``unsteady.max_inner``
     - ``50``
     - ≥ 1. Inner iterations of a step at most; a step that reaches it is reported and the run goes on (it does not
       fail the run).
   * - ``unsteady.average_start``
     - ``0``
     - Finite, ≥ 0 and < ``unsteady.stop_time`` (s). The running mean and variance of :math:`u, v, w, \theta, k`
       (cell-centred) are taken over the time after it, each step weighted by its part after this time. They are in
       every plotfile as ``u_mean`` … ``tke_mean`` and ``u_var`` … ``tke_var`` (0 before any time is averaged). A
       restart must give the checkpoint's value.
   * - ``unsteady.plot_interval``
     - ``0``
     - ≥ 0. A plotfile after every so many physical steps, counted from the run's start (a restart counts on), and
       the final one at the end. 0: the final one only. Needs ``io.write_plotfile = 1``. The plotfiles are numbered
       by the total of inner iterations, as the final one.
   * - ``unsteady.checkpoint_interval``
     - ``0``
     - ≥ 0. A checkpoint after every so many physical steps, as ``unsteady.plot_interval``, and the final one at the
       end. Needs ``io.write_checkpoint = 1``. A run that diverges writes no final checkpoint.
   * - ``unsteady.probes_file``
     - —
     - A file of probe points, one ``x y h`` per line (m; ``#`` starts a comment), :math:`h` > 0 the height above the
       ground. Each probe takes the cell of its column whose centre is nearest :math:`h` (the lower one on a tie);
       a point outside the domain or above the column's top aborts. The values of :math:`u, v, w, \theta, k` at the
       start and after each step go to ``<io.plot_file>_probes.csv``; a restart keeps its rows up to the
       checkpoint's time.

``simplec``
-----------

See :doc:`solution_algorithm`. Read with ``steady.stage = simplec`` or ``simplec_newton``.

.. list-table::
   :header-rows: 1
   :widths: 34 14 52

   * - Input
     - Default
     - Range and meaning
   * - ``simplec.pressure_reltol``
     - ``1e-3``
     - ≥ 10 machine epsilon and < 1: relative tolerance of the pressure-correction solve. Looser values can be
       unstable in stratified flow.
   * - ``simplec.alpha_p``
     - ``1``
     - > 0 and ≤ 1: pressure under-relaxation.
   * - ``simplec.rotation_relax``
     - ``0.3``
     - > 0 and ≤ 1: step of the point-implicit rotation, in inertial times. 1 often converges faster (the hill
       cases use it).
   * - ``simplec.w_vnn_cap``
     - ``0`` (off)
     - 0, or finite and ≥ 16: the cap :math:`V` on the vertical VNN number of the :math:`w` rows. Changes the path,
       not the answer. May change across a restart.
   * - ``simplec.w_vnn_cap_slope_margin``
     - ``4``
     - 0 (off), or finite and ≥ 1; refused unless ``simplec.w_vnn_cap`` > 0. The safety factor :math:`s` that lowers the
       cap on sloping ground so the pressure operator stays definite (:doc:`solution_algorithm`). 0 is the cap of
       #32, which can stall where thin cells meet steep slopes; the run then warns over terrain when the mesh
       aspect ratio is above the tested 50.78. Changes the path, not the answer. May change across a restart.
   * - ``simplec.theta``
     - ``solve``
     - ``solve`` or ``frozen``. ``frozen``: θ is the inflow profile's θ by height above the ground on every sequenced
       level and is not solved (in a canopy with ``canopy.heat_frozen = blend``, blended towards the leaves' θ,
       :doc:`canopy`); the momentum equations carry no buoyancy force, the closure keeps :math:`N^2`
       (:doc:`solution_algorithm`, prescribed potential temperature). Needs grid sequencing (``steady.sequencing = 1``,
       ``mesh.hierarchy_levels`` > 1); refused with ``steady.stage = simplec_newton``, with ``io.restart`` and with
       ``init.first_guess = column`` and a non-neutral ``init.column_stability``. With ``frozen``, ``theta`` in
       ``bc.rayleigh_fields`` and ``bc.zhi_theta_grad`` have no effect.

``solver``
----------

See :doc:`solution_algorithm`. Read when ``init.first_guess`` is not ``none`` or with the SIMPLEC stage; refused
otherwise.

.. list-table::
   :header-rows: 1
   :widths: 30 14 56

   * - Input
     - Default
     - Range and meaning
   * - ``solver.poisson_backend``
     - ``gmres_mlmg``
     - ``gmres_mlmg`` (GMRES with multigrid), ``mlmg`` (multigrid alone) or ``gmres_fft`` (ERF's operator with an
       FFT preconditioner).
   * - ``solver.poisson_reltol``
     - ``1e-10`` (``1e-5`` single)
     - ≥ 10 machine epsilon and < 1. Converged when the residual is below this fraction of the initial one. In
       single precision the rounding floor of the residual is about :math:`\varepsilon N_z^2/2`
       (:math:`\varepsilon` the machine epsilon, :math:`N_z` the vertical cells; 2e-4 on 64 layers), so the default
       1e-5 lies below it once :math:`N_z` > 13; use about :math:`4\varepsilon N_z^2` there.
   * - ``solver.poisson_abstol``
     - ``0``
     - Finite, ≥ 0. Converged when the residual is below this value.
   * - ``solver.poisson_max_iter``
     - ``200``
     - ≥ 1. GMRES iterations, or multigrid cycles with ``mlmg``; the run stops when it is reached.
   * - ``solver.poisson_restart``
     - ``50``
     - ≥ 1. GMRES restart length; not with ``mlmg``.
   * - ``solver.poisson_precond``
     - ``fft``
     - ``fft`` or ``none``; only with ``gmres_fft``.
   * - ``solver.poisson_verbose``
     - ``0``
     - ≥ 0. 1 prints each solve; 2 also prints the residual at each GMRES restart.

``io``
------

.. list-table::
   :header-rows: 1
   :widths: 24 12 64

   * - Input
     - Default
     - Range and meaning
   * - ``io.plot_file``
     - ``plt``
     - Plotfile prefix; not empty.
   * - ``io.chk_file``
     - ``chk``
     - Checkpoint prefix; not empty.
   * - ``io.write_plotfile``
     - ``1``
     - 0 or 1: write the plotfile at the end of the run.
   * - ``io.write_checkpoint``
     - ``1``
     - 0 or 1: write the checkpoint at the end of the run (and, in an unsteady run, every
       ``unsteady.checkpoint_interval`` steps).
   * - ``io.turb_plot_file``
     - none
     - Prefix of the turbulence plotfile with ERF's variable names (:doc:`output`). Needs the SIMPLEC stage.
   * - ``io.iteration_log``
     - ``table``
     - ``table`` or ``full``: the per-iteration log (:doc:`output`). ``table``: aligned columns (iteration, CFL,
       line solver, the u, v, w, tke, theta and p residuals, r). ``full``: one ``gust: it`` line per row with
       every diagnostic (the log before this input).
   * - ``io.column_guess_file``
     - ``column_guess.dat``
     - Not empty. The file the converged column's table (``z u v theta k``, full precision) is written to, on the
       I/O rank; it reads back as a ``bc.inflow_file``. A table that cannot be written gives a warning, not a stop.
       Refused without ``init.first_guess = column`` (no effect).
   * - ``io.restart``
     - none
     - Checkpoint directory to restart from (:doc:`running`). A run that would write a checkpoint into this
       directory (``io.chk_file`` and the step give the same directory, however the paths are spelled) stops.

``refavg`` (reference-averaging tool)
-------------------------------------

The separate tool ``gust_refavg`` (:doc:`references`) reads only the ``refavg`` namespace.

.. list-table::
   :header-rows: 1
   :widths: 26 12 62

   * - Input
     - Default
     - Range and meaning
   * - ``refavg.plotfiles``
     - required
     - At least 3 single-level plotfiles with the same variables, domain and box layout, in increasing time.
   * - ``refavg.output``
     - required
     - Archive directory to create; must not exist.
   * - ``refavg.case``
     - required
     - Case name, one word.
   * - ``refavg.erf_sha``
     - required
     - Full 40-character ERF commit, lower-case hexadecimal.
   * - ``refavg.machine``
     - ``local``
     - ``kestrel-cpu``, ``kestrel-gpu`` or ``local``.
   * - ``refavg.nodes``
     - ``1``
     - ≥ 1.
   * - ``refavg.ranks``
     - ``1``
     - ≥ 1. MPI ranks (CPU) or GPUs of the run.
   * - ``refavg.wall_seconds``
     - ``0``
     - Finite, ≥ 0.
   * - ``refavg.inputs_files``
     - none
     - Files to copy into the archive; each must exist, with distinct names.
   * - ``refavg.verify``
     - none
     - Instead of building, check this archive and print its metadata. The other keys are then refused.
