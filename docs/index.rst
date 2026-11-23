GUST
====

**GUST** computes steady-state wind fields in the atmospheric boundary layer (ABL) over complex terrain. It
solves the steady Reynolds-averaged Navier–Stokes (RANS) equations with a one-equation k–L turbulence model on a
terrain-following mesh. It is a stand-alone `AMReX <https://github.com/AMReX-Codes/amrex>`__ application and runs on
CPUs and GPUs. GUST stands for Grid-sequenced Unified Steady Terrain solver; it was called ablsteady until
October 2026.

A transient ABL code such as `ERF <https://github.com/erf-model/ERF>`__ reaches a steady wind field by integrating in
time for many hours of simulated flow. GUST solves for the steady state directly. It drives the residual of the
steady equations to zero with a pseudo-time iteration, and accelerates the iteration with line solvers, a good
initial field and coarse grids.

Some flows have no steady state, such as the wake of a bluff body that sheds vortices. For these GUST also runs
unsteady RANS (``steady.stage = unsteady``): it marches physical time in backward-Euler steps and converges each
step with the same pseudo-time iteration (:ref:`time-driver`).

Beyond bare terrain, GUST models forest canopies (:doc:`canopy`) and buildings (:doc:`buildings`), and refines
parts of the domain with nests (:ref:`nests-one-run`).

GUST reuses ERF components where they fit:

- the Axell–Liungman k–L turbulence model;
- the stretched terrain-following mesh with terrain smoothing;
- the terrain Poisson solver;
- the Monin–Obukhov surface layer.

The initial field comes from a mass-consistent model in the style of WindNinja, which follows
`FastWindTerrain <https://github.com/hgopalan/FastWindTerrain>`__.

.. only:: html

   This documentation is also available as a single PDF: `gust.pdf <gust.pdf>`__.

Development status
------------------

.. list-table::
   :header-rows: 1
   :widths: 50 50

   * - Capability
     - Status
   * - Pseudo-time SIMPLEC with implicit vertical lines (z-line Jacobi or line LU-SGS), a horizontal grid
       hierarchy and grid sequencing (:doc:`solution_algorithm`, :doc:`grid_hierarchy`)
     - Available; the production steady solver
   * - Jacobian-free Newton–Krylov, preconditioned by SIMPLEC (``steady.stage = simplec_newton``)
     - Available as an option. Refused with buildings, laminar flow, a frozen :math:`\theta`, the near-ground
       upwind blend, two-way nests and an immersed ground.
   * - Unsteady RANS by dual time stepping, backward Euler (``steady.stage = unsteady``)
     - Available. First order in time. Grid sequencing, the near-ground upwind blend, nests and an immersed
       ground are refused.
   * - Forest canopy: drag, turbulence and heat (:doc:`canopy`)
     - Available
   * - Buildings by immersed forcing (:doc:`buildings`)
     - Available with the SIMPLEC and unsteady stages, z-line Jacobi sweeps and one grid level. Compared with
       measurements on a cube and a cube array.
   * - Nests in one run, one-way and two-way (:ref:`nests-one-run`, :ref:`two-way-nests`)
     - Available at refinement 2 2 1 with the steady stages. Buildings are refused in a deck with ``nest.names``;
       a patch that is run separately with ``nest`` faces accepts them.
   * - Immersed terrain on a flat or hybrid mesh (:doc:`immersed_terrain`)
     - In progress: neutral, one mesh level, laterally periodic or with inflow and outflow faces. A roughness map
       and a forest canopy are accepted.
   * - Nonlinear FAS multigrid
     - Not pursued. The grid hierarchy serves grid sequencing only.

Each page states the limits of its capability. The target machine is NLR Kestrel: CPU nodes with 104 Sapphire
Rapids cores and GPU nodes with four NVIDIA H100.

Where to start
--------------

- **New users:** build the code (:doc:`building`), then work through the tutorials, starting with
  :doc:`tutorial_column`.
- **Setting up your own case:** :doc:`case_setup` (from a YAML file), :doc:`howto` and the :doc:`inputs`.
- **A flow that does not settle:** :doc:`tutorial_unsteady`, the unsteady recipe in :doc:`howto`, then
  :ref:`time-driver`.
- **The model and the numerics:** the Mathematical model and Numerical method sections below.
- **Accuracy and speed:** :doc:`results` and :doc:`scaling`.

Contents
--------

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   building
   running
   output

.. toctree::
   :maxdepth: 1
   :caption: Tutorials

   tutorial_column
   tutorial_flat_abl
   tutorial_hill
   tutorial_askervein
   tutorial_unsteady

.. toctree::
   :maxdepth: 2
   :caption: How-to guides and reference

   howto
   case_setup
   hrrr_regimes
   hrrr_outer
   inputs

.. toctree::
   :maxdepth: 2
   :caption: Mathematical model

   governing_equations
   turbulence_model
   terrain_coordinates
   boundary_conditions
   canopy
   buildings
   immersed_terrain

.. toctree::
   :maxdepth: 2
   :caption: Numerical method

   discretization
   solution_algorithm
   initial_conditions
   grid_hierarchy

.. toctree::
   :maxdepth: 2
   :caption: Verification and validation

   test_cases
   results
   scaling
   references
   bibliography

.. toctree::
   :maxdepth: 2
   :caption: Development

   testing
   ci
   extraction
   development
