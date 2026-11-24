Test cases
==========

GUST is verified on a sequence of cases of increasing complexity. Each case adds one new feature to the case
before it, so a failure points to the feature just added. A new solver or option is tested on the cases in this
order. Most cases are ERF's Canonical_RANS cases, so ERF provides a reference solution (:doc:`references`). The input
files are in ``Tests/Canonical_RANS/`` and ``Tests/Real_Terrain/``.

Steady cases
------------

.. list-table::
   :header-rows: 1
   :widths: 5 25 50 20

   * - No.
     - Case
     - What it tests
     - Test
   * - 1
     - 1-D neutral column
     - Wall model, turbulence model and vertical diffusion. The log law and the equilibrium :math:`k` are known
       exactly.
     - ``rt_steady_column``
   * - 1b
     - 1-D neutral column, immersed ground
     - The same column on a flat mesh with the ground inside a layer (:doc:`immersed_terrain`): the wall stress
       balances the force for solid fractions 0, 0.3, 0.7 and 0.9, and at 0 the column repeats the body-fitted one.
     - ``rt_immersed_column_f0``, ``_f03``, ``_f07``, ``_f09``; with ``turb.shear_stencil = log_law``:
       ``rt_log_shear_column_f0``, ``_f03``, ``_f07``, ``_f09``
   * - 1c
     - Ridge and round hill, immersed ground
     - A sloping ground immersed in a flat mesh against the terrain-following mesh (:doc:`immersed_terrain`): the
       wind along a ridge (the wall law on slopes up to 18°), across a ridge and over a round hill (crest speed and
       lee minimum).
     - ``rt_immersed_along_vertical``, ``rt_immersed_ridge_vertical``, ``rt_immersed_ridge_lusgs``;
       ``rt_immersed_hill_vertical`` (nightly)
   * - 2
     - Flat 3-D periodic (ERF ``Neutral_ABL_Flat``)
     - Horizontal directions, pressure, Coriolis forcing.
     - ``rt_steady_neutral_flat``
   * - 3
     - 2-D hill (ERF ``Neutral_Hill_2D``)
     - Terrain-following mesh, metric terms, wall distance.
     - ``rt_steady_hill2d``
   * - 4
     - 3-D hill (ERF ``Neutral_Hill_3D``)
     - Flow around as well as over a hill.
     - ``rt_steady_hill3d``
   * - 5
     - Askervein hill
     - Real terrain and field measurements (masts on lines A and AA).
     - ``rt_steady_askervein``
   * - 5b
     - Askervein, conventionally neutral
     - Case 5 with stratification aloft, Coriolis forcing and a Rayleigh layer (below).
     - ``rt_steady_askervein_cnbl``
   * - 5c
     - Askervein, ERF's setup
     - Case 5 with ERF's own inflow (power law from 210° through two faces) and Rayleigh layer on ERF's full mesh,
       compared with ERF's URANS on lines A and AA.
     - ``rt_askervein_210_erf_lines``, ``rt_askervein_210_level_history``
   * - 6
     - Bolund, Big Southern Butte
     - A steep escarpment (Bolund) and an isolated butte, with Askervein, against field measurements and the
       WindNinja results of Wagenbrenner et al. (2019). Flow separation in the lee may prevent a steady solution.
     - Comparison runs in ``Tests/Real_Terrain/WindNinja2019`` (not run by CTest); the Butte at 128 m cells with a
       log-law inflow: ``rt_sites128_butte_log_law`` (nightly)
   * - 7
     - Perdigão
     - A double ridge at 128 m cells with a target wind of 10 m/s at 100 m from the south-west, through two faces
       (``Tests/Real_Terrain/Sites_128m/Perdigao_TargetWind``). No comparison with the campaign's measurements yet.
     - ``rt_sites128_perdigao_target_wind`` (nightly), ``rt_sites128_perdigao_column``,
       ``rt_sites128_perdigao_table``

``Tests/Real_Terrain/Sites_128m`` also holds Altamont Pass (power-law inflow, ``rt_sites128_altamont_power_law``)
and the Flatirons (an inflow table from GUST's column, ``rt_sites128_flatirons_inflow_file``,
``rt_sites128_flatirons_column``, ``rt_sites128_flatirons_table``); the four nightly tests check that each case
reaches a steady state. The four sites also have hybrid-mesh decks (:doc:`immersed_terrain`), one level with
LU-SGS: ``rt_sites128_butte_log_law_hybrid``, ``rt_sites128_altamont_power_law_hybrid``,
``rt_sites128_perdigao_target_wind_hybrid`` and ``rt_sites128_flatirons_inflow_file_hybrid`` (nightly) check the
steady verdict. The column tests ``rt_sites128_perdigao_column_z0`` and ``rt_sites128_column_z0_needed``
run ``init.column_z0`` under a roughness map.

**Case 5b.** Case 5 keeps the potential temperature, gravity, Coriolis forcing at 57.2° N and the flat case's
sounding: an inversion from 700 m, 0.03 K/m above, the top raised to 1599 m and Rayleigh damping above 1200 m. The
inflow is the steady mean of a precursor column (``inputs_precursor_askervein_cnbl``), whose geostrophic wind turns
the 10 m wind onto +x with :math:`u_* = 0.6` m/s, the inflow of case 5. The wind veers with height, so the flow enters
through the west and north faces and leaves through the east and south faces. There is no ERF reference yet.

Unsteady cases
--------------

The unsteady solver (``steady.stage = unsteady``, :ref:`time-driver`) is tested on these cases.

.. list-table::
   :header-rows: 1
   :widths: 25 50 25

   * - Case
     - What it tests
     - Test
   * - 2-D hill (``Unsteady_Hill_2D``)
     - Ten steps of 60 s on 1, 2 and 4 ranks give the same steps and the same answer; a restart at a step boundary
       continues bit for bit; averages and probes.
     - ``rt_unsteady_hill2d`` (MPI builds)
   * - Laminar cylinder (``Cylinder_Laminar``)
     - At Re = 40 the unsteady run settles to the steady answer; at Re = 100 it sheds vortices
       (:doc:`buildings`, Validation).
     - Comparison runs. ``rt_cylinder_re40`` checks the steady run at Re = 40.
   * - RUSHIL Hill 5 and Askervein with ERF's setup (``URANS_Averages``)
     - The time mean of the unsteady run against the steady run of the same deck (:doc:`results`).
     - Comparison runs; ``rt_urans_mean_option`` checks the readers

Canopies, buildings, nests and immersed terrain
-----------------------------------------------

Each capability has its own cases, described on its page. Each case folder's ``README.md`` and :doc:`testing`
give the CTest names.

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Capability
     - Cases
   * - Forest canopy (:doc:`canopy`)
     - ``Canopy_Drag``, ``Homogeneous_Canopy``, ``Forest_Edge``, ``Forested_Hill``, ``Forested_Escarpment``
   * - Buildings (:doc:`buildings`)
     - ``Building_Cube``, ``Building_Courtyard``, ``Building_Hill``, ``Snyder_Lawson_Cube``,
       ``Castro_Cube_Array``, ``Cylinder_Laminar``
   * - Nests (:ref:`nests-one-run`, :ref:`two-way-nests`)
     - ``Nest_Faces``, ``Nest_Validation``, ``Nest_Shapes``, ``Nests_One_Process``, ``Nests_Two_Way``,
       ``Nests_Flat_ABL``, and ``Askervein_Nest`` and ``Columbia_Gorge`` on real terrain
   * - Immersed terrain (:doc:`immersed_terrain`)
     - ``Immersed_Column`` (case 1b above), ``Immersed_Ridge_2D``, ``Immersed_Ridge_Resolution``,
       ``Immersed_Open_Faces``, ``Immersed_Roughness_Step``, ``Immersed_Canopy``, and ``Immersed_Terrain_Cases``
       on real terrain

Transient cases
---------------

ERF's Canonical_RANS set also has two stratified flat cases. They have no steady solution, so GUST does not
use them as steady references:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Case
     - Why it is not steady
   * - ``Stable_ABL_Flat``
     - GABLS1 cools the surface throughout, and its low-level jet is an inertial oscillation.
   * - ``Convective_ABL_Flat``
     - The mixed layer deepens and warms throughout.

They can be compared through quasi-steady quantities (:math:`u_*`, jet height, boundary-layer depth).

Notes
-----

- A case with Coriolis forcing started from a uniform geostrophic wind oscillates with the inertial period
  :math:`2\pi/f` (17.5 h at :math:`f = 10^{-4}` s\ :sup:`-1`) before it settles. Its ERF reference therefore runs
  48 h (:doc:`references`).
- A case that may not be steady (case 6) is still run. If it does not meet the convergence criteria
  (:ref:`convergence-criteria`), it is reported as not steady. A steady run that stops at its iteration limit is not by
  itself a sign of an unsteady flow (:doc:`results`). The unsteady solver can follow such a flow in time, from a
  deck without grid sequencing and the upwind blend (:doc:`howto`).

Results for these cases are in :doc:`results`.
