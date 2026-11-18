Tutorial 5: An unsteady run — the 2-D hill in time
==================================================

The steady solver looks for a state that no longer changes. The unsteady solver follows the flow in physical time
instead, with the same equations and the same turbulence model. This tutorial runs the hill of
:doc:`tutorial_hill` for ten minutes of physical time in about 20 seconds on two ranks.

**You will learn to:** switch a deck to the unsteady stage, choose the time step and the inner stop, read the log
of the physical steps, and find the time averages and the probe histories.

The setup
---------

``Tests/Canonical_RANS/Unsteady_Hill_2D/inputs_unsteady_hill2d`` follows the hill deck of tutorial 3. These
lines make it unsteady:

.. code-block:: text

   # Follow the flow in time
   steady.stage = unsteady

   # z-line sweeps with a CFL cap of 4 for the inner iterations
   steady.sweep        = zjacobi
   steady.cfl_max      = 4.
   steady.zline_sweeps = 16

   # Each physical step: two orders below its first residual, at most 60 inner iterations
   unsteady.inner_drop_orders = 2
   unsteady.max_inner         = 60

The keys of the steady stop (``steady.max_iterations``, ``steady.residual_drop_orders``, the monitor and the
confidence keys) are removed: an unsteady run refuses them. ``steady.residual_tol`` (default 5e-5) is the floor of each
step's inner stop. The deck leaves the end time, the step and the outputs to the command line.

Each physical step is one backward Euler step. Inside it, the SIMPLEC iterations of the steady solver run until
the step's own residual has fallen by ``unsteady.inner_drop_orders`` orders (:ref:`physical-time-term`,
:ref:`time-driver`).

Running
-------

From a scratch directory:

.. code-block:: bash

   mpirun -n 2 ../../build-rel/gust ../../Tests/Canonical_RANS/Unsteady_Hill_2D/inputs_unsteady_hill2d \
       init.first_guess=profile \
       unsteady.stop_time=600. unsteady.dt=60. unsteady.average_start=300. \
       unsteady.probes_file=../../Tests/Canonical_RANS/Unsteady_Hill_2D/probes_hill2d.txt \
       unsteady.plot_interval=5

- ``unsteady.stop_time`` and ``unsteady.dt`` give ten steps of 60 s. ``unsteady.cfl`` in place of ``unsteady.dt``
  takes the step from the convective CFL number.
- ``unsteady.average_start`` starts the running means and variances at 300 s.
- ``unsteady.probes_file`` names a text file with one probe per line: ``x y h``, with ``h`` the height above the
  ground.
- ``unsteady.plot_interval = 5`` writes a plotfile after every fifth step, in addition to the one at the end.

The run starts from the log-law profile, so the first steps follow the adjustment of that profile to the hill.

The log
-------

.. code-block:: text

   gust: unsteady: physical steps to t = 600 s, dt 60.000000 s; each step to 2 orders below its first residual (floor steady.residual_tol 5e-05) or 60 inner iterations
   gust: unsteady: 4 probes to unsteady_hill2d_plt_probes.csv
   gust: step 1: t = 60 s, dt = 60 s, 48 inner iterations, inner residual 0.0153223
   gust: step 2: t = 120 s, dt = 60 s, 48 inner iterations, inner residual 0.0157439
   gust: step 3: t = 180 s, dt = 60 s, 21 inner iterations, inner residual 0.016622
   ...
   gust: step 10: t = 600 s, dt = 60 s, 19 inner iterations, inner residual 0.0178188
   gust: unsteady: 10 physical steps to t = 600 s, 248 inner iterations (24.8 per step), 0 steps stopped at unsteady.max_inner
   gust: unsteady run ended after 248 iterations (248 over the run)
   gust: wrote unsteady_hill2d_plt00248
   gust: wrote chk00248

There is one line per physical step. The first two steps need 48 inner iterations and the later ones about 19,
since the flow changes less from step to step. The summary counts the steps that stopped at
``unsteady.max_inner`` before reaching the drop; here there are none. A run with many such steps needs a smaller
step or a larger ``unsteady.max_inner``.

The outputs
-----------

- **Plotfiles.** ``unsteady_hill2d_plt00153`` (after step 5) and ``unsteady_hill2d_plt00248`` (at the end) are
  numbered by the total of inner iterations and carry the physical time. Besides the fields of a steady run they
  hold the running means and variances from 300 s on: ``u_mean`` … ``tke_mean`` and ``u_var`` … ``tke_var``.
- **Probes.** ``unsteady_hill2d_plt_probes.csv`` has a header that gives each probe's cell and one row per step
  with :math:`u, v, w, \theta, k` at each probe:

  .. code-block:: text

     # probe 1: 0 20 50 -> 64 0 3 49.236138712241655
     t,p1_u,p1_v,p1_w,p1_theta,p1_tke,...
     0,6.4952,...
     60,6.5378,...
     120,6.5654,...
     ...
     600,6.3627,...

  (values shortened to four decimals). At the hill top, 50 m above the ground, :math:`u` changes by about 0.03 m/s
  per step at the end: after ten minutes the flow is still adjusting from the first guess.
- **Checkpoint.** ``chk00248`` holds the state, the time, the step count and the averages. A later run with
  ``io.restart = chk00248`` and a later ``unsteady.stop_time`` continues it.

Things to try
-------------

- **A longer run.** Raise ``unsteady.stop_time``. This hill has a steady state, so the changes per step die out
  and the variances over a late window become small.
- **A start from the steady state.** Run tutorial 3's z-line deck (``inputs_steady_hill2d``, the same mesh) first,
  then give its checkpoint as ``io.restart`` in place of ``init.first_guess``. The unsteady run then starts at
  :math:`t = 0` from that state (the test ``rt_unsteady_hill2d_from_steady`` does this from an early checkpoint).
- **A smaller step.** Backward Euler is first order and damps oscillations. For a flow that sheds vortices, halve
  the step and compare the amplitudes (:doc:`buildings`, Validation, has the numbers of a laminar cylinder).

What an unsteady run does not accept (grid sequencing, the near-ground upwind blend, nests) is listed in
:doc:`howto` (Run unsteady RANS). :doc:`results` compares the time mean of unsteady runs with the steady runs of
the same decks.
