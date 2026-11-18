Running a simulation
====================

GUST takes one input file and, optionally, extra inputs on the command line:

.. code-block:: bash

   mpirun -n 2 build-rel/gust inputs_file [key=value ...]

For a first run, follow :doc:`tutorial_column`.

What a run does
---------------

#. **Check the inputs** (:doc:`inputs`). If any check fails, the run lists every error and stops.
#. **Build the mesh** and print its quality report (:doc:`terrain_coordinates`).
#. **Build the initial field** (:doc:`initial_conditions`), or read it from a checkpoint (``io.restart``).
#. **Solve**, as ``steady.stage`` says (:doc:`solution_algorithm`):

   - ``simplec``: the SIMPLEC stage, until the convergence criteria are met or ``steady.max_iterations`` is
     reached. One row of the iteration log is printed per iteration for the first ten, then every tenth
     (:doc:`output`).
   - ``simplec_newton``: SIMPLEC, then the Newton–Krylov stage once the switch rules allow it, to the same
     criteria.
   - ``unsteady``: physical time steps from 0 to ``unsteady.stop_time``, each converged by SIMPLEC iterations
     (:ref:`time-driver`). One line is printed per step. The steady criteria do not apply.

#. **Write the output**: a plotfile, the turbulence plotfile if ``io.turb_plot_file`` is set, and a checkpoint
   (:doc:`output`). An unsteady run can also write plotfiles and checkpoints during the run.
#. **Run the nests**, if the deck names any (``nest.names``): each nest runs in the same process as a finer patch
   of the outer domain's checkpoint, and with ``nest.feedback`` other than ``none`` a joint run of the outer domain
   and the nests follows (:ref:`nests-one-run`, :ref:`two-way-nests`). The nests run only if the outer domain
   reaches its steady verdict; otherwise the run stops after the outer domain's output. A joint run writes every
   domain's output a second time when it ends.

With the default ``steady.stage = none`` the run skips step 4: it writes the mesh and initial field and stops.
This is useful for checking a new mesh or initial field.

Command-line inputs
-------------------

Inputs after the file name are added to the file's inputs. A key that the file already sets cannot be given again
on the command line: a repeated key is refused.

Restarting
----------

Set ``io.restart`` to a checkpoint directory (:doc:`howto` has a step-by-step recipe).

- The mesh and terrain inputs and the floating-point precision must match the checkpoint exactly. The box layout
  (``mesh.max_grid_size``) and the number of ranks may differ. The node heights may differ from this build's by up to
  64 ulp of the top level's height (a checkpoint from another compiler; about 1.7e-11 m on a 1000 m deep mesh in
  double precision, 8 mm in single); the log says so.
- ``init.*`` keys cannot be combined with ``io.restart``: the checkpoint holds the state.
- With grid sequencing, a restart continues level 0 only.
- A restart with the SIMPLEC stage needs a checkpoint written by the stage. The checkpoint holds the solver history
  and the turbulence model's state, so the continued run reproduces an uninterrupted one on the same ranks and box
  layout, bit for bit, over flat ground and over terrain.
- The run stops if it would write a checkpoint with the name of the one it read.
- An unsteady run may start from a steady checkpoint written by ``steady.stage = simplec`` (its time starts at 0)
  or from an unsteady checkpoint (it continues that run's time, step count and averages). A
  ``steady.stage = simplec`` run may start from an unsteady checkpoint; it ignores the averages
  (:ref:`time-driver`, :ref:`unsteady-outputs`). A ``simplec_newton`` checkpoint continues only with
  ``simplec_newton``.
- A restart must give the checkpoint's canopy inputs (:doc:`canopy`), buildings (:doc:`buildings`), ``nest`` face
  inputs (:doc:`boundary_conditions`) and terrain representation (:doc:`immersed_terrain`).

AMReX runtime flags
-------------------

AMReX reads its own ``amrex.*`` flags. The tests use these, which match the GPU CI settings:

.. code-block:: text

   amrex.the_arena_init_size=0 amrex.init_snan=1 amrex.fpe_trap_invalid=1 amrex.abort_on_unused_inputs=1

``amrex.abort_on_unused_inputs=1`` makes a misspelled ``amrex.*`` flag fail the run when it finishes (after the
output is written). AppleClang builds do not trap
floating-point exceptions; use a Linux GCC build for that check.
