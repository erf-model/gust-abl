Continuous integration
======================

Every pull request runs two workflows.

``.github/workflows/ci.yml``
----------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Job
     - What it does
   * - CPU (GCC Release, GCC Debug with ``AMReX_TESTING``, Clang Release)
     - Builds with Open MPI, runs ``ctest -L "unit|smoke"`` (4-rank tests on 2 cores with ``--oversubscribe``),
       and checks that the tests left the working tree clean. Linux GCC traps floating-point exceptions.
   * - Single precision
     - Compiles only. The steady stop needs double precision.
   * - CUDA sm_90
     - Compiles only, for the H100. Hosted runners have no GPU.
   * - Windows MSVC Release
     - Builds with MSVC, no MPI and FFTW from vcpkg, and runs ``ctest -L "unit|smoke"``. MSVC does not trap
       floating-point exceptions; that check stays with Linux GCC.
   * - Windows CUDA sm_90
     - Compiles only, with ``cl`` and ``nvcc`` (CUDA 12.6, Windows SDK 10.0.22621.0).
   * - HIP gfx90a
     - Compiles only, with ROCm 6.2.4 and rocFFT, for AMD GPUs (MI250X).
   * - SYCL
     - Compiles only, with Intel ``icpx`` and oneMKL, for Intel GPUs.
   * - ``extern_erf`` drift check
     - Runs the drift check's self-tests and the check itself against a clone of ERF, and lists upstream ERF
       changes to copied files (:doc:`extraction`).
   * - Spelling and whitespace
     - ``codespell`` (settings in ``.codespellrc``) and no trailing whitespace.

``.github/workflows/docs.yml``
------------------------------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Job
     - What it does
   * - Build documentation
     - ``make -C docs html`` with warnings as errors, and uploads the HTML as a workflow artifact.
   * - Link check
     - ``make -C docs linkcheck``. It reaches third-party sites, so a failure can come from outside this
       repository.

``.github/workflows/static.yml``
--------------------------------

On every push to ``main`` (or by hand), the ``deploy`` job builds the documentation as above, with warnings as
errors, and publishes ``docs/_build/html`` to GitHub Pages. Only the built HTML is uploaded, not the repository.

Kestrel
-------

The hosted runners have no GPU, so the CUDA, HIP and SYCL jobs above only compile. ``tools/kestrel/nightly_gpu.sh``
runs the smoke tests on Kestrel H100s:

.. code-block:: bash

   tools/kestrel/nightly_gpu.sh run      # inside a GPU allocation (salloc --gpus=4)
   tools/kestrel/nightly_gpu.sh submit   # the same as a batch job (GUST_ACCOUNT=<allocation>)
   tools/kestrel/nightly_gpu.sh status   # the result of the latest run

The job clones ``erf-model/gust-abl`` at ``main`` (``GUST_NIGHTLY_REF`` for another branch) into
``GUST_NIGHTLY_WORK`` (default ``/scratch/$USER/gust-nightly``), builds it with ``build_kestrel.sh gpu``, sets one
rank per GPU (``GUST_MPIEXEC_PREFLAGS=--gpus-per-task=1``) and runs ``ctest -L smoke``. Before the tests it checks
that every smoke run of ``gust`` carries ``amrex.the_arena_init_size=0 amrex.init_snan=1
amrex.fpe_trap_invalid=1``. The smoke tests include the 2-D hill restart chain over terrain
(``rt_steady_hill2d_restart_*``) and the 2-D hill with Coriolis run twice, compared at zero tolerance
(``rt_coriolis_repeat_*``). Each run writes ``logs/<date>_<sha>/`` (``info``, ``build.log``, ``ctest.log``,
``junit.xml``, ``status``); the newest 30 are kept (``GUST_NIGHTLY_KEEP``). ``GUST_NIGHTLY_MAIL`` mails a failed
batch job.

The job is started by ``scrontab`` or by ``cron`` on a login node, not from GitHub Actions: a self-hosted runner on
Kestrel would need the HPC center's approval. Results are read on Kestrel (``status``) and are not reported to pull
requests.

The regression tests labelled ``physics`` (``ctest -L physics``) and the performance tests will run nightly on
Kestrel CPU nodes and H100s once they exist.
