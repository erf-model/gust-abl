Initial conditions
==================

The steady solver starts from an initial field. A good initial field is already close to the solution, so the solver
needs fewer iterations. GUST builds the initial field in two steps:

#. **Initial profile.** A wind profile is placed on the mesh by height above the local ground.
#. **Mass-consistent adjustment.** The field is corrected so that it conserves mass, which bends the flow over and
   around the terrain.

The method is that of mass-consistent wind models such as WindNinja (Forthofer et al. 2014) and follows
FastWindTerrain. No FastWindTerrain code is copied (it uses a Cartesian immersed-boundary mesh), but every parameter
has the same meaning. The code is in ``Source/Initializer`` and ``Source/ColumnGuess``; the inputs are in the
``init`` namespace (:doc:`inputs`).

Options
-------

``init.first_guess`` selects the initial profile:

.. list-table::
   :header-rows: 1
   :widths: 18 52 30

   * - Value
     - Initial profile
     - Recommended for
   * - ``none`` (default)
     - Fluid at rest.
     - Tests only.
   * - ``profile``
     - The inflow profile of :doc:`boundary_conditions` (log law, power law or file).
     - Cases without Coriolis forcing; flat cases with an inversion.
   * - ``column``
     - The steady solution of a single vertical column, computed first (below).
     - Terrain cases with Coriolis forcing.
   * - ``nest``
     - A parent run's checkpoint (``bc.nest_parent``) interpolated onto a patch that refines it r r 1, r = ``bc.nest_ratio`` (below).
     - A finer patch inside a coarser run.

``init.projection`` (default 1) switches the mass-consistent adjustment on or off. With 0, the solver's own pressure
correction removes the divergence in its first iteration.

Initial profile from the inflow profile
---------------------------------------

With ``init.first_guess = profile``:

- :math:`u_0, v_0` are the inflow profile at each face's height above the ground,
  :math:`h = \max(z_\mathrm{face} - z_\mathrm{ground}, h_\mathrm{min})`, where :math:`h_\mathrm{min}` is
  ``init.z_agl_min``.
- :math:`w_0 = 0`, except on the ground faces, where :math:`w_0 = u_0 z_x + v_0 z_y` (no flow through the ground).
- :math:`\theta` is the sounding ``bc.theta_profile`` at the absolute height, and :math:`p = 0`.
- :math:`k` is the profile's equilibrium value. Above :math:`z_i` = ``init.k_decay_height`` (0: no decay) it decays:

  .. math::

     k = k_\mathrm{prof}(h)\, e^{-(h - z_i)/z_i} \quad (h > z_i), \qquad k = k_\mathrm{prof}(h) \quad (h \le z_i).

Initial profile from a single column
------------------------------------

With Coriolis forcing, the steady wind turns with height (the Ekman spiral). The inflow profile has no such turning,
so the solver must first spend many iterations on the slow inertial adjustment. ``init.first_guess = column`` avoids
this: it first solves the horizontally uniform problem of the case, a single vertical column, and uses its steady
profile as the initial profile.

The column has the vertical grid of the case over flat ground. It solves for the horizontal wind, written as the
complex number :math:`u + iv`, and for :math:`k`:

.. math::

   \frac{\partial (u + iv)}{\partial t} = \frac{\partial}{\partial z}\Bigl(\nu\,\frac{\partial (u + iv)}{\partial z}\Bigr)
   - i f\,\bigl((u + iv) - (u_g + i v_g)\bigr) + F_x + i F_y, \qquad
   \frac{\partial k}{\partial t} = \frac{\partial}{\partial z}\Bigl(\gamma_k\,\frac{\partial k}{\partial z}\Bigr)
   + P_s + P_b - \varepsilon .

- In the neutral column :math:`\theta` is held at the sounding (the stratified column below solves for it).
- :math:`\nu`, :math:`\gamma_k`, :math:`P_s`, :math:`P_b`, :math:`\varepsilon` come from the turbulence model of
  :doc:`turbulence_model`.
- The ground of the neutral column uses the neutral log law; the top and the Rayleigh layer are those of the case.
  Its roughness is ``bc.surface_z0``; with a roughness map (``bc.surface_z0_file``) the column takes one value,
  ``init.column_z0``, which must then be given.

**Forest canopy.** When the run has a canopy (:doc:`canopy`) the column carries it: the layer profile :math:`f_k`
of the level's drag field (its mean over the domain's columns) enters the momentum equation as the drag
:math:`-f_k |W| W`, linearised along the wind about each step's start (the derivative diagonal of the 3-D rows; a
lagged :math:`|W|` would hand the balance speed back and forth at the column's large steps), and with ``canopy.tke``
the :math:`k` equation gains the wake production :math:`\beta_p f |W|^3`, the short-circuit sink
:math:`\beta_d f |W| k` on the Patankar diagonal and the closure's length cap :math:`L \le \alpha/f`, the same terms
as the 3-D stage's on a flat column. The first guess and an inflow column made from it then carry the forest: on
``Tests/Canonical_RANS/Canopy_Drag`` with a uniform 20 m canopy (LAI 3) the 3-D run starts from the column with
first residuals of :math:`10^{-9}` against 1 from the profile, and reaches the verdict in 1082 iterations against
2046 (``rt_canopy_column``). The column is horizontally uniform and sits on the flat levels, so the canopy must be
the same in every column and the terrain flat: on a terrain-following mesh even a canopy of uniform height gives
each column its own layer profile, and a canopy whose :math:`f` differs between columns (a forest on part of the
domain) stops the run, naming ``init.first_guess`` (use ``profile`` there, or an inflow table). Deep inside a dense
uniform canopy the column's balances are closed forms, :math:`|W| = \sqrt{F/f}` for a pressure-gradient force
:math:`F` and :math:`\beta_p |W|^3 = \beta_d |W| k + C_{\mu 0}^3 k^{3/2}/\alpha`, which the unit tests
``ColumnCanopy.*`` pin.

The column is marched to steady state with implicit (backward Euler) steps that grow from 10 s to :math:`10^5` s. It
has converged when the largest change of :math:`u, v` and :math:`k` (and :math:`\theta` in a stratified column) over a
step, relative to their scale, is below
``init.column_tol``. This typically takes 100–200 steps and a few hundredths of a second. If it has not converged
after ``init.column_max_iterations`` steps, the run stops.

The converged :math:`u, v, k` are then placed on the mesh by height above the local ground, exactly as the inflow
profile is; :math:`\theta` stays the sounding (the stratified column below can set it). The table is written to
``io.column_guess_file`` (default ``column_guess.dat`` in the run directory), in full precision, so it reads back as a
``bc.inflow_file`` exactly.
An inflow face keeps the case's inflow profile; the solver reconciles the two.

**Requirements.** The column option needs:

- a horizontal forcing: ``forcing.geostrophic = 1`` or a non-zero ``forcing.pressure_gradient``;
- the SIMPLEC stage (``steady.stage = simplec`` or ``simplec_newton``);
- a ``surface_layer`` ground with the default wall :math:`k` (``bc.surface_zref`` is not supported), neutral
  unless ``init.column_stability = heat_flux`` (then ``bc.surface_stability = heat_flux``) or ``mol`` with a
  non-zero ``init.column_mol`` (``bc.surface_stability = heat_flux`` allowed: the run takes the column's flux);
- with a non-neutral steady column, a Rayleigh layer on :math:`\theta` (below; the time-marching column needs none);
- a ``slip_wall``, ``geostrophic`` or ``no_slip_wall`` top;
- ``bc.rayleigh_fields`` either containing both ``u`` and ``v`` or neither;
- with a pressure gradient but no Coriolis force: a neutral sounding or a Rayleigh layer on ``u`` and ``v`` (the
  steady column only);
- with a forest canopy (:doc:`canopy`): flat terrain and a canopy that is the same in every column (below).

``init.k_decay_height`` cannot be used with it. Otherwise the run stops with a message naming
``init.first_guess``.

.. _column-stability:

**Stratified column.** With a non-neutral surface (``init.column_stability``, the forms of
amrTerrain's 1-D column solver) the column also solves for :math:`\theta`,

.. math::

   \frac{\partial \theta}{\partial t} = \frac{\partial}{\partial z}\Bigl(\gamma_\theta\,\frac{\partial \theta}{\partial z}\Bigr)
   - c_r(z)\,(\theta - \theta_p),

with the deck's gradient ``bc.zhi_theta_grad`` at the top and, at the ground, Monin–Obukhov similarity for a
surface temperature :math:`\theta_s` (ERF's Businger–Dyer :math:`\psi_m`, :math:`\psi_h`, or QNSE's on the stable
side with ``bc.surface_stable_functions = qnse``, :ref:`stable functions <stable-functions>`):

.. math::

   u_* = \frac{\kappa U_1}{\ln\frac{z_1 + z_0}{z_0} - \psi_m(\zeta)}, \qquad
   q = -\frac{\kappa u_*\,(\theta_1 - \theta_s)}{\ln\frac{z_1 + z_0}{z_0} - \psi_h(\zeta)}, \qquad
   \zeta = \frac{z_1}{L}, \quad \frac{1}{L} = -\frac{\kappa g q}{u_*^3 \theta_0},

- :math:`U_1`, :math:`\theta_1`: the first cell's speed and potential temperature at height :math:`z_1`;
- :math:`\zeta` from the bulk Richardson number, as ERF's surface-temperature solve:
  :math:`\zeta = \mathrm{Ri}_b\,(\ln - \psi_m)^2/(\ln - \psi_h)`,
  :math:`\mathrm{Ri}_b = z_1 g (\theta_1 - \theta_s)/(U_1^2 \theta_0)` held in
  [``init.column_rib_min``, ``init.column_rib_max``] (default ±4, ERF's), solved by ERF's fixed point
  under-relaxed by a half. Businger–Dyer's linear stable :math:`\psi` has a root only for
  :math:`\mathrm{Ri}_b < 0.2` (:math:`\zeta = \mathrm{Ri}_b \ln/(1 - 5\,\mathrm{Ri}_b)`), so with them the bound
  is at most 0.19 (:math:`\zeta = 3.8\ln`); QNSE's :math:`\phi` levels off and needs no further cap. The flux takes
  the full :math:`\theta_1 - \theta_s`, as in ERF and WRF;
- :math:`q`: the kinematic surface heat flux (K m/s, upward > 0), implicit in :math:`\theta_1`;
- the wall :math:`k` the surface layer's (Axell and Liungman 2001, Eq. 16):
  :math:`k_1 = (u_*^3 + \kappa B z_1)^{2/3}/C_{\mu 0}^2`, :math:`B = g q/\theta_0` for an upward flux only.

``surface_temperature`` fixes :math:`\theta_s` (``init.column_surface_theta``). ``mol`` and ``heat_flux`` fix
:math:`1/L` (``init.column_mol``) or :math:`q` (``bc.surface_heat_flux``), and the column finds the :math:`\theta_s`
that gives it by the secant method. The search is on :math:`1/L` for ``mol``, so it reaches every :math:`L` the
column allows, also those past its largest downward flux. :math:`\zeta` is taken at the first cell centre,
:math:`z_1/L`, as in the 3-D wall model.

**Heat budget.** At steady state the surface flux has to leave (or enter) the column somewhere. The top passes
almost nothing: its flux is the deck's gradient times a diffusivity at the floor. The Rayleigh layer on
:math:`\theta` is therefore the column's heat source or sink, and a non-neutral column requires it
(``bc.rayleigh_zdamp`` > 0 and ``theta`` in ``bc.rayleigh_fields``). The most stable column the deck allows (the
smallest :math:`L`, the largest downward flux) depends on that layer's height and strength, on the wind and on the closure's length cap
``turb.max_geom_lscale``: on ``Tests/Canonical_RANS/Column_Target`` the smallest :math:`L` is below 500 m with a
30 m cap and about 700 m with the default 100 m. A request beyond it stops the run, naming the smallest
:math:`L` (or the largest flux) reached. The numbers that follow are for a 30 m cap. On the unit test's 64 m
column (Rayleigh layer above 48 m) the smallest :math:`L` is 414 m; on a 120 m deck with the Bolund vertical grid the largest downward flux is −0.029 K m/s with the
Rayleigh layer above 80 m and −0.022 K m/s above 100 m.

The log prints the flux, :math:`L` and :math:`\theta_s` reached, and the column's total steps and columns. With
``mol`` and ``bc.surface_stability = heat_flux`` the 3-D run takes the column's flux as ``bc.surface_heat_flux``.
The wall :math:`k` is the surface layer's, :math:`(u_*^3 + \kappa B z_1)^{2/3}/C_{\mu 0}^2`; amrTerrain's column
uses the additive form :math:`(u_*^2 + (\kappa B z_1)^{2/3})/C_{\mu}^2`, up to :math:`2^{1/3}` larger for an
upward flux.

**Target wind.** ``init.column_target_u``, ``init.column_target_v`` at ``init.column_target_height`` (amrTerrain's
met-mast forcing): the geostrophic wind needed to keep a given wind at a given height under each stratification.
After each converged column the geostrophic wind and the column's Rayleigh reference wind are multiplied by
:math:`R = U_\mathrm{target}/U(z_t)` (complex: scaled and turned), the pressure gradient by :math:`R\,|R|` (the
stress goes with the square of the speed), and the column is solved again, until :math:`|R - 1| < 10^{-4}`.
:math:`U(z_t)` is MOST's profile below the first cell centre and linear in :math:`\ln(z + z_0)` above. The 3-D run
then takes that forcing. It reaches the target only if nothing holds its wind to the deck's inflow profile, so a
target is refused with a Rayleigh layer on :math:`u, v`, a sponge or an inflow face (for an inflow case, run the
column on a periodic deck and give its ``io.column_guess_file`` table as ``bc.inflow_file``), and with a stratified column
over a neutral 3-D surface.

.. _column-marching:

**Time-marching column.** ``init.column_time`` > 0 marches the same column in physical time instead: backward Euler
steps of ``init.column_dt`` (the last step ends at ``init.column_time``), :math:`\nu_t` taken from the closure at
each step without the steady column's relaxation, and no convergence test. It is the form of ERF's near-1-D column
with a stratified inflow: a column heated or cooled at a fixed rate has no steady state, so it needs no
Rayleigh layer on :math:`\theta`, and the table is its state at ``init.column_time``. The surface is fixed: the
neutral wall law (``mol`` with ``init.column_mol = 0``, :math:`\theta` held), ``surface_temperature``, or
``heat_flux``, where the flux :math:`q` = ``bc.surface_heat_flux`` enters the first cell as given (explicit, no search
for :math:`\theta_s`) and MOST gives :math:`u_*` and :math:`L` from it, as ERF's flux surface:

.. math::

   u_* = \frac{\kappa U_1}{\ln\frac{z_1 + z_0}{z_0} - \psi_m(\zeta)}, \qquad
   \zeta = -\frac{\kappa g q z_1}{\theta_0 u_*^3},

solved by the fixed point under-relaxed by a half (from the previous step's :math:`\zeta`, then from 0 if that
fails, to :math:`10^{-10}`),
:math:`\theta_* = -q/u_*`. A downward flux too large for the wind has no root and stops the run, naming
``bc.surface_heat_flux`` (at any step, and at the final state). With a slip top, ``bc.zhi_theta_grad = 0`` and no Rayleigh layer on :math:`\theta` the
column's heat content :math:`\sum \theta\,\Delta z` then grows by :math:`q t` to round-off (the unit test
``ColumnMarching.SurfaceFluxBudget``). The target wind and a non-zero ``init.column_mol`` need the steady column; the
bulk Richardson bounds are refused with ``heat_flux`` (no bulk Richardson number is solved). The log and the table
header give the surface of the final state; with ``surface_temperature`` the flux printed is the one at
``init.column_time`` only. The 3-D run is steady and has no such time: a marched :math:`\theta` with a fixed flux and
no heat sink is not a steady state of a periodic deck. The table is meant as the ``bc.inflow_file`` of a deck with an
inflow face (with ``simplec.theta = frozen``), or as a first guess alone (the example deck
stops after one iteration, ``steady.max_iterations = 1``).
On four stratified columns (240 × 10 m levels, 9 h, ``Tests/Canonical_RANS/Column_Marching``) the
column is within 9 % of ERF's bulk Richardson number between 46 and 73 m and within 1 % of its speed at 300 m, the
same for steps of 1, 10 and 60 s; a 9 h column of 1 s steps takes about 2 s. The march is first order in time and
stayed stable up to 3600 s steps; steps up to about 300 s keep the bulk Richardson number within 0.002 of the 1 s
column, while 3600 s steps damp the inertial oscillation (U(300 m) 1–3 % low).

**The first cell's k.** By default (``init.column_wall_k = fixed``) the marching column holds the first cell's
:math:`k` at the surface layer's wall value, as the steady column and the 3-D run do:

.. math::

   k_1 = \frac{\left(u_*^3 + \kappa B z_1\right)^{2/3}}{C_{\mu 0}^2}, \qquad B = \max\left(\frac{g q}{\theta_0}, 0\right)

(Axell and Liungman 2001, Eq. 16). ``init.column_wall_k = solved`` solves the first cell's :math:`k` with the
equation of the other cells instead: no :math:`k` flux through the ground, the shear production from the mean of
the wall law's gradient at the first centre (:math:`u_*/(\kappa (z_1 + z_0))` when neutral, MOST's otherwise) and the
layer difference above, the buoyancy production or destruction, the closure's dissipation. The layer difference is
smaller than the wall law's gradient, so in a neutral surface layer the solved :math:`k_1` settles below the
log-layer value :math:`u_*^2/C_{\mu 0}^2`: at 0.77 of it for :math:`z_1/z_0 = 5` and 0.74 for :math:`z_1/z_0 \ge 30`
(a discrete constant-stress layer; 0.76 in ``ColumnMarching.SolvedWallK``), with a :math:`k` minimum at the first
cell. The closer agreement with ERF's first-cell :math:`k` below is therefore partly this discretisation, not
evidence that the solved value is the more physical one. In the four columns
the solved first-cell :math:`k` is 1.2–1.3 times ERF's (the wall value is about 1.9 times), and the bulk Richardson
number moves 1–2 % of its value towards ERF's (6.5–8 % below it, from 8–9 %). The speed at 300 m and :math:`\theta`
change by at most 0.01 m/s and 0.01 K. The 3-D run keeps the wall value either way.

**Restart.** The checkpoint holds the state, not the adjusted forcing or the column's flux, and ``init.*`` keys are
refused on a restart. The log lists the inputs a restart deck needs in their place (``forcing.ug``,
``forcing.vg``, ``forcing.pressure_gradient``, ``bc.surface_heat_flux``), with all their digits.

On the Bolund case 4 deck with :math:`L` = 250 m and 7.58 m/s at 5 m, the column gives :math:`u_*` = 0.509 m/s
(0.5085 m/s measured). That agreement mainly reflects MOST at 5 m with the chosen :math:`z_0` and target. The flux,
:math:`q` = −0.039 K m/s, is sustained by the Rayleigh layer and has no measured counterpart.

The column's :math:`\theta` goes into the ``io.column_guess_file`` table (a ``bc.inflow_file`` for a run with the same
stratification). The 3-D first guess takes the column's :math:`\theta`, by height above the ground, when the 3-D run
shares the column's surface (``bc.surface_stability = heat_flux``); otherwise :math:`\theta` is the sounding.

**When to use it.** On the 2-D and 3-D hills the column start halves the iterations (:doc:`results`). On the flat
``Neutral_ABL_Flat`` case it does not converge: the column's :math:`k` at the sharp inversion base differs from the
3-D solution, and a slow oscillation of :math:`k` develops there. Use ``profile`` for that case.

Initial field from a parent run
-------------------------------

With ``init.first_guess = nest`` the run is a *patch*: a finer mesh inside a coarser *parent* run. The initial field
is the parent's checkpoint (``bc.nest_parent``) interpolated onto the patch, then the mass-consistent adjustment
below. The code is in ``Source/NestParent``.

**Placement.** The patch refines the parent by :math:`r` = ``bc.nest_ratio`` (2, the default, 4, 8 or 16) in x and y
and keeps the parent's z-levels (refinement :math:`r\,r\,1`):

- each patch cell edge lies on a parent cell face, to :math:`10^{-4}` of a parent cell;
- the patch's z inputs (``mesh.n_cell`` z, ``mesh.prob_lo`` z, ``mesh.prob_hi`` z, ``mesh.dz0``,
  ``mesh.stretch_ratio``, ``mesh.zlevels_file``) give the parent's z-levels bit for bit (the checkpoint stores them
  as hexadecimal floats);
- the patch keeps 2 parent cells from every non-periodic face of the parent, so the interpolation never reads the
  parent's boundary ghost cells. Along a periodic direction of the parent the patch may reach the parent's faces,
  and its ghost cells then read the periodic neighbours; a patch periodic in a direction must span the parent's full
  width there, from the same ``mesh.prob_lo``.

Any other placement stops the run, naming ``bc.nest_parent`` and listing each fault. The patch's terrain is its own;
it should be the parent's terrain on the finer cells.

**Interpolation.** Each field is interpolated at its own staggered point, bilinearly in x and y on each
:math:`\zeta`-level, with no vertical interpolation. For a patch cell centre :math:`i` (index from the patch's low
corner) and a parent offset :math:`o` (the parent cell under the patch's first cell), the parent position is

.. math::

   X = \frac{i + 1/2}{r} - \frac12 + o \quad\text{(cell centres)}, \qquad
   X = \frac{i}{r} + o \quad\text{(faces and nodes)},

and the value is :math:`(1 - t)\,q_I + t\,q_{I+1}` with :math:`I = \lfloor X \rfloor`, :math:`t = X - I`, in each
horizontal direction (at a face on a parent face the code takes :math:`I - 1` with :math:`t = 1`; the value is the
same). :math:`t` is a multiple of :math:`1/(2r)`; with :math:`r` a power of two, :math:`t` and :math:`1 - t` are exact
in binary and the weights are mirror-symmetric, so the values are reproducible bit for bit.
``u`` uses the face axis in x and the centre axis in y; ``v`` the reverse; ``w``,
:math:`\theta`, :math:`k` and :math:`p` the centre axes. A linear field is reproduced exactly. The parent's cells
are read on the patch's own box layout, so the result does not depend on the number of ranks.

**Values per level, not per height.** Where the patch's ground :math:`z_s` differs from the parent's interpolated
ground by :math:`\Delta z_s`, a patch point sits higher or lower than the parent point it takes its value from by about
:math:`\Delta z_s (1 - \zeta/H)`. For :math:`u, v, k`, which follow the height above the ground, this matters little.
:math:`\theta` keeps the parent's value, so in a sounding with lapse rate :math:`\Gamma` it is off by about
:math:`\Gamma\,\Delta z_s (1 - \zeta/H)` near the ground: 0.01 K in ``Nest_Patch_Hill`` (the same smooth hill on both
meshes), of order 0.1–0.3 K for a ridge that the patch resolves tens of metres higher in an inversion of 0.01 K/m.

**Ground flux.** After the interpolation, :math:`w` on the ground faces is set so that no mass crosses the ground with
the patch's own slopes (:math:`\Omega = 0`), as for ``profile``: the adjustment keeps the ground flux fixed, so the
first guess must carry it.

**Checkpoint format.** The parent checkpoint must be format 2, which stores the z-levels. A format-1 checkpoint (from
an earlier version) still restarts, but as a nest parent it stops the run; rerun the parent to write format 2.

The patch's faces use its own ``bc.*`` boundary types. With ``bc.<face> = nest`` a lateral face takes the parent's
field (:doc:`boundary_conditions`, *Nest faces*). An ``inflow`` face, a Rayleigh damping layer and a side sponge
use the ``bc.inflow_*`` profile and the sounding, not the parent's field. ``init.z_agl_min``, ``init.k_decay_height`` and ``init.obrien`` act on a profile and are refused with
``nest``. Grid
sequencing (``steady.sequencing = 1``) is refused with ``nest``: the patch starts on its own mesh from the parent.

**Frozen θ.** With ``simplec.theta = frozen`` the patch holds the parent's :math:`\theta`, interpolated as above, for
the whole run (:doc:`solution_algorithm`, prescribed potential temperature). The parent may have solved
:math:`\theta` or held it. The log gives the range and the mean of the held field. The ground ghost cells of
:math:`\theta` are extended from the patch's own first two cells, as on the sequenced levels. An ``inflow`` face of
the patch still takes its :math:`\theta` from ``bc.theta_profile``; with nest faces on all lateral faces every
lateral ghost value is the parent's, and ``bc.theta_profile`` (which the deck must still give) is not used.
``canopy.heat_frozen = blend`` is refused: the parent's :math:`\theta` may already carry the leaves' heat.
``steady.stage = unsteady`` is refused (not tested). ``Tests/Canonical_RANS/Nest_Faces`` (``rt_nest_faces_frozen``)
checks that :math:`\theta` does not change over the run and equals the parent's, interpolated.

What the patch does not take from the parent:

- **The buoyancy force.** If the parent solved :math:`\theta`, its velocity and :math:`p'` balance a force the
  patch does not carry, and the patch departs from the parent, more as :math:`U/(N h)` falls
  (:doc:`solution_algorithm`, Limits, gives the test's numbers). A parent that also held :math:`\theta` has no
  such mismatch.
- **The parent's ground.** :math:`\theta` is carried with the :math:`\zeta`-level. Where the patch's ground
  differs from the parent's by :math:`\delta h` in a domain of height :math:`H`, :math:`N^2` changes by about
  :math:`\delta h/(H - h)` (0.23 % at most in the test; 5 % for a crest 50 m sharper in a 1000 m domain).
- **The first cell's** :math:`N^2`. It is the held field's own gradient between the first two cells (the ground
  extension), not the wall-law value of a parent that solved :math:`\theta`.
- **The surface stability.** ``bc.surface_stability`` and ``bc.surface_heat_flux`` are the patch deck's own, and
  are not compared with the parent's: give the patch the parent's values. A parent that solved :math:`\theta`
  with a fixed heat flux has no steady :math:`\theta`; the patch holds the :math:`\theta` of its last iteration.

.. _nests-one-run:

Nests in one run
----------------

A deck can name nests: finer meshes that run in the same process after the outer domain. Each nest gives its name and
its two horizontal corners only (:doc:`inputs`, ``nest``):

.. code-block:: text

   nest.names   = hill edge
   nest.hill.lo = 400 300      # x y (m)
   nest.hill.hi = 1200 900
   nest.edge.lo = 1200 300
   nest.edge.hi = 1500 600

The outer domain runs to its end and writes its checkpoint. Then each nest runs in turn, in the order of
``nest.names``, as a patch of that checkpoint. A nest has no deck of its own. It takes the outer deck and replaces
what a patch must have (``Source/Nests.H``):

- the mesh: its corners in x and y, twice the outer cells there (refinement 2 2 1), the outer z-levels and
  ``mesh.max_grid_size``, one grid level; periodic only in a periodic outer direction that the nest spans whole;
- the faces: ``bc.xlo``, ``bc.xhi``, ``bc.ylo`` and ``bc.yhi = nest`` (``periodic`` where the nest is periodic), with
  ``bc.nest_parent`` the outer run's last checkpoint; the ground and the top as the outer domain's; no side sponge; the
  outer deck's Rayleigh layer without u, v and k (on those it would relax towards the inflow profile, not the outer
  field), so with the default ``bc.rayleigh_fields`` a layer on w and theta; with ``nest.damping = outer`` the outer
  deck's layer as it is; no relaxation zone;
- the start: ``init.first_guess = nest``; the outer deck's ``init.projection``, ``init.alpha_h`` and ``init.alpha_v``
  also apply to each nest's correction, while its other ``init.*`` keys and ``io.column_guess_file`` set up the outer
  first guess only; no grid sequencing; no ``steady.probes_file`` and no ``io.restart``;
- the outputs: ``io.plot_file``, ``io.chk_file`` and ``io.turb_plot_file`` with ``_<name>`` appended, for example
  ``plt_hill00096``.

The rest is the outer domain's. The terrain is the outer deck's, and a ``terrain.type = agnesi`` hill stays at the
outer domain's centre, so the nest's ground equals the outer ground on the nodes they share. The terrain mapping must
then give the same heights above the ground too: ``terrain.smoothing = 0`` (BTF) and ``2`` (Sullivan) set each node's
height from its own ground and z-level, but STF (``1``, the default) smooths each level over the mesh's whole ground,
so with nests over terrain other than ``flat`` it is refused. ERF also uses only BTF on its refined levels. The forcing is the outer
run's final forcing: ``forcing.ug``, ``forcing.vg``, the horizontal ``forcing.pressure_gradient`` and
``bc.surface_heat_flux`` after any adjustment made by a column first guess.

If the outer domain does not reach steady (``steady.max_iterations``, ``steady.residual_tol``,
``steady.residual_drop_orders``), GUST stops after the outer run and no nest runs; the message names the outer domain's
last checkpoint. If a nest does not reach steady, the other nests still run; then GUST stops with an error that names
each such nest and its last checkpoint.

The placement rules are those of a patch above: corners on the outer cell faces, at least 2 outer cells inside each
non-periodic outer face. A corner within 1e-4 cells of an outer face is moved onto it. Nests may touch but not
overlap, and a nest may not cover the whole outer domain in x and y (it would have no nest face). Every nest's inputs,
its box layout (for example ``steady.line_overlap`` against its narrowest box) and its boundary inputs are checked
before the outer domain starts, so a fault in a nest stops the run before any work is done.

With nests the outer deck needs ``steady.stage = simplec`` or ``simplec_newton`` and ``io.write_checkpoint = 1``.
Buildings are refused. With ``simplec.theta = frozen`` each nest holds the θ it starts from, the outer domain's
interpolated.

A nest's result is that of the same patch run on its own: ``Tests/Canonical_RANS/Nests_One_Process`` compares the two
bit for bit, and the outer domain with and without nests. It also checks that each nest's ground equals the outer
ground on their shared nodes, and that the nest's u, averaged to the outer cells, stays close to the outer u.

Limits of this one-way form:

- The nests do not feed back into the outer domain: each nest takes the outer field on its faces, and the outer field
  does not change. ``nest.feedback = defect`` or ``relax`` adds the feedback (:ref:`two-way-nests`,
  :ref:`two-way-relax`).
- Two nests that touch each take the outer values on their shared face, but each corrects its own net flux, so their
  normal velocities there differ.
- The outer Rayleigh layer damps u and v (relaxing towards the inflow profile); a nest's damps w and theta only, unless ``nest.damping = outer``. Near
  the top the nest's u and v may therefore differ from the outer domain's.
- A nest's convergence is measured against its own first residuals, not the outer domain's. A nest starts close to
  its steady state, so ``steady.residual_drop_orders`` may take longer to reach than ``steady.residual_tol``.
- The Rayleigh layer relaxes the nest's theta towards the deck's profile, as the outer domain's does.
- The check of each nest before the outer run uses the deck's forcing; the nest itself runs with the outer run's final
  forcing (they differ only after a column first guess adjusts it).
- A one-way nest cannot be restarted within the run. To continue one, run it as a patch on its own (a deck like
  ``Tests/Canonical_RANS/Nests_One_Process/inputs_nest``) with ``io.restart``. A joint run (``nest.feedback =
  defect`` or ``relax``) restarts from its joint record (:doc:`boundary_conditions`, Two-way nests, Restart).
- A nest is not checked against the outer sponge or the outflow face beyond the 2-cell margin; a nest placed in the
  sponge takes a damped field on its faces.
- ``terrain.nest_blend_width`` is not used by the nests (their ground already equals the outer ground on their
  faces); nor are ``steady.probes_file`` and ``io.restart``.

Mass-consistent adjustment
--------------------------

The initial profile :math:`\mathbf u_0` does not conserve mass over terrain. The adjustment finds the field
:math:`\mathbf u` closest to :math:`\mathbf u_0` that does (Sherman 1978):

.. math::

   \nabla\cdot(\Sigma\nabla\lambda) = \nabla\cdot\mathbf u_0, \qquad
   \mathbf u = \mathbf u_0 - \Sigma\nabla\lambda, \qquad
   \Sigma = \mathrm{diag}(\alpha_h^2, \alpha_h^2, \alpha_v^2) .

- :math:`\lambda`: a Lagrange multiplier (a potential), with :math:`\lambda = 0` on inflow and outflow faces and no
  flux through walls, the ground and the top.
- :math:`\alpha_h, \alpha_v`: horizontal and vertical transmissivity (``init.alpha_h``, ``init.alpha_v``).

The equation is solved by the Poisson solver of :doc:`solution_algorithm`, terrain cross terms included.

**Choosing α.** Only the ratio :math:`\alpha_v/\alpha_h` matters. A smaller ratio allows less vertical motion, so
the flow goes around steep terrain rather than over it, as in stable air. In Sherman's notation
:math:`\alpha_v/\alpha_h` is :math:`\alpha_1/\alpha_2`, WindNinja's stability parameter. Ratios of 0.1 and above are
practical; smaller ratios need more Poisson iterations.

With :math:`\alpha_v \ne \alpha_h` this one solve always uses the ``gmres_fft`` back end, whatever
``solver.poisson_backend`` says, and the log says so. The SIMPLEC stage keeps ``solver.poisson_backend``.

O'Brien correction
------------------

``init.obrien = 1`` (off by default) first sets the vertical flux :math:`\Omega_0` from column continuity, integrating
from the ground up, and removes the remaining flux at the top with a weight that grows with height (O'Brien 1970):

.. math::

   a_z\Omega_{k+1/2} = a_z\Omega_{k-1/2} - J_k D^h_k\,\Delta\zeta, \qquad
   f_k = \frac{z_k - z_\mathrm{ground}}{z_\mathrm{top} - z_\mathrm{ground}} ,

where :math:`D^h` is the horizontal part of the divergence and the top residual is removed with the weight
:math:`f_k^2`. It is meant for initial profiles taken from soundings or other models.

The resulting field
-------------------

The initial field is divergence-free to the Poisson tolerance, and :math:`\theta`, :math:`k`, :math:`p` are
consistent with the sounding. It is a starting point, not a RANS solution: it does not satisfy the momentum
equations.

A restart (``io.restart``) cannot be combined with ``init.first_guess``: the checkpoint already holds the state.
