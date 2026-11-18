Turbulence model
================

RANS equations describe the mean flow. The effect of turbulent eddies on the mean flow appears as extra mixing,
which is modelled by a turbulent viscosity :math:`\nu_t`. A turbulence model computes :math:`\nu_t` from the flow.

GUST uses the one-equation k–L model of Axell and Liungman (2001), referred to below as AL01. The model needs
two quantities:

- the turbulent kinetic energy :math:`k`, how energetic the eddies are, from its own transport equation
  (:eq:`rans`);
- the length scale :math:`L`, how large the eddies are, from an algebraic formula. Near the ground the eddies are
  small; stable stratification makes them smaller still.

The model's point relations are ERF's (``extern_erf/ERF_RANSClosure.H``). The inputs are in the ``turb`` namespace
(:doc:`inputs`).

Eddy viscosity and TKE sources
------------------------------

.. math::

   \nu_t = c_\mu\,k^{1/2}L, \qquad \nu_t' = c_\mu'\,k^{1/2}L, \qquad
   P_s = \nu_t S^2, \qquad P_b = -\nu_t' N^2, \qquad \varepsilon = C_{\mu 0}^3\,\frac{k^{3/2}}{L} .

- :math:`S^2 = 2S_{ij}S_{ij}`: the squared strain rate.
- :math:`N^2 = (g/\theta_0)\,\partial\theta/\partial z`: the squared Brunt–Väisälä frequency.
- :math:`c_\mu`, :math:`c_\mu'`: stability functions (below). The turbulent Prandtl number is
  :math:`c_\mu/c_\mu'`, equal to 1 in neutral flow.
- :math:`C_{\mu 0} = 0.5562` (``turb.Cmu0``) and :math:`\kappa = 0.41` (``turb.kappa``).
- :math:`k` is bounded below by ``turb.tke_floor`` (at least machine epsilon).

Length scale
------------

The length scale starts from the geometric length

.. math::

   l_g = \kappa\,(d + z_0),

capped harmonically at ``turb.max_geom_lscale``. Here :math:`\kappa` is the von Kármán constant, :math:`d` the wall
distance and :math:`z_0` the roughness length. Stratification modifies it:

.. math::

   \frac{1}{L^2} = \frac{1}{l_g^2} + \frac{N^2}{C_b^2 k}\quad (N^2 > 0), \qquad
   L = l_g\left(1 - \frac{C_{\mu 0}^6}{C_b^2}R_t\right)^{1/2}\quad (N^2 < 0),

with the turbulent Richardson number and the constant :math:`C_b` (``turb.Cb``, 0.35)

.. math::

   R_t = \frac{L^2 N^2}{C_{\mu 0}^6 k} .

As in ERF, the unstable length is evaluated once, from :math:`l_g`.

Stability functions
-------------------

From AL01, Eqs. 31–32:

.. math::

   c_\mu = \frac{C_{\mu 0} + 0.108\,R_t}{1 + 0.308\,R_t + 0.00837\,R_t^2}, \qquad
   c_\mu' = \frac{C_{\mu 0}}{1 + 0.277\,R_t} .

The poles of :math:`c_\mu` and :math:`c_\mu'` lie at :math:`R_t = -3.5987` and :math:`-3.6101`. Below
``turb.Rt_crit``, :math:`R_t` is smoothly mapped towards ``turb.Rt_min`` (the map of Burchard and Petersen) so that
it never reaches them; ``turb.Rt_min`` must therefore lie above :math:`-3.598`. Near the pole :math:`c_\mu` grows
without bound (1.5 at the default :math:`-3`, 7.3 at :math:`-3.5`, 78 at :math:`-3.59`).

Neutral wall layer
------------------

Because :math:`\varepsilon` contains :math:`C_{\mu 0}^3`, the model reproduces the neutral log law exactly when
:math:`l_g = \kappa(z + z_0)`, that is, where the cap ``turb.max_geom_lscale`` has no effect (set it very large to
remove it, as tutorial 1 does). In the wall layer

.. math::

   u = \frac{u_*}{\kappa}\ln\frac{z + z_0}{z_0}, \qquad k = \frac{u_*^2}{C_{\mu 0}^2}, \qquad P_s = \varepsilon

at every height, with :math:`u_*` the friction velocity.

Wall distance
-------------

``turb.wall_distance`` selects the wall distance :math:`d` to the ground. With buildings and ``buildings.wall_tke = on``
the closure takes the nearer of the ground and the nearest building wall, each with its own roughness
(:doc:`buildings`).

- ``terrain_height`` (default): the height of the cell centre above the ground of its column.
- ``tucker``: the differential wall distance of Tucker (2003),

  .. math::

     \nabla^2\phi = -1,\quad \phi = 0 \text{ on the ground}, \qquad
     d = -|\nabla\phi| + \sqrt{|\nabla\phi|^2 + 2\phi} .

  On flat ground :math:`d` is the height. Over terrain it approximates the true distance; the error is largest at
  convex crests (:doc:`results`).

Laminar flow
------------

``turb.closure = none`` turns the closure off, for laminar flow and verification cases:

.. math::

   \nu = \nu_\mathrm{mol}, \qquad \gamma_\theta = \alpha_\mathrm{mol}, \qquad \gamma_k = \nu_\mathrm{mol},
   \qquad \nu_t = P_s = P_b = \varepsilon = 0,

- :math:`\nu_\mathrm{mol}` (``turb.nu_molecular``, required > 0) and :math:`\alpha_\mathrm{mol}`
  (``turb.alpha_molecular``): the molecular viscosity and heat diffusivity.

:math:`k` is then a passive scalar: it is transported and diffused but has no sources, and nothing reads it. The
fields are constant, ghost cells included, so a laminar run costs no closure evaluation. The inputs of the k–L
closure, and the options built on it (a canopy, the surface layer, the column first guess, the Newton correction of
the SIMPLEC stage, the buildings' wall laws and wall k), are refused with ``none`` (:doc:`inputs`, ``turb``); so is
``turb.face_mean``, since every face mean of a constant gives the same value. The upwind blend
(``convection.hybrid_sigma_max``) needs ``turb.realizability``, a k–L input, so a laminar run has none. The
laminar cylinder (``Tests/Canonical_RANS/Cylinder_Laminar``, :doc:`buildings`, Validation) uses it.

Model options
-------------

- ``turb.face_mean``: how the cell diffusivities are averaged to the faces (:doc:`discretization`).
- ``turb.nut_relax``: under-relaxation of :math:`\nu_t` between iterations; ``turb.nut_relax_end`` ends it below a
  relative residual (:ref:`nut-relax-end`).
- ``turb.frozen_switch``: holds at the floor, during the Newton stage, the cells whose :math:`k` was at the floor
  when Newton started (:doc:`solution_algorithm`). The stable/unstable branch and the :math:`R_t` smoothing are not
  frozen: the length and :math:`R_t` are continuously differentiable across them (below).
- ``turb.ri_gate``: the stratification gate (below).
- ``turb.shear_stencil``: the vertical shear in :math:`S^2` near the ground (:ref:`log-shear`).

.. _log-shear:

Log-consistent shear near the ground
------------------------------------

The shear production is :math:`\nu_t S^2`, and :math:`S^2` takes the vertical gradients of :math:`u` and
:math:`v` by a central difference over two layers. Near the ground the wind follows the log law, and the central
difference of a logarithm is larger than its gradient. For :math:`u = (u_*/\kappa)\ln((d + z_0)/z_0)`, with
:math:`d_k` the height of the centre of cell :math:`k` above the ground,

.. math::

   \frac{u_{k+1} - u_{k-1}}{d_{k+1} - d_{k-1}} = \frac{\partial u}{\partial z}\bigg|_{d_k}\,
   \frac{(d_k + z_0)\,\ln\dfrac{d_{k+1} + z_0}{d_{k-1} + z_0}}{d_{k+1} - d_{k-1}} .

On uniform layers with :math:`z_0` small against the layer the factor is 1.21 in the second cell, 1.06 in the third
and 1.03 in the fourth, so :math:`S^2` is 46 %, 12 % and 6 % high there. The wall cell is not affected: its
:math:`k` is the wall value.

With ``turb.shear_stencil = log_law`` the vertical gradients of :math:`u` and :math:`v` in :math:`S^2` are

.. math::

   \frac{\partial u}{\partial z}\bigg|_{d_k} = \frac{u_{k+1} - u_{k-1}}{(d_k + z_0)\,\ln\dfrac{d_{k+1} + z_0}{d_{k-1} + z_0}}

in every cell between the wall cell and the top cell (``Source/LogShear.cpp``). It is exact for a log profile. Away
from the ground it tends to the central difference: on uniform layers of thickness :math:`\Delta` the two differ by
the factor :math:`1 + \Delta^2/(3 d^2)`. :math:`z_0` is the column's roughness length (``bc.surface_z0`` or the map
of ``bc.surface_z0_file``). The other gradients in :math:`S^2`, the momentum equations and the wall cell are as
before. The default is ``central``, so existing results do not change.

**Measured.** A neutral column driven by a geostrophic wind of 10 m/s (:math:`f_c` = 1e-4 1/s, layers of
15.625 m, :math:`z_0` = 0.1 m, ``turb.max_geom_lscale = 30``, ``turb.face_mean = logarithmic``, a slip top near
2 km), against the same column on
grids 16 and 32 times finer, extrapolated. Differences from that reference in the four cells above the wall cell
(body-fitted) or above the cut cell (immersed ground with solid fraction :math:`f`, :doc:`immersed_terrain`):

.. list-table::
   :header-rows: 1
   :widths: 22 20 20 19 19

   * - Ground
     - :math:`k`, ``central``
     - :math:`k`, ``log_law``
     - Speed, ``central``
     - Speed, ``log_law``
   * - Body-fitted
     - +4.1 to +17.3 %
     - −1.1 to −0.5 %
     - −0.8 to −0.5 %
     - −0.4 to −0.1 %
   * - Immersed, :math:`f` = 0.3
     - +4.7 to +20.4 %
     - −0.8 to −0.4 %
     - −1.0 to −0.6 %
     - −0.6 to −0.1 %
   * - Immersed, :math:`f` = 0.7
     - +5.6 to +25.5 %
     - −1.0 to +1.0 %
     - −1.2 to −0.5 %
     - −0.8 to −0.1 %
   * - Immersed, :math:`f` = 0.9
     - +6.2 to +28.0 %
     - −1.3 to +1.9 %
     - −1.3 to −0.2 %
     - −0.8 to 0.0 %
   * - Body-fitted, ``turb.face_mean = arithmetic`` (the default)
     - +4.7 to +13.7 %
     - −3.2 to +0.1 %
     - −1.4 to −1.3 %
     - −1.1 to −0.5 %

With the arithmetic face mean the cell above the wall cell keeps −3.2 % in :math:`k` (−1.2 % in the next): the
arithmetic mean of a diffusivity that grows linearly from the ground is not the one that carries a log profile's
stress. Every run reaches the steady verdict (4711 to 4785 iterations with ``log_law``). These columns are not
CTests, and their decks and scripts are not in the repository: the reference is a separate fine-grid solution. The pressure-driven columns of
``Tests/Canonical_RANS/Immersed_Column`` are (``rt_log_shear_column_*``, :doc:`immersed_terrain`).

**Limits.**

- The cut cell's speed above an immersed ground. With the length capped at 30 m the constant-stress profile is not
  a pure logarithm, :math:`u = u_*[\ln((d + z_0)/z_0)/\kappa + d/l_\mathrm{max}]`, while the wall law reads the
  reference cell's speed as logarithmic. This explains the sign and the growth with :math:`f`, not the size: the
  cell above also sits about 2 % under that profile. In the geostrophic column the cut cell's speed against the
  log law with the run's own :math:`u_*` is
  +0.2, −0.4 and −0.5 % with ``central`` and −1.1, −3.6 and −8.2 % with ``log_law`` (:math:`f` = 0.3, 0.7, 0.9;
  body-fitted wall cell: +0.9 and +0.1 %). Without the cap (``turb.max_geom_lscale = 1e12``) it is within 0.2 %
  with ``log_law`` (:doc:`immersed_terrain`, Verification). The cells above the cut cell are as in the table.
- It needs the height above the ground and its roughness length: ``bc.zlo = surface_layer`` and
  ``turb.wall_distance = terrain_height``. With buildings (``buildings.type = immersed_forcing``) the run stops:
  the closure's distance there is the distance to the nearest wall.
- The column first guess (``init.first_guess = column``) and an inflow table made from it keep the central
  difference: that column forms its own shear. The first guess and the inflow k next to the ground are then those
  of ``central``, and the 3-D solve adjusts them.
- Profiles that are not logarithmic. The formula scales the central difference by a factor below 1 on uniform
  layers: 0.83 to 0.94 in the second cell (:math:`z_0/\Delta` from 0 to 1), at least 0.94 in the third. A linear
  profile's gradient is then 17 % low in the second cell, and a stable log-linear profile's (Obukhov length
  100 m) 9 % low where the central difference is 9 % high. The buoyancy frequency keeps its central difference,
  so the gradient Richardson number of the second cell (also read by ``turb.ri_gate``) rises with ``log_law``.
  In the first full cell above a cut cell the factor falls with :math:`f` (0.57 at :math:`f` = 0.995), and on
  layers stretched by 1.3 it passes 1 a few cells up (1.006 in the sixth).
- Measured on flat, neutral columns only. It is not part of the recommended setting below. No terrain, stratified
  or canopy case has been run with it; the run is accepted there, and the profile near the ground is then not
  logarithmic.

Smoothness at neutral
---------------------

With :math:`y = l_g^2 N^2/(C_b^2 k)`, the stable length is :math:`l_g(1 + y)^{-1/2}` and the unstable one
:math:`l_g(1 - y)^{1/2}` (before the :math:`R_t` smoothing). Both are :math:`l_g(1 - y/2) + O(y^2)`, so :math:`L`,
:math:`R_t` and :math:`\nu_t` have the same value and the same first derivative in :math:`N^2` on both sides of
:math:`N^2 = 0`; only the second derivative jumps. The smoothing map of :math:`R_t` has slope 1 on both sides of
``turb.Rt_crit``. The kink of the closure is the floor :math:`k_f = \max(k, k_\mathrm{floor})`. (The neutral band
:math:`|N^2| \le` machine epsilon, where :math:`L = l_g`, is a jump in :math:`L` of relative size :math:`y/2` at
its edge: 8e-10 at :math:`k` = 1e-3, but close to 1 at :math:`k` = machine epsilon, where :math:`\nu_t` is
negligible.) Near the floor the closure is smooth but steep: at :math:`k` = 1e-6 the :math:`N^2` scale of
:math:`L`, :math:`C_b^2 k / l_g^2`, is about 4e-10 s⁻², close to the change of :math:`N^2` a finite-difference
Jacobian makes, which then misses :math:`\partial L/\partial N^2` by about 30 %.

Stratification gate
-------------------

In a nearly neutral layer, under an inversion for example, :math:`N^2` can be round-off or
numerical noise of either sign. With ``turb.ri_gate`` :math:`= Ri_0 > 0` the length :math:`L` and :math:`R_t`
(and through :math:`R_t` the stability functions and the heat diffusivity) use

.. math::

   N^2_\mathrm{eff} = N^2\,\frac{Ri^2}{Ri^2 + Ri_0^2} = N^2\,\frac{N^4}{N^4 + Ri_0^2 S^4}, \qquad Ri = \frac{N^2}{S^2},

with :math:`S^2 = 2 S_{ij}S_{ij}` (:math:`(\partial u/\partial z)^2` in simple shear, so :math:`Ri` is the gradient
Richardson number). Where the shear dominates (:math:`|Ri| \ll Ri_0`) the stratification is suppressed by
:math:`(Ri/Ri_0)^2`; in an inversion (:math:`|Ri| \gg Ri_0`) it is kept to within :math:`(Ri_0/Ri)^2`. The gate
keeps the sign of :math:`N^2`, is monotone in :math:`N^2` (its slope is at most 1.125, at :math:`|Ri| = \sqrt{3}\,
Ri_0`), is infinitely differentiable wherever :math:`S^2 > 0` (at :math:`N^2 = S^2 = 0` it is continuous but not
differentiable), and is exact at :math:`N^2 = 0`, so a neutral run is unchanged bit for bit. The buoyancy
production :math:`P_b = -\nu_t' N^2` keeps the ungated :math:`N^2`, so the buoyancy sink of :math:`k` stays that of
the heat flux :math:`\nu_t'\,\partial\theta/\partial z`; the buoyancy of the momentum equation is not gated
either.

It is a smooth, local form of the neutral switch of Kynema-SGF's ``KLAxell`` (``lengthscale_switch``: below
``meso_sponge_start``, when the deck's surface heat flux is below 1e-5, the neutral length and :math:`R_t = 0`).
That switch is a jump in value and depends on a global input: with KLAxell's defaults (``Cb_stable`` 0.25,
:math:`\lambda` = 30 m) at :math:`N^2` = 1e-4 s⁻², :math:`k` = 0.1 m²/s² and :math:`z` = 100 m it changes :math:`L`
from 7.2 m to 17.3 m. The gate changes stratified answers: in the stable
surface layer :math:`Ri` is about 0.05–0.2, so :math:`Ri_0` = 0.01 lowers :math:`N^2` there by 3.8 % (at
:math:`Ri` = 0.05) to 0.25 % (at 0.2).
It is off by default (:doc:`results`).

.. _nut-cap:

Eddy-viscosity cap for separated flow
-------------------------------------

A one-equation closure with an algebraic length keeps :math:`\nu_t` large in the shear layer
that leaves a crest, so the flow reattaches too soon behind steep terrain (no separation behind the Bolund
escarpment or the Ishihara–Zhou ridge with the default closure). With ``turb.realizability`` :math:`= C_r > 0` the
new eddy viscosity is capped before the relaxation of ``turb.nut_relax``:

.. math::

   \nu_t = w\,\min\!\Bigl(\nu_t^\mathrm{new},\ C_r\,\frac{k}{S}\Bigr) + (1 - w)\,\nu_t^\mathrm{old},
   \qquad S = \sqrt{2 S_{ij} S_{ij}},

with :math:`w` = ``turb.nut_relax`` and :math:`\nu_t^\mathrm{old}` the stored value of the last iteration. The form
is that of the SST limiter (Bradshaw's :math:`-\overline{uw}/k = a_1 = 0.31`). Durbin's (1996) time-scale bound is
:math:`\nu_t \le \alpha k/(\sqrt{6}\,|S|)` with :math:`|S|^2 = S_{ij}S_{ij}`; in the :math:`S` above it is
:math:`C_r = \alpha/\sqrt{3} \approx 0.346` for :math:`\alpha = 0.6` (0.245 = :math:`\alpha/\sqrt{6}` would be the same
formula with Durbin's :math:`|S|` taken as this :math:`S`). Unlike SST it acts everywhere, not only inside the
boundary layer. With the closure's relations, :math:`\nu_t S/k = C_{\mu 0}^2 \sqrt{P/\varepsilon}` in neutral flow, so the
cap acts where :math:`P/\varepsilon > (C_r/C_{\mu 0}^2)^2`. :math:`C_r = 0.309` is 0.1 % below :math:`C_{\mu 0}^2 = 0.30936`
for the default :math:`C_{\mu 0}`: the cap is active in an equilibrium surface layer too, where it lowers
:math:`\nu_t` by 0–3 % (most near the wall, where :math:`P/\varepsilon` is 1.04–1.05 in the first cells), and it limits
:math:`\nu_t` strongly where the strain far exceeds equilibrium, as in the shear layer leaving a crest. Applying the
cap before the relaxation lets :math:`w < 1` damp the cap's changes between iterations: an imposed two-iteration
oscillation of the cap passes with the factor :math:`w/(2 - w)` (1/3 at :math:`w = 0.5`) instead of 1. Without the
option the cap is the largest real number and the closure is unchanged bit for bit. The derivative the Newton stage
reads follows the cap (:math:`\mathrm{d}\,\mathrm{cap}/\mathrm{d}k = \mathrm{cap}/k`); the column first guess
(``init.first_guess = column``) does not apply it. Around buildings, whose edges fix their separation, the cap fades
out within ``buildings.limiter_zone`` building heights (:doc:`buildings`).

.. _nut-relax-end:

Ending the relaxation of :math:`\nu_t`
--------------------------------------

The steady states are the same for every ``turb.nut_relax``: at a fixed point
:math:`\nu_t^\mathrm{old} = \nu_t^\mathrm{new}`. Whether the iteration reaches one is not. On the Columbia Gorge outer
domain (``Tests/Real_Terrain/Columbia_Gorge``, 200 m cells) the run with ``turb.nut_relax = 0.5`` stops falling at a
relative residual of :math:`3.9 \times 10^{-3}`; with ``turb.nut_relax = 1`` it reaches the steady verdict after 164
iterations. The stall is a period-2 cycle, not a fixed point: :math:`k`, :math:`\nu_t` and the k residual change sign
every iteration and repeat every second one, so the log, which prints every tenth iteration, shows norms frozen to
seven digits. The cycle sits in a few cells of the first layers on the steep south face of the Columbia Hills (half of
the k residual's squared norm in one cell, 99 % in 12), where the cap holds; there :math:`k` alternates between 2.07
and 2.15 m²/s² and :math:`\nu_t` between 2.63 and 2.71 m²/s; no :math:`k` is clipped. The cap
:math:`C_r k/S` that sets the new :math:`\nu_t` is built from the :math:`k` at the start of the iteration, and that
:math:`k` is low when :math:`S` is high: the two swing in anti-phase (:math:`k` by ±1.9 %, :math:`S` by ±2.8 %), so
the cap swings by ±4.7 %. A cap that alternates by :math:`\pm\delta` gives a relaxed :math:`\nu_t` that alternates by
:math:`\pm\delta\,w/(2 - w)`, ±1.6 % with :math:`w = 0.5`, as observed. The cycle therefore runs through the momentum
equations (the strain) as well as through :math:`k`.

**What a one-cell model shows.** With :math:`S` held fixed, the k update alone has the gain :math:`A`; through
:math:`\nu_t` it gains the production feedback :math:`B`; :math:`h` is the stored :math:`\nu_t`, scaled by
:math:`k/(m\,\nu_t)`:

.. math::

   k^{n+1} = (A + wB)\,k^n + (1 - w)\,B\,h^n, \qquad h^{n+1} = w\,k^n + (1 - w)\,h^n,

- :math:`A = 1 - a/D`: the k update with :math:`\nu_t` held, :math:`a` the derivative of the k residual in
  :math:`k` and :math:`D` the diagonal of the update (pseudo-time, diffusion and Patankar terms, a rate);
- :math:`B = P\,m/(k\,D) \ge 0`: the change of production through :math:`\nu_t`, with :math:`P` the shear production
  and :math:`m = \partial\nu_t^\mathrm{new}/\partial k\cdot k/\nu_t^\mathrm{new}` (:math:`m = 1` where the cap holds);
- :math:`w` = ``turb.nut_relax``.

The iteration converges for :math:`A^* < A < 1 - B`, with

.. math::

   A^* = -1 - \frac{wB}{2 - w}

(:math:`-1 - B` at :math:`w = 1`, :math:`-1 - B/3` at :math:`w = 0.5`; the range is empty for
:math:`B \ge (2 - w)/(1 - w)`). At :math:`A = A^*` an eigenvalue passes through :math:`-1`: a two-iteration cycle. The
upper bound does not depend on :math:`w`. So :math:`w < 1` can only narrow the range against a cycle. For the
Patankar update of :math:`k` in neutral flow, however, :math:`A \ge -1/2`, so this model alone does not cycle; with
:math:`S` added as a third variable, random parameter samples gave no case that is stable at :math:`w = 1` and cycles
through an eigenvalue below :math:`-1` at :math:`w = 0.5`. **What it means:** the cycle needs the coupling of
:math:`k`, :math:`\nu_t` and the strain between neighbouring cells; why it persists with :math:`w = 0.5` and not
with :math:`w = 1` is not explained here. The fixed points are the same for every :math:`w`.

**The option.** With ``turb.nut_relax_end`` :math:`= r_e > 0` the relaxation ends, on each grid-sequencing level, the
first time the driver's relative residual :math:`r` (:doc:`solution_algorithm`) after a SIMPLEC iteration falls
below :math:`r_e` (Newton iterations of ``steady.stage = simplec_newton`` do not end it); from the next
iteration on :math:`w = 1` and :math:`\nu_t` is the new value bit for bit. It does not start again if :math:`r` rises.
The run prints one line when it ends. A restart replays the recorded :math:`r` values with the restart's
``turb.nut_relax_end``, so with the same value it continues with the same decision. A nest's relaxation ended in its
one-way run stays ended in the joint run. ``steady.stage = unsteady`` refuses the option (the end would fall within a
physical step). Once the relaxation has ended the Newton stage (``steady.stage = simplec_newton``) no longer hands back for
it (:doc:`solution_algorithm`). On the Columbia Gorge outer domain :math:`r_e = 0.01` ends the relaxation at the first
level-0 iteration below 0.01 (between iterations 50 and 60, where the convection blend also freezes) and the run reaches the steady
verdict after 164 iterations. The coarse levels stop at :math:`r` just below 0.01, so the relaxation ends on their
last iteration and they are unchanged.

**Limits.** After the end, :math:`\nu_t` is not damped between iterations: a cap that still alternates at the end
passes to :math:`\nu_t` in full instead of by :math:`w/(2 - w)` (three times more with :math:`w = 0.5`). On the
Columbia Gorge the run converges after the end, but nothing guarantees this in general. The rough Ishihara–Zhou
ridge with the recommended setting is a counter-example: it reaches the steady verdict after 406 iterations without
the option; with :math:`r_e = 0.01` the relaxation ends at iteration 175, :math:`r` falls to
:math:`5 \times 10^{-4}` while the :math:`k` residual grows from iteration 220 on, and from about iteration 690
:math:`r` oscillates between 0.03 and 0.13 with no steady verdict after 20000 iterations. The option is therefore not
part of the recommended setting; use it on a run that stalls in a two-iteration cycle. A run that stalls above
:math:`r_e` is not helped; the option does not detect a cycle. The default 0 keeps the relaxation to the end of the
run, so existing decks are unchanged bit for bit.

.. _hybrid-convection:

Near-ground blend towards upwind
--------------------------------

With the cap, the second-order limited convection can leave a steady run at a fixed,
unconverged state (the rough Ishihara–Zhou ridge stalls at a relative residual of 0.17 with
``van_leer`` and ``van_albada``; first-order upwind converges). With ``convection.hybrid_sigma_max``
:math:`= \sigma_\mathrm{max} > 0` the momentum and θ convection of the residual become

.. math::

   C = C_\mathrm{hi} + \sigma\,(C_\mathrm{up} - C_\mathrm{hi}), \qquad
   \sigma_t = \sigma_\mathrm{max}\,s\!\Bigl(\frac{a - 0.8}{0.2}\Bigr)\,
   s\!\Bigl(\frac{1.5 H_h - d}{0.5 H_h}\Bigr), \qquad a = \frac{\nu_t^\mathrm{new}}{C_r k / S},

- :math:`C_\mathrm{hi}`: the limited scheme of ``convection.limiter``; :math:`C_\mathrm{up}`: first-order upwind
  (the scheme of the implicit side);
- :math:`s(t) = t^2(3 - 2t)` on :math:`[0, 1]` (0 below, 1 above): a smooth step;
- :math:`a`: the ratio of the closure's eddy viscosity to the cap (:math:`a \ge 1` where the cap acts);
- :math:`d`: the height above the ground (the closure's wall distance without buildings; with them the ground's
  distance alone, :doc:`buildings`); :math:`H_h` = ``convection.hybrid_height``.

:math:`\sigma` is a cell field, interpolated to the faces of :math:`u, v, w`, updated once per iteration as
:math:`\sigma \leftarrow \sigma + \lambda(\sigma_t - \sigma)` with :math:`\lambda` = ``convection.hybrid_relax``, and
frozen once the relative residual falls below ``convection.hybrid_freeze`` (after at least :math:`2/\lambda` updates
on each grid level), so a converged state satisfies one fixed scheme. The blend is of the convection operators,
not of the face fluxes: it is not conservative where :math:`\sigma` varies (in a 1-D test with :math:`\sigma` rising
from 0 to 0.1 over a few cells, a net source of 0.2–0.4 % of the flux jump; mass is not affected). :math:`\sigma` is
not checkpointed, so ``io.restart`` with the blend aborts.

**Where it acts.** Because :math:`C_r` is at the equilibrium value, :math:`a` is 1.00–1.03 in an equilibrium surface
layer, so :math:`\sigma` reaches :math:`\sigma_\mathrm{max}` over most of the near-ground layer below :math:`H_h`, on the
windward side as well as in the lee; inside a recirculation, where :math:`k` is carried in (:math:`P < \varepsilon`),
:math:`a` can fall below 0.8 and the blend switches off. With :math:`\sigma_\mathrm{max} = 0.1` the added horizontal
numerical diffusion near the ground is of the order of the turbulent viscosity there (about 5 m²/s at 10 m on
Askervein) and changes the approach-flow speed-up by about 0.2 %. A gate that acts only above equilibrium (for
example :math:`a > 1.05`) is not in this version.

Recommended setting
~~~~~~~~~~~~~~~~~~~

For terrain with separated lees: ``turb.realizability = 0.309``, ``turb.nut_relax = 0.5``,
``convection.hybrid_sigma_max = 0.1``, ``convection.hybrid_height`` = 2.5 terrain heights and
``turb.max_geom_lscale = 100`` (on Big Southern Butte the length is uncapped with a moving lid, so that the inflow
does not drift over its 18.7 km fetch). ``turb.nut_relax_end`` (:ref:`nut-relax-end`) is not part of it. :math:`C_r = 0.309` is :math:`C_{\mu 0}^2` for the default ``turb.Cmu0``
only. On the benchmarks (``Tests/Real_Terrain/WindNinja2019``, ``IshiharaZhou``,
``RUSHIL``) every case reaches the steady verdict with it:

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - Case
     - Default closure
     - Recommended setting
   * - Ishihara–Zhou ridge, smooth (tunnel reattachment 150–200 m)
     - no reversed flow; U rms 0.80 m/s
     - reversed 42–170 m; 0.30 m/s
   * - Ishihara–Zhou ridge, rough (250–300 m)
     - reversed 59–122 m; 0.70 m/s
     - reversed 47–173 m; 0.58 m/s
   * - RUSHIL Hill 3 (tunnel 0.27a–2.17a at 0.1 H)
     - 0.64a–1.27a; 0.50 m/s
     - 0.58a–1.45a; 0.48 m/s
   * - RUSHIL Hill 5 / Hill 8
     - attached, 0.23 m/s / attached, 0.23 m/s
     - attached, 0.21 m/s / attached, 0.20 m/s
   * - Askervein, all masts (RMSE)
     - 1.24 m/s
     - 0.97 m/s
   * - Bolund cases 1 / 3 / 4, all masts (RMSE)
     - 1.46 / 0.66 / 0.89 m/s
     - 1.25 / 0.64 / 1.23 m/s

Bolund case 4 (flow from the east, the lee over the escarpment) is slower than measured in the lee with the
cap. On Big Southern Butte the lee separates only on the unsmoothed terrain with a 4 m first layer
(``Tests/Real_Terrain/WindNinja2019/inputs_bsb_srtm``: lee RMSE 2.58 m/s against 3.53 m/s on the paper-like mesh,
with 68 m first cells and smoothed terrain).

``Tests/Paper/JWEIA2026`` runs all the variants of the paper on this setting (without and with the sensor, the
limit without the blending, first-order upwinding, the sensitivity to :math:`C_r`, :math:`\sigma_\mathrm{max}`
and :math:`H_h`, and the grids) and draws its figures.

The comparison uses data made openly available by N. Wagenbrenner (windninja-tools, the Big Southern Butte
database and Wagenbrenner et al. 2019), T. Ishihara and T. Zhou (the ridge wind-tunnel data) and ERCOFTAC (the
RUSHIL data of Khurshudyan, Snyder and Nekrasov); we thank them. The case READMEs list the sources and licences.
