Forest canopy
=============

A forest slows the wind inside and above it. GUST represents a forest as a drag field: each cell holds the product
of a drag coefficient and the leaf area density of the trees in it. This page describes how that field is built
from the ``canopy`` inputs (:doc:`inputs`).

**Status.** The drag field is built, checked, written to the plotfile and acts on the momentum equations (`Drag on
the flow`_) and on the turbulence model (`Canopy turbulence`_).

The geometry follows ERF's ``ForestDrag`` (``Source/SourceTerms/ERF_ForestDrag.cpp`` in ERF): the same patch file,
the same two profile shapes and the same ``laimax`` range check (the copied ``ERF_ForestUtils.H``). The differences
are listed at the end of the page.

Leaf area density
-----------------

A canopy of height :math:`h` (m) and leaf area index :math:`\mathrm{LAI}` (leaf area per ground area) has the leaf
area density (m²/m³)

.. math::

   a(z) = \frac{\mathrm{LAI}}{h}\,\frac{P(z/h)}{N}, \qquad N = \int_0^1 P(\eta)\,d\eta, \qquad 0 \le z < h,

and :math:`a = 0` above :math:`h`. Here :math:`z` is the height above the local ground and :math:`\eta = z/h`. The
division by :math:`N` makes :math:`\int_0^h a\,dz = \mathrm{LAI}`. Two shapes :math:`P` are available:

- **Uniform** (``canopy.profile = uniform``, patch type 1): :math:`P = 1`, :math:`N = 1`.
- **Lalic–Mihailovic** (``canopy.profile = lalic_mihailovic``, patch type 2; Lalic and Mihailovic 2004): with
  :math:`\eta_m` = ``laimax`` the height of the largest density as a fraction of :math:`h`, and
  :math:`r = (1 - \eta_m)/(1 - \eta)`,

  .. math::

     P(\eta) = \begin{cases} r^6\,e^{6(1 - r)}, & \eta < \eta_m, \\ r^{1/2}\,e^{(1 - r)/2}, & \eta \ge \eta_m. \end{cases}

  :math:`P` rises to 1 at :math:`\eta_m` and falls to 0 at the canopy top.

The drag field of a cell
------------------------

The drag field is :math:`f = C_d\,a` (1/m), with :math:`C_d` the drag coefficient. Each cell takes the mean of
:math:`f` over its height above the column's ground, from :math:`z_b` to :math:`z_t`:

.. math::

   f = C_d\,\mathrm{LAI}\,\frac{F(z_t/h) - F(z_b/h)}{z_t - z_b}, \qquad F(\eta) = \frac{1}{N}\int_0^{\eta} P\,d\eta',

with :math:`F = 1` above the top. :math:`z_b`, :math:`z_t` and the ground are the means of the four corner nodes of
the cell's faces, the same means that give the cell volume factor :math:`J` (``detJ``). The column sum of
:math:`f\,J\,\Delta\zeta` is then :math:`C_d\,\mathrm{LAI}` to rounding, also when the canopy top lies inside a cell
and on any mesh spacing. A canopy taller than the domain loses the part above it; the log then warns.

**Why a cell mean.** On coarse layers a value taken at the cell centre can miss most of a thin canopy, or count a
cell as fully inside when its centre is. The mean keeps each column's total drag.

**Closed form.** :math:`F` is exact. With :math:`t = r(\eta)` and :math:`c = 1 - \eta_m`, :math:`d\eta = c\,dt/t^2`, and

.. math::

   \int_0^{\eta} P\,d\eta' = c\,[E(c) - E(t)], \quad
   E(t) = e^{6(1 - t)}\left(\frac{t^4}{6} + \frac{4t^3}{6^2} + \frac{12t^2}{6^3} + \frac{24t}{6^4} + \frac{24}{6^5}\right),
   \qquad \eta \le \eta_m,

.. math::

   \int_{\eta}^{1} P\,d\eta' = c\left[\frac{2}{\sqrt{t}}\,e^{(1 - t)/2} - \sqrt{2\pi e}\,\operatorname{erfc}\sqrt{t/2}\right],
   \qquad \eta \ge \eta_m,

and :math:`N` is the sum of the two at :math:`\eta = \eta_m`. The difference :math:`F(z_t/h) - F(z_b/h)` is formed
from these integrals directly, so a cell just below the top, where :math:`F` is within rounding of 1, keeps its
small share. Everything is evaluated in double precision.

Sources
-------

**Patches** (``canopy.type = patches``). ERF's forest file: one round patch per line,
``type x y height diameter cd lai laimax``. A column belongs to a patch when its centre is within half the diameter
of the patch centre. In a periodic direction the distance is taken to the nearest periodic image, so a patch across
the domain edge continues on the other side. Overlapping patches add their drag.

.. code-block:: text

   # type  x     y     height  diameter  cd   lai  laimax
   1       1024  512   20      200       0.2  6    0.8
   2       1024  1524  25      200       0.2  6    0.7

**Maps** (``canopy.type = map``). The canopy height, the leaf area index and, optionally, the drag coefficient as
``x y value`` point files, the format of ``terrain.file`` with ``terrain.type = xyz`` and of ``bc.surface_z0_file``.
Each is interpolated to the column centres by inverse-distance weighting with ``terrain.idw_neighbors`` and
``terrain.idw_power`` (:doc:`terrain_coordinates`). A height, an LAI or a :math:`C_d` of 0 means no canopy in
that column, so the maps may hold zeros where there is no forest. ``canopy.cd`` is the one value of every column and
must be > 0 (a run without a canopy sets ``canopy.type = none``). ``canopy.profile`` and ``canopy.laimax`` apply to
every column.

**Forest edges.** The interpolation blurs an edge over the points' spacing, so a 25 m forest next to open ground
gives intermediate heights there. The density in a column is the interpolated LAI over the interpolated height:
with the same points in both maps it stays near LAI/h across the edge, but with different points (or a large LAI
where the height is small) it can be several times the interior value. Give both maps on the same points.

.. code-block:: text

   canopy.type        = map
   canopy.height_file = canopy_height.xyz
   canopy.lai_file    = canopy_lai.xyz
   canopy.cd          = 0.2
   canopy.profile     = lalic_mihailovic
   canopy.laimax      = 0.8

Ghost cells and coarse levels
-----------------------------

A patch wider than a periodic direction of the domain would reach a column through two images and is refused.
Each cell checks every patch, so the cost grows with cells times patches (cells above the tallest patch skip the
check); it is built once per run.

The field has one ghost column in x and y. In a periodic direction it holds the periodic image; at other faces it
repeats the nearest domain column, as the roughness map does.

On the coarse levels of the grid hierarchy (:doc:`grid_hierarchy`) each layer takes

.. math::

   f_c = \frac{1}{r^2 J_c}\sum_{\text{fine cells}} f\,J,

the sum over the :math:`r \times r` fine cells under the coarse cell in the same layer. The integral of :math:`f`
over each coarse cell's volume then equals the integral over the fine cells it covers, so each coarse column keeps
the total drag. The coarse ground is smoothed, so :math:`J_c` can differ from the mean of the fine :math:`J`;
dividing by :math:`J_c` keeps the integral rather than the mean. The coarse canopy follows the layers, so its height
scales with :math:`J_c/J`: under a strongly smoothed ground it is spread over more or fewer layers than on level 0.
Only the coarse levels' solution, the first guess of level 0, sees this.

Drag on the flow
----------------

Each velocity face takes the momentum sink of ERF (``ERF_MakeMomSources.cpp``):

.. math::

   S_c = -f_c\,|U|\,u_c, \qquad |U| = \sqrt{u^2 + v^2 + w^2},

with :math:`u_c` the face's own component (:math:`c = x, y, z`), :math:`f_c` the mean of :math:`f` over the two cells
that share the face, and the other two components the means of their four faces around it. A :math:`z`-face on the
ground or the top of the domain takes the one cell inside; those rows are the boundary rows of :math:`w`, which the
boundary conditions set.

**Implicit diagonal.** In a dense canopy the drag rate :math:`f\,|U|` can exceed the pseudo-time rate
:math:`1/\Delta\tau` of the empty domain (:math:`C_d` = 0.2 and LAI 6 over 20 m at 2 m/s give about 0.1 /s). An explicit sink would
then limit the pseudo-time step. The momentum rows therefore gain a term on their diagonal
(``canopy.drag_diagonal``):

.. math::

   \texttt{derivative}:\ f_c\left(|U| + \frac{u_c^2}{|U|}\right) = \frac{\partial (f_c |U| u_c)}{\partial u_c},
   \qquad \texttt{picard}:\ f_c\,|U|, \qquad \texttt{none}:\ 0.

Each is :math:`\ge 0`, so :math:`1/\Delta\tau` plus the diagonal stays positive. The term reaches the line solves,
the SIMPLEC :math:`a_P` and its :math:`d = a_P - \sum a_{nb}`, so the pressure correction sees the drag, and the
Newton preconditioner. It changes the iterations, not the converged answer, and it does not enter the local
pseudo-time step, as with the damping layers. ``none`` keeps ERF's explicit treatment for comparison; an explicit
drag is stable only while :math:`f\,|U|\,\Delta\tau < 1` on the component along the wind (2 across it), which a
dense canopy breaks at a CFL of 1.

**Which diagonal.** On the drag alone (the scalar model :math:`du/d\tau = -f|u|u`), ``derivative`` damps every
pseudo-time step monotonically, by :math:`1/(1 + 2 f|U|\Delta\tau)`; ``picard`` damps by
:math:`(1 - x)/(1 + x)` with :math:`x = f|U|\Delta\tau`, which changes sign for :math:`x > 1` and tends to :math:`-1`
for large steps. In the coupled solver the convection and diffusion on the same diagonal bound both. On
``Tests/Canonical_RANS/Canopy_Drag`` the two took the same iterations within 5 % for LAI 1 and 3, and ``picard``
2 % and 10 % more for LAI 6 and 10. ``derivative`` is the default: it is the drag's own Jacobian diagonal, so the
SIMPLEC :math:`d` and the Newton stage's linearisation agree, and it does not depend on the other terms for its
damping.

**Linearisation state.** A solver computes the diagonal once per step, at the state it linearises about, before it
updates any component (``CanopyRows``): the SIMPLEC predictor at the current state, the Newton preconditioner at the
state of the Newton step, not at the increment it solves for.

**Face mean.** :math:`f_c` is the plain mean of the two cells, as in ERF. On terrain a mean weighted by the cell
volumes would keep the drag's integral exactly, as the coarse-level restriction does; the difference is of the
order of the change of :math:`J` between neighbouring cells.

The drag acts on the coarse levels of the grid hierarchy with their restricted field, and in the residual the Newton
stage differentiates. The column first guess (``init.first_guess = column``) carries it too, with the canopy's
:math:`k` terms, for a horizontally uniform canopy on flat terrain (:doc:`initial_conditions`, forest canopy).

Canopy turbulence
-----------------

The drag takes mean kinetic energy from the wind at the rate :math:`f|U|^3`. Part of it becomes turbulence in the
wakes of leaves and branches, at scales far below the mesh, and that turbulence dissipates fast (the short-circuit of
the energy cascade). The leaves also limit the size of the eddies. With ``canopy.tke = on`` (the default with a
canopy) the k–L closure takes ERF's canopy terms (erf-model/ERF#4340, which follows the KLAxell canopy of Kynema-SGF)
where :math:`f > 0`:

.. math::

   S_k \mathrel{+}= P_c - D_c,\qquad P_c = \beta_p f |U|^3,\qquad D_c = \beta_d f |U| k,
   \qquad L \le \alpha / f,

with :math:`|U|` the cell-centred wind, :math:`\beta_p` = ``canopy.beta_p`` the fraction of the drag work that
becomes :math:`k`, :math:`\beta_d` = ``canopy.beta_d`` the short-circuit coefficient (Green 1992, Liu et al. 1996,
Sanz 2003) and :math:`\alpha` = ``canopy.length_alpha`` (Harman and Finnigan 2007: :math:`\alpha = 2\beta^3` with
:math:`\beta = u_*/U_h`; 0.04 is :math:`\beta \approx 0.27`, the low end of their range). The defaults are ERF's
(1, 4, 0.04).

**Local balance.** Where the canopy terms and the dissipation dominate (no shear production or transport), the
capped dissipation :math:`\varepsilon = C_{\mu0}^3 k^{3/2} f/\alpha` is also proportional to :math:`f`, and the
balance :math:`P_c = D_c + \varepsilon` gives :math:`k = \kappa_c |U|^2` with

.. math::

   \beta_d\,\kappa_c + \frac{C_{\mu0}^3}{\alpha}\,\kappa_c^{3/2} = \beta_p,

independent of :math:`f` and :math:`|U|`: :math:`\kappa_c = 0.1728` for the defaults, so the cap and the sink act
together as an effective short-circuit coefficient :math:`\beta_p/\kappa_c \approx 5.8` (the uncapped k–ε values
are 4 to 5.1). Without the cap (``canopy.length_alpha = 0``) and where :math:`\varepsilon` is small,
:math:`\kappa_c \to \beta_p/\beta_d`. Near the canopy top the shear layer above adds production and transport, so
:math:`k` there is larger than :math:`\kappa_c |U|^2`.

**The length cap** acts on the turbulent length before the dissipation :math:`\varepsilon = C_{\mu0}^3
k^{3/2}/L`, the Richardson number, the stability functions and :math:`\nu_t` are formed from it (ERF's
``limited_length``, the copied ``ERF_CanopyTurbulence.H``). It acts only where :math:`f > 0` and :math:`\alpha < L
f`, so near the ground, where :math:`L = \kappa (d + z_0)` is already short, the closure is unchanged. Whether a
cell is capped is one of the closure's branch decisions: a full freeze of the decisions (a Jacobian by differences)
holds it with the others, while the Newton stage (``turb.frozen_switch``) holds only the :math:`k`-floor decisions
and takes the cap, like the stratification branches, from the state. A checkpoint of fully frozen decisions records
the cap's :math:`\alpha`, and a restart with another (``canopy.tke = off`` included) is refused.

**Implicit side.** The sink :math:`D_c` is linear in :math:`k`, so the diagonal of the :math:`k` rows gains
:math:`\beta_d f |U|` exactly (the point-implicit Patankar term, :math:`D_c/k`) and the exact source derivative the
Newton stage reads loses it. The capped length raises :math:`\varepsilon`, which is already on that diagonal as
:math:`\varepsilon/k`. The production :math:`P_c` does not depend on :math:`k` and stays an explicit source; its
dependence on the wind is lagged, as the shear production's. So neither term limits the pseudo-time step, and
:math:`k` stays positive.

**Steady form.** ERF integrates :math:`P_c` and :math:`D_c` together and exactly over each time stage (an
exponential update), so that an explicit step never removes more :math:`k` than there is, and it needs
``erf.implicit_tke_dissipation`` with the cap. GUST's diagonal does both jobs in pseudo-time: the point update
:math:`k^+ = (k + \Delta\tau P)/(1 + \Delta\tau\,a)`, with :math:`a` the diagonal, stays positive for any
:math:`\Delta\tau`. GUST's fixed point is ERF's in the limit :math:`\Delta t \to 0`; at a finite step ERF's
exponential form weights the canopy pair by :math:`\varphi_1(\beta_d f |U| \Delta t)` against the other terms, so
the two agree exactly only where the canopy pair is alone.

**Ground roughness.** Under the canopy the ground's :math:`z_0` should be that of the soil and litter (0.01 to
0.1 m). A forest-sized :math:`z_0` (about :math:`h/10`) counts the forest twice; the run warns when a column with
canopy has a ground :math:`z_0` of 0.5 m or more.

``canopy.tke = off`` keeps the drag only, the behaviour of ERF before #4340 and of GUST before this term (#75). With
the default ``on``, a run with a canopy gives other answers than it did with #75. These inputs are not part of the
checkpoint's canopy hash: like the ``turb.*`` inputs, they may change across a restart, which then converges to the
new inputs' answer (except the frozen-decision case above).

Canopy heat
-----------

With ``canopy.heat = fixed`` the leaves exchange sensible heat with the air at a fixed leaf temperature. This is the
steady form of ERF's fixed-leaf-temperature canopy biophysics (``Source/SourceTerms/ERF_CanopyBiophysics.H`` in ERF,
erf-model/ERF#3897, from the ICLASS canopy model). The :math:`\theta` equation gains

.. math::

   S_\theta = \lambda\,(\theta_l - \theta),\qquad
   \lambda = \frac{2a}{r_b},\qquad
   r_b = C_b\sqrt{\frac{l_w}{\max(|U|, U_{\min})}}.

- :math:`a`: the leaf area density (1/m), the cell means of the drag field with :math:`C_d = 1` (ERF's frontal area
  density). A patch with :math:`C_d = 0` has no drag but still has leaves.
- :math:`\theta_l`: ``canopy.leaf_theta``, the leaves' potential temperature (K).
- :math:`|U|`: the cell-centred wind, from the means of the cell's faces.
- :math:`r_b`: the leaf boundary-layer resistance (s/m), with ERF's constants :math:`C_b` = 100 s\ :sup:`1/2`/m,
  leaf width :math:`l_w` = 0.05 m and :math:`U_{\min}` = 0.01 m/s. The factor 2 counts both sides of a leaf.

ERF adds :math:`2\rho c_p (T_l - T)\,a / (r_b\,c_p\,\Pi)` to :math:`\rho\theta`, with :math:`T = \Pi\theta`; the
density, :math:`c_p` and the Exner function :math:`\Pi` cancel. ERF's transpiration term (with moisture) is not here:
GUST carries no water vapour.

**Implicit side.** :math:`\lambda` does not depend on :math:`\theta`, so the :math:`\theta` rows' diagonal gains
:math:`\lambda`, the exact derivative. Its dependence on the wind is lagged, as the drag's is. :math:`\lambda \ge 0`,
so the diagonal stays positive, and the source pulls :math:`\theta` towards :math:`\theta_l`. The coarse levels of
grid sequencing restrict :math:`a` as they restrict :math:`f`, keeping each layer's integral.

**With** ``simplec.theta = frozen``. :math:`\theta` is then data, with no equation for the source to act in
(:doc:`solution_algorithm`). ``canopy.heat_frozen`` chooses what happens:

- ``skip`` (default): the leaf heat is left out; the run warns at start-up.
- ``blend``: in the canopy :math:`\theta` is the steady point balance of the leaf exchange and of turbulent mixing
  towards the frozen profile :math:`\theta_f`, the mixing taken as a relaxation at the eddy turnover rate
  :math:`\varepsilon/k`:

  .. math::

     0 = \lambda(\theta_l - \theta) + \frac{\varepsilon}{k}(\theta_f - \theta)
     \quad\Rightarrow\quad
     \theta = \theta_f + w\,(\theta_l - \theta_f),\qquad
     w = \frac{\lambda k}{\varepsilon + \lambda k} = \frac{\mathrm{Da}}{1 + \mathrm{Da}},\qquad
     \mathrm{Da} = \frac{\lambda k}{\varepsilon}.

  :math:`k` (at ``turb.tke_floor`` where below it, as :math:`\varepsilon` takes it) and :math:`\varepsilon` are the
  closure's at the iteration's state; :math:`\varepsilon` does not include the canopy's short-circuit sink, which
  is the wake turbulence's cascade, not the renewal of the canopy air. The blend is set once per SIMPLEC iteration,
  after the pressure correction. :math:`w = 0` where :math:`a = 0`, so :math:`\theta` outside the canopy
  is :math:`\theta_f` bit for bit; :math:`0 \le w < 1`, and :math:`w` grows smoothly with the foliage and falls with
  the turbulence. There is no new input, but the rate :math:`\varepsilon/k` is a modelling choice: in a canopy it
  lies between the flushing rate :math:`u_*/h` (about 4 times slower) and :math:`\nu_t/L^2` (about 3 times
  faster), and close to the inverse Lagrangian time scale :math:`u_*/(0.3h)` (Raupach 1989), so :math:`\mathrm{Da}`
  is uncertain by about a factor of 3. With leaves cooler than the frozen air the loop
  :math:`\theta \to N^2 \to k \to w` is monotone (a cooler canopy damps :math:`k`, which raises :math:`w`); with
  thin cells and weak shear at the canopy top it can have a second fixed point with :math:`k` near its floor. The
  blended :math:`\theta` enters the closure's :math:`N^2` like the
  frozen field; buoyancy stays out of the momentum equations. The log's :math:`\theta` residual stays 0, as
  with ``frozen`` alone. A uniform profile has a jump in :math:`a` at the canopy top, so :math:`\theta`
  can jump there too; the Lalic–Mihailovic profile tapers to the top.

**Leaves warmer than the air, θ solved.** On ``Tests/Canonical_RANS/Forested_Escarpment`` (neutral 300 K air, the
LAI 5 forest, 2 ranks, default CFL cap 8) the iterations to the verdict are:

.. list-table::
   :header-rows: 1

   * - Leaves
     - 298 K
     - 300 K
     - 300.25 K
     - 300.5 K
     - 301 K
     - 301.5 K
     - 302 K
   * - Iterations
     - 201
     - 249
     - 240
     - 239
     - 225
     - 204
     - not steady (600; and 1500 at cap 2)

Without the heat the case takes 256. At 302 K :math:`r` cycles between 0.3 and 0.4, and a lower CFL cap does not
change that: the onset lies between 1.5 and 2 K above the air. The heated canopy air is unstably stratified, with
buoyancy in the momentum equations; in the top canopy cells the closure's :math:`R_t` reaches its lower smoothing
bound, so the closure cannot mix more, and a steady state, if one exists, is not one pseudo-time can reach. For warm leaves use ``simplec.theta = frozen`` with
``canopy.heat_frozen = blend``, which converged in 212 iterations with leaves at 304 K on the stable deck of the same
case.

**Costs and limits.** One cell kernel per residual and one per linearisation. The plotfile carries the leaf area
density as ``canopy_lad`` when ``canopy.heat`` is on. The column first guess
(``init.first_guess = column``) does not carry the leaf heat. The checkpoint's canopy hash does not include these
inputs: like ``canopy.tke``, they may change across a restart.

Steep terrain and the boundaries
--------------------------------

**CFL cap and the step floor.** Before the local step floor (``steady.dtau_floor``, :doc:`solution_algorithm`) a
dense canopy needed a lower ``steady.cfl_max`` than the empty domain: on ``Tests/Canonical_RANS/Forested_Escarpment``
(a 200 m tanh escarpment, steepest slope 22°, with the WINSENT canopy: uniform, 25 m, LAI 5, :math:`C_d` = 0.15, on
the slope and the plateau's edge; 8 m first layer stretched by 1.12; LU-SGS with grid sequencing over three levels, 2
ranks) the iteration cycled at caps of 4 and above and was marginal at 3 (#77). The cause is not the terrain: the same
forest on flat ground and on an 11° ramp cycles the same way at the default cap. It is a cell at rest. Under the
forest, 250 m inside its windward edge and again at its lee edge, the wind in the lowest layers slows to a few
centimetres per second and reverses. There the convective and diffusive rates of the local pseudo-time step vanish, so
:math:`\Delta\tau` is unbounded (:math:`\mathrm{CFL}/10^{-12}` s), and the drag's own diagonal
:math:`f(|U| + u_c^2/|U|)` vanishes with the wind. The momentum rows' :math:`a_P - \sum a_{nb}` then tends to zero,
the SIMPLEC :math:`d = a_P - \sum a_{nb}` and with it the bound on the pressure correction's velocity update in
the cell vanish, and the drag reacts to the overshoot at the next step: the :math:`w` and :math:`p` residuals lead a
cycle of about 40 iterations between :math:`r` 0.1 and 0.5. This is the observation; a 2-D SIMPLEC model with the
same terms does not reproduce the cycle, so what closes the loop in the 3-D solver is not identified
(:doc:`solution_algorithm`).

The step floor removes it: a cell at rest now steps at most as the level's fastest cell does, so :math:`a_P` keeps
at least the free stream's :math:`1/\Delta\tau`. From the inflow profile to the verdict (:math:`r` = 1e-4), 2 ranks:

.. list-table::
   :header-rows: 1
   :widths: 30 20 18 14 18

   * - Case, ``steady.cfl_max`` = 8
     - ``dtau_floor`` = 0 (the old step)
     - **``dtau_floor`` = 0.1 (default)**
     - ``dtau_floor`` = 0.3 / 1
     - floor 0, ``cfl_max`` = 2
   * - Forested escarpment (22°)
     - not steady after 3000 (a cycle)
     - **256**
     - 165 / 101
     - 185
   * - The same forest on an 11° ramp
     - not steady after 500
     - **192**
     - 161 / 108
     - —
   * - The same forest on flat ground
     - not steady after 500 (bursts to :math:`r` = 35)
     - **182**
     - 124 / 109
     - —
   * - The escarpment without the forest
     - 82
     - **82**
     - 84 / 70
     - —
   * - Leafless (LAI 0.7), escarpment
     - 109
     - —
     - —
     - —

Measured on the way (each on the flat and the escarpment cases at the default cap, the floor off): a fixed CFL of 3.5
or 4 cycles, so the SER controller is not the cause; a floor on the drag diagonal's speed, :math:`\max(|U|, U_{\min})`
with :math:`U_{\min}` from 0.25 to 4 m/s, also converges (280 to 540 iterations, with bursts; it bounds the same
:math:`d`, in the canopy cells only, and the fixed point is unchanged to :math:`10^{-11}`); a floor at the face's
balance speed :math:`\sqrt{|G|/f}` from the residual does not (the residual is small in the cell at rest); a bound
:math:`d \le 1/(\epsilon a_P)` does not (:math:`a_P` itself is small there); the rate floor at 0.03 to 1 times the
free stream's converges in every case and at 0.3 to 1 fastest, at 0.01 not. The floor is the general form and takes
the canopy decks back to the default cap; its default 0.1 is the largest value that leaves the empty cases' paths
close to what they were (:doc:`solution_algorithm`). ``steady.dtau_floor = 0`` restores the old step, and
``Tests/Canonical_RANS/Forested_Escarpment`` then needs ``steady.cfl_max`` 2 to 2.5 as before.

**Newton stage.** On ``Tests/Canonical_RANS/Canopy_Drag`` (the two patches, ``canopy.tke = off``, the verdict at
1e-8) ``steady.stage = simplec_newton`` switches to Newton once :math:`r` has been below 1e-4 for 20 iterations and
hands back after failed steps, with the canopy and without it (#75). Measured here, two ranks: the switch at
iteration 1548 (with the step floor at 1; 1859 without); the first Newton step converges its GMRES in 90 iterations with
the SIMPLEC preconditioner (the drag rows at the step's state) but the line search accepts only :math:`\lambda` =
1/16 and the residual falls by 1 %; the second step fails after five attempts, the continuation then halves the CFL
each step and the stage hands back; the run is steady after 4180 iterations (7955 without the floor, 5710 with
SIMPLEC alone). Without the floor GMRES also stalls at its 2000-iteration limit once the continuation's CFL passes
16. The same happens with the geostrophic forcing replaced by a pressure gradient (4950 iterations), so the cause is
neither the drag rows, the rotation nor the preconditioner's convergence at the first step: the Newton direction
itself is poor on this case. Which term makes it so (the closure's floor and realizability switches in the dense
patch, where :math:`k` sits at its floor, or the wall law's 0.1 m/s speed floor under the trees) is not found yet.
Without the canopy SIMPLEC alone does not reach 1e-8 on this periodic geostrophic deck in 10000 iterations (an
inertial adjustment), as noted in #75.

**Inflow faces and damping layers.** The inflow profile and the damping layers' targets do not feel the drag. The run
warns, with the count, when canopy cells lie in the first cell column inside an inflow face, or inside the sponge or
the Rayleigh layer (at the cell centre). The inflow check covers the first column only; keep the forest a few cells
away from inflow faces, as a terrain blend that fades the forest out does.

Output and restart
------------------

The run log gives one line for level 0, and with grid sequencing one for each coarse level as it starts: the
largest :math:`f`, the number of columns with canopy, and the mean and largest column integral of :math:`f`, which
is :math:`C_d\,\mathrm{LAI}`. A canopy that covers no column, for example a file in another coordinate frame, is
flagged with a warning. The plotfile carries level 0's :math:`f` as ``canopy_drag`` (:doc:`output`); the coarse
levels' fields are not written.

In a single-precision build the heights above the ground come from single-precision node heights: at 3000 m above
sea level they are rounded to about 2e-4 m, so a 0.1 m layer's :math:`f` is within about 0.2 %. The column integral
is still exact, since the layers telescope.

The checkpoint holds the SHA-256 of the canopy inputs: the options and the contents of the files read, byte for
byte. A restart with other canopy inputs, with a canopy where the checkpoint has none, or without one where it has
one, stops and names ``io.restart``; it is checked before the field is built. Renaming a file with the same contents
is accepted; editing a comment in it counts as a change. ``canopy.drag_diagonal`` is not part of the hash: it changes
the iterations, not the fixed point, so it may change across a restart.

Differences from ERF
--------------------

- ERF takes :math:`a(z)` at the cell centre; GUST takes its mean over the cell. On fine layers the two agree. On
  layers of :math:`h/4` or thicker the centre value's column integral departs from :math:`C_d\,\mathrm{LAI}`: for
  :math:`h` = 20 m and :math:`\eta_m` = 0.8 it is 1.11 LAI with 5 m layers, 1.43 with 10 m and 0.15 with 15 m. The
  mean gives LAI exactly.
- ERF normalises the Lalic–Mihailovic profile with a 100-point sum; GUST uses the closed form. ERF's :math:`N` is
  within 0.8 % of the exact value for :math:`\eta_m \le 0.97` (0.73 % at 0, 0.01 % at the default 0.8); it departs
  as :math:`\eta_m` approaches 1 (−10 % at 0.99).
- ERF reads gridded maps from NetCDF on a uniform grid and interpolates them bilinearly; GUST reads point files and
  uses inverse-distance weighting, as for the terrain and the roughness map, and needs no NetCDF.
- In a periodic direction GUST uses the nearest periodic image of a patch; ERF uses the plain distance.
- The leaf heat (``canopy.heat``) uses the cell mean of the leaf area density, as the drag does, so drag and heat
  see the same leaves; ERF uses the centre value of its frontal area density. ERF's transpiration term (with
  moisture) is not ported. The constants of the leaf resistance are ERF's (:math:`C_b` = 100): its
  :math:`r_b` is about 3 times below a laminar flat plate's per side, at the high-conductance end of field values.
