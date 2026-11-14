Boundary conditions
===================

The computational domain is a box with six faces: west and east (``xlo``, ``xhi``), south and north (``ylo``,
``yhi``), the ground (``zlo``) and the top (``zhi``). Each face gets one boundary type. The typical ABL setup is
periodic sides, a rough ground described by a wall model, and a slip top. The boundary types are ERF's. The code is
in ``Source/BoundaryConditions`` (the fills and Dirichlet rows), ``Source/BoundaryInputs`` (the ``bc`` reader and its
checks), ``Source/BoundaryDamping`` (the damping layers), ``Source/InflowProfile`` and ``Source/SurfaceLayer``; the
inputs are in the ``bc``
namespace (:doc:`inputs`).

Boundary types
--------------

.. list-table::
   :header-rows: 1
   :widths: 18 14 68

   * - Type
     - Faces
     - Condition
   * - ``periodic``
     - x, y
     - The flow leaving one side enters the opposite side. Used where ``mesh.is_periodic = 1``.
   * - ``inflow``
     - x, y
     - :math:`u, v, \theta, k` from the inflow profile at their height above the ground.
   * - ``outflow``
     - x, y
     - Zero normal gradient for every field; :math:`p = 0`.
   * - ``nest``
     - x, y
     - The parent run's field on the face of a nested patch (`Nest faces`_).
   * - ``slip_wall``
     - all
     - No flow through the face; no stress along it. At the top, :math:`\theta` may have a given gradient
       (``bc.zhi_theta_grad``).
   * - ``no_slip_wall``
     - all
     - The velocity equals the wall velocity (0, or a moving lid set by ``bc.<face>_wall_velocity``).
   * - ``surface_layer``
     - ground
     - Monin–Obukhov wall model (below). The default ground.
   * - ``geostrophic``
     - top
     - :math:`u, v` held at the geostrophic wind.

The defaults are periodic sides, a ``surface_layer`` ground and a ``slip_wall`` top, as in ERF's canonical ABL
cases.

**No flow through the ground and the top.** On a terrain-following mesh the velocity normal to the ground is the
contravariant velocity :math:`\Omega` (:doc:`terrain_coordinates`), not :math:`w`. The ground and top values of
:math:`w` are therefore

.. math::

   w_\mathrm{face} = z_x\,\bar u + z_y\,\bar v ,

which gives :math:`\Omega = 0`. Over flat ground :math:`w_\mathrm{face} = 0`.

**Outflow faces: the transport fluxes of the ghost cells.** The ``outflow`` fill copies :math:`u, v, w, \theta, k`
along the grid line into the ghost cells outside the face. The ghost node columns take the face node's height, so a
ghost cell is flat across the face (along the face it keeps the face line's slope). Convection uses the transport
fluxes :math:`F_x = a_x u`, :math:`F_y = a_y v` and :math:`\Omega = w - z_x \bar u - z_y \bar v` (MassFlux.H), and
the control volume of the face-point velocity reaches half a cell into the ghost cell, so it reads the ghost cell's
fluxes. By default (``bc.outflow_transport = mirror``) the ghost cell is given the edge cell's fluxes as a mirror
image through the face:

.. math::

   \Omega_{e+m} = \Omega_{e-m+1}, \qquad F_{t,\,e+m} = F_{t,\,e-m+1}, \qquad
   F_{n,\,f+m} = 2\,F_{n,\,f} - F_{n,\,f-m}, \qquad m = 1, 2, \dots

- :math:`e`: the last cell inside the face; :math:`e + m`: the :math:`m`-th ghost cell outside it (counted outwards,
  so towards lower indices on a low face), at the same level and position along the face; :math:`e - m + 1`: its
  mirror image inside.
- :math:`F_t`: the tangential flux (:math:`F_y` on an x face, :math:`F_x` on a y face), at cell positions along the
  normal; :math:`F_n`: the normal flux, on the faces :math:`f` (the boundary face) and :math:`f \pm m`.
- :math:`a_x, a_y`: the area factors of the x- and y-faces (Mesh.H); :math:`z_x, z_y`: the slopes of a z-face, from
  its four nodes;
  :math:`\bar u, \bar v`: the means of the four :math:`u` (:math:`v`) points around it.

Each ghost cell then has its edge cell's divergence (in index space), whatever the ground does along or across the
face, and the face-point control volume, half in each, has the edge cell's net flux: for a mass-consistent field the
face-point rows are like interior rows. The velocity ghosts are not changed, so the
strain rate, the diffusion and the closure at the face see the same values as before.

With ``bc.outflow_transport = from_velocity`` (the only choice up to v1.0.0) the fluxes are computed from the copied
velocities and the ghost cell's own geometry. These terms then differ from the edge cell's:

- :math:`\Omega_g - \Omega_e \approx z_{x,e}\,\bar u_e` on an x face (plus the difference of the along-face slopes
  times :math:`\bar v`): a vertical flux where the face cuts a slope.
- :math:`F_{n,f+1} \approx F_{n,f-1}` (the copied normal velocity), so the ghost cell's normal divergence is minus the
  edge cell's.
- In 3-D the tangential flux, through the flat ghost cell's area factors (a smaller term).

The face-point control volume then has a net flux the implicit diagonal (:math:`a_P = \sum a_\mathrm{nb}`) does not
hold. Where it is an outflow (a face cutting a rising slope) the true diagonal is larger, the face-point rows
overshoot and the iteration alternates between two states (``Tests/Canonical_RANS/Outflow_Steep_Terrain``: relative
residual 0.176 with ``from_velocity``, steady in 148 iterations with ``mirror``). Mirroring :math:`\Omega` alone
removes that flux but leaves the normal-flux term: on ``Aspect_Ratio_Hill_2D`` (dx/dz1 = 100, a slope of 5e-4 at the
face) the spanwise velocity residual then decayed slowly at the face column (546 iterations, against 181 with
``from_velocity`` and 161 with the full mirror). The cause of that slow mode was not isolated; the full mirror,
which makes the face-point control volume consistent, does not show it.

The two choices differ over flat ground too (the normal flux), by little: on ``Askervein_210`` the 10 m speeds change
by at most 7e-6 (relative). With backflow through the face (or an edge cell that slows the flow by more than two
thirds) the extrapolated normal flux turns the outer half of the face-point control volume into an inflow; the rows
stay well posed (upwind weights), but a zero-gradient outflow is not meant for backflow with either choice. A restart from a checkpoint written with v1.0.0 continues with ``mirror``; give
``bc.outflow_transport = from_velocity`` to continue the old run unchanged.

Wall model
----------

The first cell above the ground is usually much thicker than the thin layer where the wind drops to zero. A wall
model bridges this layer: it gives the stress at the ground from the wind in the first cell, using Monin–Obukhov
similarity theory (MOST). In each column,

.. math::

   u_* = \frac{\kappa\,U}{\max\left(\ln\dfrac{d + z_0}{z_0} - \psi_m(d/L_O),\; m\right)}, \qquad
   \frac{1}{L_O} = -\frac{\kappa\,g\,q}{u_*^3\,\theta_0} .

- :math:`u_*`: friction velocity.
- :math:`U = \max(|\mathbf U_t|, 0.1~\mathrm{m/s})`: wind speed parallel to the ground at the reference height
  (``bc.surface_zref``; by default the first cell centre).
  :math:`\mathbf U_t = \mathbf U - (\mathbf U\cdot\mathbf n)\,\mathbf n` is the cell-centre wind with its
  component along the ground normal :math:`\mathbf n` removed, interpolated in :math:`\ln(d + z_0)` between two
  centres. In the cell on the ground the centre wind is taken parallel to the ground,
  :math:`(u_c,\, v_c,\, z_x u_c + z_y v_c)`: the no-flux rule (:math:`\Omega = 0`) applied to the cell's own centre
  :math:`u_c, v_c`, which is exact for a flow along a uniform slope; the cell's :math:`w` faces are not read. The ghost
  fill writes the ground :math:`w` face and the :math:`u, v` ghosts after the wall model, so the wall model reads no
  value the fill writes, and a second fill changes nothing. Over flat ground the vertical part is projected out
  either way. (Before audit A2 the vertical part was the mean of the ground and top :math:`w` faces; the ground face
  is written by the fill from the wall law's ghosts, so the fill was not a fixed point.) A reference height between
  the last two cell centres reads the top :math:`w` face, which the fill writes too; over a top that is not flat the
  fill is then a fixed point only to that face's change.
- :math:`d = h/|A|`: distance from the ground along its normal, with :math:`|A| = \sqrt{1 + z_x^2 + z_y^2}`.
- :math:`z_0`: roughness length (``bc.surface_z0``, or per column from ``bc.surface_z0_file``,
  :ref:`roughness map <roughness-map>`).
- :math:`q`: kinematic surface heat flux (``bc.surface_heat_flux``).
- :math:`L_O`: Obukhov length.
- :math:`\psi_m`: ERF's Businger–Dyer stability function, or QNSE's on the stable side (below).
- :math:`m`: a floor on the denominator (below).

With ``bc.surface_stability = heat_flux``, :math:`u_*` and :math:`L_O` are solved together by Newton's method to
round-off. (ERF uses a relaxed fixed point to 1e-5; a converged :math:`u_*` keeps the residual smooth, which a
finite-difference Jacobian needs.)

**Denominator floor.** Without a floor, strong convection at low wind gives :math:`u_* < 0`.
``bc.surface_denominator_floor`` selects the floor:

- ``log_law`` (default): :math:`m = \min(\ln((d + z_0)/z_0),\, 1)`. In neutral flow the floor never acts, so
  :math:`u_*` is the log-law value at every height.
- ``one``: :math:`m = 1`. In neutral flow this lowers :math:`u_*` when the reference point is below
  :math:`(e - 1)\,z_0 \approx 1.72\,z_0`, for example a 1 m first layer over :math:`z_0 = 0.5` m.

The two give identical results where :math:`\ln((d + z_0)/z_0) \ge 1`.

.. _stable-functions:

**Stable functions.** ``bc.surface_stable_functions`` selects the stable side
(:math:`\zeta = z/L > 0`) of the flux-gradient functions :math:`\phi_m, \phi_h` and their integrals
:math:`\psi = \int_0^\zeta (1 - \phi(x))/x\,dx`. Both use one form,

.. math::

   \phi(\zeta) = \frac{1 + a\,\zeta}{1 + b\,\zeta}, \qquad
   \psi(\zeta) = -\frac{a - b}{b}\,\ln(1 + b\,\zeta) .

- ``businger_dyer`` (default): :math:`a = 5`, :math:`b = 0`, so :math:`\phi = 1 + 5\zeta` and
  :math:`\psi = -5\zeta` (Dyer 1974, ERF's functions).
- ``qnse``: the quasi-normal scale elimination functions (Sukoriansky, Galperin and Perov 2005), as ERF's
  ``enable_qnse_stable_functions``: :math:`a_m, b_m` (``bc.surface_qnse_am``, ``_bm``) for momentum and
  :math:`a_h, b_h` (``_ah``, ``_bh``) for heat, ERF's defaults 2.5 and 0.2. :math:`\phi` tends to
  :math:`a/b = 12.5` in very stable air, so the surface stress does not collapse as fast as with
  Businger–Dyer's linear :math:`\phi`.

The unstable side is Businger–Dyer's with either choice. The selection applies to the wall model's
:math:`u_*` (``bc.surface_stability = heat_flux``), the ground gradients of the ghost cells, and the stratified
column (:ref:`stratified column <column-stability>`).

The wall model replaces the diffusive fluxes through the ground:

.. list-table::
   :header-rows: 1
   :widths: 18 82

   * - Field
     - Ground flux
   * - :math:`u, v`
     - :math:`u_*^2\,|A|\,\mathbf U_t/U`, along the ground wind
   * - :math:`\theta`
     - :math:`-q\,|A|`
   * - :math:`k`
     - 0; the first cell relaxes to :math:`k_w = (u_*^3 + \kappa B d_0)^{2/3}/C_{\mu 0}^2` (AL01 Eq. 16),
       with :math:`d_0` the first-cell-centre distance and :math:`B = gq/\theta_0` for :math:`q > 0`. With
       ``bc.surface_wall_k = dirichlet`` (default) this replaces the first cell's :math:`k` equation; with
       ``source`` it is added to it.

The ghost cells below the ground continue the MOST profile.

.. _roughness-map:

**Roughness map.** ``bc.surface_z0_file`` replaces the single :math:`z_0` with one value per ground column, for a
surface such as water beside land. The file holds ``x y z0`` points (the format of ``terrain.file``). Each ground
column takes the inverse-distance-weighted mean of :math:`\ln z_0` over the ``terrain.idw_neighbors`` nearest
points (weight :math:`r^{-p}`, :math:`p` = ``terrain.idw_power``) at its cell centre:

.. math::

   z_{0,ij} = z_{0,\mathrm{ref}}\,\exp\!\left(\frac{\sum_n w_n \ln(z_{0,n}/z_{0,\mathrm{ref}})}{\sum_n w_n}\right),
   \qquad w_n = r_n^{-p},

- :math:`z_{0,n}`: the file's value at point :math:`n`, at distance :math:`r_n` from the cell centre;
- :math:`z_{0,\mathrm{ref}}`: the file's largest value (a file of one value then gives that value exactly).

The log law depends on :math:`z_0` through :math:`\ln z_0`, so the mean is taken on that scale; an arithmetic mean
would let a few rough points dominate a smooth neighbour. A point within :math:`10^{-6}` m of the centre gives its own
value. Ghost columns take the column they image (periodic image, or the nearest domain column). Every reader of
:math:`z_0` on the ground uses the column's value: the wall model above, the ghost fill's wall-law gradient, and the
turbulence model's geometric length :math:`\kappa(d + z_0)`. The inflow profile keeps ``bc.inflow_z0``; where it
differs from the map next to the inflow face, an internal boundary layer grows from that face (over water
:math:`z_0 = 0.0003` m below an inflow of :math:`z_0 = 0.01` m, about 9 m deep after 200 m by Elliott's 1958 estimate),
so the inflow should take the upstream surface's :math:`z_0` and a matching :math:`u_*`.

On the coarse levels of grid sequencing (``mesh.hierarchy_levels`` > 1) a coarse column covering
:math:`r \times r` level-0 columns takes their geometric mean,
:math:`z_{0,c} = z_{0,\mathrm{ref}} \exp\big(\tfrac{1}{r^2}\sum \ln(z_{0}/z_{0,\mathrm{ref}})\big)`. This is a
simple aggregation on the log-law scale, not a blending-height drag average (Mason 1988): in a column half water
and half land it gives about 20–30 % less drag at a 0.5–1 m first centre than the area mean of the two drags. The
coarse levels only seed level 0, which uses the map itself. The run prints the minimum, maximum and geometric mean of each level's
map, and the turbulence plotfile (``io.turb_plot_file``) holds it as ``z0`` on the first cells. A restart reads the
file the deck names again: the checkpoint does not record the map. ``init.first_guess = column`` needs a single
:math:`z_0` and refuses a map.

``Tests/Real_Terrain/WindNinja2019/bolund_z0_map.py`` writes the Bolund map from the Bolund terrain points: 0.0003 m over water (terrain at or
below the 0.75 m water level) and 0.015 m over land, with land roughness east of :math:`x = 327` m, the values and
rules of the Bolund blind comparison (Bechmann et al. 2011).

Inflow profile
--------------

The inflow profile gives the wind, :math:`\theta` and :math:`k` as functions of the height :math:`z` above the
ground. It is used on inflow faces, in the damping layers and for the initial field (:doc:`initial_conditions`). The
wind speed :math:`S(z)` points along ``bc.inflow_direction``:

.. math::

   \text{log\_law:}\quad S = \frac{u_*}{\kappa}\ln\frac{z + z_0}{z_0}, \qquad
   \text{power\_law:}\quad S = S_\mathrm{ref}\left(\frac{z}{z_\mathrm{ref}}\right)^\alpha .

A ``file`` profile tabulates ``z u v theta [k]``. For the analytic profiles :math:`\theta` is the sounding
``bc.theta_profile``.

**Turbulence of the profile.** Without a :math:`k` column, :math:`k` and :math:`\nu_t` are in equilibrium with the
wind's shear: shear production equals dissipation in the neutral k–L closure (:doc:`turbulence_model`),

.. math::

   u_{*l} = l\,\frac{dS}{dz}, \qquad k = \frac{u_{*l}^2}{C_{\mu 0}^2}, \qquad \nu_t = l\,u_{*l}, \qquad
   \frac{1}{l} = \frac{1}{\kappa (z + z_0)} + \frac{1}{l_\mathrm{max}} .

- :math:`u_{*l}`: the local friction velocity; :math:`l`: the length; :math:`l_\mathrm{max}`: ``turb.max_geom_lscale``.
- ``bc.inflow_tke_length = closure`` (the default) uses this :math:`l`, the closure's own geometric length, so the
  inflow :math:`k` is the closure's equilibrium. ``unbounded`` drops the cap (:math:`l = \kappa(z + z_0)`); the log
  law then has :math:`u_{*l} = u_*` and :math:`k = u_*^2/C_{\mu 0}^2` at every height.
- **Why the cap.** Without it, :math:`\kappa z\,dS/dz` grows with height on a power law, and so does :math:`k`: on
  ERF's Askervein inflow it reaches about 9 m²/s² at 316 m, where ERF's flow has 0.31; with the cap, 0.31 at 316 m.
- A :math:`k` column is used as given, and ``bc.inflow_tke_length`` is refused (it would change nothing).
- The equilibrium is the neutral one, also above an inversion, and the height is the height above the ground (the
  closure's wall distance with ``turb.wall_distance = terrain_height``). The roughness is ``bc.inflow_z0``, which
  should equal ``bc.surface_z0``, or with ``bc.surface_z0_file`` the map's :math:`z_0` next to the inflow face.
- A log-law wind with the capped length is not a momentum equilibrium (its stress falls with height as
  :math:`(l/\kappa(z + z_0))^2`); with the unbounded length the :math:`k` equation was out of balance instead. The
  domain adjusts either way; the capped form at least balances the :math:`k` equation.

Damping layers
--------------

Damping layers absorb waves before they reflect from the top or the sides. Both are off by default. They relax a
field :math:`q` towards the inflow profile :math:`q_\mathrm{ref}` (0 for :math:`w`) by adding a term to the
residual:

.. math::

   \text{Rayleigh (top):}\quad R \mathrel{+}= c_r \sin^2\!\left(\frac{\pi}{2}\,\frac{z - z_d}{z_\mathrm{top} - z_d}\right)
   (q - q_\mathrm{ref}),
   \qquad
   \text{sponge (sides):}\quad R \mathrel{+}= c_s\left(1 - \frac{s}{w_s}\right)^2 (q - q_\mathrm{ref}) .

- :math:`z_d`: start height of the Rayleigh layer (``bc.rayleigh_zdamp``); :math:`z_\mathrm{top}`: the top level of the
  mesh, which is not ``mesh.prob_hi`` z on a stretched grid (:doc:`terrain_coordinates`).
- :math:`c_r`: Rayleigh coefficient (``bc.rayleigh_dampcoef``); ``bc.rayleigh_fields`` selects the fields.
- :math:`s`: distance to an inflow or outflow face; :math:`w_s` is the sponge width (``bc.sponge_width``) and
  :math:`c_s` its strength (``bc.sponge_strength``).
- :math:`\theta` is relaxed to the sounding at the absolute height by default (``bc.rayleigh_height``).

In stratified flow over a hill, damp :math:`\theta` together with :math:`w`. Gravity waves exchange energy between
the two, and on ERF's 2-D hill damping :math:`w` alone let the pseudo-time iteration grow.

Nest faces
----------

A *patch* is a finer run inside a finished, coarser *parent* run (:doc:`initial_conditions`). A lateral face of the
patch with ``bc.<face> = nest`` takes its values from the parent's checkpoint (``bc.nest_parent``). The nesting is one
way: the parent is not changed. The code is in ``Source/NestFaces``.

On a nest face:

- the normal velocity on the face is the parent's (a Dirichlet row);
- the two ghost layers outside it hold the parent's :math:`u, v, w, \theta` and :math:`k`, interpolated at the ghost
  points themselves (bilinear per :math:`\zeta`-level, as the first guess);
- :math:`p'` has zero gradient through the face.

One face may hold parts where the parent's flow enters and parts where it leaves. A patch with a nest face must start
from the parent's field (``init.first_guess = nest``) or from a restart, and is refused with grid sequencing. A nest
face may be combined with ``periodic`` and ``outflow`` faces; with an ``inflow`` face only if there is an
``outflow`` face as well. A lateral wall beside a nest face is refused: the parent's flow crosses that face, and the
first guess is the parent's field.

**Net flux.** With zero-gradient :math:`p'` on every face, the pressure correction can remove a divergence only if the
flux through the boundary has no net. The parent conserves mass on its own cells, not on the patch's faces after
interpolation and with the patch's own ground. With :math:`V_f` the parent's normal velocity at face point :math:`f`,
:math:`s_f = -1` on a low face and :math:`+1` on a high face, and :math:`W_f = a_f / \Delta x_d` (:math:`a_f` the flux
area of the face in the terrain-following metric, :math:`\Delta x_d` the cell size along the normal),

.. math::

   \mathrm{net} = \sum_f s_f W_f V_f, \qquad \delta = \frac{\mathrm{net}}{\sum_f W_f}, \qquad
   V_f \leftarrow V_f - s_f\,\delta .

Every nest face point's outward velocity drops by the same :math:`\delta`, so each face's flux changes in proportion to
its area. The two ghost layers of the normal velocity outside the face drop by :math:`\delta` with it, so the gradient
across the face stays the parent's. The correction is applied once, at set-up, when no face of the patch is an
``outflow`` (an outflow face's :math:`p' = 0` carries the net instead). The sums are formed in double precision in
every build. The log reports the net as a fraction of the flux entering, and :math:`\delta`. On the hill of
``Nest_Faces`` the net is 6.3e-6 of the flux entering and :math:`\delta` is 1.5e-5 m/s.

The correction changes the parent's flow on the faces by :math:`\delta`. The log prints a warning when the net is
above 1e-3 of the flux entering. The usual cause is a patch ground that differs from the parent's on the nest faces:
with the same hill 100 m off the patch's centre the net was 2.4e-2 of the flux entering (:math:`\delta` = 5.7e-2 m/s),
and 3.0e-5 with ``terrain.nest_blend_width = 250``.

**Relaxation zone** (``bc.nest_relax_width``, off by default). In the :math:`N` cells inside each nest face, where the
parent's flow enters, :math:`u, v, w, \theta` and :math:`k` are relaxed towards the parent's field
:math:`q_\mathrm{parent}` (Davies 1976):

.. math::

   R \mathrel{+}= c\,(q - q_\mathrm{parent}), \qquad
   c = \left(1 - \frac{s}{N}\right) f\,\frac{\max(0, U_\mathrm{in})}{N\,\Delta x_d} .

- :math:`s`: distance from the face in cells (:math:`m + 1/2` for cell :math:`m` from the face, :math:`m` for a normal
  velocity point), so the weight falls linearly from 1 at the face to 0 at the zone's inner end;
- :math:`f`: ``bc.nest_relax_factor``; :math:`U_\mathrm{in}`: the parent's normal velocity on the face in the point's
  own row, positive into the patch. Where the parent's flow leaves, :math:`c = 0`;
- :math:`c` is added to the implicit diagonal. It does not depend on the solution, so the diagonal is the term's exact
  derivative;
- in a corner the two faces' rates add.

With :math:`f = 1` the rate at the face is one over the time the parent's flow takes to cross the zone. With width 0
no term is added and the residual is bit for bit that of the hard faces alone.

*What the zone does and does not do.* A difference from the parent that enters the zone at the face is carried inwards
while it decays at the rate :math:`c`. Integrated across the zone, the share that leaves the inner end is about
:math:`e^{-f/2}`: 61 % with :math:`f = 1`, 0.7 % with :math:`f = 10`. The zone therefore holds the patch to the
*parent's* (coarser) solution over its width. It is a robustness setting, for a patch whose iteration does not settle
at hard faces; it does not make the patch more accurate. On the hill of ``Nest_Faces`` every setting below reached
the steady verdict (73 to 109 iterations), and the difference from a run of the whole parent domain at the patch's
resolution grew with :math:`N f`:

.. list-table::
   :header-rows: 1
   :widths: 40 30 30

   * - Setting
     - rms :math:`|u - u_\mathrm{fine}|` (m/s)
     - max (m/s)
   * - **no zone, no blend (default)**
     - 1.37e-3
     - 1.00e-2
   * - :math:`N = 4,\ f = 1`
     - 1.41e-3
     - 1.01e-2
   * - :math:`N = 4,\ f = 100`
     - 2.70e-3
     - 1.35e-2
   * - :math:`N = 8,\ f = 100`
     - 5.34e-3
     - 3.10e-2
   * - blend :math:`W = 100` m
     - 1.55e-3
     - 1.57e-2
   * - blend :math:`W = 250` m
     - 1.75e-3
     - 1.82e-2

In the cell on the ground the zone's term is added to the row that holds :math:`k` at its wall value, as the side
sponge's is (the row's residual and diagonal both gain :math:`c`); the value the row holds moves towards the parent's
:math:`k` there by the share :math:`c / (1 + c)` of their difference.

**Ground blend** (``terrain.nest_blend_width``, off by default). The patch's terrain is its own. Where it differs from
the parent's at a nest face, the parent's values are placed at other heights above the ground than they were computed
at. With a width :math:`W > 0` the patch's ground :math:`h` becomes

.. math::

   h \leftarrow (1 - S)\,h_\mathrm{parent} + S\,h, \qquad S = t^2 (3 - 2t), \qquad t = \min(d / W, 1),

with :math:`d` the distance to the nearest nest face and :math:`h_\mathrm{parent}` the parent's ground interpolated
bilinearly at the node: the parent's ground on the face, the patch's own from :math:`W` inwards, with a continuous
slope between. The log reports the largest change of the ground. The blend adds a slope of up to
:math:`1.5\,|h - h_\mathrm{parent}| / W` to the ground, so a narrow width over a large difference steepens it. Like the
zone, the blend moves the patch towards the parent (table above); use it when the grounds differ on the faces
(the warning on the net flux), not otherwise. With a terrain mapping that smooths the levels (``terrain.smoothing`` 1 or 2) the node heights above the
ground can still differ from the parent's at the face, because the smoothing sees each mesh's own ground.

**Restart.** The checkpoint holds a record of the nest settings: the types of the four lateral faces, the parent's
step and cells, a digest of the parent's content (its ``Header`` and the headers of its field and height files, which
hold every box's minimum and maximum), the patch's place in the parent, the zone's width and factor and the blend
width. A restart with other settings, or with another parent, stops, naming ``io.restart``. The parent's path is not
part of the record: the parent may be moved.

**Limits.**

- Results on different numbers of ranks or box layouts agree to the steady tolerance, not bit for bit (the line sweeps
  are local to a box).
- Where the parent's flow leaves, a nest face still holds the parent's normal velocity and ghost values. With matched
  ground this over-specified part showed no layer in the tests; a patch whose flow differs from the parent's at such a
  part (a wake the parent does not resolve) meets the parent's values there. An ``outflow`` face is the alternative.
- The Rayleigh layer and the side sponge relax :math:`u, v` and :math:`k` towards the ``bc.inflow_*`` profile, not the
  parent's field. With a nest face and no ``inflow`` face the Rayleigh layer may damp :math:`w` and :math:`\theta`
  only (``bc.rayleigh_fields = w`` or ``w theta``) and the sponge is refused.
- Buildings (``buildings.type``) are refused with a nest face.
- A canopy (``canopy.type``) on the patch: the parent's values on the nest faces feel the drag only if the parent run
  had the same canopy. This is not checked, and the warning for canopy cells beside an ``inflow`` face does not count
  nest faces. Give the parent the patch's canopy, or keep the forest a few cells inside the nest faces.
- The nest faces are frozen: with ``steady.stage = unsteady`` the patch's boundary values do not change in time. The
  tests cover ``steady.stage = simplec``; ``simplec_newton`` and ``unsteady`` use the same rows, fills and diagonal and
  are accepted, but have no nest test yet. A patch is solved on level 0 only (grid sequencing is refused).

Pressure
--------

The pressure correction is periodic where the mesh is periodic, Dirichlet (:math:`p = 0`) on ``outflow`` faces and
Neumann on all other faces (``nest`` faces included). The types follow ERF's ``get_fft_bc``.

Ghost cells
-----------

The solution has two layers of ghost cells outside each face. With interior values :math:`c_0, c_1` and ghost
values :math:`g_1, g_2` along the face normal:

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Rule
     - Ghost values
   * - even
     - :math:`g_1 = c_0,\ g_2 = c_1`
   * - odd
     - :math:`g_1 = -c_0,\ g_2 = -c_1` (0 on the face)
   * - value :math:`V`
     - :math:`g_1 = 2V - c_0,\ g_2 = 2V - c_1`
   * - profile
     - :math:`g = P(h)` at the ghost's height above the ground
   * - gradient :math:`G`
     - the centred difference at :math:`c_0` equals :math:`G`

The directions are filled in the order x, y, z, so corner ghosts get values too. Each direction fills every field
at both faces on all boxes of a rank in one launch (audit A9); on a domain no wider than the two ghost layers it is two
launches, the low face first. The values are those of one launch per field, face, box and layer.
Every copy of a ghost cell then holds one value. In its outermost lateral ghost layer a box lacks some of what the
rules read (the node heights have one lateral ghost layer, the wall law's gradient two), so after the fill each cell
outside the domain takes the value of the box whose column it is in, and the ghost images of cells inside the domain
are refreshed from the neighbouring boxes (audit A3). No valid cell changes. The results then do not depend on the box
layout. Before, they did where an inflow face crosses terrain at a box seam and ``turb.realizability`` is on: the
cap's strain in the first lateral ghost layer read inflow face points computed from the clamped heights.
