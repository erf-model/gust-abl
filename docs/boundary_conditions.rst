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
  either way. (With the mean of the ground and top :math:`w` faces as the vertical part the fill would not be a
  fixed point: the ground face is written by the fill from the wall law's ghosts.) A reference height between
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

.. _surface-temperature:

**Surface temperature.** ``bc.surface_stability = surface_temperature`` gives the ground's temperature
:math:`\theta_s` instead of the flux: one value (``bc.surface_temperature``) or a map of ``x y theta_s`` points
(``bc.surface_temperature_file``). Each column then finds its flux from the difference between the air at the
reference point and the ground, so the flux falls where the wind falls. The column's bulk Richardson number fixes
:math:`\zeta = d/L_O`:

.. math::

   \mathrm{Ri}_b = \frac{g\,d\,(\theta_r - \theta_s)}{\theta_0\,U^2}, \qquad
   \zeta - \mathrm{Ri}_b\,\frac{c_m(\zeta)^2}{c_h(\zeta)} = 0, \qquad
   c_m = \max\left(\ln\frac{d + z_0}{z_0} - \psi_m(\zeta),\, m\right), \quad
   c_h = \max\left(\ln\frac{d + z_0}{z_0} - \psi_h(\zeta),\, m\right),

.. math::

   u_* = \frac{\kappa\,U}{c_m}, \qquad \theta_* = \frac{\kappa\,(\theta_r - \theta_s)}{c_h}, \qquad
   q = -u_*\,\theta_* = c\,(\theta_s - \theta_r), \qquad c = \frac{\kappa\,u_*}{c_h} .

- :math:`\theta_r`: :math:`\theta` at the reference point, interpolated as the wind (the first cell's own
  :math:`\theta` by default).
- :math:`\theta_s`: the ground temperature, from the value or the map.
- :math:`\theta_0` (``forcing.theta0``, default 300 K) and :math:`g` (``forcing.gravity``): the buoyancy's reference
  temperature and gravity.
- :math:`\psi_m, \psi_h`: Businger–Dyer (``qnse`` is refused with this ground condition); :math:`m`, :math:`U`,
  :math:`d` and :math:`z_0` as above.
- :math:`c`: the exchange velocity (m/s) of the heat flux.

GUST's :math:`\theta` is the virtual potential temperature :math:`\theta_v`, so :math:`\theta_s` is the ground's
:math:`\theta_v` too (``tools/hrrr_outer.py`` writes :math:`\theta_{v,s} = T_\mathrm{skin}\,(10^5/p_s)^{0.2857}
(1 + 0.61\,q_{2\mathrm{m}})` from HRRR's skin temperature).

:math:`\mathrm{Ri}_b` is bounded to :math:`[-4, 0.19]` before the solve. The lower bound is ERF's. The upper one
comes from Businger–Dyer's stable side, :math:`c_m = c_h = \ln((d + z_0)/z_0) + 5\zeta`: it gives
:math:`\mathrm{Ri}_b = \zeta c_h/c_m^2 = \zeta/(\ln((d + z_0)/z_0) + 5\zeta) < 0.2` for every :math:`\zeta`, so
without the bound a calm, cooled column has no root and no flux. With the bound it keeps a small downward flux.
The root is bracketed in :math:`[-100, 0]` or :math:`[0, 100]` (the side of :math:`\mathrm{Ri}_b`'s sign) and found
by Newton's method with bisection steps, to round-off, as the ``heat_flux`` solve. ERF's ``surface_temp`` mode uses the
stability functions of Jiménez et al. instead; on the stable side the two differ by up to about 10 % at low
wind (0.1 % on the unstable side). Thermal roughness is not modelled: :math:`z_{0h} = z_0`. Over a radiometric skin
temperature this overestimates the flux (by a factor of about 1.7 for :math:`kB^{-1} = 2` at 5 m/s and a 6 K
difference), and the bound 0.19 lies far beyond the range Businger–Dyer was fitted for (:math:`\zeta \lesssim 1`),
so strongly stable columns are only bounded, not modelled. A bracket without a root (only
possible at the bounds) gives its far end. At zero wind :math:`U` is still floored at 0.1 m/s, so the flux stays
finite.

The flux enters every place where ``heat_flux`` uses :math:`q`: the ground flux of :math:`\theta` (now
:math:`-q|A|` with each column's own :math:`q`), the ghost cells' :math:`\theta` gradient, :math:`L_O` and
:math:`k_w`. In the implicit rows of :math:`\theta` the flux is linearised in the first cell's :math:`\theta`: the
diagonal gains :math:`c\,|A|\,(1 - f)/(J\,\Delta\zeta)`, the derivative of the flux with :math:`c` held (a Picard step),
so a warm ground pulls the air towards :math:`\theta_s` without overshoot. :math:`f` is the weight of the cell above in
the reference point's interpolation: 0 at the default reference point (the first cell centre). With a
``bc.surface_zref`` above the second cell centre the first cell's :math:`\theta` does not enter the flux and the
diagonal gains nothing.

With ``simplec.theta = frozen`` or ``frozen_buoyant`` :math:`\theta` does not change, so :math:`\theta_r - \theta_s` is
fixed; the flux then changes only through :math:`U` and sets :math:`u_*`, :math:`L_O` and :math:`k_w`. The heat the
ground gives is not carried into :math:`\theta`.

On the coarse levels of grid sequencing each level reads the map at its own cell centres. A restart reads the file the
deck names again: the checkpoint does not record the map. ``init.first_guess = column`` allows a neutral column only
(the column does not take the ground temperature); ``abl.theta = surface_flux`` accepts this ground condition.
``Tests/Canonical_RANS/Surface_Temperature_Flat`` checks the wall model's :math:`u_*` against the formulas above in a
3-D run with stable and unstable columns.

The map's interpolation is the inverse-distance-weighted mean of :math:`\theta_s` itself over the
``terrain.idw_neighbors`` nearest points (weight :math:`r^{-p}`, :math:`p` = ``terrain.idw_power``); a point within
:math:`10^{-6}` m of a cell centre gives its own value, and ghost columns take the column they image.

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
file the deck names again: the checkpoint does not record the map. ``init.first_guess = column`` runs its
column on one :math:`z_0`: with a map it needs ``init.column_z0``, and the 3-D run keeps the map.

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
  ERF's Askervein inflow it reaches about 9 m²/s² at 316 m, where ERF's flow has 0.31; with a 30 m cap (ERF's default), 0.31 at 316 m; with the default 100 m, about 1.7.
- A :math:`k` column is used as given, and ``bc.inflow_tke_length`` is refused (it would change nothing).
- The equilibrium is the neutral one, also above an inversion, and the height is the height above the ground (the
  closure's wall distance with ``turb.wall_distance = terrain_height``). The roughness is ``bc.inflow_z0``, which
  should equal ``bc.surface_z0``, or with ``bc.surface_z0_file`` the map's :math:`z_0` next to the inflow face.
- A log-law wind with the capped length is not a momentum equilibrium (its stress falls with height as
  :math:`(l/\kappa(z + z_0))^2`); with the unbounded length the :math:`k` equation was out of balance instead. The
  domain adjusts either way; the capped form at least balances the :math:`k` equation.

.. _inflow-column:

**The run's own column** (``bc.inflow_profile = column``). With ``init.first_guess = column`` the run solves its
horizontally uniform column at start-up (:doc:`initial_conditions`), including the target-wind search
(``init.column_target_u``, ``_v``, ``_height``) that adjusts the geostrophic wind and the pressure gradient. With
``bc.inflow_profile = column`` the converged column is also the profile of the inflow faces, of the side sponge
and of the :math:`\theta` Rayleigh layer, and the run keeps the adjusted forcing. The log prints that forcing
(from ``Tests/Canonical_RANS/Column_Inflow``):

.. code-block:: text

   gust: inflow faces: the converged column (bc.inflow_profile = column), 33 rows, z0 0.05 m, with the target wind's forcing (kept in the checkpoint):
     forcing.ug = 9.26889497
     forcing.vg = 5.13447251
     forcing.pressure_gradient = 0 0 0

- **What it replaces.** Two runs: a periodic column run that writes its table (``io.column_guess_file``) and logs
  the forcing, then the 3-D run with ``bc.inflow_profile = file`` on that table and the forcing copied into
  ``forcing.ug``, ``forcing.vg``. The single run gives the same table to the last bit and the same forcing
  (``Tests/Canonical_RANS/Column_Inflow``). The log prints the forcing in 9 digits, so a two-run result that copies
  it differs from the single run's by about :math:`10^{-11}` (relative); the single run keeps the full precision.
- **The column's start.** The column starts from the analytic defaults (the log law with its default :math:`u_*` and
  :math:`z_0`, along x), as a periodic column deck without ``bc.inflow_*`` keys does.
- **Roughness.** The faces take the column's roughness: ``init.column_z0`` under a roughness map, else
  ``bc.surface_z0``. ``bc.inflow_z0`` is refused, so the column and the faces cannot differ. The column always has
  :math:`k`, so ``bc.inflow_tke_length`` is refused too. Under a roughness map the faces keep the column's
  :math:`z_0`, not the map's value at the face; :math:`z_0` reaches the faces only through the diagnostic
  :math:`\nu_t`, since the table holds :math:`k`. Choose ``init.column_z0`` near the map's values upwind.
- **Restart.** The checkpoint holds the column, the faces' :math:`z_0` (and whether a roughness map was given), the
  forcing, and the Coriolis setting with ``forcing.latitude`` and ``forcing.rotation_period`` (a text file
  ``InflowColumn`` in the checkpoint directory, every value in ``%.17g``). A restart with
  ``bc.inflow_profile = column`` takes the column and :math:`z_0` from it. When the column had a target wind, the
  record's forcing replaces the restart deck's ``forcing.ug``, ``forcing.vg`` and ``forcing.pressure_gradient``
  (x, y). Without a target wind, the deck's forcing must equal the record's, or the restart aborts naming the key.
  A record with a geostrophic wind needs Coriolis with geostrophic forcing in the restart deck. The restart deck must
  keep ``forcing.coriolis``, ``forcing.latitude`` and ``forcing.rotation_period``, and, without a map, the same
  ``bc.surface_z0`` and no ``bc.surface_z0_file``; a column that ran under a map needs ``bc.surface_z0_file`` (its
  :math:`z_0` cannot be compared: ``init.column_z0`` is refused on a restart). Each difference aborts, naming the
  key. The :math:`\theta` sounding (``bc.theta_ref_*``) is not compared. A restart with no
  stage reads no ``bc.*`` keys; its checkpoint carries the record of the one it read. A checkpoint without the
  record is refused under ``bc.inflow_profile = column``.
- **Refused.** A Rayleigh layer on ``u``, ``v`` or ``tke`` (also on a restart): the column's own layer would relax
  towards the profile it is solving for. Nests (``nest.names``): their record of the outer run's faces does not hold
  the column.
- **Limits** (as for ``bc.inflow_profile = file``). Below the first cell centre the faces interpolate linearly from
  the table's :math:`z = 0` row, not the column's log law. The target wind is measured with the column's log
  interpolation between cell centres, while the faces interpolate the table linearly: a face cell between the
  column's centres sees up to about 0.5 % less wind (50 m cells, target at 100 m). Nothing refuses a column whose
  wind, turned by Coriolis, leaves through a face set to ``inflow``; the log warns, per face, with the rows and
  heights (``gust: warning: bc.ylo = inflow: the column's wind leaves through this face in ...``). A wind that leaves
  aloft only can still converge; one that leaves through most of a face does not. A stratified column's surface heat
  flux is not in the record; a restart takes it from its deck.
- **Drive.** With ``abl.drive = geostrophic_inflow`` (or ``geostrophic_column``), a target wind or a restart of a
  column with one stands for ``forcing.ug``, ``forcing.vg``: no placeholder value is needed.

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
way: the parent is not changed (for two-way nests, ``nest.feedback`` other than ``none``, see :ref:`two-way-nests`
and :ref:`two-way-relax`). The code is in ``Source/NestFaces``.

On a nest face:

- the normal velocity on the face is the parent's (a Dirichlet row);
- the two ghost layers outside it hold the parent's :math:`u, v, w, \theta` and :math:`k`, interpolated at the ghost
  points themselves (bilinear per :math:`\zeta`-level, as the first guess);
- :math:`p'` has zero gradient through the face.

One face may hold parts where the parent's flow enters and parts where it leaves. A patch with a nest face must start
from the parent's field (``init.first_guess = nest``) or from a restart, and is refused with grid sequencing. A nest
face may be combined with ``periodic`` and ``outflow`` faces; with an ``inflow`` face only if there is an
``outflow`` face as well. A lateral wall beside a nest face is refused: the parent's flow crosses that face, and the
first guess is the parent's field. With ``simplec.theta = frozen`` the patch holds the parent's :math:`\theta`, as
interpolated at the start (:doc:`solution_algorithm`, prescribed potential temperature); the nest faces' ghost layers
hold the parent's :math:`\theta` as before.

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
with the hill 100 m off the parent's (measured with the Agnesi hill at the patch's own centre; the code places it
at the parent's) the net was 2.4e-2 of the flux entering (:math:`\delta` = 5.7e-2 m/s), and 3.0e-5 with
``terrain.nest_blend_width = 250``. The centre is the parent's domain centre, so a patch of a patch places the hill
at the centre of its own parent, which is itself a patch: in such a chain keep every patch centred on the outermost
domain, or use ``terrain.type = xyz``.

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
(the warning on the net flux), not otherwise. With STF (``terrain.smoothing = 1``) the node heights above the
ground can still differ from the parent's at the face, because the smoothing sees each mesh's own ground; BTF and
Sullivan set each node from its own ground.

**Buildings in a patch** (``buildings.type``). The patch builds its buildings from ``buildings.file`` on
its own cells, as a run of the whole domain on the same cells would: the same solid fractions, fixed faces and wall
cells (unit test ``NestBuildings.PatchHasTheWholeDomainsFractionsAndWalls``). The parent's values on the nest faces
know nothing of the buildings, so every building tile (a point of the height map with a height above 0 and the
rectangle it owns, :doc:`buildings`) must keep clear of the nest faces by

.. math::

   g_d = 2\,\Delta x_d + \max\left(N\,\Delta x_d,\; W\right)

for the faces normal to :math:`d`, with :math:`\Delta x_d` the patch's cell size, :math:`N` = ``bc.nest_relax_width``
and :math:`W` = ``terrain.nest_blend_width``. The zone and the blend hold the patch to the
parent's field and ground, which have no building; the two cells beyond them keep a building's wall cells and their
reference cells out of both, and away from the face's fixed row and the two ghost layers the parent fills. A tile within :math:`g_d` of a nest face stops the
run, naming the tile's x and y extents, its top, the gap and the face (up to eight tiles, then a count). A tile edge on
the gap's inner edge is clear. Along a face, a tile beside the patch in a non-periodic direction does not count. The
check uses the tiles, so it does not depend on the box layout or the rank count. Buildings need z-line sweeps
(``steady.sweep = zjacobi``), as in a whole-domain run.

The gap is a minimum for the solver, not for accuracy. Around a 30 m cube (``Tests/Canonical_RANS/Nest_Shapes``) the
patch came within 0.5 % of the uniform fine run near the cube only with the nest faces 85 m (2.8 cube widths) from its
walls; see *Where to put the nest faces* below. For a patch around buildings, keep the faces about three building
widths from the walls.

**Restart.** The checkpoint holds a record of the nest settings: the types of the four lateral faces, the parent's
step and cells, a digest of the parent's content (its ``Header`` and the headers of its field and height files, which
hold every box's minimum and maximum), the patch's place in the parent, the zone's width and factor and the blend
width, and the refinement ratio (``bc.nest_ratio``, recorded only when it is not 2, so checkpoints written before the
key keep their record). A restart with other settings, or with another parent, stops, naming ``io.restart``. A
restart with another ratio stops earlier: the patch's cells fit only one ratio, and the placement check names
``bc.nest_parent``. The parent's path is not
part of the record: the parent may be moved.

**Where to put the nest faces** (measured). A patch was compared with a uniform run of
the whole domain on the patch's cells, on a 150 m round Agnesi hill of half-width :math:`a` = 250 m
(``Tests/Canonical_RANS/Nest_Validation``) and on Askervein (``Tests/Real_Terrain/Askervein_Nest``). :math:`D` is the
half-width of a square patch centred on the hill. On the hill the difference is the RMS of the velocity-vector
difference over the cells within :math:`a` of the centre and below 100 m, as a share of the speed scale
:math:`(u_*/0.4)\ln(h/z_0)` = 7.31 m/s (the inflow speed at the hill's height :math:`h` is 7.14 m/s, with the solver's
:math:`\kappa` = 0.41); *lee* is a square of half-width :math:`a` centred :math:`2a` downwind; *faces* is the two rings of
cells beside the nest faces. One nest from a 50 m parent against the 25 m run:

.. list-table::
   :header-rows: 1
   :widths: 28 12 12 12 12 12 12

   * - Run
     - Cells
     - Iterations
     - Hill (%)
     - Lee (%)
     - Faces (%)
     - Faces, largest (%)
   * - 25 m, whole domain (the reference)
     - 196 608
     - 91
     -
     -
     -
     -
   * - 50 m, whole domain (the parent)
     - 49 152
     - 107
     - 0.700
     - 0.541
     - 0.133
     - 1.826
   * - Patch, :math:`D = 1.6a`
     - 16 384
     - 124
     - 0.082
     - 0.268
     - 0.126
     - 1.211
   * - **Patch,** :math:`D = 2.4a`
     - 36 864
     - 131
     - 0.019
     - 0.247
     - 0.125
     - 1.391
   * - Patch, :math:`D = 3.2a`
     - 65 536
     - 125
     - 0.011
     - 0.102
     - 0.128
     - 1.551
   * - Patch, :math:`D = 4a`
     - 102 400
     - 122
     - 0.008
     - 0.031
     - 0.131
     - 1.690

(The lee region reaches :math:`2.9a` downwind, so the patches of :math:`D = 1.6a` and :math:`2.4a` hold only part of
it, ending on a nest face; over those parts the parent is 0.440 and 0.508 %. The parent's other numbers are over the
regions of the :math:`D = 4a` patch.)

- Over the hill the patch is within 0.02 % of the uniform fine run once the nest faces are :math:`2.4a` from the
  hill's centre, against 0.70 % for the parent.
- Beside the nest faces the patch keeps the parent's difference from the fine run (0.13 % RMS; 1.2 to 1.7 % at the
  largest point, where the hill's wake leaves through a face).
- *Rule (from this one case: attached flow, one wind direction,* :math:`a` *= five parent cells):* keep every nest
  face at least one hill half-width beyond the region whose results are used. Faces :math:`0.3a` to :math:`0.6a`
  beyond the region gave 0.08 to 0.10 %; faces :math:`1.1a` to :math:`1.4a` beyond gave 0.02 to 0.03 %. Whether the
  distance scales with the hill or with the parent's cells is not separated by this case.
- A chain of three nests from a 100 m run (50, 25 and 12.5 m cells; :math:`D = 4a`, :math:`3.2a`, :math:`2.4a`) was
  within 0.044 % of the uniform 12.5 m run over the hill (the 100 m run: 3.11 %), with 250 880 cells in the four runs
  against 786 432. No single 12.5 m nest from a uniform 25 m parent was run, so what the chain adds is not separated
  from the finer cells.
- The relaxation zone (3 and 5 cells, factor 1, at :math:`D = 2.4a`) changes none of these numbers by more than
  0.001 %.
- To the steady verdict a patch takes more iterations (122 to 131) than the uniform fine run started from the inflow
  profile (91). The verdict is relative to each run's own first residuals, which the parent's field lowers (u: 0.015
  against 0.37), and in these neutral runs the patch's last equation to pass is θ at its floor. The patch reaches
  the fine run's final momentum residuals between iterations 70 and 80, and the two converge at the same rate. The
  saving is in the cells.
- On 1, 2 and 4 ranks (:math:`D = 2.4a` in four boxes) the patch takes 127 iterations and the fields agree to
  7e-14 of the speed scale.

On Askervein (40 m parent, 20 m patches centred on the hilltop, ``terrain.nest_blend_width = 200``; the masts of
lines A and AA lie within 870 m of the hilltop in x and in y) the 10 m speeds at the 30 masts differ from the uniform 20 m run by
0.076 %, 0.035 % and 0.015 % RMS for :math:`D` = 1000, 1400 and 2000 m (180 000, 352 800 and 720 000 cells against
1 620 000), and by 1.21 % for the parent. At :math:`D` = 1000 m one mast (line A, 130 m from a face) lies in the ground
blend, where the patch's ground differs from the reference's by up to 0.05 m; without it the RMS is 0.077 %. All
numbers are from macOS arm64, double precision, four ranks.

Around five canonical shapes (``Tests/Canonical_RANS/Nest_Shapes``) a patch at 2 2 1
was compared with the uniform fine run at three distances :math:`D` from the shape's centre to the nest faces
(velocity RMS over the shape, % of the shape's speed scale; the parent in brackets):

.. list-table::
   :header-rows: 1

   * - Shape
     - :math:`D`
     - Patch
   * - 30 × 30 × 20 m cube (Building_Cube), 5 m cells
     - 40, 70, **100 m**
     - 4.60, 1.39, **0.45 %** (18.4 %)
   * - round Agnesi hill, h = a = 150 m, 25 m cells
     - 300, 600, **900 m**
     - 0.57, 0.13, **0.04 %** (3.4 %)
   * - tanh escarpment, 100 m high, 25 m cells
     - 300, 600, **900 m**
     - 0.019, 0.005, **0.004 %** (0.53 %)
   * - plate 4 m long, 1 m above the ground, 0.5 m cells
     - 4, 8, 12 m
     - 0.49, 0.18, 0.08 % (2.5 %)
   * - 2-D Agnesi ridge, h = a = 150 m (slope 33°), 25 m cells
     - 300, 600, 900 m
     - 5.65, 8.14, 8.61 % (8.50 %)

The patches took about as many iterations as the uniform fine run, and on 1, 2 and 4 ranks (four boxes) the fields
agree to 2.3e-10 of the speed scale or better. Behind the ridge the flow separates, and the steady problem on these
cells has more than one converged solution: on the same domain, with the same inflow and outflow faces and cells, a
run from the log-law profile ends its reversed flow near the ground at x = 2312.5 m, like the uniform run, and a run
from the parent's field (``init.first_guess = nest``) at 2187.5 m, like the patches (8.56 % against 0.18 %). A patch
always starts from its parent's field, so where the steady state is not unique it can settle on another one than the
fine run. The README of the case gives the runs.

**Refinement ratio.** ``bc.nest_ratio`` (2, 4, 8 or 16; default 2) sets the patch's refinement of the parent in x and
y; z is not refined. With ratio :math:`r` the patch's cells are the parent's divided by :math:`r`, and its faces lie
on the parent's cell faces. The interpolation weights are multiples of :math:`1/(2r)`; for a power of two they are
exact in binary and mirror-symmetric. A field linear in x and y is reproduced to round-off (unit tests at 2, 4, 8 and
16). Other ratios are refused. Instead of one large ratio, a patch can be the parent of a finer patch (a chain); each
link is an ordinary one-way nest. Around the cube from a 20 m parent, against the 5 m uniform run (velocity RMS over
the cube's region, % of the speed scale; four boxes, two ranks):

.. list-table::
   :header-rows: 1

   * - Run
     - Cells
     - RMS
   * - 20 m parent (whole domain)
     - 6 912
     - 42.0 %
   * - 5 m patch, ratio 4, half-widths 70 m × 60 m
     - 16 128
     - 3.97 %
   * - 10 m patch, ratio 2 (90 m × 100 m), then 5 m patch, ratio 2 (70 m × 60 m)
     - 8 640 + 16 128
     - **2.60 %**
   * - 5 m uniform (reference)
     - 110 592
     - --

The ratio-4 patch and the chain differ by 1.40 %. The chain's 10 m link gives the 5 m patch's faces a finer parent
field, and the faces are where the two differ most. The largest difference from the reference in the 10 m ring of
cells along the nest faces is 32 % for the ratio-4 patch, whose faces take the 20 m parent's values, and 15.5 % for
the chain (31 % for the parent, in its 20 m ring). In this case (one shape, one patch size, ratio 4) the chain was
closer to the reference, at the cost of one more run and 54 % more cells (8 640 + 16 128). A ratio-4 patch with a
wider gap between its faces and the cube was not run. The ratio-4 patch on 1, 2 and 4 ranks gives the same field to
:math:`10^{-8}` of the speed scale (CTests ``rt_nest_shapes_cube_ratio4*`` and ``rt_nest_shapes_cube_chain2``).

**Solar farm.** Six rows of flat panels, 4 m long in x, 0.1 m thick and 1.5 m above the ground,
10 m apart and spanning the periodic width, on a 160 m × 2 m domain with 0.02 m layers up to 3.2 m. The 1 m parent
is refined to 0.25 m in two ways: a chain (a 0.5 m middle patch of half-width 64 m, then 0.25 m patches of it) and
ratio 4 straight from the parent, each at :math:`D` = 32, 40 and 48 m from the array's centre (the outer panel edges
are :math:`D` − 27 m from the faces). Velocity RMS over the array, % of the speed scale, against the 0.25 m uniform
run:

.. list-table::
   :header-rows: 1

   * - Run
     - :math:`D` = 32, 40, 48 m
   * - 1 m parent
     - 1.845 %
   * - 0.5 m middle patch (:math:`D` = 64 m)
     - 0.910 %
   * - 0.25 m chain
     - 0.092, 0.036, **0.011 %**
   * - 0.25 m ratio 4
     - 0.196, 0.082, 0.042 %

The chain and the ratio-4 patch differ by 0.287, 0.117 and 0.052 % at the three :math:`D`; the 4 4 1 patches took
10 to 16 % more iterations. Against the reference the chain is closer at every :math:`D` (0.011 against 0.042 % at 48 m), as around the cube. On 1, 2 and 4 ranks (four boxes) the chain's
patch at :math:`D` = 40 m gives the same field to 1.4e-15 of the speed scale. The solar runs use
``simplec.pressure_reltol`` = 1e-6. At the default 1e-3 the k residual rose by factors of 3 or more after iteration
500, from the first layer behind the array: there the ground's k row, with the weight :math:`u_*/d_0` (about
35 s\ :sup:`-1`, :math:`d_0` = 0.01 m), shows the change of the first cell's velocity between iterations, which the
pressure solve's error at 1e-3 makes non-monotone. The fields of the two tolerances differ by at most 3.5e-8 of the
speed scale in velocity, the iterations to the steady verdict hardly change, and the 0.25 m patches take 3 to 4
times as long per iteration (the case's README gives the numbers). The 0.25 m reference
(1.1 million cells) takes about five hours (2574 iterations) on eight cores of a laptop, so the solar farm is run by hand, not in the
nightly study.

What these comparisons do not cover: the reference is the same solver with the same closure, so they measure the
nesting error, not accuracy; neutral flow, refinement 2 2 1 (twice in a chain for the solar farm; 4 4 1 for the cube
and the solar farm above) with the same layers in every run, one wind direction per case; the patches are centred on
the shape and away from the domain's own boundaries (the solar farm's middle patch by 16 m);
the zone was tried with factor 1 at one patch size.

**Limits.**

- Results on different numbers of ranks or box layouts agree to the steady tolerance, not bit for bit (the line sweeps
  are local to a box).
- Where the parent's flow leaves, a nest face still holds the parent's normal velocity and ghost values. With matched
  ground this over-specified part showed no layer in the tests; a patch whose flow differs from the parent's at such a
  part (a wake the parent does not resolve) meets the parent's values there. An ``outflow`` face is the alternative.
- The Rayleigh layer and the side sponge relax :math:`u, v` and :math:`k` towards the ``bc.inflow_*`` profile, not the
  parent's field. With a nest face and no ``inflow`` face the Rayleigh layer may damp :math:`w` and :math:`\theta`
  only (``bc.rayleigh_fields = w`` or ``w theta``) and the sponge is refused.
- Buildings in a patch must keep clear of its nest faces (above); a building cut by a nest face is not supported.
- A patch starts from its parent's field. Where the steady problem has more than one solution (separated flow behind
  a steep ridge, above), the patch can settle on another one than a uniform fine run from a profile.
- A canopy (``canopy.type``) on the patch: the parent's values on the nest faces feel the drag only if the parent run
  had the same canopy. This is not checked, and the warning for canopy cells beside an ``inflow`` face does not count
  nest faces. Give the parent the patch's canopy, or keep the forest a few cells inside the nest faces.
- The nest faces are frozen: with ``steady.stage = unsteady`` the patch's boundary values do not change in time. The
  tests cover ``steady.stage = simplec``; ``simplec_newton`` and ``unsteady`` use the same rows, fills and diagonal and
  are accepted, but have no nest test yet. A patch is solved on level 0 only (grid sequencing is refused).

.. _two-way-nests:

Two-way nests (``nest.feedback = defect``)
------------------------------------------

With ``nest.names`` in the outer deck (:doc:`inputs`, ``nest``) the nests run one way after the outer domain. With
``nest.feedback = defect`` the run goes on: the outer domain and all its nests stay in memory and iterate together,
one *exchange* per iteration (``Source/JointRun``, ``Source/NestDefect``).

The use in mind is one coarse domain over a region with several sites of interest, each with its own fine nest
(``nest.names = a b c``). The nests run together in the same exchange; each one's flow reaches the coarse domain
around it, and through it the other nests' faces.

Each exchange:

1. **Up.** Each nest's state :math:`q_f` is restricted onto the outer cells under it (its *footprint*); the outer
   state with these values on the footprint is :math:`q^*`. The outer residual at :math:`q^*`, without any
   source, gives the defect source

   .. math::

      S = \chi_s\, R_c(q^*), \qquad q^* = q_c \text{ with } I q_f \text{ on the footprint}.

2. The outer domain takes one iteration of :math:`R_c(q_c) - S = 0`.
3. **Down.** Each nest's parent images and nest faces are taken from the outer state now.
4. Each nest takes one iteration.

With :math:`n` nests, step 1 restricts the cell fields, :math:`u` and :math:`v` of every nest first, then fills the
outer ghost cells, then forms :math:`w` of every nest, then fills them again (``restrict_patches``): two outer ghost
fills per exchange for any number of nests, against :math:`2n + 1` if each nest were restricted in full before the
next. :math:`w` reads the outer :math:`u` and :math:`v`, so on a face that two touching footprints share it sees the
face's final flux, and the outer :math:`\Omega` on every footprint is the restricted nest :math:`\Omega` (unit test
``TouchingNestsKeepTheRestrictedOmegaOnBothFootprints``). The outer node heights are copied into :math:`q^*` once, as
they do not change during the joint run.

Such a shared face is a lateral face of both nests. Each nest's face takes the outer flow interpolated, minus that
nest's own net-flux shift (``Source/NestFaces``), so the two restrictions of the face differ slightly. The face keeps
the restriction of the nest listed later in ``nest.names``, and the outer cells next to it on the other nest's side
keep a continuity mismatch of :math:`a_x \Delta u` (no order can make both sides divergence-free). With each nest
restricted in full before the next, the earlier nest's column next to the face would also have :math:`\Omega` off by about
:math:`z_x \Delta u / 2`. On the two-way hill deck with two touching nests (``crest`` and ``lee``) the defect run is
the same bit for bit in both orders: the changed column is outside the stencil of the source's region there. With
``nest.feedback = relax`` (factor 4) or with ``turb.realizability = 0.309`` it reads that column, and the steady
fields differ by at most 1.7e-6 m/s and 1.0e-8 m/s in :math:`u`, with the same number of exchanges.

Here

- :math:`R_c`: the outer residual of :math:`u, v, w, \theta` and :math:`k` (with the stage's upwind blend);
- :math:`q_c`: the outer state; :math:`q_f`: a nest's state;
- :math:`I`: the restriction at 2 2 1. Cell values are the :math:`J`-weighted mean of the nest cells; :math:`u` and
  :math:`v` are the nest's face fluxes summed over the outer face's area; :math:`w` comes from the restricted
  vertical flux :math:`\Omega`. Both areas are the open part of the face, :math:`a_d (1 - \beta_d)`
  (``Mesh::flux_area``, :doc:`buildings`), as in the divergence and the pressure correction; an outer face closed in
  full gets 0 (unit test ``RestrictionWeightsTheOpenFaceAreas``). Without buildings the open part is the whole face.
  Joint nests are still refused with buildings. The restricted field is divergence-free on the outer cells when the
  nest's field is divergence-free on its own, so continuity gets no source. The pressure is the :math:`J`-weighted mean of the nest cells,
  shifted by one constant per nest so that its mean over the footprint is the outer pressure's mean there (the two
  pressures have different levels). The shift uses plain means, not :math:`J`-weighted ones. The level does not
  reach :math:`S`: the momentum rows of :math:`\chi_s` see the pressure only through differences inside the
  footprint, where every value carries the same constant (unit test ``ThePressureLevelDoesNotReachTheSource``);
- :math:`\chi_s`: 1 on the outer points at least 2 outer cells inside the footprint's lateral faces (all
  :math:`z`), 0 elsewhere. The 2-cell ring between keeps the nest's own face values out of the cells of :math:`S`.
  The rows of :math:`u` (and :math:`v`) on the edge faces of :math:`\chi_s` still read the nest's face two faces
  away through the upwind reconstruction of their momentum flux.

At a joint fixed point :math:`R_c(q_c) = 0` outside :math:`\chi_s` and :math:`R_c(q_c) = R_c(q^*)` inside it, so the
outer domain carries the nests' solution on :math:`\chi_s` and its own equations see it through their stencils. On the
ring :math:`q^*` takes :math:`I q_f` too, so the rows of :math:`\chi_s` see the nest's state alone at :math:`q^*`;
at the joint steady state the outer state on the ring is its own, close to :math:`I q_f` but not equal to it, and
the rows of :math:`\chi_s` within two cells of the ring see that difference, so the outer state differs from
:math:`I q_f` there too, less so further in. Taking
:math:`q^*` = :math:`I q_f` on :math:`\chi_s` only (the ring at the outer state) makes :math:`q_c = I q_f` an exact
fixed point, but on the Columbia Gorge case the exchanges then settled into a cycle near a residual of 1e-4 and did
not reach steady in 600 exchanges; with :math:`q^*` over the whole footprint the same case is steady in 276. This is the
local defect correction of Ferket & Reusken (1996) and Anthonissen et al. (2003), applied to the steady RANS
residual. :math:`S` is constant during the outer iteration and changes only at an exchange. The rows a boundary
condition or the wall law replaces replace :math:`S` too.

**The nests' coordinate surfaces.** In the joint run a nest's coordinate surfaces follow the outer domain's, with
the nest's own ground detail decaying with height (a SLEVE-like coordinate, Schär et al. 2002):

.. math::

   z_k = z_{p,k} + (h - h_p)\, \frac{\sinh\big((H - \zeta_k)/s\big)}{\sinh(H/s)}, \qquad
   s = \max\Big(r\,\Delta x,\; \frac{2 \max|h - h_p|}{m}\Big),

with :math:`z_{p,k}` the outer domain's surface :math:`k` interpolated to the nest's node (bicubic, Catmull-Rom), :math:`h` the
nest's ground, :math:`h_p = z_{p,0}`, :math:`\zeta_k` the z-levels, :math:`H` the top level, :math:`r\,\Delta x` the
outer domain's smaller horizontal cell size and :math:`m` the smallest slope :math:`\partial z_p / \partial \zeta`
of the outer surfaces over the nest. The ground stays the nest's and the top is the outer domain's. The layers keep
:math:`\partial z / \partial \zeta \ge m - \tfrac{m}{2}\coth(H/s)`, at least :math:`0.48\,m` when :math:`H \ge 2s`;
a nest whose ground stands far above the outer one's can still invert them. With
:math:`s = 2\max|h - h_p|/m` the bound :math:`m - \max|h - h_p| \coth(H/s)/s` is positive only while
:math:`\max|h - h_p| < m H / (2\,\mathrm{atanh}(1/2)) \approx 0.91\, m H`; a larger :math:`s` (:math:`r\,\Delta x`)
only helps. The bound is sufficient, not necessary, so the levels are checked layer by layer once they are set: a
layer :math:`z_{k+1} - z_k \le 0` at a node stops the run with an error that names the nest (``nest.<name>``), the
node's position, the layer's thickness, :math:`m`, :math:`\max|h - h_p|`, :math:`H`, :math:`s` and the bound, and
suggests raising ``mesh.prob_hi`` z or refining the outer mesh. This check comes before the general mesh check, whose
message names no nest. **Why:** the restriction takes
:math:`w` from the vertical flux :math:`\Omega` summed over the nest faces. With the nest's own terrain-fitted
surfaces its finer ground detail decays only as :math:`1 - \zeta/H` (BTF), and between the outer nodes the nest's
surfaces aloft are not the outer domain's. The restricted :math:`\Omega` then carries the slope-velocity covariance
of that detail as vertical motion: on the Columbia Gorge case the outer :math:`w` error outside the footprint grew to
0.035 m/s (one way 0.026) and downstream to 0.053 m/s (one way 0.031). With the surfaces above, the nest's surfaces
lie on the outer domain's within a few :math:`s` of the ground, and the errors are 0.025 and 0.026 m/s. The cost is
over the footprint: there the outer :math:`w` error grows from 0.044 m/s one way to 0.068 m/s, because the outer
cells see the nest's flow on surfaces that bend faster than their own. The log gives
:math:`\max|h - h_p|`, :math:`s` and the largest move of a node. A one-way nest keeps its own terrain-fitted surfaces.

**Pseudo-time step.** When the joint run starts, each domain's CFL is held at no less than the value its one-way run
ended at. The first exchange raises every domain's residual above its one-way reference, where the SER law would cut
the CFL; with the pseudo-time step the SIMPLEC denominators and the Poisson coefficients shrink, and the pressure
correction grows in proportion. A floor at ``steady.cfl0`` instead kept the nest's step a ``cfl_max / cfl0`` fraction
of the outer domain's, and the exchanges settled into a cycle that did not reach steady.

**Stopping.** The run is steady when every domain's driver is steady after the same exchange (each domain's ``steady.*``
criteria, with :math:`S` and the faces as they are then) and the feedback has settled. Each exchange's change from the
last,

.. math::

   \delta_n = \max_{q \in \{u, v, w, \theta, k\}} \frac{\max|F_q^{n} - F_q^{n-1}|}{\max(\sigma_q, \phi_q)},

is tested over the outer domain's ``steady.monitor_window`` exchanges :math:`W` (one exchange is one iteration) as a
monitor is over that many iterations (:ref:`convergence-criteria`):

.. math::

   \sum_{m=n-W+1}^{n} \delta_m \le \epsilon, \qquad
   \frac{d_2\,\rho}{1 - \rho} \le \epsilon, \quad \rho = \frac{d_2}{d_1},

with :math:`\epsilon` the outer domain's ``steady.monitor_rel_tol`` and :math:`d_1, d_2` the sums of :math:`\delta`
over the window's two halves (the second test is skipped when both are below :math:`10^{-3}\epsilon`, and fails when
:math:`\rho \ge 1`). Here :math:`n` is the exchange, :math:`F` the defect source :math:`S` (maxima over its region),
:math:`\sigma_q` the largest :math:`|S_q^n|` and :math:`\phi_q` the outer domain's
``steady.residual_tol / steady.monitor_rel_tol`` times the residual scale of :math:`q` (:math:`U^2/H`,
:math:`U\theta_s/H`, :math:`U k_s/H`; :ref:`convergence-criteria`), so a change of ``steady.residual_tol`` times the
scale passes at any size of :math:`S`. With ``nest.feedback = relax``, :math:`F` is the target :math:`q^*` where the
rate is positive, :math:`\sigma` the largest :math:`|u^*|, |v^*|, |w^*|` for the three velocity components, 0 for
:math:`\theta` and the largest :math:`k^*` for :math:`k`, and :math:`\phi` the fixed scales :math:`U`,
:math:`\theta_s`, :math:`k_s` of :ref:`convergence-criteria` (10 m/s, 10 K, 0.1 m\ :sup:`2`/s\ :sup:`2`; not inputs), so :math:`\theta` is measured in kelvin, not
relative to its ~300 K. :math:`\delta` is infinite at the first exchange, and an infinite :math:`\delta` empties the
window. **Why:** :math:`S` is a constant of each step, so a domain's own verdict does not see how much :math:`S` moved
between exchanges. The sum bounds the feedback's change over the window, as a monitor's change is bounded over its
window; the extrapolation keeps a slow drift from passing on a small change per exchange. The floor :math:`\phi`
keeps a source below what the residual test resolves (:math:`\theta` in a neutral run) from holding the run. The run
also stops when one domain diverges, or after ``steady.max_iterations`` exchanges; when every domain was steady at the
last exchange but the feedback had not settled, the log says so. The joint run starts once the outer domain is steady,
also when a nest's one-way run was not; GUST stops with an error after the joint run, naming every domain that is not
steady. The log gives, for the first ten exchanges and then every tenth, the largest :math:`|S|` per equation,
:math:`\delta_n` (``S change``; ``target change`` with ``relax``) with the window's sum, and the largest change of a
nest's face velocities. A domain name longer than 24 characters shows in these lines as its first 12 and last 9 characters
around ``...``; the per-domain steady line after the run gives it in full.

**Restart.** A joint run that stops after ``steady.max_iterations`` exchanges, not steady, writes every domain's
plotfile and checkpoint, and a joint record, the file ``Joint``, in the outer domain's checkpoint. The record holds
the feedback mode, the exchanges taken, each driver's CFL floor (above), the stop test's window of changes
:math:`\delta_n` (above), ``nest.feedback_factor`` (0 with ``defect``) and each nest's checkpoint; beside it the outer checkpoint holds the feedback of the last
exchange (``Level_0/joint_feedback_*``: :math:`S` with ``defect``, the target :math:`q^*` with ``relax``). The error
message names that checkpoint. To go on, rerun the same deck with ``io.restart`` set to it, without the ``init.*``
keys, and with a larger ``steady.max_iterations``; the joint run counts its exchanges over the run, so the value is
the total. Every domain then reads its checkpoint and its CFL floor and joins the joint run at once, with no one-way
iterations. The log says ``joint run: restarted after exchange N``.

- ``nest.names`` must list the record's nests in the same order. Each nest's ``lo`` and ``hi`` must be the ones
  its checkpoint was written with (its node heights and nest-face record are checked against them).
- The record names each nest's checkpoint as written, relative to the directory the run started in. Restart from
  the same directory: a recorded checkpoint that is not there stops the run before anything runs, naming
  ``io.restart`` and the missing paths.
- The restart takes ``nest.feedback_factor`` from the deck or its default, not from the record. With the same
  ``nest.feedback`` and another factor the run goes on with the deck's value, the window starts empty and the log
  gives both values. With another ``nest.feedback`` the factor is the new one's, given or default (from ``relax`` to
  ``relax_consistent`` without the key, 0.5 becomes 1), and the log names the two feedbacks.
- A record written before the factor was added (``GUST_JOINT 1`` or ``2``) holds none. With the same
  ``nest.feedback``, other than ``defect``, the log says so. A record of version 2 then restores its window whatever
  the deck's factor: if that is not the recorded run's, the window holds another factor's changes and the run may
  stop early, so give the factor the run had. A ``relax`` run without the key from the time its default was 1 needs
  ``nest.feedback_factor = 1``.
- A GUST build from before ``GUST_JOINT 3`` stops on a record of version 3, also with ``nest.feedback = none``.
- A restart needs the precision the checkpoint was written in; the other one stops, naming ``io.restart``.
- ``convection.hybrid_sigma_max`` (the blend) is refused with ``io.restart``, as for any restart: the blend is not
  in the checkpoint. A joint run with the blend therefore restarts without it, and its error message says so.
- On the same number of ranks a restart continues the run bit for bit: the regression ``rt_nests_joint_restart``
  stops a run after 150 exchanges, restarts it, and compares the last plotfiles and checkpoints byte for byte. Two
  things make this hold. A joint nest's boxes go to the ranks in a fixed order (a space-filling curve in rank order),
  not to the ranks with the least memory in use, which depends on what the outer domain holds when the nest is
  built. And after each exchange the nest's ghost cells are filled again from its new faces, as the first step of a
  restart fills them; the driver's local pseudo-time step reads them before the iteration does. Both change a joint
  run without a restart at round-off (the rank order changes the order of reductions; the second fill changes the
  local step), and the extra fill costs one ghost fill per nest per exchange.
- ``steady.max_iterations`` also caps each nest's one-way run. A nest whose one-way run stops there without reaching
  steady joins the joint run as it is (the log notes it). A run stopped and restarted then repeats a straight run with
  the larger value only when the stop point is above every one-way run's iterations.
- A restart restores the window and the last feedback, so its first exchange measures :math:`\delta_n` and tests
  the window as the run without the restart does, and it stops after the same exchange. A record written before the
  window was added (``GUST_JOINT 1``) still restarts, as does a restart that changes ``nest.feedback`` or
  ``nest.feedback_factor``; there the window starts empty (the log notes it), so the restart cannot stop before it is
  full again, ``steady.monitor_window`` exchanges or more later. From a version 1 record, with nothing changed, it
  never stops earlier than the run it continues and may stop later. With another feedback or factor the iteration
  changes: the run may take fewer or more exchanges, and except for a change of ``inject``'s factor alone its steady
  state may differ.
- ``nest.feedback`` may change at a restart (the log notes it). ``nest.feedback_flux`` is not in the record and may
  change at a restart without a note: the ring correction is formed again at every exchange, so nothing of it is
  restored. ``nest.damping`` is not in the record either; the deck must keep it, since it sets the nests' Rayleigh
  layer and so the joint fixed point. ``steady.max_iterations`` at or below the recorded
  exchanges, a record that names other nests, and a record that cannot be read stop the run, naming
  ``io.restart``.
- A checkpoint without a record, such as the outer domain's when its one-way run did not reach steady, restarts that
  one-way run; the nests then start from its last checkpoint, as in a first run. A joint checkpoint restarted with
  ``nest.feedback = none`` does the same (the log notes it); the outer domain's history then already holds the joint
  run's iterations.
- The checkpoints of a joint run can start a one-way run (``init.first_guess``, a parent for ``bc.nest_parent``). A
  joint nest's checkpoint does not restart as a stand-alone patch with ``bc.nest_parent``: its nest-face record
  leaves out the parent's step and content, which change at every exchange, so the check against a parent fails.

**Damping.** By default a nest damps only :math:`w` and :math:`\theta` (above). Over each footprint the joint run
carries the nests' equations, so an outer Rayleigh term on :math:`u, v` or :math:`k` would act outside the footprints
and not inside them: a jump in the forcing at each footprint's edge above ``bc.rayleigh_zdamp``.
The joint run (``nest.feedback`` other than ``none``) therefore needs ``bc.rayleigh_fields = w`` or ``w theta``, and stops with an
error otherwise, unless ``nest.damping = outer``: the nests then keep the outer deck's ``bc.rayleigh_fields`` and damp
:math:`u, v` and :math:`k` towards the same ``bc.inflow_*`` profile as the outer domain, so the forcing is the same on
both sides of each footprint's edge. On a flat, horizontally uniform ABL with an inversion
(``Tests/Canonical_RANS/Nests_Flat_ABL``, below) the joint run needs it: with the nests damping :math:`w` and
:math:`\theta` only, rolls grow in the nest below the inversion (4.5e-3 m/s departures from uniform in the buoyant
case) or the joint run does not reach steady (neutral case, 10000 exchanges). A nest has no side
sponge either: with ``bc.sponge_width > 0`` the run prints a note, and the nests are best placed outside the sponge.
Nor does a nest of ``nest.names`` have a relaxation zone or a ground blend (Nest faces, above):
each nest's ``bc.nest_relax_width`` is 0 (``Source/Nests.cpp``), and ``bc.nest_relax_width`` or
``terrain.nest_blend_width`` in the outer deck stops the run, since the outer domain has no ``nest`` face. A joint
nest's ground is its own up to the faces, and its coordinate surfaces follow the outer domain's aloft (above).

**Flat ABL.** ``Tests/Canonical_RANS/Nests_Flat_ABL`` checks the nests where the answer is known: on
a doubly periodic flat outer domain that starts from its steady column, the steady state is the column, the same at
every :math:`(x, y)`, and a nest at 2 2 1 must return the outer domain's level means. Three decks: neutral (case S
of Moeng and Sullivan 1994, as rerun by Churchfield et al. 2010), buoyant (case B) and stable (the set-up of GABLS1,
Beare et al. 2006, with a fixed surface flux of -0.0025 K m/s). S and B freeze :math:`\theta` at the column's
(``simplec.theta = frozen``; with :math:`\theta` solved, entrainment at the inversion warms the mixed layer and the
case has no steady state) and damp :math:`u, v` and :math:`w` towards the column's table from below the inversion (450 m and 700 m) with
``nest.damping = outer``. The stable deck solves :math:`\theta` and damps :math:`w` and :math:`\theta` only. On two
ranks the nest's level means are the outer domain's to 1.1e-7 m/s (S), 2e-8 m/s (B) and 3.2e-5 m/s (stable) one
way, and to 8.5e-10, 6e-10 and 6.5e-5 m/s after the joint run. In the stable case the outer domain's own stop
(``steady.residual_tol = 5e-5``) leaves its mean 0.066 m/s (near 316 m) and 0.013 K from the state it reaches
alone after 20000 iterations; the joint run, which iterates on, moves it by the same amounts. The joint stable run
also leaves a roll in :math:`u` aloft, where :math:`k` is about 0 and :math:`u` and :math:`v` are not damped (4.0e-4
m/s from uniform in both domains; the outer domain alone stays uniform to 2e-6): a :math:`u(y)` in geostrophic
balance is a neutral mode there, and its size depends on the iteration history. The outer domain's level means
differ from the column's by up to 0.15 m/s (S), 0.005 m/s (B) and 0.09 m/s (stable): the column and the 3-D stage
do not discretize the capped, stratified column alike. The case's ``README.md`` has the
details.

.. _conservative-faces:

**Conservative faces (``nest.feedback_flux = nest``).** With the default, ``outer``, the outer row just outside
:math:`\chi_s` (the *ring*) keeps its own equation, and its face on the edge of :math:`\chi_s` carries the outer
scheme's flux of :math:`q^*`. The nest's own flux through that face can differ, most on a lee slope (Limits, below).
With ``nest`` the ring's equation takes the nest's flux there, as in the conservative local defect correction of
Anthonissen et al. (2003) (``Source/NestFaceFlux``):

.. math::

   R_c(q_c) - S = 0 \text{ on the ring}, \qquad
   S = \frac{s}{J\,\Delta\xi} \Big[ (G_c - G_f) + w \big( (G_c^{up} - G_c) - (G_f^{up} - G_f) \big) \Big],

with

- :math:`G_c`: the outer flux of :math:`q^*` through the face (convective minus diffusive, as the residual forms
  them, after the wall law and the building closures);
- :math:`G_f`: the nest's flux of :math:`q_f` through the nest faces the outer face covers, summed and divided by
  their number (the mean flux density; on a staggered face of :math:`u` or :math:`v` the nest's faces at the two
  ends weigh half). Where the outer face is nodal along its normal (:math:`u` in :math:`x`, :math:`v` in
  :math:`y`) it lies inside a nest control volume, and :math:`G_f` is the mean of the two nest faces around it, an
  interpolation, not a sum;
- :math:`G^{up}`: the same fluxes with the upwind limiter, and :math:`w` the ring point's blend weight
  (``convection.hybrid_sigma_max``; :math:`w = 0` without the blend). :math:`k` has no blend (bounded upwind);
- :math:`s`: +1 on the ring below :math:`\chi_s` along the face's normal, −1 above; :math:`J` the ring point's
  Jacobian, :math:`\Delta\xi` the outer cell size along the normal.

The ring's flux through the face is then :math:`G_f + [G_c(q_c) - G_c(q^*)]`: the nest's flux plus the outer flux of
what :math:`q_c` and :math:`q^*` still differ by there (:math:`q^*` covers the whole footprint, the ring and the
margin keep the outer equations). It is the nest's flux only to that difference; a one-dimensional model of the
exchange (advection and diffusion, ratio 2) gives the nest's flux to round-off, but in two or three dimensions only
the sum over the edge of :math:`\chi_s` has to match. Only the convective and diffusive fluxes are corrected: the
ring's pressure gradient, continuity and the source terms (buoyancy, production of :math:`k`, the wall law) stay the
outer domain's, so the momentum through the face is not the nest's. The blend term matches the nest's only where the
nest's own blend weight next to the face is the ring's :math:`w`. The joint log lines end with the ring's largest
correction (``ring max |C| u ...``), so a correction that grows or oscillates shows. :math:`S` on the ring is formed
at every exchange, with :math:`S` on :math:`\chi_s`, from the outer residual at :math:`q^*` and one more residual of
each nest (the cost: one nest residual per exchange). Only the lateral faces of :math:`\chi_s`, in the directions a
nest has faces; the ratio must be even (2 for named nests). With any other ``nest.feedback`` the key is
refused. Results (RMS of :math:`u` against the uniform 25 m run, joint over one way,
``Tests/Canonical_RANS/Nests_Two_Way``, 2 ranks; the deck as it is, and the recommended setting in brackets):

.. list-table::
   :header-rows: 1

   * - Nests
     - Nest, ``outer``
     - Nest, ``nest``
     - Outer domain outside the footprints, ``outer``
     - Outer domain outside the footprints, ``nest``
   * - ``hill`` (faces 250 m and more from the crest)
     - 0.25 (0.30)
     - 0.26 (0.34)
     - 0.34 (0.31)
     - 0.33 (0.30)
   * - ``crest`` 300 to 800 m (downstream face on the lee slope)
     - 1.54 (1.68)
     - **1.03** (1.08)
     - 1.59 (1.77)
     - **0.90** (0.93)
   * - ``crest`` 300 to 900 m
     - 3.31 (3.47)
     - **1.22** (1.37)
     - 2.74 (2.95)
     - **0.96** (1.10)
   * - ``crest`` 300 to 1000 m
     - 0.40 (0.77)
     - **0.38** (0.57)
     - 0.66 (0.54)
     - **0.50** (0.46)
   * - ``crest`` to 1000 m and ``lee`` (touching)
     - 0.42 / 0.61 (0.80 / 0.50)
     - **0.37 / 0.43** (0.56 / 0.41)
     - 0.80 (0.70)
     - **0.51** (0.48)
   * - ``crest`` to 900 m and ``lee`` from 900 m
     - 3.37 / 2.38 (3.49 / 2.36)
     - **1.23 / 0.65** (1.38 / 0.69)
     - 2.19 (2.19)
     - **0.85** (0.92)

Where the nest's face is on the lee slope the joint run is now within 1.4 times the one-way nest's error, and the
outer domain outside the footprints within 1.1 times its one-way error; it was 3.5 and 3.0 times. On the hill, whose
faces are far from the crest, the nest's :math:`u` error is 3 % larger than with ``outer`` (2.03e-3 to 2.10e-3; 13 %
with the recommended setting). Each exchange costs one more nest residual; the number of exchanges changed by at most
8 % (wall times were not measured on an idle machine). The regression ``rt_nests_feedback_flux`` checks the crest to
900 m with both values and the hill with ``nest``. Two more layouts, by hand (2 ranks; :math:`u`, :math:`w`, :math:`k`
of the nest, joint over one way, ``outer`` → ``nest``): the crest nest ending at 950 m, where the flow near the ground
is slowest behind the hill, 2.04 → 0.65, 2.20 → 0.72, 2.09 → 0.78 (189 and 188 exchanges); the crest to 900 m on a
hill 150 m high, 1.96 → 0.62, 2.41 → 0.87, 2.02 → 0.69 (249 and 255 exchanges). Neither reference separates behind the
hill, so a face that carries reversed flow is not tested. In all four runs the ring's largest correction settles
without growing or oscillating.

**Limits.**

- With ``turb.nut_relax`` below 1 the outer :math:`\nu_t` blends with its own history, so :math:`R_c(q^*)` uses that
  history, not one of :math:`q^*`'s own. At a steady state the history equals :math:`\nu_t` of the outer state, so
  the steady state is the same; whether the iteration reaches it can depend on ``turb.nut_relax``
  (:ref:`nut-relax-end`).
- ``steady.probes_file`` rows stop at the end of the one-way run: the joint run writes none.
- The outer domain's turbulence plotfile after the joint run shows the residual with :math:`S` included.
- ``steady.stage = simplec_newton`` is refused: its residual is not the one the source is formed with.
- Two-way is not always closer to the fine run than one-way (Harris & Durran 2010): the outer domain's error
  outside the footprints can grow where the nest's solution differs from the outer one's. Compare with
  ``nest.feedback = none`` on a new site. A nest face on a lee slope is one such place: on the steep hill of
  ``Tests/Canonical_RANS/Nests_Two_Way`` a nest ending 100 m behind the crest (two thirds of the half-width) has a
  joint error 3.3 times its one-way one, and one ending 200 m or more behind it 0.4 times or less. A likely cause:
  the flux leaving :math:`\chi_s` is the outer scheme's flux of :math:`I q_f`, not the nest's own face flux, and
  where the gradients are strong the two differ. A 1-D model problem (second-order upwind, ratio 2) shows the same
  pattern, and a conservative variant (the summed nest flux on that face) reduced it there. That variant is
  ``nest.feedback_flux = nest`` (Conservative faces, above); on the same nest it brings the joint error to 1.2 times
  the one-way one.

**Where to put a joint nest.** ``tools/nest_placement.py`` checks a nest from the terrain file and the outer mesh
alone, before any run (``--nest X0 Y0 X1 Y1``), or searches for one of a given size (``--recommend WIDTH HEIGHT``).
It reports the share of high-misfit columns a nest covers (where the outer cells miss the ground) and, for each
face, the channels, ridges and steps it cuts. For each wind direction a face is inflow or outflow when the wind
is within 60° of its normal, and a side face otherwise. An inflow or outflow face is flagged when it cuts a
feature the outer cells cannot follow (local misfit above 5 m at 200 m cells) at more than 20° to the face
normal; the angle is measured against the face normal, not the wind. ``--direction`` takes one or more wind
directions or ``all`` (every 22.5°); a nest is clear only when it is clear for each. The thresholds scale with
the outer cell size only (the misfit is the outer mesh's error in the ground height); the tool prints the
values it used. The nest corners, the domain and the ``--recommend`` sizes must be whole outer cells. What the
runs support so far:

- The angle and misfit rules were fitted on four nests of one channel site (200 m outer cells, the wind normal
  to the faces): the nests at 25° and 28° raised the outer cells above 0.25 m/s from 12 to 663 and 817 in a
  joint run, and those at 14° and 0° to 26 and 12. They have not been tested on an independent case.
- A channel bending near the outflow face did not harm a joint run (12 to 12 cells), so a bend is a note, not a
  flag.
- A nest 4 outer cells from the outer inflow face did not converge in a joint run on the Columbia Gorge; the
  search keeps 10 outer cells from the outer faces by default, and a nest closer than that gets a note.
- The step rule was tested on an escarpment 60 m high, with 50 m outer cells (misfit at most 1.8 m) and two nests
  whose faces cut the step at 35° and 55° (``Tests/Canonical_RANS/Nest_Placement_Runs``). For winds
  along an axis the verdicts match the runs: the two flagged nests make the joint run worse than one way (up to
  1.35 times in the nest and 1.75 outside it at 55°), and the two clear ones do not (within 1.02).
  ``nest.feedback_flux = nest`` brings the flagged nests back to within 1.02 (nest) and 1.08 (outside) of one way.
  For an oblique wind (from 225°, every face inflow or outflow) the tool flags both nests, but neither joint run
  is worse than one way: that wind blows within 10° of the step's contour lines, and the tool measures the angle
  against the face normal, not the wind. The flag is then on the safe side, and the thresholds are unchanged.
  Other cell sizes have not been tested by a run. With the default step rise (0.125 dx, 25 m) no
  9.6 km x 7.2 km position on the Columbia Gorge is clear for wind from 270°; with ``--step-rise 100``, 3 are.

When no position is clear, the tool says so, recommends none and exits with status 2. Compare
``nest.feedback = none`` and ``defect`` on the best-covering flagged position, or try another size. A checked
nest with a flagged face exits with status 1, and an input error with status 3.

**Results.** RMS of :math:`u` against a run of the whole domain at the nest's resolution, one way / joint.
"One way" is the first pass of the defect run (nest levels as in the joint run) for the canonical hill, and a
separate ``nest.feedback = none`` run (the nest's own surfaces) for Askervein and the Columbia Gorge:

.. list-table::
   :header-rows: 1

   * - Case
     - Nest
     - Outer, outside the footprint
     - Cost of the joint run
   * - Steep hill (``Tests/Canonical_RANS/Nests_Two_Way``; 50 m / 25 m, 2 ranks)
     - 8.23e-3 / **2.03e-3** m/s
     - 1.44e-2 / **4.89e-3** m/s
     - 208 exchanges, 5.3 s
   * - Askervein (40 m / 20 m nest of 2 km x 2 km, 4 ranks)
     - 3.37e-3 / 3.37e-3 m/s
     - 4.62e-3 / **3.21e-3** m/s
     - 156 exchanges, 45 s
   * - Columbia Gorge (``Tests/Real_Terrain/Columbia_Gorge``; 200 m / 100 m, 4 ranks)
     - 0.125 / **0.122** m/s
     - 0.105 / **0.093** m/s
     - 276 exchanges, 260 s

On Askervein the joint run moves the outer domain towards the 20 m run outside the footprint in every field
(:math:`v` 9.37e-3 to 6.93e-3, :math:`k` 3.73e-3 to 3.12e-3 RMS) and leaves the nest's :math:`u` where it was. In the
nest :math:`w` improves (2.87e-3 to 1.48e-3) and :math:`v` and :math:`k` get worse (5.68e-3 to 6.03e-3 and 1.56e-3 to
2.76e-3). There the one-way nest is already close to the 20 m run: its 10 m speed at the masts of lines A and AA
differs by 0.04 % RMS, the joint nest's by 0.11 %; the outer domain's moves from 1.18 % to 1.05 %. No nest cell in
the lowest 200 m differs from the 20 m run by more than 0.25 m/s in either run (largest 0.23 and 0.14 m/s).

Why the nest's :math:`u` does not move: the joint run changes two things, the nest's levels (they
follow the outer domain's surfaces aloft) and the feedback. The first pass of the joint run (one way, with the joint
levels) separates them. The levels move by 0.10 m RMS (at most 1.3 m, largest at 130 m above the ground), and alone
raise the nest's RMS error in :math:`u` from 3.37e-3 to 3.92e-3, in :math:`v` from 5.68e-3 to 6.72e-3 and in :math:`k`
from 1.56e-3 to 2.80e-3 (:math:`w` 2.87e-3 to 2.84e-3). The feedback then lowers all four, to 3.37e-3, 6.03e-3,
1.48e-3 and 2.76e-3. :math:`u` ends where it began because the two cancel. The feedback acts near the faces: within 4
nest cells of them it lowers :math:`u` by 31 %, :math:`v` by 32 % and :math:`w` by 62 %; 8 to 40 cells in, :math:`u`
and :math:`v` rise by up to 14 %. The :math:`k` error added by the levels lies in the lowest three levels. All these
errors are at most 0.2 % of the fields' RMS (:math:`u` 8.8 m/s, :math:`v` 15.2 m/s, :math:`k` 1.66 m²/s²): the nest
and the 20 m run share the cell size, so the comparison measures the nest's faces and levels only. With
``nest.feedback_flux = nest`` (:ref:`Conservative faces <conservative-faces>`; 156 exchanges, as before) the nest's
errors fall to 2.98e-3 (:math:`u`), 4.74e-3 (:math:`v`), 1.21e-3 (:math:`w`) and 2.65e-3 (:math:`k`): below the
one-way run's for :math:`u`, :math:`v` and :math:`w`, while the :math:`k` error the levels add stays. Outside the
footprint the outer domain moves a little further towards the 20 m run (:math:`u` 3.19e-3, :math:`v` 6.73e-3,
:math:`w` 2.44e-3, :math:`k` 3.05e-3).

The joint run costs more than the one-way run: on Askervein the joint run takes 45 s after the one-way runs, and the
uniform 20 m run 85 s (4 ranks). Over the footprint the outer
:math:`w` follows the restricted contravariant flux, not a cell mean of the nest's :math:`w`; on Askervein's slopes
the two differ by 0.05 m/s RMS, so read :math:`w` there from the nest. On a gentler hill (60 m, half-width 250 m, the
hill of ``Nests_One_Process``) where the outer cells already resolve the flow, the nest's :math:`u` error moved from
1.36e-3 to 1.70e-3 m/s (``check_nests_two_way.py`` of ``Nests_Two_Way`` with ``terrain.hill_height = 60`` and
``terrain.hill_half_width = 250`` in a copy of its deck).

.. _two-way-relax:

Two-way nests by relaxation (``nest.feedback = relax``)
-------------------------------------------------------

``nest.feedback = relax`` runs the same joint run as ``defect`` (:ref:`two-way-nests`: the one-way runs, then one
exchange per iteration, the same coordinate surfaces, pseudo-time step, stopping rule, damping rule and refusals)
with another way up (``Source/NestRelax``). Each exchange, :math:`q^*` is formed as for ``defect``, and the outer
domain takes one iteration of

.. math::

   R_c(q_c) + c\,(q_c - q^*) = 0, \qquad c = w(s)\, \frac{f\, U}{\Delta x}, \qquad
   w(s) = \min\Big(1, \max\big(0, \frac{s}{M}\big)\Big),

for :math:`u, v, w, \theta` and :math:`k`, with

- :math:`f`: ``nest.feedback_factor`` (default 0.5; see The default, below);
- :math:`U`: the mean over the footprint's cells of the horizontal speed of :math:`I q_f`,
  :math:`\sqrt{u_c^2 + v_c^2}` with :math:`u_c, v_c` the means of the cell's two faces; one value per nest, set at
  each exchange and printed at the first;
- :math:`\Delta x`: the outer domain's smaller horizontal cell size;
- :math:`s`: a point's distance from the nearest lateral face of the footprint, in outer cells (a cell
  :math:`i - L + \tfrac12` from a face at :math:`L`, a face point :math:`i - L`); a nest that spans a periodic
  direction has no faces in it;
- :math:`M = 2`: the defect source's margin. :math:`c` is 0 on the footprint's faces and outside it and rises to its
  full value 2 cells in. Where footprints overlap the larger rate holds.

**Why this rate:** :math:`U / \Delta x` is the :math:`x`-advection part of an outer cell's diagonal, so at
:math:`f = 1` the relaxation holds the outer domain towards the nest about as strongly as its own advection along the
flow holds it to its neighbours. The full diagonal also has the :math:`y` and vertical advection and the diffusion;
near the ground the diffusion dominates it, and there the relaxation is weaker relative to the cell's own terms. At the fixed point :math:`q_c - q^* = -R_c(q_c)/c`: the outer domain under a nest follows it more
closely as :math:`f` grows, but not exactly at any finite :math:`f` (with ``defect``, :math:`q_c = I q_f` on
:math:`\chi_s` up to the ring's effect).

:math:`c` does not depend on :math:`q_c`, and it is added to the diagonal of each equation's implicit rows, so the
diagonal is the source's exact derivative. No row is fixed and continuity gets no source. The pressure correction
acts on every outer cell, but its coefficient there is :math:`1/d` with :math:`d = 1/\Delta\tau + c` (SIMPLEC, as
with canopy drag), so under a nest it is smaller by :math:`1/(1 + f\,w\,U\Delta\tau/\Delta x)`: about 0.1 at
:math:`f = 1` and under 0.01 at :math:`f = 16` for a step of 8 cells. The rows a boundary condition or the wall law replaces (the ground
and top :math:`w` faces, the first cells' wall-:math:`k` rows) replace the source too. The kernels have the form of
the nest relaxation zone's (``Source/NestFaces``), kept apart from it. The log gives, for the first ten exchanges and
then every tenth, the largest and RMS :math:`|q_c - q^*|` where :math:`c > 0`, per equation.

**The pressure under a nest.** The defect source carries the nest's pressure into the outer momentum rows through
:math:`R_c(q^*)`. The relaxation has no pressure term: the outer pressure under a nest is the outer domain's own, set
by its continuity equation with the relaxation in the momentum rows. Under the nest the relaxation force takes the place of
the outer domain's own momentum balance: the target has no nest pressure gradient or transport in it. As :math:`f`
grows the outer velocity under the nest is held at :math:`I q_f` whatever the outer pressure is (the pressure
correction's coefficient falls too, above), while the nest still takes its faces from the outer velocity. The flow
rate through the footprint is then set by a loop: the nest returns the flow rate it is given, and the relaxation
imposes it on the outer domain. A two-dimensional linear model of the coupling gives this loop a
gain of about :math:`cA/(cA + D)`, with :math:`A` the footprint's area and :math:`D` the outer domain's friction,
which tends to 1 as :math:`f` grows; the exchanges grow as :math:`1/(1 - \text{gain})`. On the canonical hill this lets both domains drift: the nest and the
outer domain under it grow faster than the fine run, more so at larger :math:`f` (Results).

**Results.** ``Tests/Canonical_RANS/Nests_Two_Way`` (``check_nests_two_way.py --relax``; 50 m / 25 m, 2 ranks), RMS of
:math:`u` in m/s: the outer domain minus the restricted nest where :math:`c > 0` (the log's last value), and the nest
and the outer domain outside the footprint against the uniform 25 m run:

.. list-table::
   :header-rows: 1

   * - Run
     - Exchanges
     - Outer − nest
     - Nest
     - Outer, outside the footprint
   * - One way
     -
     -
     - 8.23e-3
     - 1.44e-2
   * - ``defect``
     - 208
     -
     - **2.03e-3**
     - **4.89e-3**
   * - ``relax``, :math:`f` = 0.25
     - 150
     - 1.33e-2
     - 2.32e-3
     - 5.90e-3
   * - ``relax``, :math:`f` = 0.5 (default)
     - 158
     - 1.07e-2
     - 2.35e-3
     - 5.02e-3
   * - ``relax``, :math:`f` = 1
     - 160
     - 8.04e-3
     - 3.38e-3
     - 5.03e-3
   * - ``relax``, :math:`f` = 2
     - 159
     - 5.52e-3
     - 5.60e-3
     - 6.14e-3
   * - ``relax``, :math:`f` = 4
     - 177
     - 3.49e-3
     - 9.94e-3
     - 9.54e-3
   * - ``relax``, :math:`f` = 16
     - 508
     - 1.50e-3
     - 3.27e-2
     - 2.94e-2

The rows at :math:`f` = 0.25 and 2 are from #116; the others are from the code of #130.

The difference between the two domains under the nest falls with :math:`f`, from about :math:`f^{-0.3}` at small
:math:`f` to :math:`f^{-0.6}` at :math:`f \ge 1`, rather than :math:`1/f`: the nest it is measured against changes
with :math:`f` too, and :math:`c\,|q_c - q^*|` (the outer residual at the fixed point) grows sevenfold over the
table. Against one way, both domains are closer
to the fine run at every :math:`f \le 2`. As :math:`f` grows from 0.25 the nest moves away from the fine run, and the
outer domain does so above :math:`f = 1`; at :math:`f \ge 4` the nest is further from the fine run than one way.
The mean :math:`u` error over the nest is :math:`-9\times10^{-5}` m/s at
:math:`f = 0.25` and :math:`+2.0\times10^{-2}` m/s at :math:`f = 16`, largest in the lowest layers. Leaving
:math:`k` out of the source does not change this (nest 3.26e-2 at :math:`f = 16`).

On the Columbia Gorge (``Tests/Real_Terrain/Columbia_Gorge``, steep real terrain, 4 ranks), the RMS of :math:`u`
(m/s) against the uniform 100 m run:

.. list-table::
   :header-rows: 1

   * - Run
     - Exchanges
     - Nest
     - Outer outside
     - Footprint
     - Downstream
   * - one way
     -
     - 0.125
     - 0.105
     - 0.204
     - 0.152
   * - ``defect``
     - 276
     - 0.122
     - 0.093
     - 0.164
     - 0.117
   * - :math:`f = 0.5`
     - 174
     - 0.126
     - 0.092
     - 0.166
     - 0.116
   * - :math:`f = 1`
     - 180
     - 0.128
     - 0.092
     - 0.167
     - 0.115
   * - :math:`f = 4`
     - 439
     - 0.147
     - 0.099
     - 0.180
     - 0.133

At :math:`f` = 0.5 and 1 the outer domain, the footprint and the region downstream are as close to the fine run as
with ``defect``, in about two thirds of the exchanges; the nest's :math:`u` is slightly further from it than one
way (0.128 against 0.125). The footprint :math:`w` is 0.045 at :math:`f = 0.5` (``defect`` 0.068, one way 0.044).
At :math:`f = 4` every region moves away from the fine run and the exchanges grow, as on the hill.

**The default.** :math:`f = 1` is the rate argument above: the relaxation as strong as the outer cell's own
advection, and was the default until #130. A sweep of :math:`f` = 0.1, 0.25, 0.5, 0.75, 1 and 2
on four nests (the hill, the lee-face nest, the Columbia Gorge and Askervein, with the recommended setting;
``Tests/Canonical_RANS/Nests_Two_Way/README.md``, Choice of the relaxation factor) had every run steady. At
:math:`f = 0.5` the nest's :math:`u` was as close to the fine run as at 1 or closer on all four (0.35 against 0.51 of
the one-way error on the hill, 0.88 against 1.04 on the lee nest, 0.90 against 1.05 on Askervein, 0.99 against 1.01
on the Gorge), the outer :math:`u` as close or closer, and Askervein took 143 exchanges against 191 (the other
three differ by less than the run-to-run spread of about 5 %). The costs: on the hill the outer :math:`w` and
:math:`k` outside the footprint are 8 to 9 % further from the fine run than at 1. At 0.25 the errors are
similar on average, but on the hill the outer :math:`u` is 0.37 of the one-way error (0.32 at 0.5) and the outer
:math:`k` 14 % further than at 1; at 0.1 the hill's errors grow and the Gorge takes 212 exchanges. The default is 0.5; ``nest.feedback_factor = 1`` keeps the old answer (``rt_nests_relax`` runs both).
The relaxation's target has no nest pressure gradient (above); a target that carries one could move the best
factor.

**Limits.**

- At a finite factor the outer domain under a nest differs from the restricted nest by :math:`O(R_c/c)` (the table).
- Some fields stay further from the fine run than one way at every factor: on the lee-face nest :math:`w` and
  :math:`k` over the nest (1.13 and 1.05 of one way at :math:`f = 0.5`, 1.43 and 1.31 at 1), on Askervein :math:`k`
  over the nest (1.69; 1.77 with ``defect``).
- A large factor slows the joint run: at :math:`f = 16` it took 508 exchanges, against 160 at :math:`f = 1`. The
  cause is the coupling loop's gain tending to 1 (above), not stiff rows: :math:`c` on the diagonal damps each outer
  step more as :math:`f` grows.
- The pressure under a nest is the outer domain's own (above).
- Two nests that touch each have :math:`w = 0` on their shared face, so a strip about :math:`2M` outer cells wide
  between them is relaxed weakly, although the nests cover it (the defect source's region has the same gap).
- The limits of ``defect`` that concern the joint run (no probe rows after the one-way run,
  ``simplec_newton`` refused, :math:`\nu_t` relaxation) hold for ``relax`` too.

.. _two-way-relax-consistent:

Pressure-consistent relaxation (``nest.feedback = relax_consistent``)
---------------------------------------------------------------------

``nest.feedback = relax_consistent`` is ``relax`` (:ref:`two-way-relax`: the same joint run, rate, weight, margin,
kernels' form and refusals) with the outer residual of the target taken out of the source. The outer domain takes one
iteration of

.. math::

   R_c(q_c) + c\,(q_c - q^*) - w(s)\, R_c(q^*) = 0,

with :math:`c`, :math:`w(s)` and :math:`q^*` as for ``relax`` and :math:`R_c(q^*)` the outer residual of the
restricted nests, evaluated as for ``defect`` (the same function, with the defect and relaxation sources left out),
for :math:`u, v, w, \theta` and :math:`k`.
Its ``nest.feedback_factor`` defaults to 1, not to ``relax``'s 0.5: the sweep behind ``relax``'s default
(:ref:`two-way-relax`, The default) did not cover ``relax_consistent``.

**What it changes.** At :math:`q_c = q^*` the left side is :math:`(1 - w)\,R_c(q^*)`: zero where :math:`w = 1`.
Pointwise, the restricted nests then satisfy the outer momentum, :math:`\theta` and :math:`k` rows there at every
factor, pressure gradient included; with ``relax`` the fixed point is :math:`q_c - q^* = -R_c(q_c)/c` instead. On the
margin (:math:`0 < w < 1`) a fraction :math:`(1 - w)` of :math:`R_c(q^*)` remains, so :math:`q_c` differs from
:math:`q^*` on the margin and, much less, on the rows with :math:`w = 1` that read it. In a 1-D linear model
the largest :math:`|q_c - q^*|` where :math:`w = 1` is :math:`5\times10^{-6}` at :math:`f = 1`, against
:math:`1.4\times10^{-2}` with ``relax`` and :math:`1.7\times10^{-4}` with ``defect``. The outer pressure under a nest
mostly follows the nest's: in a 2-D Stokes model the range of :math:`p - p^*` over the inner box is 0.03 to 0.04 (model units),
against 0.18 with ``relax`` and 0.20 with no feedback. The rest is set at the margin and does
not decay inward (the pressure is elliptic and continuity has no source).

**Why a value of** ``nest.feedback``: the term changes what the joint run's fixed point is, as the choice between
``defect`` and ``relax`` does; it is not a parameter of the relaxation's strength. A separate switch would read with
``defect``, ``inject`` and ``none``, where it has no meaning, and would need its own refusals.
``nest.feedback_factor`` keeps its meaning (:math:`f` in :math:`c`), and ``nest.feedback_flux = nest`` is refused, as
with ``relax``.

**Implementation.** :math:`w R_c(q^*)` is a constant of the outer iteration: it is formed once per exchange, when
:math:`q^*` is set, and subtracted in each residual evaluation by its own kernel (``Source/NestRelax.cpp``). The
diagonal is that of ``relax`` (the term does not depend on :math:`q_c`), so the pseudo-time iteration's amplification
factor with :math:`q^*` held is unchanged. The exchange loop's gain is not: in the 1-D model its spectral radius is
0.12 against 0.05 with ``relax`` at :math:`f = 0.25`, 0.20 against 0.13 at :math:`f = 1` and 0.45 against 0.43 at
:math:`f = 16`. The loop is stable with the domains stepped in turn, as the joint run does (target, outer step, nest
faces from the new outer, nest steps); stepped together from the same old states it can diverge at low cell Péclet
numbers, so that order is required. The cost is one more outer residual evaluation per exchange. The weights are cached at the first
exchange; the footprints may not change during the run. A restart rebuilds the term from :math:`q^*`. The log's
start line says ``with the weighted residual of the target``.

**Results.** On the canonical hill (``check_nests_two_way.py --relax-consistent``; the pinned deck, 50 m / 25 m,
2 ranks), RMS of :math:`u` in m/s as in :ref:`two-way-relax`:

.. list-table::
   :header-rows: 1

   * - Run
     - Exchanges
     - Outer − nest
     - Nest
     - Outer, outside the footprint
   * - One way
     -
     -
     - 8.23e-3
     - 1.44e-2
   * - ``defect``
     - 208
     -
     - **2.03e-3**
     - **4.89e-3**
   * - ``relax``, :math:`f` = 1
     - 160
     - 8.04e-3
     - 3.38e-3
     - 5.03e-3
   * - ``relax_consistent``, :math:`f` = 1
     - 169
     - 1.32e-3
     - 4.24e-3
     - 5.18e-3
   * - ``relax_consistent``, :math:`f` = 4
     - 174
     - 1.22e-3
     - 1.07e-2
     - 1.01e-2

The outer domain under the nest is about six times closer to the restricted nest than with ``relax`` at the same
factor, as the fixed point predicts. Both domains are slightly further from the fine run than with ``relax`` at
:math:`f = 1` (mean nest :math:`u` error :math:`+2.3\times10^{-3}` m/s), and further again at :math:`f = 4`
(:math:`+6.5\times10^{-3}` m/s): the coupling loop of :ref:`two-way-relax` (the nest returns the flow rate it is
given) is not changed by the term, and holding the outer domain to the nest more exactly lets the pair drift together.

On the Columbia Gorge (``Tests/Real_Terrain/Columbia_Gorge``, 4 ranks) and Askervein
(``Tests/Real_Terrain/Askervein_Nest``, the ``two_way`` layout, 4 ranks), both with the recommended setting, RMS of
:math:`u` (m/s) against each case's uniform fine run (the nest's cells everywhere; the Gorge's reference, one-way
and ``defect`` runs are those of #130):

.. list-table::
   :header-rows: 1

   * - Case
     - Run
     - Exchanges
     - Outer − nest
     - Nest
     - Outer, outside the footprint
   * - Gorge
     - one way
     -
     -
     - 0.122
     - 0.102
   * - Gorge
     - ``defect``
     - 279
     -
     - **0.120**
     - **0.088**
   * - Gorge
     - ``relax``, :math:`f` = 1
     - 184
     - 4.25e-2
     - 0.123
     - 0.088
   * - Gorge
     - ``relax_consistent``, :math:`f` = 1
     - 198
     - **1.42e-2**
     - 0.121
     - 0.091
   * - Askervein
     - one way
     -
     -
     - 3.37e-3
     - 4.62e-3
   * - Askervein
     - ``defect``
     - 156
     -
     - **3.37e-3**
     - **3.21e-3**
   * - Askervein
     - ``relax``, :math:`f` = 1
     - 191
     - 6.48e-3
     - 3.53e-3
     - 3.75e-3
   * - Askervein
     - ``relax_consistent``, :math:`f` = 1
     - 235
     - **1.31e-3**
     - 4.02e-3
     - 3.74e-3

The outer domain follows the nest three times (Gorge) and five times (Askervein) more closely than with ``relax``.
Against the fine run the changes are small and of both signs: on the Gorge the nest is slightly closer (0.121
against 0.123) and the outer domain slightly further (0.091 against 0.088); on Askervein the nest is further
(4.02e-3 against 3.53e-3, 1.19 times one way) and the outer domain the same. The joint run takes 8 % (Gorge) and
23 % (Askervein) more exchanges than ``relax``, fewer than ``defect`` on the Gorge and more on Askervein.

**Limits.**

- The target satisfies the outer rows only where :math:`w = 1`; on the margin :math:`(1 - w) R_c(q^*)` remains,
  and the pressure under the nest keeps part of the outer domain's.
- The joint run's order of steps is required (above).
- It does not bring either domain consistently closer to a fine run than ``relax`` (hill, Gorge, Askervein);
  ``defect`` remains closest on all three.
- One more outer residual evaluation per exchange.
- The footprints are fixed for the run (an abort names the change).
- The limits of ``relax`` (two touching nests, :math:`\nu_t` relaxation, no probe rows after the one-way run,
  ``simplec_newton`` refused) hold too.

.. _two-way-inject:

Two-way nests by injection (``nest.feedback = inject``)
-------------------------------------------------------

``nest.feedback = inject`` runs the same joint run as ``defect`` (:ref:`two-way-nests`: the one-way runs, then one
exchange per iteration, the same coordinate surfaces, pseudo-time step, stopping rule, damping rule, refusals and
restart) with a third way up (``Source/NestInject``). Each exchange, :math:`q^*` is formed as for ``defect``, and the
outer state on the defect source's region is moved towards it:

.. math::

   q_c \leftarrow q_c + \omega\,(q^* - q_c) \quad \text{on } \chi_s,

for :math:`u, v, w, p, \theta` and :math:`k`, with

- :math:`q_c`: the outer state;
- :math:`q^*`: the restricted nests :math:`I q_f` on their footprints (with the restriction's pressure level, as for
  ``defect``), the outer state elsewhere;
- :math:`\omega`: ``nest.feedback_factor`` (default 0.5; :math:`0 < \omega \le 1`), the fraction of the way taken
  per exchange;
- :math:`\chi_s`: the defect source's region, the outer points at least 2 outer cells inside the footprint's lateral
  faces, all heights (cells, and faces from L to H + 1 along their normal, so the region's boundary faces are moved
  too). The ring between :math:`\chi_s` and the footprint's faces keeps the outer domain's own state and equations,
  as with ``defect`` and ``relax``: the nest's faces take their values from the outer domain there.

**The outer iteration holds** :math:`\chi_s` **as data.** The outer residual's rows (momentum, continuity,
:math:`\theta` and :math:`k`) are zero on :math:`\chi_s`, so the residual norms and the stopping rule see only the
points outside it. The SIMPLEC step resets the velocity on :math:`\chi_s` after the momentum predictor and every
field at the end of the iteration. The pressure correction leaves :math:`\chi_s` out of its Poisson solve, as it does
buildings' solid cells: a held cell's row is :math:`p' = 0`, no correction flux crosses a held face, and the mean
removed from the right-hand side is the :math:`J`-weighted mean over the other cells. The ring cells beside
:math:`\chi_s` therefore see the injected boundary fluxes as given (a Dirichlet boundary of the outer solve).
:math:`w` on the ground and top faces stays the boundary condition's (on the terrain it follows :math:`u`).

**The fixed point.** At steady :math:`q_c = q^*` on :math:`\chi_s` whatever :math:`\omega` is: :math:`\omega`
changes the path, not the answer. After a steady joint run the outer state is injected once in full
(:math:`\omega = 1`), so the outer output under a nest is :math:`I q_f` of the final nest. A run stopped before
steady (``steady.max_iterations``) writes the iterate, so its restart continues the straight run bit for bit.

**Why** :math:`\omega < 1`. With :math:`\omega = 1` the outer domain and the nest are two regions coupled by
Dirichlet data, overlapping by 2 outer cells, each taking one pseudo-time step per exchange. On the canonical hill
neither 1 nor 0.75 reached steady (table below). A partial move per exchange damps that coupling. Smaller :math:`\omega`
reaches the same state in more exchanges.

**Results.** ``Tests/Canonical_RANS/Nests_Two_Way`` (``check_nests_two_way.py --inject``; 50 m / 25 m, 2 ranks),
RMS of :math:`u` in m/s against the uniform 25 m run, over the nest and over the outer cells outside the footprint:

.. list-table::
   :header-rows: 1

   * - Run
     - Exchanges
     - Nest
     - Outer, outside the footprint
   * - One way
     -
     - 8.23e-3
     - 1.44e-2
   * - ``defect``
     - 208
     - **2.03e-3**
     - **4.89e-3**
   * - ``relax``, :math:`f` = 1
     - 160
     - 3.38e-3
     - 5.03e-3
   * - ``inject``, :math:`\omega` = 0.25
     - 602
     - 4.22e-3
     - 5.41e-3
   * - ``inject``, :math:`\omega` = 0.5
     - 344
     - 4.22e-3
     - 5.41e-3
   * - ``inject``, :math:`\omega` = 0.6
     - 300
     - 4.22e-3
     - 5.41e-3
   * - ``inject``, :math:`\omega` = 0.75
     - not steady after 3000
     -
     -
   * - ``inject``, :math:`\omega` = 1
     - not steady after 600
     -
     -

(0.25, 0.6, 0.75 and 1 are not in the test.) Up to :math:`\omega = 0.6` the three runs reach the same state to
the digits shown, in fewer exchanges as :math:`\omega` grows. At 0.75 the difference between the outer state and
the restricted nest falls for about 8 exchanges and then grows until it settles into an oscillation (nest-face
change about 1 m/s per exchange, RMS of :math:`v` under the nest about 6e-2 m/s); at 1 the oscillation is larger
(about 0.28 m/s). The default 0.5 is within a factor of 1.5 of that limit here; the limit depends on the case.

Both domains are closer to the fine run than one way: the nest's error is about half the one-way one and the outer
domain's about 0.38 of it. On this hill ``inject`` is further from the fine run than ``defect`` (nest 4.22e-3
against 2.03e-3, outer 5.41e-3 against 4.89e-3) and than ``relax`` at :math:`f = 1`, and takes more exchanges than
either.

**Terrain.** ``Tests/Real_Terrain/Columbia_Gorge`` (200 m outer, 100 m nest, 4 ranks, ``--relax 1 --inject 0.5``)
and ``Tests/Real_Terrain/Askervein_Nest`` (40 m / 20 m, ``two_way relax inject``), with the recommended setting.
RMS against the uniform fine run, in m/s; "outer" is the outer domain outside the footprint:

.. list-table::
   :header-rows: 1

   * - Case, run
     - Exchanges
     - Nest :math:`u`
     - Nest :math:`v`
     - Outer :math:`u`
     - Outer :math:`v`
   * - Gorge, one way
     -
     - 0.122
     - 0.033
     - 0.102
     - 0.034
   * - Gorge, ``defect``
     - 279
     - **0.120**
     - **0.030**
     - **0.088**
     - **0.031**
   * - Gorge, ``relax``, :math:`f` = 1
     - 184
     - 0.123
     - **0.030**
     - **0.088**
     - **0.031**
   * - Gorge, ``inject``, :math:`\omega` = 0.5
     - 555
     - 0.128
     - 0.049
     - 0.093
     - 0.036
   * - Askervein, one way
     -
     - 3.37e-3
     - **5.68e-3**
     - 4.62e-3
     - 9.37e-3
   * - Askervein, ``defect``
     - 156
     - **3.37e-3**
     - 6.03e-3
     - **3.21e-3**
     - **6.93e-3**
   * - Askervein, ``relax``, :math:`f` = 1
     - 191
     - 3.53e-3
     - 6.33e-3
     - 3.75e-3
     - 7.56e-3
   * - Askervein, ``inject``, :math:`\omega` = 0.5
     - 2153
     - 9.92e-3
     - 3.45e-2
     - 1.37e-2
     - 2.46e-2

Both terrain runs reach steady at the default :math:`\omega`. On both, ``inject`` is further from the fine run than
``defect`` and ``relax`` and takes the most exchanges. On the Gorge its outer domain is still closer to the fine run
than one way, its nest :math:`u` is 0.006 m/s further than one way, and its nest :math:`v` 0.016 m/s further. On
Askervein, where every run is within 0.035 m/s of the fine run, ``inject`` is the furthest of the four in
:math:`u` and :math:`v`, nest and outer; at the masts its 10 m speed is within 0.22 % of the fine run (``defect``
0.11 %), and the RMS of its speed-up against the observations, on either mast line, differs from the ``defect``
run's by at most 0.004. Counts of cells in the lowest 200 m more than 1 m/s from the Gorge's fine run: footprint 96
(``defect`` 111, ``relax`` 93, one way 134), nest 3 (6, 4, 3).

**Limits.**

- :math:`\omega` = 0.75 and 1 (plain injection) did not reach steady on the canonical hill; values above the
  default are accepted but may not converge, and a joint run that stops reports it (``steady.max_iterations``).
- With ``turb.nut_relax`` < 1 the closure's :math:`\nu_t` on :math:`\chi_s` blends with the outer domain's own
  history, as with ``defect`` and ``relax``; it enters the ring's diffusive fluxes on :math:`\chi_s`'s boundary faces.
- The pressure correction's masks need a GMRES back end of the Poisson solve, as buildings do.
- ``steady.stage = simplec_newton`` is refused, as with ``defect`` and ``relax``.
- Cost: an outer iteration with ``inject`` took 1.9 to 4.9 times the time of one with ``defect`` (hill 2.0,
  Askervein 1.9, Gorge 4.9; not profiled), and the joint runs took 2.2 (hill), 20 (Askervein) and 5.9
  (Gorge) times as long. The runs shared the machine with other work; the ratios are indicative.

**Buildings on terrain.** The masked Poisson solve used by buildings' solid cells (#86) shifts
:math:`\phi` by its mean after the solve on level ground only. On terrain the open faces beside masked cells read those cells'
:math:`\phi` in their cross terms, so a constant is not a null vector of the masked operator, and the shift would leave a
residual next to the masked cells. There the solve keeps GMRES's :math:`\phi` when cells are masked (test
``MaskedPoisson.EnclosedBlockOnABumpSolvesItsSystem``). The shift
applies only when the system is singular (no pressure boundary); every building deck in ``Tests`` has an outflow
face, and ``building_hill``, the one on terrain, gives the same plotfiles with and without the change.

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
at both faces on all boxes of a rank in one launch; on a domain no wider than the two ghost layers it is two
launches, the low face first. The values are those of one launch per field, face, box and layer.
Every copy of a ghost cell then holds one value. In its outermost lateral ghost layer a box lacks some of what the
rules read (the node heights have one lateral ghost layer, the wall law's gradient two), so after the fill each cell
outside the domain takes the value of the box whose column it is in, and the ghost images of cells inside the domain
are refreshed from the neighbouring boxes. No valid cell changes. The results then do not depend on the box
layout. Without the refresh they would where an inflow face crosses terrain at a box seam and ``turb.realizability`` is on: the
cap's strain in the first lateral ghost layer would read inflow face points computed from the clamped heights.
