Forest canopy
=============

A forest slows the wind inside and above it. GUST represents a forest as a drag field: each cell holds the product
of a drag coefficient and the leaf area density of the trees in it. This page describes how that field is built
from the ``canopy`` inputs (:doc:`inputs`).

**Status.** The drag field is built, checked and written to the plotfile. It does not act on the flow yet: the drag
in the momentum equations and the canopy terms of the turbulence model follow in later changes (Phase 8 of the
development tracker). Until then a forest can still be represented by its roughness alone, with
``bc.surface_z0_file``.

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
is accepted; editing a comment in it counts as a change.

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
