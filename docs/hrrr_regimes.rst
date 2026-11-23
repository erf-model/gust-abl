Stratification class of an HRRR hour
====================================

In GUST runs at the AWAKEN site with an HRRR outer level, the difference between the GUST and HRRR winds changed
from hour to hour: five hours were run, and each was in a different part of the day (a night with a low-level jet,
a morning, an afternoon, an evening and a later night). ``tools/hrrr_regime.py`` puts each hour into one of nine
stratification classes, with three modifiers, so that the runs and their settings can be compared class by class
and over many hours.

The script reads HRRR only; it runs no GUST. Its classes and diagnostics need numpy. Downloading an hour also needs
`Herbie <https://github.com/blaylockbk/Herbie>`_, xarray and cfgrib.

Commands
--------

.. code-block:: bash

   # download hours (one-hour forecasts) and write a column for each
   python3 tools/hrrr_regime.py fetch 2023-05-10T05:00Z 2023-05-10T06:00Z --lat 36.362 --lon -97.405 \
       --save grib --out-dir columns --jobs 4
   # or make a column from a crop saved earlier (.npz with the same keys)
   python3 tools/hrrr_regime.py column crop.npz --valid 2023-05-10T05:00Z --lat 36.362 --lon -97.405 --out col.json
   # the class of one hour, with its neighbours for the transitions and the tendency
   python3 tools/hrrr_regime.py classify col.json --before prev.json --after next.json
   # many hours: neighbours are found by valid time; CSV of the diagnostics and counts per class and month
   python3 tools/hrrr_regime.py catalogue columns/*.json --jobs 4 --out catalogue.csv
   # the default thresholds, as JSON to edit and pass with --thresholds
   python3 tools/hrrr_regime.py thresholds

``fetch`` downloads byte-range subsets of the native-level (``wrfnat``) and surface (``wrfsfc``) files: wind,
height, temperature, humidity, pressure and cloud water on the first 20 hybrid levels (to about 6.4 km), and the
sensible and latent heat flux, :math:`u_*`, :math:`z_0`, the PBL height, the ground height, the 2 m temperature and
humidity and the surface pressure; about 190 MB an hour. ``--fxx`` (default 1, at least 1) is the forecast hour; the
analysis (f00) is refused, as for the outer level itself. ``--jobs`` downloads or classifies hours in that many
processes; ``catalogue`` gives the same results for any ``--jobs``.

``classify --relief H`` replaces the relief of the column (m) in the Froude number, for a GUST domain whose terrain
differs from HRRR's smoothed ground. ``classify`` prints the class, the modifiers, the rule that decided the class,
notes (for example that no neighbouring hour was given) and the diagnostics, as JSON; a value that is not finite is
written as the string ``"inf"``, ``"-inf"`` or ``"nan"``. ``catalogue`` prints a line per hour and the hours per
class and month, and ``--out`` writes the same as CSV (with the rule and a 0/1 jet column). An hour whose diagnostics
cannot be worked out is listed as refused and the exit status is 1. Input errors exit with status 2 and name the
input: a threshold, a file that cannot be read, a column with a missing or out-of-range field, two columns with the
same valid time. ``fetch`` and ``column`` check ``--lat`` (21 to 53), ``--lon`` (−135 to −60, or 225 to 300),
``--half-km`` (> 0) and ``--grad-half-km`` (at least ``--half-km``) before any download or read.

The column
----------

A column is the mean over a square of half-width ``--half-km`` (default 12.5 km, 81 HRRR cells) around the point:
the height of each hybrid level above HRRR's ground, the earth-relative wind (HRRR's grid-relative wind turned by
:math:`\sin(38.5^\circ)(\lambda + 97.5^\circ)`, :math:`\lambda` the longitude), the virtual potential temperature
:math:`\theta_v = T (10^5/p)^{0.2857} (1 + 0.61 q)`, and the largest cloud water mixing ratio in the square at each
level. The horizontal gradient of :math:`\theta_v` is a plane fit over a larger square of half-width
``--grad-half-km`` (default 40 km), turned from the sloping hybrid level to constant height. ``fetch`` crops
:math:`\sqrt{2}` times that half-width plus three cells, so the square fits in the crop at any turn of the HRRR
grid against north; a square that reaches the edge of the crop or of the HRRR domain is refused:

.. math::

   \nabla \theta_v \big|_z = \nabla \theta_v \big|_k - \frac{\partial \theta_v}{\partial z} \nabla z \big|_k ,

with :math:`\nabla \cdot |_k` the gradient along level :math:`k` and :math:`z` the level's height. The surface
fields are means over the small square; the relief :math:`H` is the range of HRRR's ground height in it.

Diagnostics
-----------

The virtual heat flux at the ground, from HRRR's sensible heat flux :math:`H_s` and latent heat flux :math:`E`
(both W/m², upward positive), and the Obukhov length:

.. math::

   q_v = \frac{H_s}{\rho c_p}(1 + 0.61 q) + 0.61\,\theta\,\frac{E}{\rho L_v}, \qquad
   L = -\frac{u_*^3\,\theta_v}{\kappa\, g\, q_v},

with :math:`\rho = p_s / (R_d T_v)` from the surface pressure :math:`p_s` and the 2 m virtual temperature
:math:`T_v = T (1 + 0.61 q)`, :math:`T` and :math:`q` the 2 m temperature and specific humidity, :math:`\theta`
and :math:`\theta_v` their potential temperatures, :math:`c_p` = 1005 J/(kg K), :math:`L_v` = 2.5e6 J/kg,
:math:`R_d` = 287.04 J/(kg K), :math:`g` = 9.81 m/s², :math:`\kappa` = 0.41 (the default of ``turb.kappa``) and
:math:`u_*` HRRR's friction velocity. :math:`z_i` is HRRR's PBL height. :math:`H_s` is taken as
:math:`\rho c_p \overline{w'\theta'}`; if HRRR's flux is :math:`\rho c_p \overline{w'T'}` instead, the factor
:math:`(10^5/p_s)^{0.2857}` left out is 0.7 % at AWAKEN and about 10 % at a 700 hPa site. :math:`L` varies as
:math:`u_*^3`, so 10 % in :math:`u_*` is 30 % in :math:`L`: for a weak flux the classes use :math:`\mathrm{Ri}_b`
instead (below). The other diagnostics:

.. math::

   \mathrm{Ri}_b = \frac{g}{\bar\theta_v} \frac{\Delta\theta_v\, \Delta z}{\Delta u^2 + \Delta v^2}, \qquad
   N^2 = \frac{g}{\bar\theta_v} \frac{\Delta\theta_v}{\Delta z}, \qquad
   \frac{\partial \mathbf{V}_g}{\partial z} = \frac{g}{f\,\theta_v}\, \mathbf{k} \times \nabla\theta_v \big|_z , \qquad
   \mathrm{Fr} = \frac{U}{N H}.

:math:`\mathrm{Ri}_b` is a layer value between the lowest level (about 10 m) and ``rib_top``, with :math:`\Delta u`
from the lowest level, not from no wind at the ground; a ground-referenced value is 2 to 6 times smaller at night.
:math:`N` is taken between :math:`\max(z_i, z_1)` and that height plus ``n_depth``, :math:`z_1` the lowest level;
it is a diagnostic only (no class uses it) and is NaN when the layer passes the column top. The geostrophic shear is the thermal wind (:math:`f` the Coriolis parameter,
:math:`\mathbf{k}` the vertical unit vector), averaged as a vector over :math:`z_i` to :math:`z_i` +
``thermal_depth``, where the temperature field is closer to balance than inside the boundary layer. In
:math:`\mathrm{Fr}`, :math:`U` and :math:`N` are taken over the lowest level to :math:`\max(H,` ``froude_min_depth``
:math:`, z_2)`, :math:`z_2` the second level; with no relief, or :math:`N^2 \le 0`, Fr is infinite.

The jet test takes the largest speed at or below ``jet_nose_max`` as the nose. The nose is not the lowest level.
The fall-off is the nose speed minus the smallest speed between the nose and ``jet_depth`` above it. A jet needs a
fall-off of at least ``jet_drop`` and at least ``jet_drop_frac`` of the nose speed.

The tendency is the change over an hour of the wind averaged over the ground (no wind) to ``tendency_top``: from
the hour before to the hour after when both are given (divided by 2 h), else from the one given.

Classes
-------

The first rule that holds gives the class.

.. list-table::
   :header-rows: 1
   :widths: 8 22 70

   * - Class
     - Regime
     - Rule
   * - C9
     - Cloud or moist
     - Cloud water above ``cloud_mixing_ratio`` anywhere in the square at or below ``cloud_top`` (the domain top).
       GUST is dry, so these hours are set aside.
   * - C7
     - Evening transition
     - :math:`q_v > 0` in the earlier hour and :math:`q_v \le 0` in the later hour, with :math:`|q_v| >`
       ``flux_eps`` in at least one of them, for the pair before–now or now–after.
   * - C8
     - Morning transition
     - The other way round.
   * - C1
     - Convective
     - :math:`q_v >` ``flux_eps``, :math:`-L \le` ``neutral_length`` and :math:`-z_i/L >` ``convective_ratio``.
   * - C2
     - Weakly unstable
     - :math:`q_v >` ``flux_eps``, :math:`-L \le` ``neutral_length`` and :math:`-z_i/L \le` ``convective_ratio``.
   * - C5
     - Low-level jet
     - :math:`q_v < 0`, :math:`\mathrm{Ri}_b > 0` and a jet, whatever :math:`L`. A night jet over a well-mixed
       surface layer, or over a surface where latent heat offsets most of the sensible cooling, has a large
       :math:`L` and would otherwise be neutral.
   * - C3, C6
     - Weak flux
     - :math:`|q_v| \le` ``flux_eps``: C6 when :math:`\mathrm{Ri}_b >` ``rib_crit``, else C3. :math:`L` is not
       used, since it is ill-conditioned when the flux is near zero.
   * - C3
     - Neutral
     - :math:`|L| >` ``neutral_length``.
   * - C6
     - Very stable
     - :math:`z_i/L >` ``very_stable_ratio`` or :math:`\mathrm{Ri}_b >` ``rib_crit``.
   * - C4
     - Weakly stable
     - The other stable hours.

Modifiers, added to any class: **Fr** when :math:`\mathrm{Fr} <` ``froude``, **B** (baroclinic) when the
geostrophic shear is above ``geo_shear``, **T** (unsteady) when the tendency is above ``tendency``. Without a
neighbouring hour there is no transition or T test, and the result says so.

Thresholds
----------

The defaults are starting values, not taken from the literature. They are to be set from a catalogue of many hours
and checked against published criteria. A JSON (or, with PyYAML, YAML) file given with ``--thresholds`` replaces
any of them; an unknown name, a value that is not a finite number or a value out of range stops the script and
names the threshold.

.. list-table::
   :header-rows: 1
   :widths: 22 12 16 50

   * - Name
     - Default
     - Range
     - Meaning
   * - ``flux_eps``
     - 0.005 K m/s
     - :math:`\ge 0`
     - :math:`|q_v|` up to this is a weak flux: the class comes from :math:`\mathrm{Ri}_b` (C3, C6), and a sign
       change between two such hours is not a transition (C7, C8); about 6 W/m²
   * - ``neutral_length``
     - 500 m
     - :math:`> 0`
     - :math:`|L|` above this is neutral
   * - ``convective_ratio``
     - 10
     - :math:`> 0`
     - :math:`-z_i/L` above this is convective
   * - ``very_stable_ratio``
     - 1
     - :math:`> 0`
     - :math:`z_i/L` above this is very stable
   * - ``rib_top``
     - 100 m
     - :math:`> 0`
     - Top of the layer for :math:`\mathrm{Ri}_b`
   * - ``rib_crit``
     - 0.25
     - :math:`> 0`
     - :math:`\mathrm{Ri}_b` above this is very stable
   * - ``n_depth``
     - 500 m
     - :math:`> 0`
     - Depth above :math:`z_i` for :math:`N` (a diagnostic only)
   * - ``jet_nose_max``
     - 1000 m
     - :math:`> 0`
     - Highest nose
   * - ``jet_depth``
     - 1000 m
     - :math:`> 0`
     - Depth above the nose searched for the fall-off
   * - ``jet_drop``
     - 2 m/s
     - :math:`> 0`
     - Least fall-off
   * - ``jet_drop_frac``
     - 0.1
     - (0, 1]
     - Least fall-off as a fraction of the nose speed
   * - ``cloud_top``
     - 3000 m
     - :math:`> 0`
     - Domain top: cloud above it is not counted
   * - ``cloud_mixing_ratio``
     - 1e-5 kg/kg
     - :math:`> 0`
     - Cloud water above this is a cloud. In a YAML file write ``1.0e-5``: YAML 1.1 reads ``1e-5`` as a string
   * - ``thermal_depth``
     - 1000 m
     - :math:`> 0`
     - Depth above :math:`z_i` for the geostrophic shear
   * - ``geo_shear``
     - 2 m/s per km
     - :math:`> 0`
     - Geostrophic shear above this sets B
   * - ``tendency_top``
     - 1000 m
     - :math:`> 0`
     - Top of the layer for the tendency
   * - ``tendency``
     - 1 m/s per h
     - :math:`> 0`
     - Tendency above this sets T
   * - ``froude_min_depth``
     - 100 m
     - :math:`> 0`
     - Least depth of the layer for :math:`U` and :math:`N` in Fr
   * - ``froude``
     - 1
     - :math:`> 0`
     - Fr below this sets Fr

A height a class or modifier needs that lies outside the column (for example :math:`z_i` + ``thermal_depth`` above the
20th level) stops the hour with a message naming the threshold; in a catalogue the hour is listed as refused.

AWAKEN hours
------------

``Tests/Infrastructure/HRRR_Regimes`` holds the columns of 16 hours of 10–11 May 2023 at the AWAKEN site
(36.362 N, 97.405 W) and pins their classes (CTest ``hrrr_regime_awaken``). With the default thresholds:

.. list-table::
   :header-rows: 1
   :widths: 18 10 72

   * - Hours (UTC)
     - Class
     - Note
   * - 10 May 04, 05, 06
     - C5
     - Jet nose 14–17 m/s at 455–866 m, fall-off 3.0–3.5 m/s; :math:`L` = 24.6 km, 770 and 361 m. At 04 UTC
       latent heat offsets most of the sensible cooling (:math:`q_v` = −4e-4 K m/s)
   * - 11
     - C6
     - :math:`z_i/L` = 4.7
   * - 12, 13
     - C8
     - The surface flux turns upward between them
   * - 19, 20, 21, 23
     - C1
     - :math:`-z_i/L` = 31 to 12
   * - 11 May 00, 01
     - C7
     - The surface flux turns downward between them
   * - 02
     - C4
     - :math:`z_i/L` = 0.84, :math:`\mathrm{Ri}_b` = 0.04
   * - 08, 09
     - C6
     - :math:`z_i/L` = 1.8
   * - 10
     - C5
     - A jet with 2.6 m/s fall-off; :math:`q_v` = −7e-4 K m/s and :math:`u_*` = 0.026 m/s

B and T are each set on 13 of the 16 hours with the starting thresholds. The thermal wind is 0.25 to 9.4 m/s per
km above :math:`z_i` and the tendency 0.8 to 3.5 m/s per hour. Some classes sit close to a threshold (the jet
fall-off at 11 May 09 UTC is 1.99 m/s against ``jet_drop`` = 2); the case README lists them. The pins are a
regression check of the script, not evidence that the thresholds are right.

Limits
------

- The classes are those of the box mean at one point. A domain with very different ground in it (a coast, a
  ridge) may hold several regimes in one hour.
- :math:`u_*` and :math:`z_i` are HRRR's; no surface-layer solve is made.
- The transitions are found from the surface flux alone, an hour apart. A transition that starts and ends within
  the hour is not seen.
- The thermal wind above :math:`z_i` changes by up to 9 m/s per km from one night hour to the next at AWAKEN, and
  its direction often differs from the box-mean shear over the same layer: in the residual layer the column
  gradient carries mesoscale structure as well as the synoptic baroclinicity. A ``geo_shear`` of 2 m/s per km is
  about 0.5 K per 100 km. The column does not store the slope correction, so its size cannot be read back.
- The thresholds are starting values; a catalogue of a year at a site needs about 190 MB of downloads an hour
  (about 1.7 TB a year), which this script does not reduce.

Data and acknowledgements
-------------------------

The High-Resolution Rapid Refresh (HRRR) is NOAA's model, developed at the Global Systems Laboratory and run
operationally by the National Centers for Environmental Prediction. The script reads its forecasts from the NOAA
Open Data Dissemination archive on AWS (https://registry.opendata.aws/noaa-hrrr-pds/) through Herbie, by Brian
Blaylock. We thank NOAA and the HRRR developers, Amazon Web Services and Brian Blaylock. Only the 16 small columns of
``Tests/Infrastructure/HRRR_Regimes`` are redistributed in this repository.
