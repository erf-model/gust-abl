"""Wind helpers: directions, the faces the wind enters and leaves through, and
the inflow profile's parameters.

Two direction conventions meet here.

* The case description uses the meteorological one: the direction the wind
  comes FROM, in degrees clockwise from north (270 = from the west).
* GUST's ``bc.inflow_direction`` is the direction the wind blows TOWARDS, in
  degrees counter-clockwise from the x axis (docs/howto.rst).

With x east and y north, towards = 270 - from (mod 360).

The log law is GUST's inflow profile (docs/boundary_conditions.rst),

    U(z) = u* / kappa  ln((z + z0) / z0),

z the height above the ground, so u* = kappa U_ref / ln((z_ref + z0) / z0).
"""
import math

KAPPA = 0.41          # turb.kappa's default
EARTH_ROTATION = 2.0 * math.pi / 86400.0   # rad/s: GUST's f uses forcing.rotation_period (default 86400 s)
DRAG_LAW_A = 1.8      # the neutral geostrophic drag law's constants (Troen and Petersen, 1989,
DRAG_LAW_B = 4.5      # European Wind Atlas, section 8.1)
ALIGNED_TOL = 1.0e-9  # |cos| or |sin| below this: the wind is along an axis


def towards_from_meteorological(direction_from):
    """``bc.inflow_direction`` (degrees, in [0, 360)) of a meteorological direction."""
    return (270.0 - float(direction_from)) % 360.0


def meteorological_from_towards(towards):
    """The inverse of :func:`towards_from_meteorological`."""
    return (270.0 - float(towards)) % 360.0


def unit_vector(towards):
    """(cos, sin) of the flow direction, exactly 0 or +-1 along an axis."""
    c = math.cos(math.radians(towards))
    s = math.sin(math.radians(towards))
    if abs(c) < ALIGNED_TOL:
        c, s = 0.0, math.copysign(1.0, s)
    if abs(s) < ALIGNED_TOL:
        s, c = 0.0, math.copysign(1.0, c)
    return c, s


def lateral_faces(towards, side_faces="slip_wall"):
    """The four lateral boundary types for a flow direction.

    A face the wind enters through is ``inflow`` and one it leaves through is
    ``outflow``. When the wind is along an axis, the two faces parallel to it
    take ``side_faces`` (``slip_wall`` or ``periodic``).
    Returns a dict with the keys xlo, xhi, ylo, yhi.
    """
    if side_faces not in ("slip_wall", "periodic"):
        raise ValueError("side_faces must be slip_wall or periodic")
    c, s = unit_vector(towards)
    faces = {}
    for lo, hi, comp in (("xlo", "xhi", c), ("ylo", "yhi", s)):
        if comp > 0.0:
            faces[lo], faces[hi] = "inflow", "outflow"
        elif comp < 0.0:
            faces[lo], faces[hi] = "outflow", "inflow"
        else:
            faces[lo] = faces[hi] = side_faces
    return faces


def friction_velocity(speed, height, z0, kappa=KAPPA):
    """u* of the log law that gives ``speed`` at ``height`` above the ground."""
    return kappa * speed / math.log((height + z0) / z0)


def log_law_speed(ustar, height, z0, kappa=KAPPA):
    """The log law's speed at a height above the ground."""
    return ustar / kappa * math.log((height + z0) / z0)


def coriolis_parameter(latitude, rotation_rate=EARTH_ROTATION):
    """f = 2 Omega sin(latitude), 1/s (negative in the southern hemisphere)."""
    return 2.0 * rotation_rate * math.sin(math.radians(latitude))


def geostrophic_speed(ustar, z0, latitude, kappa=KAPPA):
    """The neutral geostrophic drag law's wind speed above a surface of roughness z0,

        G = u* / kappa  sqrt((ln(u* / (|f| z0)) - A)^2 + B^2),   A = 1.8, B = 4.5.

    It is the first value the column run starts from when the case gives a wind at
    a height: GUST's column then turns and scales it until that wind is met.
    """
    f = abs(coriolis_parameter(latitude))
    if not f > 0.0:
        raise ValueError("the drag law needs a latitude away from the equator")
    return ustar / kappa * math.hypot(math.log(ustar / (f * z0)) - DRAG_LAW_A, DRAG_LAW_B)


def sweep_directions(spec):
    """The directions of a sweep: a list, or start/stop/step (stop included)."""
    if isinstance(spec, dict):
        start, stop, step = float(spec["start"]), float(spec["stop"]), float(spec["step"])
        count = int(math.floor((stop - start) / step + 1.0e-9)) + 1
        return [start + i * step for i in range(count)]
    return [float(d) for d in spec]


def direction_label(direction_from):
    """A folder-name label to 0.001 degree: 222 -> '222', 22.5 -> '022p5', 359.9999 -> '000'."""
    thousandths = int(round(float(direction_from) * 1000.0)) % 360000
    whole, frac = divmod(thousandths, 1000)
    if frac == 0:
        return "%03d" % whole
    return "%03dp%s" % (whole, ("%03d" % frac).rstrip("0"))
