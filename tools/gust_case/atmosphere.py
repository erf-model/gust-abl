"""Atmosphere helpers: the reference sounding and the Rayleigh layer.

The sounding of an inversion (amrTerrain's inversionHeight, inversionWidth,
inversionStrength and lapseRate) is piecewise linear:

    theta = T0                                 0 <= z <= base
    theta = T0 + strength (z - base)/thickness base <= z <= base + thickness
    theta = T0 + strength + lapse (z - base - thickness)   above

GUST interpolates ``bc.theta_profile`` linearly between its heights, so the
four corner points give this profile exactly.
"""
from .errors import CaseError


def inversion_sounding(reference_temperature, top, base=None, thickness=None, strength=None, lapse_rate=None):
    """Heights and theta of the corner points, from the ground to ``top``."""
    t0 = float(reference_temperature)
    lapse = float(lapse_rate or 0.0)
    if base is None:
        if lapse == 0.0:
            return [0.0], [t0]
        return [0.0, top], [t0, t0 + lapse * top]
    thickness = float(thickness)
    strength = float(strength)
    if base + thickness >= top:
        raise CaseError("atmosphere.inversion: base + thickness = %g m reaches the top of the domain (domain.top = %g m)"
                        % (base + thickness, top))
    heights = [0.0, base, base + thickness, top]
    theta = [t0, t0, t0 + strength, t0 + strength + lapse * (top - base - thickness)]
    return heights, theta


def check_sounding(heights, theta, top):
    """Refuse a sounding GUST would refuse: heights from 0, strictly increasing, equal lengths."""
    if len(heights) != len(theta):
        raise CaseError("atmosphere.sounding: heights has %d values and theta %d; they must be equal in number"
                        % (len(heights), len(theta)))
    if heights[0] != 0.0:
        raise CaseError("atmosphere.sounding.heights: the first height must be 0 (it is %g)" % heights[0])
    for a, b in zip(heights[:-1], heights[1:]):
        if b <= a:
            raise CaseError("atmosphere.sounding.heights: must increase strictly (%g then %g)" % (a, b))
    if heights[-1] < top:
        raise CaseError("atmosphere.sounding.heights: the last height, %g m, is below the top of the domain (%g m)"
                        % (heights[-1], top))


def is_uniform(theta):
    return all(t == theta[0] for t in theta)


def never_falls(theta):
    return all(b >= a for a, b in zip(theta[:-1], theta[1:]))


def gradient_at_top(heights, theta):
    """d theta / dz of the sounding's last segment (0 for a single point)."""
    if len(heights) < 2:
        return 0.0
    return (theta[-1] - theta[-2]) / (heights[-1] - heights[-2])
