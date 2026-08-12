"""HEALPix spherical projection: (z, a) ↔ (t, u)."""

from __future__ import annotations

import math

from ..constants import PI, PI_2, PI_4, PI_8
from ..types import TU, V3, ZA
from ..utils import square
from .spherical import za2vec


def sigma(z: float) -> float:
    """Sigma function for the polar cap projection (equal-area)."""
    if z < 0:
        return -sigma(-z)
    return 2 - math.sqrt(3 * (1 - z))


def za2tu(z: float, a: float) -> TU:
    """HEALPix forward spherical projection: sphere → 2D plane."""
    if abs(z) <= 2 / 3:
        t = a
        u = 3 * PI_8 * z
        return TU(t=t, u=u)

    p_t = a % PI_2
    sigma_z = sigma(z)
    t = a - (abs(sigma_z) - 1) * (p_t - PI_4)
    u = PI_4 * sigma_z
    return TU(t=t, u=u)


def tu2za(t: float, u: float) -> ZA:
    """HEALPix inverse spherical projection: 2D plane → sphere."""
    abs_u = abs(u)

    if abs_u >= PI_2:
        return ZA(z=math.copysign(1.0, u), a=0.0)

    if abs_u <= PI_4:
        z = (8 / (3 * PI)) * u
        return ZA(z=z, a=t)

    t_t = t % PI_2
    a = t - ((abs_u - PI_4) / (abs_u - PI_2)) * (t_t - PI_4)
    z = math.copysign(1.0, u) * (1 - (1 / 3) * square(2 - (4 * abs_u) / PI))
    return ZA(z=z, a=a)


def tu2vec(t: float, u: float) -> V3:
    """Convert projection coordinates directly to 3D vector."""
    za = tu2za(t, u)
    return za2vec(za.z, za.a)
