"""Conversions between 3D vectors and spherical coordinates."""

from __future__ import annotations

import math

from ..constants import PI2
from ..types import AngularCoords, V3, ZA


def vec2za(X: float, Y: float, z: float) -> ZA:
    """Convert 3D Cartesian coordinates to normalized spherical (z, a)."""
    r2 = X * X + Y * Y

    if r2 == 0:
        return ZA(z=-1.0 if z < 0 else 1.0, a=0.0)

    a = (math.atan2(Y, X) + PI2) % PI2
    z = z / math.sqrt(z * z + r2)
    return ZA(z=z, a=a)


def za2vec(z: float, a: float) -> V3:
    """Convert spherical (z, a) to 3D Cartesian unit vector."""
    sin_theta = math.sqrt(1 - z * z)
    X = sin_theta * math.cos(a)
    Y = sin_theta * math.sin(a)
    return (X, Y, z)


def ang2vec(theta: float, phi: float) -> V3:
    """Convert angular (theta, phi) to 3D unit vector."""
    z = math.cos(theta)
    return za2vec(z, phi)


def vec2ang(v: V3) -> AngularCoords:
    """Convert 3D unit vector to angular (theta, phi)."""
    za = vec2za(v[0], v[1], v[2])
    return AngularCoords(theta=math.acos(za.z), phi=za.a)


def angularDistance(a: V3, b: V3) -> float:
    """Angular distance between two unit vectors in radians."""
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    dist2 = dx * dx + dy * dy + dz * dz
    return 2 * math.asin(math.sqrt(dist2) / 2)
