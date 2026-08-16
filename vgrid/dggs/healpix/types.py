"""Core type definitions for HEALPix."""

from __future__ import annotations

from typing import NamedTuple

# 3D vector on the unit sphere: [X, Y, Z]
# X = sin(theta)*cos(phi), Y = sin(theta)*sin(phi), Z = cos(theta)
V3 = tuple[float, float, float]

# Geographic coordinates: [lon, lat] in degrees
LonLat = tuple[float, float]

# Bounding box: [minLon, minLat, maxLon, maxLat] in degrees
BBox = tuple[float, float, float, float]


class FXY(NamedTuple):
    """Pixel coordinates within the HEALPix grid.

    - f: Base pixel index {0..11}
    - x: Index along the north-east direction [0, nside)
    - y: Index along the north-west direction [0, nside)
    """

    f: int
    x: int
    y: int


class ZA(NamedTuple):
    """Spherical coordinates in (z, a) form.

    - z: cos(colatitude) in [-1, 1]
    - a: azimuth angle in [0, 2π)
    """

    z: float
    a: float


class AngularCoords(NamedTuple):
    """Angular coordinates (theta, phi).

    - theta: Colatitude [0, π]
    - phi: Longitude [0, 2π)
    """

    theta: float
    phi: float


class TU(NamedTuple):
    """HEALPix projection coordinates.

    - t: longitude-like coordinate [0, 2π)
    - u: latitude-like coordinate [-π/2, π/2]
    """

    t: float
    u: float


class OrderPix(NamedTuple):
    """Order and nested pixel index pair (UNIQ scheme)."""

    order: int
    ipix: int
