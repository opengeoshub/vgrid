"""High-level API for position ↔ pixel index conversions."""

from __future__ import annotations

import math

from ..coordinates.projection import tu2za, za2tu
from ..coordinates.spherical import vec2za, za2vec
from ..pixel.fxy import fxy2tu, tu2fxy
from ..pixel.geometry import fxyCorners, fxySubpixel
from ..schemes.conversion import nest2ring, ring2nest
from ..schemes.nested import fxy2nest, nest2fxy
from ..types import AngularCoords, V3


def vec2PixNest(nside: int, v: V3) -> int:
    """Convert 3D vector to NESTED pixel index."""
    za = vec2za(v[0], v[1], v[2])
    tu = za2tu(za.z, za.a)
    fxy = tu2fxy(nside, tu.t, tu.u)
    return fxy2nest(nside, fxy.f, fxy.x, fxy.y)


def vec2PixRing(nside: int, v: V3) -> int:
    """Convert 3D vector to RING pixel index."""
    return nest2ring(nside, vec2PixNest(nside, v))


def ang2PixNest(nside: int, theta: float, phi: float) -> int:
    """Convert angular position to NESTED pixel index."""
    z = math.cos(theta)
    tu = za2tu(z, phi)
    fxy = tu2fxy(nside, tu.t, tu.u)
    return fxy2nest(nside, fxy.f, fxy.x, fxy.y)


def ang2PixRing(nside: int, theta: float, phi: float) -> int:
    """Convert angular position to RING pixel index."""
    return nest2ring(nside, ang2PixNest(nside, theta, phi))


def pix2VecNest(nside: int, ipix: int) -> V3:
    """Convert NESTED pixel index to 3D vector (pixel center)."""
    fxy = nest2fxy(nside, ipix)
    tu = fxy2tu(nside, fxy.f, fxy.x, fxy.y)
    za = tu2za(tu.t, tu.u)
    return za2vec(za.z, za.a)


def pix2VecRing(nside: int, ipix: int) -> V3:
    """Convert RING pixel index to 3D vector (pixel center)."""
    return pix2VecNest(nside, ring2nest(nside, ipix))


def pix2AngNest(nside: int, ipix: int) -> AngularCoords:
    """Convert NESTED pixel index to angular coordinates (pixel center)."""
    fxy = nest2fxy(nside, ipix)
    tu = fxy2tu(nside, fxy.f, fxy.x, fxy.y)
    za = tu2za(tu.t, tu.u)
    return AngularCoords(theta=math.acos(za.z), phi=za.a)


def pix2AngRing(nside: int, ipix: int) -> AngularCoords:
    """Convert RING pixel index to angular coordinates (pixel center)."""
    return pix2AngNest(nside, ring2nest(nside, ipix))


def cornersNest(nside: int, ipix: int) -> list[V3]:
    """Return the 4 corner vertices of a NESTED pixel."""
    fxy = nest2fxy(nside, ipix)
    return fxyCorners(nside, fxy.f, fxy.x, fxy.y)


def cornersRing(nside: int, ipix: int) -> list[V3]:
    """Return the 4 corner vertices of a RING pixel."""
    return cornersNest(nside, ring2nest(nside, ipix))


def pixcoord2VecNest(nside: int, ipix: int, ne: float, nw: float) -> V3:
    """Return 3D vector for a position within a NESTED pixel."""
    fxy = nest2fxy(nside, ipix)
    return fxySubpixel(nside, fxy.f, fxy.x, fxy.y, ne, nw)


def pixcoord2VecRing(nside: int, ipix: int, ne: float, nw: float) -> V3:
    """Return 3D vector for a position within a RING pixel."""
    return pixcoord2VecNest(nside, ring2nest(nside, ipix), ne, nw)
