"""Pixel geometry functions: corners, sub-pixel positions, pixel radius."""

from __future__ import annotations

from ..constants import PI_4
from ..coordinates.projection import tu2vec, tu2za
from ..coordinates.spherical import angularDistance, za2vec
from ..types import V3
from .fxy import fxy2tu


def fxyCorners(nside: int, f: int, x: int, y: int) -> list[V3]:
    """Return the 4 corner vertices: [north, west, south, east]."""
    tu = fxy2tu(nside, f, x, y)
    d = PI_4 / nside

    corners: list[V3] = []
    for tt, uu in ((0, d), (-d, 0), (0, -d), (d, 0)):
        za = tu2za(tu.t + tt, tu.u + uu)
        corners.append(za2vec(za.z, za.a))
    return corners


def fxySubpixel(
    nside: int, f: int, x: int, y: int, ne: float, nw: float
) -> V3:
    """Return 3D vector for a sub-pixel position within (f, x, y)."""
    tu = fxy2tu(nside, f, x, y)
    d = PI_4 / nside
    za = tu2za(tu.t + d * (ne - nw), tu.u + d * (ne + nw - 1))
    return za2vec(za.z, za.a)


def maxPixelRadius(nside: int) -> float:
    """Maximum angular radius of a pixel (center to farthest corner)."""
    unit = PI_4 / nside
    return angularDistance(
        tu2vec(unit, nside * unit),
        tu2vec(unit, (nside + 1) * unit),
    )
