"""Conversion between NESTED and RING schemes."""

from __future__ import annotations

from .nested import fxy2nest, nest2fxy
from .ring import fxy2ring, ring2fxy


def nest2ring(nside: int, ipix: int) -> int:
    """Convert NESTED pixel index to RING pixel index."""
    fxy = nest2fxy(nside, ipix)
    return fxy2ring(nside, fxy.f, fxy.x, fxy.y)


def ring2nest(nside: int, ipix: int) -> int:
    """Convert RING pixel index to NESTED pixel index."""
    if nside == 1:
        return ipix
    fxy = ring2fxy(nside, ipix)
    return fxy2nest(nside, fxy.f, fxy.x, fxy.y)
