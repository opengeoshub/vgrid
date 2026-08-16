"""Functions for converting between HEALPix resolution parameters."""

from __future__ import annotations

import math

from .constants import PI
from .utils import ilog2


def order2nside(order: int) -> int:
    """Convert HEALPix order (resolution level) to nside (grid divisions).

    nside = 2^order. Order must be in [0, 29].
    """
    if order < 0 or order > 29:
        raise ValueError(f"order must be between 0 and 29, got {order}")
    return 1 << order


def nside2order(nside: int) -> int:
    """Convert nside (grid divisions) back to order (resolution level)."""
    if nside <= 0 or not isinstance(nside, int):
        raise ValueError(f"nside must be a positive integer, got {nside}")
    if (nside & (nside - 1)) != 0:
        raise ValueError(f"nside must be a power of 2, got {nside}")
    return ilog2(nside)


def nside2npix(nside: int) -> int:
    """Total number of pixels at a given resolution: 12 * nside^2."""
    return 12 * nside * nside


def nside2pixarea(nside: int) -> float:
    """Solid angle of each pixel in steradians: π / (3 * nside^2)."""
    return PI / (3 * nside * nside)


def nside2resol(nside: int) -> float:
    """Approximate angular size of pixels in radians: sqrt(π/3) / nside."""
    return math.sqrt(PI / 3) / nside
