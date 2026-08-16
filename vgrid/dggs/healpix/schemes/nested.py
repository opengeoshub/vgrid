"""NESTED numbering scheme implementation."""

from __future__ import annotations

from ..types import FXY
from ..utils import assert_

BITS = 26
MAX_COORD = (1 << BITS) - 1


def spread1By1(n: int) -> int:
    """Spread 32-bit value: insert 0 between each bit."""
    n = n & 0xFFFFFFFF
    n = (n | (n << 16)) & 0x0000FFFF0000FFFF
    n = (n | (n << 8)) & 0x00FF00FF00FF00FF
    n = (n | (n << 4)) & 0x0F0F0F0F0F0F0F0F
    n = (n | (n << 2)) & 0x3333333333333333
    n = (n | (n << 1)) & 0x5555555555555555
    return n


def compact1By1(n: int) -> int:
    """Compact 64-bit value: extract even bits."""
    n = n & 0x5555555555555555
    n = (n | (n >> 1)) & 0x3333333333333333
    n = (n | (n >> 2)) & 0x0F0F0F0F0F0F0F0F
    n = (n | (n >> 4)) & 0x00FF00FF00FF00FF
    n = (n | (n >> 8)) & 0x0000FFFF0000FFFF
    n = (n | (n >> 16)) & 0x00000000FFFFFFFF
    return n


def bitCombine(x: int, y: int) -> int:
    """Interleave bits of x and y (Morton code / Z-order curve)."""
    assert_(x <= MAX_COORD and x >= 0, "x must fit in 26 bits")
    assert_(y <= MAX_COORD and y >= 0, "y must fit in 26 bits")
    return spread1By1(x) | (spread1By1(y) << 1)


def bitDecombine(p: int) -> tuple[int, int]:
    """De-interleave bits to recover x and y coordinates."""
    return compact1By1(p), compact1By1(p >> 1)


def fxy2nest(nside: int, f: int, x: int, y: int) -> int:
    """Convert pixel coordinates to NESTED pixel index."""
    return f * nside * nside + bitCombine(x, y)


def nest2fxy(nside: int, ipix: int) -> FXY:
    """Extract pixel coordinates from NESTED pixel index."""
    nside2 = nside * nside
    f = ipix // nside2
    k = ipix % nside2
    x, y = bitDecombine(k)
    return FXY(f=f, x=x, y=y)
