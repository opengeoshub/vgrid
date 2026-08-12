"""Conversions between projection coordinates (t, u) and pixel coordinates (f, x, y)."""

from __future__ import annotations

import math

from ..constants import PI_2, PI_4
from ..types import FXY, TU
from ..utils import clip, wrap


def tu2fpq(t: float, u: float) -> tuple[int, float, float]:
    """Convert projection coordinates to base pixel and fractional position.

    Returns (f, p, q) where f is base pixel and (p, q) are fractional coords [0, 1).
    """
    t /= PI_4
    u /= PI_4

    t = wrap(t, 8)
    t += -4
    u += 5

    pp = clip((u + t) / 2, 0, 5)
    PP = int(math.floor(pp))

    qq = clip((u - t) / 2, 3 - PP, 6 - PP)
    QQ = int(math.floor(qq))

    V = 5 - (PP + QQ)

    if V < 0:
        return 0, 1.0, 1.0

    H = PP - QQ + 4
    f = 4 * V + ((H >> 1) % 4)

    p = pp % 1
    q = qq % 1
    return f, p, q


def tu2fxy(nside: int, t: float, u: float) -> FXY:
    """Convert projection coordinates to base pixel and integer pixel indices."""
    f, p, q = tu2fpq(t, u)
    x = int(clip(math.floor(nside * p), 0, nside - 1))
    y = int(clip(math.floor(nside * q), 0, nside - 1))
    return FXY(f=f, x=x, y=y)


def fxy2tu(nside: int, f: int, x: int, y: int) -> TU:
    """Convert pixel coordinates back to projection coordinates (pixel center)."""
    f_row = f // 4
    f1 = f_row + 2
    f2 = 2 * (f % 4) - (f_row % 2) + 1

    v = x + y
    h = x - y

    i = f1 * nside - v - 1
    k = f2 * nside + h + 8 * nside

    t = (k / nside) * PI_4
    u = PI_2 - (i / nside) * PI_4
    return TU(t=t, u=u)


def fxyEqual(a: FXY, b: FXY) -> bool:
    """Compare two FXY pixel coordinates for equality."""
    return a.x == b.x and a.y == b.y and a.f == b.f


def rightNextPixel(nside: int, fxy: FXY) -> FXY:
    """Return the next pixel to the right (east) along a ring."""
    f, x, y = fxy.f, fxy.x, fxy.y

    x += 1

    if x == nside:
        row = f // 4
        if row == 0:
            f = (f + 1) % 4
            x = y
            y = nside
        elif row == 1:
            f = f - 4
            x = 0
        elif row == 2:
            f = 4 + ((f + 1) % 4)
            x = 0

    y -= 1

    if y == -1:
        row = f // 4
        if row == 0:
            f = 4 + ((f + 1) % 4)
            y = nside - 1
        elif row == 1:
            f = f + 4
            y = nside - 1
        elif row == 2:
            f = 8 + ((f + 1) % 4)
            y = x - 1
            x = 0

    return FXY(f=f, x=x, y=y)
