"""RING numbering scheme implementation."""

from __future__ import annotations

import math

from ..types import FXY


def fxy2ring(nside: int, f: int, x: int, y: int) -> int:
    """Convert pixel coordinates to RING pixel index."""
    f_row = f // 4
    f1 = f_row + 2
    v = x + y
    i = f1 * nside - v - 1

    if i < nside:
        f_col = f % 4
        return 2 * i * (i - 1) + i * f_col + nside - y - 1

    if i < 3 * nside:
        h = x - y
        f2 = 2 * (f % 4) - (f_row % 2) + 1
        k = (f2 * nside + h + 8 * nside) % (8 * nside)
        offset = 2 * nside * (nside - 1)
        return offset + (i - nside) * 4 * nside + (k >> 1)

    i_i = 4 * nside - i
    i_f_col = 3 - (f % 4)
    j = 4 * i_i - i_i * i_f_col - y
    i_j = 4 * i_i - j + 1
    return 12 * nside * nside - 2 * i_i * (i_i - 1) - i_j


def ring2fxy(nside: int, ipix: int) -> FXY:
    """Extract pixel coordinates from RING pixel index."""
    polar_lim = 2 * nside * (nside - 1)

    if ipix < polar_lim:
        i = int(math.floor((math.sqrt(1 + 2 * ipix) + 1) / 2))
        j = ipix - 2 * i * (i - 1)
        f = j // i
        k = j % i
        x = nside - i + k
        y = nside - 1 - k
        return FXY(f=f, x=x, y=y)

    if ipix < polar_lim + 8 * nside * nside:
        k = ipix - polar_lim
        ring = 4 * nside
        i = nside - (k // ring)
        s = 1 if i % 2 == 0 else 0
        j = 2 * (k % ring) + s

        jj = j - 4 * nside
        ii = i + 5 * nside - 1
        pp = (ii + jj) / 2
        qq = (ii - jj) / 2

        PP = int(math.floor(pp / nside))
        QQ = int(math.floor(qq / nside))
        V = 5 - (PP + QQ)
        H = PP - QQ + 4
        f = 4 * V + ((H >> 1) % 4)

        x = int(pp % nside)
        y = int(qq % nside)
        return FXY(f=f, x=x, y=y)

    p = 12 * nside * nside - ipix - 1
    i = int(math.floor((math.sqrt(1 + 2 * p) + 1) / 2))
    j = p - 2 * i * (i - 1)
    f = 11 - (j // i)
    k = j % i
    x = i - k - 1
    y = k
    return FXY(f=f, x=x, y=y)
