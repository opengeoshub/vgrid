"""Disc query functions for finding pixels within a circular region."""

from __future__ import annotations

import math
from collections.abc import Callable

from ..constants import PI, PI2, PI_2, PI_4
from ..coordinates.projection import tu2za, za2tu
from ..coordinates.spherical import angularDistance, vec2za
from ..lookup.lookup import pix2VecNest
from ..pixel.fxy import fxyEqual, rightNextPixel, tu2fxy
from ..pixel.geometry import maxPixelRadius
from ..schemes.conversion import nest2ring
from ..schemes.nested import fxy2nest
from ..types import V3
from ..utils import square, wrap


def queryDiscInclusiveNest(
    nside: int,
    v: V3,
    radius: float,
    cb: Callable[[int], None],
) -> None:
    """Find all NESTED pixels that overlap with a disc on the sphere.

    Inclusive means pixels are included if any part overlaps the disc.
    Radius must be < π/2.
    """
    if radius > PI_2:
        raise ValueError("query_disc: radius must be < PI/2")

    pixrad = maxPixelRadius(nside)
    d = PI_4 / nside

    za0 = vec2za(v[0], v[1], v[2])
    z0, a0 = za0.z, za0.a
    sin_t = math.sqrt(1 - z0 * z0)

    cos_r = math.cos(radius)
    sin_r = math.sin(radius)

    z1 = z0 * cos_r + sin_t * sin_r
    z2 = z0 * cos_r - sin_t * sin_r

    u1 = za2tu(z1, 0).u
    u2 = za2tu(z2, 0).u

    cover_north_pole = sin_t * cos_r - z0 * sin_r < 0
    cover_south_pole = sin_t * cos_r + z0 * sin_r < 0

    i1 = int(math.floor((PI_2 - u1) / d))
    i2 = int(math.floor((PI_2 - u2) / d + 1))

    if cover_north_pole:
        i1 += 1
        for i in range(1, i1 + 1):
            _walkRing(nside, i, cb)
        i1 += 1
    if i1 == 0:
        _walkRing(nside, 1, cb)
        i1 = 2

    if cover_south_pole:
        i2 -= 1
        for i in range(i2, 4 * nside):
            _walkRing(nside, i, cb)
        i2 -= 1
    if i2 == 4 * nside:
        _walkRing(nside, 4 * nside - 1, cb)
        i2 = 4 * nside - 2

    theta = math.acos(z0)
    for i in range(i1, i2 + 1):

        def _filter(ipix: int, _i: int = i) -> None:
            if angularDistance(pix2VecNest(nside, ipix), v) <= radius + pixrad:
                cb(ipix)

        _walkRingAround(nside, i, a0, theta, radius + pixrad, _filter)


def queryDiscInclusiveRing(
    nside: int,
    v: V3,
    radius: float,
    cb_ring: Callable[[int], None],
) -> None:
    """Find all RING pixels that overlap with a disc."""

    def _cb(ipix: int) -> None:
        cb_ring(nest2ring(nside, ipix))

    queryDiscInclusiveNest(nside, v, radius, _cb)


def _walkRing(nside: int, i: int, cb: Callable[[int], None]) -> None:
    """Iterate all pixels in a specific ring."""
    u = PI_4 * (2 - i / nside)
    t = PI_4 * (1 + (1 - (i % 2)) / nside)

    begin = tu2fxy(nside, t, u)
    s = begin

    while True:
        cb(fxy2nest(nside, s.f, s.x, s.y))
        s = rightNextPixel(nside, s)
        if fxyEqual(s, begin):
            break


def _walkRingAround(
    nside: int,
    i: int,
    a0: float,
    theta: float,
    r: float,
    cb: Callable[[int], None],
) -> None:
    """Iterate pixels in a ring near a specific azimuth."""
    if theta < r or theta + r > PI:
        _walkRing(nside, i, cb)
        return

    u = PI_4 * (2 - i / nside)
    z = tu2za(PI_4, u).z

    st = math.sin(theta)
    ct = math.cos(theta)
    sr = math.sin(r)
    cr = math.cos(r)

    w = math.atan2(
        math.sqrt(-square(z - ct * cr) / (square(st) * sr * sr) + 1) * sr,
        (-z * ct + cr) / st,
    )

    if w >= PI:
        _walkRing(nside, i, cb)
        return

    t1 = _centerT(nside, i, za2tu(z, wrap(a0 - w, PI2)).t)
    t2 = _centerT(nside, i, za2tu(z, wrap(a0 + w, PI2)).t)

    begin = tu2fxy(nside, t1, u)
    end = rightNextPixel(nside, tu2fxy(nside, t2, u))

    s = begin
    while not fxyEqual(s, end):
        cb(fxy2nest(nside, s.f, s.x, s.y))
        s = rightNextPixel(nside, s)


def _centerT(nside: int, i: int, t: float) -> float:
    """Snap a t-coordinate to the center of the nearest pixel in ring i."""
    d = PI_4 / nside
    t /= d
    # Match JS bitwise truncation toward zero for positive values
    t_i = int(t)
    t = (((t_i + (i % 2)) >> 1) << 1) + 1 - (i % 2)
    t *= d
    return t
