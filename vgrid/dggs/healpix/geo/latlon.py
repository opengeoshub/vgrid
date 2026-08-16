"""Geographic coordinate conversions (latitude/longitude)."""

from __future__ import annotations

from ..constants import PI2, PI_2
from ..coordinates.spherical import vec2ang
from ..lookup.lookup import (
    ang2PixNest,
    ang2PixRing,
    cornersNest,
    cornersRing,
    pix2AngNest,
    pix2AngRing,
)
from ..types import AngularCoords, LonLat, V3
from ..utils import deg2Rad, rad2Deg


def pix2LonLatNest(nside: int, ipix: int) -> LonLat:
    """Convert a NESTED pixel index to latitude/longitude (pixel center)."""
    ang = pix2AngNest(nside, ipix)
    return _ang2LonLat(ang.theta, ang.phi)


def pix2LonLatRing(nside: int, ipix: int) -> LonLat:
    """Convert a RING pixel index to latitude/longitude (pixel center)."""
    ang = pix2AngRing(nside, ipix)
    return _ang2LonLat(ang.theta, ang.phi)


def lonLat2PixNest(nside: int, lon: float, lat: float) -> int:
    """Convert latitude/longitude to NESTED pixel index."""
    ang = _lonLat2Ang(lon, lat)
    return ang2PixNest(nside, ang.theta, ang.phi)


def lonLat2PixRing(nside: int, lon: float, lat: float) -> int:
    """Convert latitude/longitude to RING pixel index."""
    ang = _lonLat2Ang(lon, lat)
    return ang2PixRing(nside, ang.theta, ang.phi)


def cornersNestLonLat(nside: int, ipix: int) -> list[LonLat]:
    """Get the 4 corner coordinates of a NESTED pixel in lon/lat."""
    corners = cornersNest(nside, ipix)
    center_lon, _ = pix2LonLatNest(nside, ipix)
    return _unwrapCornerLons([vec2LonLat(v) for v in corners], center_lon)


def cornersRingLonLat(nside: int, ipix: int) -> list[LonLat]:
    """Get the 4 corner coordinates of a RING pixel in lon/lat."""
    corners = cornersRing(nside, ipix)
    center_lon, _ = pix2LonLatRing(nside, ipix)
    return _unwrapCornerLons([vec2LonLat(v) for v in corners], center_lon)


def _unwrapCornerLons(corners: list[LonLat], ref_lon: float) -> list[LonLat]:
    """Unwrap corner longitudes so they form a continuous polygon."""
    result: list[LonLat] = []
    for lon, lat in corners:
        if abs(lat) >= 90 - 1e-6:
            lon = ref_lon
        while lon - ref_lon > 180:
            lon -= 360
        while lon - ref_lon < -180:
            lon += 360
        result.append((lon, lat))
    return result


def vec2LonLat(v: V3) -> LonLat:
    """Convert 3D vector to latitude/longitude."""
    ang = vec2ang(v)
    return _ang2LonLat(ang.theta, ang.phi)


def _ang2LonLat(theta: float, phi: float) -> LonLat:
    lat = rad2Deg(PI_2 - theta)
    lon = rad2Deg(phi)
    while lon > 180:
        lon -= 360
    return (lon, lat)


def _lonLat2Ang(lon: float, lat: float) -> AngularCoords:
    theta = deg2Rad(90 - lat)
    phi = deg2Rad(lon)
    while phi < 0:
        phi += PI2
    return AngularCoords(theta=theta, phi=phi)
