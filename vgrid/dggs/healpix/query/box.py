"""Bounding box queries and conversions between HEALPix cells and geographic boxes."""

from __future__ import annotations

import math
from collections.abc import Callable

from ..constants import PI, PI_2, PI_4
from ..coordinates.projection import tu2za, za2tu
from ..geo.latlon import cornersNestLonLat, pix2LonLatNest
from ..lookup.lookup import cornersNest, pix2VecNest
from ..pixel.fxy import fxy2tu, fxyEqual, rightNextPixel, tu2fxy
from ..pixel.geometry import maxPixelRadius
from ..schemes.conversion import nest2ring, ring2nest
from ..schemes.nested import fxy2nest
from ..schemes.ring import ring2fxy
from ..types import BBox
from ..utils import deg2Rad, rad2Deg


def queryBoxInclusiveNest(nside: int, bbox: BBox) -> list[int]:
    """Find all NESTED pixels that intersect a lat/lon bounding box."""
    return _queryBoxInclusiveImpl(nside, bbox, as_ring=False)


def queryBoxInclusiveRing(nside: int, bbox: BBox) -> list[int]:
    """Find all RING pixels that intersect a lat/lon bounding box."""
    return _queryBoxInclusiveImpl(nside, bbox, as_ring=True)


def cellsToBoundingBox(
    nside: int,
    cells: list[int],
    numberingScheme: str,
) -> BBox:
    """Compute the tightest lon/lat bounding box containing HEALPix cells.

    numberingScheme must be 'nest' or 'ring'.
    """
    if not cells:
        raise ValueError("cellsToBoundingBox: cells must not be empty")

    d = PI_4 / nside

    max_u = -math.inf
    min_u = math.inf
    max_u_ring_cell = -1
    min_u_ring_cell = -1

    az_values: list[float] = []
    az_cell_refs: list[int] = []

    for cell in cells:
        ring_cell = cell if numberingScheme == "ring" else nest2ring(nside, cell)
        fxy = ring2fxy(nside, ring_cell)
        tu = fxy2tu(nside, fxy.f, fxy.x, fxy.y)

        if tu.u > max_u:
            max_u = tu.u
            max_u_ring_cell = ring_cell
        if tu.u < min_u:
            min_u = tu.u
            min_u_ring_cell = ring_cell

        for dt, du in ((d, 0), (-d, 0), (0, d), (0, -d)):
            za = tu2za(tu.t + dt, tu.u + du)
            if abs(za.z) < 1 - 1e-12:
                az_values.append(_normalizeAngleRad(za.a))
                az_cell_refs.append(ring_cell)

    start_cell, end_cell, span = _arcBoundaryCells(az_values, az_cell_refs, 2 * PI)

    boundary_cells: set[int] = {min_u_ring_cell, max_u_ring_cell}
    if span < 2 * PI - 1e-10:
        boundary_cells.add(start_cell)
        boundary_cells.add(end_cell)

    min_z = math.inf
    max_z = -math.inf
    corner_azimuths: list[float] = []

    for ring_cell in boundary_cells:
        nest_cell = ring2nest(nside, ring_cell)
        corners = cornersNest(nside, nest_cell)
        center = pix2VecNest(nside, nest_cell)
        center_az = _normalizeAngleRad(math.atan2(center[1], center[0]))

        for cx, cy, cz in corners:
            z = max(-1.0, min(1.0, cz))
            min_z = min(min_z, z)
            max_z = max(max_z, z)
            if abs(z) >= 1 - 1e-12:
                az = center_az
            else:
                az = _normalizeAngleRad(math.atan2(cy, cx))
            corner_azimuths.append(az)

    min_lat = rad2Deg(math.asin(min_z))
    max_lat = rad2Deg(math.asin(max_z))

    if span >= 2 * PI - 1e-10:
        return (-180.0, min_lat, 180.0, max_lat)

    start_az, end_az, _ = _minimalCoveringArc(corner_azimuths, 2 * PI)

    return (
        _normalizeLon(rad2Deg(start_az)),
        min_lat,
        _normalizeLon(rad2Deg(end_az)),
        max_lat,
    )


def _normalizeLon(lon: float) -> float:
    """Normalize a longitude to (-180, 180]."""
    while lon > 180:
        lon -= 360
    while lon <= -180:
        lon += 360
    return lon


def _normalizeAngleRad(a: float) -> float:
    """Normalize an angle in radians to [0, 2π)."""
    period = 2 * PI
    while a < 0:
        a += period
    while a >= period:
        a -= period
    return a


def _normalizeDelta(d: float) -> float:
    """Normalize a longitude difference to (-180, 180]."""
    while d > 180:
        d -= 360
    while d <= -180:
        d += 360
    return d


def _arcBoundaryCells(
    azimuths: list[float],
    cells: list[int],
    period: float,
) -> tuple[int, int, float]:
    """Identify cells that bound the minimal longitude arc."""
    if not azimuths:
        cell0 = cells[0] if cells else -1
        return cell0, cell0, 0.0

    idx = list(range(len(azimuths)))
    idx.sort(key=lambda i: azimuths[i])

    max_gap = -math.inf
    gap_after = 0

    for i in range(len(idx)):
        cur = azimuths[idx[i]]
        if i == len(idx) - 1:
            nxt = azimuths[idx[0]] + period
        else:
            nxt = azimuths[idx[i + 1]]
        gap = nxt - cur
        if gap > max_gap:
            max_gap = gap
            gap_after = i

    start_idx = (gap_after + 1) % len(idx)
    return cells[idx[start_idx]], cells[idx[gap_after]], period - max_gap


def _minimalCoveringArc(
    values: list[float],
    period: float,
) -> tuple[float, float, float]:
    """Find the smallest arc on a circle that contains all given angles."""
    if len(values) == 1:
        value = values[0]
        return value, value, 0.0

    sorted_vals = sorted(values)
    max_gap = -math.inf
    max_gap_idx = 0

    for i in range(len(sorted_vals)):
        cur = sorted_vals[i]
        if i == len(sorted_vals) - 1:
            nxt = sorted_vals[0] + period
        else:
            nxt = sorted_vals[i + 1]
        gap = nxt - cur
        if gap > max_gap:
            max_gap = gap
            max_gap_idx = i

    start = sorted_vals[(max_gap_idx + 1) % len(sorted_vals)]
    end = sorted_vals[max_gap_idx]
    if end < start:
        end += period

    return start, end, end - start


def _queryBoxInclusiveImpl(nside: int, bbox: BBox, as_ring: bool) -> list[int]:
    """Internal implementation for both NESTED and RING box queries."""
    raw_min_lon, min_lat, raw_max_lon, max_lat = bbox
    results: list[int] = []

    min_lon = _normalizeLon(raw_min_lon)
    max_lon = _normalizeLon(raw_max_lon)

    crosses_antimeridian = min_lon > max_lon
    lon_span = (
        360 - min_lon + max_lon if crosses_antimeridian else max_lon - min_lon
    )

    if lon_span == 0 and raw_min_lon != raw_max_lon:
        lon_span = 360

    half_lon_span = lon_span / 2
    ref_lon = min_lon + half_lon_span

    pix_rad_deg = maxPixelRadius(nside) * (180 / PI)
    expanded_min_lat = max(-90.0, min_lat - pix_rad_deg)
    expanded_max_lat = min(90.0, max_lat + pix_rad_deg)
    expanded_half_lon_span = half_lon_span + pix_rad_deg

    d = PI_4 / nside
    z_top = math.sin(deg2Rad(expanded_max_lat))
    z_bot = math.sin(deg2Rad(expanded_min_lat))
    u_top = za2tu(z_top, 0).u
    u_bot = za2tu(z_bot, 0).u
    ring_min = max(1, int(math.floor((PI_2 - u_top) / d)))
    ring_max = min(4 * nside - 1, int(math.floor((PI_2 - u_bot) / d + 1)))

    for ring in range(ring_min, ring_max + 1):

        def _on_pixel(ipix_nest: int) -> None:
            if _pixelIntersectsBox(
                nside, ipix_nest, ref_lon, half_lon_span, min_lat, max_lat
            ):
                results.append(
                    nest2ring(nside, ipix_nest) if as_ring else ipix_nest
                )

        _walkRingFiltered(
            nside, ring, ref_lon, expanded_half_lon_span, _on_pixel
        )

    return results


def _walkRingFiltered(
    nside: int,
    ring: int,
    ref_lon: float,
    half_lon_span: float,
    cb: Callable[[int], None],
) -> None:
    """Walk pixels in a ring, filtering by longitude range."""
    u = PI_4 * (2 - ring / nside)
    t = PI_4 * (1 + (1 - (ring % 2)) / nside)

    begin = tu2fxy(nside, t, u)
    s = begin

    while True:
        ipix = fxy2nest(nside, s.f, s.x, s.y)
        lon, _ = pix2LonLatNest(nside, ipix)
        delta = _normalizeDelta(lon - ref_lon)
        if abs(delta) <= half_lon_span:
            cb(ipix)
        s = rightNextPixel(nside, s)
        if fxyEqual(s, begin):
            break


def _pixelIntersectsBox(
    nside: int,
    ipix: int,
    ref_lon: float,
    half_lon_span: float,
    min_lat: float,
    max_lat: float,
) -> bool:
    """Check if a pixel intersects a bounding box using corner geometry."""
    corners = cornersNestLonLat(nside, ipix)

    pix_min_lat = math.inf
    pix_max_lat = -math.inf
    pix_min_delta = math.inf
    pix_max_delta = -math.inf

    for lon, lat in corners:
        if lat < pix_min_lat:
            pix_min_lat = lat
        if lat > pix_max_lat:
            pix_max_lat = lat
        delta = _normalizeDelta(lon - ref_lon)
        if delta < pix_min_delta:
            pix_min_delta = delta
        if delta > pix_max_delta:
            pix_max_delta = delta

    if pix_max_lat < min_lat or pix_min_lat > max_lat:
        return False

    return pix_max_delta >= -half_lon_span and pix_min_delta <= half_lon_span
