"""
World Geographic Reference System (GEOREF), adapted from GeographicLib's Georef.

Copyright (c) Charles Karney (2015-2024) <karney@alum.mit.edu> and licensed under the
MIT/X11 License.  For more information, see https://geographiclib.sourceforge.io/
"""

from __future__ import annotations

import math
import sys
from typing import Tuple

digits_ = "0123456789"
lontile_ = "ABCDEFGHJKLMNPQRSTUVWXYZ"
lattile_ = "ABCDEFGHJKLM"
degrees_ = "ABCDEFGHJKLMNPQ"
tile_ = 15
lonorig_ = -180
latorig_ = -90
base_ = 10
baselen_ = 4
maxprec_ = 11
maxlen_ = 2 * maxprec_ + 4
_qd = 90.0  # Math::qd


def _encode_out_len(prec: int) -> int:
    """GEOREF string length: prec 0 → 2, prec 1 → 4, prec ≥ 2 → 2*prec + 4."""
    if prec <= 0:
        return 2
    if prec == 1:
        return 4
    return 2 * prec + 4


def _decode_prec_from_len(leng: int) -> int:
    """Precision implied by string length (inverse of :func:`_encode_out_len`)."""
    if leng == 2:
        return 0
    if leng == 4:
        return 1
    if leng < 8 or (leng - 4) % 2 != 0:
        raise GeorefException(f"Invalid georef length: {leng}")
    prec = (leng - 4) // 2
    if prec < 2 or prec > maxprec_:
        raise GeorefException(f"Invalid georef length: {leng}")
    return prec


class GeorefException(Exception):
    pass


def _ang_normalize(lon: float) -> float:
    """Longitude in degrees to [-180, 180), so tile indexing stays in range (e.g. lon=-180)."""
    x = (lon + 180.0) % 360.0
    if x < 0:
        x += 360.0
    return x - 180.0


def _lookup(alphabet: str, ch: str) -> int:
    k = alphabet.find(ch)
    return k


def _find_first_not_of(s: str, charset: str, start: int = 0) -> int:
    for i in range(start, len(s)):
        if s[i] not in charset:
            return i
    return -1


def encode(lat: float, lon: float, prec: int) -> str:
    """
    Convert geographic coordinates to a GEOREF string (GeographicLib Forward).

    If lat or lon is NaN, returns ``\"INVALID\"``.
    """
    if abs(lat) > _qd:
        raise GeorefException(
            f"Latitude {lat} not in [-{_qd}, {_qd}]"
        )
    if math.isnan(lat) or math.isnan(lon):
        return "INVALID"
    lon = _ang_normalize(lon)
    if lat == _qd:
        lat = lat * (1 - sys.float_info.epsilon / 2)
    prec = max(-1, min(maxprec_, int(prec)))
    m = 60000000000
    x = int(math.floor(lon * m)) - lonorig_ * m
    y = int(math.floor(lat * m)) - latorig_ * m
    ilon = x // m
    ilat = y // m
    out_len = _encode_out_len(prec)
    georef1 = [""] * out_len
    georef1[0] = lontile_[ilon // tile_]
    georef1[1] = lattile_[ilat // tile_]
    if prec >= 1:
        georef1[2] = degrees_[ilon % tile_]
        georef1[3] = degrees_[ilat % tile_]
        if prec >= 2:
            x -= m * ilon
            y -= m * ilat
            d = int(pow(base_, maxprec_ - prec))
            x //= d
            y //= d
            c = prec
            while c:
                c -= 1
                georef1[baselen_ + c] = digits_[x % base_]
                x //= base_
                georef1[baselen_ + c + prec] = digits_[y % base_]
                y //= base_
    return "".join(georef1)


def decode(georef: str, centerp: bool = False) -> Tuple[float, float, int]:
    """
    Convert a GEOREF string to latitude, longitude, and precision (GeographicLib Reverse).

    Letter case is ignored. If the string begins with ``\"INV\"`` (case-insensitive),
    ``lat`` and ``lon`` are NaN (GeographicLib sets NaN and leaves prec unchanged; here
    ``prec`` is set to ``-1``).
    """
    if georef is None:
        raise GeorefException("Invalid Georef string: None")
    s = georef.strip()
    leng = len(s)
    if leng >= 3 and s[0:3].upper() == "INV":
        return float("nan"), float("nan"), -1
    if leng < baselen_ - 2:
        raise GeorefException(
            f"Georef must start with at least 2 letters: {georef!r}"
        )
    u = s.upper()
    prec1 = _decode_prec_from_len(leng)
    k = _lookup(lontile_, u[0])
    if k < 0:
        raise GeorefException(f"Bad longitude tile letter in georef: {georef!r}")
    lon1 = k + lonorig_ / tile_
    k = _lookup(lattile_, u[1])
    if k < 0:
        raise GeorefException(f"Bad latitude tile letter in georef: {georef!r}")
    lat1 = k + latorig_ / tile_
    unit = 1.0
    if leng > 2:
        unit *= tile_
        k = _lookup(degrees_, u[2])
        if k < 0:
            raise GeorefException(
                f"Bad longitude degree letter in georef: {georef!r}"
            )
        lon1 = lon1 * tile_ + k
        if leng < 4:
            raise GeorefException(
                f"Missing latitude degree letter in georef: {georef!r}"
            )
        k = _lookup(degrees_, u[3])
        if k < 0:
            raise GeorefException(
                f"Bad latitude degree letter in georef: {georef!r}"
            )
        lat1 = lat1 * tile_ + k
        if prec1 >= 2:
            if _find_first_not_of(u, digits_, baselen_) != -1:
                raise GeorefException(
                    f"Non digits in trailing portion of georef: {u[baselen_:]!r}"
                )
            if (leng - baselen_) % 2:
                raise GeorefException(
                    f"Georef must end with an even number of digits: {u[baselen_:]!r}"
                )
            if prec1 > maxprec_:
                raise GeorefException(
                    f"More than {2 * maxprec_} digits in georef: {u[baselen_:]!r}"
                )
            for i in range(prec1):
                m = base_ if i else 6
                unit *= m
                x = _lookup(digits_, u[baselen_ + i])
                y = _lookup(digits_, u[baselen_ + i + prec1])
                if not (i or (x < m and y < m)):
                    raise GeorefException(
                        "Minutes terms in georef must be less than 60 "
                        f"{u[baselen_:]!r}"
                    )
                lon1 = m * lon1 + x
                lat1 = m * lat1 + y
    if centerp:
        unit *= 2
        lat1 = 2 * lat1 + 1
        lon1 = 2 * lon1 + 1
    lat = (tile_ * lat1) / unit
    lon = (tile_ * lon1) / unit
    return lat, lon, prec1

## Added by Vgrid
def georefcell(georef_id: str):
    """Decode GEOREF to center and bounding box (Vgrid helper)."""
    center_lat, center_lon, resolution = decode(georef_id, True)
    if resolution == 0:
        grid_size = 15.0
    elif resolution == 1:
        grid_size = 1.0
    elif resolution > 1:
        grid_size = 1 / (6 * 10 ** (resolution - 1))
    else:
        raise ValueError(f"Invalid resolution: {resolution}")

    min_lon = float(int(center_lon // grid_size) * grid_size)
    max_lon = min_lon + grid_size
    min_lat = float(int(center_lat // grid_size) * grid_size)
    max_lat = min_lat + grid_size

    return center_lat, center_lon, min_lat, min_lon, max_lat, max_lon, resolution
