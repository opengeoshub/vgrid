"""General utility functions used throughout the library."""

from __future__ import annotations

import math

from .constants import DEG2RAD, RAD2DEG


def rad2deg(radians: float) -> float:
    """Convert radians to degrees."""
    return radians * RAD2DEG


def deg2rad(degrees: float) -> float:
    """Convert degrees to radians."""
    return degrees * DEG2RAD


def square(x: float) -> float:
    """Return x squared."""
    return x * x


def clip(value: float, min_val: float, max_val: float) -> float:
    """Clamp value to [min_val, max_val]."""
    if value < min_val:
        return min_val
    if value > max_val:
        return max_val
    return value


def wrap(value: float, period: float) -> float:
    """Wrap value to [0, period).

    Handles negative values correctly (unlike a naive modulo in some languages).
    """
    return value % period


def ilog2(x: int) -> int:
    """Integer log base 2 (floor)."""
    return int(math.floor(math.log2(x)))


def assert_(condition: bool, message: str | None = None) -> None:
    """Debug assertion that raises AssertionError on failure."""
    if not condition:
        raise AssertionError(message or "assertion failed")


# Keep camelCase aliases matching the TypeScript API.
rad2Deg = rad2deg
deg2Rad = deg2rad
