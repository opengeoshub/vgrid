"""Hierarchical pixel operations: parent/child relationships."""

from __future__ import annotations

import math

from ..schemes.conversion import nest2ring, ring2nest


def nestParent(ipix: int) -> int:
    """Get the parent pixel index in the NESTED scheme."""
    return ipix >> 2


def nestChildren(ipix: int) -> tuple[int, int, int, int]:
    """Get the 4 child pixel indices in the NESTED scheme."""
    base = ipix << 2
    return (base, base + 1, base + 2, base + 3)


def nestAncestor(ipix: int, levels: int) -> int:
    """Get the parent pixel at a specific ancestor level."""
    return ipix >> (2 * levels)


def nestDescendants(ipix: int, levels: int) -> list[int]:
    """Get all descendants at a specific depth."""
    if levels <= 0:
        return [ipix]

    count = 1 << (2 * levels)
    base = ipix << (2 * levels)
    return [base + i for i in range(count)]


def isNestAncestor(
    ancestor: int,
    ancestorNside: int,
    descendant: int,
    descendantNside: int,
) -> bool:
    """Check if one nested pixel is an ancestor of another."""
    if ancestorNside >= descendantNside:
        return False

    ratio = descendantNside / ancestorNside
    levels = math.log2(ratio)

    if not levels.is_integer():
        return False

    return (descendant >> (2 * int(levels))) == ancestor


def ringParent(nside: int, ipix: int) -> int:
    """Get the parent pixel index in the RING scheme."""
    parent_nside = nside >> 1
    nest_idx = ring2nest(nside, ipix)
    return nest2ring(parent_nside, nestParent(nest_idx))


def ringChildren(nside: int, ipix: int) -> tuple[int, int, int, int]:
    """Get the 4 child pixel indices in the RING scheme."""
    child_nside = nside << 1
    nest_idx = ring2nest(nside, ipix)
    nest_kids = nestChildren(nest_idx)
    return tuple(nest2ring(child_nside, k) for k in nest_kids)  # type: ignore[return-value]


def ringAncestor(nside: int, ipix: int, levels: int) -> int:
    """Get the ancestor pixel at a specific level in the RING scheme."""
    ancestor_nside = nside >> levels
    nest_idx = ring2nest(nside, ipix)
    return nest2ring(ancestor_nside, nestAncestor(nest_idx, levels))


def ringDescendants(nside: int, ipix: int, levels: int) -> list[int]:
    """Get all descendants at a specific depth in the RING scheme."""
    if levels <= 0:
        return [ipix]

    desc_nside = nside << levels
    nest_idx = ring2nest(nside, ipix)
    return [nest2ring(desc_nside, d) for d in nestDescendants(nest_idx, levels)]


def isRingAncestor(
    ancestor: int,
    ancestorNside: int,
    descendant: int,
    descendantNside: int,
) -> bool:
    """Check if one ring pixel is an ancestor of another."""
    if ancestorNside >= descendantNside:
        return False

    ratio = descendantNside / ancestorNside
    levels = math.log2(ratio)
    if not levels.is_integer():
        return False

    nest_anc = ring2nest(ancestorNside, ancestor)
    nest_desc = ring2nest(descendantNside, descendant)
    return isNestAncestor(nest_anc, ancestorNside, nest_desc, descendantNside)
