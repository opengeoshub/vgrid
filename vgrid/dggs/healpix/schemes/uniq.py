"""UNIQ (Unique Identifier) scheme implementation."""

from __future__ import annotations

from ..types import OrderPix


def orderpix2uniq(order: int, ipix: int) -> int:
    """Pack (order, ipix_nested) into a single unique integer."""
    return 4 * ((1 << (2 * order)) - 1) + ipix


def uniq2orderpix(uniq: int) -> OrderPix:
    """Unpack a unique identifier into (order, ipix_nested)."""
    if uniq < 0:
        raise ValueError(f"uniq must be non-negative, got {uniq}")

    order = 0
    l = (uniq >> 2) + 1
    while l >= 4:
        l >>= 2
        order += 1

    ipix = uniq - (((1 << (2 * order)) - 1) << 2)
    return OrderPix(order=order, ipix=ipix)


def uniqParent(uniq: int) -> int:
    """Get the parent uniq value (coarser resolution)."""
    order, ipix = uniq2orderpix(uniq)
    if order == 0:
        raise ValueError("Base pixels (order 0) have no parent")
    return orderpix2uniq(order - 1, ipix >> 2)


def uniqChildren(uniq: int) -> tuple[int, int, int, int]:
    """Get the 4 child uniq values (finer resolution)."""
    order, ipix = uniq2orderpix(uniq)
    child_order = order + 1
    base = ipix << 2
    return (
        orderpix2uniq(child_order, base),
        orderpix2uniq(child_order, base + 1),
        orderpix2uniq(child_order, base + 2),
        orderpix2uniq(child_order, base + 3),
    )
