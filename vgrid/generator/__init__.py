"""
Generator module for vgrid.

This module provides functions to generate discrete global grid systems (DGGS)
for various coordinate systems and geographic areas.
"""

# Main grid generation functions
from .h3grid import h3grid, h3grid_cli
from .s2grid import s2grid, s2grid_cli
from .a5grid import a5grid, a5grid_cli
from .rhealpixgrid import rhealpixgrid, rhealpixgrid_cli
from .isea4tgrid import isea4tgrid, isea4tgrid_cli
from .isea3hgrid import isea3hgrid, isea3hgrid_cli
from .easegrid import easegrid, easegrid_cli
from .dggridgen import dggridgen, dggridgen_cli
from .dggalgen import dggalgen, dggalgen_cli
from .qtmgrid import qtmgrid, qtmgrid_cli
from .olcgrid import olcgrid, olcgrid_cli
from .geohashgrid import geohashgrid, geohashgrid_cli
from .georefgrid import georefgrid, georefgrid_cli
from .mgrsgrid import mgrsgrid, mgrsgrid_cli
from .tilecodegrid import tilecodegrid, tilecodegrid_cli
from .quadkeygrid import quadkeygrid, quadkeygrid_cli
from .maidenheadgrid import maidenheadgrid, maidenheadgrid_cli
from .garsgrid import garsgrid, garsgrid_cli
from .digipingrid import digipingrid, digipingrid_cli

__all__ = [
    "h3grid",
    "h3grid_cli",
    "s2grid",
    "s2grid_cli",
    "a5grid",
    "a5grid_cli",
    "rhealpixgrid",
    "rhealpixgrid_cli",
    "isea4tgrid",
    "isea4tgrid_cli",
    "isea3hgrid",
    "isea3hgrid_cli",
    "easegrid",
    "easegrid_cli",
    "dggridgen",
    "dggridgen_cli",
    "dggalgen",
    "dggalgen_cli",
    "qtmgrid",
    "qtmgrid_cli",
    "olcgrid",
    "olcgrid_cli",
    "geohashgrid",
    "geohashgrid_cli",
    "georefgrid",
    "georefgrid_cli",
    "mgrsgrid",
    "mgrsgrid_cli",
    "tilecodegrid",
    "tilecodegrid_cli",
    "quadkeygrid",
    "quadkeygrid_cli",
    "maidenheadgrid",
    "maidenheadgrid_cli",
    "garsgrid",
    "garsgrid_cli",
    "digipingrid",
    "digipingrid_cli",
]
