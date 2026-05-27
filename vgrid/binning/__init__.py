"""
Binning module for vgrid.

This module provides functions to bin and aggregate data using various
discrete global grid systems (DGGS), including statistical analysis
and data categorization.
"""

# Import all binning functions
from .h3bin import h3bin, h3bin_cli
from .s2bin import s2bin, s2bin_cli
from .a5bin import a5bin, a5bin_cli
from .rhealpixbin import rhealpixbin, rhealpixbin_cli
from .isea4tbin import isea4tbin, isea4tbin_cli
from .easebin import easebin, easebin_cli
from .dggalbin import dggalbin, dggalbin_cli
from .dggridbin import dggridbin, dggridbin_cli
from .qtmbin import qtmbin, qtmbin_cli
from .olcbin import olcbin, olcbin_cli
from .geohashbin import geohashbin, geohashbin_cli
from .georefbin import georefbin, georefbin_cli
from .garsbin import garsbin, garsbin_cli
from .maidenheadbin import maidenheadbin, maidenheadbin_cli
from .tilecodebin import tilecodebin, tilecodebin_cli
from .quadkeybin import quadkeybin, quadkeybin_cli
from .polygonbin import polygonbin, polygonbin_cli


__all__ = [
    "h3bin",
    "h3bin_cli",
    "s2bin",
    "s2bin_cli",
    "a5bin",
    "a5bin_cli",
    "rhealpixbin",
    "rhealpixbin_cli",
    "isea4tbin",
    "isea4tbin_cli",
    "easebin",
    "easebin_cli",
    "dggalbin",
    "dggalbin_cli",
    "dggridbin",
    "dggridbin_cli",
    "qtmbin",
    "qtmbin_cli",
    "olcbin",
    "olcbin_cli",
    "geohashbin",
    "geohashbin_cli",
    "georefbin",
    "georefbin_cli",
    "garsbin",
    "garsbin_cli",
    "maidenheadbin",
    "maidenheadbin_cli",
    "tilecodebin",
    "tilecodebin_cli",
    "quadkeybin",
    "quadkeybin_cli",
    "polygonbin",
    "polygonbin_cli",
]
