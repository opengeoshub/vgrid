"""
Raster to DGGS conversion functions.

This submodule provides functions to convert raster data to various
discrete global grid systems (DGGS).
"""

from .raster2h3 import raster2h3, raster2h3_cli
from .raster2s2 import raster2s2, raster2s2_cli
from .raster2a5 import raster2a5, raster2a5_cli
from .raster2healpix import raster2healpix, raster2healpix_cli
from .raster2rhealpix import raster2rhealpix, raster2rhealpix_cli
from .raster2isea4t import raster2isea4t, raster2isea4t_cli
from .raster2dggrid import raster2dggrid, raster2dggrid_cli
from .raster2dggal import raster2dggal, raster2dggal_cli
from .raster2qtm import raster2qtm, raster2qtm_cli
from .raster2olc import raster2olc, raster2olc_cli
from .raster2geohash import raster2geohash, raster2geohash_cli
from .raster2georef import raster2georef, raster2georef_cli
from .raster2gars import raster2gars, raster2gars_cli
from .raster2maidenhead import raster2maidenhead, raster2maidenhead_cli
from .raster2tilecode import raster2tilecode, raster2tilecode_cli
from .raster2quadkey import raster2quadkey, raster2quadkey_cli
from .raster2digipin import raster2digipin, raster2digipin_cli

__all__ = [
    "raster2h3",
    "raster2h3_cli",
    "raster2s2",
    "raster2s2_cli",
    "raster2a5",
    "raster2a5_cli",
    "raster2healpix",
    "raster2healpix_cli",
    "raster2rhealpix",
    "raster2rhealpix_cli",
    "raster2isea4t",
    "raster2isea4t_cli",
    "raster2dggrid",
    "raster2dggrid_cli",
    "raster2dggal",
    "raster2dggal_cli",
    "raster2qtm",
    "raster2qtm_cli",
    "raster2olc",
    "raster2olc_cli",
    "raster2geohash",
    "raster2geohash_cli",
    "raster2georef",
    "raster2georef_cli",
    "raster2gars",
    "raster2gars_cli",
    "raster2maidenhead",
    "raster2maidenhead_cli",
    "raster2tilecode",
    "raster2tilecode_cli",
    "raster2quadkey",
    "raster2quadkey_cli",
    "raster2digipin",
    "raster2digipin_cli",
]
