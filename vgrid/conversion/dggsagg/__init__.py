"""Aggregate DGGS cells by parent resolution."""

from .h3agg import h3agg, h3agg_cli, h3_agg
from .s2agg import s2agg, s2agg_cli, s2_agg
from .a5agg import a5agg, a5agg_cli, a5_agg
from .rhealpixagg import rhealpixagg, rhealpixagg_cli, rhealpix_agg
from .isea4tagg import isea4tagg, isea4tagg_cli, isea4t_agg
from .isea3hagg import isea3hagg, isea3hagg_cli, isea3h_agg
from .easeagg import easeagg, easeagg_cli, ease_agg
from .dggalagg import dggalagg, dggalagg_cli, dggal_agg
from .qtmagg import qtmagg, qtmagg_cli, qtm_agg
from .olcagg import olcagg, olcagg_cli, olc_agg
from .geohashagg import geohashagg, geohashagg_cli, geohash_agg
from .tilecodeagg import tilecodeagg, tilecodeagg_cli, tilecode_agg
from .quadkeyagg import quadkeyagg, quadkeyagg_cli, quadkey_agg
from .digipinagg import digipinagg, digipinagg_cli, digipin_agg

__all__ = [
    "h3agg",
    "h3_agg",
    "h3agg_cli",
    "s2agg",
    "s2_agg",
    "s2agg_cli",
    "a5agg",
    "a5_agg",
    "a5agg_cli",
    "rhealpixagg",
    "rhealpix_agg",
    "rhealpixagg_cli",
    "isea4tagg",
    "isea4t_agg",
    "isea4tagg_cli",
    "isea3hagg",
    "isea3h_agg",
    "isea3hagg_cli",
    "easeagg",
    "ease_agg",
    "easeagg_cli",
    "dggalagg",
    "dggal_agg",
    "dggalagg_cli",
    "qtmagg",
    "qtm_agg",
    "qtmagg_cli",
    "olcagg",
    "olc_agg",
    "olcagg_cli",
    "geohashagg",
    "geohash_agg",
    "geohashagg_cli",
    "tilecodeagg",
    "tilecode_agg",
    "tilecodeagg_cli",
    "quadkeyagg",
    "quadkey_agg",
    "quadkeyagg_cli",
    "digipinagg",
    "digipin_agg",
    "digipinagg_cli",
]
