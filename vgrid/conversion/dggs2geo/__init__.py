"""
DGGS to Geographic coordinate conversion functions.

This module provides functions to convert various discrete global grid systems (DGGS)
to shapely polygon and geojson.
"""

from .h32geo import (
    h32geo,
    h32geojson,
    h32geo_cli,
    h32geojson_cli,
)
from .s22geo import (
    s22geo,
    s22geojson,
    s22geo_cli,
    s22geojson_cli,
)
from .rhealpix2geo import (
    rhealpix2geo,
    rhealpix2geojson,
    rhealpix2geo_cli,
    rhealpix2geojson_cli,
)
from .isea4t2geo import (
    isea4t2geo,
    isea4t2geojson,
    isea4t2geo_cli,
    isea4t2geojson_cli,
)
from .isea3h2geo import (
    isea3h2geo,
    isea3h2geojson,
    isea3h2geo_cli,
    isea3h2geojson_cli,
)
from .ease2geo import (
    ease2geo,
    ease2geojson,
    ease2geo_cli,
    ease2geojson_cli,
)
from .qtm2geo import (
    qtm2geo,
    qtm2geojson,
    qtm2geo_cli,
    qtm2geojson_cli,
)
from .olc2geo import (
    olc2geo,
    olc2geojson,
    olc2geo_cli,
    olc2geojson_cli,
)
from .geohash2geo import (
    geohash2geo,
    geohash2geojson,
    geohash2geo_cli,
    geohash2geojson_cli,
)
from .georef2geo import (
    georef2geo,
    georef2geojson,
    georef2geo_cli,
    georef2geojson_cli,
)
from .mgrs2geo import (
    mgrs2geo,
    mgrs2geojson,
    mgrs2geo_cli,
    mgrs2geojson_cli,
)
from .tilecode2geo import (
    tilecode2geo,
    tilecode2geojson,
    tilecode2geo_cli,
    tilecode2geojson_cli,
)
from .quadkey2geo import (
    quadkey2geo,
    quadkey2geojson,
    quadkey2geo_cli,
    quadkey2geojson_cli,
)
from .maidenhead2geo import (
    maidenhead2geo,
    maidenhead2geojson,
    maidenhead2geo_cli,
    maidenhead2geojson_cli,
)
from .gars2geo import (
    gars2geo,
    gars2geojson,
    gars2geo_cli,
    gars2geojson_cli,
)
from .digipin2geo import (
    digipin2geo,
    digipin2geojson,
    digipin2geo_cli,
    digipin2geojson_cli,
)
from .a52geo import (
    a52geo,
    a52geojson,
    a52geo_cli,
    a52geojson_cli,
)
from .healpix2geo import (
    healpix2geo,
    healpix2geojson,
    healpix2geo_cli,
    healpix2geojson_cli,
)

from .dggal2geo import (
    dggal2geo,
    dggal2geojson,
    dggal2geo_cli,
    dggal2geojson_cli,
)
from .dggrid2geo import (
    dggrid2geo,
    dggrid2geojson,
    dggrid2geo_cli,
    dggrid2geojson_cli,
)

__all__ = [
    "h32geo",
    "h32geojson",
    "h32geo_cli",
    "h32geojson_cli",
    "s22geo",
    "s22geojson",
    "s22geo_cli",
    "s22geojson_cli",
    "rhealpix2geo",
    "rhealpix2geojson",
    "rhealpix2geo_cli",
    "rhealpix2geojson_cli",
    "isea4t2geo",
    "isea4t2geojson",
    "isea4t2geo_cli",
    "isea4t2geojson_cli",
    "isea3h2geo",
    "isea3h2geojson",
    "isea3h2geo_cli",
    "isea3h2geojson_cli",
    "ease2geo",
    "ease2geojson",
    "ease2geo_cli",
    "ease2geojson_cli",
    "qtm2geo",
    "qtm2geojson",
    "qtm2geo_cli",
    "qtm2geojson_cli",
    "olc2geo",
    "olc2geojson",
    "olc2geo_cli",
    "olc2geojson_cli",
    "geohash2geo",
    "geohash2geojson",
    "geohash2geo_cli",
    "geohash2geojson_cli",
    "georef2geo",
    "georef2geojson",
    "georef2geo_cli",
    "georef2geojson_cli",
    "mgrs2geo",
    "mgrs2geojson",
    "mgrs2geo_cli",
    "mgrs2geojson_cli",
    "tilecode2geo",
    "tilecode2geojson",
    "tilecode2geo_cli",
    "tilecode2geojson_cli",
    "quadkey2geo",
    "quadkey2geojson",
    "quadkey2geo_cli",
    "quadkey2geojson_cli",
    "maidenhead2geo",
    "maidenhead2geojson",
    "maidenhead2geo_cli",
    "maidenhead2geojson_cli",
    "gars2geo",
    "gars2geojson",
    "gars2geo_cli",
    "gars2geojson_cli",
    "digipin2geo",
    "digipin2geojson",
    "digipin2geo_cli",
    "digipin2geojson_cli",
    "a52geo",
    "a52geojson",
    "a52geo_cli",
    "a52geojson_cli",
    "healpix2geo",
    "healpix2geojson",
    "healpix2geo_cli",
    "healpix2geojson_cli",
    "dggal2geo",
    "dggal2geojson",
    "dggal2geo_cli",
    "dggal2geojson_cli",
    "dggrid2geo",
    "dggrid2geojson",
    "dggrid2geo_cli",
    "dggrid2geojson_cli",
]
