"""HEALPix (Hierarchical Equal Area isoLatitude Pixelization) Python library.

Ported from healpix-ts (https://github.com/developmentseed/healpix-ts),
itself based on michitaro/healpix.

Based on: Górski et al. (2005) http://iopscience.iop.org/article/10.1086/427976/pdf

Quick start::

    from vgrid.dggs.healpix import ang2PixNest, pix2LonLatNest, order2nside

    nside = order2nside(8)  # nside = 256
    ipix = ang2PixNest(nside, math.pi / 4, math.pi / 2)
    lon, lat = pix2LonLatNest(nside, ipix)
"""

from __future__ import annotations

from .constants import DEG2RAD, PI, PI2, PI_2, PI_4, PI_8, RAD2DEG
from .coordinates import (
    ang2vec,
    angularDistance,
    sigma,
    tu2vec,
    tu2za,
    vec2ang,
    vec2za,
    za2tu,
    za2vec,
)
from .geo import (
    cornersNestLonLat,
    cornersRingLonLat,
    lonLat2PixNest,
    lonLat2PixRing,
    pix2LonLatNest,
    pix2LonLatRing,
    vec2LonLat,
)
from .lookup import (
    ang2PixNest,
    ang2PixRing,
    cornersNest,
    cornersRing,
    pix2AngNest,
    pix2AngRing,
    pix2VecNest,
    pix2VecRing,
    pixcoord2VecNest,
    pixcoord2VecRing,
    vec2PixNest,
    vec2PixRing,
)
from .pixel import (
    fxy2tu,
    fxyCorners,
    fxyEqual,
    fxySubpixel,
    isNestAncestor,
    isRingAncestor,
    maxPixelRadius,
    nestAncestor,
    nestChildren,
    nestDescendants,
    nestParent,
    rightNextPixel,
    ringAncestor,
    ringChildren,
    ringDescendants,
    ringParent,
    tu2fpq,
    tu2fxy,
)
from .query import (
    cellsToBoundingBox,
    queryBoxInclusiveNest,
    queryBoxInclusiveRing,
    queryDiscInclusiveNest,
    queryDiscInclusiveRing,
)
from .resolution import (
    nside2npix,
    nside2order,
    nside2pixarea,
    nside2resol,
    order2nside,
)
from .schemes import (
    bitCombine,
    bitDecombine,
    fxy2nest,
    fxy2ring,
    nest2fxy,
    nest2ring,
    orderpix2uniq,
    ring2fxy,
    ring2nest,
    uniq2orderpix,
    uniqChildren,
    uniqParent,
)
from .types import (
    AngularCoords,
    BBox,
    FXY,
    LonLat,
    OrderPix,
    TU,
    V3,
    ZA,
)
from .utils import clip, deg2Rad, ilog2, rad2Deg, square, wrap

__all__ = [
    # Types
    "AngularCoords",
    "BBox",
    "FXY",
    "LonLat",
    "OrderPix",
    "TU",
    "V3",
    "ZA",
    # Constants
    "DEG2RAD",
    "PI",
    "PI2",
    "PI_2",
    "PI_4",
    "PI_8",
    "RAD2DEG",
    # Utils
    "clip",
    "deg2Rad",
    "ilog2",
    "rad2Deg",
    "square",
    "wrap",
    # Resolution
    "nside2npix",
    "nside2order",
    "nside2pixarea",
    "nside2resol",
    "order2nside",
    # Coordinates
    "ang2vec",
    "angularDistance",
    "sigma",
    "tu2vec",
    "tu2za",
    "vec2ang",
    "vec2za",
    "za2tu",
    "za2vec",
    # Pixel
    "fxy2tu",
    "fxyCorners",
    "fxyEqual",
    "fxySubpixel",
    "isNestAncestor",
    "isRingAncestor",
    "maxPixelRadius",
    "nestAncestor",
    "nestChildren",
    "nestDescendants",
    "nestParent",
    "rightNextPixel",
    "ringAncestor",
    "ringChildren",
    "ringDescendants",
    "ringParent",
    "tu2fpq",
    "tu2fxy",
    # Schemes
    "bitCombine",
    "bitDecombine",
    "fxy2nest",
    "fxy2ring",
    "nest2fxy",
    "nest2ring",
    "orderpix2uniq",
    "ring2fxy",
    "ring2nest",
    "uniq2orderpix",
    "uniqChildren",
    "uniqParent",
    # Lookup
    "ang2PixNest",
    "ang2PixRing",
    "cornersNest",
    "cornersRing",
    "pix2AngNest",
    "pix2AngRing",
    "pix2VecNest",
    "pix2VecRing",
    "pixcoord2VecNest",
    "pixcoord2VecRing",
    "vec2PixNest",
    "vec2PixRing",
    # Query
    "cellsToBoundingBox",
    "queryBoxInclusiveNest",
    "queryBoxInclusiveRing",
    "queryDiscInclusiveNest",
    "queryDiscInclusiveRing",
    # Geo
    "cornersNestLonLat",
    "cornersRingLonLat",
    "lonLat2PixNest",
    "lonLat2PixRing",
    "pix2LonLatNest",
    "pix2LonLatRing",
    "vec2LonLat",
]
