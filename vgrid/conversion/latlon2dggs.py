"""
Latitude/Longitude to DGGS Conversion Module

This module provides functions to convert latitude and longitude coordinates to various
Discrete Global Grid System (DGGS) cell identifiers. It supports multiple DGGS types
including H3, S2, A5, HEALPix, RHEALPix, ISEA4T, ISEA3H, DGGRID, DGGAL, EASE, QTM, OLC,
Geohash, GEOREF, MGRS, Tilecode, Quadkey, Maidenhead, GARS, and DIGIPIN.

Each DGGS type has its own resolution range and addressing scheme. The module includes
both programmatic functions and command-line interfaces (CLI) for each conversion type.

Each function follows the pattern: latlon2<dggs_type>(lat, lon, res=<default_res>)
and returns the corresponding cell identifier as a string or appropriate data type.

CLI Usage Examples:
    latlon2h3 10.775275567242561 106.70679737574993 13
    latlon2s2 10.775275567242561 106.70679737574993 21
    latlon2dggrid ISEA7H 10.775275567242561 106.70679737574993 13
    latlon2dggal gnosis 10.775275567242561 106.70679737574993 8
    latlon2digipin 28.6139 77.2090 10
"""

from vgrid.dggs import s2, olc, geohash, georef, mgrs, maidenhead, tilecode, qtm
from vgrid.dggs.digipin import latlon_to_digipin
from vgrid.dggs.healpix import lonLat2PixNest, order2nside, orderpix2uniq
import h3
import a5
from dggal import *

from gars_field.garsgrid import GARSGrid
from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.dggs.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

import platform
import argparse

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.dggs.eaggr.shapes.lat_long_point import LatLongPoint
    from vgrid.dggs.eaggr.enums.model import Model

    isea4t_dggs = Eaggr(Model.ISEA4T)
    isea3h_dggs = Eaggr(Model.ISEA3H)

from vgrid.utils.constants import (
    ISEA3H_RES_ACCURACY_DICT,
    ISEA4T_RES_ACCURACY_DICT,
    DGGAL_TYPES,
)

import geopandas as gpd
from dggrid4py.dggrid_runner import output_address_types

from shapely import Point
from vgrid.utils.io import (
    validate_h3_resolution,
    validate_s2_resolution,
    validate_a5_resolution,
    validate_healpix_resolution,
    validate_rhealpix_resolution,
    validate_isea4t_resolution,
    validate_isea3h_resolution,
    validate_ease_resolution,
    validate_qtm_resolution,
    validate_geohash_resolution,
    validate_georef_resolution,
    validate_mgrs_resolution,
    validate_tilecode_resolution,
    validate_quadkey_resolution,
    validate_maidenhead_resolution,
    validate_gars_resolution,
    validate_olc_resolution,
    olc_resolutions,
    validate_digipin_resolution,
    validate_dggal_type,
    validate_dggal_resolution,
    validate_dggrid_type,
    validate_dggrid_resolution,
    create_dggrid_instance,
)
from vgrid.utils.constants import DGGS_TYPES, DGGRID_TYPES
from ease_dggs.dggs.grid_addressing import geos_to_grid_ids

app = Application(appGlobals=globals())
pydggal_setup(app)


def latlon2h3(lat, lon, res):
    """
    Convert latitude and longitude to H3 cell identifier.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): H3 resolution level [0-15]

    Returns:
        str: H3 cell identifier

    Example:
        >>> latlon2h3(10.775275567242561, 106.70679737574993, 13)
        '8b194e64992ffff'
    """
    if res is None:
        res = DGGS_TYPES["h3"]["default_res"]
    res = validate_h3_resolution(res)
    h3_id = h3.latlng_to_cell(lat, lon, res)
    return h3_id


def latlon2h3_cli():
    """
    Command-line interface for latlon2h3.
    """
    min_res = DGGS_TYPES["h3"]["min_res"]
    max_res = DGGS_TYPES["h3"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to H3 code at a specific Resolution [{min_res}..{max_res}]. \
                                    Usage: latlon2h3 <lat> <lon> <res> [{min_res}..{max_res}]. \
                                    Ex: latlon2h3 10.775275567242561 106.70679737574993 13"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["h3"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    h3_id = latlon2h3(lat, lon, res)
    print(h3_id)


def latlon2s2(lat, lon, res=None):
    """
    Convert latitude and longitude to S2 cell identifier.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): S2 resolution level [0-30]

    Returns:
        str: S2 cell token

    Example:
        >>> latlon2s2(10.775275567242561, 106.70679737574993, 21)
        '1d4b9b0b8c8c8c8c'
    """
    if res is None:
        res = DGGS_TYPES["s2"]["default_res"]
    res = validate_s2_resolution(res)
    lat_lng = s2.LatLng.from_degrees(lat, lon)
    s2_cell = s2.CellId.from_lat_lng(lat_lng)  # return S2 cell at max level 30
    s2_cell_res = s2_cell.parent(res)  # get S2 cell at resolution
    s2_token = s2.CellId.to_token(
        s2_cell_res
    )  # get Cell ID Token, shorter than cell_id.id()
    return s2_token


def latlon2s2_cli():
    """
    Command-line interface for latlon2s2.
    """
    min_res = DGGS_TYPES["s2"]["min_res"]
    max_res = DGGS_TYPES["s2"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to S2 code at a specific Resolution [0..30]. \
                                     Usage: latlon2s2 <lat> <lon> <res> [0..30]. \
                                     Ex: latlon2s2 10.775275567242561 106.70679737574993 21"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["s2"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()
    res = args.res
    lat = args.lat
    lon = args.lon
    s2_token = latlon2s2(lat, lon, res)
    print(s2_token)


def latlon2a5(lat, lon, res):
    """
    Convert latitude and longitude to A5 cell identifier.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): A5 resolution level [0-29]

    Returns:
        str: A5 cell identifier in hexadecimal format

    Example:
        >>> latlon2a5(10.775275567242561, 106.70679737574993, 15)
        'a5c8b9a8b9a8b9a8'
    """
    if res is None:
        res = DGGS_TYPES["a5"]["default_res"]
    res = validate_a5_resolution(res)
    a5_id = a5.lonlat_to_cell((lon, lat), res)
    a5_hex = a5.u64_to_hex(a5_id)
    return a5_hex


def latlon2a5_cli():
    """
    Command-line interface for latlon2a5.
    """
    min_res = DGGS_TYPES["a5"]["min_res"]
    max_res = DGGS_TYPES["a5"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to A5 code at a specific Resolution [0..29]. \
                                     Usage: latlon2a5 <lat> <lon> <res> [0..29]. \
                                     Ex: latlon2a5 10.775275567242561 106.70679737574993 15"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["a5"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon

    a5_hex = latlon2a5(lat, lon, res)
    print(a5_hex)


def latlon2healpix(lat, lon, res=None):
    """
    Convert latitude and longitude to a HEALPix UNIQ cell ID.

    The UNIQ scheme packs resolution (order) and nested pixel index into a
    single integer, so downstream converters (e.g. healpix2geo) do not need a
    separate resolution argument.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): HEALPix resolution/order [0-29] (nside = 2^res)

    Returns:
        int: HEALPix UNIQ cell ID encoding both resolution and nested index

    Example:
        >>> latlon2healpix(10.775275567242561, 106.70679737574993, 10)
        9941583
    """
    if res is None:
        res = DGGS_TYPES["healpix"]["default_res"]
    res = validate_healpix_resolution(res)
    nside = order2nside(res)
    ipix = lonLat2PixNest(nside, lon, lat)
    return orderpix2uniq(res, ipix)


def latlon2healpix_cli():
    """
    Command-line interface for latlon2healpix.
    """
    min_res = DGGS_TYPES["healpix"]["min_res"]
    max_res = DGGS_TYPES["healpix"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to HEALPix UNIQ ID at a specific Resolution [0..29]. \
                                     Usage: latlon2healpix <lat> <lon> <res> [0..29]. \
                                     Ex: latlon2healpix 10.775275567242561 106.70679737574993 10"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["healpix"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    healpix_id = latlon2healpix(args.lat, args.lon, args.res)
    print(healpix_id)


def latlon2rhealpix(lat, lon, res):
    """
    Convert latitude and longitude to RHEALPix cell identifier.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): RHEALPix resolution level [0-15]

    Returns:
        str: RHEALPix cell identifier

    Example:
        >>> latlon2rhealpix(10.775275567242561, 106.70679737574993, 8)
        'N:1:8:1:2:3:4:5:6:7:8'
    """
    res = validate_rhealpix_resolution(res)
    E = WGS84_ELLIPSOID
    rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)
    point = (lon, lat)
    rhealpix_cell = rhealpix_dggs.cell_from_point(res, point, plane=False)
    rhealpix_id = str(rhealpix_cell)
    return rhealpix_id


def latlon2rhealpix_cli():
    """
    Command-line interface for latlon2rhealpix.
    """
    min_res = DGGS_TYPES["rhealpix"]["min_res"]
    max_res = DGGS_TYPES["rhealpix"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to Rhealpix code at a specific Resolution [0..15]. \
                                     Usage: latlon2rhealpix <lat> <lon> <res> [0..15]. \
                                     Ex: latlon2rhealpix 10.775275567242561 106.70679737574993 {DGGS_TYPES['rhealpix']['default_res']}"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["rhealpix"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon

    rhealpix_id = latlon2rhealpix(lat, lon, res)
    print(rhealpix_id)


def latlon2isea4t(lat, lon, res):
    """
    Convert latitude and longitude to ISEA4T cell identifier.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): ISEA4T resolution level [0-39] (Windows only)

    Returns:
        str: ISEA4T cell identifier

    Example:
        >>> latlon2isea4t(10.775275567242561, 106.70679737574993, 20)
        '0123456789abcdef012345'

    Note:
        This function is only available on Windows systems.
    """
    if res is None:
        res = DGGS_TYPES["isea4t"]["default_res"]
    res = validate_isea4t_resolution(res)
    max_accuracy = ISEA4T_RES_ACCURACY_DICT[
        39
    ]  # maximum cell_id length with 41 characters
    lat_long_point = LatLongPoint(lat, lon, max_accuracy)
    isea4t_cell_max_accuracy = isea4t_dggs.convert_point_to_dggs_cell(lat_long_point)
    cell_id_len = res + 2
    isea4t_cell = DggsCell(isea4t_cell_max_accuracy._cell_id[:cell_id_len])
    return isea4t_cell._cell_id


def latlon2isea4t_cli():
    """
    Command-line interface for latlon2isea4t.
    """
    min_res = DGGS_TYPES["isea4t"]["min_res"]
    max_res = DGGS_TYPES["isea4t"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to OpenEaggr ISEA4T code at a specific Resolution [0..39]. \
                                     Usage: latlon2isea4t <lat> <lon> <res> [0..39]. \
                                     Ex: latlon2isea4t 10.775275567242561 106.70679737574993 20"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["isea4t"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    if platform.system() == "Windows":
        isea4t_id = latlon2isea4t(lat, lon, res)
        print(isea4t_id)
    else:
        print("ISEA4T is only supported on Windows systems")


def latlon2isea3h(lat, lon, res):
    """
    Convert latitude and longitude to ISEA3H cell identifier.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): ISEA3H resolution level [0-40] (Windows only)

    Returns:
        str: ISEA3H cell identifier

    Example:
        >>> latlon2isea3h(10.775275567242561, 106.70679737574993, 27)
        '0123456789abcdef012345678'

    Note:
        This function is only available on Windows systems.
        Resolution 27 is suitable for geocoding applications.
    """
    # res: [0..40], res=27 is suitable for geocoding
    if res is None:
        res = DGGS_TYPES["isea3h"]["default_res"]
    res = validate_isea3h_resolution(res)
    accuracy = ISEA3H_RES_ACCURACY_DICT.get(res)
    lat_long_point = LatLongPoint(lat, lon, accuracy)
    isea3h_cell = isea3h_dggs.convert_point_to_dggs_cell(lat_long_point)
    return str(isea3h_cell.get_cell_id())


def latlon2isea3h_cli():
    """
    Command-line interface for latlon2isea3h.
    """
    min_res = DGGS_TYPES["isea3h"]["min_res"]
    max_res = DGGS_TYPES["isea3h"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to OpenEaggr ISEA3H code at a specific Resolution [0..40]. \
                                     Usage: latlon2isea3h <lat> <lon> <res> [0..40]. \
                                     Ex: latlon2isea3h 10.775275567242561 106.70679737574993 27"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["isea3h"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    if platform.system() == "Windows":
        isea3h_id = latlon2isea3h(lat, lon, res)
        print(isea3h_id)
    else:
        print("ISEA3H is only supported on Windows systems")


def latlon2dggrid(
    dggrid_instance, dggs_type, lat, lon, res, output_address_type="SEQNUM"
):
    """
    Convert latitude and longitude to DGGRID cell identifier.

    Args:
        dggrid_instance: DGGRID instance for processing
        dggs_type (str): DGGRID DGGS type (e.g., 'ISEA7H', 'ISEA4T')
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int, optional): Resolution level. Defaults to type-specific default.
        output_address_type (str): Output address format ('SEQNUM', 'INTERLEAVE', etc.)

    Returns:
        str: DGGRID cell identifier in specified format

    Example:
        >>> dggrid_instance = create_dggrid_instance()
        >>> latlon2dggrid(dggrid_instance, 'ISEA7H', 10.775275567242561, 106.70679737574993, 13)
        '123456789012345'
    """
    point = Point(lon, lat)
    dggs_type = validate_dggrid_type(dggs_type)
    if res is None:
        res = DGGRID_TYPES[dggs_type]["default_res"]
    res = validate_dggrid_resolution(dggs_type, res)
    geodf_points_wgs84 = gpd.GeoDataFrame([{"geometry": point}], crs="EPSG:4326")

    dggrid_cell = dggrid_instance.cells_for_geo_points(
        geodf_points_wgs84=geodf_points_wgs84,
        cell_ids_only=True,
        dggs_type=dggs_type,
        resolution=res,
        output_address_type=output_address_type,
    )
    # seqnum = dggrid_cell.loc[0, "name"]
    # if output_address_type == "SEQNUM":
    #     return seqnum
    # address_type_transform = dggrid_instance.address_transform(
    #     [seqnum],
    #     dggs_type=dggs_type,
    #     resolution=res,
    #     mixed_aperture_level=None,
    #     input_address_type="SEQNUM",
    #     output_address_type=output_address_type,
    # )
    # dggrid_id = address_type_transform.loc[0, output_address_type]
    # return dggrid_id
    return dggrid_cell.loc[0, "name"]


def latlon2dggrid_cli():
    """
    Command-line interface for latlon2dggrid.
    """
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to DGGRID cell at a specific Resolution. \
                                     Usage: latlon2dggrid <dggs_type> <lat> <lon> <res>. \
                                     Ex: latlon2dggrid ISEA7H 10.775275567242561 106.70679737574993 13"
    )
    parser.add_argument(
        "dggs_type",
        choices=DGGRID_TYPES.keys(),
        help="Select a DGGS type from the available options.",
    )

    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        help="Input Resolution (range depends on dggs_type)",
    )
    parser.add_argument(
        "output_address_type",
        choices=output_address_types,
        default="SEQNUM",
        nargs="?",  # This makes the argument optional
        help="Select an output address type from the available options.",
    )
    args = parser.parse_args()
    dggs_type = args.dggs_type
    res = args.res
    output_address_type = args.output_address_type
    dggrid_instance = create_dggrid_instance()
    dggrid_cell_id = latlon2dggrid(
        dggrid_instance, dggs_type, args.lat, args.lon, res, output_address_type
    )
    print(dggrid_cell_id)


def latlon2ease(lat, lon, res):
    """
    Convert latitude and longitude to EASE-DGGS cell identifier.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): EASE resolution level [0-6]

    Returns:
        str: EASE-DGGS cell identifier

    Example:
        >>> latlon2ease(10.775275567242561, 106.70679737574993, 3)
        '0123456789abcdef'
    """
    if res is None:
        res = DGGS_TYPES["ease"]["default_res"]
    res = validate_ease_resolution(res)
    ease_cell = geos_to_grid_ids([(lon, lat)], level=res)
    ease_id = ease_cell["result"]["data"][0]
    return ease_id


def latlon2ease_cli():
    """
    Command-line interface for latlon2isea3h.
    """
    min_res = DGGS_TYPES["ease"]["min_res"]
    max_res = DGGS_TYPES["ease"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to EASE-DGGS cell at a specific Resolution [0..6]. \
                                     Usage: latlon2ease <lat> <lon> <res> [0..6]. \
                                            Ex: latlon2ease 10.775275567242561 106.70679737574993 3"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["ease"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon

    ease_id = latlon2ease(lat, lon, res)
    print(ease_id)


def latlon2qtm(lat, lon, res):
    """
    Convert latitude and longitude to QTM cell identifier.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): QTM resolution level [1-24]

    Returns:
        str: QTM cell identifier

    Example:
        >>> latlon2qtm(10.775275567242561, 106.70679737574993, 12)
        '0123456789ab'
    """
    if res is None:
        res = DGGS_TYPES["qtm"]["default_res"]
    res = validate_qtm_resolution(res)
    qtm_id = qtm.latlon_to_qtm_id(lat, lon, res)
    return qtm_id


def latlon2qtm_cli():
    """
    Command-line interface for latlon2qtm.
    """
    min_res = DGGS_TYPES["qtm"]["min_res"]
    max_res = DGGS_TYPES["qtm"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to QTM. \
                                     Usage: latlon2qtm <lat> <lon> <res> [1..24]. \
                                     Ex: latlon2qtm 10.775275567242561 106.70679737574993 12"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["qtm"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon

    qtm_id = latlon2qtm(lat, lon, res)
    print(qtm_id)


def latlon2olc(lat, lon, res):
    """
    Convert latitude and longitude to Open Location Code (OLC/Plus Code).

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): OLC code length [2,4,6,8,10-15]

    Returns:
        str: Open Location Code

    Example:
        >>> latlon2olc(10.775275567242561, 106.70679737574993, 12)
        '6P5XQ2+2Q'
    """
    if res is None:
        res = DGGS_TYPES["olc"]["default_res"]
    res = validate_olc_resolution(res)
    olc_id = olc.encode(lat, lon, res)
    return olc_id


def latlon2olc_cli():
    """
    Command-line interface for latlon2olc.
    """
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to OLC/ Google Plus Code at a specific Code length [10..15]. \
                                     Usage: latlon2olc <lat> <lon> <res> [2,4,6,8,10..15]. \
                                     Ex: latlon2olc 10.775275567242561 106.70679737574993 12"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["olc"]["default_res"],
        choices=olc_resolutions,
        help=f"Resolution of the OLC DGGS (choose from {', '.join(map(str, olc_resolutions))})",
    )

    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon

    olc_id = latlon2olc(lat, lon, res)
    print(olc_id)


def latlon2geohash(lat, lon, res):
    """
    Convert latitude and longitude to Geohash.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): Geohash precision [1-10]

    Returns:
        str: Geohash string

    Example:
        >>> latlon2geohash(10.775275567242561, 106.70679737574993, 6)
        'w8k3x2'
    """
    if res is None:
        res = DGGS_TYPES["geohash"]["default_res"]
    res = validate_geohash_resolution(res)
    geohash_id = geohash.encode(lat, lon, res)
    return geohash_id


def latlon2geohash_cli():
    """
    Command-line interface for latlon2geohash.
    """
    min_res = DGGS_TYPES["geohash"]["min_res"]
    max_res = DGGS_TYPES["geohash"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to Geohash code at a specific resolution [1..10]. \
                                    Usage: latlon2geohash <lat> <lon> <res>[1..10]. \
                                    Ex: latlon2geohash 10.775275567242561 106.70679737574993 6"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["geohash"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    geohash_id = latlon2geohash(lat, lon, res)
    print(geohash_id)


def latlon2georef(lat, lon, res):
    """
    Convert latitude and longitude to GEOREF.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): GEOREF resolution [0-10]

    Returns:
        str: GEOREF code

    Example:
        >>> latlon2georef(10.775275567242561, 106.70679737574993, 5)
        'MK1234567890'
    """
    res = validate_georef_resolution(res)
    georef_id = georef.encode(lat, lon, res)
    return georef_id


def latlon2georef_cli():
    """
    Command-line interface for latlon2georef.
    """
    min_res = DGGS_TYPES["georef"]["min_res"]
    max_res = DGGS_TYPES["georef"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to GEOREF code at a specific resolution [0..10]. \
                                     Usage: latlon2georef <lat> <lon> <res> [0..10]. \
                                    Ex: latlon2georef 10.775275567242561 106.70679737574993 5"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["georef"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    georef_id = latlon2georef(lat, lon, res)
    print(georef_id)


def latlon2mgrs(lat, lon, res):
    """
    Convert latitude and longitude to MGRS.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): MGRS precision [0-5]

    Returns:
        str: MGRS coordinate

    Example:
        >>> latlon2mgrs(10.775275567242561, 106.70679737574993, 3)
        '48PXV123456'
    """
    if res is None:
        res = DGGS_TYPES["mgrs"]["default_res"]
    res = validate_mgrs_resolution(res)
    mgrs_cell = mgrs.toMgrs(lat, lon, res)
    return mgrs_cell


def latlon2mgrs_cli():
    """
    Command-line interface for latlon2mgrs.
    """
    min_res = DGGS_TYPES["mgrs"]["min_res"]
    max_res = DGGS_TYPES["mgrs"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to GEOREF code at a specific resolution [0..5]. \
                                     Usage: latlon2mgrs <lat> <lon> <res> [0..5]. \
                                     Ex: latlon2mgrs 10.775275567242561 106.70679737574993 3"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["mgrs"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution  [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    mgrs_id = latlon2mgrs(lat, lon, res)
    print(mgrs_id)


def latlon2tilecode(lat, lon, res):
    """
    Convert latitude and longitude to Tilecode.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): Tile zoom level [0-29]

    Returns:
        str: Tilecode identifier

    Example:
        >>> latlon2tilecode(10.775275567242561, 106.70679737574993, 15)
        '123456789012345'
    """
    if res is None:
        res = DGGS_TYPES["tilecode"]["default_res"]
    res = validate_tilecode_resolution(res)
    tilecode_id = tilecode.latlon2tilecode(lat, lon, res)
    return tilecode_id


def latlon2tilecode_cli():
    """
    Command-line interface for latlon2tilecode.
    """
    min_res = DGGS_TYPES["tilecode"]["min_res"]
    max_res = DGGS_TYPES["tilecode"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to Tile code at a specific resolution/ zoom level [0..29]. \
                                    Usage: latlon2tilecode <lat> <lon> <res> [0..29]. \
                                    Ex: latlon2tilecode 10.775275567242561 106.70679737574993 15"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["tilecode"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution/ Zoom level [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    tilecode_id = latlon2tilecode(lat, lon, res)
    print(tilecode_id)


def latlon2quadkey(lat, lon, res):
    """
    Convert latitude and longitude to Quadkey.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): Quadkey zoom level [0-29]

    Returns:
        str: Quadkey identifier

    Example:
        >>> latlon2quadkey(10.775275567242561, 106.70679737574993, 15)
        '123456789012345'
    """
    if res is None:
        res = DGGS_TYPES["quadkey"]["default_res"]
    res = validate_quadkey_resolution(res)
    quadkey_id = tilecode.latlon2quadkey(lat, lon, res)
    return quadkey_id


def latlon2quadkey_cli():
    """
    Command-line interface for latlon2tilecode.
    """
    min_res = DGGS_TYPES["quadkey"]["min_res"]
    max_res = DGGS_TYPES["quadkey"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to Quadkey at a specific resolution/ zoom level [0..29]. \
                                     Usage: latlon2quadkey <lat> <lon> <res> [0..29]. \
                                     Ex: latlon2quadkey 10.775275567242561 106.70679737574993 15"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["quadkey"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution/ Zoom level [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    quadkey_id = latlon2quadkey(lat, lon, res)
    print(quadkey_id)


def latlon2maidenhead(lat, lon, res):
    """
    Convert latitude and longitude to Maidenhead grid square.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): Maidenhead precision [1-4]

    Returns:
        str: Maidenhead grid square

    Example:
        >>> latlon2maidenhead(10.775275567242561, 106.70679737574993, 2)
        'OK12'
    """
    if res is None:
        res = DGGS_TYPES["maidenhead"]["default_res"]
    res = validate_maidenhead_resolution(res)
    maidenhead_id = maidenhead.toMaiden(lat, lon, res)
    # maidenhead_id = maidenhead.to_maiden(lat, lon, res)
    return maidenhead_id


def latlon2maidenhead_cli():
    """
    Command-line interface for latlon2maidenhead.
    """
    min_res = DGGS_TYPES["maidenhead"]["min_res"]
    max_res = DGGS_TYPES["maidenhead"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to Tile code at a specific resolution [1..4]. \
                                    Usage: latlon2maidenhead <lat> <lon> <res> [1..4]. \
                                    Ex: latlon2maidenhead 10.775275567242561 106.70679737574993 2"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["maidenhead"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    maidenhead_id = latlon2maidenhead(lat, lon, res)
    print(maidenhead_id)


def latlon2gars(lat, lon, res):
    """
    Convert latitude and longitude to GARS.

    Args:
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int): GARS resolution [1-4] (1=30min, 2=15min, 3=5min, 4=1min)

    Returns:
        str: GARS identifier

    Example:
        >>> latlon2gars(10.775275567242561, 106.70679737574993, 2)
        '123456789012345'
    """
    if res is None:
        res = DGGS_TYPES["gars"]["default_res"]
    res = validate_gars_resolution(res)
    # Convert res to minutes: 1->30, 2->15, 3->5, 4->1
    minutes_map = {1: 30, 2: 15, 3: 5, 4: 1}
    minutes = minutes_map[res]
    gars_cell = GARSGrid.from_latlon(lat, lon, minutes)
    return str(gars_cell)


def latlon2gars_cli():
    """
    Command-line interface for latlon2gars.
    """
    min_res = DGGS_TYPES["gars"]["min_res"]
    max_res = DGGS_TYPES["gars"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to GARS code at a specific resolution [1..4]. \
                                     Usage: latlon2gars <lat> <lon> <res> [1..4]. \
                                     Ex: latlon2gars 10.775275567242561 106.70679737574993 2"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["gars"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}] (1=30min, 2=15min, 3=5min, 4=1min)",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    gars_id = latlon2gars(lat, lon, res)
    print(gars_id)


def latlon2digipin(lat, lon, res=None):
    """
    Convert latitude and longitude to DIGIPIN code.

    Args:
        lat (float): Latitude in decimal degrees (must be between 2.5 and 38.5)
        lon (float): Longitude in decimal degrees (must be between 63.5 and 99.5)
        res (int): DIGIPIN resolution [1-15] (number of characters in the code)

    Returns:
        str: DIGIPIN identifier with dashes (e.g., 'F3K-492-6P96')

    Example:
        >>> latlon2digipin(28.6139, 77.2090, 10)
        'F3K-492-6P96'

    Note:
        DIGIPIN is a geocoding system for India. Coordinates outside the bounds
        (lat: 2.5-38.5, lon: 63.5-99.5) will return 'Out of Bound'.
    """
    if res is None:
        res = DGGS_TYPES["digipin"]["default_res"]
    res = validate_digipin_resolution(res)
    digipin_id = latlon_to_digipin(lat, lon, res)
    return digipin_id


def latlon2digipin_cli():
    """
    Command-line interface for latlon2digipin.
    """
    min_res = DGGS_TYPES["digipin"]["min_res"]
    max_res = DGGS_TYPES["digipin"]["max_res"]
    parser = argparse.ArgumentParser(
        description="Convert Lat, Long to DIGIPIN code at a specific resolution [1..10]. \
                                     Usage: latlon2digipin <lat> <lon> <res> [1..10]. \
                                     Ex: latlon2digipin 17.414718, 78.482992 10"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        default=DGGS_TYPES["digipin"]["default_res"],
        choices=range(min_res, max_res + 1),
        help=f"Input Resolution [{min_res}..{max_res}]",
    )
    args = parser.parse_args()

    res = args.res
    lat = args.lat
    lon = args.lon
    digipin_id = latlon2digipin(lat, lon, res)
    print(digipin_id)


def latlon2dggal(dggs_type, lat, lon, res):
    """
    Convert latitude and longitude to DGGAL zone identifier.

    Args:
        dggs_type (str): DGGAL DGGS type (e.g., 'gnosis', 'dggrid')
        lat (float): Latitude in decimal degrees
        lon (float): Longitude in decimal degrees
        res (int, optional): Resolution level. Defaults to type-specific default.

    Returns:
        str: DGGAL zone identifier

    Example:
        >>> latlon2dggal('gnosis', 10.775275567242561, 106.70679737574993, 8)
        '123456789012345'
    """
    dggs_type = validate_dggal_type(dggs_type)
    if res is None:
        res = DGGAL_TYPES[dggs_type]["default_res"]
    res = validate_dggal_resolution(dggs_type, res)
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = globals()[dggs_class_name]()
    dggal_zone = dggrs.getZoneFromWGS84Centroid(res, GeoPoint(lat, lon))
    dggal_zoneid = dggrs.getZoneTextID(dggal_zone)
    return dggal_zoneid


def latlon2dggal_cli():
    """
    Command-line interface for latlon2dggal.

    Usage: latlon2dggal <dggs_type> <lat> <lon> <res>
    """
    parser = argparse.ArgumentParser(
        description=(
            "Convert Lat, Long to DGGAL ZoneID via dgg CLI (zone). "
            "Usage: latlon2dggal <dggs_type> <lat> <lon> <res>. "
            "Ex: latlon2dggal gnosis 10.775275567242561 106.70679737574993 8"
        )
    )
    parser.add_argument(
        "dggs_type", type=str, choices=DGGAL_TYPES.keys(), help="DGGAL type"
    )
    parser.add_argument("lat", type=float, help="Input Latitude")
    parser.add_argument("lon", type=float, help="Input Longitude")
    parser.add_argument(
        "res",
        type=int,
        help="Input Resolution (range depends on dggs_type)",
    )

    args = parser.parse_args()

    # Default res=8
    zone_id = latlon2dggal(args.dggs_type, args.lat, args.lon, args.res)
    if zone_id is None:
        raise SystemExit(
            "Failed to compute DGGAL ZoneID. Ensure `dgg` is installed and on PATH."
        )
    print(zone_id)
