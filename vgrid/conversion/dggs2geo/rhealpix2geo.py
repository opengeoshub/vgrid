"""
RHEALPix to Geometry Module

This module provides functionality to convert RHEALPix (Rectified HEALPix) cell IDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    rhealpix2geo: Convert RHEALPix cell IDs to Shapely Polygons
    rhealpix2geojson: Convert RHEALPix cell IDs to GeoJSON FeatureCollection
    rhealpix2geo_cli: Command-line interface for polygon conversion
    rhealpix2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.dggs.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
from vgrid.utils.geometry import geodesic_dggs_to_feature, rhealpix_cell_to_polygon
from vgrid.utils.geometry import shift_balanced, shift_west, shift_east
from vgrid.utils.antimeridian import fix_polygon
from pyproj import Geod

geod = Geod(ellps="WGS84")
E = WGS84_ELLIPSOID
rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3)


def rhealpix2geo(rhealpix_ids, fix_antimeridian=None):
    """
    Convert RHEALPix cell IDs to Shapely geometry objects.

    Accepts a single rhealpix_id (string) or a list of rhealpix_ids. For each valid RHEALPix cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    rhealpix_ids : str or list of str
        RHEALPix cell ID(s) to convert. Can be a single string or a list of strings.
        Each ID should be a string starting with 'R' followed by numeric digits.
        Example format: "R31260335553825"
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single RHEALPix cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding RHEALPix cell.

    Examples
    --------
    >>> rhealpix2geo("R31260335553825")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> rhealpix2geo(["R31260335553825", "R31260335553826"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(rhealpix_ids, str):
        rhealpix_ids = [rhealpix_ids]
    rhealpix_polygons = []
    for rhealpix_id in rhealpix_ids:
        try:
            rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
            rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
            cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
            if fix_antimeridian == "shift" or fix_antimeridian == "shift_balanced":
                cell_polygon = shift_balanced(
                    cell_polygon, threshold_west=-149, threshold_east=149
                )
            elif fix_antimeridian == "shift_west":
                cell_polygon = shift_west(cell_polygon, threshold=-149)
            elif fix_antimeridian == "shift_east":
                cell_polygon = shift_east(cell_polygon, threshold=149)
            elif fix_antimeridian == "split":
                cell_polygon = fix_polygon(cell_polygon)
            rhealpix_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(rhealpix_polygons) == 1:
        return rhealpix_polygons[0]
    return rhealpix_polygons


def rhealpix2geo_cli():
    """
    Command-line interface for converting RHEALPix cell IDs to Shapely Polygons.

    This function provides a command-line interface that accepts multiple RHEALPix
    cell IDs as command-line arguments and returns the corresponding Shapely
    Polygon objects.

    Returns:
        list: A list of Shapely Polygon objects representing the converted cells.

    Usage:
        rhealpix2geo R31260335553825 R31260335553826

    Note:
        This function is designed to be called from the command line and will
        parse arguments using argparse. Invalid cell IDs are silently skipped.
    """
    parser = argparse.ArgumentParser(
        description="Convert Rhealpix cell ID(s) to Shapely Polygons"
    )
    parser.add_argument(
        "rhealpix",
        nargs="+",
        help="Input Rhealpix cell ID(s), e.g., rhealpix2geo R31260335553825 R31260335553826",
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=[
            "shift",
            "shift_balanced",
            "shift_west",
            "shift_east",
            "split",
            "none",
        ],
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    args = parser.parse_args()
    polys = rhealpix2geo(args.rhealpix, fix_antimeridian=args.fix_antimeridian)
    return polys


def rhealpix2geojson(rhealpix_ids, fix_antimeridian=None):
    """
    Convert RHEALPix cell IDs to GeoJSON FeatureCollection.

    Accepts a single rhealpix_id (string) or a list of rhealpix_ids. For each valid RHEALPix cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    rhealpix_ids : str or list of str
        RHEALPix cell ID(s) to convert. Can be a single string or a list of strings.
        Each ID should be a string starting with 'R' followed by numeric digits.
        Example format: "R31260335553825"
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid RHEALPix cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the RHEALPix cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> rhealpix2geojson("R31260335553825")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> rhealpix2geojson(["R31260335553825", "R31260335553826"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(rhealpix_ids, str):
        rhealpix_ids = [rhealpix_ids]
    rhealpix_features = []
    for rhealpix_id in rhealpix_ids:
        try:
            rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
            rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
            resolution = rhealpix_cell.resolution
            cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
            if fix_antimeridian == "shift" or fix_antimeridian == "shift_balanced":
                cell_polygon = shift_balanced(
                    cell_polygon, threshold_west=-128, threshold_east=160
                )
            elif fix_antimeridian == "shift_west":
                cell_polygon = shift_west(cell_polygon, threshold=-128)
            elif fix_antimeridian == "shift_east":
                cell_polygon = shift_east(cell_polygon, threshold=160)
            elif fix_antimeridian == "split":
                cell_polygon = fix_polygon(cell_polygon)
            num_edges = 4
            if rhealpix_cell.ellipsoidal_shape() == "dart":
                num_edges = 3
            rhealpix_feature = geodesic_dggs_to_feature(
                "rhealpix", rhealpix_id, resolution, cell_polygon, num_edges
            )
            rhealpix_features.append(rhealpix_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": rhealpix_features}


def rhealpix2geojson_cli():
    """
    Command-line interface for converting RHEALPix cell IDs to GeoJSON.

    This function provides a command-line interface that accepts multiple RHEALPix
    cell IDs as command-line arguments and outputs the corresponding GeoJSON
    FeatureCollection as a JSON string to stdout.

    Usage:
        rhealpix2geojson R31260335553825 R31260335553826

    Output:
        Prints a JSON string representing a GeoJSON FeatureCollection to stdout.

    Example:
        $ python -m vgrid.conversion.dggs2geo.rhealpix2geo R31260335553825
        {"type": "FeatureCollection", "features": [...]}

    Note:
        This function is designed to be called from the command line and will
        parse arguments using argparse. The GeoJSON output is formatted as a
        JSON string printed to stdout. Invalid cell IDs are silently skipped.
    """
    parser = argparse.ArgumentParser(
        description="Convert Rhealpix cell ID(s) to GeoJSON"
    )
    parser.add_argument(
        "rhealpix",
        nargs="+",
        help="Input Rhealpix cell ID(s), e.g., rhealpix2geojson R31260335553825 R31260335553826",
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=[
            "shift",
            "shift_balanced",
            "shift_west",
            "shift_east",
            "split",
            "none",
        ],
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(
        rhealpix2geojson(args.rhealpix, fix_antimeridian=args.fix_antimeridian)
    )
    print(geojson_data)
