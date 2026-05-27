"""
DGGAL to Geometry Module

This module provides functionality to convert DGGAL ZoneIDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    dggal2geo: Convert DGGAL ZoneIDs to Shapely Polygons
    dggal2geojson: Convert DGGAL ZoneIDs to GeoJSON FeatureCollection
    dggal2geo_cli: Command-line interface for polygon conversion
    dggal2geojson_cli: Command-line interface for GeoJSON conversion
"""

import argparse
import json

# Try to import dggal library, handle gracefully if import fails
from dggal import *
from vgrid.utils.geometry import geodesic_dggs_to_feature
from vgrid.utils.constants import DGGAL_TYPES
from vgrid.utils.geometry import dggal_to_geo
from vgrid.utils.io import validate_dggal_type
from vgrid.utils.antimeridian import fix_polygon

app = Application(appGlobals=globals())
pydggal_setup(app)


def dggal2geo(
    dggs_type: str, zone_ids: str, options: dict = {}, split_antimeridian=False
):
    """
    Convert DGGAL ZoneIDs to Shapely geometry objects.

    Accepts a single zone_id (string) or a list of zone_ids. For each valid DGGAL ZoneID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    dggs_type : str
        DGGAL DGGS type (e.g., "isea3h", "isea4t", "rhealpix").
    zone_ids : str or list of str
        DGGAL ZoneID(s) to convert. Can be a single string or a list of strings.
        Example format: "A4-0-A"
    options : dict, optional
        Additional options for the conversion process.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single DGGAL ZoneID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding DGGAL cell.

    Examples
    --------
    >>> dggal2geo("isea3h", "A4-0-A")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> dggal2geo("isea3h", ["A4-0-A", "A4-0-B"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(zone_ids, str):
        zone_ids = [zone_ids]
    zone_polygons = []
    for zone_id in zone_ids:
        try:
            zone_polygon = dggal_to_geo(dggs_type, zone_id, options)
            if zone_polygon and split_antimeridian:
                zone_polygon = fix_polygon(zone_polygon)
            if zone_polygon:
                zone_polygons.append(zone_polygon)
        except Exception:
            continue
    if len(zone_polygons) == 1:
        return zone_polygons[0]
    return zone_polygons


def dggal2geo_cli():
    parser = argparse.ArgumentParser(
        description=(
            "Convert DGGAL ZoneID to Shapely geometry. "
            "Usage: dggal2geo <dggs_type> <ZoneID> [ZoneID ...]"
        )
    )
    parser.add_argument(
        "dggs_type", type=str, choices=DGGAL_TYPES.keys(), help="DGGAL DGGS type"
    )
    parser.add_argument(
        "zone_id", nargs="+", help="ZoneIDs, e.g., dggal2geo isea3h A4-0-A A4-0-B"
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=True,
        help="Enable Antimeridian splitting",
    )
    args = parser.parse_args()
    polys = dggal2geo(
        args.dggs_type, args.zone_id, split_antimeridian=args.split_antimeridian
    )
    return polys


def dggal2geojson(
    dggs_type: str, zone_ids: str, options: dict = {}, split_antimeridian=False
):
    """
    Convert DGGAL ZoneIDs to GeoJSON FeatureCollection.

    Accepts a single zone_id (string) or a list of zone_ids. For each valid DGGAL ZoneID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    dggs_type : str
        DGGAL DGGS type (e.g., "isea3h", "isea4t", "rhealpix").
    zone_ids : str or list of str
        DGGAL ZoneID(s) to convert. Can be a single string or a list of strings.
        Example format: "A4-0-A"
    options : dict, optional
        Additional options for the conversion process.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid DGGAL cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the DGGAL ZoneID, resolution level, and cell metadata

    Examples
    --------
    >>> dggal2geojson("isea3h", "A4-0-A")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> dggal2geojson("isea3h", ["A4-0-A", "A4-0-B"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    dggs_type = validate_dggal_type(dggs_type)
    # Create the appropriate DGGS instance
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = globals()[dggs_class_name]()

    if isinstance(zone_ids, str):
        zone_ids = [zone_ids]
    zone_features = []

    for zone_id in zone_ids:
        try:
            zone = dggrs.getZoneFromTextID(zone_id)
            resolution = dggrs.getZoneLevel(zone)
            num_edges = dggrs.countZoneEdges(zone)
            cell_polygon = dggal2geo(
                dggs_type, zone_id, options, split_antimeridian=split_antimeridian
            )
            zone_feature = geodesic_dggs_to_feature(
                f"dggal_{dggs_type}", zone_id, resolution, cell_polygon, num_edges
            )
            zone_features.append(zone_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": zone_features}


def dggal2geojson_cli():
    """
    Command-line interface for converting DGGAL ZoneIDs to GeoJSON.

    This function provides a command-line interface that accepts multiple DGGAL
    ZoneIDs as command-line arguments and outputs the corresponding GeoJSON
    FeatureCollection as a JSON string to stdout.

    Usage:
        dggal2geojson isea3h A4-0-A A4-0-B

    Output:
        Prints a JSON string representing a GeoJSON FeatureCollection to stdout.

    Example:
        $ python -m vgrid.conversion.dggs2geo.dggal2geo isea3h A4-0-A
        {"type": "FeatureCollection", "features": [...]}

    Note:
        This function is designed to be called from the command line and will
        parse arguments using argparse. The GeoJSON output is formatted as a
        JSON string printed to stdout. Invalid cell IDs are silently skipped.
    """
    parser = argparse.ArgumentParser(description="Convert DGGAL ZoneID(s) to GeoJSON")
    parser.add_argument(
        "dggs_type", type=str, choices=DGGAL_TYPES.keys(), help="DGGAL DGGS type"
    )
    parser.add_argument(
        "zone_id",
        nargs="+",
        help="Input DGGAL ZoneID(s), e.g., dggal2geojson isea3h A4-0-A A4-0-B",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Enable Antimeridian splitting",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(
        dggal2geojson(
            args.dggs_type, args.zone_id, split_antimeridian=args.split_antimeridian
        )
    )
    print(geojson_data)  # print to stdout
