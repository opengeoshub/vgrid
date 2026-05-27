"""
GEOREF to Geometry Module

This module provides functionality to convert GEOREF (World Geographic Reference System) codes to Shapely Polygons and GeoJSON FeatureCollection.
GEOREF is a standardized system for identifying locations on the Earth's surface using a grid-based coordinate system.

Key Functions:
    georef2geo: Convert GEOREF codes to Shapely Polygons
    georef2geojson: Convert GEOREF codes to GeoJSON FeatureCollection
    georef2geo_cli: Command-line interface for polygon conversion
    georef2geojson_cli: Command-line interface for GeoJSON conversion
"""

from vgrid.dggs import georef
from shapely.geometry import Polygon
import json
import argparse
from vgrid.utils.geometry import graticule_dggs_to_feature


def georef2geo(georef_ids):
    """
    Convert GEOREF codes to Shapely geometry objects.

    Accepts a single georef_id (string) or a list of georef_ids. For each valid GEOREF code,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    georef_ids : str or list of str
        GEOREF code(s) to convert. Can be a single string or a list of strings.
        Example format: "VGBL42404651"

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single GEOREF code is provided, returns a single Shapely Polygon object.
        If a list of codes is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding GEOREF cell.

    Examples
    --------
    >>> georef2geo("VGBL42404651")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> georef2geo(["VGBL42404651", "VGBL42404652"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(georef_ids, str):
        georef_ids = [georef_ids]
    georef_polygons = []
    for georef_id in georef_ids:
        try:
            _, _, min_lat, min_lon, max_lat, max_lon, _ = georef.georefcell(georef_id)
            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat],
                ]
            )
            georef_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(georef_polygons) == 1:
        return georef_polygons[0]
    return georef_polygons


def georef2geo_cli():
    """
    Command-line interface for georef2geo supporting multiple GEOREF codes.
    """
    parser = argparse.ArgumentParser(
        description="Convert GEOREF code(s) to Shapely Polygons"
    )
    parser.add_argument(
        "georef",
        nargs="+",
        help="Input GEOREF code(s), e.g., georef2geo VGBL42404651 ...",
    )
    args = parser.parse_args()
    polys = georef2geo(args.georef)
    return polys


def georef2geojson(georef_ids):
    """
    Convert GEOREF codes to GeoJSON FeatureCollection.

    Accepts a single georef_id (string) or a list of georef_ids. For each valid GEOREF code,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    georef_ids : str or list of str
        GEOREF code(s) to convert. Can be a single string or a list of strings.
        Example format: "VGBL42404651"

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid GEOREF cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the GEOREF code, resolution level, and cell metadata

    Examples
    --------
    >>> georef2geojson("VGBL42404651")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> georef2geojson(["VGBL42404651", "VGBL42404652"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(georef_ids, str):
        georef_ids = [georef_ids]
    georef_features = []
    for georef_id in georef_ids:
        try:
            _, _, min_lat, min_lon, max_lat, max_lon, resolution = georef.georefcell(
                georef_id
            )
            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat],
                ]
            )
            georef_feature = graticule_dggs_to_feature(
                "georef", georef_id, resolution, cell_polygon
            )
            georef_features.append(georef_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": georef_features}


def georef2geojson_cli():
    """
    Command-line interface for georef2geojson supporting multiple GEOREF codes.
    """
    parser = argparse.ArgumentParser(description="Convert GEOREF code(s) to GeoJSON")
    parser.add_argument(
        "georef",
        nargs="+",
        help="Input GEOREF code(s), e.g., georef2geojson VGBL42404651 ...",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(georef2geojson(args.georef))
    print(geojson_data)
