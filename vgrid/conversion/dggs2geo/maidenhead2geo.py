"""
Maidenhead to Geometry Module

This module provides functionality to convert Maidenhead locator grid cell IDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    maidenhead2geo: Convert Maidenhead cell IDs to Shapely Polygons
    maidenhead2geojson: Convert Maidenhead cell IDs to GeoJSON FeatureCollection
    maidenhead2geo_cli: Command-line interface for polygon conversion
    maidenhead2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from shapely.geometry import Polygon
from vgrid.utils.geometry import graticule_dggs_to_feature
from vgrid.dggs import maidenhead


def maidenhead2geo(maidenhead_ids):
    """
    Convert Maidenhead cell IDs to Shapely geometry objects.

    Accepts a single maidenhead_id (string) or a list of maidenhead_ids. For each valid Maidenhead cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    maidenhead_ids : str or list of str
        Maidenhead cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "OK3046."

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single Maidenhead cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding Maidenhead cell.

    Examples
    --------
    >>> maidenhead2geo("OK3046.")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> maidenhead2geo(["OK3046.", "OK3047."])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(maidenhead_ids, str):
        maidenhead_ids = [maidenhead_ids]
    maidenhead_polygons = []
    for maidenhead_id in maidenhead_ids:
        try:
            _, _, min_lat, min_lon, max_lat, max_lon, _ = maidenhead.maidenGrid(
                maidenhead_id
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
            maidenhead_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(maidenhead_polygons) == 1:
        return maidenhead_polygons[0]
    return maidenhead_polygons


def maidenhead2geo_cli():
    """
    Command-line interface for maidenhead2geo supporting multiple Maidenhead cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert Maidenhead cell ID(s) to Shapely Polygons"
    )
    parser.add_argument(
        "maidenhead",
        nargs="+",
        help="Input Maidenhead cell ID(s), e.g., maidenhead2geo OK3046.",
    )
    args = parser.parse_args()
    polys = maidenhead2geo(args.maidenhead)
    return polys


def maidenhead2geojson(maidenhead_ids):
    """
    Convert Maidenhead cell IDs to GeoJSON FeatureCollection.

    Accepts a single maidenhead_id (string) or a list of maidenhead_ids. For each valid Maidenhead cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    maidenhead_ids : str or list of str
        Maidenhead cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "OK3046."

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid Maidenhead cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the Maidenhead cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> maidenhead2geojson("OK3046.")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> maidenhead2geojson(["OK3046.", "OK3047."])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(maidenhead_ids, str):
        maidenhead_ids = [maidenhead_ids]
    maidenhead_features = []
    for maidenhead_id in maidenhead_ids:
        try:
            cell_polygon = maidenhead2geo(maidenhead_id)
            resolution = int(len(maidenhead_id) / 2)
            maidenhead_feature = graticule_dggs_to_feature(
                "maidenhead", maidenhead_id, resolution, cell_polygon
            )
            maidenhead_features.append(maidenhead_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": maidenhead_features}


def maidenhead2geojson_cli():
    """
    Command-line interface for maidenhead2geojson supporting multiple Maidenhead cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert Maidenhead cell ID(s) to GeoJSON"
    )
    parser.add_argument(
        "maidenhead",
        nargs="+",
        help="Input Maidenhead cell ID(s), e.g., maidenhead2geojson OK3046.",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(maidenhead2geojson(args.maidenhead))
    print(geojson_data)
