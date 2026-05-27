"""
Geohash to Geometry Module

This module provides functionality to convert Geohash cell identifiers to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    geohash2geo: Convert Geohash IDs to Shapely Polygons
    geohash2geojson: Convert Geohash IDs to GeoJSON FeatureCollection
    geohash2geo_cli: Command-line interface for polygon conversion
    geohash2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from shapely.geometry import Polygon
from vgrid.dggs import geohash
from vgrid.utils.geometry import graticule_dggs_to_feature


def geohash2geo(geohash_ids):
    """
    Convert Geohash cell IDs to Shapely geometry objects.

    Accepts a single geohash_id (string) or a list of geohash_ids. For each valid Geohash cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    geohash_ids : str or list of str
        Geohash cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "w3gvk1td8"

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single Geohash cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding Geohash cell.

    Examples
    --------
    >>> geohash2geo("w3gvk1td8")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> geohash2geo(["w3gvk1td8", "w3gvk1td9"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(geohash_ids, str):
        geohash_ids = [geohash_ids]
    geohash_polygons = []
    for geohash_id in geohash_ids:
        try:
            bbox = geohash.bbox(geohash_id)
            min_lat, min_lon = bbox["s"], bbox["w"]
            max_lat, max_lon = bbox["n"], bbox["e"]
            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat],
                ]
            )
            geohash_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(geohash_polygons) == 1:
        return geohash_polygons[0]
    return geohash_polygons


def geohash2geo_cli():
    """
    Command-line interface for geohash2geo supporting multiple Geohash cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert Geohash cell ID(s) to Shapely Polygons"
    )
    parser.add_argument(
        "geohash",
        nargs="+",
        help="Input Geohash cell ID(s), e.g., geohash2geo w3gvk1td8 ...",
    )
    args = parser.parse_args()
    polys = geohash2geo(args.geohash)
    return polys


def geohash2geojson(geohash_ids):
    """
    Convert Geohash cell IDs to GeoJSON FeatureCollection.

    Accepts a single geohash_id (string) or a list of geohash_ids. For each valid Geohash cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    geohash_ids : str or list of str
        Geohash cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "w3gvk1td8"

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid Geohash cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the Geohash cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> geohash2geojson("w3gvk1td8")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> geohash2geojson(["w3gvk1td8", "w3gvk1td9"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(geohash_ids, str):
        geohash_ids = [geohash_ids]
    geohash_features = []
    for geohash_id in geohash_ids:
        try:
            cell_polygon = geohash2geo(geohash_id)
            resolution = len(geohash_id)
            geohash_feature = graticule_dggs_to_feature(
                "geohash", geohash_id, resolution, cell_polygon
            )
            geohash_features.append(geohash_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": geohash_features}


def geohash2geojson_cli():
    """
    Command-line interface for geohash2geojson supporting multiple Geohash cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert Geohash cell ID(s) to GeoJSON"
    )
    parser.add_argument(
        "geohash",
        nargs="+",
        help="Input Geohash cell ID(s), e.g., geohash2geojson w3gvk1td8 ...",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(geohash2geojson(args.geohash))
    print(geojson_data)
