"""
DIGIPIN to Geometry Module

This module provides functionality to convert DIGIPIN identifiers to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    digipin2geo: Convert DIGIPIN IDs to Shapely Polygons
    digipin2geojson: Convert DIGIPIN IDs to GeoJSON FeatureCollection
    digipin2geo_cli: Command-line interface for polygon conversion
    digipin2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from shapely.geometry import Polygon
from vgrid.dggs.digipin import digipin_to_bounds
from vgrid.utils.geometry import graticule_dggs_to_feature


def digipin2geo(digipin_ids):
    """
    Convert DIGIPIN cell IDs to Shapely geometry objects.

    Accepts a single digipin_id (string) or a list of digipin_ids. For each valid DIGIPIN cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    digipin_ids : str or list of str
        DIGIPIN cell ID(s) to convert. Can be a single string or a list of strings.
        Format: Alphanumeric code with optional dashes (e.g., 'F3K-492-6P96' or 'F3K4926P96')
        DIGIPIN codes represent locations in India (lat: 2.5-38.5, lon: 63.5-99.5)

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single DIGIPIN cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding DIGIPIN cell.

    Examples
    --------
    >>> digipin2geo("F3K")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> digipin2geo(["F3K", "39J-438-TJC7"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(digipin_ids, str):
        digipin_ids = [digipin_ids]
    digipin_polygons = []
    for digipin_id in digipin_ids:
        try:
            bounds = digipin_to_bounds(digipin_id)
            if isinstance(bounds, str):  # Error message like 'Invalid DIGIPIN'
                continue
            min_lat = bounds["minLat"]
            max_lat = bounds["maxLat"]
            min_lon = bounds["minLon"]
            max_lon = bounds["maxLon"]
            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat],
                ]
            )
            digipin_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(digipin_polygons) == 1:
        return digipin_polygons[0]
    return digipin_polygons


def digipin2geo_cli():
    """
    Command-line interface for digipin2geo supporting multiple DIGIPIN codes.
    """
    parser = argparse.ArgumentParser(
        description="Convert DIGIPIN code(s) to Shapely Polygons"
    )
    parser.add_argument(
        "digipin_id", nargs="+", help="Input DIGIPIN code(s), e.g. F3K 39J-438-TJC7"
    )
    args = parser.parse_args()
    polys = digipin2geo(args.digipin_id)
    return polys


def digipin2geojson(digipin_ids):
    """
    Convert DIGIPIN cell IDs to GeoJSON FeatureCollection.

    Accepts a single digipin_id (string) or a list of digipin_ids. For each valid DIGIPIN cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    digipin_ids : str or list of str
        DIGIPIN cell ID(s) to convert. Can be a single string or a list of strings.
        Format: Alphanumeric code with optional dashes (e.g., 'F3K-492-6P96' or 'F3K4926P96')
        DIGIPIN codes represent locations in India (lat: 2.5-38.5, lon: 63.5-99.5)

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid DIGIPIN cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the DIGIPIN cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> digipin2geojson("F3K")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> digipin2geojson(["F3K", "39J-438-TJC7"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(digipin_ids, str):
        digipin_ids = [digipin_ids]
    digipin_features = []
    for digipin_id in digipin_ids:
        try:
            bounds = digipin_to_bounds(digipin_id)
            if isinstance(bounds, str):  # Error message like 'Invalid DIGIPIN'
                continue
            min_lat = bounds["minLat"]
            max_lat = bounds["maxLat"]
            min_lon = bounds["minLon"]
            max_lon = bounds["maxLon"]
            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat],
                ]
            )
            # Calculate resolution from DIGIPIN code length (excluding dashes)
            clean_id = digipin_id.replace("-", "")
            resolution = len(clean_id)
            digipin_feature = graticule_dggs_to_feature(
                "digipin_id", digipin_id, resolution, cell_polygon
            )
            digipin_features.append(digipin_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": digipin_features}


def digipin2geojson_cli():
    """
    Command-line interface for digipin2geojson supporting multiple DIGIPIN codes.
    """
    parser = argparse.ArgumentParser(description="Convert DIGIPIN code(s) to GeoJSON")
    parser.add_argument(
        "digipin_id", nargs="+", help="Input DIGIPIN code(s), e.g. F3K 39J-438-TJC7"
    )
    args = parser.parse_args()
    geojson_data = json.dumps(digipin2geojson(args.digipin_id))
    print(geojson_data)
