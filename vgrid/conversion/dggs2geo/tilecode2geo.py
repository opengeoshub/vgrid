"""
Tilecode to Geometry Module

This module provides functionality to convert Tilecode identifiers to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    tilecode2geo: Convert Tilecode IDs to Shapely Polygons
    tilecode2geojson: Convert Tilecode IDs to GeoJSON FeatureCollection
    tilecode2geo_cli: Command-line interface for polygon conversion
    tilecode2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import re
import argparse
from shapely.geometry import Polygon
from vgrid.dggs import mercantile
from vgrid.utils.geometry import graticule_dggs_to_feature


def tilecode2geo(tilecode_ids):
    """
    Convert Tilecode cell IDs to Shapely geometry objects.

    Accepts a single tilecode_id (string) or a list of tilecode_ids. For each valid Tilecode cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    tilecode_ids : str or list of str
        Tilecode cell ID(s) to convert. Can be a single string or a list of strings.
        Format: 'z{x}x{y}y{z}' where z is zoom level and x,y are tile coordinates.
        Example format: "z0x0y0"

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single Tilecode cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding Tilecode cell.

    Examples
    --------
    >>> tilecode2geo("z0x0y0")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> tilecode2geo(["z0x0y0", "z1x1y1"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(tilecode_ids, str):
        tilecode_ids = [tilecode_ids]
    tilecode_polygons = []
    for tilecode_id in tilecode_ids:
        try:
            match = re.match(r"z(\d+)x(\d+)y(\d+)", tilecode_id)
            if not match:
                continue
            z = int(match.group(1))
            x = int(match.group(2))
            y = int(match.group(3))
            bounds = mercantile.bounds(x, y, z)
            min_lat, min_lon = bounds.south, bounds.west
            max_lat, max_lon = bounds.north, bounds.east
            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat],
                ]
            )
            tilecode_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(tilecode_polygons) == 1:
        return tilecode_polygons[0]
    return tilecode_polygons


def tilecode2geo_cli():
    """
    Command-line interface for tilecode2geo supporting multiple Tilecodes.
    """
    parser = argparse.ArgumentParser(
        description="Convert Tilecode(s) to Shapely Polygons"
    )
    parser.add_argument(
        "tilecode_id", nargs="+", help="Input Tilecode(s), e.g. z0x0y0 z1x1y1"
    )
    args = parser.parse_args()
    polys = tilecode2geo(args.tilecode_id)
    return polys


def tilecode2geojson(tilecode_ids):
    """
    Convert Tilecode cell IDs to GeoJSON FeatureCollection.

    Accepts a single tilecode_id (string) or a list of tilecode_ids. For each valid Tilecode cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    tilecode_ids : str or list of str
        Tilecode cell ID(s) to convert. Can be a single string or a list of strings.
        Format: 'z{x}x{y}y{z}' where z is zoom level and x,y are tile coordinates.
        Example format: "z0x0y0"

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid Tilecode cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the Tilecode cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> tilecode2geojson("z0x0y0")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> tilecode2geojson(["z0x0y0", "z1x1y1"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(tilecode_ids, str):
        tilecode_ids = [tilecode_ids]
    tilecode_features = []
    for tilecode_id in tilecode_ids:
        try:
            match = re.match(r"z(\d+)x(\d+)y(\d+)", tilecode_id)
            if not match:
                continue
            z = int(match.group(1))
            x = int(match.group(2))
            y = int(match.group(3))
            bounds = mercantile.bounds(x, y, z)
            cell_polygon = tilecode2geo(tilecode_id)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                cell_polygon = Polygon(
                    [
                        [min_lon, min_lat],
                        [max_lon, min_lat],
                        [max_lon, max_lat],
                        [min_lon, max_lat],
                        [min_lon, min_lat],
                    ]
                )
                resolution = z
                tilecode_feature = graticule_dggs_to_feature(
                    "tilecode_id", tilecode_id, resolution, cell_polygon
                )
                tilecode_features.append(tilecode_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": tilecode_features}


def tilecode2geojson_cli():
    """
    Command-line interface for tilecode2geojson supporting multiple Tilecodes.
    """
    parser = argparse.ArgumentParser(description="Convert Tilecode(s) to GeoJSON")
    parser.add_argument(
        "tilecode_id", nargs="+", help="Input Tilecode(s), e.g. z0x0y0 z1x1y1"
    )
    args = parser.parse_args()
    geojson_data = json.dumps(tilecode2geojson(args.tilecode_id))
    print(geojson_data)
