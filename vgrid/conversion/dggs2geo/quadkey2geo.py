"""
Quadkey to Geometry Module

This module provides functionality to convert Quadkey identifiers to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    quadkey2geo: Convert Quadkey IDs to Shapely Polygons
    quadkey2geojson: Convert Quadkey IDs to GeoJSON FeatureCollection
    quadkey2geo_cli: Command-line interface for polygon conversion
    quadkey2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from shapely.geometry import Polygon
from vgrid.dggs import mercantile
from vgrid.utils.geometry import graticule_dggs_to_feature


def quadkey2geo(quadkey_ids):
    """
    Convert Quadkey cell IDs to Shapely geometry objects.

    Accepts a single quadkey_id (string) or a list of quadkey_ids. For each valid Quadkey cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    quadkey_ids : str or list of str
        Quadkey cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "13223011131020220011133"

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single Quadkey cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding Quadkey cell.

    Examples
    --------
    >>> quadkey2geo("13223011131020220011133")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> quadkey2geo(["13223011131020220011133", "13223011131020220011134"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(quadkey_ids, str):
        quadkey_ids = [quadkey_ids]
    quadkey_polygons = []
    for quadkey_id in quadkey_ids:
        try:
            tile = mercantile.quadkey_to_tile(quadkey_id)
            z = tile.z
            x = tile.x
            y = tile.y
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
            quadkey_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(quadkey_polygons) == 1:
        return quadkey_polygons[0]
    return quadkey_polygons


def quadkey2geo_cli():
    """
    Command-line interface for quadkey2geo supporting multiple Quadkeys.
    """
    parser = argparse.ArgumentParser(
        description="Convert Quadkey(s) to Shapely Polygons"
    )
    parser.add_argument(
        "quadkey", nargs="+", help="Input Quadkey(s), e.g. 13223011131020220011133 ..."
    )
    args = parser.parse_args()
    polys = quadkey2geo(args.quadkey)
    return polys


def quadkey2geojson(quadkey_ids):
    """
    Convert Quadkey cell IDs to GeoJSON FeatureCollection.

    Accepts a single quadkey_id (string) or a list of quadkey_ids. For each valid Quadkey cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    quadkey_ids : str or list of str
        Quadkey cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "13223011131020220011133"

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid Quadkey cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the Quadkey cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> quadkey2geojson("13223011131020220011133")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> quadkey2geojson(["13223011131020220011133", "13223011131020220011134"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(quadkey_ids, str):
        quadkey_ids = [quadkey_ids]
    quadkey_features = []
    for quadkey_id in quadkey_ids:
        try:
            tile = mercantile.quadkey_to_tile(quadkey_id)
            z = tile.z
            x = tile.x
            y = tile.y
            bounds = mercantile.bounds(x, y, z)
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
                quadkey_feature = graticule_dggs_to_feature(
                    "quadkey", quadkey_id, resolution, cell_polygon
                )
                quadkey_features.append(quadkey_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": quadkey_features}


def quadkey2geojson_cli():
    """
    Command-line interface for quadkey2geojson supporting multiple Quadkeys.
    """
    parser = argparse.ArgumentParser(description="Convert Quadkey(s) to GeoJSON")
    parser.add_argument(
        "quadkey", nargs="+", help="Input Quadkey(s), e.g. 13223011131020220011133 ..."
    )
    args = parser.parse_args()
    geojson_data = json.dumps(quadkey2geojson(args.quadkey))
    print(geojson_data)
