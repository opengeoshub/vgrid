"""
EASE-DGGS to Geometry Module

This module provides functionality to convert EASE-DGGS (Equal-Area Scalable Earth
Grid) cell IDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    ease2geo: Convert EASE-DGGS cell IDs to Shapely Polygons
    ease2geojson: Convert EASE-DGGS cell IDs to GeoJSON FeatureCollection
    ease2geo_cli: Command-line interface for polygon conversion
    ease2geojson_cli: Command-line interface for GeoJSON conversion

"""

import json
import argparse
from shapely.geometry import Polygon
from ease_dggs.constants import levels_specs
from ease_dggs.dggs.grid_addressing import grid_ids_to_geos
from vgrid.utils.geometry import geodesic_dggs_to_feature, get_ease_resolution


def ease2geo(ease_ids):
    """
    Convert EASE-DGGS codes to Shapely geometry objects.

    Accepts a single ease_id (string) or a list of ease_ids. For each valid EASE-DGGS code,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    ease_ids : str or list of str
        EASE-DGGS code(s) to convert. Can be a single string or a list of strings.
        Example format: "L4.165767.02.02.20.71"

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single EASE-DGGS code is provided, returns a single Shapely Polygon object.
        If a list of codes is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding EASE-DGGS cell.

    Examples
    --------
    >>> ease2geo("L4.165767.02.02.20.71")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> ease2geo(["L4.165767.02.02.20.71", "L4.165768.02.02.20.71"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(ease_ids, str):
        ease_ids = [ease_ids]
    ease_polygons = []
    for ease_id in ease_ids:
        try:
            level = int(ease_id[1])
            level_spec = levels_specs[level]
            n_row = level_spec["n_row"]
            n_col = level_spec["n_col"]
            geo = grid_ids_to_geos([ease_id])
            center_lon, center_lat = geo["result"]["data"][0]
            cell_min_lat = center_lat - (180 / (2 * n_row))
            cell_max_lat = center_lat + (180 / (2 * n_row))
            cell_min_lon = center_lon - (360 / (2 * n_col))
            cell_max_lon = center_lon + (360 / (2 * n_col))
            cell_polygon = Polygon(
                [
                    [cell_min_lon, cell_min_lat],
                    [cell_max_lon, cell_min_lat],
                    [cell_max_lon, cell_max_lat],
                    [cell_min_lon, cell_max_lat],
                    [cell_min_lon, cell_min_lat],
                ]
            )
            ease_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(ease_polygons) == 1:
        return ease_polygons[0]
    return ease_polygons


def ease2geo_cli():
    """
    Command-line interface for ease2geo supporting multiple EASE-DGGS codes.
    """
    parser = argparse.ArgumentParser(
        description="Convert EASE-DGGS code(s) to Shapely Polygons"
    )
    parser.add_argument(
        "ease",
        nargs="+",
        help="Input EASE-DGGS code(s), e.g., ease2geo L4.165767.02.02.20.71 ...",
    )
    args = parser.parse_args()
    polys = ease2geo(args.ease)
    return polys


def ease2geojson(ease_ids):
    """
    Convert a list of EASE-DGGS codes to GeoJSON FeatureCollection.

    Accepts a single ease_id (string) or a list of ease_ids. For each valid EASE-DGGS code,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    ease_ids : str or list of str
        EASE-DGGS code(s) to convert. Can be a single string or a list of strings.
        Example format: "L4.165767.02.02.20.71"

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid EASE-DGGS cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the EASE-DGGS code, resolution level, and cell metadata

    Examples
    --------
    >>> ease2geojson("L4.165767.02.02.20.71")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> ease2geojson(["L4.165767.02.02.20.71", "L4.165768.02.02.20.71"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(ease_ids, str):
        ease_ids = [ease_ids]
    ease_features = []
    for ease_id in ease_ids:
        try:
            cell_polygon = ease2geo(ease_id)
            resolution = get_ease_resolution(ease_id)
            num_edges = 4
            ease_feature = geodesic_dggs_to_feature(
                "ease", ease_id, resolution, cell_polygon, num_edges
            )
            ease_features.append(ease_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": ease_features}


def ease2geojson_cli():
    """
    Command-line interface for ease2geojson supporting multiple EASE-DGGS codes.
    """
    parser = argparse.ArgumentParser(description="Convert EASE-DGGS code(s) to GeoJSON")
    parser.add_argument(
        "ease",
        nargs="+",
        help="Input EASE-DGGS code(s), e.g., ease2geojson L4.165767.02.02.20.71 ...",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(ease2geojson(args.ease))
    print(geojson_data)
