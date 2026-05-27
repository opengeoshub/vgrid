"""
QTM to Geometry Module

This module provides functionality to convert QTM (Ternary Triangular Mesh) cell IDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    qtm2geo: Convert QTM cell IDs to Shapely Polygons
    qtm2geojson: Convert QTM cell IDs to GeoJSON FeatureCollection
    qtm2geo_cli: Command-line interface for polygon conversion
    qtm2geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from vgrid.dggs.qtm import constructGeometry, qtm_id_to_facet
from vgrid.utils.geometry import geodesic_dggs_to_feature


def qtm2geo(qtm_ids):
    """
    Convert QTM cell IDs to Shapely geometry objects.

    Accepts a single qtm_id (string) or a list of qtm_ids. For each valid QTM cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    qtm_ids : str or list of str
        QTM cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "42012321"

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single QTM cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding QTM cell.

    Examples
    --------
    >>> qtm2geo("42012321")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> qtm2geo(["42012321", "42012322"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(qtm_ids, str):
        qtm_ids = [qtm_ids]
    qtm_polygons = []
    for qtm_id in qtm_ids:
        try:
            facet = qtm_id_to_facet(qtm_id)
            cell_polygon = constructGeometry(facet)
            qtm_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(qtm_polygons) == 1:
        return qtm_polygons[0]
    return qtm_polygons


def qtm2geo_cli():
    """
    Command-line interface for qtm2geo supporting multiple QTM cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert QTM cell ID(s) to Shapely Polygons"
    )
    parser.add_argument(
        "qtm",
        nargs="+",
        help="Input QTM cell ID(s), e.g., qtm2geo 42012321 42012322",
    )
    args = parser.parse_args()
    polys = qtm2geo(args.qtm)
    return polys


def qtm2geojson(qtm_ids):
    """
    Convert QTM cell IDs to GeoJSON FeatureCollection.

    Accepts a single qtm_id (string) or a list of qtm_ids. For each valid QTM cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    qtm_ids : str or list of str
        QTM cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "42012321"

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid QTM cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the QTM cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> qtm2geojson("42012321")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> qtm2geojson(["42012321", "42012322"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(qtm_ids, str):
        qtm_ids = [qtm_ids]
    qtm_features = []
    for qtm_id in qtm_ids:
        try:
            cell_polygon = qtm2geo(qtm_id)
            resolution = len(qtm_id)
            num_edges = 3
            qtm_feature = geodesic_dggs_to_feature(
                "qtm", qtm_id, resolution, cell_polygon, num_edges
            )
            qtm_features.append(qtm_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": qtm_features}


def qtm2geojson_cli():
    """
    Command-line interface for qtm2geojson supporting multiple QTM cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert QTM cell ID(s) to GeoJSON")
    parser.add_argument(
        "qtm",
        nargs="+",
        help="Input QTM cell ID(s), e.g., qtm2geojson 42012321 42012322",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(qtm2geojson(args.qtm))
    print(geojson_data)
