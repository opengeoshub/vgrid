"""
GARS to Geometry Module

This module provides functionality to convert GARS (Global Area Reference System) cell IDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    gars2geo: Convert GARS cell IDs to Shapely Polygons
    gars2geojson: Convert GARS cell IDs to GeoJSON FeatureCollection
    gars2geo_cli: Command-line interface for polygon conversion
    gars2geojson_cli: Command-line interface for GeoJSON conversion
"""

from gars_field import garsgrid
from shapely.geometry import Polygon
import json
import argparse
from vgrid.utils.geometry import graticule_dggs_to_feature
from pyproj import Geod

geod = Geod(ellps="WGS84")


def gars2geo(gars_ids):
    """
    Convert GARS cell IDs to Shapely geometry objects.

    Accepts a single gars_id (string) or a list of gars_ids. For each valid GARS cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    gars_ids : str or list of str
        GARS cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "574JK1918"

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single GARS cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding GARS cell.

    Examples
    --------
    >>> gars2geo("574JK1918")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> gars2geo(["574JK1918", "574JK1919"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(gars_ids, str):
        gars_ids = [gars_ids]
    gars_polygons = []
    for gars_id in gars_ids:
        try:
            gars_grid = garsgrid.GARSGrid(gars_id)
            wkt_polygon = gars_grid.polygon
            cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
            gars_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(gars_polygons) == 1:
        return gars_polygons[0]
    return gars_polygons


def gars2geo_cli():
    """
    Command-line interface for gars2geo supporting multiple GARS cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert GARS cell ID(s) to Shapely Polygons"
    )
    parser.add_argument(
        "gars",
        nargs="+",
        help="Input GARS cell ID(s), e.g., gars2geo 574JK1918 ...",
    )
    args = parser.parse_args()
    polys = gars2geo(args.gars)
    return polys


def gars2geojson(gars_ids):
    """
    Convert GARS cell IDs to GeoJSON FeatureCollection.

    Accepts a single gars_id (string) or a list of gars_ids. For each valid GARS cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    gars_ids : str or list of str
        GARS cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "574JK1918"

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid GARS cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the GARS cell ID, resolution level, and cell metadata

    Examples
    --------
    >>> gars2geojson("574JK1918")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> gars2geojson(["574JK1918", "574JK1919"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(gars_ids, str):
        gars_ids = [gars_ids]
    gars_features = []
    for gars_id in gars_ids:
        try:
            gars_grid = garsgrid.GARSGrid(gars_id)
            wkt_polygon = gars_grid.polygon
            resolution_minute = gars_grid.resolution
            resolution = 1
            if resolution_minute == 30:
                resolution = 1
            elif resolution_minute == 15:
                resolution = 2
            elif resolution_minute == 5:
                resolution = 3
            elif resolution_minute == 1:
                resolution = 4
            cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
            gars_feature = graticule_dggs_to_feature(
                "gars", gars_id, resolution, cell_polygon
            )
            gars_features.append(gars_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": gars_features}


def gars2geojson_cli():
    """
    Command-line interface for gars2geojson supporting multiple GARS cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert GARS cell ID(s) to GeoJSON")
    parser.add_argument(
        "gars",
        nargs="+",
        help="Input GARS cell ID(s), e.g., gars2geojson 574JK1918 ...",
    )
    args = parser.parse_args()
    geojson_data = json.dumps(gars2geojson(args.gars))
    print(geojson_data)
