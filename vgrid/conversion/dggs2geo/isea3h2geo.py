"""
ISEA3H to Geometry Module

This module provides functionality to convert ISEA3H (Icosahedral Snyder Equal Area Aperture 3 Hexagon) cell IDs to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    isea3h2geo: Convert ISEA3H cell IDs to Shapely Polygons
    isea3h2geojson: Convert ISEA3H cell IDs to GeoJSON FeatureCollection
    isea3h2geo_cli: Command-line interface for polygon conversion
    isea3h2geojson_cli: Command-line interface for GeoJSON conversion

Note: This module is only supported on Windows systems due to OpenEaggr dependency.
"""

import json
import argparse
from shapely.geometry import mapping
import platform

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.dggs.eaggr.enums.model import Model
    from vgrid.utils.constants import ISEA3H_ACCURACY_RES_DICT

    isea3h_dggs = Eaggr(Model.ISEA3H)

from pyproj import Geod
from vgrid.utils.geometry import (
    isea3h_cell_to_polygon,
    shift_balanced,
    shift_west,
    shift_east,
)

geod = Geod(ellps="WGS84")
from vgrid.utils.antimeridian import fix_polygon


def isea3h2geo(isea3h_ids, fix_antimeridian=None):
    """
    Convert ISEA3H cell IDs to Shapely geometry objects.

    Accepts a single isea3h_id (string) or a list of isea3h_ids. For each valid ISEA3H cell ID,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    isea3h_ids : str or list of str
        ISEA3H cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "1327916769,-55086"
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none

    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single ISEA3H cell ID is provided, returns a single Shapely Polygon object.
        If a list of IDs is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding ISEA3H cell.

    Examples
    --------
    >>> isea3h2geo("1327916769,-55086")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> isea3h2geo(["1327916769,-55086", "1327916770,-55087"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(isea3h_ids, str):
        isea3h_ids = [isea3h_ids]
    isea3h_polygons = []
    for isea3h_id in isea3h_ids:
        try:
            isea3h_cell = DggsCell(isea3h_id)
            cell_polygon = isea3h_cell_to_polygon(isea3h_cell)
            if fix_antimeridian == "shift" or fix_antimeridian == "shift_balanced":
                cell_polygon = shift_balanced(cell_polygon)
            elif fix_antimeridian == "shift_west":
                cell_polygon = shift_west(cell_polygon)
            elif fix_antimeridian == "shift_east":
                cell_polygon = shift_east(cell_polygon)
            elif fix_antimeridian == "split":
                cell_polygon = fix_polygon(cell_polygon)

            isea3h_polygons.append(cell_polygon)
        except Exception:
            continue

    if len(isea3h_polygons) == 1:
        return isea3h_polygons[0]
    return isea3h_polygons


def isea3h2geo_cli():
    """
    Command-line interface for isea3h2geo supporting multiple ISEA3H cell IDs.
    """
    parser = argparse.ArgumentParser(
        description="Convert ISEA3H cell ID(s) to Shapely Polygons"
    )
    parser.add_argument(
        "isea3h",
        nargs="+",
        help="Input ISEA3H cell ID(s), e.g., isea3h2geo 1327916769,-55086 ...",
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=[
            "shift",
            "shift_balanced",
            "shift_west",
            "shift_east",
            "split",
            "none",
        ],
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    args = parser.parse_args()
    if platform.system() == "Windows":
        polys = isea3h2geo(args.isea3h, fix_antimeridian=args.fix_antimeridian)
        return polys
    else:
        print("ISEA3H is only supported on Windows systems")


def isea3h2geojson(isea3h_ids, fix_antimeridian=None):
    """
    Convert ISEA3H cell IDs to GeoJSON FeatureCollection.

    Accepts a single isea3h_id (string) or a list of isea3h_ids. For each valid ISEA3H cell ID,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    isea3h_ids : str or list of str
        ISEA3H cell ID(s) to convert. Can be a single string or a list of strings.
        Example format: "1327916769,-55086"
    fix_antimeridian : str, optional
        Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none

    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid ISEA3H cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the ISEA3H cell ID, resolution level, center coordinates, edge length, and cell area

    Examples
    --------
    >>> isea3h2geojson("1327916769,-55086")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> isea3h2geojson(["1327916769,-55086", "1327916770,-55087"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(isea3h_ids, str):
        isea3h_ids = [isea3h_ids]
    features = []
    for isea3h_id in isea3h_ids:
        try:
            isea3h_cell = DggsCell(isea3h_id)
            cell_polygon = isea3h2geo(isea3h_id)
            if fix_antimeridian:
                cell_polygon = fix_polygon(cell_polygon)
            cell_centroid = cell_polygon.centroid
            center_lat = cell_centroid.y
            center_lon = cell_centroid.x
            cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
            cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
            isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)
            cell_accuracy = isea3h2point._accuracy
            avg_edge_len = cell_perimeter / 6
            cell_resolution = ISEA3H_ACCURACY_RES_DICT.get(cell_accuracy)
            if cell_resolution == 0:
                avg_edge_len = cell_perimeter / 3
            if cell_accuracy == 0.0:
                if round(avg_edge_len, 2) == 0.06:
                    cell_resolution = 33
                elif round(avg_edge_len, 2) == 0.03:
                    cell_resolution = 34
                elif round(avg_edge_len, 2) == 0.02:
                    cell_resolution = 35
                elif round(avg_edge_len, 2) == 0.01:
                    cell_resolution = 36
                elif round(avg_edge_len, 3) == 0.007:
                    cell_resolution = 37
                elif round(avg_edge_len, 3) == 0.004:
                    cell_resolution = 38
                elif round(avg_edge_len, 3) == 0.002:
                    cell_resolution = 39
                elif round(avg_edge_len, 3) <= 0.001:
                    cell_resolution = 40
            feature = {
                "type": "Feature",
                "geometry": mapping(cell_polygon),
                "properties": {
                    "isea3h": isea3h_id,
                    "resolution": cell_resolution,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "avg_edge_len": round(avg_edge_len, 3),
                    "cell_area": cell_area,
                },
            }
            features.append(feature)
        except Exception:
            continue
    feature_collection = {"type": "FeatureCollection", "features": features}
    return feature_collection


def isea3h2geojson_cli():
    """
    Command-line interface for isea3h2geojson supporting multiple ISEA3H cell IDs.
    """
    parser = argparse.ArgumentParser(description="Convert ISEA3H ID(s) to GeoJSON")
    parser.add_argument(
        "isea3h",
        nargs="+",
        help="Input ISEA3H cell ID(s), e.g., isea3h2geojson 1327916769,-55086 ...",
    )
    parser.add_argument(
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=[
            "shift",
            "shift_balanced",
            "shift_west",
            "shift_east",
            "split",
            "none",
        ],
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    args = parser.parse_args()
    if platform.system() == "Windows":
        geojson_data = json.dumps(
            isea3h2geojson(args.isea3h, fix_antimeridian=args.fix_antimeridian)
        )
        print(geojson_data)
    else:
        print("ISEA3H is only supported on Windows systems")
