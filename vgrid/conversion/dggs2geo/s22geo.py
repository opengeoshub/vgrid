"""
S2 to Geometry Module

This module provides functionality to convert S2 cell tokens to Shapely Polygons and GeoJSON FeatureCollection.

Key Functions:
    s22geo: Convert S2 tokens to Shapely Polygons
    s22geojson: Convert S2 tokens to GeoJSON FeatureCollection
    s22_cli: Command-line interface for polygon conversion
    s22geojson_cli: Command-line interface for GeoJSON conversion
"""

import json
import argparse
from shapely.geometry import Polygon
from vgrid.utils.geometry import geodesic_dggs_to_feature
from vgrid.dggs import s2
from vgrid.utils.antimeridian import fix_polygon
from vgrid.utils.geometry import shift_balanced, shift_west, shift_east


def s22geo(s2_tokens, fix_antimeridian=None):
    """
    Convert S2 cell tokens to Shapely geometry objects.

    Accepts a single s2_token (string) or a list of s2_tokens. For each valid S2 cell token,
    creates a Shapely Polygon representing the grid cell boundaries. Skips invalid or
    error-prone cells.

    Parameters
    ----------
    s2_tokens : str or list of str
        S2 cell token(s) to convert. Can be a single string or a list of strings.
        Example format: "31752f45cc94"
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
    Returns
    -------
    shapely.geometry.Polygon or list of shapely.geometry.Polygon
        If a single S2 cell token is provided, returns a single Shapely Polygon object.
        If a list of tokens is provided, returns a list of Shapely Polygon objects.
        Each polygon represents the boundaries of the corresponding S2 cell.

    Examples
    --------
    >>> s22geo("31752f45cc94")
    <shapely.geometry.polygon.Polygon object at ...>

    >>> s22geo(["31752f45cc94", "31752f45cc95"])
    [<shapely.geometry.polygon.Polygon object at ...>, <shapely.geometry.polygon.Polygon object at ...>]
    """
    if isinstance(s2_tokens, str):
        s2_tokens = [s2_tokens]
    s2_polygons = []
    for s2_token in s2_tokens:
        try:
            cell_id = s2.CellId.from_token(s2_token)
            cell = s2.Cell(cell_id)
            vertices = [cell.get_vertex(i) for i in range(4)]
            shapely_vertices = []
            for vertex in vertices:
                lat_lng = s2.LatLng.from_point(vertex)
                longitude = lat_lng.lng().degrees
                latitude = lat_lng.lat().degrees
                shapely_vertices.append((longitude, latitude))
            shapely_vertices.append(shapely_vertices[0])
            cell_polygon = Polygon(shapely_vertices)
            if fix_antimeridian == "shift" or fix_antimeridian == "shift_balanced":
                cell_polygon = shift_balanced(
                    cell_polygon, threshold_west=-90, threshold_east=90
                )
            elif fix_antimeridian == "shift_west":
                cell_polygon = shift_west(cell_polygon, threshold=-90)
            elif fix_antimeridian == "shift_east":
                cell_polygon = shift_east(cell_polygon, threshold=90)  # 130?
            elif fix_antimeridian == "split":
                cell_polygon = fix_polygon(cell_polygon)
                # if resolution == 0:  #
                #     cell_polygon = fix_polygon(cell_polygon)
                # elif (
                #         s2_token.startswith("00")
                #         or s2_token.startswith("09")
                #         or s2_token.startswith("14")
                #         or s2_token.startswith("04")
                #         or s2_token.startswith("19")
                #     ):
                #     cell_polygon = fix_polygon(cell_polygon)
            s2_polygons.append(cell_polygon)
        except Exception:
            continue
    if len(s2_polygons) == 1:
        return s2_polygons[0]
    return s2_polygons


def s22geo_cli():
    """
    Command-line interface for s22geo supporting multiple S2 cell tokens.
    """
    parser = argparse.ArgumentParser(
        description="Convert S2 cell token(s) to Shapely Polygons"
    )
    parser.add_argument(
        "s2",
        nargs="+",
        help="Input S2 cell token(s), e.g., s22geo 31752f45cc94 31752f45cc95",
    )
    args = parser.parse_args()
    polys = s22geo(args.s2)
    return polys


def s22geojson(s2_tokens, fix_antimeridian=None):
    """
    Convert S2 cell tokens to GeoJSON FeatureCollection.

    Accepts a single s2_token (string) or a list of s2_tokens. For each valid S2 cell token,
    creates a GeoJSON feature with polygon geometry representing the grid cell boundaries.
    Skips invalid or error-prone cells.

    Parameters
    ----------
    s2_tokens : str or list of str
        S2 cell token(s) to convert. Can be a single string or a list of strings.
        Example format: "31752f45cc94"
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
    Returns
    -------
    dict
        A GeoJSON FeatureCollection containing polygon features for each valid S2 cell.
        Each feature includes:
        - geometry: Polygon representing the cell boundaries
        - properties: Contains the S2 cell token, resolution level, and cell metadata

    Examples
    --------
    >>> s22geojson("31752f45cc94")
    {'type': 'FeatureCollection', 'features': [...]}

    >>> s22geojson(["31752f45cc94", "31752f45cc95"])
    {'type': 'FeatureCollection', 'features': [...]}
    """
    if isinstance(s2_tokens, str):
        s2_tokens = [s2_tokens]
    s2_features = []
    for s2_token in s2_tokens:
        try:
            cell_id = s2.CellId.from_token(s2_token)
            cell_polygon = s22geo(s2_token, fix_antimeridian=fix_antimeridian)
            resolution = cell_id.level()
            num_edges = 4
            s2_feature = geodesic_dggs_to_feature(
                "s2", s2_token, resolution, cell_polygon, num_edges
            )
            s2_features.append(s2_feature)
        except Exception:
            continue
    return {"type": "FeatureCollection", "features": s2_features}


def s22geojson_cli():
    """
    Command-line interface for s22geojson supporting multiple S2 cell tokens.
    """
    parser = argparse.ArgumentParser(description="Convert S2 cell token(s) to GeoJSON")
    parser.add_argument(
        "s2",
        nargs="+",
        help="Input S2 cell token(s), e.g., s22geojson 31752f45cc94 31752f45cc95",
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
    fix_antimeridian = args.fix_antimeridian
    geojson_data = json.dumps(s22geojson(args.s2, fix_antimeridian=fix_antimeridian))
    print(geojson_data)
