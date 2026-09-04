"""
Vector to QTM Module

This module provides functionality to convert vector geometries to QTM grid cells with flexible input and output formats.

Key Functions:
    point2qtm: Convert point geometries to QTM cells
    polyline2qtm: Convert line geometries to QTM cells
    polygon2qtm: Convert polygon geometries to QTM cells with spatial predicates
    geodataframe2qtm: Convert GeoDataFrame to QTM cells with topology support
    vector2qtm: Main function for converting various input formats to QTM cells
    vector2qtm_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
import math
from tqdm import tqdm
from shapely.geometry import Polygon, MultiPoint
import geopandas as gpd
from vgrid.dggs import qtm
from vgrid.dggs.qtm import QTM_INITIAL_FACETS
from vgrid.conversion.dggscompact.qtmcompact import qtmcompact
from vgrid.conversion.dggs2geo.qtm2geo import qtm2geo
from vgrid.utils.geometry import (
    check_predicate,
    shortest_point_distance,
    geodesic_dggs_to_geoseries,
)
from vgrid.stats.qtmstats import qtm_metrics
from vgrid.utils.io import (
    validate_qtm_resolution,
    process_input_data_vector,
    convert_to_output_format,
    DGGS_TYPES,
    add_verbose_argument,
    add_compact_depth_argument,
)
from vgrid.utils.constants import STRUCTURED_FORMATS, OUTPUT_FORMATS

min_res = DGGS_TYPES["qtm"]["min_res"]
max_res = DGGS_TYPES["qtm"]["max_res"]


def point2qtm(
    feature,
    resolution,
    feature_properties=None,
    include_properties=True,
):
    """
    Convert a point geometry to QTM grid cells.

    Converts point or multipoint geometries to QTM grid cells at the specified resolution.
    Each point is assigned to its containing QTM cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to QTM cells.
    resolution : int
        QTM resolution level [1..24].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable QTM compact mode (not used for points).
    include_properties : bool, optional
        Whether to include properties in output.

    Returns
    -------
    list of dict
        List of dictionaries representing QTM cells containing the point(s).
        Each dictionary contains QTM cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2qtm(point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2qtm(points, 8)
    >>> len(cells)
    2

    Returns:
        list: List of GeoJSON feature dictionaries representing QTM cells containing the point
    """
    qtm_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        latitude = point.y
        longitude = point.x
        qtm_id = qtm.latlon_to_qtm_id(latitude, longitude, resolution)
        cell_polygon = qtm2geo(qtm_id)
        if cell_polygon:
            num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "qtm", qtm_id, resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            qtm_rows.append(row)
    return qtm_rows


def polyline2qtm(
    feature,
    resolution,
    feature_properties=None,
    include_properties=True,
):
    """
    Convert line geometries (LineString, MultiLineString) to QTM grid cells.

    Args:
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Line geometry to convert
        resolution (int): QTM resolution [1..24]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for lines)
        compact (bool, optional): Enable QTM compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2qtm)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing QTM cells intersecting the line
    """
    qtm_rows = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []
    for polyline in polylines:
        levelFacets = {}
        QTMID = {}
        for lvl in range(resolution):
            levelFacets[lvl] = []
            QTMID[lvl] = []
            if lvl == 0:
                for i, facet in enumerate(QTM_INITIAL_FACETS):
                    QTMID[0].append(str(i + 1))
                    levelFacets[0].append(facet)
                    facet_geom = qtm.constructGeometry(facet)
                    if Polygon(facet_geom).intersects(polyline) and resolution == 1:
                        qtm_id = QTMID[0][i]
                        num_edges = 3
                        row = geodesic_dggs_to_geoseries(
                            "qtm", qtm_id, resolution, facet_geom, num_edges
                        )
                        if include_properties and feature_properties:
                            row.update(feature_properties)
                        qtm_rows.append(row)
                        return qtm_rows
            else:
                for i, pf in enumerate(levelFacets[lvl - 1]):
                    subdivided_facets = qtm.divideFacet(pf)
                    for j, subfacet in enumerate(subdivided_facets):
                        subfacet_geom = qtm.constructGeometry(subfacet)
                        if Polygon(subfacet_geom).intersects(polyline):
                            new_id = QTMID[lvl - 1][i] + str(j)
                            QTMID[lvl].append(new_id)
                            levelFacets[lvl].append(subfacet)
                            if lvl == resolution - 1:
                                num_edges = 3
                                row = geodesic_dggs_to_geoseries(
                                    "qtm", new_id, resolution, subfacet_geom, num_edges
                                )
                                if include_properties and feature_properties:
                                    row.update(feature_properties)
                                qtm_rows.append(row)
    return qtm_rows


def polygon2qtm(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    depth=-1,
    include_properties=True,
    verbose=True,
):
    """
    Convert polygon geometries (Polygon, MultiPolygon) to QTM grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): QTM resolution [1..24]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable QTM compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2qtm)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing QTM cells based on predicate
    """
    qtm_rows = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        levelFacets = {}
        QTMID = {}
        for lvl in range(resolution):
            levelFacets[lvl] = []
            QTMID[lvl] = []
            if lvl == 0:
                for i, facet in enumerate(QTM_INITIAL_FACETS):
                    QTMID[0].append(str(i + 1))
                    levelFacets[0].append(facet)
                    facet_geom = qtm.constructGeometry(facet)
                    if Polygon(facet_geom).intersects(polygon) and resolution == 1:
                        qtm_id = QTMID[0][i]
                        num_edges = 3
                        row = geodesic_dggs_to_geoseries(
                            "qtm", qtm_id, resolution, facet_geom, num_edges
                        )
                        if include_properties and feature_properties:
                            row.update(feature_properties)
                        qtm_rows.append(row)
                        return qtm_rows
            else:
                for i, pf in enumerate(levelFacets[lvl - 1]):
                    subdivided_facets = qtm.divideFacet(pf)
                    for j, subfacet in enumerate(subdivided_facets):
                        subfacet_geom = qtm.constructGeometry(subfacet)
                        if Polygon(subfacet_geom).intersects(polygon):
                            new_id = QTMID[lvl - 1][i] + str(j)
                            QTMID[lvl].append(new_id)
                            levelFacets[lvl].append(subfacet)
                            if lvl == resolution - 1:
                                if not check_predicate(
                                    Polygon(subfacet_geom), polygon, predicate
                                ):
                                    continue
                                num_edges = 3
                                row = geodesic_dggs_to_geoseries(
                                    "qtm", new_id, resolution, subfacet_geom, num_edges
                                )
                                if include_properties and feature_properties:
                                    row.update(feature_properties)
                                qtm_rows.append(row)
    if compact and qtm_rows:
        temp_gdf = gpd.GeoDataFrame(qtm_rows, geometry="geometry", crs="EPSG:4326")
        compacted_gdf = qtmcompact(
            temp_gdf, qtm_id="qtm", output_format="gpd", verbose=verbose, depth=depth
        )
        if compacted_gdf is not None:
            qtm_rows = compacted_gdf.to_dict("records")
    return qtm_rows


def geodataframe2qtm(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    depth=-1,
    topology=False,
    include_properties=True,
    verbose=True,
):
    """
    Convert a GeoDataFrame to QTM grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): QTM resolution [1..24]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable QTM compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint QTM cells
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with QTM grid cells
    """

    # Process topology for points and multipoints if enabled
    if topology:
        resolution = max_res
        # Collect all points for topology preservation
        points_list = []
        for _, row in gdf.iterrows():
            geom = row.geometry
            if geom is None:
                continue
            if geom.geom_type == "Point":
                points_list.append(geom)
            elif geom.geom_type == "MultiPoint":
                points_list.extend(list(geom.geoms))

        if points_list:
            all_points = MultiPoint(points_list)

            # Calculate the shortest distance between all points
            shortest_distance = shortest_point_distance(all_points)

            # Find resolution where QTM cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint QTM cells
            if shortest_distance > 0:
                for res in range(
                    min_res, max_res + 1
                ):  # QTM resolution range is [1..24]
                    _, avg_edge_length, _, _ = qtm_metrics(res)
                    # Use a factor to ensure sufficient separation (triangle diameter is ~2x edge length)
                    triangle_diameter = (2 * avg_edge_length * math.sqrt(3)) / 3
                    if triangle_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_qtm_resolution(resolution)

    qtm_rows = []

    for _, row in tqdm(gdf.iterrows(), desc="Processing features", total=len(gdf), disable=not verbose):
        geom = row.geometry
        if geom is None:
            continue

        props = row.to_dict()
        if "geometry" in props:
            del props["geometry"]

        if not include_properties:
            props = {}

        if geom.geom_type == "Point" or geom.geom_type == "MultiPoint":
            qtm_rows.extend(
                point2qtm(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    include_properties=include_properties,
                )
            )

        elif geom.geom_type in ("LineString", "MultiLineString"):
            qtm_rows.extend(
                polyline2qtm(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            qtm_rows.extend(
                polygon2qtm(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    depth=depth,
                    include_properties=include_properties,
                    verbose=verbose,
                )
            )
    return gpd.GeoDataFrame(qtm_rows, geometry="geometry", crs="EPSG:4326")


# --- Main vector2qtm function ---
def vector2qtm(
    vector_data,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    output_format='gpd',
    include_properties=True,
    verbose=True,
    depth=-1,
    **kwargs,
):
    """
    Convert vector data to QTM grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): QTM resolution [1..24]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable QTM compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint QTM cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_qtm_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2qtm(
        gdf, resolution, predicate, compact, depth, topology, include_properties,
        verbose=verbose,
    )

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2qtm"
        else:
            output_name = "qtm"
    return convert_to_output_format(result, output_format, output_name)


def vector2qtm_cli():
    """
    Command-line interface for vector2qtm conversion.
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to QTM grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        metavar=f"[{min_res}..{max_res}]",
        help=f"QTM resolution [{min_res}..{max_res}] (coarsest={min_res}, finest={max_res})",
    )
    parser.add_argument(
        "-p",
        "--predicate",
        choices=["intersect", "within", "centroid_within", "largest_overlap"],
        help="Spatial predicate: intersect, within, centroid_within, largest_overlap for polygons",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable QTM compact mode for polygons",
    )
    parser.add_argument(
        "-t", "--topology", action="store_true", help="Enable topology preserving mode"
    )
    parser.add_argument(
        "-np",
        "-no-props",
        dest="include_properties",
        action="store_false",
        help="Do not include original feature properties.",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format (default: gpd).",
    )

    add_compact_depth_argument(parser)
    add_verbose_argument(parser)
    args = parser.parse_args()
    args.resolution = validate_qtm_resolution(args.resolution)
    try:
        result = vector2qtm(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            depth=args.depth,
            topology=args.topology,
            output_format=args.output_format,
            include_properties=args.include_properties,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
        # For file outputs, the utility prints the saved path
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2qtm_cli()
