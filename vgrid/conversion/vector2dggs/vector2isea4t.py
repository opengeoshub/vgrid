"""
Vector to ISEA4T Module

This module provides functionality to convert vector geometries to ISEA4T grid cells with flexible input and output formats.

Key Functions:
    point2isea4t: Convert point geometries to ISEA4T cells
    polyline2isea4t: Convert line geometries to ISEA4T cells
    polygon2isea4t: Convert polygon geometries to ISEA4T cells with spatial predicates
    geodataframe2isea4t: Convert GeoDataFrame to ISEA4T cells with topology support
    vector2isea4t: Main function for converting various input formats to ISEA4T cells
    vector2isea4t_cli: Command-line interface for batch processing

Note: This module is only supported on Windows systems due to OpenEaggr dependency.
"""

import sys
import os
import math
import argparse
import geopandas as gpd
from tqdm import tqdm
from shapely.geometry import box, MultiPoint
from vgrid.utils.geometry import check_predicate, shortest_point_distance
from vgrid.utils.io import (
    validate_isea4t_resolution,
    process_input_data_vector,
    convert_to_output_format,
    DGGS_TYPES,
    add_verbose_argument,
)
import platform
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS

min_res = DGGS_TYPES["isea4t"]["min_res"]
max_res = DGGS_TYPES["isea4t"]["max_res"]

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.enums.model import Model
    from vgrid.dggs.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.dggs.eaggr.shapes.lat_long_point import LatLongPoint
    from vgrid.generator.isea4tgrid import get_isea4t_children_cells_within_bbox
    from vgrid.utils.constants import ISEA4T_RES_ACCURACY_DICT
    from vgrid.utils.geometry import geodesic_dggs_to_geoseries
    from vgrid.conversion.dggscompact.isea4tcompact import isea4t_compact
    from vgrid.conversion.dggs2geo.isea4t2geo import isea4t2geo
    from vgrid.stats.isea4tstats import isea4t_metrics

    isea4t_dggs = Eaggr(Model.ISEA4T)


def point2isea4t(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    fix_antimeridian=None,
):
    """
    Convert a point geometry to ISEA4T grid cells.

    Converts point or multipoint geometries to ISEA4T grid cells at the specified resolution.
    Each point is assigned to its containing ISEA4T cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to ISEA4T cells.
    resolution : int
        ISEA4T resolution level [0..30].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable ISEA4T compact mode (not used for points).
    topology : bool, optional
        Enable topology preserving mode (handled by geodataframe2isea4t).
    include_properties : bool, optional
        Whether to include properties in output.
    fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none

    Returns
    -------
    list of dict
        List of dictionaries representing ISEA4T cells containing the point(s).
        Each dictionary contains ISEA4T cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2isea4t(point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2isea4t(points, 8)
    >>> len(cells)
    2
    """
    isea4t_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        accuracy = ISEA4T_RES_ACCURACY_DICT.get(resolution)
        lat_long_point = LatLongPoint(point.y, point.x, accuracy)
        isea4t_cell = isea4t_dggs.convert_point_to_dggs_cell(lat_long_point)
        isea4t_id = isea4t_cell.get_cell_id()
        cell_polygon = isea4t2geo(isea4t_id, fix_antimeridian=fix_antimeridian)
        num_edges = 3
        row = geodesic_dggs_to_geoseries(
            "isea4t", isea4t_id, resolution, cell_polygon, num_edges
        )
        if include_properties and feature_properties:
            row.update(feature_properties)
        isea4t_rows.append(row)
    return isea4t_rows


def polyline2isea4t(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    fix_antimeridian=None,
):
    """
    Convert a polyline geometry to ISEA4T grid cells.

    Args:
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Polyline geometry to convert
        resolution (int): ISEA4T resolution level [0..30]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for polylines)
        compact (bool, optional): Enable ISEA4T compact mode (not used for polylines)
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2isea4t)
        include_properties (bool, optional): Whether to include properties in output
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none

    Returns:
        list: List of dictionaries representing ISEA4T cells intersecting the polyline

    Example:
        >>> from shapely.geometry import LineString
        >>> line = LineString([(-122.4194, 37.7749), (-122.4000, 37.7800)])
        >>> cells = polyline2isea4t(line, 10, {"name": "route"})
        >>> len(cells) > 0
        True
    """
    isea4t_rows = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        accuracy = ISEA4T_RES_ACCURACY_DICT.get(resolution)
        bounding_box = box(*polyline.bounds)
        bounding_box_wkt = bounding_box.wkt
        shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(
            bounding_box_wkt, ShapeStringFormat.WKT, accuracy
        )
        shape = shapes[0]
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_child_cells = get_isea4t_children_cells_within_bbox(
            bounding_cell.get_cell_id(), bounding_box, resolution
        )
        for child in bounding_child_cells:
            isea4t_id = child
            cell_polygon = isea4t2geo(isea4t_id, fix_antimeridian=fix_antimeridian)
            if cell_polygon.intersects(polyline):
                num_edges = 3
                cell_resolution = len(isea4t_id) - 2
                row = geodesic_dggs_to_geoseries(
                    "isea4t", isea4t_id, cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                isea4t_rows.append(row)
    return isea4t_rows


def polygon2isea4t(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    fix_antimeridian=None,
):
    """
    Convert a polygon geometry to ISEA4T grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): ISEA4T resolution level [0..30]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable ISEA4T compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2isea4t)
        include_properties (bool, optional): Whether to include properties in output
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none

    Returns:
        list: List of dictionaries representing ISEA4T cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])
        >>> cells = polygon2isea4t(poly, 10, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """
    isea4t_rows = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        accuracy = ISEA4T_RES_ACCURACY_DICT.get(resolution)
        bounding_box = box(*polygon.bounds)
        bounding_box_wkt = bounding_box.wkt
        shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(
            bounding_box_wkt, ShapeStringFormat.WKT, accuracy
        )
        shape = shapes[0]
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_child_cells = get_isea4t_children_cells_within_bbox(
            bounding_cell.get_cell_id(), bounding_box, resolution
        )
        for child in bounding_child_cells:
            isea4t_id = child
            cell_polygon = isea4t2geo(isea4t_id, fix_antimeridian=fix_antimeridian)
            if check_predicate(cell_polygon, polygon, predicate):
                num_edges = 3
                cell_resolution = len(isea4t_id) - 2
                row = geodesic_dggs_to_geoseries(
                    "isea4t", isea4t_id, cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                isea4t_rows.append(row)

        # Compact mode: apply to isea4t_rows after predicate check
        if compact:
            # Extract cell IDs from isea4t_rows
            cells_to_process = [row.get("isea4t") for row in isea4t_rows]
            # Apply compact
            cells_to_process = isea4t_compact(cells_to_process)
            # Rebuild isea4t_rows with compacted cells
            isea4t_rows = []
            for cell_id in cells_to_process:
                cell_polygon = isea4t2geo(cell_id, fix_antimeridian=fix_antimeridian)
                num_edges = 3
                cell_resolution = len(cell_id) - 2
                row = geodesic_dggs_to_geoseries(
                    "isea4t", cell_id, cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                isea4t_rows.append(row)
    return isea4t_rows


def geodataframe2isea4t(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    fix_antimeridian=None,
    verbose=True,
):
    """
    Convert a GeoDataFrame to ISEA4T grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): ISEA4T resolution level [0..30]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable ISEA4T compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint ISEA4T cells
        include_properties (bool, optional): Whether to include properties in output
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with ISEA4T grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2isea4t(gdf, 10)
        >>> len(result) > 0
        True
    """

    # Process topology for points and multipoints if enabled
    if topology:
        estimated_resolution = max_res
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

            # Find resolution where ISEA4T cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint ISEA4T cells
            if shortest_distance > 0:
                for res in range(min_res, max_res + 1):
                    _, avg_edge_length, _, _ = isea4t_metrics(res)
                    cell_diameter = (2 * avg_edge_length * math.sqrt(3)) / 3
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_isea4t_resolution(resolution)

    isea4t_rows = []

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
            isea4t_rows.extend(
                point2isea4t(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    topology=topology,  # Topology already processed above
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )

        elif geom.geom_type in ("LineString", "MultiLineString"):
            isea4t_rows.extend(
                polyline2isea4t(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            isea4t_rows.extend(
                polygon2isea4t(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )

    if isea4t_rows:
        gdf = gpd.GeoDataFrame(isea4t_rows, geometry="geometry", crs="EPSG:4326")
    else:
        gdf = gpd.GeoDataFrame(columns=["geometry"], crs="EPSG:4326")
    return gdf


def vector2isea4t(
    vector_data,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    output_format="gpd",
    include_properties=True,
    fix_antimeridian=None,
    verbose=True,
    **kwargs,
):
    """
    Convert vector data to ISEA4T grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): ISEA4T resolution level [0..30]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable ISEA4T compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint ISEA4T cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2isea4t("data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_isea4t_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2isea4t(
        gdf,
        resolution,
        predicate,
        compact,
        topology,
        include_properties,
        fix_antimeridian,
        verbose=verbose,
    )

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2isea4t"
        else:
            output_name = "isea4t"
    return convert_to_output_format(result, output_format, output_name)


def vector2isea4t_cli():
    """
    Command-line interface for vector2isea4t conversion.

    This function provides a command-line interface for converting vector data to ISEA4T grid cells.
    It supports various input formats and output formats, with options for resolution control,
    spatial predicates, compact mode, and topology preservation.

    Usage:
        python vector2isea4t.py -i input.geojson -r 10 -f geojson
        python vector2isea4t.py -i input.shp -r 8 -p intersect -c -t
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to ISEA4T grid cells."
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"ISEA4T resolution [{min_res}..{max_res}]. Required when topology=False, auto-calculated when topology=True",
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
        help="Enable ISEA4T compact mode for polygons",
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
        choices=OUTPUT_FORMATS,
        help="Output format (default: gpd).",
        default="gpd",
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
    add_verbose_argument(parser)
    args = parser.parse_args()

    # Allow running on all platforms
    if platform.system() == "Windows":
        try:
            result = vector2isea4t(
                vector_data=args.input,
                resolution=args.resolution,
                predicate=args.predicate,
                compact=args.compact,
                topology=args.topology,
                output_format=args.output_format,
                include_properties=args.include_properties,
                fix_antimeridian=args.fix_antimeridian,
                verbose=args.verbose,
            )
            if args.output_format in STRUCTURED_FORMATS:
                print(result)
            # For file outputs, the utility prints the saved path
        except Exception as e:
            print(f"Error: {str(e)}", file=sys.stderr)
            sys.exit(1)
    else:
        print(
            "ISEA4T conversion is only supported on Windows systems.", file=sys.stderr
        )
        sys.exit(1)


if __name__ == "__main__":
    vector2isea4t_cli()
