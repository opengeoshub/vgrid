"""
Vector to S2 Module

This module provides functionality to convert vector geometries to S2 grid cells with flexible input and output formats.

Key Functions:
    point2s2: Convert point geometries to S2 cells
    polyline2s2: Convert line geometries to S2 cells
    polygon2s2: Convert polygon geometries to S2 cells with spatial predicates
    geodataframe2s2: Convert GeoDataFrame to S2 cells with topology support
    vector2s2: Main function for converting various input formats to S2 cells
    vector2s2_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
from collections import deque
from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import LineString, MultiPoint, Point
from vgrid.dggs import s2
from vgrid.utils.geometry import (
    geodesic_dggs_to_geoseries,
    check_predicate,
    shortest_point_distance,
    strip_duplicate_and_collinear_vertices,
)
from vgrid.conversion.dggs2geo.s22geo import s22geo
from vgrid.conversion.latlon2dggs import latlon2s2
from vgrid.conversion.dggscompact.s2compact import s2_compact
from vgrid.stats.s2stats import s2_metrics
from vgrid.utils.io import (
    validate_s2_resolution,
    process_input_data_vector,
    convert_to_output_format,
    add_verbose_argument,
    add_compact_depth_argument,
)
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.io import DGGS_TYPES

min_res = DGGS_TYPES["s2"]["min_res"]
max_res = DGGS_TYPES["s2"]["max_res"]


def point2s2(
    feature,
    resolution,
    feature_properties=None,
    include_properties=True,
    fix_antimeridian=None,
):
    """
    Convert a point geometry to S2 grid cells.

    Converts point or multipoint geometries to S2 grid cells at the specified resolution.
    Each point is assigned to its containing S2 cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to S2 cells.
    resolution : int
        S2 resolution level [0..30].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable S2 compact mode (not used for points).
    include_properties : bool, optional
        Whether to include properties in output.

    Returns
    -------
    list of dict
        List of dictionaries representing S2 cells containing the point(s).
        Each dictionary contains S2 cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2s2(point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2s2(points, 8)
    >>> len(cells)
    2
    """
    s2_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        lat_lng = s2.LatLng.from_degrees(point.y, point.x)
        cell_id_max_res = s2.CellId.from_lat_lng(lat_lng)
        cell_id = cell_id_max_res.parent(resolution)
        s2_cell = s2.Cell(cell_id)
        cell_token = s2.CellId.to_token(s2_cell.id())
        if s2_cell:
            cell_polygon = s22geo(cell_token, fix_antimeridian=fix_antimeridian)
            cell_resolution = cell_id.level()
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "s2", cell_token, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            s2_rows.append(row)
    return s2_rows


def _s2_segment_cells(
    resolution,
    start_xy,
    end_xy,
    fix_antimeridian=None,
):
    """Walk S2 cells along one polyline segment via intersecting edge neighbors."""
    start_x, start_y = start_xy
    end_x, end_y = end_xy
    segment = LineString([(start_x, start_y), (end_x, end_y)])
    start_id = latlon2s2(start_y, start_x, resolution)
    if segment.is_empty or segment.length == 0:
        return [start_id]

    end_id = latlon2s2(end_y, end_x, resolution)
    end_pt = Point(end_x, end_y)

    start_poly = s22geo(start_id, fix_antimeridian=fix_antimeridian)
    if start_poly is None:
        return [end_id] if start_id != end_id else [start_id]

    if start_id == end_id or start_poly.covers(end_pt):
        ordered = [start_id]
        if end_id != start_id:
            ordered.append(end_id)
        return ordered

    ordered = []
    visited = set()
    queue = deque([start_id])
    reached_endpoint = False

    while queue:
        cell_id = queue.popleft()
        if cell_id in visited:
            continue
        visited.add(cell_id)

        cell_polygon = s22geo(cell_id, fix_antimeridian=fix_antimeridian)
        if cell_polygon is None or cell_polygon.is_empty:
            continue
        if cell_id != start_id and not cell_polygon.intersects(segment):
            continue

        if not ordered or ordered[-1] != cell_id:
            ordered.append(cell_id)

        if cell_id == end_id or cell_polygon.covers(end_pt):
            reached_endpoint = True
            continue

        cell = s2.CellId.from_token(cell_id)
        for neighbor in cell.get_edge_neighbors():
            neighbor_id = neighbor.to_token()
            if neighbor_id in visited:
                continue
            neighbor_polygon = s22geo(
                neighbor_id, fix_antimeridian=fix_antimeridian
            )
            if neighbor_polygon is None or neighbor_polygon.is_empty:
                continue
            if neighbor_polygon.intersects(segment):
                queue.append(neighbor_id)

    if not reached_endpoint and end_id not in visited:
        ordered.append(end_id)
    return ordered


def polyline2s2(
    feature,
    resolution,
    feature_properties=None,
    include_properties=True,
    all_polylines=None,
    fix_antimeridian=None,
):
    """
    Convert each polyline to S2 cells by walking intersecting neighbors.

    For each LineString or MultiLineString:
    - Remove consecutive duplicate vertices and collinear intermediate vertices
    - Split the cleaned polyline into segments
    - Start from the cell containing each segment start
    - If that cell already contains the endpoint, continue to the next segment
    - Otherwise expand edge neighbors, keep those that intersect the segment,
      and repeat until the endpoint cell is reached
    """
    s2_rows = []
    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        coords = strip_duplicate_and_collinear_vertices(polyline)
        if len(coords) < 2:
            continue

        ordered_cells = []
        for i in range(len(coords) - 1):
            segment_cells = _s2_segment_cells(
                resolution,
                coords[i],
                coords[i + 1],
                fix_antimeridian=fix_antimeridian,
            )
            for cell_id in segment_cells:
                if ordered_cells and ordered_cells[-1] == cell_id:
                    continue
                ordered_cells.append(cell_id)

        for cell_id in ordered_cells:
            cell_polygon = s22geo(cell_id, fix_antimeridian=fix_antimeridian)
            if cell_polygon is None:
                continue
            cell = s2.CellId.from_token(cell_id)
            cell_resolution = cell.level()
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "s2", cell_id, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            s2_rows.append(row)
    return s2_rows


def polygon2s2(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    depth=-1,
    include_properties=True,
    all_polygons=None,
    fix_antimeridian=None,
    verbose=True,
):
    """
    Convert a polygon geometry to S2 grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): S2 resolution level [0..30]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable S2 compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2s2)
        include_properties (bool, optional): Whether to include properties in output
        all_polygons (list, optional): List of all polygons for topology preservation (not used in this function)

    Returns:
        list: List of dictionaries representing S2 cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])
        >>> cells = polygon2s2(poly, 10, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """
    s2_rows = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        polygon_rows = []
        min_lng, min_lat, max_lng, max_lat = polygon.bounds
        level = resolution
        coverer = s2.RegionCoverer()
        coverer.min_level = level
        coverer.max_level = level
        region = s2.LatLngRect(
            s2.LatLng.from_degrees(min_lat, min_lng),
            s2.LatLng.from_degrees(max_lat, max_lng),
        )
        covering = coverer.get_covering(region)
        cell_ids = covering
        for cell_id in cell_ids:
            cell_polygon = s22geo(cell_id.to_token(), fix_antimeridian=fix_antimeridian)
            if not check_predicate(cell_polygon, polygon, predicate):
                continue
            cell_token = s2.CellId.to_token(cell_id)
            cell_resolution = cell_id.level()
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "s2", cell_token, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            polygon_rows.append(row)

        if compact and polygon_rows:
            tokens = [row.get("s2") for row in polygon_rows if row.get("s2")]
            compact_tokens = s2_compact(tokens, depth=depth, verbose=verbose)
            compact_rows = []
            for cell_token in compact_tokens:
                cell_polygon = s22geo(
                    cell_token, fix_antimeridian=fix_antimeridian
                )
                compact_cell = s2.CellId.from_token(cell_token)
                cell_resolution = compact_cell.level()
                num_edges = 4
                row = geodesic_dggs_to_geoseries(
                    "s2", cell_token, cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                compact_rows.append(row)
            polygon_rows = compact_rows

        s2_rows.extend(polygon_rows)

    return s2_rows


def geodataframe2s2(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    depth=-1,
    topology=False,
    include_properties=True,
    fix_antimeridian=None,
    verbose=True,
):
    """
    Convert a GeoDataFrame to S2 grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): S2 resolution level [0..30]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable S2 compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint S2 cells
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with S2 grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2s2(gdf, 10)
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

            # Find resolution where S2 cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint S2 cells
            if shortest_distance > 0:
                for res in range(min_res, max_res + 1):
                    _, avg_edge_length, _, _ = s2_metrics(res)
                    cell_diameter = (
                        avg_edge_length * 2
                    )  # math.sqrt(2) # to bypass the big aperture
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_s2_resolution(resolution)

    s2_rows = []

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
            s2_rows.extend(
                point2s2(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )

        elif geom.geom_type in ("LineString", "MultiLineString"):
            s2_rows.extend(
                polyline2s2(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            s2_rows.extend(
                polygon2s2(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    depth=depth,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                    verbose=verbose,
                )
            )
    return gpd.GeoDataFrame(s2_rows, geometry="geometry", crs="EPSG:4326")


# --- Main vector2s2 function ---
def vector2s2(
    vector_data,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    output_format='gpd',
    include_properties=True,
    fix_antimeridian=None,
    verbose=True,
    depth=-1,
    **kwargs,
):
    """
    Convert vector data to S2 grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): S2 resolution level [0..30]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable S2 compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint S2 cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2s2("data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_s2_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2s2(
        gdf,
        resolution,
        predicate,
        compact,
        depth,
        topology,
        include_properties,
        fix_antimeridian,
        verbose=verbose,
    )

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2s2"
        else:
            output_name = "s2"
    return convert_to_output_format(result, output_format, output_name)


# --- CLI ---
def vector2s2_cli():
    """
    Command-line interface for vector2s2 conversion.

    This function provides a command-line interface for converting vector data to S2 grid cells.
    It supports various input formats and output formats, with options for resolution control,
    spatial predicates, compact mode, and topology preservation.

    Usage:
        python vector2s2.py -i input.geojson -r 10 -f geojson
        python vector2s2.py -i input.shp -r 8 -p intersect -c -t
    """
    parser = argparse.ArgumentParser(description="Convert vector data to S2 grid cells")
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"S2 resolution [{min_res}..{max_res}]. Required when topology=False, auto-calculated when topology=True",
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
        help="Enable S2 compact mode for polygons",
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
    add_compact_depth_argument(parser)
    add_verbose_argument(parser)
    args = parser.parse_args()

    try:
        result = vector2s2(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            depth=args.depth,
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


if __name__ == "__main__":
    vector2s2_cli()
