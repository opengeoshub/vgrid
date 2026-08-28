"""
Vector to A5 Module

This module provides functionality to convert vector geometries to A5 grid cells with flexible input and output formats.

Key Functions:
    point2a5: Convert point geometries to A5 cells
    polyline2a5: Convert line geometries to A5 cells
    polygon2a5: Convert polygon geometries to A5 cells with spatial predicates
    geodataframe2a5: Convert GeoDataFrame to A5 cells with topology support
    vector2a5: Main function for converting various input formats to A5 cells
    vector2a5_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
import json
from collections import deque
from tqdm import tqdm
from pyproj import Geod
import geopandas as gpd
from shapely.geometry import MultiPoint
import a5
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.geometry import (
    check_predicate,
    shortest_point_distance,
)
from vgrid.utils.io import (
    validate_a5_resolution,
    process_input_data_vector,
    convert_to_output_format,
    DGGS_TYPES,
    add_verbose_argument,
)

from vgrid.conversion.latlon2dggs import latlon2a5
from vgrid.conversion.dggs2geo.a52geo import a52geo, a52geo_u64
from vgrid.conversion.dggscompact.a5compact import a5compact
from vgrid.stats.a5stats import a5_metrics
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS

geod = Geod(ellps="WGS84")
min_res = DGGS_TYPES["a5"]["min_res"]
max_res = DGGS_TYPES["a5"]["max_res"]


def point2a5(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    options=None,
    split_antimeridian=False,
):
    """
    Convert a point geometry to A5 grid cells.

    Converts point or multipoint geometries to A5 grid cells at the specified resolution.
    Each point is assigned to its containing A5 cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to A5 cells.
    resolution : int
        A5 resolution level [0..28].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable A5 compact mode (not used for points).
    topology : bool, optional
        Enable topology preserving mode (handled by geodataframe2a5).
    include_properties : bool, optional
        Whether to include properties in output.
    options : dict, optional
        Options for a52geo.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    Returns
    -------
    list of dict
        List of dictionaries representing A5 cells containing the point(s).
        Each dictionary contains A5 cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2a5(point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2a5(points, 8)
    >>> len(cells)
    2
    """
    a5_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []
    for point in points:
        a5_hex = latlon2a5(point.y, point.x, resolution)
        cell_polygon = a52geo(a5_hex, options, split_antimeridian=split_antimeridian)
        cell_resolution = a5.get_resolution(a5.hex_to_u64(a5_hex))
        num_edges = 5
        if cell_resolution == 1:
            num_edges = 3
        row = geodesic_dggs_to_geoseries(
            "a5", a5_hex, cell_resolution, cell_polygon, num_edges
        )

        # Add properties if requested
        if include_properties and feature_properties:
            row.update(feature_properties)

        a5_rows.append(row)
    return a5_rows


def polyline2a5(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    options=None,
    split_antimeridian=False,
):
    """
    Convert each polyline to an A5 path along its vertices.

    For each LineString or MultiLineString, this function:
    - Collects vertices as (lon, lat) waypoints
    - Uses a5.line_string_to_cells(waypoints, resolution) to trace cells
      along great-circle arcs between consecutive waypoints
    - Converts those cells back to geometries and appends them to a5_rows

    Parameters
    ----------
    feature : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Polyline geometry to convert.
    resolution : int
        A5 resolution level [0..28].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate (not used for polylines).
    compact : bool, optional
        Enable A5 compact mode (not used for polylines).
    topology : bool, optional
        Enable topology preserving mode (handled by geodataframe2a5).
    include_properties : bool, optional
        Whether to include properties in output.
    options : dict, optional
        Options for a52geo.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.

    Returns
    -------
    list of dict
        List of dictionaries representing A5 cells along the polyline(s).

    Examples
    --------
    >>> from shapely.geometry import LineString
    >>> line = LineString([(-122.4194, 37.7749), (-122.4000, 37.7800)])
    >>> cells = polyline2a5(line, 10, {"name": "route"})
    >>> len(cells) > 0
    True
    """
    a5_rows = []

    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        coords = list(polyline.coords)
        if len(coords) < 2:
            continue

        waypoints = [(lon, lat) for lon, lat in coords]

        try:
            ordered_cell_ids = a5.line_string_to_cells(waypoints, resolution)
        except Exception:
            ordered_cell_ids = []

        for cell_id in ordered_cell_ids:
            cell_hex = a5.u64_to_hex(cell_id)
            cell_polygon = a52geo(
                cell_hex, options, split_antimeridian=split_antimeridian
            )
            cell_resolution = a5.get_resolution(cell_id)
            num_edges = 5
            if cell_resolution == 1:
                num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "a5", cell_hex, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            a5_rows.append(row)

    return a5_rows


def polygon2a5(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    options=None,
    split_antimeridian=False,
    verbose=True,
):
    """
    Convert a polygon geometry to A5 grid cells.

    Discovery phase: pick a seed cell from the polygon's representative point,
    then BFS outward (grid disk radius 1) while cells intersect the polygon.
    Candidate cells are finally filtered by `predicate` against the input polygon.
    """
    a5_rows = []
    if feature.geom_type in ("Polygon",):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon",):
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        if polygon is None or polygon.is_empty:
            continue

        rep_pt = polygon.representative_point()
        seed_cell_id = a5.lonlat_to_cell((rep_pt.x, rep_pt.y), resolution)
        seed_cell_resolution = a5.get_resolution(seed_cell_id)
        seed_cell_polygon = a52geo_u64(
            seed_cell_id,
            options=options,
            split_antimeridian=split_antimeridian,
        )
        if seed_cell_polygon.contains(polygon):
            num_edges = 5
            if seed_cell_resolution == 1:
                num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "a5", seed_cell_id, seed_cell_resolution, seed_cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            a5_rows.append(row)
            return a5_rows
        else:
            intersecting_cells = {}  # {cell_u64: cell_polygon}
            covered_cells = set()
            queue = deque([seed_cell_id])

            while queue:
                current_cell_id = queue.popleft()
                if current_cell_id in covered_cells:
                    continue
                covered_cells.add(current_cell_id)

                cell_polygon = a52geo_u64(
                    current_cell_id,
                    options=options,
                    split_antimeridian=split_antimeridian,
                )
                if cell_polygon is None or cell_polygon.is_empty:
                    continue

                if cell_polygon.intersects(polygon):
                    intersecting_cells[current_cell_id] = cell_polygon
                    neighbors = a5.uncompact(
                        a5.grid_disk_vertex(current_cell_id, 1), resolution
                    )
                    for neighbor_id in neighbors:
                        if neighbor_id not in covered_cells:
                            queue.append(neighbor_id)

            for cell_id, cell_polygon in tqdm(
                intersecting_cells.items(), desc="Generating A5 cells", unit=" cells", disable=not verbose
            ):
                if check_predicate(cell_polygon, polygon, predicate):
                    cell_hex = a5.u64_to_hex(cell_id)
                    cell_resolution = a5.get_resolution(cell_id)
                    num_edges = 5
                    if cell_resolution == 1:
                        num_edges = 3
                    row = geodesic_dggs_to_geoseries(
                        "a5", cell_hex, cell_resolution, cell_polygon, num_edges
                    )
                    if include_properties and feature_properties:
                        row.update(feature_properties)
                    a5_rows.append(row)

            # Apply compact mode if enabled
            if compact and a5_rows:
                temp_gdf = gpd.GeoDataFrame(
                    a5_rows, geometry="geometry", crs="EPSG:4326"
                )
                compacted_gdf = a5compact(temp_gdf, a5_hex="a5", output_format="gpd")
                if compacted_gdf is not None:
                    a5_rows = compacted_gdf.to_dict("records")

    return a5_rows


def polygon2a5_new(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    options=None,
    split_antimeridian=False,
):
    """
    Convert a polygon geometry to A5 grid cells.

    Uses ``a5.polygon_to_cells`` (dense boundary sampling + flood fill, with
    center-point containment) to discover all A5 cells covering each ring,
    then applies the requested spatial predicate against each cell polygon,
    with optional compaction at the end.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): A5 resolution [0..28]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable A5 compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2a5)
        include_properties (bool, optional): Whether to include properties in output
        options (dict, optional): Options for a52geo.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.

    Returns:
        list: List of dictionaries representing A5 cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])
        >>> cells = polygon2a5(poly, 10, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """
    a5_rows = []
    if feature.geom_type == "Polygon":
        polygons = [feature]
    elif feature.geom_type == "MultiPolygon":
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        if polygon is None or polygon.is_empty:
            continue

        # a5.polygon_to_cells expects an unclosed [(lon, lat), ...] ring
        ring = [(lon, lat) for lon, lat in polygon.exterior.coords[:-1]]
        if len(ring) < 3:
            continue

        try:
            compacted_cells = a5.polygon_to_cells(ring, resolution)
        except Exception:
            compacted_cells = []
        if not compacted_cells:
            continue

        # Expand to the target resolution so each cell can be predicate-filtered
        try:
            candidate_cells = a5.uncompact(compacted_cells, resolution)
        except Exception:
            candidate_cells = list(compacted_cells)

        # First collect cells that pass the predicate check
        filtered_cells = []
        for cell_id in candidate_cells:
            cell_polygon = a52geo_u64(
                cell_id,
                options=options,
                split_antimeridian=split_antimeridian,
            )
            if cell_polygon is None or cell_polygon.is_empty:
                continue
            if not check_predicate(cell_polygon, polygon, predicate):
                continue
            filtered_cells.append(cell_id)

        # Apply compact after predicate check
        if compact and filtered_cells:
            try:
                filtered_cells = a5.compact(filtered_cells)
            except Exception:
                pass

        # Convert filtered/compacted cells to rows
        for cell_id in filtered_cells:
            cell_hex = a5.u64_to_hex(cell_id)
            cell_polygon = a52geo_u64(
                cell_id,
                options=options,
                split_antimeridian=split_antimeridian,
            )
            cell_resolution = a5.get_resolution(cell_id)
            num_edges = 5
            if cell_resolution == 1:
                num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "a5", cell_hex, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            a5_rows.append(row)

    return a5_rows


def geodataframe2a5(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    options=None,
    split_antimeridian=False,
    verbose=True,
):
    """
    Convert a GeoDataFrame to A5 grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): A5 resolution level [0..28]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable A5 compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint A5 cells
        include_properties (bool, optional): Whether to include properties in output
        options (dict, optional): Options for a52geo.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with A5 grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2a5(gdf, 10)
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

            # Find resolution where A5 cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint A5 cells
            if shortest_distance > 0:
                for res in range(min_res, max_res + 1):
                    _, avg_edge_length, _, _ = a5_metrics(res)
                    cell_diameter = avg_edge_length * 1.4
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_a5_resolution(resolution)

    a5_rows = []

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
            a5_rows.extend(
                point2a5(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    topology=topology,  # Topology already processed above
                    include_properties=include_properties,
                    options=options,
                    split_antimeridian=split_antimeridian,
                )
            )

        elif geom.geom_type in ("LineString", "MultiLineString"):
            a5_rows.extend(
                polyline2a5(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                    options=options,
                    split_antimeridian=split_antimeridian,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            a5_rows.extend(
                polygon2a5(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                    options=options,
                    split_antimeridian=split_antimeridian,
                    verbose=verbose,
                )
            )
            # polygon2a5_new only supports predicate "centroid_within"
    return gpd.GeoDataFrame(a5_rows, geometry="geometry", crs="EPSG:4326")


# --- Main vector2a5 function ---
def vector2a5(
    vector_data,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    output_format="gpd",
    include_properties=True,
    options=None,
    split_antimeridian=False,
    verbose=True,
    **kwargs,
):
    """
    Convert vector data to A5 grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): A5 resolution level [0..28]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable A5 compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint A5 cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        options (dict, optional): Options for a52geo.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2a5("data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_a5_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2a5(
        gdf,
        resolution,
        predicate,
        compact,
        topology,
        include_properties,
        options,
        split_antimeridian=split_antimeridian,
        verbose=verbose,
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2a5"
        else:
            output_name = "a5"
    return convert_to_output_format(result, output_format, output_name)


# --- CLI ---
def vector2a5_cli():
    """
    Command-line interface for vector2a5 conversion.

    This function provides a command-line interface for converting vector data to A5 grid cells.
    It supports various input formats and output formats, with options for resolution control,
    spatial predicates, compact mode, and topology preservation.

    Usage:
        python vector2a5.py -i input.geojson -r 10 -f geojson
        python vector2a5.py -i input.shp -r 8 -p intersect -c -t
    """
    parser = argparse.ArgumentParser(description="Convert vector data to A5 grid cells")
    parser.add_argument("-i", "--input", help="Input file path, URL")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"A5 resolution [{min_res}..{max_res}]. Required when topology=False, auto-calculated when topology=True",
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
        help="Enable A5 compact mode for polygons",
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
        help="Output format (default: gpd).",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Apply antimeridian fixing to the resulting polygons",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to a52geo. "
        "Example: '{\"segments\": 1000}'",
    )

    add_verbose_argument(parser)
    args = parser.parse_args()

    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}", file=sys.stderr)
            sys.exit(1)

    try:
        result = vector2a5(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            output_format=args.output_format,
            include_properties=args.include_properties,
            options=options,
            split_antimeridian=args.split_antimeridian,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
        # For file outputs, the utility prints the saved path
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2a5_cli()
