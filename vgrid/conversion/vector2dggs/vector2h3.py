"""
Vector to H3 Module

This module provides functionality to convert vector geometries to H3 grid cells with flexible input and output formats.

Key Functions:
    point2h3: Convert point geometries to H3 cells
    polyline2h3: Convert line geometries to H3 cells
    polygon2h3: Convert polygon geometries to H3 cells with spatial predicates
    geodataframe2h3: Convert GeoDataFrame to H3 cells with topology support
    vector2h3: Main function for converting various input formats to H3 cells
    vector2h3_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
from tqdm import tqdm
import geopandas as gpd
from shapely.geometry import box, MultiPoint
from vgrid.conversion.dggs2geo.h32geo import h32geo
import h3
from vgrid.utils.geometry import (
    geodesic_buffer,
    check_predicate,
    shortest_point_distance,
    geodesic_dggs_to_geoseries,
)
from vgrid.utils.io import (
    process_input_data_vector,
    convert_to_output_format,
    DGGS_TYPES,
)
from vgrid.utils.io import validate_h3_resolution
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS

min_res = DGGS_TYPES["h3"]["min_res"]
max_res = DGGS_TYPES["h3"]["max_res"]


# Function to generate grid for Point
# --- Replace geojson feature output with geoseries dict output ---
def point2h3(
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
    Convert a point geometry to H3 grid cells.

    Converts point or multipoint geometries to H3 grid cells at the specified resolution.
    Each point is assigned to its containing H3 cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to H3 cells.
    resolution : int
        H3 resolution level [0..15].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable H3 compact mode (not used for points).
    topology : bool, optional
        Enable topology preserving mode (handled by geodataframe2h3).
    include_properties : bool, optional
        Whether to include properties in output.

    Returns
    -------
    list of dict
        List of dictionaries representing H3 cells containing the point(s).
        Each dictionary contains H3 cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2h3(point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2h3(points, 8)
    >>> len(cells)
    2
    """

    h3_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []
    for point in points:
        h3_id = h3.latlng_to_cell(point.y, point.x, resolution)
        cell_polygon = h32geo(h3_id, fix_antimeridian=fix_antimeridian)
        cell_resolution = h3.get_resolution(h3_id)
        num_edges = 6
        if h3.is_pentagon(h3_id):
            num_edges = 5
        row = geodesic_dggs_to_geoseries(
            "h3", h3_id, cell_resolution, cell_polygon, num_edges
        )

        # Add properties if requested
        if include_properties and feature_properties:
            row.update(feature_properties)

        h3_rows.append(row)
    return h3_rows


def polyline2h3(
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
    Convert each polyline to an H3 path using all vertices.

    For each LineString or MultiLineString, this function:
    - Iterates over every consecutive vertex pair
    - Converts each pair to H3 cells
    - Uses h3.grid_path_cells(start_cell, end_cell) per segment
    - Merges segment paths into a single ordered, de-duplicated list of cells
    - Converts those cells back to geometries and appends them to h3_rows
    """

    h3_rows = []

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

        ordered_cells = []
        for i in range(len(coords) - 1):
            start_x, start_y = coords[i]
            end_x, end_y = coords[i + 1]

            start_cell = h3.latlng_to_cell(start_y, start_x, resolution)
            end_cell = h3.latlng_to_cell(end_y, end_x, resolution)

            try:
                segment_cells = h3.grid_path_cells(start_cell, end_cell)
            except Exception:
                segment_cells = []

            for cell_id in segment_cells:
                # Preserve each segment path exactly; only drop immediate duplicates
                # at segment boundaries.
                if ordered_cells and ordered_cells[-1] == cell_id:
                    continue
                ordered_cells.append(cell_id)

        for cell_id in ordered_cells:
            cell_polygon = h32geo(cell_id, fix_antimeridian=fix_antimeridian)
            cell_resolution = h3.get_resolution(cell_id)
            num_edges = 6
            if h3.is_pentagon(cell_id):
                num_edges = 5

            row = geodesic_dggs_to_geoseries(
                "h3", cell_id, cell_resolution, cell_polygon, num_edges
            )

            if include_properties and feature_properties:
                row.update(feature_properties)

            h3_rows.append(row)

    return h3_rows


# --- Polygon ---
def polygon2h3(
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
    Convert a polygon geometry to H3 grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): H3 resolution [0..15]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable H3 compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2h3)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of dictionaries representing H3 cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])
        >>> cells = polygon2h3(poly, 10, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """

    h3_rows = []
    if feature.geom_type == "Polygon":
        polygons = [feature]
    elif feature.geom_type == "MultiPolygon":
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        bbox = box(*polygon.bounds)
        bbox_cells = h3.geo_to_cells(bbox, resolution)

        # First collect cells that pass the predicate check
        filtered_cells = []
        for bbox_cell in bbox_cells:
            cell_polygon = h32geo(bbox_cell, fix_antimeridian=fix_antimeridian)
            # Use the check_predicate function to determine if we should keep this cell
            if not check_predicate(cell_polygon, polygon, predicate):
                continue  # Skip non-matching cells
            filtered_cells.append(bbox_cell)

        # Apply compact after predicate check
        if compact:
            filtered_cells = h3.compact_cells(filtered_cells)

        # Convert filtered/compacted cells to rows
        for cell_id in filtered_cells:
            cell_polygon = h32geo(cell_id, fix_antimeridian=fix_antimeridian)
            cell_resolution = h3.get_resolution(cell_id)
            num_edges = 6
            if h3.is_pentagon(cell_id):
                num_edges = 5
            row = geodesic_dggs_to_geoseries(
                "h3", cell_id, cell_resolution, cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            h3_rows.append(row)

    return h3_rows


def geodataframe2h3(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    fix_antimeridian=None,
):
    """
    Convert a GeoDataFrame to H3 grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): H3 resolution [0..15]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable H3 compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint H3 cells
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with H3 grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2h3(gdf, 10)
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

            # Find resolution where H3 cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint H3 cells
            if shortest_distance > 0:
                for res in range(min_res, max_res + 1):
                    avg_edge_length = h3.average_hexagon_edge_length(res, unit="m")
                    # Use a factor to ensure sufficient separation (cell diameter is ~2x edge length)
                    cell_diameter = avg_edge_length * 2
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_h3_resolution(resolution)

    h3_rows = []

    for _, row in tqdm(gdf.iterrows(), desc="Processing features", total=len(gdf)):
        geom = row.geometry
        if geom is None:
            continue

        props = row.to_dict()
        if "geometry" in props:
            del props["geometry"]

        if not include_properties:
            props = {}

        if geom.geom_type == "Point" or geom.geom_type == "MultiPoint":
            h3_rows.extend(
                point2h3(
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
            h3_rows.extend(
                polyline2h3(
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
            h3_rows.extend(
                polygon2h3(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )
    return gpd.GeoDataFrame(h3_rows, geometry="geometry", crs="EPSG:4326")


def vector2h3(
    vector_data,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    output_format="gpd",
    include_properties=True,
    fix_antimeridian=None,
    **kwargs,
):
    """
    Convert vector data to H3 grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): H3 resolution [0..15]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable H3 compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint H3 cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2h3("data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_h3_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2h3(
        gdf,
        resolution,
        predicate,
        compact,
        topology,
        include_properties,
        fix_antimeridian=fix_antimeridian,
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2h3"
        else:
            output_name = "h3"
    return convert_to_output_format(result, output_format, output_name)


def vector2h3_cli():
    """
    Command-line interface for vector2h3 conversion.

    This function provides a command-line interface for converting vector data to H3 grid cells.
    It supports various input formats and output formats, with options for resolution control,
    spatial predicates, compact mode, and topology preservation.

    Usage:
        python vector2h3.py -i input.geojson -r 10 -f geojson
        python vector2h3.py -i input.shp -r 8 -p intersect -c -t
    """
    parser = argparse.ArgumentParser(description="Convert vector data to H3 grid cells")
    parser.add_argument("-i", "--input", help="Input file path, URL")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"H3 resolution [{min_res}..{max_res}]. Required when topology=False, auto-calculated when topology=True",
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
        help="Enable H3 compact mode for polygons",
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
    try:
        result = vector2h3(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            output_format=args.output_format,
            include_properties=args.include_properties,
            fix_antimeridian=fix_antimeridian,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
        # For file outputs, the utility prints the saved path
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2h3_cli()
