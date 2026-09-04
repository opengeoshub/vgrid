"""
Vector to RHEALPix Module

This module provides functionality to convert vector geometries to RHEALPix grid cells with flexible input and output formats.

Key Functions:
    point2rhealpix: Convert point geometries to RHEALPix cells
    polyline2rhealpix: Convert line geometries to RHEALPix cells
    polygon2rhealpix: Convert polygon geometries to RHEALPix cells with spatial predicates
    geodataframe2rhealpix: Convert GeoDataFrame to RHEALPix cells with topology support
    vector2rhealpix: Main function for converting various input formats to RHEALPix cells
    vector2rhealpix_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
import math
from collections import deque
from shapely.geometry import MultiPoint
import geopandas as gpd
from tqdm import tqdm
from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.dggs.rhealpixdggs.rhp_wrappers import linetrace
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.conversion.dggscompact.rhealpixcompact import rhealpix_compact
from vgrid.conversion.dggs2geo.rhealpix2geo import rhealpix2geo
from vgrid.utils.geometry import check_predicate
from vgrid.stats.rhealpixstats import rhealpix_metrics
from vgrid.utils.geometry import (
    shortest_point_distance,
)
from vgrid.utils.io import (
    validate_rhealpix_resolution,
    process_input_data_vector,
    convert_to_output_format,
    add_verbose_argument,
    add_compact_depth_argument,
)
from vgrid.utils.constants import STRUCTURED_FORMATS, OUTPUT_FORMATS
from vgrid.utils.io import DGGS_TYPES

min_res = DGGS_TYPES["rhealpix"]["min_res"]
max_res = DGGS_TYPES["rhealpix"]["max_res"]
rhealpix_dggs = RHEALPixDGGS()


def point2rhealpix(
    feature,
    resolution,
    feature_properties=None,
    include_properties=True,
    fix_antimeridian=None,
):
    """
    Convert a point geometry to RHEALPix grid cells.

    Converts point or multipoint geometries to RHEALPix grid cells at the specified resolution.
    Each point is assigned to its containing RHEALPix cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to RHEALPix cells.
    resolution : int
        RHEALPix resolution level [0..30].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable RHEALPix compact mode (not used for points).
    include_properties : bool, optional
        Whether to include properties in output.
    fix_antimeridian : str, optional
        Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        Defaults to None when omitted.

    Returns
    -------
    list of dict
        List of dictionaries representing RHEALPix cells containing the point(s).
        Each dictionary contains RHEALPix cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2rhealpix(point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2rhealpix(points, 8)
    >>> len(cells)
    2
    """
    rhealpix_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        seed_cell = rhealpix_dggs.cell_from_point(
            resolution, (point.x, point.y), plane=False
        )

        seed_cell_id = str(seed_cell)
        seed_cell_polygon = rhealpix2geo(
            seed_cell_id, fix_antimeridian=fix_antimeridian
        )
        if seed_cell_polygon:
            num_edges = 4
            if seed_cell.ellipsoidal_shape() == "dart":
                num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "rhealpix", seed_cell_id, resolution, seed_cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            rhealpix_rows.append(row)
    return rhealpix_rows


def _rhealpix_cell_from_id(cell_id):
    """Return an rHEALPix Cell from its string id (e.g. ``N45``)."""
    rhealpix_uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
    return rhealpix_dggs.cell(rhealpix_uids)


def _rhealpix_rows_from_cell_ids(
    cell_ids,
    resolution,
    feature_properties,
    include_properties,
    fix_antimeridian,
):
    """Build geoseries rows from rHEALPix cell id strings (polyline / linetrace)."""
    rows = []
    for cell_id in cell_ids:
        cell_polygon = rhealpix2geo(cell_id, fix_antimeridian=fix_antimeridian)
        if cell_polygon is None or cell_polygon.is_empty:
            continue
        cell = _rhealpix_cell_from_id(cell_id)
        num_edges = 4
        if cell.ellipsoidal_shape() == "dart":
            num_edges = 3
        row = geodesic_dggs_to_geoseries(
            "rhealpix", cell_id, resolution, cell_polygon, num_edges
        )
        if include_properties and feature_properties:
            row.update(feature_properties)
        rows.append(row)
    return rows


def polyline2rhealpix(
    feature,
    resolution,
    feature_properties=None,
    include_properties=True,
    fix_antimeridian=None,
):
    """
    Convert a polyline geometry to rHEALPix grid cells using ``linetrace``.

    Cells are those touched by the line at the requested resolution
    (see ``rhp_wrappers.linetrace``).

    Parameters
    ----------
    feature : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Polyline geometry to convert.
    resolution : int
        rHEALPix resolution level [0..15].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Not used for polylines (kept for API compatibility).
    compact : bool, optional
        Not used for polylines (kept for API compatibility).
    topology : bool, optional
        Handled by geodataframe2rhealpix.
    include_properties : bool, optional
        Whether to include properties in output.
    fix_antimeridian : str, optional
        Antimeridian fixing method.

    Returns
    -------
    list of dict
        rHEALPix cell rows along the polyline.

    Examples
    --------
    >>> from shapely.geometry import LineString
    >>> line = LineString([(-122.4194, 37.7749), (-122.4000, 37.7800)])
    >>> cells = polyline2rhealpix(line, 10, {"name": "route"})
    >>> len(cells) > 0
    True
    """
    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []

    seen_ids = set()
    ordered_cell_ids = []
    for polyline in polylines:
        if polyline.is_empty or polyline.length == 0:
            continue
        traced = linetrace(polyline, resolution, plane=False, dggs=rhealpix_dggs)
        if not traced:
            continue
        for cell_id in traced:
            if cell_id in seen_ids:
                continue
            seen_ids.add(cell_id)
            ordered_cell_ids.append(cell_id)

    return _rhealpix_rows_from_cell_ids(
        ordered_cell_ids,
        resolution,
        feature_properties,
        include_properties,
        fix_antimeridian,
    )


def polygon2rhealpix(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    depth=-1,
    include_properties=True,
    fix_antimeridian=None,
    verbose=True,
):
    """
    Convert a polygon geometry to rHEALPix grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): rHEALPix resolution level [0..30]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable rHEALPix compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2rhealpix)
        include_properties (bool, optional): Whether to include properties in output
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        Defaults to None when omitted.

    Returns:
        list: List of dictionaries representing rHEALPix cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])
        >>> cells = polygon2rhealpix(poly, 10, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """
    rhealpix_rows = []
    polygons = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)

    for polygon in polygons:
        rep_pt = polygon.representative_point()
        seed_point = (rep_pt.x, rep_pt.y)
        seed_cell = rhealpix_dggs.cell_from_point(resolution, seed_point, plane=False)
        seed_cell_id = str(seed_cell)
        seed_cell_polygon = rhealpix2geo(
            seed_cell_id, fix_antimeridian=fix_antimeridian
        )
        if seed_cell_polygon.contains(polygon):
            num_edges = 4
            if seed_cell.ellipsoidal_shape() == "dart":
                num_edges = 3
            cell_resolution = resolution
            row = geodesic_dggs_to_geoseries(
                "rhealpix", seed_cell_id, cell_resolution, seed_cell_polygon, num_edges
            )
            if include_properties and feature_properties:
                row.update(feature_properties)
            rhealpix_rows.append(row)
            return rhealpix_rows
        else:
            covered_cells = set()
            queue = deque([seed_cell])  # Use deque for BFS

            while queue:
                current_cell = queue.popleft()  # BFS: FIFO
                current_cell_id = str(current_cell)
                if current_cell_id in covered_cells:
                    continue
                covered_cells.add(current_cell_id)

                cell_polygon = rhealpix2geo(
                    current_cell_id, fix_antimeridian=fix_antimeridian
                )

                if not cell_polygon.intersects(polygon):
                    continue

                neighbors = current_cell.neighbors(plane=False)
                for _, neighbor in neighbors.items():
                    neighbor_id = str(neighbor)
                    if neighbor_id not in covered_cells:
                        queue.append(neighbor)

            for cell_id in covered_cells:
                cell_polygon = rhealpix2geo(cell_id, fix_antimeridian=fix_antimeridian)
                rhealpix_uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
                rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
                cell_resolution = rhealpix_cell.resolution

                if not check_predicate(cell_polygon, polygon, predicate):
                    continue

                num_edges = 4
                if rhealpix_cell.ellipsoidal_shape() == "dart":
                    num_edges = 3
                row = geodesic_dggs_to_geoseries(
                    "rhealpix", cell_id, cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                rhealpix_rows.append(row)

            # Compact mode: apply to rhealpix_rows after predicate check
            if compact:
                # Extract cell IDs from rhealpix_rows
                cells_to_process = [row.get("rhealpix") for row in rhealpix_rows]
                # Apply compact
                cells_to_process = rhealpix_compact(cells_to_process, depth=depth, verbose=verbose)
                # Rebuild rhealpix_rows with compacted cells
                rhealpix_rows = []
                for cell_id in cells_to_process:
                    cell_polygon = rhealpix2geo(
                        cell_id, fix_antimeridian=fix_antimeridian
                    )
                    rhealpix_uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
                    rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
                    cell_resolution = rhealpix_cell.resolution

                    num_edges = 4
                    if rhealpix_cell.ellipsoidal_shape() == "dart":
                        num_edges = 3
                    row = geodesic_dggs_to_geoseries(
                        "rhealpix", cell_id, cell_resolution, cell_polygon, num_edges
                    )
                    if include_properties and feature_properties:
                        row.update(feature_properties)
                    rhealpix_rows.append(row)

    return rhealpix_rows


def geodataframe2rhealpix(
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
    Convert a GeoDataFrame to rHEALPix grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): rHEALPix resolution level [0..30]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable rHEALPix compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint rHEALPix cells
        include_properties (bool, optional): Whether to include properties in output
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        Defaults to None when omitted.

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with rHEALPix grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2rhealpix(gdf, 10)
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

            # Find resolution where rHEALPix cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint rHEALPix cells
            if shortest_distance > 0:
                for res in range(min_res, max_res + 1):
                    _, avg_edge_length, _, _ = rhealpix_metrics(res)
                    cell_diameter = avg_edge_length * math.sqrt(2)
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_rhealpix_resolution(resolution)

    rhealpix_rows = []

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
            rhealpix_rows.extend(
                point2rhealpix(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )

        elif geom.geom_type in ("LineString", "MultiLineString"):
            rhealpix_rows.extend(
                polyline2rhealpix(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            rhealpix_rows.extend(
                polygon2rhealpix(
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
            #   void using native rhp polyfill because it only supports "within" predicate
    return gpd.GeoDataFrame(rhealpix_rows, geometry="geometry", crs="EPSG:4326")


def vector2rhealpix(
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
    Convert vector data to rHEALPix grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): rHEALPix resolution level [0..30]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable rHEALPix compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint rHEALPix cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        Defaults to None when omitted.
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2rhealpix("data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_rhealpix_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2rhealpix(
        gdf,
        resolution,
        predicate,
        compact,
        depth,
        topology,
        include_properties,
        fix_antimeridian=fix_antimeridian,
        verbose=verbose,
    )

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2rhealpix"
        else:
            output_name = "rhealpix"
    return convert_to_output_format(result, output_format, output_name)


def vector2rhealpix_cli():
    """
    Command-line interface for vector2rhealpix conversion.

    This function provides a command-line interface for converting vector data to rHEALPix grid cells.
    It supports various input formats and output formats, with options for resolution control,
    spatial predicates, compact mode, and topology preservation.

    Usage:
        python vector2rhealpix.py -i input.geojson -r 10 -f geojson
        python vector2rhealpix.py -i input.shp -r 8 -p intersect -c -t
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to rHEALPix grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"rHEALPix resolution [{min_res}..{max_res}]. Required when topology=False, auto-calculated when topology=True",
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
        help="Enable rHEALPix compact mode for polygons",
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
        "--fix-antimeridian",
        type=str,
        choices=[
            "shift",
            "shift_balanced",
            "shift_west",
            "shift_east",
            "split",
            "none",
        ],
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    add_compact_depth_argument(parser)
    add_verbose_argument(parser)
    args = parser.parse_args()

    try:
        result = vector2rhealpix(
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
    vector2rhealpix_cli()
