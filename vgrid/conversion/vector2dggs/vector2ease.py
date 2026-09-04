"""
Vector to EASE Module

This module provides functionality to convert vector geometries to EASE grid cells with flexible input and output formats.

Key Functions:
    point2ease: Convert point geometries to EASE cells
    polyline2ease: Convert line geometries to EASE cells
    polygon2ease: Convert polygon geometries to EASE cells with spatial predicates
    geodataframe2ease: Convert GeoDataFrame to EASE cells with topology support
    vector2ease: Main function for converting various input formats to EASE cells
    vector2ease_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
import math
from tqdm import tqdm
from shapely.geometry import box, MultiPoint
import geopandas as gpd
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from ease_dggs.constants import levels_specs, geo_crs, ease_crs
from ease_dggs.dggs.grid_addressing import geo_polygon_to_grid_ids
from vgrid.conversion.dggscompact.easecompact import ease_compact
from vgrid.utils.geometry import (
    check_predicate,
    shortest_point_distance,
    get_ease_resolution,
)
from vgrid.utils.io import (
    validate_ease_resolution,
    process_input_data_vector,
    convert_to_output_format,
    DGGS_TYPES,
    add_verbose_argument,
    add_compact_depth_argument,
)
from vgrid.conversion.dggs2geo.ease2geo import ease2geo
from vgrid.conversion.latlon2dggs import latlon2ease
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS

min_res = DGGS_TYPES["ease"]["min_res"]
max_res = DGGS_TYPES["ease"]["max_res"]


def point2ease(
    feature,
    resolution,
    feature_properties=None,
    include_properties=True,
):
    """
    Convert a point geometry to EASE grid cells.

    Converts point or multipoint geometries to EASE grid cells at the specified resolution.
    Each point is assigned to its containing EASE cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to EASE cells.
    resolution : int
        EASE resolution level [0..6].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable EASE compact mode (not used for points).
    include_properties : bool, optional
        Whether to include properties in output.

    Returns
    -------
    list of dict
        List of dictionaries representing EASE cells containing the point(s).
        Each dictionary contains EASE cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2ease(point, 4, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2ease(points, 3)
    >>> len(cells)
    2
    """
    ease_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []
    for point in points:
        ease_id = latlon2ease(point.y, point.x, resolution)
        cell_polygon = ease2geo(ease_id)
        num_edges = 4
        row = geodesic_dggs_to_geoseries(
            "ease", ease_id, int(ease_id[1]), cell_polygon, num_edges
        )
        if include_properties and feature_properties:
            row.update(feature_properties)
        ease_rows.append(row)
    return ease_rows


def polyline2ease(
    feature,
    resolution,
    feature_properties=None,
    include_properties=True,
):
    """
    Convert line geometries (LineString, MultiLineString) to EASE grid cells.

    Args:
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Line geometry to convert
        resolution (int): EASE resolution level [0..6]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for lines)
        compact (bool, optional): Enable EASE compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2ease)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing EASE cells intersecting the line
    """
    ease_rows = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        poly_bbox = box(*polyline.bounds)
        polygon_bbox_wkt = poly_bbox.wkt
        cells_bbox = geo_polygon_to_grid_ids(
            polygon_bbox_wkt,
            resolution,
            geo_crs,
            ease_crs,
            levels_specs,
            return_centroids=True,
            wkt_geom=True,
        )
        ease_ids = cells_bbox["result"]["data"]
        for ease_id in ease_ids:
            cell_resolution = int(ease_id[1])
            # Use ease2geo to get the cell geometry
            cell_polygon = ease2geo(ease_id)
            if cell_polygon and cell_polygon.intersects(polyline):
                num_edges = 4
                row = geodesic_dggs_to_geoseries(
                    "ease", str(ease_id), cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                ease_rows.append(row)
    return ease_rows


def polygon2ease(
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
    Convert polygon geometries (Polygon, MultiPolygon) to EASE grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): EASE resolution level [0..6]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable EASE compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2ease)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of GeoJSON feature dictionaries representing EASE cells based on predicate
    """
    ease_rows = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []
    for polygon in polygons:
        poly_bbox = box(*polygon.bounds)
        polygon_bbox_wkt = poly_bbox.wkt
        cells_bbox = geo_polygon_to_grid_ids(
            polygon_bbox_wkt,
            resolution,
            geo_crs,
            ease_crs,
            levels_specs,
            return_centroids=True,
            wkt_geom=True,
        )
        ease_ids = cells_bbox["result"]["data"]
        if not ease_ids:
            continue
        polygon_ease_rows = []
        for ease_id in ease_ids:
            cell_resolution = int(ease_id[1])
            # Use ease2geo to get the cell geometry
            cell_polygon = ease2geo(ease_id)
            if cell_polygon and check_predicate(cell_polygon, polygon, predicate):
                num_edges = 4
                row = geodesic_dggs_to_geoseries(
                    "ease", str(ease_id), cell_resolution, cell_polygon, num_edges
                )
                if feature_properties:
                    row.update(feature_properties)
                polygon_ease_rows.append(row)

        # Compact mode: apply to polygon_ease_rows after predicate check
        if compact:
            # Extract cell IDs from polygon_ease_rows
            cells_to_process = [row.get("ease") for row in polygon_ease_rows]
            # Apply compact
            cells_to_process = ease_compact(cells_to_process, depth=depth, verbose=verbose)
            # Rebuild polygon_ease_rows with compacted cells
            polygon_ease_rows = []
            for cell_id in cells_to_process:
                cell_polygon = ease2geo(cell_id)
                cell_resolution = get_ease_resolution(cell_id)

                # No need to re-check predicate for parent cells from compact mode
                # if not check_predicate(cell_polygon, polygon, predicate):
                #     continue

                num_edges = 4
                row = geodesic_dggs_to_geoseries(
                    "ease", cell_id, cell_resolution, cell_polygon, num_edges
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                polygon_ease_rows.append(row)

        ease_rows.extend(polygon_ease_rows)
    return ease_rows


def geodataframe2ease(
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
    Convert a GeoDataFrame to EASE grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): EASE resolution level [0..6]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable EASE compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint EASE cells
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with EASE grid cells
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

            # Find resolution where EASE cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint EASE cells
            if shortest_distance > 0:
                for res in range(
                    min_res, max_res + 1
                ):  # EASE resolution range is [0..6]
                    if res in levels_specs:
                        cell_width = levels_specs[res]["x_length"]
                        # Use a factor to ensure sufficient separation (cell diagonal is ~1.4x cell width)
                        cell_diagonal = cell_width * math.sqrt(2)
                        if cell_diagonal < shortest_distance:
                            estimated_resolution = res
                            break

        resolution = estimated_resolution

    resolution = validate_ease_resolution(resolution)

    ease_rows = []

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
            ease_rows.extend(
                point2ease(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    include_properties=include_properties,
                )
            )

        elif geom.geom_type in ("LineString", "MultiLineString"):
            ease_rows.extend(
                polyline2ease(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            ease_rows.extend(
                polygon2ease(
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
    return gpd.GeoDataFrame(ease_rows, geometry="geometry", crs="EPSG:4326")


# --- Main vector2ease function ---
def vector2ease(
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
    Convert vector data to EASE grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): EASE resolution level [0..6]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable EASE compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint EASE cells
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
        resolution = validate_ease_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2ease(
        gdf, resolution, predicate, compact, depth, topology, include_properties,
        verbose=verbose,
    )

    output_name = kwargs.get("output_name", None)
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2ease"
        else:
            output_name = "ease"
    return convert_to_output_format(result, output_format, output_name)


# --- CLI ---
def vector2ease_cli():
    """
    Command-line interface for vector2ease conversion.
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to EASE grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"EASE resolution [{min_res}..{max_res}] (0=coarsest, {max_res}=finest)",
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
        help="Enable EASE compact mode for polygons",
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
    add_compact_depth_argument(parser)
    add_verbose_argument(parser)
    args = parser.parse_args()
    args.resolution = validate_ease_resolution(args.resolution)
    output_name = None
    try:
        result = vector2ease(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            topology=args.topology,
            compact=args.compact,
            depth=args.depth,
            output_format=args.output_format,
            output_name=output_name,
            include_properties=args.include_properties,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2ease_cli()
