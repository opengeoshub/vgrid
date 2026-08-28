"""
Vector to Quadkey Module

This module provides functionality to convert vector geometries to Quadkey grid cells with flexible input and output formats.

Key Functions:
    point2quadkey: Convert point geometries to Quadkey cells
    polyline2quadkey: Convert line geometries to Quadkey cells
    polygon2quadkey: Convert polygon geometries to Quadkey cells with spatial predicates
    geodataframe2quadkey: Convert GeoDataFrame to Quadkey cells with topology support
    vector2quadkey: Main function for converting various input formats to Quadkey cells
    vector2quadkey_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
from tqdm import tqdm
from shapely.geometry import Polygon, MultiPoint
import geopandas as gpd
from vgrid.dggs import tilecode
from vgrid.dggs import mercantile
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.conversion.dggscompact.quadkeycompact import quadkeycompact
from vgrid.stats.quadkeystats import quadkey_metrics
from vgrid.utils.geometry import (
    check_predicate,
    shortest_point_distance,
)
from math import sqrt
from vgrid.utils.io import (
    validate_quadkey_resolution,
    process_input_data_vector,
    convert_to_output_format,
    DGGS_TYPES,
    add_verbose_argument,
)
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS

min_res = DGGS_TYPES["quadkey"]["min_res"]
max_res = DGGS_TYPES["quadkey"]["max_res"]


def point2quadkey(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a point geometry to Quadkey grid cells.

    Converts point or multipoint geometries to Quadkey grid cells at the specified resolution.
    Each point is assigned to its containing Quadkey cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to Quadkey cells.
    resolution : int
        Quadkey resolution level [0..29].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable Quadkey compact mode (not used for points).
    topology : bool, optional
        Enable topology preserving mode (handled by geodataframe2quadkey).
    include_properties : bool, optional
        Whether to include properties in output.

    Returns
    -------
    list of dict
        List of dictionaries representing Quadkey cells containing the point(s).
        Each dictionary contains Quadkey cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2quadkey(point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2quadkey(points, 8)
    >>> len(cells)
    2
    """
    quadkey_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        quadkey_id = tilecode.latlon2quadkey(point.y, point.x, resolution)
        quadkey_cell = mercantile.tile(point.x, point.y, resolution)
        bounds = mercantile.bounds(quadkey_cell)
        if bounds:
            min_lat, min_lon = bounds.south, bounds.west
            max_lat, max_lon = bounds.north, bounds.east
            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat],
                ]
            )
            quadkey_row = graticule_dggs_to_geoseries(
                "quadkey", quadkey_id, resolution, cell_polygon
            )
            if include_properties and feature_properties:
                quadkey_row.update(feature_properties)
            quadkey_rows.append(quadkey_row)
    return quadkey_rows


def polyline2quadkey(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a polyline geometry to Quadkey grid cells.

    Args:
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Polyline geometry to convert
        resolution (int): Quadkey resolution level [0..29]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for polylines)
        compact (bool, optional): Enable Quadkey compact mode (not used for polylines)
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2quadkey)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of dictionaries representing Quadkey cells intersecting the polyline

    Example:
        >>> from shapely.geometry import LineString
        >>> line = LineString([(-122.4194, 37.7749), (-122.4000, 37.7800)])
        >>> cells = polyline2quadkey(line, 10, {"name": "route"})
        >>> len(cells) > 0
        True
    """
    quadkey_rows = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        min_lon, min_lat, max_lon, max_lat = polyline.bounds
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        for tile in tiles:
            z, x, y = tile.z, tile.x, tile.y
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                quadkey_id = mercantile.quadkey(tile)
                cell_polygon = Polygon(
                    [
                        [min_lon, min_lat],
                        [max_lon, min_lat],
                        [max_lon, max_lat],
                        [min_lon, max_lat],
                        [min_lon, min_lat],
                    ]
                )
                if cell_polygon.intersects(polyline):
                    quadkey_row = graticule_dggs_to_geoseries(
                        "quadkey", quadkey_id, resolution, cell_polygon
                    )
                    if include_properties and feature_properties:
                        quadkey_row.update(feature_properties)
                    quadkey_rows.append(quadkey_row)

    return quadkey_rows


def polygon2quadkey(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a polygon geometry to Quadkey grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): Quadkey resolution level [0..29]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable Quadkey compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2quadkey)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of dictionaries representing Quadkey cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])
        >>> cells = polygon2quadkey(poly, 10, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """
    quadkey_rows = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        for tile in tiles:
            z, x, y = tile.z, tile.x, tile.y
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                quadkey_id = mercantile.quadkey(tile)
                cell_polygon = Polygon(
                    [
                        [min_lon, min_lat],
                        [max_lon, min_lat],
                        [max_lon, max_lat],
                        [min_lon, max_lat],
                        [min_lon, min_lat],
                    ]
                )
                if check_predicate(cell_polygon, polygon, predicate):
                    quadkey_row = graticule_dggs_to_geoseries(
                        "quadkey", quadkey_id, resolution, cell_polygon
                    )
                    if include_properties and feature_properties:
                        quadkey_row.update(feature_properties)
                    quadkey_rows.append(quadkey_row)

    # Apply compact mode if enabled
    if compact and quadkey_rows:
        # Create a GeoDataFrame from the current results
        temp_gdf = gpd.GeoDataFrame(quadkey_rows, geometry="geometry", crs="EPSG:4326")

        # Use quadkeycompact function directly
        compacted_gdf = quadkeycompact(
            temp_gdf, quadkey_id="quadkey", output_format="gpd"
        )

        if compacted_gdf is not None:
            # Convert back to list of dictionaries
            quadkey_rows = compacted_gdf.to_dict("records")
        # If compaction failed, keep original results

    return quadkey_rows


def geodataframe2quadkey(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    verbose=True,
):
    """
    Convert a GeoDataFrame to Quadkey grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): Quadkey resolution level [0..29]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable Quadkey compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint Quadkey cells
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with Quadkey grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2quadkey(gdf, 10)
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

            # Find resolution where Quadkey cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint Quadkey cells
            if shortest_distance > 0:
                for res in range(min_res, max_res + 1):
                    _, avg_edge_len, _, _ = quadkey_metrics(res, unit="m")
                    cell_diameter = avg_edge_len * sqrt(2)
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_quadkey_resolution(resolution)

    quadkey_rows = []

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
            quadkey_rows.extend(
                point2quadkey(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    topology=topology,  # Topology already processed above
                    include_properties=include_properties,
                )
            )

        elif geom.geom_type in ("LineString", "MultiLineString"):
            quadkey_rows.extend(
                polyline2quadkey(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            quadkey_rows.extend(
                polygon2quadkey(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                )
            )
    return gpd.GeoDataFrame(quadkey_rows, geometry="geometry", crs="EPSG:4326")


# --- Main vector2quadkey function ---
def vector2quadkey(
    vector_data,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    output_format="gpd",
    include_properties=True,
    verbose=True,
    **kwargs,
):
    """
    Convert vector data to Quadkey grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): Quadkey resolution level [0..29]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable Quadkey compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint Quadkey cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2quadkey("data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_quadkey_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2quadkey(
        gdf, resolution, predicate, compact, topology, include_properties,
        verbose=verbose,
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2quadkey"
        else:
            output_name = "quadkey"
    return convert_to_output_format(result, output_format, output_name)


# --- CLI ---
def vector2quadkey_cli():
    """
    Command-line interface for vector2quadkey conversion.

    This function provides a command-line interface for converting vector data to Quadkey grid cells.
    It supports various input formats and output formats, with options for resolution control,
    spatial predicates, compact mode, and topology preservation.

    Usage:
        python vector2quadkey.py -i input.geojson -r 10 -f geojson
        python vector2quadkey.py -i input.shp -r 8 -p intersect -c -t
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to Quadkey grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"Quadkey resolution [{min_res}..{max_res}]. Required when topology=False, auto-calculated when topology=True",
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
        help="Enable Quadkey compact mode for polygons",
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
    add_verbose_argument(parser)
    args = parser.parse_args()

    try:
        result = vector2quadkey(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
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
    vector2quadkey_cli()
