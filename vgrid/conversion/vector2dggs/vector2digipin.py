"""
Vector to DIGIPIN Module

This module provides functionality to convert vector geometries to DIGIPIN grid cells with flexible input and output formats.
DIGIPIN is a hierarchical geocoding system for India that divides geographic areas into a 4x4 grid recursively.

Key Functions:
    point2digipin: Convert point geometries to DIGIPIN cells
    polyline2digipin: Convert line geometries to DIGIPIN cells
    polygon2digipin: Convert polygon geometries to DIGIPIN cells with spatial predicates
    geodataframe2digipin: Convert GeoDataFrame to DIGIPIN cells with topology support
    vector2digipin: Main function for converting various input formats to DIGIPIN cells
    vector2digipin_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
from math import sqrt
from tqdm import tqdm
from shapely.geometry import MultiPoint
import geopandas as gpd
from vgrid.dggs.digipin import BOUNDS
from vgrid.conversion.latlon2dggs import latlon2digipin
from vgrid.conversion.dggs2geo.digipin2geo import digipin2geo
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.conversion.dggscompact.digipincompact import digipincompact
from vgrid.stats.digipinstats import digipin_metrics
from vgrid.utils.geometry import (
    check_predicate,
    shortest_point_distance,
)
from vgrid.utils.io import (
    validate_digipin_resolution,
    process_input_data_vector,
    convert_to_output_format,
    DGGS_TYPES,
    add_verbose_argument,
)
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS

min_res = DGGS_TYPES["digipin"]["min_res"]
max_res = DGGS_TYPES["digipin"]["max_res"]


def point2digipin(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a point geometry to DIGIPIN grid cells.

    Converts point or multipoint geometries to DIGIPIN grid cells at the specified resolution.
    Each point is assigned to its containing DIGIPIN cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to DIGIPIN cells.
    resolution : int
        DIGIPIN resolution level [1..10].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable DIGIPIN compact mode (not used for points).
    topology : bool, optional
        Enable topology preserving mode (handled by geodataframe2digipin).
    include_properties : bool, optional
        Whether to include properties in output.

    Returns
    -------
    list of dict
        List of dictionaries representing DIGIPIN cells containing the point(s).
        Each dictionary contains DIGIPIN cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(77.2090, 28.6139)  # Delhi, India
    >>> cells = point2digipin(point, 10, {"name": "Delhi"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(77.2090, 28.6139), (72.8777, 19.0760)])  # Delhi, Mumbai
    >>> cells = point2digipin(points, 8)
    >>> len(cells)
    2
    """
    digipin_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        digipin_id = latlon2digipin(point.y, point.x, resolution)

        # Skip if point is out of bounds
        if digipin_id == "Out of Bound":
            continue

        cell_polygon = digipin2geo(digipin_id)
        if isinstance(cell_polygon, str):
            continue

        digipin_row = graticule_dggs_to_geoseries(
            "digipin", digipin_id, resolution, cell_polygon
        )
        if include_properties and feature_properties:
            digipin_row.update(feature_properties)
        digipin_rows.append(digipin_row)
    return digipin_rows


def polyline2digipin(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a polyline geometry to DIGIPIN grid cells.

    Args:
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Polyline geometry to convert
        resolution (int): DIGIPIN resolution level [1..10]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for polylines)
        compact (bool, optional): Enable DIGIPIN compact mode (not used for polylines)
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2digipin)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of dictionaries representing DIGIPIN cells intersecting the polyline

    Example:
        >>> from shapely.geometry import LineString
        >>> line = LineString([(77.2090, 28.6139), (77.2200, 28.6200)])  # Delhi area
        >>> cells = polyline2digipin(line, 10, {"name": "route"})
        >>> len(cells) > 0
        True
    """

    digipin_rows = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        min_lon, min_lat, max_lon, max_lat = polyline.bounds

        # Constrain to DIGIPIN bounds (India region)
        min_lat = max(min_lat, BOUNDS["minLat"])
        min_lon = max(min_lon, BOUNDS["minLon"])
        max_lat = min(max_lat, BOUNDS["maxLat"])
        max_lon = min(max_lon, BOUNDS["maxLon"])

        # Calculate sampling density based on resolution (following digipingrid.py approach)
        # Each level divides the cell by 4 (2x2 grid)
        base_width = 9.0  # degrees at resolution 1
        factor = 0.25 ** (resolution - 1)  # each level divides by 4
        sample_width = base_width * factor

        seen_cells = set()

        # Sample points across the bounding box
        lon = min_lon
        while lon <= max_lon:
            lat = min_lat
            while lat <= max_lat:
                try:
                    # Get DIGIPIN code for this point at the specified resolution
                    digipin_id = latlon2digipin(lat, lon, resolution)

                    if digipin_id == "Out of Bound":
                        lat += sample_width
                        continue

                    if digipin_id in seen_cells:
                        lat += sample_width
                        continue

                    seen_cells.add(digipin_id)

                    # Get the bounds for this DIGIPIN cell
                    cell_polygon = digipin2geo(digipin_id)

                    if isinstance(cell_polygon, str):  # Error like 'Invalid DIGIPIN'
                        lat += sample_width
                        continue

                    # Check if cell intersects with polyline
                    if cell_polygon.intersects(polyline):
                        digipin_row = graticule_dggs_to_geoseries(
                            "digipin", digipin_id, resolution, cell_polygon
                        )
                        if include_properties and feature_properties:
                            digipin_row.update(feature_properties)
                        digipin_rows.append(digipin_row)

                except Exception:
                    # Skip cells with errors
                    pass

                lat += sample_width
            lon += sample_width
    return digipin_rows


def polygon2digipin(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    verbose=True,
):
    """
    Convert a polygon geometry to DIGIPIN grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): DIGIPIN resolution level [1..10]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable DIGIPIN compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2digipin)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of dictionaries representing DIGIPIN cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(77.1, 28.5), (77.3, 28.5), (77.3, 28.7), (77.1, 28.7)])  # Delhi area
        >>> cells = polygon2digipin(poly, 10, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """
    from vgrid.dggs.digipin import BOUNDS

    digipin_rows = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        min_lon, min_lat, max_lon, max_lat = polygon.bounds

        # Constrain to DIGIPIN bounds (India region)
        min_lat = max(min_lat, BOUNDS["minLat"])
        min_lon = max(min_lon, BOUNDS["minLon"])
        max_lat = min(max_lat, BOUNDS["maxLat"])
        max_lon = min(max_lon, BOUNDS["maxLon"])

        # Calculate sampling density based on resolution (following digipingrid.py approach)
        # Each level divides the cell by 4 (2x2 grid)
        base_width = 9.0  # degrees at resolution 1
        factor = 0.25 ** (resolution - 1)  # each level divides by 4
        sample_width = base_width * factor

        seen_cells = set()

        # Sample points across the bounding box
        lon = min_lon
        while lon <= max_lon:
            lat = min_lat
            while lat <= max_lat:
                try:
                    # Get DIGIPIN code for this point at the specified resolution
                    digipin_id = latlon2digipin(lat, lon, resolution)

                    if digipin_id == "Out of Bound":
                        lat += sample_width
                        continue

                    if digipin_id in seen_cells:
                        lat += sample_width
                        continue

                    seen_cells.add(digipin_id)

                    # Get the bounds for this DIGIPIN cell
                    cell_polygon = digipin2geo(digipin_id)

                    if isinstance(cell_polygon, str):  # Error like 'Invalid DIGIPIN'
                        lat += sample_width
                        continue
                    # Check spatial predicate
                    if check_predicate(cell_polygon, polygon, predicate):
                        digipin_row = graticule_dggs_to_geoseries(
                            "digipin", digipin_id, resolution, cell_polygon
                        )
                        if include_properties and feature_properties:
                            digipin_row.update(feature_properties)
                        digipin_rows.append(digipin_row)

                except Exception:
                    # Skip cells with errors
                    pass

                lat += sample_width
            lon += sample_width

    # Apply compact mode if enabled
    if compact and digipin_rows:
        # Create a GeoDataFrame from the current results
        temp_gdf = gpd.GeoDataFrame(digipin_rows, geometry="geometry", crs="EPSG:4326")

        # Use digipincompact function directly
        compacted_gdf = digipincompact(
            temp_gdf, digipin_id="digipin", output_format="gpd", verbose=verbose)

        if compacted_gdf is not None:
            # Convert back to list of dictionaries
            digipin_rows = compacted_gdf.to_dict("records")
        # If compaction failed, keep original results

    return digipin_rows


def geodataframe2digipin(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    verbose=True,
):
    """
    Convert a GeoDataFrame to DIGIPIN grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): DIGIPIN resolution level [1..15]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable DIGIPIN compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint DIGIPIN cells
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with DIGIPIN grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['Delhi'],
        ...     'geometry': [Point(77.2090, 28.6139)]
        ... })
        >>> result = geodataframe2digipin(gdf, 10)
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

            # Find resolution where DIGIPIN cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint DIGIPIN cells
            if shortest_distance > 0:
                for res in range(min_res, max_res + 1):
                    _, avg_edge_len, _, _ = digipin_metrics(res, unit="m")
                    cell_diameter = avg_edge_len * sqrt(2)
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break
        resolution = estimated_resolution

    resolution = validate_digipin_resolution(resolution)

    digipin_rows = []

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
            digipin_rows.extend(
                point2digipin(
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
            digipin_rows.extend(
                polyline2digipin(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            digipin_rows.extend(
                polygon2digipin(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                    verbose=verbose,
                )
            )
    return gpd.GeoDataFrame(digipin_rows, geometry="geometry", crs="EPSG:4326")


# --- Main vector2digipin function ---
def vector2digipin(
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
    Convert vector data to DIGIPIN grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): DIGIPIN resolution level [1..15]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable DIGIPIN compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint DIGIPIN cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2digipin("data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_digipin_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2digipin(
        gdf, resolution, predicate, compact, topology, include_properties,
        verbose=verbose,
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2digipin"
        else:
            output_name = "digipin"
    return convert_to_output_format(result, output_format, output_name)


# --- CLI ---
def vector2digipin_cli():
    """
    Command-line interface for vector2digipin conversion.

    This function provides a command-line interface for converting vector data to DIGIPIN grid cells.
    It supports various input formats and output formats, with options for resolution control,
    spatial predicates, compact mode, and topology preservation.

    Usage:
        python vector2digipin.py -i input.geojson -r 10 -f geojson
        python vector2digipin.py -i input.shp -r 8 -p intersect -c -t
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to DIGIPIN grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"DIGIPIN resolution [{min_res}..{max_res}]. Required when topology=False, auto-calculated when topology=True",
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
        help="Enable DIGIPIN compact mode for polygons",
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
        result = vector2digipin(
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
    vector2digipin_cli()
