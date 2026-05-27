"""
Vector to OLC Module

This module provides functionality to convert vector geometries to OLC grid cells with flexible input and output formats.

Key Functions:
    point2olc: Convert point geometries to OLC cells
    polyline2olc: Convert line geometries to OLC cells
    polygon2olc: Convert polygon geometries to OLC cells with spatial predicates
    geodataframe2olc: Convert GeoDataFrame to OLC cells with topology support
    vector2olc: Main function for converting various input formats to OLC cells
    vector2olc_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse
import math
from tqdm import tqdm
from shapely.geometry import MultiPoint
import geopandas as gpd
from vgrid.dggs import olc
from vgrid.generator.olcgrid import olc_grid, olc_refine_cell
from vgrid.utils.geometry import (
    check_predicate,
    shortest_point_distance,
)
from vgrid.stats.olcstats import olc_metrics
from vgrid.utils.io import (
    validate_olc_resolution,
    convert_to_output_format,
    process_input_data_vector,
    DGGS_TYPES,
    olc_resolutions,
)
from vgrid.utils.constants import STRUCTURED_FORMATS, OUTPUT_FORMATS
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.conversion.dggs2geo.olc2geo import olc2geo
from vgrid.conversion.dggscompact.olccompact import olccompact

min_res = DGGS_TYPES["olc"]["min_res"]
max_res = DGGS_TYPES["olc"]["max_res"]


def point2olc(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a point geometry to OLC grid cells.

    Converts point or multipoint geometries to OLC grid cells at the specified resolution.
    Each point is assigned to its containing OLC cell.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to OLC cells.
    resolution : int
        OLC resolution level [2,4,6,8,10,11,12,13,14,15].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable OLC compact mode (not used for points).
    topology : bool, optional
        Enable topology preserving mode (handled by geodataframe2olc).
    include_properties : bool, optional
        Whether to include properties in output.

    Returns
    -------
    list of dict
        List of dictionaries representing OLC cells containing the point(s).
        Each dictionary contains OLC cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2olc(point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2olc(points, 8)
    >>> len(cells)
    2
    """
    olc_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        olc_id = olc.encode(point.y, point.x, resolution)
        cell_polygon = olc2geo(olc_id)
        if cell_polygon:
            olc_row = graticule_dggs_to_geoseries(
                "olc", olc_id, resolution, cell_polygon
            )
            if include_properties and feature_properties:
                olc_row.update(feature_properties)
            olc_rows.append(olc_row)
    return olc_rows


def polyline2olc(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a polyline geometry to OLC grid cells.

    Args:
        feature (shapely.geometry.LineString or shapely.geometry.MultiLineString): Polyline geometry to convert
        resolution (int): OLC resolution level [2,4,6,8,10,11,12,13,14,15]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply (not used for polylines)
        compact (bool, optional): Enable OLC compact mode (not used for polylines)
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2olc)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of dictionaries representing OLC cells intersecting the polyline

    Example:
        >>> from shapely.geometry import LineString
        >>> line = LineString([(-122.4194, 37.7749), (-122.4000, 37.7800)])
        >>> cells = polyline2olc(line, 10, {"name": "route"})
        >>> len(cells) > 0
        True
    """
    olc_rows = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        base_resolution = 2
        base_cells = olc_grid(base_resolution, verbose=False)
        seed_cells = []
        for idx, base_cell in base_cells.iterrows():
            base_cell_poly = base_cell["geometry"]
            if polyline.intersects(base_cell_poly):
                seed_cells.append(base_cell)
        refined_features = []
        for seed_cell in seed_cells:
            seed_cell_poly = seed_cell["geometry"]
            if seed_cell_poly.contains(polyline) and resolution == base_resolution:
                refined_features.append(seed_cell)
            else:
                refined_features.extend(
                    olc_refine_cell(
                        seed_cell_poly.bounds, base_resolution, resolution, polyline
                    )
                )
        # refined_features may be a mix of GeoDataFrame rows and dicts from refine_cell
        # Normalize all to dicts for downstream processing
        normalized_features = []
        for feat in refined_features:
            if isinstance(feat, dict):
                normalized_features.append(feat)
            else:
                # Convert GeoDataFrame row to dict
                d = dict(feat)
                d["geometry"] = feat["geometry"]
                normalized_features.append(d)
        resolution_features = [
            refined_feature
            for refined_feature in normalized_features
            if refined_feature["resolution"] == resolution
        ]
        seen_olc_codes = set()
        for resolution_feature in resolution_features:
            olc_id = resolution_feature["olc"]
            if olc_id not in seen_olc_codes:
                cell_polygon = olc2geo(olc_id)
                olc_row = graticule_dggs_to_geoseries(
                    "olc", olc_id, resolution, cell_polygon
                )
                if include_properties and feature_properties:
                    olc_row.update(feature_properties)
                olc_rows.append(olc_row)
                seen_olc_codes.add(olc_id)
    return olc_rows


def polygon2olc(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a polygon geometry to OLC grid cells.

    Args:
        feature (shapely.geometry.Polygon or shapely.geometry.MultiPolygon): Polygon geometry to convert
        resolution (int): OLC resolution level [2,4,6,8,10,11,12,13,14,15]
        feature_properties (dict, optional): Properties to include in output features
        predicate (str, optional): Spatial predicate to apply ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable OLC compact mode to reduce cell count
        topology (bool, optional): Enable topology preserving mode (handled by geodataframe2olc)
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        list: List of dictionaries representing OLC cells based on predicate

    Example:
        >>> from shapely.geometry import Polygon
        >>> poly = Polygon([(-122.5, 37.7), (-122.3, 37.7), (-122.3, 37.9), (-122.5, 37.9)])
        >>> cells = polygon2olc(poly, 10, {"name": "area"}, predicate="intersect")
        >>> len(cells) > 0
        True
    """
    olc_rows = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        base_resolution = 2
        base_cells = olc_grid(base_resolution, verbose=False)
        seed_cells = []
        for idx, base_cell in base_cells.iterrows():
            base_cell_poly = base_cell["geometry"]
            if polygon.intersects(base_cell_poly):
                seed_cells.append(base_cell)
        refined_features = []
        for seed_cell in seed_cells:
            seed_cell_poly = seed_cell["geometry"]
            if seed_cell_poly.contains(polygon) and resolution == base_resolution:
                refined_features.append(seed_cell)
            else:
                refined_features.extend(
                    olc_refine_cell(
                        seed_cell_poly.bounds, base_resolution, resolution, polygon
                    )
                )
        # refined_features may be a mix of GeoDataFrame rows and dicts from refine_cell
        # Normalize all to dicts for downstream processing
        normalized_features = []
        for feat in refined_features:
            if isinstance(feat, dict):
                normalized_features.append(feat)
            else:
                # Convert GeoDataFrame row to dict
                d = dict(feat)
                d["geometry"] = feat["geometry"]
                normalized_features.append(d)
        resolution_features = [
            refined_feature
            for refined_feature in normalized_features
            if refined_feature["resolution"] == resolution
        ]
        seen_olc_codes = set()
        for resolution_feature in resolution_features:
            olc_id = resolution_feature["olc"]
            if olc_id not in seen_olc_codes:
                cell_geom = olc2geo(olc_id)
                if not check_predicate(cell_geom, polygon, predicate):
                    continue
                olc_row = graticule_dggs_to_geoseries(
                    "olc", olc_id, resolution, cell_geom
                )
                if include_properties and feature_properties:
                    olc_row.update(feature_properties)
                olc_rows.append(olc_row)
                seen_olc_codes.add(olc_id)

    # Apply compact mode if enabled
    if compact and olc_rows:
        # Create a GeoDataFrame from the current results
        temp_gdf = gpd.GeoDataFrame(olc_rows, geometry="geometry", crs="EPSG:4326")

        # Use olccompact function directly
        compacted_gdf = olccompact(temp_gdf, olc_id="olc", output_format="gpd")

        if compacted_gdf is not None:
            # Convert back to list of dictionaries
            olc_rows = compacted_gdf.to_dict("records")
        # If compaction failed, keep original results

    return olc_rows


def geodataframe2olc(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
):
    """
    Convert a GeoDataFrame to OLC grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        resolution (int, optional): OLC resolution level [2,4,6,8,10,11,12,13,14,15]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable OLC compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint OLC cells
        include_properties (bool, optional): Whether to include properties in output

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with OLC grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2olc(gdf, 10)
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

            # Find resolution where OLC cell size is smaller than shortest distance
            # This ensures disjoint points have disjoint OLC cells
            if shortest_distance > 0:
                for res in olc_resolutions:  # OLC valid resolutions
                    _, avg_edge_length, _, _ = olc_metrics(res)
                    cell_diameter = avg_edge_length * math.sqrt(2) 
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_olc_resolution(resolution)

    olc_rows = []

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
            olc_rows.extend(
                point2olc(
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
            olc_rows.extend(
                polyline2olc(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            olc_rows.extend(
                polygon2olc(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                )
            )
    return gpd.GeoDataFrame(olc_rows, geometry="geometry", crs="EPSG:4326")


# --- Main vector2olc function ---
def vector2olc(
    vector_data,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    output_format="gpd",
    include_properties=True,
    **kwargs,
):
    """
    Convert vector data to OLC grid cells from various input formats.

    Args:
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): OLC resolution level [2,4,6,8,10,11,12,13,14,15]. Required when topology=False, auto-calculated when topology=True
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable OLC compact mode for polygons
        topology (bool, optional): Enable topology preserving mode to ensure disjoint features have disjoint OLC cells
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        include_properties (bool, optional): Whether to include properties in output
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2olc("data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    # Validate resolution parameter
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    # Only validate resolution if it's not None
    if resolution is not None:
        resolution = validate_olc_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2olc(
        gdf, resolution, predicate, compact, topology, include_properties
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2olc"
        else:
            output_name = "olc"
    return convert_to_output_format(result, output_format, output_name)


# --- CLI ---
def vector2olc_cli():
    """
    Command-line interface for vector2olc conversion.

    This function provides a command-line interface for converting vector data to OLC grid cells.
    It supports various input formats and output formats, with options for resolution control,
    spatial predicates, compact mode, and topology preservation.

    Usage:
        python vector2olc.py -i input.geojson -r 10 -f geojson
        python vector2olc.py -i input.shp -r 8 -p intersect -c -t
    """
    parser = argparse.ArgumentParser(
        description="Convert vector data to OLC grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=olc_resolutions,
        help=f"OLC resolution {olc_resolutions}. Required when topology=False, auto-calculated when topology=True",
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
        help="Enable OLC compact mode for polygons",
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
    args = parser.parse_args()

    try:
        result = vector2olc(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            output_format=args.output_format,
            include_properties=args.include_properties,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
        # For file outputs, the utility prints the saved path
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2olc_cli()
