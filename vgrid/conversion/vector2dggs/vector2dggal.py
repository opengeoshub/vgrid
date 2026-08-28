"""
Vector to DGGAL Module

This module provides functionality to convert vector geometries to DGGAL grid cells with flexible input and output formats.

Key Functions:
    point2dggal: Convert point geometries to DGGAL cells
    polyline2dggal: Convert line geometries to DGGAL cells
    polygon2dggal: Convert polygon geometries to DGGAL cells with spatial predicates
    vector2dggal: Main function for converting various input formats to DGGAL cells
    vector2dggal_cli: Command-line interface for batch processing
"""

from __future__ import annotations

import argparse
import math
import sys
import geopandas as gpd
from shapely.geometry import MultiPoint
from vgrid.utils.io import (
    process_input_data_vector,
    convert_to_output_format,
    validate_dggal_resolution,
    add_verbose_argument,
)
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS, DGGAL_TYPES
from vgrid.utils.geometry import (
    check_predicate,
    geodesic_dggs_to_geoseries,
    shortest_point_distance,
)
from vgrid.stats.dggalstats import dggal_metrics
from vgrid.conversion.latlon2dggs import latlon2dggal
from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo
from vgrid.conversion.dggscompact.dggalcompact import dggalcompact
from tqdm import tqdm

from dggal import *

app = Application(appGlobals=globals())
pydggal_setup(app)

_DGGAL_HEX_TYPES = {
    "isea3h",
    "isea7h",
    "isea7h_z7",
    "ivea3h",
    "ivea7h",
    "ivea7h_z7",
    "rtea3h",
    "rtea7h",
    "rtea7h_z7",
}


def _dggal_cell_diameter(dggs_type: str, resolution: int) -> float:
    """Approximate cell diameter (m) from ``dggal_metrics``."""
    _, avg_edge_len, _, _ = dggal_metrics(dggs_type, resolution, unit="m")
    if dggs_type in _DGGAL_HEX_TYPES:
        return avg_edge_len * 2
    return avg_edge_len * math.sqrt(2)


def point2dggal(
    dggs_type: str,
    feature=None,
    resolution=None,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    split_antimeridian=False,
):
    """
    Convert a point geometry to DGGAL grid cells.

    Converts point or multipoint geometries to DGGAL grid cells at the specified resolution.
    Each point is assigned to its containing DGGAL cell.

    Parameters
    ----------
    dggs_type : str
        DGGAL DGGS type (e.g., "isea3h", "isea4t", "isea7h", "isea9h").
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to DGGAL cells.
    resolution : int
        DGGAL resolution level.
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable DGGAL compact mode (not used for points).
    topology : bool, optional
        Enable topology preserving mode (handled by ``geodataframe2dggal``).
    include_properties : bool, optional
        Whether to include properties in output.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.

    Returns
    -------
    list of dict
        List of dictionaries representing DGGAL cells containing the point(s).
        Each dictionary contains DGGAL cell properties and geometry.

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> cells = point2dggal("isea3h", point, 10, {"name": "SF"})
    >>> len(cells)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> cells = point2dggal("isea4t", points, 8)
    >>> len(cells)
    2
    """
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = globals()[dggs_class_name]()

    dggal_rows = []
    if feature.geom_type in ("Point"):
        points = [feature]
    elif feature.geom_type in ("MultiPoint"):
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        zone_id = latlon2dggal(dggs_type, point.y, point.x, resolution)
        zone = dggrs.getZoneFromTextID(zone_id)
        cell_resolution = dggrs.getZoneLevel(zone)
        num_edges = dggrs.countZoneEdges(zone)
        cell_polygon = dggal2geo(
            dggs_type, zone_id, split_antimeridian=split_antimeridian
        )
        row = geodesic_dggs_to_geoseries(
            f"dggal_{dggs_type}", zone_id, cell_resolution, cell_polygon, num_edges
        )
        # Add properties if requested
        if include_properties and feature_properties:
            row.update(feature_properties)

        dggal_rows.append(row)

    return dggal_rows


def polyline2dggal(
    dggs_type: str,
    feature=None,
    resolution=None,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    split_antimeridian=False,
):
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = globals()[dggs_class_name]()

    dggal_rows = []
    if feature.geom_type in ("LineString"):
        polylines = [feature]
    elif feature.geom_type in ("MultiLineString"):
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        try:
            min_lon, min_lat, max_lon, max_lat = polyline.bounds
            ll = GeoPoint(min_lat, min_lon)
            ur = GeoPoint(max_lat, max_lon)
            geo_extent = GeoExtent(ll, ur)
            zones = dggrs.listZones(resolution, geo_extent)
            for zone in zones:
                zone_id = dggrs.getZoneTextID(zone)
                cell_polygon = dggal2geo(
                    dggs_type, zone_id, split_antimeridian=split_antimeridian
                )
                if not check_predicate(cell_polygon, polyline, "intersects"):
                    continue
                cell_resolution = dggrs.getZoneLevel(zone)
                num_edges = dggrs.countZoneEdges(zone)
                row = geodesic_dggs_to_geoseries(
                    f"dggal_{dggs_type}",
                    zone_id,
                    cell_resolution,
                    cell_polygon,
                    num_edges,
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                dggal_rows.append(row)
        except Exception:
            pass

    return dggal_rows


def polygon2dggal(
    dggs_type: str | None = None,
    feature=None,
    resolution=None,
    feature_properties=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    split_antimeridian=False,
):
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = globals()[dggs_class_name]()

    dggal_rows = []
    if feature.geom_type in ("Polygon"):
        polygons = [feature]
    elif feature.geom_type in ("MultiPolygon"):
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        try:
            min_lon, min_lat, max_lon, max_lat = polygon.bounds
            ll = GeoPoint(min_lat, min_lon)
            ur = GeoPoint(max_lat, max_lon)
            geo_extent = GeoExtent(ll, ur)
            zones = dggrs.listZones(resolution, geo_extent)
            for zone in zones:
                zone_id = dggrs.getZoneTextID(zone)
                cell_polygon = dggal2geo(
                    dggs_type, zone_id, split_antimeridian=split_antimeridian
                )
                if not check_predicate(cell_polygon, polygon, predicate):
                    continue
                cell_resolution = dggrs.getZoneLevel(zone)
                num_edges = dggrs.countZoneEdges(zone)
                row = geodesic_dggs_to_geoseries(
                    f"dggal_{dggs_type}",
                    zone_id,
                    cell_resolution,
                    cell_polygon,
                    num_edges,
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                dggal_rows.append(row)
        except Exception:
            pass

    if compact and dggal_rows:
        # Create a GeoDataFrame from the current results
        temp_gdf = gpd.GeoDataFrame(dggal_rows, geometry="geometry", crs="EPSG:4326")
        # Use a5compact function directly
        compacted_gdf = dggalcompact(dggs_type, temp_gdf, output_format="gpd")

        if compacted_gdf is not None:
            # Convert back to list of dictionaries
            dggal_rows = compacted_gdf.to_dict("records")
        # If compaction failed, keep original results

    return dggal_rows


def geodataframe2dggal(
    dggs_type: str,
    gdf,
    resolution=None,
    predicate: str | None = None,
    compact: bool = False,
    topology: bool = False,
    include_properties: bool = True,
    split_antimeridian: bool = False,
    verbose=True,
):
    """
    Convert a GeoDataFrame to DGGAL grid cells.

    Args:
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        dggs_type (str): One of DGGAL_TYPES
        resolution (int, optional): Integer resolution. Required when topology=False,
            auto-calculated from point spacing when topology=True.
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable DGGAL compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode to ensure disjoint
            points have disjoint DGGAL cells
        include_properties (bool, optional): Whether to include properties in output
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with DGGAL grid cells

    Example:
        >>> import geopandas as gpd
        >>> from shapely.geometry import Point
        >>> gdf = gpd.GeoDataFrame({
        ...     'name': ['San Francisco'],
        ...     'geometry': [Point(-122.4194, 37.7749)]
        ... })
        >>> result = geodataframe2dggal("isea3h", gdf, 10)
        >>> len(result) > 0
        True
    """
    dggs_min_res = int(DGGAL_TYPES[dggs_type]["min_res"])
    dggs_max_res = int(DGGAL_TYPES[dggs_type]["max_res"])

    # Process topology for points and multipoints if enabled
    if topology:
        estimated_resolution = dggs_max_res
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
            shortest_distance = shortest_point_distance(all_points)

            if shortest_distance > 0:
                for res in range(dggs_min_res, dggs_max_res + 1):
                    cell_diameter = _dggal_cell_diameter(dggs_type, res)
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        resolution = estimated_resolution

    resolution = validate_dggal_resolution(dggs_type, resolution)

    dggal_rows = []
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
            dggal_rows.extend(
                point2dggal(
                    dggs_type,
                    geom,
                    resolution,
                    props,
                    predicate,
                    compact,
                    topology=topology,  # Topology already processed above
                    include_properties=include_properties,
                    split_antimeridian=split_antimeridian,
                )
            )

        elif geom.geom_type in ("LineString", "MultiLineString"):
            dggal_rows.extend(
                polyline2dggal(
                    dggs_type,
                    geom,
                    resolution,
                    props,
                    predicate,
                    compact,
                    include_properties,
                    split_antimeridian=split_antimeridian,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            dggal_rows.extend(
                polygon2dggal(
                    dggs_type,
                    geom,
                    resolution,
                    props,
                    predicate,
                    compact,
                    include_properties,
                    split_antimeridian=split_antimeridian,
                )
            )
    return gpd.GeoDataFrame(dggal_rows, geometry="geometry", crs="EPSG:4326")


def vector2dggal(
    dggs_type: str,
    vector_data,
    resolution=None,
    predicate: str | None = None,
    compact: bool = False,
    topology: bool = False,
    include_properties: bool = True,
    output_format: str = "gpd",
    split_antimeridian: bool = False,
    verbose=True,
    **kwargs,
):
    """
    Convert vector data to DGGAL grid cells for a given type and resolution.

    Args:
        dggs_type (str): One of DGGAL_TYPES
        vector_data (str, geopandas.GeoDataFrame, pandas.DataFrame, dict, or list): Input vector data
        resolution (int, optional): Integer resolution. Required when topology=False,
            auto-calculated when topology=True.
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable DGGAL compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode to ensure disjoint points
            have disjoint DGGAL cells
        include_properties (bool, optional): Whether to include properties in output
        output_format (str, optional): Output format ('gpd', 'geojson', 'csv', 'shapefile', 'gpkg', 'parquet', 'geoparquet')
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        geopandas.GeoDataFrame, dict, or str: Output in the specified format. If output_format is a file-based format,
        the output will be saved to a file in the current directory with a default name based on the input.
        Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Example:
        >>> result = vector2dggal("isea3h", "data/points.geojson", resolution=10, output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    if resolution is not None:
        resolution = validate_dggal_resolution(dggs_type, resolution)

    gdf_input = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2dggal(
        dggs_type,
        gdf_input,
        resolution,
        predicate,
        compact,
        topology,
        include_properties,
        split_antimeridian=split_antimeridian,
        verbose=verbose,
    )

    # Return or export
    output_name = None
    if output_format in OUTPUT_FORMATS:
        res_label = resolution
        if res_label is None and not result.empty and "resolution" in result.columns:
            res_label = result["resolution"].iloc[0]
        output_name = f"{dggs_type}_grid_{res_label}"
    return convert_to_output_format(result, output_format, output_name)


def vector2dggal_cli():
    parser = argparse.ArgumentParser(
        description="Convert vector data to DGGAL grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path or URL")
    parser.add_argument(
        "-dggs",
        dest="dggs_type",
        type=str,
        required=True,
        choices=DGGAL_TYPES.keys(),
        help="DGGAL DGGS type",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        help="Resolution (integer). Required when topology=False, auto-calculated when topology=True",
    )
    parser.add_argument(
        "-p",
        "--predicate",
        choices=["intersect", "within", "centroid_within", "largest_overlap"],
        help="Spatial predicate for polygon conversion",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Use compact grid generation",
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
    add_verbose_argument(parser)
    args = parser.parse_args()

    try:
        result = vector2dggal(
            dggs_type=args.dggs_type,
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            include_properties=args.include_properties,
            output_format=args.output_format,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2dggal_cli()
