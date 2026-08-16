"""
Vector to GEOREF Module

Convert vector geometries to World Geographic Reference System (GEOREF) cells,
mirroring :mod:`vector2geohash` but using :func:`~vgrid.conversion.latlon2dggs.latlon2georef`,
:func:`~vgrid.conversion.dggs2geo.georef2geo.georef2geo`, and the same bbox lattice as
:func:`~vgrid.generator.georefgrid.georef_grid` (``numpy.arange`` over bounds at
:data:`~vgrid.utils.constants.GEOREF_RESOLUTION_DEGREES` for the target resolution).
"""

import argparse
import os
import sys
from math import sqrt

import geopandas as gpd
import numpy as np
from shapely.geometry import MultiPoint
from tqdm import tqdm

from vgrid.conversion.dggs2geo.georef2geo import georef2geo
from vgrid.conversion.latlon2dggs import latlon2georef
from vgrid.stats.georefstats import georef_metrics
from vgrid.utils.constants import (
    DGGS_TYPES,
    GEOREF_RESOLUTION_DEGREES,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
)
from vgrid.utils.geometry import (
    check_predicate,
    graticule_dggs_to_geoseries,
    shortest_point_distance,
)
from vgrid.utils.io import (
    convert_to_output_format,
    process_input_data_vector,
    validate_georef_resolution,
)

min_res = DGGS_TYPES["georef"]["min_res"]
max_res = DGGS_TYPES["georef"]["max_res"]


def point2georef(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Convert point or multipoint geometries to GEOREF cells at ``resolution``."""
    rows = []
    if feature.geom_type == "Point":
        points = [feature]
    elif feature.geom_type == "MultiPoint":
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        georef_id = latlon2georef(point.y, point.x, resolution)
        cell_polygon = georef2geo(georef_id)
        row = graticule_dggs_to_geoseries("georef", georef_id, resolution, cell_polygon)
        if include_properties and feature_properties:
            row.update(feature_properties)
        rows.append(row)
    return rows


def polyline2georef(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Collect GEOREF cells at ``resolution`` that intersect the line geometry."""
    rows = []
    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []

    for polyline in polylines:
        min_lon, min_lat, max_lon, max_lat = polyline.bounds
        resolution_degrees = GEOREF_RESOLUTION_DEGREES.get(resolution)
        longitudes = np.arange(min_lon, max_lon, resolution_degrees)
        latitudes = np.arange(min_lat, max_lat, resolution_degrees)
        for lon in longitudes:
            for lat in latitudes:
                georef_id = latlon2georef(lat, lon, resolution)
                cell_polygon = georef2geo(georef_id)
                if cell_polygon is not None and cell_polygon.intersects(polyline):
                    row = graticule_dggs_to_geoseries(
                        "georef", georef_id, resolution, cell_polygon
                    )
                    if include_properties and feature_properties:
                        row.update(feature_properties)
                    rows.append(row)
    return rows


def polygon2georef(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Collect GEOREF cells at ``resolution`` using ``predicate`` against the polygon."""
    rows = []
    if feature.geom_type == "Polygon":
        polygons = [feature]
    elif feature.geom_type == "MultiPolygon":
        polygons = list(feature.geoms)
    else:
        return []

    for polygon in polygons:
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        resolution_degrees = GEOREF_RESOLUTION_DEGREES.get(resolution)
        longitudes = np.arange(min_lon, max_lon, resolution_degrees)
        latitudes = np.arange(min_lat, max_lat, resolution_degrees)
        for lon in longitudes:
            for lat in latitudes:
                georef_id = latlon2georef(lat, lon, resolution)
                cell_polygon = georef2geo(georef_id)
                if cell_polygon is not None and check_predicate(
                    cell_polygon, polygon, predicate
                ):
                    row = graticule_dggs_to_geoseries(
                        "georef", georef_id, resolution, cell_polygon
                    )
                    if include_properties and feature_properties:
                        row.update(feature_properties)
                    rows.append(row)
    return rows


def geodataframe2georef(
    gdf,
    resolution=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Convert a GeoDataFrame to GEOREF cells."""
    if topology:
        estimated_resolution = max_res
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
                for res in range(min_res, max_res + 1):
                    _, avg_edge_len, _, _ = georef_metrics(res, unit="m")
                    cell_diameter_m = avg_edge_len * sqrt(2)
                    if cell_diameter_m < shortest_distance:
                        estimated_resolution = res
                        break
        resolution = estimated_resolution

    resolution = validate_georef_resolution(resolution)

    geom_col = gdf.geometry.name
    georef_rows = []
    for _, row in tqdm(gdf.iterrows(), desc="Processing features", total=len(gdf)):
        geom = row.geometry
        if geom is None:
            continue

        props = row.to_dict()
        if geom_col in props:
            del props[geom_col]

        if not include_properties:
            props = {}

        if geom.geom_type in ("Point", "MultiPoint"):
            georef_rows.extend(
                point2georef(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("LineString", "MultiLineString"):
            georef_rows.extend(
                polyline2georef(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            georef_rows.extend(
                polygon2georef(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
    if not georef_rows:
        return gpd.GeoDataFrame(
            {"geometry": gpd.GeoSeries([], crs="EPSG:4326")},
            geometry="geometry",
            crs="EPSG:4326",
        )
    return gpd.GeoDataFrame(georef_rows, geometry="geometry", crs="EPSG:4326")


def vector2georef(
    vector_data,
    resolution=None,
    predicate=None,
    topology=False,
    output_format="gpd",
    include_properties=True,
    **kwargs,
):
    """
    Convert vector data to GEOREF grid cells.

    Parameters mirror :func:`vector2geohash` except GEOREF has no compact mode.
    ``resolution`` is in ``[0..10]`` (see :data:`~vgrid.utils.constants.DGGS_TYPES`).
    """
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    if resolution is not None:
        resolution = validate_georef_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2georef(
        gdf, resolution, predicate, topology, include_properties
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2georef"
        else:
            output_name = "georef"
    return convert_to_output_format(result, output_format, output_name)


def vector2georef_cli():
    """CLI for :func:`vector2georef`."""
    parser = argparse.ArgumentParser(
        description="Convert vector data to GEOREF grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"GEOREF resolution [{min_res}..{max_res}]. Required when topology=False.",
    )
    parser.add_argument(
        "-p",
        "--predicate",
        choices=["intersect", "within", "centroid_within", "largest_overlap"],
        help="Spatial predicate for polygons",
    )
    parser.add_argument(
        "-t",
        "--topology",
        action="store_true",
        help="Topology-preserving resolution for disjoint points",
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
        result = vector2georef(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            topology=args.topology,
            output_format=args.output_format,
            include_properties=args.include_properties,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2georef_cli()
