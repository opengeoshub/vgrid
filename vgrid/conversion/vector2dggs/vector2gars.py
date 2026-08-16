"""
Vector to GARS Module

Convert vector geometries to GARS (Global Area Reference System) cells using
:func:`~vgrid.conversion.latlon2dggs.latlon2gars`, :func:`~vgrid.conversion.dggs2geo.gars2geo.gars2geo`,
and the same bbox lattice as :func:`~vgrid.generator.garsgrid.gars_grid`.
"""

import argparse
import os
import sys
from math import sqrt

import geopandas as gpd
import numpy as np
from shapely.geometry import MultiPoint
from tqdm import tqdm

from vgrid.conversion.dggs2geo.gars2geo import gars2geo
from vgrid.conversion.latlon2dggs import latlon2gars
from vgrid.stats.garsstats import gars_metrics
from vgrid.utils.constants import (
    DGGS_TYPES,
    GARS_RESOLUTION_MINUTES,
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
    validate_gars_resolution,
)

min_res = DGGS_TYPES["gars"]["min_res"]
max_res = DGGS_TYPES["gars"]["max_res"]


def point2gars(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Convert point or multipoint geometries to GARS cells at ``resolution``."""
    rows = []
    if feature.geom_type == "Point":
        points = [feature]
    elif feature.geom_type == "MultiPoint":
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        gars_id = latlon2gars(point.y, point.x, resolution)
        cell_polygon = gars2geo(gars_id)
        row = graticule_dggs_to_geoseries("gars", gars_id, resolution, cell_polygon)
        if include_properties and feature_properties:
            row.update(feature_properties)
        rows.append(row)
    return rows


def polyline2gars(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Collect GARS cells at ``resolution`` that intersect the line geometry."""
    rows = []
    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []

    resolution_minutes = GARS_RESOLUTION_MINUTES.get(resolution)
    resolution_degrees = resolution_minutes / 60.0

    for polyline in polylines:
        min_lon, min_lat, max_lon, max_lat = polyline.bounds
        longitudes = np.arange(min_lon, max_lon, resolution_degrees)
        latitudes = np.arange(min_lat, max_lat, resolution_degrees)
        for lon in longitudes:
            for lat in latitudes:
                gars_id = latlon2gars(lat, lon, resolution)
                cell_polygon = gars2geo(gars_id)
                if not cell_polygon.intersects(polyline):
                    continue
                row = graticule_dggs_to_geoseries(
                    "gars", gars_id, resolution, cell_polygon
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                rows.append(row)
    return rows


def polygon2gars(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Collect GARS cells at ``resolution`` using ``predicate`` against the polygon."""
    rows = []
    if feature.geom_type == "Polygon":
        polygons = [feature]
    elif feature.geom_type == "MultiPolygon":
        polygons = list(feature.geoms)
    else:
        return []

    resolution_minutes = GARS_RESOLUTION_MINUTES.get(resolution)
    resolution_degrees = resolution_minutes / 60.0

    for polygon in polygons:
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        longitudes = np.arange(min_lon, max_lon, resolution_degrees)
        latitudes = np.arange(min_lat, max_lat, resolution_degrees)
        seen = set()
        for lon in longitudes:
            for lat in latitudes:
                gars_id = latlon2gars(lat, lon, resolution)
                cell_polygon = gars2geo(gars_id)
                if not check_predicate(cell_polygon, polygon, predicate):
                    continue
                seen.add(gars_id)
                row = graticule_dggs_to_geoseries(
                    "gars", gars_id, resolution, cell_polygon
                )
                if include_properties and feature_properties:
                    row.update(feature_properties)
                rows.append(row)
    return rows


def geodataframe2gars(
    gdf,
    resolution=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Convert a GeoDataFrame to GARS cells."""
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
                    _, avg_edge_len, _, _ = gars_metrics(res, unit="m")
                    cell_diameter_m = avg_edge_len * sqrt(2)
                    if cell_diameter_m < shortest_distance:
                        estimated_resolution = res
                        break
        resolution = estimated_resolution

    resolution = validate_gars_resolution(resolution)

    geom_col = gdf.geometry.name
    gars_rows = []
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
            gars_rows.extend(
                point2gars(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("LineString", "MultiLineString"):
            gars_rows.extend(
                polyline2gars(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            gars_rows.extend(
                polygon2gars(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
    if not gars_rows:
        return gpd.GeoDataFrame(
            {"geometry": gpd.GeoSeries([], crs="EPSG:4326")},
            geometry="geometry",
            crs="EPSG:4326",
        )
    return gpd.GeoDataFrame(gars_rows, geometry="geometry", crs="EPSG:4326")


def vector2gars(
    vector_data,
    resolution=None,
    predicate=None,
    topology=False,
    output_format="gpd",
    include_properties=True,
    **kwargs,
):
    """
    Convert vector data to GARS grid cells.

    ``resolution`` is ``1..4`` (30′, 15′, 5′, 1′ cells). There is no compact mode.
    """
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    if resolution is not None:
        resolution = validate_gars_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2gars(gdf, resolution, predicate, topology, include_properties)
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2gars"
        else:
            output_name = "gars"
    return convert_to_output_format(result, output_format, output_name)


def vector2gars_cli():
    """CLI for :func:`vector2gars`."""
    parser = argparse.ArgumentParser(
        description="Convert vector data to GARS grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=(
            f"GARS resolution [{min_res}..{max_res}] "
            "(1=30min, 2=15min, 3=5min, 4=1min). Required when topology=False."
        ),
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
        result = vector2gars(
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
    vector2gars_cli()
