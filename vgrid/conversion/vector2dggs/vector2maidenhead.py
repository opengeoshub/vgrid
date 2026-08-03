"""
Vector to Maidenhead Module

Convert vector geometries to Maidenhead locator cells using
:func:`~vgrid.conversion.latlon2dggs.latlon2maidenhead`,
:func:`~vgrid.conversion.dggs2geo.maidenhead2geo.maidenhead2geo`, and for lines and polygons
the same bbox-relative grid walk as :func:`~vgrid.generator.maidenheadgrid.maidenhead_grid_within_bbox`
(indexed from ``(-90, -180)`` with resolution-dependent cell sizes).
"""

import argparse
import math
import os
import sys

import geopandas as gpd
from shapely.geometry import MultiPoint
from tqdm import tqdm

from vgrid.conversion.dggs2geo.maidenhead2geo import maidenhead2geo
from vgrid.conversion.latlon2dggs import latlon2maidenhead
from vgrid.stats.maidenheadstats import maidenhead_metrics
from vgrid.dggs import maidenhead
from vgrid.utils.constants import (
    DGGS_TYPES,
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
    validate_maidenhead_resolution,
)

min_res = DGGS_TYPES["maidenhead"]["min_res"]
max_res = DGGS_TYPES["maidenhead"]["max_res"]


def _maidenhead_lon_lat_width(resolution):
    """Cell width in degrees (lon, lat) — matches :mod:`vgrid.generator.maidenheadgrid`."""
    if resolution == 1:
        return 20.0, 10.0
    if resolution == 2:
        return 2.0, 1.0
    if resolution == 3:
        return 0.083333, 0.041666
    if resolution == 4:
        return 0.008333, 0.004167
    raise ValueError("Unsupported resolution")


def _maidenhead_cell_records_for_bbox(resolution, bbox):
    """
    Build graticule records for all Maidenhead cells whose axis-aligned bounds overlap
    ``bbox`` (``[min_lon, min_lat, max_lon, max_lat]``), same indexing as
    ``maidenhead_grid_within_bbox`` (no tqdm).
    """
    resolution = validate_maidenhead_resolution(resolution)
    lon_width, lat_width = _maidenhead_lon_lat_width(resolution)
    min_lon, min_lat, max_lon, max_lat = bbox
    base_lat, base_lon = -90, -180
    start_x = math.floor((min_lon - base_lon) / lon_width)
    end_x = math.floor((max_lon - base_lon) / lon_width)
    start_y = math.floor((min_lat - base_lat) / lat_width)
    end_y = math.floor((max_lat - base_lat) / lat_width)

    records = []
    for x in range(start_x, end_x + 1):
        for y in range(start_y, end_y + 1):
            cell_min_lon = base_lon + x * lon_width
            cell_max_lon = cell_min_lon + lon_width
            cell_min_lat = base_lat + y * lat_width
            cell_max_lat = cell_min_lat + lat_width

            if (
                cell_max_lon < min_lon
                or cell_min_lon > max_lon
                or cell_max_lat < min_lat
                or cell_min_lat > max_lat
            ):
                continue

            cell_center_lat = (cell_min_lat + cell_max_lat) / 2
            cell_center_lon = (cell_min_lon + cell_max_lon) / 2
            maidenhead_id = maidenhead.toMaiden(
                cell_center_lat, cell_center_lon, resolution
            )
            cell_polygon = maidenhead2geo(maidenhead_id)
            records.append(
                graticule_dggs_to_geoseries(
                    "maidenhead", maidenhead_id, resolution, cell_polygon
                )
            )
    return records


def point2maidenhead(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Convert point or multipoint geometries to Maidenhead cells at ``resolution``."""
    rows = []
    if feature.geom_type == "Point":
        points = [feature]
    elif feature.geom_type == "MultiPoint":
        points = list(feature.geoms)
    else:
        return []

    for point in points:
        maidenhead_id = latlon2maidenhead(point.y, point.x, resolution)
        cell_polygon = maidenhead2geo(maidenhead_id)
        row = graticule_dggs_to_geoseries(
            "maidenhead", maidenhead_id, resolution, cell_polygon
        )
        if include_properties and feature_properties:
            row.update(feature_properties)
        rows.append(row)
    return rows


def polyline2maidenhead(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Collect cells from the bbox Maidenhead lattice that intersect the line."""
    rows = []
    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []

    seen = set()
    for polyline in polylines:
        bbox = list(polyline.bounds)
        for record in _maidenhead_cell_records_for_bbox(resolution, bbox):
            mid = record["maidenhead"]
            if mid in seen:
                continue
            cell_geom = record["geometry"]
            if not cell_geom.intersects(polyline):
                continue
            seen.add(mid)
            row = dict(record)
            if include_properties and feature_properties:
                row = {**row, **feature_properties}
            rows.append(row)
    return rows


def polygon2maidenhead(
    feature,
    resolution,
    feature_properties=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Collect cells from the bbox Maidenhead lattice, filtered by ``predicate``."""
    rows = []
    if feature.geom_type == "Polygon":
        polygons = [feature]
    elif feature.geom_type == "MultiPolygon":
        polygons = list(feature.geoms)
    else:
        return []

    seen = set()
    for polygon in polygons:
        bbox = list(polygon.bounds)
        for record in _maidenhead_cell_records_for_bbox(resolution, bbox):
            mid = record["maidenhead"]
            if mid in seen:
                continue
            cell_geom = record["geometry"]
            if not check_predicate(cell_geom, polygon, predicate):
                continue
            seen.add(mid)
            row = dict(record)
            if include_properties and feature_properties:
                row = {**row, **feature_properties}
            rows.append(row)
    return rows


def geodataframe2maidenhead(
    gdf,
    resolution=None,
    predicate=None,
    topology=False,
    include_properties=True,
):
    """Convert a GeoDataFrame to Maidenhead cells."""
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
                    _, avg_edge_len, _, _ = maidenhead_metrics(res, unit="m")
                    cell_diameter_m = avg_edge_len * math.sqrt(2)
                    if cell_diameter_m < shortest_distance:
                        estimated_resolution = res
                        break
        resolution = estimated_resolution

    resolution = validate_maidenhead_resolution(resolution)

    geom_col = gdf.geometry.name
    maidenhead_rows = []
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
            maidenhead_rows.extend(
                point2maidenhead(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("LineString", "MultiLineString"):
            maidenhead_rows.extend(
                polyline2maidenhead(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            maidenhead_rows.extend(
                polygon2maidenhead(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    topology=topology,
                    include_properties=include_properties,
                )
            )
    if not maidenhead_rows:
        return gpd.GeoDataFrame(
            {"geometry": gpd.GeoSeries([], crs="EPSG:4326")},
            geometry="geometry",
            crs="EPSG:4326",
        )
    return gpd.GeoDataFrame(maidenhead_rows, geometry="geometry", crs="EPSG:4326")


def vector2maidenhead(
    vector_data,
    resolution=None,
    predicate=None,
    topology=False,
    output_format="gpd",
    include_properties=True,
    **kwargs,
):
    """
    Convert vector data to Maidenhead grid cells.

    ``resolution`` is ``1..4``. Lines and polygons walk the same bbox-relative
    Maidenhead lattice as ``maidenhead_grid_within_bbox`` (implemented locally).
    """
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    if resolution is not None:
        resolution = validate_maidenhead_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2maidenhead(
        gdf, resolution, predicate, topology, include_properties
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2maidenhead"
        else:
            output_name = "maidenhead"
    return convert_to_output_format(result, output_format, output_name)


def vector2maidenhead_cli():
    """CLI for :func:`vector2maidenhead`."""
    parser = argparse.ArgumentParser(
        description="Convert vector data to Maidenhead grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=f"Maidenhead resolution [{min_res}..{max_res}]. Required when topology=False.",
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
        result = vector2maidenhead(
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
    vector2maidenhead_cli()
