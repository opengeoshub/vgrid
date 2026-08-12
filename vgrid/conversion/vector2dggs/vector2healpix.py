"""
Vector to HEALPix Module

This module provides functionality to convert vector geometries to HEALPix grid
cells (UNIQ scheme) with flexible input and output formats.

Key Functions:
    point2healpix: Convert point geometries to HEALPix cells
    polyline2healpix: Convert line geometries to HEALPix cells
    polygon2healpix: Convert polygon geometries to HEALPix cells with spatial predicates
    geodataframe2healpix: Convert GeoDataFrame to HEALPix cells with topology support
    vector2healpix: Main function for converting various input formats to HEALPix cells
    vector2healpix_cli: Command-line interface for batch processing
"""

import sys
import os
import argparse

import geopandas as gpd
from shapely.geometry import MultiPoint
from tqdm import tqdm

from vgrid.conversion.latlon2dggs import latlon2healpix
from vgrid.conversion.dggs2geo.healpix2geo import healpix2geo
from vgrid.conversion.dggscompact.healpixcompact import healpix_compact
from vgrid.dggs.healpix import (
    nside2resol,
    order2nside,
    orderpix2uniq,
    queryBoxInclusiveNest,
    uniq2orderpix,
)
from vgrid.utils.geometry import (
    check_predicate,
    geodesic_dggs_to_geoseries,
    shortest_point_distance,
)
from vgrid.utils.io import (
    convert_to_output_format,
    process_input_data_vector,
    validate_healpix_resolution,
)
from vgrid.utils.constants import DGGS_TYPES, OUTPUT_FORMATS, STRUCTURED_FORMATS

min_res = DGGS_TYPES["healpix"]["min_res"]
max_res = DGGS_TYPES["healpix"]["max_res"]
EARTH_RADIUS_M = 6371008.8


def _empty_healpix_gdf():
    return gpd.GeoDataFrame(
        columns=[
            "healpix",
            "resolution",
            "center_lat",
            "center_lon",
            "avg_edge_len",
            "cell_area",
            "cell_perimeter",
            "geometry",
        ],
        geometry="geometry",
        crs="EPSG:4326",
    )


def _healpix_avg_edge_m(resolution: int) -> float:
    """Approximate average HEALPix cell edge length in meters."""
    return EARTH_RADIUS_M * nside2resol(order2nside(resolution))


def _nested_to_uniq(resolution: int, nested_ids):
    """Convert nested pixel indices at a fixed resolution to UNIQ IDs."""
    return [orderpix2uniq(resolution, int(ipix)) for ipix in nested_ids]


def _bbox_cells(geometry, resolution: int):
    """Return UNIQ cell IDs covering a geometry's bounding box."""
    min_lon, min_lat, max_lon, max_lat = geometry.bounds
    nside = order2nside(resolution)
    nested_ids = queryBoxInclusiveNest(
        nside, (min_lon, min_lat, max_lon, max_lat)
    )
    return _nested_to_uniq(resolution, nested_ids)


def _rows_from_uniq_ids(
    uniq_ids,
    feature_properties=None,
    include_properties=True,
    fix_antimeridian=None,
):
    """Build geoseries rows from HEALPix UNIQ IDs."""
    rows = []
    for uniq_id in uniq_ids:
        cell_polygon = healpix2geo(uniq_id, fix_antimeridian=fix_antimeridian)
        if cell_polygon is None or getattr(cell_polygon, "is_empty", False):
            continue
        order, _ipix = uniq2orderpix(int(uniq_id))
        row = geodesic_dggs_to_geoseries(
            "healpix", uniq_id, order, cell_polygon, num_edges=4
        )
        if include_properties and feature_properties:
            row.update(feature_properties)
        rows.append(row)
    return rows


def point2healpix(
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
    Convert a point geometry to HEALPix grid cells.

    Parameters
    ----------
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert.
    resolution : int
        HEALPix resolution/order [0..29].
    feature_properties : dict, optional
        Properties to include in output features.
    predicate : str, optional
        Not used for points (API compatibility).
    compact : bool, optional
        Not used for points (API compatibility).
    topology : bool, optional
        Handled by geodataframe2healpix.
    include_properties : bool, optional
        Whether to include properties in output.
    fix_antimeridian : str, optional
        Antimeridian fixing method.

    Returns
    -------
    list of dict
        HEALPix cell rows containing the point(s).
    """
    if feature.geom_type == "Point":
        points = [feature]
    elif feature.geom_type == "MultiPoint":
        points = list(feature.geoms)
    else:
        return []

    healpix_rows = []
    for point in points:
        healpix_id = latlon2healpix(point.y, point.x, resolution)
        cell_polygon = healpix2geo(healpix_id, fix_antimeridian=fix_antimeridian)
        order, _ipix = uniq2orderpix(healpix_id)
        row = geodesic_dggs_to_geoseries(
            "healpix", healpix_id, order, cell_polygon, num_edges=4
        )
        if include_properties and feature_properties:
            row.update(feature_properties)
        healpix_rows.append(row)
    return healpix_rows


def polyline2healpix(
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
    Convert a polyline geometry to HEALPix grid cells.

    Uses ``queryBoxInclusiveNest`` on the line bbox, then keeps cells that
    intersect the polyline.
    """
    if feature.geom_type == "LineString":
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        polylines = list(feature.geoms)
    else:
        return []

    healpix_rows = []
    for polyline in polylines:
        if polyline is None or polyline.is_empty:
            continue
        candidate_ids = _bbox_cells(polyline, resolution)
        intersecting_ids = []
        for uniq_id in candidate_ids:
            cell_polygon = healpix2geo(uniq_id, fix_antimeridian=fix_antimeridian)
            if cell_polygon is None or cell_polygon.is_empty:
                continue
            if cell_polygon.intersects(polyline):
                intersecting_ids.append(uniq_id)

        healpix_rows.extend(
            _rows_from_uniq_ids(
                intersecting_ids,
                feature_properties=feature_properties,
                include_properties=include_properties,
                fix_antimeridian=fix_antimeridian,
            )
        )
    return healpix_rows


def polygon2healpix(
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
    Convert a polygon geometry to HEALPix grid cells.

    Uses ``queryBoxInclusiveNest`` on the polygon bbox as candidates, then
    filters with ``check_predicate``. Optional UNIQ compaction is applied after
    filtering.
    """
    if feature.geom_type == "Polygon":
        polygons = [feature]
    elif feature.geom_type == "MultiPolygon":
        polygons = list(feature.geoms)
    else:
        return []

    healpix_rows = []
    for polygon in polygons:
        if polygon is None or polygon.is_empty:
            continue

        candidate_ids = _bbox_cells(polygon, resolution)
        if not candidate_ids:
            continue

        filtered_ids = []
        for uniq_id in tqdm(
            candidate_ids, desc="Generating HEALPix cells", unit=" cells"
        ):
            cell_polygon = healpix2geo(uniq_id, fix_antimeridian=fix_antimeridian)
            if cell_polygon is None or cell_polygon.is_empty:
                continue
            if not check_predicate(cell_polygon, polygon, predicate):
                continue
            filtered_ids.append(uniq_id)

        if compact and filtered_ids:
            filtered_ids = healpix_compact(filtered_ids)

        healpix_rows.extend(
            _rows_from_uniq_ids(
                filtered_ids,
                feature_properties=feature_properties,
                include_properties=include_properties,
                fix_antimeridian=fix_antimeridian,
            )
        )
    return healpix_rows


def geodataframe2healpix(
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    fix_antimeridian=None,
):
    """
    Convert a GeoDataFrame to HEALPix grid cells.

    When ``topology=True``, resolution is estimated so disjoint points receive
    disjoint HEALPix cells.
    """
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
                    cell_diameter = _healpix_avg_edge_m(res) * 2
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break
        resolution = estimated_resolution

    resolution = validate_healpix_resolution(resolution)

    healpix_rows = []
    for _, row in tqdm(gdf.iterrows(), desc="Processing features", total=len(gdf)):
        geom = row.geometry
        if geom is None:
            continue

        props = row.to_dict()
        if "geometry" in props:
            del props["geometry"]
        if not include_properties:
            props = {}

        if geom.geom_type in ("Point", "MultiPoint"):
            healpix_rows.extend(
                point2healpix(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    topology=topology,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )
        elif geom.geom_type in ("LineString", "MultiLineString"):
            healpix_rows.extend(
                polyline2healpix(
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
            healpix_rows.extend(
                polygon2healpix(
                    feature=geom,
                    resolution=resolution,
                    feature_properties=props,
                    predicate=predicate,
                    compact=compact,
                    include_properties=include_properties,
                    fix_antimeridian=fix_antimeridian,
                )
            )

    if not healpix_rows:
        return _empty_healpix_gdf()
    return gpd.GeoDataFrame(healpix_rows, geometry="geometry", crs="EPSG:4326")


def vector2healpix(
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
    Convert vector data to HEALPix grid cells from various input formats.

    Parameters
    ----------
    vector_data : str, GeoDataFrame, DataFrame, dict, or list
        Input vector data.
    resolution : int, optional
        HEALPix resolution/order [0..29]. Required when topology=False.
    predicate : str, optional
        Spatial predicate for polygons.
    compact : bool, optional
        Enable UNIQ compact mode for polygons.
    topology : bool, optional
        Auto-estimate resolution for disjoint points.
    output_format : str, optional
        Output format (gpd, geojson, csv, ...).
    include_properties : bool, optional
        Whether to include original feature properties.
    fix_antimeridian : str, optional
        Antimeridian fixing method.
    **kwargs
        Passed to ``process_input_data_vector``.
    """
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")

    if resolution is not None:
        resolution = validate_healpix_resolution(resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2healpix(
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
            output_name = f"{base}2healpix"
        else:
            output_name = "healpix"
    return convert_to_output_format(result, output_format, output_name)


def vector2healpix_cli():
    """Command-line interface for vector2healpix conversion."""
    parser = argparse.ArgumentParser(
        description="Convert vector data to HEALPix grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path, URL")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=range(min_res, max_res + 1),
        help=(
            f"HEALPix resolution [{min_res}..{max_res}]. "
            "Required when topology=False, auto-calculated when topology=True"
        ),
    )
    parser.add_argument(
        "-p",
        "--predicate",
        choices=["intersect", "within", "centroid_within", "largest_overlap"],
        help="Spatial predicate for polygons",
    )
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable HEALPix UNIQ compact mode for polygons",
    )
    parser.add_argument(
        "-t",
        "--topology",
        action="store_true",
        help="Enable topology preserving mode",
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
        help="Antimeridian fixing method",
    )
    args = parser.parse_args()

    try:
        result = vector2healpix(
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            output_format=args.output_format,
            include_properties=args.include_properties,
            fix_antimeridian=args.fix_antimeridian,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2healpix_cli()
