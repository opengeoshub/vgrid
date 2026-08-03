"""
Vector to DGGRID Module

This module provides functionality to convert vector geometries to DGGRID grid cells with flexible input and output formats.

Key Functions:
    point2dggrid: Convert point geometries to DGGRID cells
    polyline2dggrid: Convert line geometries to DGGRID cells
    polygon2dggrid: Convert polygon geometries to DGGRID cells with spatial predicates
    geodataframe2dggrid: Convert GeoDataFrame to DGGRID cells
    vector2dggrid: Main function for converting various input formats to DGGRID cells
    vector2dggrid_cli: Command-line interface for batch processing
"""

import argparse
import json
import math
import sys
import os
from tqdm import tqdm
from pyproj import Geod
import geopandas as gpd
import pandas as pd
from shapely.geometry import box
from vgrid.utils.io import (
    validate_dggrid_type,
    validate_dggrid_resolution,
    create_dggrid_instance,
)
from vgrid.conversion.latlon2dggs import latlon2dggrid
from vgrid.conversion.dggs2geo.dggrid2geo import dggrid2geo
from vgrid.utils.geometry import (
    _metric_crs,
    check_predicate,
    dggrid_num_edges,
)
from vgrid.stats.dggridstats import dggrid_metrics, dggridstats
from dggrid4py.dggrid_runner import output_address_types
from vgrid.utils.io import process_input_data_vector, convert_to_output_format
from vgrid.utils.constants import DGGRID_TYPES, OUTPUT_FORMATS, STRUCTURED_FORMATS

geod = Geod(ellps="WGS84")


def _dggrid_cell_diameter(
    dggrid_instance,
    dggs_type: str,
    resolution: int,
    grid_stats: pd.DataFrame | None = None,
) -> float:
    """Approximate cell diameter (m) from ``dggrid_metrics``."""
    _, avg_edge_len, _, _ = dggrid_metrics(
        dggrid_instance,
        dggs_type,
        resolution,
        unit="m",
        grid_stats=grid_stats,
    )
    if dggrid_num_edges(dggs_type) == 6:
        return avg_edge_len * 2
    return avg_edge_len * math.sqrt(2)


def _shortest_point_distance_sjoin(points_gdf: gpd.GeoDataFrame) -> float:
    """
    Minimum distance (meters) between distinct points via ``sjoin_nearest``.

    Each point is joined to its nearest *other* point (self-matches excluded).
    Geographic inputs are reprojected to a local metric CRS first.
    """
    if len(points_gdf) < 2:
        return 0.0

    metric_crs = _metric_crs(points_gdf)
    if (
        metric_crs is not None
        and points_gdf.crs is not None
        and points_gdf.crs.is_geographic
    ):
        pts = points_gdf.to_crs(metric_crs)
    else:
        pts = points_gdf

    shortest_distance = float("inf")
    for idx in pts.index:
        query = pts.loc[[idx]]
        others = pts.drop(index=idx)
        joined = gpd.sjoin_nearest(query, others, how="inner", distance_col="distance")
        if not joined.empty:
            shortest_distance = min(shortest_distance, joined["distance"].min())

    return 0.0 if shortest_distance == float("inf") else float(shortest_distance)


# Function to generate grid for Point
def point2dggrid(
    dggrid_instance,
    dggs_type,
    feature,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    feature_properties=None,
    output_address_type="SEQNUM",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    """
    Convert a point geometry to DGGRID grid cells.

    Converts point or multipoint geometries to DGGRID grid cells at the specified resolution.
    Each point is assigned to its containing DGGRID cell.

    Parameters
    ----------
    dggrid_instance : object
        DGGRID instance for grid operations.
    dggs_type : str
        DGGRID DGGS type (e.g., "isea4h", "fuller").
    feature : shapely.geometry.Point or shapely.geometry.MultiPoint
        Point geometry to convert to DGGRID cells.
    resolution : int
        DGGRID resolution level.
    predicate : str, optional
        Spatial predicate to apply (not used for points).
    compact : bool, optional
        Enable DGGRID compact mode (not used for points).
    topology : bool, optional
        Enable topology preserving mode (handled by ``geodataframe2dggrid``).
    include_properties : bool, optional
        Whether to include properties in output.
    feature_properties : dict, optional
        Properties to include in output features.
    output_address_type : str, optional
        Output address type (e.g., "SEQNUM", "Q2DI", "Q2DD"). Defaults to "SEQNUM".
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.

    Returns
    -------
    geopandas.GeoDataFrame
        GeoDataFrame containing DGGRID cells with the point(s).

    Examples
    --------
    >>> from shapely.geometry import Point
    >>> point = Point(-122.4194, 37.7749)  # San Francisco
    >>> gdf = point2dggrid(dggrid_instance, "isea4h", point, 10)
    >>> len(gdf)
    1

    >>> from shapely.geometry import MultiPoint
    >>> points = MultiPoint([(-122.4194, 37.7749), (-74.0060, 40.7128)])
    >>> gdf = point2dggrid(dggrid_instance, "fuller", points, 8)
    >>> len(gdf)
    2
    """
    dggs_type = validate_dggrid_type(dggs_type)
    resolution = validate_dggrid_resolution(dggs_type, resolution)

    if feature.geom_type == "Point":
        points = [feature]
    elif feature.geom_type == "MultiPoint":
        points = list(feature.geoms)
    else:
        return gpd.GeoDataFrame(columns=["geometry"], crs="EPSG:4326")

    seqnums = []
    for point in points:
        seqnum = latlon2dggrid(
            dggrid_instance,
            dggs_type,
            float(point.y),
            float(point.x),
            resolution,
            output_address_type,
        )
        seqnums.append(seqnum)

    # Build polygons from SEQNUM ids
    gdf = dggrid2geo(
        dggrid_instance,
        dggs_type,
        seqnums,
        resolution,
        output_address_type,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        options=options,
    )
    if include_properties and feature_properties:
        for key, value in feature_properties.items():
            gdf[key] = value
    return gdf


# Function to generate grid for Polyline
def polyline2dggrid(
    dggrid_instance,
    dggs_type,
    feature,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    feature_properties=None,
    output_address_type="SEQNUM ",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    """
    Generate DGGRID cells intersecting with a LineString or MultiLineString geometry.

    Args:
        dggrid_instance: DGGRIDv7 instance for grid operations.
        dggs_type (str): Type of DGGS (e.g., ISEA4H, FULLER, etc.).
        res (int): Resolution for the DGGRID.
        address_type (str): Address type for the output grid cells.
        geometry (shapely.geometry.LineString or MultiLineString): Input geometry.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
    Returns:
        geopandas.GeoDataFrame: GeoDataFrame containing DGGRID cells intersecting with the input geometry.
    """
    # Initialize an empty list to store filtered grid cells
    merged_grids = []

    # Check the geometry type
    if feature.geom_type == "LineString":
        # Handle single LineString
        polylines = [feature]
    elif feature.geom_type == "MultiLineString":
        # Handle MultiLineString: process each line separately
        polylines = list(feature.geoms)

    # Process each polyline
    for polyline in polylines:
        # Get bounding box for the current polyline
        bounding_box = box(*polyline.bounds)

        # Generate grid cells for the bounding box
        kwargs = {
            "split_dateline": split_antimeridian,
            "output_address_type": output_address_type,
        }
        if options:
            kwargs.update(options)
        dggrid_gdf = dggrid_instance.grid_cell_polygons_for_extent(
            dggs_type,
            resolution,
            clip_geom=bounding_box,
            **kwargs,
        )

        # Keep only grid cells that match predicate (defaults to intersects)
        dggrid_gdf = dggrid_gdf[dggrid_gdf.intersects(polyline)]

        try:
            if output_address_type != "SEQNUM":

                def address_transform(
                    dggrid_seqnum, dggal_type, resolution, address_type
                ):
                    address_type_transform = dggrid_instance.address_transform(
                        [dggrid_seqnum],
                        dggs_type=dggs_type,
                        resolution=resolution,
                        mixed_aperture_level=None,
                        input_address_type="SEQNUM",
                        output_address_type=output_address_type,
                    )
                    return address_type_transform.loc[0, address_type]

                dggrid_gdf["name"] = dggrid_gdf["name"].astype(str)
                dggrid_gdf["name"] = dggrid_gdf["name"].apply(
                    lambda val: address_transform(
                        val, dggs_type, resolution, output_address_type
                    )
                )
                dggrid_gdf = dggrid_gdf.rename(
                    columns={"name": output_address_type.lower()}
                )
            else:
                dggrid_gdf = dggrid_gdf.rename(columns={"name": "seqnum"})

        except Exception:
            pass
        # Append the filtered GeoDataFrame to the list
        if include_properties and feature_properties and not dggrid_gdf.empty:
            for key, value in feature_properties.items():
                dggrid_gdf[key] = value
        merged_grids.append(dggrid_gdf)

    # Merge all filtered grids into one GeoDataFrame
    if merged_grids:
        final_grid = gpd.GeoDataFrame(
            pd.concat(merged_grids, ignore_index=True), crs=merged_grids[0].crs
        )
    else:
        final_grid = gpd.GeoDataFrame(columns=["geometry"], crs="EPSG:4326")

    if split_antimeridian:
        if aggregate:
            final_grid = final_grid.dissolve(by="global_id", as_index=False)
    return final_grid


def polygon2dggrid(
    dggrid_instance,
    dggs_type,
    feature,
    resolution,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    feature_properties=None,
    output_address_type="SEQNUM",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    """
    Generate DGGRID cells intersecting with a given polygon or multipolygon geometry.

    Args:
        dggrid_instance: DGGRIDv7 instance for grid operations.
        dggs_type (str): Type of DGGS (e.g., ISEA4H, FULLER, etc.).
        res (int): Resolution for the DGGRID.
        address_type (str): Address type for the output grid cells.
        geometry (shapely.geometry.Polygon or MultiPolygon): Input geometry.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
    Returns:
        geopandas.GeoDataFrame: GeoDataFrame containing DGGRID cells intersecting with the input geometry.
    """
    # Initialize an empty list to store filtered grid cells
    merged_grids = []

    # Check the geometry type
    if feature.geom_type == "Polygon":
        # Handle single Polygon
        polygons = [feature]
    elif feature.geom_type == "MultiPolygon":
        # Handle MultiPolygon: process each polygon separately
        polygons = list(feature.geoms)  # Use .geoms to get components of MultiPolygon

    # Process each polygon
    for polygon in polygons:
        # Get bounding box for the current polygon
        bounding_box = box(*polygon.bounds)

        # Generate grid cells for the bounding box
        kwargs = {
            "split_dateline": split_antimeridian,
            "output_address_type": output_address_type,
        }
        if options:
            kwargs.update(options)
        dggrid_gdf = dggrid_instance.grid_cell_polygons_for_extent(
            dggs_type,
            resolution,
            clip_geom=bounding_box,
            **kwargs,
        )

        # Keep only grid cells that satisfy predicate (defaults to intersects)
        if predicate:
            dggrid_gdf = dggrid_gdf[
                dggrid_gdf.geometry.apply(
                    lambda cell: check_predicate(cell, polygon, predicate)
                )
            ]
        else:
            dggrid_gdf = dggrid_gdf[dggrid_gdf.intersects(polygon)]
        try:
            if output_address_type != "SEQNUM":

                def address_transform(
                    dggrid_seqnum, dggs_type, resolution, address_type
                ):
                    address_type_transform = dggrid_instance.address_transform(
                        [dggrid_seqnum],
                        dggs_type=dggs_type,
                        resolution=resolution,
                        mixed_aperture_level=None,
                        input_address_type="SEQNUM",
                        output_address_type=output_address_type,
                    )
                    return address_type_transform.loc[0, address_type]

                dggrid_gdf["name"] = dggrid_gdf["name"].astype(str)
                dggrid_gdf["name"] = dggrid_gdf["name"].apply(
                    lambda val: address_transform(
                        val, dggs_type, resolution, output_address_type
                    )
                )
                dggrid_gdf = dggrid_gdf.rename(
                    columns={"name": output_address_type.lower()}
                )
            else:
                dggrid_gdf = dggrid_gdf.rename(columns={"name": "seqnum"})

        except Exception:
            pass

        # Append the filtered GeoDataFrame to the list
        if include_properties and feature_properties and not dggrid_gdf.empty:
            for key, value in feature_properties.items():
                dggrid_gdf[key] = value
        merged_grids.append(dggrid_gdf)

    # Merge all filtered grids into one GeoDataFrame
    if merged_grids:
        final_grid = gpd.GeoDataFrame(
            pd.concat(merged_grids, ignore_index=True), crs=merged_grids[0].crs
        )
    else:
        final_grid = gpd.GeoDataFrame(columns=["geometry"], crs="EPSG:4326")
    if split_antimeridian:
        if aggregate:
            final_grid = final_grid.dissolve(by="global_id", as_index=False)
    return final_grid


def geodataframe2dggrid(
    dggrid_instance,
    dggs_type,
    gdf,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    output_address_type="SEQNUM",
    split_antimeridian=False,
    aggregate=False,
    options=None,
):
    """
    Convert a GeoDataFrame to DGGRID grid cells.

    Args:
        dggrid_instance: DGGRIDv7 instance for grid operations
        gdf (geopandas.GeoDataFrame): GeoDataFrame to convert
        dggs_type (str): One of DGGRID_TYPES
        resolution (int, optional): Integer resolution. Required when topology=False,
            auto-calculated when topology=True.
        predicate (str, optional): Spatial predicate to apply for polygons ('intersect', 'within', 'centroid_within', 'largest_overlap')
        compact (bool, optional): Enable compact mode for polygons and lines
        topology (bool, optional): Enable topology preserving mode to ensure disjoint
            points have disjoint DGGRID cells
        include_properties (bool, optional): Whether to include properties in output
        output_address_type (str, optional): Output address type (SEQNUM, Q2DI, Q2DD, etc.)
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.

    Returns:
        geopandas.GeoDataFrame: GeoDataFrame with DGGRID grid cells
    """
    dggs_type = validate_dggrid_type(dggs_type)

    if topology:
        min_res = int(DGGRID_TYPES[dggs_type]["min_res"])
        max_res = int(DGGRID_TYPES[dggs_type]["max_res"])
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

        if len(points_list) >= 2:
            points_gdf = gpd.GeoDataFrame(
                geometry=points_list, crs=gdf.crs or "EPSG:4326"
            )
            shortest_distance = _shortest_point_distance_sjoin(points_gdf)

            if shortest_distance > 0:
                grid_stats = dggridstats(dggrid_instance, dggs_type, unit="m")
                for res in range(min_res, max_res + 1):
                    try:
                        cell_diameter = _dggrid_cell_diameter(
                            dggrid_instance,
                            dggs_type,
                            res,
                            grid_stats=grid_stats,
                        )
                    except ValueError:
                        continue
                    if cell_diameter < shortest_distance:
                        estimated_resolution = res
                        break

        if estimated_resolution is None:
            estimated_resolution = int(DGGRID_TYPES[dggs_type]["default_res"])
        resolution = estimated_resolution

    resolution = validate_dggrid_resolution(dggs_type, resolution)

    # Build GeoDataFrames per geometry type and concatenate for performance
    dggrid_rows = []
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
            gdf_result = point2dggrid(
                dggrid_instance,
                dggs_type,
                geom,
                resolution,
                predicate=predicate,
                compact=compact,
                topology=False,  # topology already processed above
                include_properties=include_properties,
                feature_properties=props,
                output_address_type=output_address_type,
                split_antimeridian=split_antimeridian,
                aggregate=aggregate,
                options=options,
            )
            if not gdf_result.empty:
                dggrid_rows.append(gdf_result)

        elif geom.geom_type in ("LineString", "MultiLineString"):
            gdf_result = polyline2dggrid(
                dggrid_instance,
                dggs_type,
                geom,
                resolution,
                predicate,
                compact,
                topology,
                include_properties,
                props,
                output_address_type,
                split_antimeridian=split_antimeridian,
                aggregate=aggregate,
                options=options,
            )
            if not gdf_result.empty:
                dggrid_rows.append(gdf_result)

        elif geom.geom_type in ("Polygon", "MultiPolygon"):
            gdf_result = polygon2dggrid(
                dggrid_instance,
                dggs_type,
                geom,
                resolution,
                predicate,
                compact,
                topology,
                include_properties,
                props,
                output_address_type,
                split_antimeridian=split_antimeridian,
                aggregate=aggregate,
                options=options,
            )
            if not gdf_result.empty:
                dggrid_rows.append(gdf_result)

    if dggrid_rows:
        final_grid = gpd.GeoDataFrame(
            pd.concat(dggrid_rows, ignore_index=True), crs=dggrid_rows[0].crs
        )
    else:
        final_grid = gpd.GeoDataFrame(columns=["geometry"], crs="EPSG:4326")

    return final_grid


def vector2dggrid(
    dggrid_instance,
    dggs_type,
    vector_data,
    resolution=None,
    predicate=None,
    compact=False,
    topology=False,
    include_properties=True,
    output_address_type="SEQNUM",
    output_format="gpd",
    split_antimeridian=False,
    aggregate=False,
    options=None,
    **kwargs,
):
    """
    Convert vector data to DGGRID grid cells from various input formats.
    If output_format is a file-based format (csv, geojson, shapefile, gpkg, parquet, geoparquet),
    the output will be saved to a file in the current directory with a default name based on the input.
    Otherwise, returns a Python object (GeoDataFrame, dict, etc.) depending on output_format.

    Args:
        data: Input data (file path, URL, GeoDataFrame, GeoJSON, etc.)
        dggrid_instance: DGGRIDv7 instance for grid operations.
        dggs_type: DGGS type (e.g., ISEA4H, FULLER, etc.)
        resolution: Resolution for the DGGRID
        address_type: Output address type (default: SEQNUM)
        output_format: Output format (gpd, geojson, csv, etc.)
        include_properties: Whether to include original feature properties
        split_antimeridian: When True, apply antimeridian fixing to the resulting polygons.
        **kwargs: Additional arguments passed to process_input_data_vector

    Returns:
        GeoDataFrame or file path depending on output_format
    """
    dggs_type = validate_dggrid_type(dggs_type)
    if not topology and resolution is None:
        raise ValueError("resolution parameter is required when topology=False")
    if resolution is not None:
        resolution = validate_dggrid_resolution(dggs_type, resolution)

    gdf = process_input_data_vector(vector_data, **kwargs)
    result = geodataframe2dggrid(
        dggrid_instance,
        dggs_type,
        gdf,
        resolution,
        predicate,
        compact,
        topology,
        include_properties,
        output_address_type,
        split_antimeridian=split_antimeridian,
        aggregate=aggregate,
        options=options,
    )

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(vector_data, str):
            base = os.path.splitext(os.path.basename(vector_data))[0]
            output_name = f"{base}2dggrid_{dggs_type}_{resolution}"
        else:
            output_name = f"dggrid_{dggs_type}_{resolution}"

    return convert_to_output_format(result, output_format, output_name)
    # return result


def vector2dggrid_cli():
    parser = argparse.ArgumentParser(
        description="Convert vector data to DGGRID grid cells"
    )
    parser.add_argument("-i", "--input", help="Input file path or URL")
    parser.add_argument(
        "-dggs",
        dest="dggs_type",
        type=str,
        required=True,
        choices=DGGRID_TYPES.keys(),
        help="DGGRID DGGS type",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Resolution (integer)",
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
        "-a",
        "--output_address_type",
        choices=output_address_types,
        default="SEQNUM",
        nargs="?",  # This makes the argument optional
        help="Select an output address type from the available options.",
    )

    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Apply antimeridian fixing to the resulting polygons",
    )
    parser.add_argument(
        "-aggregate",
        "--aggregate",
        action="store_true",
        help="Aggregate the resulting polygons",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to grid_cell_polygons_for_extent or grid_cell_polygons_from_cellids. "
        "Example: '{\"densification\": 2}'",
    )
    args = parser.parse_args()
    dggrid_instance = create_dggrid_instance()

    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}", file=sys.stderr)
            sys.exit(1)

    try:
        result = vector2dggrid(
            dggrid_instance=dggrid_instance,
            dggs_type=args.dggs_type,
            vector_data=args.input,
            resolution=args.resolution,
            predicate=args.predicate,
            compact=args.compact,
            topology=args.topology,
            include_properties=args.include_properties,
            output_address_type=args.output_address_type,
            output_format=args.output_format,
            split_antimeridian=args.split_antimeridian,
            aggregate=args.aggregate,
            options=options,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    vector2dggrid_cli()
