"""
Maidenhead Grid Binning Module

Bins point data into Maidenhead locator cells using
:func:`~vgrid.generator.maidenheadgrid.maidenhead_grid_within_bbox` over the points'
bounding box and the same aggregation pattern as :mod:`geohashbin`.

Key Functions:
- maidenhead_bin(): Core binning with spatial join and aggregation
- maidenheadbin(): User-facing function with multiple output formats
- maidenheadbin_cli(): Command-line interface
"""

import argparse
import os

import geopandas as gpd

from vgrid.generator.maidenheadgrid import maidenhead_grid_within_bbox
from vgrid.utils.constants import (
    DGGS_TYPES,
    OUTPUT_FORMATS,
    STATS_OPTIONS,
    STRUCTURED_FORMATS,
)
from vgrid.utils.io import (
    convert_to_output_format,
    process_input_data_bin,
    validate_maidenhead_resolution,
    aggregate_joined,
)

min_res = DGGS_TYPES["maidenhead"]["min_res"]
max_res = DGGS_TYPES["maidenhead"]["max_res"]
default_res = DGGS_TYPES["maidenhead"]["default_res"]


def maidenhead_bin(
    data,
    resolution,
    stats="count",
    category_col=None,
    numeric_col=None,
    lat_col="lat",
    lon_col="lon",
    **kwargs,
):
    """
    Bin point data into Maidenhead cells and compute statistics (grid + spatial join + groupby).

    Returns a GeoDataFrame with Maidenhead cell stats and geometry (EPSG:4326).
    """
    resolution = validate_maidenhead_resolution(int(resolution))

    if stats != "count" and not numeric_col:
        raise ValueError("A numeric_col is required for statistics other than 'count'")

    points_gdf = process_input_data_bin(
        data, lat_col=lat_col, lon_col=lon_col, **kwargs
    )
    if not points_gdf.empty:
        points_gdf = points_gdf[
            points_gdf.geometry.geom_type.isin(["Point", "MultiPoint"])
        ].copy()
        if "MultiPoint" in set(points_gdf.geometry.geom_type.unique()):
            points_gdf = points_gdf.explode(index_parts=False, ignore_index=True)

    minx, miny, maxx, maxy = points_gdf.total_bounds
    id_col = "maidenhead"
    grid_gdf = maidenhead_grid_within_bbox(resolution, bbox=(minx, miny, maxx, maxy))

    join_cols = []
    if category_col and category_col in points_gdf.columns:
        join_cols.append(category_col)
    if stats != "count" and numeric_col:
        if numeric_col not in points_gdf.columns:
            raise ValueError(f"numeric_col '{numeric_col}' not found in input data")
        join_cols.append(numeric_col)
    left = points_gdf[[c for c in ["geometry", *join_cols] if c is not None]]
    joined = gpd.sjoin(
        left, grid_gdf[[id_col, "geometry"]], how="inner", predicate="within"
    )

    grouped = aggregate_joined(
        joined, id_col, stats=stats, category_col=category_col, numeric_col=numeric_col
    )
    grouped = grouped.reset_index()

    out = grid_gdf.merge(grouped, on=id_col, how="inner")
    if "resolution" not in out.columns:
        out["resolution"] = resolution
    result_gdf = gpd.GeoDataFrame(
        out, geometry="geometry", crs=grid_gdf.crs or "EPSG:4326"
    )
    return result_gdf


def maidenheadbin(
    data,
    resolution,
    stats="count",
    category_col=None,
    numeric_col=None,
    output_format="gpd",
    **kwargs,
):
    resolution = validate_maidenhead_resolution(resolution)
    if stats != "count" and not numeric_col:
        raise ValueError("A numeric_col is required for statistics other than 'count'")
    result_gdf = maidenhead_bin(
        data, resolution, stats, category_col, numeric_col, **kwargs
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(data, str):
            base = os.path.splitext(os.path.basename(data))[0]
            output_name = f"{base}_maidenheadbin_{resolution}"
        else:
            output_name = f"maidenheadbin_{resolution}"
    return convert_to_output_format(result_gdf, output_format, output_name)


def maidenheadbin_cli():
    parser = argparse.ArgumentParser(
        description="Binning point data to Maidenhead DGGS"
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input data: GeoJSON file path, URL, or other vector file formats",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=default_res,
        help=f"Maidenhead resolution [{min_res}..{max_res}]",
    )
    parser.add_argument(
        "-stats",
        "--statistics",
        choices=STATS_OPTIONS,
        default="count",
        help="Statistic option",
    )
    parser.add_argument(
        "-category",
        "--category",
        dest="category_col",
        required=False,
        help="Optional category field for grouping",
    )
    parser.add_argument(
        "-numeric_col",
        "--numeric_col",
        dest="numeric_col",
        required=False,
        help="Numeric field to compute statistics (required if stats != 'count')",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        required=False,
        default="gpd",
        choices=OUTPUT_FORMATS,
    )
    args = parser.parse_args()
    try:
        result = maidenheadbin(
            data=args.input,
            resolution=args.resolution,
            stats=args.statistics,
            category_col=args.category_col,
            numeric_col=args.numeric_col,
            output_format=args.output_format,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    maidenheadbin_cli()
