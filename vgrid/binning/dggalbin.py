"""
DGGAL Grid Binning Module

Bins point data into DGGAL (Discrete Global Grids with Adaptive Localization) cells and computes various statistics for multiple grid types including ISEA3H, ISEA9R, IVEA3H, IVEA9R, RTEA3H, RTEA9R, and rHEALPix.

Key Functions:
- dggal_bin(): Core binning function with spatial joins and aggregation
- dggalbin(): Main user-facing function with multiple input/output formats
- dggalbin_cli(): Command-line interface for binning functionality
"""

import os
import argparse
import geopandas as gpd
from vgrid.generator.dggalgen import dggalgen
from vgrid.utils.io import (
    process_input_data_bin,
    convert_to_output_format,
    validate_dggal_resolution,
    aggregate_joined,
)
from vgrid.utils.constants import (
    STATS_OPTIONS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    DGGAL_TYPES,
)


def dggal_bin(
    dggs_type: str,
    data,
    resolution: int,
    stats: str = "count",
    category_col: str | None = None,
    numeric_col: str | None = None,
    lat_col: str = "lat",
    lon_col: str = "lon",
    split_antimeridian: bool = False,
    **kwargs,
):
    """
    Bin point data into DGGAL grid cells and compute statistics using a single
    grid generation + spatial join, followed by pandas groupby aggregation.

    This avoids per-point subprocess calls and is significantly faster.

    Returns a GeoDataFrame with DGGAL cell stats and geometry.
    """

    resolution = validate_dggal_resolution(dggs_type, int(resolution))

    if stats != "count" and not numeric_col:
        raise ValueError(
            "A numeric_col is required for statistics other than 'count'"
        )

    # 1) Normalize input to GeoDataFrame of points
    points_gdf = process_input_data_bin(
        data, lat_col=lat_col, lon_col=lon_col, **kwargs
    )
    # Keep only points and multipoints; ignore others
    if not points_gdf.empty:
        points_gdf = points_gdf[
            points_gdf.geometry.geom_type.isin(["Point", "MultiPoint"])
        ].copy()
        if "MultiPoint" in set(points_gdf.geometry.geom_type.unique()):
            points_gdf = points_gdf.explode(index_parts=False, ignore_index=True)

    minx, miny, maxx, maxy = points_gdf.total_bounds  # lon/lat order
    bbox = (minx, miny, maxx, maxy)
    id_col = f"dggal_{dggs_type}"
    grid_gdf = dggalgen(
        dggs_type=dggs_type,
        resolution=resolution,
        output_format="gpd",
        bbox=bbox,
        split_antimeridian=split_antimeridian,
    )
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
    result_gdf = gpd.GeoDataFrame(out, geometry="geometry", crs=grid_gdf.crs or "EPSG:4326")
    return result_gdf


def dggalbin(
    dggs_type: str,
    data,
    resolution: int,
    stats: str = "count",
    category_col: str | None = None,
    numeric_col: str | None = None,
    output_format: str = "gpd",
    split_antimeridian: bool = False,
    **kwargs,
):
    """
    Bin point data into DGGAL grid cells and compute statistics from various input formats.
    """
    result_gdf = dggal_bin(
        dggs_type=dggs_type,
        data=data,
        resolution=resolution,
        stats=stats,
        category_col=category_col,
        numeric_col=numeric_col,
        split_antimeridian=split_antimeridian,
        **kwargs,
    )

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(data, str):
            base = os.path.splitext(os.path.basename(data))[0]
            output_name = f"{base}_{dggs_type}bin_{resolution}"
        else:
            output_name = f"{dggs_type}bin_{resolution}"
    return convert_to_output_format(result_gdf, output_format, output_name)


def dggalbin_cli():
    """Command-line interface for DGGAL binning."""
    parser = argparse.ArgumentParser(description="Binning point data to DGGAL DGGS")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input data: GeoJSON file path, URL, or other vector file formats",
    )
    parser.add_argument(
        "-t",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGAL_TYPES.keys(),
        help="DGGAL type",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Resolution (integer)",
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
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Apply antimeridian fixing to the resulting polygons",
    )

    args = parser.parse_args()

    try:
        result = dggalbin(
            dggs_type=args.dggs_type,
            data=args.input,
            resolution=args.resolution,
            stats=args.statistics,
            category_col=args.category_col,
            numeric_col=args.numeric_col,
            output_format=args.output_format,
            split_antimeridian=args.split_antimeridian,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    dggalbin_cli()
