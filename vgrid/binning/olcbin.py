"""
OLC Grid Binning Module

Bins point data into OLC (Open Location Code) grid cells and computes various statistics using human-readable location codes for global coverage.

Key Functions:
- olc_bin(): Core binning function with spatial joins and aggregation
- olcbin(): Main user-facing function with multiple input/output formats
- olcbin_cli(): Command-line interface for binning functionality
"""

import argparse
import geopandas as gpd
from vgrid.utils.io import (
    process_input_data_bin,
    convert_to_output_format,
    validate_olc_resolution,
    aggregate_joined,
    olc_resolutions,
)
from vgrid.utils.constants import STATS_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS


def olc_bin(
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
    Bin point data into OLC grid cells using grid generation + spatial join and
    aggregate with pandas groupby. Supports custom stats (range, variety, minority,
    majority). Only Point/MultiPoint geometries are considered.
    """
    resolution = validate_olc_resolution(resolution)
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

    # Generate OLC grid covering the points' bounding box
    minx, miny, maxx, maxy = points_gdf.total_bounds
    id_col = "olc"
    from vgrid.generator.olcgrid import olc_grid_within_bbox

    grid_gdf = olc_grid_within_bbox(
        resolution=resolution, bbox=(minx, miny, maxx, maxy)
    )

    # Spatial join points -> cells with only needed columns
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

    # Aggregate
    grouped = aggregate_joined(
        joined, id_col, stats=stats, category_col=category_col, numeric_col=numeric_col
    )
    grouped = grouped.reset_index()

    # Join back to grid and return GeoDataFrame
    out = grid_gdf.merge(grouped, on=id_col, how="inner")
    if "resolution" not in out.columns:
        out["resolution"] = resolution
    result_gdf = gpd.GeoDataFrame(out, geometry="geometry", crs=grid_gdf.crs or "EPSG:4326")
    return result_gdf


def olcbin(
    data,
    resolution,
    stats="count",
    category_col=None,
    numeric_col=None,
    output_format="gpd",
    **kwargs,
):
    resolution = validate_olc_resolution(resolution)
    if stats != "count" and not numeric_col:
        raise ValueError(
            "A numeric_col is required for statistics other than 'count'"
        )
    result_gdf = olc_bin(data, resolution, stats, category_col, numeric_col, **kwargs)
    output_name = None
    if output_format in OUTPUT_FORMATS:
        import os

        if isinstance(data, str):
            base = os.path.splitext(os.path.basename(data))[0]
            output_name = f"{base}_olcbin_{resolution}"
        else:
            output_name = f"olcbin_{resolution}"
    return convert_to_output_format(result_gdf, output_format, output_name)


def olcbin_cli():
    parser = argparse.ArgumentParser(description="Binning point data to OLC DGGS")
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
        choices=olc_resolutions,
        default=8,
        help=f"Resolution of the OLC DGGS (choose from {', '.join(map(str, olc_resolutions))})",
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
    # Removed -o/--output; output is saved in CWD with predefined name
    parser.add_argument(
        "-f",
        "--output_format",
        required=False,
        default="gpd",
        choices=OUTPUT_FORMATS,
    )
    args = parser.parse_args()
    try:
        result = olcbin(
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
    olcbin_cli()
