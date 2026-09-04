"""
Geohash Grid Binning Module

Bins point data into Geohash grid cells and computes various statistics using hierarchical geocoding system with alphanumeric identifiers.

Key Functions:
- geohash_bin(): Core binning function with spatial joins and aggregation
- geohashbin(): Main user-facing function with multiple input/output formats
- geohashbin_cli(): Command-line interface for binning functionality
"""

import argparse
import os
import geopandas as gpd
from vgrid.generator.geohashgrid import geohash_grid_within_bbox
from vgrid.utils.io import (
    process_input_data_bin,
    convert_to_output_format,
    validate_geohash_resolution,
    aggregate_joined,
    validate_agg,
)
from vgrid.utils.constants import AGG_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS


def geohash_bin(
    data,
    resolution,
    agg="count",
    category_col=None,
    numeric_col=None,
    lat_col="lat",
    lon_col="lon",
    verbose=True,
    **kwargs,
):
    """
    Bin point data into Geohash grid cells and compute statistics using a single
    grid generation + spatial join, followed by pandas groupby aggregation.

    Returns a GeoDataFrame with Geohash cell stats and geometry.
    """
    resolution = validate_geohash_resolution(int(resolution))

    if not validate_agg(agg):
        raise ValueError(f"Invalid aggregation '{agg}'")
    if agg != "count" and not numeric_col:
        raise ValueError("A numeric_col is required for statistics other than 'count'")

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

    # 2) Generate Geohash grid covering the points' bounding box
    minx, miny, maxx, maxy = points_gdf.total_bounds  # lon/lat order
    id_col = "geohash"
    grid_gdf = geohash_grid_within_bbox(
        resolution=resolution, bbox=(minx, miny, maxx, maxy), verbose=verbose
    )

    # 3) Spatial join points -> cells with only needed columns
    join_cols = []
    if category_col and category_col in points_gdf.columns:
        join_cols.append(category_col)
    if agg != "count" and numeric_col:
        if numeric_col not in points_gdf.columns:
            raise ValueError(f"numeric_col '{numeric_col}' not found in input data")
        join_cols.append(numeric_col)
    left = points_gdf[[c for c in ["geometry", *join_cols] if c is not None]]
    joined = gpd.sjoin(
        left, grid_gdf[[id_col, "geometry"]], how="inner", predicate="within"
    )

    # 4) Aggregate
    grouped = aggregate_joined(
        joined, id_col, agg=agg, category_col=category_col, numeric_col=numeric_col
    )
    grouped = grouped.reset_index()

    # 5) Join back to grid and return GeoDataFrame
    out = grid_gdf.merge(grouped, on=id_col, how="inner")
    if "resolution" not in out.columns:
        out["resolution"] = resolution
    result_gdf = gpd.GeoDataFrame(
        out, geometry="geometry", crs=grid_gdf.crs or "EPSG:4326"
    )
    return result_gdf


def geohashbin(
    data,
    resolution,
    agg="count",
    category_col=None,
    numeric_col=None,
    output_format="gpd",
    verbose=True,
    **kwargs,
):
    resolution = validate_geohash_resolution(resolution)
    if not validate_agg(agg):
        raise ValueError(f"Invalid aggregation '{agg}'")
    if agg != "count" and not numeric_col:
        raise ValueError("A numeric_col is required for statistics other than 'count'")
    result_gdf = geohash_bin(
        data, resolution, agg, category_col, numeric_col, verbose=verbose, **kwargs
    )
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(data, str):
            base = os.path.splitext(os.path.basename(data))[0]
            output_name = f"{base}_geohashbin_{resolution}"
        else:
            output_name = f"geohashbin_{resolution}"
    return convert_to_output_format(result_gdf, output_format, output_name)


def geohashbin_cli():
    parser = argparse.ArgumentParser(description="Binning point data to Geohash DGGS")
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
        default=6,
        help="Resolution of the grid [1..10]",
    )
    parser.add_argument(
        "-agg",
        "--agg",
        choices=AGG_OPTIONS,
        default="count",
        help="Aggregation option",
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
        help="Numeric field to compute statistics (required if agg != 'count')",
    )
    # Removed -o/--output; output is saved in CWD with predefined name
    parser.add_argument(
        "-f",
        "--output_format",
        required=False,
        default="gpd",
        choices=OUTPUT_FORMATS,
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show progress bar (default: True). Use --no-verbose to hide it.",
    )

    args = parser.parse_args()
    try:
        result = geohashbin(
            data=args.input,
            resolution=args.resolution,
            agg=args.agg,
            category_col=args.category_col,
            numeric_col=args.numeric_col,
            output_format=args.output_format,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    geohashbin_cli()
