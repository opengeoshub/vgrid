"""
S2 Grid Binning Module

Bins point data into S2 spherical grid cells and computes various statistics using Google's hierarchical grid system.

Key Functions:
- s2_bin(): Core binning function with spatial joins and aggregation
- s2bin(): Main user-facing function with multiple input/output formats
- s2bin_cli(): Command-line interface for binning functionality
"""

import argparse
import geopandas as gpd
from vgrid.utils.io import (
    process_input_data_bin,
    convert_to_output_format,
    validate_s2_resolution,
    aggregate_joined,
    validate_agg,
)
from vgrid.utils.constants import (
    AGG_OPTIONS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
    FIX_ANTIMERIDIAN_CHOICES,
)


def s2_bin(
    data,
    resolution,
    agg="count",
    category_col=None,
    numeric_col=None,
    lat_col="lat",
    lon_col="lon",
    fix_antimeridian=None,
    verbose=True,
    **kwargs,
):
    """
    Grid + spatial join + groupby approach for S2 binning (like a5bin).
    """
    resolution = validate_s2_resolution(resolution)
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

    # Generate S2 grid covering the points' bounding box
    minx, miny, maxx, maxy = points_gdf.total_bounds
    id_col = "s2"
    from vgrid.generator.s2grid import s2_grid

    grid_gdf = s2_grid(
        resolution=resolution,
        bbox=(minx, miny, maxx, maxy),
        fix_antimeridian=fix_antimeridian,
        verbose=verbose,
    )

    # Spatial join points -> cells with only needed columns
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

    # Aggregate
    grouped = aggregate_joined(
        joined, id_col, agg=agg, category_col=category_col, numeric_col=numeric_col
    )
    grouped = grouped.reset_index()

    # Join back to grid and return GeoDataFrame
    out = grid_gdf.merge(grouped, on=id_col, how="inner")
    if "resolution" not in out.columns:
        out["resolution"] = resolution
    result_gdf = gpd.GeoDataFrame(
        out, geometry="geometry", crs=grid_gdf.crs or "EPSG:4326"
    )
    return result_gdf


def s2bin(
    data,
    resolution,
    agg="count",
    category_col=None,
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
    **kwargs,
):
    resolution = validate_s2_resolution(resolution)
    if not validate_agg(agg):
        raise ValueError(f"Invalid aggregation '{agg}'")
    if agg != "count" and not numeric_col:
        raise ValueError("A numeric_col is required for statistics other than 'count'")
    # Process input data and bin
    result_gdf = s2_bin(
        data,
        resolution,
        agg,
        category_col,
        numeric_col,
        fix_antimeridian=fix_antimeridian,
        verbose=verbose,
        **kwargs,
    )
    # Convert to output output_format if specified
    output_name = None
    if output_format in OUTPUT_FORMATS:
        import os

        if isinstance(data, str):
            base = os.path.splitext(os.path.basename(data))[0]
            output_name = f"{base}_s2bin_{resolution}"
        else:
            output_name = f"s2bin_{resolution}"
    return convert_to_output_format(result_gdf, output_format, output_name)


def s2bin_cli():
    """
    Command-line interface for s2bin conversion.

    This function provides a command-line interface for binning point data to S2 grid cells.
    It parses command-line arguments and calls the main s2bin function.

    Usage:
            python s2bin.py -i input.shp -r 10 -agg count -f geojson

    CLI Arguments:
            -i, --input: Input file path, URL, or other vector file formats
            -r, --resolution: S2 resolution [0..30]
            -agg, --agg: Statistic to compute (count, min, max, sum, mean, median, std, var, range, minority, majority, variety)
            -category, --category: Optional category field for grouping
            -numeric_col, --numeric_col: Numeric field to compute statistics (required if agg != 'count')
            -f, --output_format: Output format (geojson, gpkg, parquet, csv, shapefile)

    Example:
            >>> # Bin shapefile to S2 cells at resolution 10 with count statistics
            >>> # python s2bin.py -i cities.shp -r 10 -agg count -f geojson
    """
    parser = argparse.ArgumentParser(description="Binning point data to S2 DGGS")
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
        default=13,
        help="Resolution of the grid [0..30]",
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
        "-fix",
        "--fix_antimeridian",
        type=str,
        choices=FIX_ANTIMERIDIAN_CHOICES,
        default=None,
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
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
        # Use the s2bin function
        result = s2bin(
            data=args.input,
            resolution=args.resolution,
            agg=args.agg,
            category_col=args.category_col,
            numeric_col=args.numeric_col,
            output_format=args.output_format,
            fix_antimeridian=args.fix_antimeridian,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except Exception as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    s2bin_cli()
