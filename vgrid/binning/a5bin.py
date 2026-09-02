"""
A5 Grid Binning Module

Bins point data into A5 (Adaptive 5) grid cells and computes various statistics using hierarchical geospatial indexing.

Key Functions:
- a5_bin(): Core binning function with spatial joins and aggregation
- a5bin(): Main user-facing function with multiple input/output formats
- a5bin_cli(): Command-line interface for binning functionality
"""

import os
import argparse
import json
import geopandas as gpd
from vgrid.generator.a5grid import a5_grid
from vgrid.utils.io import (
    process_input_data_bin,
    convert_to_output_format,
    validate_a5_resolution,
    aggregate_joined,
    validate_agg,
)
from vgrid.utils.constants import AGG_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS


def a5_bin(
    data,
    resolution,
    agg="count",
    category_col=None,
    numeric_col=None,
    lat_col="lat",
    lon_col="lon",
    options=None,
    split_antimeridian=False,
    verbose=True,
    **kwargs,
):
    """
    Bin point data into A5 grid cells and compute statistics using a single
    grid generation + spatial join, followed by pandas groupby aggregation.

    Returns a GeoDataFrame with A5 cell stats and geometry.
    options : dict, optional
        Options for a52geo.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
    """
    resolution = validate_a5_resolution(int(resolution))

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

    # 2) Generate A5 grid covering the points' bounding box
    minx, miny, maxx, maxy = points_gdf.total_bounds  # lon/lat order
    id_col = "a5"
    grid_gdf = a5_grid(
        resolution=resolution,
        bbox=(minx, miny, maxx, maxy),
        options=options,
        split_antimeridian=split_antimeridian,
        verbose=verbose,
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

    # 5) Join stats back to grid (keep per-cell metrics from a5_grid / geodesic_dggs_metrics)
    out = grid_gdf.merge(grouped, on=id_col, how="inner")
    if "resolution" not in out.columns:
        out["resolution"] = resolution
    result_gdf = gpd.GeoDataFrame(
        out, geometry="geometry", crs=grid_gdf.crs or "EPSG:4326"
    )
    return result_gdf


def a5bin(
    data,
    resolution,
    agg="count",
    category_col=None,
    numeric_col=None,
    output_format="gpd",
    options=None,
    split_antimeridian=False,
    verbose=True,
    **kwargs,
):
    """
    Bin point data into A5 grid cells and compute statistics from various input formats.

    This is the main function that handles binning of point data to A5 grid cells.
    It supports multiple input formats including file paths, URLs, DataFrames, GeoDataFrames,
    GeoJSON dictionaries, and lists of features.

    Args:
        data: Input data in one of the following formats:
            - File path (str): Path to vector file (shapefile, GeoJSON, etc.)
            - URL (str): URL to vector data
            - pandas.DataFrame: DataFrame with lat/lon columns
            - geopandas.GeoDataFrame: GeoDataFrame with point geometries
            - dict: GeoJSON dictionary
            - list: List of GeoJSON feature dictionaries
        resolution (int): A5 resolution level [0..29] (0=coarsest, 29=finest)
        agg (str): Statistic to compute:
            - 'count': Count of points in each cell
            - 'sum': Sum of field values
            - 'min': Minimum field value
            - 'max': Maximum field value
            - 'mean': Mean field value
            - 'median': Median field value
            - 'std': Standard deviation of field values
            - 'var': Variance of field values
            - 'range': Range of field values
            - 'minority': Least frequent value
            - 'majority': Most frequent value
            - 'variety': Number of unique values
        category_col (str, optional): Category column for grouping statistics. When provided,
            statistics are computed separately for each category value.
        numeric_col (str, optional): Numeric field to compute statistics (required if agg != 'count')
        output_format (str, optional): Output format. Options include:
            - 'gpd', 'geopandas', 'gdf', 'geodataframe': Return GeoDataFrame
            - 'geojson_dict', 'json_dict': Return GeoJSON dictionary
            - 'geojson', 'json': Save as GeoJSON file or return string
            - 'csv': Save as CSV file or return string
            - 'shp', 'shapefile': Save as shapefile
            - 'gpkg', 'geopackage': Save as GeoPackage
            - 'parquet', 'geoparquet': Save as Parquet file
            - None: Return list of dictionaries
        options : dict, optional
            Options for a52geo.
        split_antimeridian : bool, optional
            When True, apply antimeridian fixing to the resulting polygons.
        **kwargs: Additional arguments passed to geopandas read functions (e.g., lat_col, lon_col)

    Returns:
        Various types depending on output_format:
        - GeoDataFrame: When output_format is 'gpd', 'geopandas', 'gdf', 'geodataframe'
        - dict: When output_format is 'geojson_dict', 'json_dict', or None
        - str: When output_format is 'geojson', 'json', or 'csv' (returns data as string)
        - str: File path when output_format is a file-based format (geojson, csv, shp, gpkg, parquet)

    Raises:
        ValueError: If input data type is not supported, conversion fails, or required parameters are missing
        TypeError: If resolution is not an integer

    Example:
        >>> # Bin from file with count statistics
        >>> result = a5bin("cities.shp", 10, "count")

        >>> # Bin from GeoDataFrame with mean statistics
        >>> import geopandas as gpd
        >>> gdf = gpd.read_file("cities.shp")
        >>> result = a5bin(gdf, 10, "mean", numeric_col="population")

        >>> # Bin from GeoJSON dict with category grouping
        >>> geojson = {"type": "FeatureCollection", "features": [...]}
        >>> result = a5bin(geojson, 10, "sum", numeric_col="value", category_col="type")

        >>> # Save output as GeoJSON file
        >>> result = a5bin("points.csv", 8, "count", output_format="geojson")
        >>> print(f"Output saved to: {result}")
    """

    if not validate_agg(agg):
        raise ValueError(f"Invalid aggregation '{agg}'")
    if agg != "count" and not numeric_col:
        raise ValueError("A numeric_col is required for statistics other than 'count'")

    # Process input data and bin
    result_gdf = a5_bin(
        data,
        resolution,
        agg,
        category_col,
        numeric_col,
        options=options,
        split_antimeridian=split_antimeridian,
        verbose=verbose,
        **kwargs,
    )

    # Convert to output output_format if specified
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(data, str):
            base = os.path.splitext(os.path.basename(data))[0]
            output_name = f"{base}_a5bin_{resolution}"
        else:
            output_name = f"a5bin_{resolution}"
    return convert_to_output_format(result_gdf, output_format, output_name)


def a5bin_cli():
    """
    Command-line interface for a5bin conversion.

    This function provides a command-line interface for binning point data to A5 grid cells.
    It parses command-line arguments and calls the main a5bin function.

    Usage:
        python a5bin.py -i input.shp -r 10 -agg count -f geojson -o output.geojson

    Arguments:
        -i, --input: Input file path, URL, or other vector file formats
        -r, --resolution: A5 resolution [0..29]
        -agg, --agg: Statistic to compute (count, min, max, sum, mean, median, std, var, range, minority, majority, variety)
        -category, --category: Optional category field for grouping
        -numeric_col, --numeric_col: Numeric field to compute statistics (required if agg != 'count')
        -split, --split_antimeridian: Apply antimeridian fixing to the resulting polygons
        -o, --output: Output file path (optional, will auto-generate if not provided)
        -f, --output_format: Output output_format (geojson, gpkg, parquet, csv, shapefile)

    Example:
        >>> # Bin shapefile to A5 cells at resolution 10 with count statistics
        >>> # python a5bin.py -i cities.shp -r 10 -agg count -f geojson
    """
    parser = argparse.ArgumentParser(description="Binning point data to A5 DGGS")
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
        help="Resolution of the grid [0..29]",
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
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to a52geo. "
        "Example: '{\"segments\": 1000}'",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Apply antimeridian fixing to the resulting polygons",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show progress bar (default: True). Use --no-verbose to hide it.",
    )

    args = parser.parse_args()

    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return

    try:
        # Use the a5bin function
        result = a5bin(
            data=args.input,
            resolution=args.resolution,
            agg=args.agg,
            category_col=args.category_col,
            numeric_col=args.numeric_col,
            output_format=args.output_format,
            options=options,
            split_antimeridian=args.split_antimeridian,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
        # Print notification is now handled in convert_to_output_format
    except Exception as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    a5bin_cli()
