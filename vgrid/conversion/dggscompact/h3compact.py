"""
H3 Compact Module

This module provides functionality to compact and expand H3 cells with flexible input and output formats.

Key Functions:
    h3compact: Compact a set of H3 cells to their minimal covering set
    h3expand: Expand (uncompact) H3 cells to a specified resolution
    h3compact_cli: Command-line interface for compaction
    h3expand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
import h3
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    process_input_data_compact,
    convert_to_output_format,
    validate_h3_resolution,
)
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.conversion.dggs2geo.h32geo import h32geo


def h3compact(
    input_data,
    h3_id=None,
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Compact H3 cells to their minimal covering set.

    Compacts a set of H3 cells using the H3 library's built-in compaction functionality.
    Supports flexible input and output formats.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing H3 cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of H3 cell IDs
    h3_id : str, optional
        Name of the column containing H3 cell IDs. Defaults to "h3".
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The compacted H3 cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = h3compact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = h3compact(["8e65b56628e0d07", "8e65b56628e0d08"])

    >>> # Compact to GeoJSON file
    >>> result = h3compact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if h3_id is None:
        h3_id = "h3"
    gdf = process_input_data_compact(input_data, h3_id)
    h3_ids = gdf[h3_id].drop_duplicates().tolist()
    if not h3_ids:
        print(f"No H3 IDs found in <{h3_id}> field.")
        return
    try:
        h3_ids_compact = h3.compact_cells(h3_ids)
    except Exception:
        h3_ids_compact = (
            h3_ids  # to handle "Input cells must all share the same resolution."
        )
    if not h3_ids_compact:
        return None
    # Build output GeoDataFrame
    rows = []
    for h3_id_compact in h3_ids_compact:
        try:
            cell_polygon = h32geo(h3_id_compact, fix_antimeridian=fix_antimeridian)
            cell_resolution = h3.get_resolution(h3_id_compact)
            num_edges = 6
            if h3.is_pentagon(h3_id_compact):
                num_edges = 5
            row = geodesic_dggs_to_geoseries(
                "h3", h3_id_compact, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    # If output_format is file-based, set ouput_name as just the filename in current directory
    ouput_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            ouput_name = f"{base}_h3_compacted"
        else:
            ouput_name = "h3_compacted"

    return convert_to_output_format(out_gdf, output_format, ouput_name)


def h3compact_cli():
    """
    Command-line interface for h3compact with flexible input/output.
    """
    parser = argparse.ArgumentParser(description="H3 Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input H3 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="H3 ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
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
        help="Enable Antimeridian fixing",
    )

    args = parser.parse_args()
    fix_antimeridian = args.fix_antimeridian
    input_data = args.input
    cellid = args.cellid
    output_format = args.output_format

    result = h3compact(
        input_data,
        h3_id=cellid,
        output_format=output_format,
        fix_antimeridian=fix_antimeridian,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


def h3expand(
    input_data,
    resolution,
    h3_id=None,
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Expand (uncompact) H3 cells to a target resolution.

    Expands H3 cells to their children at the specified resolution using the H3 library's
    uncompact functionality. The target resolution must be greater than or equal to the
    maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing H3 cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of H3 cell IDs
    resolution : int
        Target H3 resolution to expand the cells to. Must be >= maximum input resolution.
    h3_id : str, optional
        Name of the column containing H3 cell IDs. Defaults to "h3".
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The expanded H3 cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = h3expand("cells.geojson", resolution=5)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = h3expand(["8e65b56628e0d07"], resolution=5)

    >>> # Expand to GeoJSON file
    >>> result = h3expand("cells.geojson", resolution=5, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if h3_id is None:
        h3_id = "h3"
    resolution = validate_h3_resolution(resolution)
    gdf = process_input_data_compact(input_data, h3_id)
    h3_ids = gdf[h3_id].drop_duplicates().tolist()
    if not h3_ids:
        print(f"No H3 IDs found in <{h3_id}> field.")
        return
    try:
        max_res = max(h3.get_resolution(hid) for hid in h3_ids)
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None
        h3_ids_expand = h3.uncompact_cells(h3_ids, resolution)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your H3 ID field and resolution."
        )
    if not h3_ids_expand:
        return None
    # Build output GeoDataFrame
    rows = []
    for h3_id_expand in h3_ids_expand:
        try:
            cell_polygon = h32geo(h3_id_expand, fix_antimeridian=fix_antimeridian)
            cell_resolution = h3.get_resolution(h3_id_expand)
            num_edges = 6
            if h3.is_pentagon(h3_id_expand):
                num_edges = 5
            row = geodesic_dggs_to_geoseries(
                "h3", h3_id_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    # If output_format is file-based, set ouput_name as just the filename in current directory
    ouput_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            ouput_name = f"{base}_h3_expanded"
        else:
            ouput_name = "h3_expanded"

    return convert_to_output_format(out_gdf, output_format, ouput_name)


def h3expand_cli():
    """
    Command-line interface for h3expand with flexible input/output.
    """
    parser = argparse.ArgumentParser(description="H3 Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input H3 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target H3 resolution to expand to (must be greater than input cells)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="H3 ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
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
        help="Enable Antimeridian fixing",
    )

    args = parser.parse_args()
    input_data = args.input
    resolution = args.resolution
    cellid = args.cellid
    output_format = args.output_format

    result = h3expand(
        input_data,
        resolution,
        h3_id=cellid,
        output_format=output_format,
        fix_antimeridian=args.fix_antimeridian,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)
