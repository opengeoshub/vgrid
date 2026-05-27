"""
S2 Compact Module

This module provides functionality to compact and expand S2 cells with flexible input and output formats.

Key Functions:
    s2compact: Compact a set of S2 cells to their minimal covering set
    s2expand: Expand (uncompact) S2 cells to a specified resolution
    s2compact_cli: Command-line interface for compaction
    s2expand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from vgrid.dggs import s2
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    process_input_data_compact,
    convert_to_output_format,
    validate_s2_resolution,
)
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.conversion.dggs2geo.s22geo import s22geo


def s2compact(
    input_data,
    s2_token="s2",
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Compact S2 cells to their minimal covering set.

    Compacts a set of S2 cells using the S2 library's CellUnion functionality.
    Supports flexible input and output formats.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing S2 cell tokens. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of S2 cell tokens
    s2_token : str, default "s2"
        Name of the column containing S2 cell tokens.
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
        The compacted S2 cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = s2compact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = s2compact(["31752f45cc94", "31752f45cc95"])

    >>> # Compact to GeoJSON file
    >>> result = s2compact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    gdf = process_input_data_compact(input_data, s2_token)
    s2_tokens = gdf[s2_token].drop_duplicates().tolist()
    if not s2_tokens:
        print(f"No S2 tokens found in <{s2_token}> field.")
        return
    try:
        s2_cells = [s2.CellId.from_token(token) for token in s2_tokens]
        s2_cells = list(set(s2_cells))
        if not s2_cells:
            print(f"No valid S2 tokens found in <{s2_token}> field.")
            return
        covering = s2.CellUnion(s2_cells)
        covering.normalize()
        s2_tokens_compact = [cell_id.to_token() for cell_id in covering.cell_ids()]
    except Exception:
        raise Exception("Compact cells failed. Please check your S2 ID field.")
    if not s2_tokens_compact:
        return None
    # Build output GeoDataFrame
    rows = []
    for s2_token_compact in s2_tokens_compact:
        try:
            cell_polygon = s22geo(s2_token_compact, fix_antimeridian=fix_antimeridian)
            cell_resolution = s2.CellId.from_token(s2_token_compact).level()
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "s2", s2_token_compact, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    # If output_format is file-based, set output_name as just the filename in current directory
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_s2_compacted"
        else:
            output_name = "s2_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def s2compact_cli():
    """
    Command-line interface for s2compact with flexible input/output.
    """
    parser = argparse.ArgumentParser(description="S2 Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input S2 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="S2 ID field")
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
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    args = parser.parse_args()
    input_data = args.input
    cellid = args.cellid
    output_format = args.output_format
    fix_antimeridian = args.fix_antimeridian
    result = s2compact(
        input_data,
        s2_token=cellid,
        output_format=output_format,
        fix_antimeridian=fix_antimeridian,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


def s2expand(
    input_data,
    resolution,
    s2_token="s2",
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Expand (uncompact) S2 cells to a target resolution.

    Expands S2 cells to their children at the specified resolution. The target resolution
    must be greater than or equal to the maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing S2 cell tokens. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of S2 cell tokens
    resolution : int
        Target S2 resolution to expand the cells to. Must be >= maximum input resolution.
    s2_token : str, default "s2"
        Name of the column containing S2 cell tokens.
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
        The expanded S2 cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = s2expand("cells.geojson", resolution=10)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = s2expand(["31752f45cc94"], resolution=10)

    >>> # Expand to GeoJSON file
    >>> result = s2expand("cells.geojson", resolution=10, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    resolution = validate_s2_resolution(resolution)
    gdf = process_input_data_compact(input_data, s2_token)
    s2_tokens = gdf[s2_token].drop_duplicates().tolist()
    if not s2_tokens:
        print(f"No S2 tokens found in <{s2_token}> field.")
        return
    try:
        s2_cells = [s2.CellId.from_token(token) for token in s2_tokens]
        s2_cells = list(set(s2_cells))
        if not s2_cells:
            print(f"No valid S2 tokens found in <{s2_token}> field.")
            return
        max_res = max(cell.level() for cell in s2_cells)
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None
        # Expand each cell to the target resolution
        expanded_cells = []
        for cell in s2_cells:
            if cell.level() >= resolution:
                expanded_cells.append(cell)
            else:
                expanded_cells.extend(cell.children(resolution))
        s2_tokens_expand = [cell_id.to_token() for cell_id in expanded_cells]
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your S2 ID field and resolution."
        )
    if not s2_tokens_expand:
        return None
    # Build output GeoDataFrame
    rows = []
    for s2_token_expand in s2_tokens_expand:
        try:
            cell_polygon = s22geo(s2_token_expand, fix_antimeridian=fix_antimeridian)
            cell_resolution = resolution
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "s2", s2_token_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    # If output_format is file-based, set output_name as just the filename in current directory
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_s2_expanded"
        else:
            output_name = "s2_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def s2expand_cli():
    """
    Command-line interface for s2expand with flexible input/output.
    """
    parser = argparse.ArgumentParser(description="S2 Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input S2 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target S2 resolution to expand to (must be greater than input cells)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="S2 Token field")
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
        help="Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none",
    )
    args = parser.parse_args()
    input_data = args.input
    resolution = args.resolution
    cellid = args.cellid
    output_format = args.output_format
    fix_antimeridian = args.fix_antimeridian
    result = s2expand(
        input_data,
        resolution,
        s2_token=cellid,
        output_format=output_format,
        fix_antimeridian=fix_antimeridian,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)
