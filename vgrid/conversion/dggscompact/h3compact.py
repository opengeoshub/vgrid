"""
H3 Compact Module

This module provides functionality to compact and expand H3 cells with flexible input and output formats.

Key Functions:
    h3_compact: Compact a list of H3 IDs with an optional parent depth
    h3compact: Compact a set of H3 cells to their covering set
    h3expand: Expand (uncompact) H3 cells to a specified resolution
    h3compact_cli: Command-line interface for compaction
    h3expand_cli: Command-line interface for expansion
"""

import os
import argparse

import geopandas as gpd
import h3
from tqdm import tqdm
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    aggregate_values,
    compact_cells,
    convert_to_output_format,
    prepare_compact_bags,
    process_input_data_compact,
    validate_h3_resolution,
)
from vgrid.utils.constants import AGG_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.conversion.dggs2geo.h32geo import h32geo


def h3_compact(h3_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of H3 cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    h3_ids : list of str
        H3 cell IDs to compact. Mixed resolutions are allowed.
    depth : int, default -1
        How many parent levels to climb:
        - ``0``: do nothing (return the unique input cells)
        - ``-1``: compact as far as possible (same result as ``h3.compact_cells``
          when all inputs share a resolution)
        - ``1``: replace complete sibling sets with their direct parent
        - ``2``: then compact those parents (grandparents), and so on
    bags : dict of list, optional
        Per-cell lists of original values. When a complete child set is replaced
        by its parent, child lists are concatenated onto the parent. Mutated
        in place so remaining keys match the compacted IDs.
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

    Returns
    -------
    list of str
        Sorted compacted H3 cell IDs.
    """
    def parent_fn(h3_id):
        cell_res = h3.get_resolution(h3_id)
        if cell_res <= 0:
            return None
        return h3.cell_to_parent(h3_id, cell_res - 1)

    return compact_cells(
        h3_ids,
        parent_fn,
        h3.cell_to_children,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting H3",
    )


def h3compact(
    input_data,
    h3_id=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
):
    """
    Compact H3 cells to their covering set at a given parent depth.

    Compacts a set of H3 cells by replacing complete sets of children with their
    parent cells. Unlike ``h3.compact_cells``, mixed input resolutions are allowed
    and ``depth`` limits how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg`` (same options as ``h3bin``). If ``agg`` is
    ``"count"``, ``numeric_col`` is ignored and the output ``count`` is the
    number of original input cells in each compacted cell.

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
    depth : int, default -1
        Compaction depth: ``0`` leaves cells unchanged, ``-1`` compact as far as
        possible, ``1`` merges to the direct parent, ``2`` to the grandparent, etc.
    agg : str, default "count"
        Aggregation applied to original child values when cells compact into a
        parent. Same options as ``h3bin`` (``count``, ``min``, ``max``, ``sum``,
        ``mean``, ``median``, ``std``, ``var``, ``range``, ``minority``,
        ``majority``, ``variety``).
    numeric_col : str, optional
        Numeric field to aggregate. Required when ``agg`` is not ``"count"``;
        ignored when ``agg`` is ``"count"``.
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

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

    >>> # Compact only one parent level
    >>> result = h3compact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = h3compact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = h3compact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if h3_id is None:
        h3_id = "h3"
    bags, agg_col = prepare_compact_bags(
        input_data,
        h3_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="H3 cells",
    )
    if bags is None:
        print(f"No H3 IDs found in <{h3_id}> field.")
        return

    h3_ids_compact = h3_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not h3_ids_compact:
        return None

    rows = []
    for h3_id_compact in tqdm(
        h3_ids_compact,
        desc="Building H3 compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = h32geo(h3_id_compact, fix_antimeridian=fix_antimeridian)
            cell_resolution = h3.get_resolution(h3_id_compact)
            num_edges = 6
            if h3.is_pentagon(h3_id_compact):
                num_edges = 5
            row = geodesic_dggs_to_geoseries(
                "h3", h3_id_compact, cell_resolution, cell_polygon, num_edges
            )
            row[agg_col] = aggregate_values(bags.get(h3_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

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
    parser.add_argument(
        "-d",
        "--depth",
        type=int,
        default=-1,
        help="Compaction depth: 0 = no-op, -1 = compact fully (default), "
        "1 = direct parent, 2 = grandparent, ...",
    )
    parser.add_argument(
        "-agg",
        "--agg",
        choices=AGG_OPTIONS,
        default="count",
        help="Aggregation option",
    )
    parser.add_argument(
        "-numeric_col",
        "--numeric_col",
        dest="numeric_col",
        required=False,
        help="Numeric field to aggregate (required if agg != 'count')",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show progress bar (default: True). Use --no-verbose to hide it.",
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
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


def h3expand(
    input_data,
    resolution,
    h3_id=None,
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
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
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

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
    rows = []
    for h3_id_expand in tqdm(
        h3_ids_expand,
        desc="Building H3 expand",
        unit=" cells",
        disable=not verbose,
    ):
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
    parser.add_argument(
        "-v",
        "--verbose",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Show progress bar (default: True). Use --no-verbose to hide it.",
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
        verbose=args.verbose,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)
