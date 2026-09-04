"""
A5 Compact Module

This module provides functionality to compact and expand A5 cells with flexible input and output formats.

Key Functions:
    a5_compact: Compact a list of A5 hex IDs with an optional parent depth
    a5compact: Compact a set of A5 cells to their covering set
    a5expand: Expand (uncompact) A5 cells to a target resolution or by depth
    a5compact_cli: Command-line interface for compaction
    a5expand_cli: Command-line interface for expansion
"""

import os
import argparse
import json
import geopandas as gpd
import a5
from tqdm import tqdm
from vgrid.conversion.dggs2geo.a52geo import a52geo
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    add_verbose_argument,
    aggregate_values,
    compact_cells,
    convert_to_output_format,
    prepare_compact_bags,
    process_input_data_compact,
    validate_dggs_compact_depth,
    validate_dggs_expand_depth,
    validate_dggs_expand_resolution,
)
from vgrid.utils.constants import AGG_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS


def _a5_parent(hex_id):
    u = a5.hex_to_u64(hex_id)
    if a5.get_resolution(u) <= 0:
        return None
    return a5.u64_to_hex(a5.cell_to_parent(u))


def _a5_children(hex_id):
    u = a5.hex_to_u64(hex_id)
    return {a5.u64_to_hex(c) for c in a5.cell_to_children(u)}


def a5_compact(a5_hexes, depth=-1, bags=None, verbose=True):
    """
    Compact a list of A5 hex cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    a5_hexes : list of str
        A5 hex string cell IDs to compact. Mixed resolutions are allowed.
    depth : int, default -1
        How many parent levels to climb (``-1`` to ``max_res``):
        - ``0``: do nothing (return the unique input cells)
        - ``-1``: compact as far as possible
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
        Sorted compacted A5 hex string cell IDs.
    """
    depth = validate_dggs_compact_depth("a5", depth)
    return compact_cells(
        a5_hexes,
        _a5_parent,
        _a5_children,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting A5",
    )


def a5_expand(a5_hexes, resolution=None, depth=None, verbose=True):
    """
    Expand A5 hex strings to a target resolution, or by a relative child depth.

    When ``resolution`` is set, ``depth`` is ignored and all cells are uncompacted
    to that absolute resolution. When only ``depth`` is set, ``resolution`` is
    ignored and each cell (at any resolution) is expanded to its descendants
    ``depth`` levels down (``1`` = direct children, ``2`` = grandchildren, and
    so on).

    Parameters
    ----------
    a5_hexes : list of str
        List of A5 hex string cell IDs. Mixed resolutions are allowed when
        expanding by depth.
    resolution : int, optional
        Target A5 resolution to expand all cells to. When set, ``depth`` is
        ignored.
    depth : int, optional
        Relative expansion depth (``1 <= depth <= max_res``). Used when
        ``resolution`` is not set. ``1`` expands each cell to its direct
        children; ``2`` to the next level, and so on.

    Returns
    -------
    list of str
        List of expanded A5 hex string cell IDs.

    Examples
    --------
    >>> hexes = ["8e65b56628e0d07"]
    >>> expanded = a5_expand(hexes, resolution=5)
    >>> children = a5_expand(hexes, depth=1)
    >>> grandchildren = a5_expand(hexes, depth=2)
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("a5", resolution)
        a5_u64s = [a5.hex_to_u64(a5_hex) for a5_hex in a5_hexes]
        a5_u64s_expand = a5.core.compact.uncompact(a5_u64s, resolution)
        return [a5.u64_to_hex(u64) for u64 in a5_u64s_expand]

    if depth is None:
        raise ValueError("Either resolution or depth must be specified.")
    depth = validate_dggs_expand_depth("a5", depth)
    a5_hexes_expand = []
    for a5_hex in tqdm(a5_hexes, desc="Expanding A5", unit=" cells", disable=not verbose):
        try:
            u = a5.hex_to_u64(a5_hex)
            child_res = a5.get_resolution(u) + depth
            a5_hexes_expand.extend(
                a5.u64_to_hex(c) for c in a5.cell_to_children(u, child_res)
            )
        except Exception:
            continue
    return a5_hexes_expand


def a5compact(
    input_data,
    a5_hex=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    options=None,
    split_antimeridian=False,
    verbose=True,
):
    """
    Compact A5 cells to their covering set at a given parent depth.

    Compacts a set of A5 cells by replacing complete sets of children with their
    parent cells. Mixed input resolutions are allowed and ``depth`` limits how
    far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg``. If ``agg`` is ``"count"``, ``numeric_col`` is
    ignored and the output ``count`` is the number of original input cells in
    each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing A5 cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of A5 cell IDs
    a5_hex : str, optional
        Name of the column containing A5 cell IDs. Defaults to "a5".
    depth : int, default -1
        Compaction depth (``-1`` to ``max_res``): ``0`` leaves cells unchanged,
        ``-1`` compact as far as possible, ``1`` merges to the direct parent,
        ``2`` to the grandparent, etc.
    agg : str, default "count"
        Aggregation applied to original child values when cells compact into a
        parent. Same options as DGGS binning (``count``, ``min``, ``max``,
        ``sum``, ``mean``, ``median``, ``std``, ``var``, ``range``,
        ``minority``, ``majority``, ``variety``).
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
    options : dict, optional
        Options for a52geo.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The compacted A5 cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = a5compact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = a5compact(["8e65b56628e0d07", "8e65b56628e0d08"])

    >>> # Compact only one parent level
    >>> result = a5compact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = a5compact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = a5compact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not a5_hex:
        a5_hex = "a5"
    bags, agg_col = prepare_compact_bags(
        input_data,
        a5_hex,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="A5 cells",
    )
    if bags is None:
        print(f"No A5 IDs found in <{a5_hex}> field.")
        return

    a5_hexes_compact = a5_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not a5_hexes_compact:
        return None
    rows = []
    for a5_hex_compact in tqdm(
        a5_hexes_compact,
        desc="Building A5 compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = a52geo(
                a5_hex_compact, options, split_antimeridian=split_antimeridian
            )
            cell_resolution = a5.get_resolution(a5.hex_to_u64(a5_hex_compact))
            num_edges = 5  # A5 cells are pentagons
            if cell_resolution == 1:
                num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "a5", a5_hex_compact, cell_resolution, cell_polygon, num_edges
            )
            row[agg_col] = aggregate_values(bags.get(a5_hex_compact, []), agg)
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    ouput_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            ouput_name = f"{base}_a5_compacted"
        else:
            ouput_name = "a5_compacted"
    return convert_to_output_format(out_gdf, output_format, ouput_name)


def a5compact_cli():
    """
    Command-line interface for a5compact with flexible input/output.
    """
    parser = argparse.ArgumentParser(description="A5 Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input A5 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="A5 Hex field")
    parser.add_argument(
        "-f", "--output_format", type=str, default="gpd", choices=OUTPUT_FORMATS
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Enable Antimeridian splitting",
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
        "-d",
        "--depth",
        type=int,
        default=-1,
        help="Compaction depth [-1, A5 max_res]: 0 = no-op, -1 = compact fully "
        "(default), 1 = direct parent, 2 = grandparent, ...",
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
    input_data = args.input
    cellid = args.cellid
    output_format = args.output_format
    split_antimeridian = args.split_antimeridian

    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return

    result = a5compact(
        input_data,
        a5_hex=cellid,
        output_format=output_format,
        options=options,
        split_antimeridian=split_antimeridian,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


def a5expand(
    input_data,
    resolution=None,
    a5_hex=None,
    output_format="gpd",
    options=None,
    split_antimeridian=False,
    verbose=True,
    depth=None,
):
    """
    Expand (uncompact) A5 cells to a target resolution or by a relative depth.

    When ``resolution`` is set, ``depth`` is ignored and cells are expanded to
    that absolute resolution (must be >= the maximum input resolution). When
    only ``depth`` is set, ``resolution`` is ignored: mixed-resolution input is
    allowed and each cell is expanded to its descendants ``depth`` levels down
    (``1`` = direct children, ``2`` = grandchildren, and so on).

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing A5 cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of A5 cell IDs
    resolution : int, optional
        Target A5 resolution to expand the cells to. Must be >= maximum input
        resolution. When set, ``depth`` is ignored.
    a5_hex : str, optional
        Name of the column containing A5 cell IDs. Defaults to "a5".
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path
    options : dict, optional
        Options for a52geo.
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    depth : int, optional
        Relative expansion depth (``1 <= depth <= max_res``). Used when
        ``resolution`` is not set. Each input cell is expanded ``depth``
        levels: ``1`` = direct children, ``2`` = grandchildren, and so on.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The expanded A5 cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = a5expand("cells.geojson", resolution=5)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = a5expand(["8e65b56628e0d07"], resolution=5)

    >>> # Expand mixed-resolution cells by relative depth
    >>> result = a5expand(cells, depth=1)
    >>> result = a5expand(cells, depth=2)

    >>> # Expand to GeoJSON file
    >>> result = a5expand("cells.geojson", resolution=5, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if a5_hex is None:
        a5_hex = "a5"
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("a5", resolution)
    elif depth is not None:
        depth = validate_dggs_expand_depth("a5", depth)
    else:
        raise ValueError("Either resolution or depth must be specified.")
    gdf = process_input_data_compact(input_data, a5_hex)
    a5_hexes = gdf[a5_hex].drop_duplicates().tolist()
    if not a5_hexes:
        print(f"No A5 Hexes found in <{a5_hex}> field.")
        return
    try:
        if resolution is not None:
            max_res = max(
                a5.get_resolution(a5.hex_to_u64(a5_hex)) for a5_hex in a5_hexes
            )
            if resolution < max_res:
                print(f"Target expand resolution ({resolution}) must >= {max_res}.")
                return None
            a5_hexes_expand = a5_expand(a5_hexes, resolution=resolution, verbose=verbose)
        else:
            a5_hexes_expand = a5_expand(a5_hexes, depth=depth, verbose=verbose)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your A5 ID field, resolution, or depth."
        )
    if not a5_hexes_expand:
        return None
    rows = []
    for a5_hex_expand in tqdm(
        a5_hexes_expand,
        desc="Building A5 expand",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = a52geo(
                a5_hex_expand, options, split_antimeridian=split_antimeridian
            )
            cell_resolution = a5.get_resolution(a5.hex_to_u64(a5_hex_expand))
            num_edges = 5  # A5 cells are pentagons
            if cell_resolution == 1:
                num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "a5", a5_hex_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    ouput_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            ouput_name = f"{base}_a5_expanded"
        else:
            ouput_name = "a5_expanded"
    return convert_to_output_format(out_gdf, output_format, ouput_name)


def a5expand_cli():
    """
    Command-line interface for a5expand with flexible input/output.
    """
    parser = argparse.ArgumentParser(description="A5 Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input A5 (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-r",
        "--resolution",
        type=int,
        help="Target A5 resolution to expand to (must be >= maximum input resolution). "
        "Ignores --depth.",
    )
    mode.add_argument(
        "-d",
        "--depth",
        type=int,
        help="Expand each cell by this many child levels (1 = direct children, "
        "2 = grandchildren, ...; 1 <= depth <= A5 max_res). "
        "Mixed input resolutions are allowed. Ignores --resolution.",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="A5 Hex field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
    )
    parser.add_argument(
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Enable Antimeridian splitting",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to a52geo. "
        "Example: '{\"segments\": 1000}'",
    )
    add_verbose_argument(parser)
    args = parser.parse_args()
    input_data = args.input
    resolution = args.resolution
    cellid = args.cellid
    output_format = args.output_format
    split_antimeridian = args.split_antimeridian

    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return

    result = a5expand(
        input_data,
        resolution=resolution,
        a5_hex=cellid,
        output_format=output_format,
        options=options,
        split_antimeridian=split_antimeridian,
        depth=args.depth,
        verbose=args.verbose,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)
