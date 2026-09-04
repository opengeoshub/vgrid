"""
S2 Compact Module

This module provides functionality to compact and expand S2 cells with flexible input and output formats.

Key Functions:
    s2_compact: Compact a list of S2 tokens with an optional parent depth
    s2compact: Compact a set of S2 cells to their covering set
    s2expand: Expand (uncompact) S2 cells to a target resolution or by depth
    s2compact_cli: Command-line interface for compaction
    s2expand_cli: Command-line interface for expansion
"""

import os
import argparse

import geopandas as gpd
from tqdm import tqdm
from vgrid.dggs import s2
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
from vgrid.conversion.dggs2geo.s22geo import s22geo


def _s2_parent(token):
    cid = s2.CellId.from_token(str(token))
    if cid.level() <= 0:
        return None
    return cid.parent().to_token()


def _s2_children(token):
    cid = s2.CellId.from_token(str(token))
    return {c.to_token() for c in cid.children()}


def s2_compact(s2_tokens, depth=-1, bags=None, verbose=True):
    """
    Compact a list of S2 cell tokens by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    s2_tokens : list of str
        S2 cell tokens to compact. Mixed resolutions are allowed.
    depth : int, default -1
        How many parent levels to climb:
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
        Sorted compacted S2 cell tokens.
    """
    depth = validate_dggs_compact_depth("s2", depth)
    return compact_cells(
        s2_tokens,
        _s2_parent,
        _s2_children,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting S2",
    )


def s2compact(
    input_data,
    s2_token="s2",
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
):
    """
    Compact S2 cells to their covering set at a given parent depth.

    Compacts a set of S2 cells by replacing complete sets of children with their
    parent cells. Mixed input resolutions are allowed and ``depth`` limits how
    far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg``. If ``agg`` is ``"count"``, ``numeric_col`` is
    ignored and the output ``count`` is the number of original input cells in
    each compacted cell.

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
    depth : int, default -1
        Compaction depth: ``0`` leaves cells unchanged, ``-1`` compact as far as
        possible, ``1`` merges to the direct parent, ``2`` to the grandparent, etc.
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
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

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

    >>> # Compact only one parent level
    >>> result = s2compact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = s2compact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = s2compact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not s2_token:
        s2_token = "s2"
    bags, agg_col = prepare_compact_bags(
        input_data,
        s2_token,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="S2 cells",
    )
    if bags is None:
        print(f"No S2 tokens found in <{s2_token}> field.")
        return

    s2_tokens_compact = s2_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not s2_tokens_compact:
        return None

    rows = []
    for s2_token_compact in tqdm(
        s2_tokens_compact,
        desc="Building S2 compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = s22geo(s2_token_compact, fix_antimeridian=fix_antimeridian)
            cell_resolution = s2.CellId.from_token(s2_token_compact).level()
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "s2", s2_token_compact, cell_resolution, cell_polygon, num_edges
            )
            row[agg_col] = aggregate_values(bags.get(s2_token_compact, []), agg)
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

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
    input_data = args.input
    cellid = args.cellid
    output_format = args.output_format
    fix_antimeridian = args.fix_antimeridian
    result = s2compact(
        input_data,
        s2_token=cellid,
        output_format=output_format,
        fix_antimeridian=fix_antimeridian,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


def s2_expand(s2_tokens, resolution=None, depth=None, verbose=True):
    """
    Expand S2 cell tokens to a target resolution, or by a relative child depth.

    When ``resolution`` is set, ``depth`` is ignored and all cells are uncompacted
    to that absolute resolution. When only ``depth`` is set, ``resolution`` is
    ignored and each cell (at any resolution) is expanded to its descendants
    ``depth`` levels down (``1`` = direct children, ``2`` = grandchildren, and
    so on).
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("s2", resolution)
        s2_cells = list({s2.CellId.from_token(token) for token in s2_tokens})
        covering = s2.CellUnion(s2_cells, raw=False)
        expanded_cells = covering.denormalize(resolution, 1)
        return [cell_id.to_token() for cell_id in expanded_cells]

    if depth is None:
        raise ValueError("Either resolution or depth must be specified.")
    depth = validate_dggs_expand_depth("s2", depth)
    expanded = []
    for token in tqdm(s2_tokens, desc="Expanding S2", unit=" cells", disable=not verbose):
        try:
            cid = s2.CellId.from_token(str(token))
            expanded.extend(c.to_token() for c in cid.children(cid.level() + depth))
        except Exception:
            continue
    return expanded


def s2expand(
    input_data,
    resolution=None,
    s2_token="s2",
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
    depth=None,
):
    """
    Expand (uncompact) S2 cells to a target resolution or by a relative depth.

    When ``resolution`` is set, ``depth`` is ignored and cells are expanded to
    that absolute resolution (must be >= the maximum input resolution). When
    only ``depth`` is set, ``resolution`` is ignored: mixed-resolution input is
    allowed and each cell is expanded to its descendants ``depth`` levels down
    (``1`` = direct children, ``2`` = grandchildren, and so on).

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing S2 cell tokens. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of S2 cell tokens
    resolution : int, optional
        Target S2 resolution to expand the cells to. Must be >= maximum input
        resolution. When set, ``depth`` is ignored.
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
    depth : int, optional
        Relative expansion depth (``1 <= depth <= max_res``). Used when
        ``resolution`` is not set. Each input cell is expanded ``depth``
        levels: ``1`` = direct children, ``2`` = grandchildren, and so on.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The expanded S2 cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> result = s2expand("cells.geojson", resolution=10)
    >>> result = s2expand(["31752f45cc94"], resolution=10)
    >>> result = s2expand(cells, depth=1)
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("s2", resolution)
    elif depth is not None:
        depth = validate_dggs_expand_depth("s2", depth)
    else:
        raise ValueError("Either resolution or depth must be specified.")
    gdf = process_input_data_compact(input_data, s2_token)
    s2_tokens = gdf[s2_token].drop_duplicates().tolist()
    if not s2_tokens:
        print(f"No S2 tokens found in <{s2_token}> field.")
        return
    try:
        if resolution is not None:
            s2_cells = [s2.CellId.from_token(token) for token in s2_tokens]
            if not s2_cells:
                print(f"No valid S2 tokens found in <{s2_token}> field.")
                return
            max_res = max(cell.level() for cell in s2_cells)
            if resolution < max_res:
                print(f"Target expand resolution ({resolution}) must >= {max_res}.")
                return None
            s2_tokens_expand = s2_expand(s2_tokens, resolution=resolution, verbose=verbose)
        else:
            s2_tokens_expand = s2_expand(s2_tokens, depth=depth, verbose=verbose)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your S2 ID field, resolution, or depth."
        )
    if not s2_tokens_expand:
        return None
    rows = []
    for s2_token_expand in tqdm(
        s2_tokens_expand,
        desc="Building S2 expand",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = s22geo(s2_token_expand, fix_antimeridian=fix_antimeridian)
            cell_resolution = s2.CellId.from_token(s2_token_expand).level()
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "s2", s2_token_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

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
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-r",
        "--resolution",
        type=int,
        help="Target S2 resolution to expand to (must be >= maximum input resolution). "
        "Ignores --depth.",
    )
    mode.add_argument(
        "-d",
        "--depth",
        type=int,
        help="Expand each cell by this many child levels (1 = direct children, "
        "2 = grandchildren, ...; 1 <= depth <= S2 max_res). "
        "Mixed input resolutions are allowed. Ignores --resolution.",
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
    add_verbose_argument(parser)
    args = parser.parse_args()
    result = s2expand(
        args.input,
        resolution=args.resolution,
        s2_token=args.cellid,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        depth=args.depth,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)
