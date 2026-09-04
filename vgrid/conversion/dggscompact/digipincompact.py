"""
Digipin Compact Module

This module provides functionality to compact and expand DIGIPIN cells with flexible input and output formats.

Key Functions:
    digipin_compact: Compact a list of DIGIPIN IDs with an optional parent depth
    digipincompact: Compact a set of DIGIPIN cells to their covering set
    digipinexpand: Expand (uncompact) DIGIPIN cells to a target resolution or by depth
    digipincompact_cli: Command-line interface for compaction
    digipinexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from tqdm import tqdm
from vgrid.utils.geometry import graticule_dggs_to_geoseries
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
from vgrid.dggs.digipin import digipin_parent, digipin_children, digipin_resolution
from vgrid.conversion.dggs2geo.digipin2geo import digipin2geo


def digipin_compact(digipin_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of DIGIPIN cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    digipin_ids : list of str
        DIGIPIN cell IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted DIGIPIN cell IDs.
    """
    depth = validate_dggs_compact_depth("digipin", depth)

    def parent_fn(digipin_id):
        parent = digipin_parent(digipin_id)
        if parent == "Invalid DIGIPIN":
            return None
        return parent

    def children_fn(parent):
        parent_resolution = digipin_resolution(parent)
        if isinstance(parent_resolution, str):
            raise ValueError("Invalid DIGIPIN resolution")
        return digipin_children(parent, parent_resolution + 1)

    return compact_cells(
        digipin_ids,
        parent_fn,
        children_fn,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting DIGIPIN",
    )


def digipincompact(
    input_data,
    digipin_id="digipin",
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    verbose=True,
):
    """
    Compact DIGIPIN cells to their covering set at a given parent depth.

    Compacts a set of DIGIPIN cells by replacing complete sets of children with
    their parent cells. ``depth`` limits how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg`` (same options as ``h3bin``). If ``agg`` is
    ``"count"``, ``numeric_col`` is ignored and the output ``count`` is the
    number of original input cells in each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing DIGIPIN cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of DIGIPIN cell IDs
    digipin_id : str, default "digipin"
        Name of the column containing DIGIPIN cell IDs.
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
        The compacted DIGIPIN cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = digipincompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = digipincompact(["F3K-F", "F3K-C", "F3K-9", "F3K-8"])

    >>> # Compact only one parent level
    >>> result = digipincompact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = digipincompact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = digipincompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not digipin_id:
        digipin_id = "digipin"

    bags, agg_col = prepare_compact_bags(
        input_data,
        digipin_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="DIGIPIN cells",
    )
    if bags is None:
        print(f"No DIGIPIN IDs found in <{digipin_id}> field.")
        return

    digipin_ids_compact = digipin_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not digipin_ids_compact:
        return None

    rows = []
    for digipin_id_compact in tqdm(
        digipin_ids_compact,
        desc="Building DIGIPIN compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = digipin2geo(digipin_id_compact)
            cell_resolution = digipin_resolution(digipin_id_compact)
            if isinstance(cell_resolution, str):
                continue  # Skip invalid resolutions
            row = graticule_dggs_to_geoseries(
                "digipin", digipin_id_compact, cell_resolution, cell_polygon
            )
            row[agg_col] = aggregate_values(bags.get(digipin_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_digipin_compacted"
        else:
            output_name = "digipin_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def digipincompact_cli():
    """Command-line interface for DIGIPIN compaction."""
    parser = argparse.ArgumentParser(description="DIGIPIN Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input DIGIPIN (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="DIGIPIN ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
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

    result = digipincompact(
        input_data,
        digipin_id=cellid,
        output_format=output_format,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def digipin_expand(digipin_ids, resolution=None, depth=None, verbose=True):
    """
    Expand DIGIPIN cell IDs to a target resolution, or by a relative child depth.

    When ``resolution`` is set, ``depth`` is ignored and all cells are expanded
    to that absolute resolution. When only ``depth`` is set, ``resolution`` is
    ignored and each cell is expanded ``depth`` levels down (``1`` = direct
    children, ``2`` = grandchildren, and so on).
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("digipin", resolution)
        expand_cells = []
        for digipin_id in tqdm(digipin_ids, desc="Expanding DIGIPIN", unit=" cells", disable=not verbose):
            current_resolution = digipin_resolution(digipin_id)
            if isinstance(current_resolution, str):
                raise ValueError("Invalid DIGIPIN format.")
            if current_resolution >= resolution:
                expand_cells.append(digipin_id)
            else:
                expand_cells.extend(digipin_children(digipin_id, resolution))
        return expand_cells

    if depth is None:
        raise ValueError("Either resolution or depth must be specified.")
    depth = validate_dggs_expand_depth("digipin", depth)
    expand_cells = []
    for digipin_id in tqdm(digipin_ids, desc="Expanding DIGIPIN", unit=" cells", disable=not verbose):
        try:
            current_resolution = digipin_resolution(digipin_id)
            if isinstance(current_resolution, str):
                continue
            expand_cells.extend(
                digipin_children(digipin_id, current_resolution + depth)
            )
        except Exception:
            continue
    return expand_cells


def digipinexpand(
    input_data,
    resolution=None,
    digipin_id="digipin",
    output_format="gpd",
    verbose=True,
    depth=None,
):
    """
    Expand (uncompact) DIGIPIN cells to a target resolution or by a relative depth.

    When ``resolution`` is set, ``depth`` is ignored and cells are expanded to
    that absolute resolution (must be >= the maximum input resolution). When
    only ``depth`` is set, ``resolution`` is ignored: mixed-resolution input is
    allowed and each cell is expanded to its descendants ``depth`` levels down.
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("digipin", resolution)
    elif depth is not None:
        depth = validate_dggs_expand_depth("digipin", depth)
    else:
        raise ValueError("Either resolution or depth must be specified.")

    gdf = process_input_data_compact(input_data, digipin_id)
    digipin_ids = gdf[digipin_id].drop_duplicates().tolist()

    if not digipin_ids:
        print(f"No DIGIPIN IDs found in <{digipin_id}> field.")
        return

    try:
        if resolution is not None:
            max_res = max(digipin_resolution(tid) for tid in digipin_ids)
            if isinstance(max_res, str):
                raise ValueError("Invalid DIGIPIN format.")
            if resolution < max_res:
                print(f"Target expand resolution ({resolution}) must >= {max_res}.")
                return None
            digipin_ids_expand = digipin_expand(digipin_ids, resolution=resolution, verbose=verbose)
        else:
            digipin_ids_expand = digipin_expand(digipin_ids, depth=depth, verbose=verbose)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your DIGIPIN ID field, resolution, or depth."
        )

    if not digipin_ids_expand:
        return None

    rows = []
    for digipin_id_expand in tqdm(
        digipin_ids_expand,
        desc="Building DIGIPIN expand",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = digipin2geo(digipin_id_expand)
            cell_resolution = digipin_resolution(digipin_id_expand)
            row = graticule_dggs_to_geoseries(
                "digipin", digipin_id_expand, cell_resolution, cell_polygon
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_digipin_expanded"
        else:
            output_name = "digipin_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def digipinexpand_cli():
    """Command-line interface for DIGIPIN expansion."""
    parser = argparse.ArgumentParser(description="DIGIPIN Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input DIGIPIN (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-r",
        "--resolution",
        type=int,
        help="Target DIGIPIN resolution to expand to (must be >= maximum input resolution). "
        "Ignores --depth.",
    )
    mode.add_argument(
        "-d",
        "--depth",
        type=int,
        help="Expand each cell by this many child levels (1 = direct children, "
        "2 = grandchildren, ...; 1 <= depth <= DIGIPIN max_res). "
        "Mixed input resolutions are allowed. Ignores --resolution.",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="DIGIPIN ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
    )

    add_verbose_argument(parser)
    args = parser.parse_args()
    result = digipinexpand(
        args.input,
        resolution=args.resolution,
        digipin_id=args.cellid,
        output_format=args.output_format,
        depth=args.depth,
        verbose=args.verbose,
    )

    if args.output_format in STRUCTURED_FORMATS:
        print(result)
