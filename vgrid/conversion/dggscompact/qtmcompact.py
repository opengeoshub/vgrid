"""
QTM Compact Module

This module provides functionality to compact and expand QTM cells with flexible input and output formats.

Key Functions:
    qtm_compact: Compact a list of QTM IDs with an optional parent depth
    qtmcompact: Compact a set of QTM cells to their covering set
    qtmexpand: Expand (uncompact) QTM cells to a target resolution or by depth
    qtmcompact_cli: Command-line interface for compaction
    qtmexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from tqdm import tqdm

from vgrid.conversion.dggs2geo.qtm2geo import qtm2geo
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
from vgrid.dggs import qtm


# --- QTM Compaction/Expansion Logic ---
def get_qtm_resolution(qtm_id):
    """Get the resolution of a QTM cell ID."""
    try:
        return len(qtm_id)
    except Exception as e:
        raise ValueError(f"Invalid QTM ID <{qtm_id}> : {e}")


def qtm_compact(qtm_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of QTM cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    qtm_ids : list of str
        QTM cell IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted QTM cell IDs.
    """
    depth = validate_dggs_compact_depth("qtm", depth)

    def parent_fn(qtm_id):
        parent = qtm.qtm_parent(qtm_id)
        if not parent or parent == qtm_id:
            return None
        return parent

    def children_fn(parent):
        return qtm.qtm_children(parent, len(parent) + 1)

    return compact_cells(
        qtm_ids,
        parent_fn,
        children_fn,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting QTM",
    )


def qtmcompact(
    input_data,
    qtm_id="qtm",
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    verbose=True,
):
    """
    Compact QTM cells to their covering set at a given parent depth.

    Compacts a set of QTM cells by replacing complete sets of children with
    their parent cells. ``depth`` limits how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg`` (same options as ``h3bin``). If ``agg`` is
    ``"count"``, ``numeric_col`` is ignored and the output ``count`` is the
    number of original input cells in each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing QTM cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of QTM cell IDs
    qtm_id : str, default "qtm"
        Name of the column containing QTM cell IDs.
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
        The compacted QTM cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = qtmcompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = qtmcompact(["A0", "A1", "A2", "A3"])

    >>> # Compact only one parent level
    >>> result = qtmcompact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = qtmcompact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = qtmcompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not qtm_id:
        qtm_id = "qtm"

    bags, agg_col = prepare_compact_bags(
        input_data,
        qtm_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="QTM cells",
    )
    if bags is None:
        print(f"No QTM IDs found in <{qtm_id}> field.")
        return

    qtm_ids_compact = qtm_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not qtm_ids_compact:
        return None

    rows = []
    for qtm_id_compact in tqdm(
        qtm_ids_compact,
        desc="Building QTM compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = qtm2geo(qtm_id_compact)
            cell_resolution = get_qtm_resolution(qtm_id_compact)
            num_edges = 3  # QTM cells are triangular
            row = geodesic_dggs_to_geoseries(
                "qtm", qtm_id_compact, cell_resolution, cell_polygon, num_edges
            )
            row[agg_col] = aggregate_values(bags.get(qtm_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_qtm_compacted"
        else:
            output_name = "qtm_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def qtmcompact_cli():
    """Command-line interface for QTM compaction."""
    parser = argparse.ArgumentParser(description="QTM Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input QTM (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="QTM ID field")
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

    result = qtmcompact(
        input_data,
        qtm_id=cellid,
        output_format=output_format,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def qtm_expand(qtm_ids, resolution=None, depth=None, verbose=True):
    """
    Expand QTM cell IDs to a target resolution, or by a relative child depth.

    When ``resolution`` is set, ``depth`` is ignored and all cells are expanded
    to that absolute resolution. When only ``depth`` is set, ``resolution`` is
    ignored and each cell is expanded ``depth`` levels down (``1`` = direct
    children, ``2`` = grandchildren, and so on).
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("qtm", resolution)
        expand_cells = []
        for qtm_id in tqdm(qtm_ids, desc="Expanding QTM", unit=" cells", disable=not verbose):
            cell_resolution = len(qtm_id)
            if cell_resolution >= resolution:
                expand_cells.append(qtm_id)
            else:
                expand_cells.extend(qtm.qtm_children(qtm_id, resolution))
        return expand_cells

    if depth is None:
        raise ValueError("Either resolution or depth must be specified.")
    depth = validate_dggs_expand_depth("qtm", depth)
    expand_cells = []
    for qtm_id in tqdm(qtm_ids, desc="Expanding QTM", unit=" cells", disable=not verbose):
        try:
            expand_cells.extend(qtm.qtm_children(qtm_id, len(qtm_id) + depth))
        except Exception:
            continue
    return expand_cells


def qtmexpand(
    input_data,
    resolution=None,
    qtm_id="qtm",
    output_format="gpd",
    verbose=True,
    depth=None,
):
    """
    Expand (uncompact) QTM cells to a target resolution or by a relative depth.

    When ``resolution`` is set, ``depth`` is ignored and cells are expanded to
    that absolute resolution (must be >= the maximum input resolution). When
    only ``depth`` is set, ``resolution`` is ignored: mixed-resolution input is
    allowed and each cell is expanded to its descendants ``depth`` levels down.
    """
    if qtm_id is None:
        qtm_id = "qtm"
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("qtm", resolution)
    elif depth is not None:
        depth = validate_dggs_expand_depth("qtm", depth)
    else:
        raise ValueError("Either resolution or depth must be specified.")

    gdf = process_input_data_compact(input_data, qtm_id)
    qtm_ids = gdf[qtm_id].drop_duplicates().tolist()

    if not qtm_ids:
        print(f"No QTM IDs found in <{qtm_id}> field.")
        return

    try:
        if resolution is not None:
            max_res = max(len(qid) for qid in qtm_ids)
            if resolution < max_res:
                print(f"Target expand resolution ({resolution}) must >= {max_res}.")
                return None
            qtm_ids_expand = qtm_expand(qtm_ids, resolution=resolution, verbose=verbose)
        else:
            qtm_ids_expand = qtm_expand(qtm_ids, depth=depth, verbose=verbose)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your QTM ID field, resolution, or depth."
        )

    if not qtm_ids_expand:
        return None

    rows = []
    for qtm_id_expand in tqdm(
        qtm_ids_expand,
        desc="Building QTM expand",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = qtm2geo(qtm_id_expand)
            cell_resolution = len(qtm_id_expand)
            num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "qtm", qtm_id_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_qtm_expanded"
        else:
            output_name = "qtm_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def qtmexpand_cli():
    """Command-line interface for QTM expansion."""
    parser = argparse.ArgumentParser(description="QTM Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input QTM (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-r",
        "--resolution",
        type=int,
        help="Target QTM resolution to expand to (must be >= maximum input resolution). "
        "Ignores --depth.",
    )
    mode.add_argument(
        "-d",
        "--depth",
        type=int,
        help="Expand each cell by this many child levels (1 = direct children, "
        "2 = grandchildren, ...; 1 <= depth <= QTM max_res). "
        "Mixed input resolutions are allowed. Ignores --resolution.",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="QTM ID field")
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
    result = qtmexpand(
        args.input,
        resolution=args.resolution,
        qtm_id=args.cellid,
        output_format=args.output_format,
        depth=args.depth,
        verbose=args.verbose,
    )

    if args.output_format in STRUCTURED_FORMATS:
        print(result)
