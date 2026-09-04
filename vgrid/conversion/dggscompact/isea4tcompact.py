"""
ISEA4T Compact Module

This module provides functionality to compact and expand ISEA4T cells with flexible input and output formats.

Key Functions:
    isea4tcompact: Compact a set of ISEA4T cells to their minimal covering set
    isea4texpand: Expand (uncompact) ISEA4T cells to a target resolution or by depth
    isea4tcompact_cli: Command-line interface for compaction
    isea4texpand_cli: Command-line interface for expansion

Note: This module is only supported on Windows systems due to OpenEaggr dependency.
"""

import os
import argparse
import geopandas as gpd
import platform

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.dggs.eaggr.enums.model import Model

    isea4t_dggs = Eaggr(Model.ISEA4T)

from tqdm import tqdm
from vgrid.conversion.dggs2geo.isea4t2geo import isea4t2geo
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


# --- ISEA4T Compaction/Expansion Logic ---
def get_isea4t_resolution(isea4t_id):
    return len(isea4t_id) - 2


def get_isea4t_cell_children(isea4t_cell, resolution):
    """Recursively expands a DGGS cell until all children reach the desired resolution."""
    cell_id = isea4t_cell.get_cell_id()
    cell_resolution = len(cell_id) - 2

    if cell_resolution >= resolution:
        return [
            isea4t_cell
        ]  # Base case: return the cell if it meets/exceeds resolution

    expanded_cells = []
    children = isea4t_dggs.get_dggs_cell_children(isea4t_cell)

    for child in children:
        expanded_cells.extend(get_isea4t_cell_children(child, resolution))

    return expanded_cells


def isea4t_compact(isea4t_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of ISEA4T cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    isea4t_ids : list of str
        ISEA4T cell IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted ISEA4T cell IDs.
    """
    depth = validate_dggs_compact_depth("isea4t", depth)

    def parent_fn(cell_id):
        if len(cell_id) > 2:
            return cell_id[:-1]
        return None

    def children_fn(parent):
        return {
            child.get_cell_id()
            for child in isea4t_dggs.get_dggs_cell_children(DggsCell(parent))
        }

    return compact_cells(
        isea4t_ids,
        parent_fn,
        children_fn,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting ISEA4T",
    )


def isea4t_expand(isea4t_ids, resolution=None, depth=None, verbose=True):
    """
    Expand ISEA4T cells to a target resolution, or by a relative child depth.

    When ``resolution`` is set, ``depth`` is ignored and all cells are expanded
    to that absolute resolution. When only ``depth`` is set, ``resolution`` is
    ignored and each cell is expanded ``depth`` levels down (``1`` = direct
    children, ``2`` = grandchildren, and so on).

    Returns cell objects (callers typically map ``.get_cell_id()``).
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("isea4t", resolution)
        expand_cells = []
        for isea4t_id in tqdm(isea4t_ids, desc="Expanding ISEA4T", unit=" cells", disable=not verbose):
            isea4t_cell = DggsCell(isea4t_id)
            expand_cells.extend(get_isea4t_cell_children(isea4t_cell, resolution))
        return expand_cells

    if depth is None:
        raise ValueError("Either resolution or depth must be specified.")
    depth = validate_dggs_expand_depth("isea4t", depth)
    expand_cells = []
    for isea4t_id in tqdm(isea4t_ids, desc="Expanding ISEA4T", unit=" cells", disable=not verbose):
        try:
            current = get_isea4t_resolution(isea4t_id)
            expand_cells.extend(
                get_isea4t_cell_children(DggsCell(isea4t_id), current + depth)
            )
        except Exception:
            continue
    return expand_cells


def isea4tcompact(
    input_data,
    isea4t_id=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
):
    """
    Compact ISEA4T cells to their covering set at a given parent depth.

    Compacts a set of ISEA4T cells by replacing complete sets of children with their
    parent cells. Mixed input resolutions are allowed and ``depth`` limits how far
    up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg``. If ``agg`` is ``"count"``, ``numeric_col`` is
    ignored and the output ``count`` is the number of original input cells in
    each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing ISEA4T cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of ISEA4T cell IDs
    isea4t_id : str, optional
        Name of the column containing ISEA4T cell IDs. Defaults to "isea4t".
    depth : int, default -1
        Compaction depth: ``0`` leaves cells unchanged, ``-1`` compact as far as
        possible, ``1`` merges to the direct parent, ``2`` to the grandparent, etc.
    agg : str, default "count"
        Aggregation applied to original child values when cells compact into a
        parent (``count``, ``min``, ``max``, ``sum``, ``mean``, ``median``,
        ``std``, ``var``, ``range``, ``minority``, ``majority``, ``variety``).
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
    fix_antimeridian : str, optional
        Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        Defaults to None when omitted.
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The compacted ISEA4T cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = isea4tcompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = isea4tcompact(["A0", "A1", "A2", "A3"])

    >>> # Compact only one parent level
    >>> result = isea4tcompact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = isea4tcompact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = isea4tcompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not isea4t_id:
        isea4t_id = "isea4t"
    bags, agg_col = prepare_compact_bags(
        input_data,
        isea4t_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="ISEA4T cells",
    )
    if bags is None:
        print(f"No ISEA4T isea4t_ids found in <{isea4t_id}> field.")
        return
    isea4t_ids_compact = isea4t_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not isea4t_ids_compact:
        return None
    rows = []
    for isea4t_id_compact in tqdm(
        isea4t_ids_compact,
        desc="Building ISEA4T compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = isea4t2geo(
                isea4t_id_compact, fix_antimeridian=fix_antimeridian
            )
            cell_resolution = get_isea4t_resolution(isea4t_id_compact)
            num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "isea4t", isea4t_id_compact, cell_resolution, cell_polygon, num_edges
            )
            row[agg_col] = aggregate_values(bags.get(isea4t_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    ouput_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            ouput_name = f"{base}_isea4t_compacted"
        else:
            ouput_name = "isea4t_compacted"
    return convert_to_output_format(out_gdf, output_format, ouput_name)


def isea4tcompact_cli():
    parser = argparse.ArgumentParser(description="ISEA4T Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input ISEA4T (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="ISEA4T ID field")
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
    result = isea4tcompact(
        input_data,
        isea4t_id=cellid,
        output_format=output_format,
        fix_antimeridian=args.fix_antimeridian,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


def isea4texpand(
    input_data,
    resolution=None,
    isea4t_id=None,
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
    depth=None,
):
    """
    Expand (uncompact) ISEA4T cells to a target resolution or by a relative depth.

    When ``resolution`` is set, ``depth`` is ignored and cells are expanded to
    that absolute resolution (must be >= the maximum input resolution). When
    only ``depth`` is set, ``resolution`` is ignored: mixed-resolution input is
    allowed and each cell is expanded to its descendants ``depth`` levels down.
    """
    if isea4t_id is None:
        isea4t_id = "isea4t"
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("isea4t", resolution)
    elif depth is not None:
        depth = validate_dggs_expand_depth("isea4t", depth)
    else:
        raise ValueError("Either resolution or depth must be specified.")
    gdf = process_input_data_compact(input_data, isea4t_id)
    isea4t_ids = gdf[isea4t_id].drop_duplicates().tolist()
    if not isea4t_ids:
        print(f"No ISEA4T IDs found in <{isea4t_id}> field.")
        return
    try:
        if resolution is not None:
            max_res = max(get_isea4t_resolution(cid) for cid in isea4t_ids)
            if resolution < max_res:
                print(f"Target expand resolution ({resolution}) must >= {max_res}.")
                return None
            isea4t_cells_expand = isea4t_expand(isea4t_ids, resolution=resolution, verbose=verbose)
        else:
            isea4t_cells_expand = isea4t_expand(isea4t_ids, depth=depth, verbose=verbose)
        isea4t_ids_expand = [c.get_cell_id() for c in isea4t_cells_expand]
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your ISEA4T ID field, resolution, or depth."
        )
    if not isea4t_ids_expand:
        return None
    rows = []
    for isea4t_id_expand in tqdm(
        isea4t_ids_expand,
        desc="Building ISEA4T expand",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = isea4t2geo(
                isea4t_id_expand, fix_antimeridian=fix_antimeridian
            )
            cell_resolution = get_isea4t_resolution(isea4t_id_expand)
            num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "isea4t", isea4t_id_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    ouput_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            ouput_name = f"{base}_isea4t_expanded"
        else:
            ouput_name = "isea4t_expanded"
    return convert_to_output_format(out_gdf, output_format, ouput_name)


def isea4texpand_cli():
    parser = argparse.ArgumentParser(description="ISEA4T Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input ISEA4T (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-r",
        "--resolution",
        type=int,
        help="Target ISEA4T resolution to expand to (must be >= maximum input resolution). "
        "Ignores --depth.",
    )
    mode.add_argument(
        "-d",
        "--depth",
        type=int,
        help="Expand each cell by this many child levels (1 = direct children, "
        "2 = grandchildren, ...; 1 <= depth <= ISEA4T max_res). "
        "Mixed input resolutions are allowed. Ignores --resolution.",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="ISEA4T ID field")
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
    if platform.system() == "Windows":
        result = isea4texpand(
            args.input,
            resolution=args.resolution,
            isea4t_id=args.cellid,
            output_format=args.output_format,
            fix_antimeridian=args.fix_antimeridian,
            depth=args.depth,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    else:
        print("ISEA4T is only supported on Windows systems")
