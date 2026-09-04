"""
EASE Compact Module

This module provides functionality to compact and expand EASE cells with flexible input and output formats.

Key Functions:
    ease_compact: Compact a list of EASE IDs with an optional parent depth
    easecompact: Compact a set of EASE cells to their covering set
    easeexpand: Expand (uncompact) EASE cells to a target resolution or by depth
    easecompact_cli: Command-line interface for compaction
    easeexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
import re
from tqdm import tqdm
from vgrid.conversion.dggs2geo.ease2geo import ease2geo
from vgrid.utils.geometry import geodesic_dggs_to_geoseries, get_ease_resolution
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
from ease_dggs.dggs.hierarchy import _parent_to_children


def _ease_parent(ease_id):
    match = re.match(r"L(\d+)\.(.+)", ease_id)
    if not match:
        return None
    resolution = int(match.group(1))
    if resolution == 0:
        return None
    return f"L{resolution - 1}." + ".".join(match.group(2).split(".")[:-1])


def _ease_children(parent):
    match = re.match(r"L(\d+)\..+", parent)
    resolution = int(match.group(1))
    return set(_parent_to_children(parent, resolution + 1))


def ease_compact(ease_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of EASE cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    ease_ids : list of str
        List of EASE cell IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted EASE cell IDs.

    Examples
    --------
    >>> ease_ids = ["L4.165767.02.02.20.71", "L4.165767.02.02.20.72"]
    >>> compacted = ease_compact(ease_ids)
    >>> print(f"Compacted {len(ease_ids)} cells to {len(compacted)} cells")
    """
    depth = validate_dggs_compact_depth("ease", depth)
    return compact_cells(
        ease_ids,
        _ease_parent,
        _ease_children,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting EASE",
    )


def easecompact(
    input_data,
    ease_id=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    verbose=True,
):
    """
    Compact EASE cells to their covering set at a given parent depth.

    Compacts a set of EASE cells by replacing complete sets of children with
    their parent cells. Mixed input resolutions are allowed and ``depth`` limits
    how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg``. If ``agg`` is ``"count"``, ``numeric_col`` is
    ignored and the output ``count`` is the number of original input cells in
    each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing EASE cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of EASE cell IDs
    ease_id : str, optional
        Name of the column containing EASE cell IDs. Defaults to "ease".
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
        The compacted EASE cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = easecompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = easecompact(["L4.165767.02.02.20.71", "L4.165767.02.02.20.72"])

    >>> # Compact only one parent level
    >>> result = easecompact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = easecompact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = easecompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not ease_id:
        ease_id = "ease"

    bags, agg_col = prepare_compact_bags(
        input_data,
        ease_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="EASE cells",
    )
    if bags is None:
        print(f"No EASE IDs found in <{ease_id}> field.")
        return

    ease_ids_compact = ease_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not ease_ids_compact:
        return None

    rows = []
    for ease_id_compact in tqdm(
        ease_ids_compact,
        desc="Building EASE compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = ease2geo(ease_id_compact)
            cell_resolution = get_ease_resolution(ease_id_compact)
            num_edges = 4  # EASE cells are rectangular
            row = geodesic_dggs_to_geoseries(
                "ease", ease_id_compact, cell_resolution, cell_polygon, num_edges
            )
            row[agg_col] = aggregate_values(bags.get(ease_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_ease_compacted"
        else:
            output_name = "ease_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def easecompact_cli():
    """Command-line interface for EASE compaction."""
    parser = argparse.ArgumentParser(description="EASE Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input EASE (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="EASE ID field")
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

    result = easecompact(
        input_data,
        ease_id=cellid,
        output_format=output_format,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def ease_expand(ease_ids, resolution=None, depth=None, verbose=True):
    """
    Expand EASE cell IDs to a target resolution, or by a relative child depth.

    When ``resolution`` is set, ``depth`` is ignored. When only ``depth`` is
    set, each cell is expanded ``depth`` levels down (``1`` = direct children,
    ``2`` = grandchildren, and so on).
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("ease", resolution)
        uncompacted_cells = []
        for ease_id in tqdm(ease_ids, desc="Expanding EASE", unit=" cells", disable=not verbose):
            ease_resolution = int(ease_id[1])
            if ease_resolution >= resolution:
                uncompacted_cells.append(ease_id)
            else:
                uncompacted_cells.extend(
                    _parent_to_children(ease_id, ease_resolution + 1)
                )
        return uncompacted_cells

    if depth is None:
        raise ValueError("Either resolution or depth must be specified.")
    depth = validate_dggs_expand_depth("ease", depth)
    cells = list(ease_ids)
    for _ in range(depth):
        nxt = []
        for ease_id in tqdm(cells, desc="Expanding EASE", unit=" cells", disable=not verbose):
            try:
                match = re.match(r"L(\d+)\..+", ease_id)
                res = int(match.group(1))
                nxt.extend(_parent_to_children(ease_id, res + 1))
            except Exception:
                continue
        cells = nxt
    return cells


def easeexpand(
    input_data,
    resolution=None,
    ease_id=None,
    output_format="gpd",
    verbose=True,
    depth=None,
):
    """
    Expand (uncompact) EASE cells to a target resolution or by a relative depth.

    When ``resolution`` is set, ``depth`` is ignored and cells are expanded to
    that absolute resolution (must be >= the maximum input resolution). When
    only ``depth`` is set, ``resolution`` is ignored: mixed-resolution input is
    allowed and each cell is expanded to its descendants ``depth`` levels down.
    """
    if ease_id is None:
        ease_id = "ease"
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("ease", resolution)
    elif depth is not None:
        depth = validate_dggs_expand_depth("ease", depth)
    else:
        raise ValueError("Either resolution or depth must be specified.")

    gdf = process_input_data_compact(input_data, ease_id)
    ease_ids = gdf[ease_id].drop_duplicates().tolist()

    if not ease_ids:
        print(f"No EASE IDs found in <{ease_id}> field.")
        return

    try:
        if resolution is not None:
            max_res = max(int(eid[1]) for eid in ease_ids)
            if resolution < max_res:
                print(f"Target expand resolution ({resolution}) must >= {max_res}.")
                return None
            ease_ids_expand = ease_expand(ease_ids, resolution=resolution, verbose=verbose)
        else:
            ease_ids_expand = ease_expand(ease_ids, depth=depth, verbose=verbose)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your EASE ID field, resolution, or depth."
        )

    if not ease_ids_expand:
        return None

    rows = []
    for ease_id_expand in tqdm(
        ease_ids_expand,
        desc="Building EASE expand",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = ease2geo(ease_id_expand)
            cell_resolution = get_ease_resolution(ease_id_expand)
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "ease", ease_id_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_ease_expanded"
        else:
            output_name = "ease_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def easeexpand_cli():
    """Command-line interface for EASE expansion."""
    parser = argparse.ArgumentParser(description="EASE Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input EASE (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-r",
        "--resolution",
        type=int,
        help="Target EASE resolution to expand to (must be >= maximum input resolution). "
        "Ignores --depth.",
    )
    mode.add_argument(
        "-d",
        "--depth",
        type=int,
        help="Expand each cell by this many child levels (1 = direct children, "
        "2 = grandchildren, ...; 1 <= depth <= EASE max_res). "
        "Mixed input resolutions are allowed. Ignores --resolution.",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="EASE ID field")
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
    result = easeexpand(
        args.input,
        resolution=args.resolution,
        ease_id=args.cellid,
        output_format=args.output_format,
        depth=args.depth,
        verbose=args.verbose,
    )

    if args.output_format in STRUCTURED_FORMATS:
        print(result)
