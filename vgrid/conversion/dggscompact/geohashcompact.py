"""
Geohash Compact Module

This module provides functionality to compact and expand Geohash cells with flexible input and output formats.

Key Functions:
    geohash_compact: Compact a list of Geohash IDs with an optional parent depth
    geohashcompact: Compact a set of Geohash cells to their covering set
    geohashexpand: Expand (uncompact) Geohash cells to a target resolution or by depth
    geohashcompact_cli: Command-line interface for compaction
    geohashexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from tqdm import tqdm

from vgrid.conversion.dggs2geo.geohash2geo import geohash2geo
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
from vgrid.dggs import geohash


# --- Geohash Compaction/Expansion Logic ---
def get_geohash_resolution(geohash_id):
    """Get the resolution of a Geohash cell ID."""
    return len(geohash_id)


def geohash_compact(geohash_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of Geohash cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    geohash_ids : list of str
        Geohash cell IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted Geohash cell IDs.
    """
    depth = validate_dggs_compact_depth("geohash", depth)

    def parent_fn(geohash_id):
        if len(geohash_id) <= 1:
            return None
        return geohash.geohash_parent(geohash_id)

    def children_fn(parent):
        return geohash.geohash_children(parent, len(parent) + 1)

    return compact_cells(
        geohash_ids,
        parent_fn,
        children_fn,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting Geohash",
    )


def geohashcompact(
    input_data,
    geohash_id=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    verbose=True,
):
    """
    Compact Geohash cells to their covering set at a given parent depth.

    Compacts a set of Geohash cells by replacing complete sets of children with
    their parent cells. ``depth`` limits how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg`` (same options as ``h3bin``). If ``agg`` is
    ``"count"``, ``numeric_col`` is ignored and the output ``count`` is the
    number of original input cells in each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing Geohash cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of Geohash cell IDs
    geohash_id : str, optional
        Name of the column containing Geohash cell IDs. Defaults to "geohash".
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
        The compacted Geohash cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = geohashcompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = geohashcompact(["w3gvk1td8", "w3gvk1td9"])

    >>> # Compact only one parent level
    >>> result = geohashcompact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = geohashcompact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = geohashcompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not geohash_id:
        geohash_id = "geohash"

    bags, agg_col = prepare_compact_bags(
        input_data,
        geohash_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="Geohash cells",
    )
    if bags is None:
        print(f"No Geohash IDs found in <{geohash_id}> field.")
        return

    geohash_ids_compact = geohash_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not geohash_ids_compact:
        return None

    rows = []
    for geohash_id_compact in tqdm(
        geohash_ids_compact,
        desc="Building Geohash compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = geohash2geo(geohash_id_compact)
            cell_resolution = get_geohash_resolution(geohash_id_compact)
            row = graticule_dggs_to_geoseries(
                "geohash", geohash_id_compact, cell_resolution, cell_polygon
            )
            row[agg_col] = aggregate_values(bags.get(geohash_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_geohash_compacted"
        else:
            output_name = "geohash_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def geohashcompact_cli():
    """Command-line interface for Geohash compaction."""
    parser = argparse.ArgumentParser(description="Geohash Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input Geohash (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="Geohash ID field")
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

    result = geohashcompact(
        input_data,
        geohash_id=cellid,
        output_format=output_format,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def geohash_expand(geohash_ids, resolution=None, depth=None, verbose=True):
    """
    Expand Geohash cell IDs to a target resolution, or by a relative child depth.

    When ``resolution`` is set, ``depth`` is ignored and all cells are expanded
    to that absolute resolution. When only ``depth`` is set, ``resolution`` is
    ignored and each cell is expanded ``depth`` levels down (``1`` = direct
    children, ``2`` = grandchildren, and so on).
    """
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("geohash", resolution)
        expand_cells = []
        for geohash_id in tqdm(geohash_ids, desc="Expanding Geohash", unit=" cells", disable=not verbose):
            cell_resolution = len(geohash_id)
            if cell_resolution >= resolution:
                expand_cells.append(geohash_id)
            else:
                expand_cells.extend(
                    geohash.geohash_children(geohash_id, resolution)
                )
        return expand_cells

    if depth is None:
        raise ValueError("Either resolution or depth must be specified.")
    depth = validate_dggs_expand_depth("geohash", depth)
    expand_cells = []
    for geohash_id in tqdm(geohash_ids, desc="Expanding Geohash", unit=" cells", disable=not verbose):
        try:
            expand_cells.extend(
                geohash.geohash_children(geohash_id, len(geohash_id) + depth)
            )
        except Exception:
            continue
    return expand_cells


def geohashexpand(
    input_data,
    resolution=None,
    geohash_id=None,
    output_format="gpd",
    verbose=True,
    depth=None,
):
    """
    Expand (uncompact) Geohash cells to a target resolution or by a relative depth.

    When ``resolution`` is set, ``depth`` is ignored and cells are expanded to
    that absolute resolution (must be >= the maximum input resolution). When
    only ``depth`` is set, ``resolution`` is ignored: mixed-resolution input is
    allowed and each cell is expanded to its descendants ``depth`` levels down.
    """
    if geohash_id is None:
        geohash_id = "geohash"
    if resolution is not None:
        resolution = validate_dggs_expand_resolution("geohash", resolution)
    elif depth is not None:
        depth = validate_dggs_expand_depth("geohash", depth)
    else:
        raise ValueError("Either resolution or depth must be specified.")

    gdf = process_input_data_compact(input_data, geohash_id)
    geohash_ids = gdf[geohash_id].drop_duplicates().tolist()

    if not geohash_ids:
        print(f"No Geohash IDs found in <{geohash_id}> field.")
        return

    try:
        if resolution is not None:
            max_res = max(len(gid) for gid in geohash_ids)
            if resolution < max_res:
                print(f"Target expand resolution ({resolution}) must >= {max_res}.")
                return None
            geohash_ids_expand = geohash_expand(geohash_ids, resolution=resolution, verbose=verbose)
        else:
            geohash_ids_expand = geohash_expand(geohash_ids, depth=depth, verbose=verbose)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your Geohash ID field, resolution, or depth."
        )

    if not geohash_ids_expand:
        return None

    rows = []
    for geohash_id_expand in tqdm(
        geohash_ids_expand,
        desc="Building Geohash expand",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = geohash2geo(geohash_id_expand)
            cell_resolution = len(geohash_id_expand)
            row = graticule_dggs_to_geoseries(
                "geohash", geohash_id_expand, cell_resolution, cell_polygon
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_geohash_expanded"
        else:
            output_name = "geohash_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def geohashexpand_cli():
    """Command-line interface for Geohash expansion."""
    parser = argparse.ArgumentParser(description="Geohash Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input Geohash (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-r",
        "--resolution",
        type=int,
        help="Target Geohash resolution to expand to (must be >= maximum input resolution). "
        "Ignores --depth.",
    )
    mode.add_argument(
        "-d",
        "--depth",
        type=int,
        help="Expand each cell by this many child levels (1 = direct children, "
        "2 = grandchildren, ...; 1 <= depth <= Geohash max_res). "
        "Mixed input resolutions are allowed. Ignores --resolution.",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="Geohash ID field")
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
    result = geohashexpand(
        args.input,
        resolution=args.resolution,
        geohash_id=args.cellid,
        output_format=args.output_format,
        depth=args.depth,
        verbose=args.verbose,
    )

    if args.output_format in STRUCTURED_FORMATS:
        print(result)
