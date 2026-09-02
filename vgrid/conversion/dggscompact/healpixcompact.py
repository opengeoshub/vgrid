"""
HEALPix Compact Module

This module provides functionality to compact and expand HEALPix UNIQ cells with
flexible input and output formats.

Key Functions:
    healpix_compact: Compact a list of HEALPix UNIQ IDs with an optional parent depth
    healpixcompact: Compact a set of HEALPix cells to their covering set
    healpixexpand: Expand (uncompact) HEALPix cells to a target resolution
    healpixcompact_cli: Command-line interface for compaction
    healpixexpand_cli: Command-line interface for expansion
"""

import os
import argparse

import geopandas as gpd
from tqdm import tqdm

from vgrid.conversion.dggs2geo.healpix2geo import healpix2geo
from vgrid.dggs.healpix import (
    nestDescendants,
    orderpix2uniq,
    uniq2orderpix,
    uniqChildren,
)
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    aggregate_values,
    compact_cells,
    convert_to_output_format,
    prepare_compact_bags,
    process_input_data_compact,
    validate_healpix_resolution,
)
from vgrid.utils.constants import AGG_OPTIONS, DGGS_TYPES, OUTPUT_FORMATS, STRUCTURED_FORMATS

min_res = DGGS_TYPES["healpix"]["min_res"]
max_res = DGGS_TYPES["healpix"]["max_res"]


def _healpix_parent(uid):
    uid = int(uid)
    order, ipix = uniq2orderpix(uid)
    if order == 0:
        return None
    return orderpix2uniq(order - 1, ipix >> 2)


def _healpix_children(uid):
    return set(uniqChildren(int(uid)))


def healpix_compact(healpix_ids, depth=-1, bags=None, verbose=True):
    """
    Compact HEALPix UNIQ IDs by replacing complete sibling sets with their parent.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    healpix_ids : list of int or str
        HEALPix UNIQ cell IDs to compact. Mixed resolutions are allowed.
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
    list of int
        Sorted compacted HEALPix UNIQ cell IDs.
    """
    return compact_cells(
        [int(uid) for uid in healpix_ids],
        _healpix_parent,
        _healpix_children,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting HEALPix",
    )


def healpix_expand(healpix_ids, resolution):
    """
    Expand HEALPix UNIQ IDs to a target resolution.

    Parameters
    ----------
    healpix_ids : list of int or str
        HEALPix UNIQ cell IDs.
    resolution : int
        Target HEALPix resolution/order. Must be >= each input cell's order.

    Returns
    -------
    list of int
        Expanded HEALPix UNIQ cell IDs at the target resolution.
    """
    resolution = validate_healpix_resolution(resolution)
    expanded = []
    for uid in healpix_ids:
        order, ipix = uniq2orderpix(int(uid))
        if resolution < order:
            raise ValueError(
                f"Target resolution ({resolution}) must be >= cell order ({order})"
            )
        levels = resolution - order
        for child_ipix in nestDescendants(ipix, levels):
            expanded.append(orderpix2uniq(resolution, child_ipix))
    return expanded


def healpixcompact(
    input_data,
    healpix_id=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
):
    """
    Compact HEALPix cells to their covering set at a given parent depth.

    Compacts a set of HEALPix cells by replacing complete sets of children with
    their parent cells. Mixed input resolutions are allowed and ``depth`` limits
    how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg``. If ``agg`` is ``"count"``, ``numeric_col`` is
    ignored and the output ``count`` is the number of original input cells in
    each compacted cell.

    Parameters
    ----------
    input_data : str, dict, GeoDataFrame, or list
        Input containing HEALPix UNIQ cell IDs.
    healpix_id : str, optional
        Column name for HEALPix IDs. Defaults to ``"healpix"``.
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
        Output format.
    fix_antimeridian : str, optional
        Antimeridian fixing method.
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

    Returns
    -------
    GeoDataFrame, str, dict, or None
        Compacted HEALPix cells in the requested format.
    """
    if healpix_id is None:
        healpix_id = "healpix"
    bags, agg_col = prepare_compact_bags(
        input_data,
        healpix_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="HEALPix cells",
        id_cast=int,
    )
    if bags is None:
        print(f"No HEALPix IDs found in <{healpix_id}> field.")
        return

    healpix_ids_compact = healpix_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not healpix_ids_compact:
        return None

    rows = []
    for uniq_id in tqdm(
        healpix_ids_compact,
        desc="Building HEALPix compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = healpix2geo(uniq_id, fix_antimeridian=fix_antimeridian)
            order, _ipix = uniq2orderpix(int(uniq_id))
            row = geodesic_dggs_to_geoseries(
                "healpix", uniq_id, order, cell_polygon, num_edges=4
            )
            row[agg_col] = aggregate_values(bags.get(uniq_id, []), agg)
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_healpix_compacted"
        else:
            output_name = "healpix_compacted"
    return convert_to_output_format(out_gdf, output_format, output_name)


def healpixcompact_cli():
    """Command-line interface for healpixcompact."""
    parser = argparse.ArgumentParser(description="HEALPix Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input HEALPix (GeoJSON, Shapefile, CSV, Parquet, or GeoDataFrame)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="HEALPix ID field")
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
        help="Antimeridian fixing method",
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
    result = healpixcompact(
        args.input,
        healpix_id=args.cellid,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


def healpixexpand(
    input_data,
    resolution,
    healpix_id=None,
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Expand (uncompact) HEALPix cells to a target resolution.

    Parameters
    ----------
    input_data : str, dict, GeoDataFrame, or list
        Input containing HEALPix UNIQ cell IDs.
    resolution : int
        Target resolution/order. Must be >= maximum input cell order.
    healpix_id : str, optional
        Column name for HEALPix IDs. Defaults to ``"healpix"``.
    output_format : str, default "gpd"
        Output format.
    fix_antimeridian : str, optional
        Antimeridian fixing method.

    Returns
    -------
    GeoDataFrame, str, dict, or None
        Expanded HEALPix cells in the requested format.
    """
    if healpix_id is None:
        healpix_id = "healpix"
    resolution = validate_healpix_resolution(resolution)
    gdf = process_input_data_compact(input_data, healpix_id)
    healpix_ids = gdf[healpix_id].drop_duplicates().tolist()
    if not healpix_ids:
        print(f"No HEALPix IDs found in <{healpix_id}> field.")
        return
    try:
        max_order = max(uniq2orderpix(int(uid)).order for uid in healpix_ids)
        if resolution < max_order:
            print(f"Target expand resolution ({resolution}) must >= {max_order}.")
            return None
        healpix_ids_expand = healpix_expand(healpix_ids, resolution)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your HEALPix ID field and resolution."
        )
    if not healpix_ids_expand:
        return None

    rows = []
    for uniq_id in healpix_ids_expand:
        try:
            cell_polygon = healpix2geo(uniq_id, fix_antimeridian=fix_antimeridian)
            order, _ipix = uniq2orderpix(int(uniq_id))
            row = geodesic_dggs_to_geoseries(
                "healpix", uniq_id, order, cell_polygon, num_edges=4
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_healpix_expanded"
        else:
            output_name = "healpix_expanded"
    return convert_to_output_format(out_gdf, output_format, output_name)


def healpixexpand_cli():
    """Command-line interface for healpixexpand."""
    parser = argparse.ArgumentParser(description="HEALPix Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input HEALPix (GeoJSON, Shapefile, CSV, Parquet, or GeoDataFrame)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        choices=range(min_res, max_res + 1),
        help=f"Target HEALPix resolution [{min_res}..{max_res}]",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="HEALPix ID field")
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
        help="Antimeridian fixing method",
    )
    args = parser.parse_args()
    result = healpixexpand(
        args.input,
        args.resolution,
        healpix_id=args.cellid,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)


if __name__ == "__main__":
    healpixcompact_cli()
