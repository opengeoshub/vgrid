"""
RHEALPix Compact Module

This module provides functionality to compact and expand RHEALPix cells with flexible input and output formats.

Key Functions:
    rhealpix_compact: Compact a list of RHEALPix IDs with an optional parent depth
    rhealpixcompact: Compact a set of RHEALPix cells to their covering set
    rhealpixexpand: Expand (uncompact) a set of RHEALPix cells to a target resolution
    rhealpixcompact_cli: Command-line interface for compaction
    rhealpixexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from tqdm import tqdm
from vgrid.dggs.rhealpixdggs.dggs import WGS84_003 as rhealpix_dggs
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    aggregate_values,
    compact_cells,
    convert_to_output_format,
    prepare_compact_bags,
    process_input_data_compact,
    validate_rhealpix_resolution,
)
from vgrid.utils.constants import AGG_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.conversion.dggs2geo.rhealpix2geo import rhealpix2geo


def _rhealpix_parent(rid):
    if len(rid) <= 1:
        return None
    return rid[:-1]


def _rhealpix_children(parent):
    parent_uids = (parent[0],) + tuple(map(int, parent[1:]))
    parent_cell = rhealpix_dggs.cell(parent_uids)
    return {str(subcell) for subcell in parent_cell.subcells()}


def rhealpix_compact(rhealpix_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of RHEALPix cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    rhealpix_ids : list of str
        List of RHEALPix cell IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted RHEALPix cell IDs.

    Examples
    --------
    >>> rhealpix_ids = ["A0", "A1", "A2", "A3"]
    >>> compacted = rhealpix_compact(rhealpix_ids)
    >>> print(f"Compacted {len(rhealpix_ids)} cells to {len(compacted)} cells")
    """
    return compact_cells(
        rhealpix_ids,
        _rhealpix_parent,
        _rhealpix_children,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting rHEALPix",
    )


def rhealpix_expand(rhealpix_ids, resolution):
    """
    Expand a list of RHEALPix cells to the target resolution.

    Takes RHEALPix cells and expands them to their children at the specified resolution.

    Parameters
    ----------
    rhealpix_ids : list of str
        List of RHEALPix cell IDs to expand.
    resolution : int
        Target resolution to expand the cells to.

    Returns
    -------
    list of str
        List of expanded RHEALPix cell IDs at the target resolution.

    Examples
    --------
    >>> rhealpix_ids = ["A0"]
    >>> expanded = rhealpix_expand(rhealpix_ids, 3)
    >>> print(f"Expanded to {len(expanded)} cells at resolution 3")
    """
    expand_cells = []
    for rhealpix_id in rhealpix_ids:
        rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        cell_resolution = rhealpix_cell.resolution

        if cell_resolution >= resolution:
            expand_cells.append(rhealpix_cell)
        else:
            expand_cells.extend(
                rhealpix_cell.subcells(resolution)
            )  # Expand to the target level
    return expand_cells


def get_rhealpix_resolution(rhealpix_id):
    try:
        rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
        rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
        return rhealpix_cell.resolution
    except Exception as e:
        raise ValueError(f"Invalid cell ID <{rhealpix_id}>: {e}")


def rhealpixcompact(
    input_data,
    rhealpix_id="rhealpix",
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    fix_antimeridian=None,
    verbose=True,
):
    """
    Compact RHEALPix cells to their covering set at a given parent depth.

    Compacts a set of RHEALPix cells by replacing complete sets of children with
    their parent cells. Mixed input resolutions are allowed and ``depth`` limits
    how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg``. If ``agg`` is ``"count"``, ``numeric_col`` is
    ignored and the output ``count`` is the number of original input cells in
    each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing RHEALPix cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of RHEALPix cell IDs
    rhealpix_id : str, default "rhealpix"
        Name of the column containing RHEALPix cell IDs.
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
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The compacted RHEALPix cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = rhealpixcompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = rhealpixcompact(["A0", "A1", "A2", "A3"])

    >>> # Compact only one parent level
    >>> result = rhealpixcompact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = rhealpixcompact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = rhealpixcompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not rhealpix_id:
        rhealpix_id = "rhealpix"
    bags, agg_col = prepare_compact_bags(
        input_data,
        rhealpix_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="rHEALPix cells",
    )
    if bags is None:
        print(f"No rHEALPix tokens found in <{rhealpix_id}> field.")
        return

    rhealpix_tokens_compact = rhealpix_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not rhealpix_tokens_compact:
        return None
    rows = []
    for rhealpix_token_compact in tqdm(
        rhealpix_tokens_compact,
        desc="Building rHEALPix compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = rhealpix2geo(
                rhealpix_token_compact, fix_antimeridian=fix_antimeridian
            )
            rhealpix_uids = (rhealpix_token_compact[0],) + tuple(
                map(int, rhealpix_token_compact[1:])
            )
            rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
            cell_resolution = rhealpix_cell.resolution
            num_edges = 4
            if rhealpix_cell.ellipsoidal_shape() == "dart":
                num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "rhealpix",
                rhealpix_token_compact,
                cell_resolution,
                cell_polygon,
                num_edges,
            )
            row[agg_col] = aggregate_values(bags.get(rhealpix_token_compact, []), agg)
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_rhealpix_compacted"
        else:
            output_name = "rhealpix_compacted"
    return convert_to_output_format(out_gdf, output_format, output_name)


def rhealpixcompact_cli():
    parser = argparse.ArgumentParser(description="rHEALPix Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input rHEALPix (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="rHEALPix ID field")
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
    result = rhealpixcompact(
        input_data,
        rhealpix_id=cellid,
        output_format=output_format,
        fix_antimeridian=fix_antimeridian,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


def rhealpixexpand(
    input_data,
    resolution,
    rhealpix_id="rhealpix",
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Expand (uncompact) RHEALPix cells to a target resolution.

    Expands RHEALPix cells to their children at the specified resolution. The target resolution
    must be greater than or equal to the maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing RHEALPix cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of RHEALPix cell IDs
    resolution : int
        Target RHEALPix resolution to expand the cells to. Must be >= maximum input resolution.
    rhealpix_id : str, default "rhealpix"
        Name of the column containing RHEALPix cell IDs.
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path
    fix_antimeridian : Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The expanded RHEALPix cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = rhealpixexpand("cells.geojson", resolution=3)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = rhealpixexpand(["A0"], resolution=3)

    >>> # Expand to GeoJSON file
    >>> result = rhealpixexpand("cells.geojson", resolution=3, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    resolution = validate_rhealpix_resolution(resolution)
    gdf = process_input_data_compact(input_data, rhealpix_id)
    rhealpix_ids = sorted(gdf[rhealpix_id].drop_duplicates().tolist())
    if not rhealpix_ids:
        print(f"No rHEALPix tokens found in <{rhealpix_id}> field.")
        return
    try:
        # Get max resolution in input
        max_res = max(get_rhealpix_resolution(token) for token in rhealpix_ids)
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None
        expanded_cells = rhealpix_expand(rhealpix_ids, resolution)
        rhealpix_tokens_expand = [str(cell) for cell in expanded_cells]
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your rHEALPix ID field and resolution."
        )
    if not rhealpix_tokens_expand:
        return None
    rows = []
    for rhealpix_token_expand in rhealpix_tokens_expand:
        try:
            cell_polygon = rhealpix2geo(
                rhealpix_token_expand, fix_antimeridian=fix_antimeridian
            )
            rhealpix_uids = (rhealpix_token_expand[0],) + tuple(
                map(int, rhealpix_token_expand[1:])
            )
            rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
            cell_resolution = resolution
            num_edges = 4
            if rhealpix_cell.ellipsoidal_shape() == "dart":
                num_edges = 3
            row = geodesic_dggs_to_geoseries(
                "rhealpix",
                rhealpix_token_expand,
                cell_resolution,
                cell_polygon,
                num_edges,
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")
    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_rhealpix_expanded"
        else:
            output_name = "rhealpix_expanded"
    return convert_to_output_format(out_gdf, output_format, output_name)


def rhealpixexpand_cli():
    parser = argparse.ArgumentParser(description="rHEALPix Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input rHEALPix (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target rHEALPix resolution to expand to (must be greater than input cells)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="rHEALPix ID field")
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
    result = rhealpixexpand(
        input_data,
        resolution,
        rhealpix_id=cellid,
        output_format=output_format,
        fix_antimeridian=fix_antimeridian,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)
