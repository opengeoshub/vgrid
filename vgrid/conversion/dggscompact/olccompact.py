"""
OLC Compact Module

This module provides functionality to compact and expand OLC cells with flexible input and output formats.

Key Functions:
    olc_compact: Compact a list of OLC IDs with an optional parent depth
    olccompact: Compact a set of OLC cells to their covering set
    olcexpand: Expand (uncompact) a set of OLC cells to a target resolution
    olccompact_cli: Command-line interface for compaction
    olcexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from tqdm import tqdm

from vgrid.conversion.dggs2geo.olc2geo import olc2geo
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import (
    aggregate_values,
    compact_cells,
    convert_to_output_format,
    prepare_compact_bags,
    process_input_data_compact,
)
from vgrid.utils.constants import AGG_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.dggs import olc


# --- OLC Compaction/Expansion Logic ---
def get_olc_resolution(olc_id):
    """Get the resolution of an OLC cell ID."""
    try:
        coord = olc.decode(olc_id)
        return coord.codeLength
    except Exception as e:
        raise ValueError(f"Invalid OLC ID <{olc_id}> : {e}")


def olc_compact(olc_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of OLC cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    olc_ids : list of str
        OLC cell IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted OLC cell IDs.
    """

    def parent_fn(olc_id):
        return olc.olc_parent(olc_id)

    def children_fn(parent):
        coord = olc.decode(parent)
        if coord.codeLength <= 10:
            next_res = coord.codeLength + 2
        else:
            next_res = coord.codeLength + 1
        return olc.olc_children(parent, next_res)

    return compact_cells(
        olc_ids,
        parent_fn,
        children_fn,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting OLC",
    )


def olccompact(
    input_data,
    olc_id=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    verbose=True,
):
    """
    Compact OLC cells to their covering set at a given parent depth.

    Compacts a set of OLC cells by replacing complete sets of children with
    their parent cells. ``depth`` limits how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg`` (same options as ``h3bin``). If ``agg`` is
    ``"count"``, ``numeric_col`` is ignored and the output ``count`` is the
    number of original input cells in each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing OLC cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of OLC cell IDs
    olc_id : str, optional
        Name of the column containing OLC cell IDs. Defaults to "olc".
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
        The compacted OLC cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = olccompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = olccompact(["7P28QPG4+4P7", "7P28QPG4+4P8"])

    >>> # Compact only one parent level
    >>> result = olccompact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = olccompact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = olccompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not olc_id:
        olc_id = "olc"

    bags, agg_col = prepare_compact_bags(
        input_data,
        olc_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="OLC cells",
    )
    if bags is None:
        print(f"No OLC IDs found in <{olc_id}> field.")
        return

    olc_ids_compact = olc_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not olc_ids_compact:
        return None

    rows = []
    for olc_id_compact in tqdm(
        olc_ids_compact,
        desc="Building OLC compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = olc2geo(olc_id_compact)
            cell_resolution = get_olc_resolution(olc_id_compact)
            row = graticule_dggs_to_geoseries(
                "olc", olc_id_compact, cell_resolution, cell_polygon
            )
            row[agg_col] = aggregate_values(bags.get(olc_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_olc_compacted"
        else:
            output_name = "olc_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def olccompact_cli():
    """Command-line interface for OLC compaction."""
    parser = argparse.ArgumentParser(description="OLC Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input OLC (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="OLC ID field")
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

    result = olccompact(
        input_data,
        olc_id=cellid,
        output_format=output_format,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def olc_expand(olc_ids, resolution):
    """
    Expand a list of OLC cells to the target resolution.

    Takes OLC cells and expands them to their children at the specified resolution.

    Parameters
    ----------
    olc_ids : list of str
        List of OLC cell IDs to expand.
    resolution : int
        Target resolution to expand the cells to.

    Returns
    -------
    list of str
        List of expanded OLC cell IDs at the target resolution.

    Examples
    --------
    >>> olc_ids = ["7P28QPG4+4P7"]
    >>> expanded = olc_expand(olc_ids, 8)
    >>> print(f"Expanded to {len(expanded)} cells at resolution 8")
    """
    expand_cells = []
    for olc_id in olc_ids:
        coord = olc.decode(olc_id)
        cell_resolution = coord.codeLength
        if cell_resolution >= resolution:
            expand_cells.append(olc_id)
        else:
            expand_cells.extend(
                olc.olc_children(olc_id, resolution)
            )  # Expand to the target level
    return expand_cells


def olcexpand(
    input_data,
    resolution,
    olc_id=None,
    output_format="gpd",
):
    """
    Expand (uncompact) OLC cells to a target resolution.

    Expands OLC cells to their children at the specified resolution. The target resolution
    must be greater than or equal to the maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing OLC cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of OLC cell IDs
    resolution : int
        Target OLC resolution to expand the cells to. Must be >= maximum input resolution.
    olc_id : str, optional
        Name of the column containing OLC cell IDs. Defaults to "olc".
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The expanded OLC cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = olcexpand("cells.geojson", resolution=8)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = olcexpand(["7P28QPG4+4P7"], resolution=8)

    >>> # Expand to GeoJSON file
    >>> result = olcexpand("cells.geojson", resolution=8, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if olc_id is None:
        olc_id = "olc"

    gdf = process_input_data_compact(input_data, olc_id)
    olc_ids = gdf[olc_id].drop_duplicates().tolist()

    if not olc_ids:
        print(f"No OLC IDs found in <{olc_id}> field.")
        return

    try:
        max_res = max(olc.decode(olc_id).codeLength for olc_id in olc_ids)
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None

        olc_ids_expand = olc_expand(olc_ids, resolution)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your OLC ID field and resolution."
        )

    if not olc_ids_expand:
        return None

    rows = []
    for olc_id_expand in olc_ids_expand:
        try:
            cell_polygon = olc2geo(olc_id_expand)
            cell_resolution = resolution
            row = graticule_dggs_to_geoseries(
                "olc", olc_id_expand, cell_resolution, cell_polygon
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_olc_expanded"
        else:
            output_name = "olc_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def olcexpand_cli():
    """Command-line interface for OLC expansion."""
    parser = argparse.ArgumentParser(description="OLC Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input OLC (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target OLC resolution to expand to (must be greater than input cells)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="OLC ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default="gpd",
        choices=OUTPUT_FORMATS,
        help="Output format",
    )

    args = parser.parse_args()
    input_data = args.input
    resolution = args.resolution
    cellid = args.cellid
    output_format = args.output_format

    result = olcexpand(
        input_data,
        resolution,
        olc_id=cellid,
        output_format=output_format,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)
