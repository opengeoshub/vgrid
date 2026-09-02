"""
Tilecode Compact Module

This module provides functionality to compact and expand Tilecode cells with flexible input and output formats.

Key Functions:
    tilecode_compact: Compact a list of Tilecode IDs with an optional parent depth
    tilecodecompact: Compact a set of Tilecode cells to their covering set
    tilecodeexpand: Expand (uncompact) a set of Tilecode cells to a target resolution
    tilecodecompact_cli: Command-line interface for compaction
    tilecodeexpand_cli: Command-line interface for expansion
"""

import os
import re
import argparse
import geopandas as gpd
from tqdm import tqdm
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import (
    aggregate_values,
    compact_cells,
    convert_to_output_format,
    prepare_compact_bags,
    process_input_data_compact,
)
from vgrid.utils.constants import AGG_OPTIONS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.dggs import tilecode
from vgrid.dggs.tilecode import tilecode_resolution
from vgrid.conversion.dggs2geo.tilecode2geo import tilecode2geo


def tilecode_compact(tilecode_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of Tilecode cell IDs by replacing complete child sets with parents.

    Groups cells by their immediate parent and replaces a parent when every child
    is present. Repeats until ``depth`` parent levels have been applied, or until
    no further compaction is possible.

    Parameters
    ----------
    tilecode_ids : list of str
        Tilecode cell IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted Tilecode cell IDs.
    """

    def parent_fn(tilecode_id):
        if not re.match(r"z(\d+)x(\d+)y(\d+)", tilecode_id):
            return None
        return tilecode.tilecode_parent(tilecode_id)

    def children_fn(parent):
        match = re.match(r"z(\d+)x(\d+)y(\d+)", parent)
        parent_res = int(match.group(1))
        return tilecode.tilecode_children(parent, parent_res + 1)

    return compact_cells(
        tilecode_ids,
        parent_fn,
        children_fn,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting Tilecode",
    )


def tilecodecompact(
    input_data,
    tilecode_id="tilecode",
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    verbose=True,
):
    """
    Compact Tilecode cells to their covering set at a given parent depth.

    Compacts a set of Tilecode cells by replacing complete sets of children with
    their parent cells. ``depth`` limits how far up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg`` (same options as ``h3bin``). If ``agg`` is
    ``"count"``, ``numeric_col`` is ignored and the output ``count`` is the
    number of original input cells in each compacted cell.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing Tilecode cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of Tilecode cell IDs
    tilecode_id : str, default "tilecode"
        Name of the column containing Tilecode cell IDs.
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
        The compacted Tilecode cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = tilecodecompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = tilecodecompact(["z3x1y1", "z3x1y2", "z3x2y1", "z3x2y2"])

    >>> # Compact only one parent level
    >>> result = tilecodecompact(cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = tilecodecompact(cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = tilecodecompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not tilecode_id:
        tilecode_id = "tilecode"

    bags, agg_col = prepare_compact_bags(
        input_data,
        tilecode_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="Tilecode cells",
    )
    if bags is None:
        print(f"No Tilecode IDs found in <{tilecode_id}> field.")
        return

    tilecode_ids_compact = tilecode_compact(
        list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )
    if not tilecode_ids_compact:
        return None

    rows = []
    for tilecode_id_compact in tqdm(
        tilecode_ids_compact,
        desc="Building Tilecode compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            cell_polygon = tilecode2geo(tilecode_id_compact)
            cell_resolution = tilecode_resolution(tilecode_id_compact)
            row = graticule_dggs_to_geoseries(
                "tilecode", tilecode_id_compact, cell_resolution, cell_polygon
            )
            row[agg_col] = aggregate_values(bags.get(tilecode_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_tilecode_compacted"
        else:
            output_name = "tilecode_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def tilecodecompact_cli():
    """Command-line interface for Tilecode compaction."""
    parser = argparse.ArgumentParser(description="Tilecode Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input Tilecode (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="Tilecode ID field")
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

    result = tilecodecompact(
        input_data,
        tilecode_id=cellid,
        output_format=output_format,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def tilecode_expand(tilecode_ids, resolution):
    """
    Expand a list of Tilecode cells to the target resolution.

    Takes Tilecode cells and expands them to their children at the specified resolution.

    Parameters
    ----------
    tilecode_ids : list of str
        List of Tilecode cell IDs to expand.
    resolution : int
        Target resolution to expand the cells to.

    Returns
    -------
    list of str
        List of expanded Tilecode cell IDs at the target resolution.

    Examples
    --------
    >>> tilecode_ids = ["z3x1y1"]
    >>> expanded = tilecode_expand(tilecode_ids, 5)
    >>> print(f"Expanded to {len(expanded)} cells at resolution 5")
    """
    expand_cells = []
    for tilecode_id in tilecode_ids:
        match = re.match(r"z(\d+)x(\d+)y(\d+)", tilecode_id)
        if not match:
            raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")
        cell_resolution = int(match.group(1))

        if cell_resolution >= resolution:
            expand_cells.append(tilecode_id)
        else:
            expand_cells.extend(
                tilecode.tilecode_children(tilecode_id, resolution)
            )  # Expand to the target level
    return expand_cells


def tilecodeexpand(
    input_data,
    resolution,
    tilecode_id="tilecode",
    output_format="gpd",
):
    """
    Expand (uncompact) Tilecode cells to a target resolution.

    Expands Tilecode cells to their children at the specified resolution. The target resolution
    must be greater than or equal to the maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing Tilecode cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of Tilecode cell IDs
    resolution : int
        Target Tilecode resolution to expand the cells to. Must be >= maximum input resolution.
    tilecode_id : str, default "tilecode"
        Name of the column containing Tilecode cell IDs.
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
        The expanded Tilecode cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = tilecodeexpand("cells.geojson", resolution=5)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = tilecodeexpand(["z3x1y1"], resolution=5)

    >>> # Expand to GeoJSON file
    >>> result = tilecodeexpand("cells.geojson", resolution=5, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """

    gdf = process_input_data_compact(input_data, tilecode_id)
    tilecode_ids = gdf[tilecode_id].drop_duplicates().tolist()

    if not tilecode_ids:
        print(f"No Tilecode IDs found in <{tilecode_id}> field.")
        return

    try:
        max_res = max(
            int(re.match(r"z(\d+)x(\d+)y(\d+)", tid).group(1)) for tid in tilecode_ids
        )
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None

        tilecode_ids_expand = tilecode_expand(tilecode_ids, resolution)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your Tilecode ID field and resolution."
        )

    if not tilecode_ids_expand:
        return None

    rows = []
    for tilecode_id_expand in tilecode_ids_expand:
        try:
            cell_polygon = tilecode2geo(tilecode_id_expand)
            cell_resolution = resolution
            row = graticule_dggs_to_geoseries(
                "tilecode", tilecode_id_expand, cell_resolution, cell_polygon
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_tilecode_expanded"
        else:
            output_name = "tilecode_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def tilecodeexpand_cli():
    """Command-line interface for Tilecode expansion."""
    parser = argparse.ArgumentParser(description="Tilecode Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input Tilecode (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target Tilecode resolution to expand to (must be greater than input cells)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="Tilecode ID field")
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

    result = tilecodeexpand(
        input_data,
        resolution,
        tilecode_id=cellid,
        output_format=output_format,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)
