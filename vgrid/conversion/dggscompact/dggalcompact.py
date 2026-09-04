"""
DGGAL Compact Module

This module provides functionality to compact and expand DGGAL cells with flexible input and output formats.

Key Functions:
    dggalcompact: Compact a set of DGGAL cells to their minimal covering set
    dggalexpand: Expand (uncompact) DGGAL cells to a target resolution or by depth
    dggalcompact_cli: Command-line interface for compaction
    dggalexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from tqdm import tqdm
from vgrid.conversion.dggs2geo.dggal2geo import dggal2geo
from vgrid.utils.io import (
    add_verbose_argument,
    aggregate_values,
    compact_cells,
    convert_to_output_format,
    prepare_compact_bags,
    process_input_data_compact,
    validate_dggal_type,
    validate_dggs_compact_depth,
    validate_dggs_expand_depth,
    validate_dggs_expand_resolution,
)
from vgrid.utils.constants import AGG_OPTIONS, DGGAL_TYPES, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries

from dggal import *

# Initialize dggal application
app = Application(appGlobals=globals())
pydggal_setup(app)
# --- DGGAL Compaction/Expansion Logic ---


def dggal_compact(dggs_type, zone_ids, depth=-1, bags=None, verbose=True):
    """
    Compact a list of DGGAL cell IDs by replacing complete child sets with parents.

    A zone may have multiple parents. Groups cells by every parent and replaces a
    parent when every child is present. Repeats until ``depth`` parent levels have
    been applied, or until no further compaction is possible.

    Parameters
    ----------
    dggs_type : str
        DGGAL DGGS type (e.g., "isea3h", "isea4t", "rhealpix").
    zone_ids : list of str
        DGGAL zone IDs to compact. Mixed resolutions are allowed.
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
        Sorted compacted DGGAL zone IDs.
    """
    dggs_type = validate_dggal_type(dggs_type)
    depth = validate_dggs_compact_depth(
        dggs_type, depth, max_res=int(DGGAL_TYPES[dggs_type]["max_res"])
    )
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = getattr(dggal, dggs_class_name)()

    def parent_fn(zone_id):
        zone = dggrs.getZoneFromTextID(zone_id)
        if dggrs.getZoneLevel(zone) <= 0:
            return None
        return [dggrs.getZoneTextID(p) for p in dggrs.getZoneParents(zone)]

    def children_fn(parent_zone_id):
        parent_zone = dggrs.getZoneFromTextID(parent_zone_id)
        return {dggrs.getZoneTextID(z) for z in dggrs.getZoneChildren(parent_zone)}

    return compact_cells(
        zone_ids,
        parent_fn,
        children_fn,
        depth=depth,
        bags=bags,
        verbose=verbose,
        desc="Compacting DGGAL",
    )


def dggalcompact(
    dggs_type,
    input_data,
    zone_id=None,
    depth=-1,
    agg="count",
    numeric_col=None,
    output_format="gpd",
    split_antimeridian=False,
    verbose=True,
):
    """
    Compact DGGAL cells to their covering set at a given parent depth.

    Compacts a set of DGGAL cells by replacing complete sets of children with their
    parent cells. Mixed input resolutions are allowed and ``depth`` limits how far
    up the hierarchy to merge.

    When a complete sibling set is replaced by its parent, original child values
    are combined with ``agg``. If ``agg`` is ``"count"``, ``numeric_col`` is
    ignored and the output ``count`` is the number of original input cells in
    each compacted cell.

    Parameters
    ----------
    dggs_type : str
        DGGAL DGGS type (e.g., "isea3h", "isea4t", "rhealpix").
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing DGGAL zone IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of DGGAL zone IDs
    zone_id : str, optional
        Name of the column containing DGGAL zone IDs. Defaults to "dggal_{dggs_type}".
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
    split_antimeridian : bool, optional
        When True, apply antimeridian fixing to the resulting polygons.
        Defaults to False when None or omitted.
    verbose : bool, default True
        Show tqdm progress bars. Use ``False`` to hide them.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The compacted DGGAL cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = dggalcompact("isea3h", "cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = dggalcompact("isea3h", ["A0", "A1", "A2", "A3"])

    >>> # Compact only one parent level
    >>> result = dggalcompact("isea3h", cells, depth=1)

    >>> # Mean of a numeric field on compacted parents
    >>> result = dggalcompact("isea3h", cells, agg="mean", numeric_col="value")

    >>> # Compact to GeoJSON file
    >>> result = dggalcompact("isea3h", "cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    dggs_type = validate_dggal_type(dggs_type)
    if not zone_id:
        zone_id = f"dggal_{dggs_type}"

    bags, agg_col = prepare_compact_bags(
        input_data,
        zone_id,
        agg=agg,
        numeric_col=numeric_col,
        verbose=verbose,
        label="DGGAL cells",
    )
    if bags is None:
        print(f"No DGGAL IDs found in <{zone_id}> field.")
        return

    # Create the appropriate DGGS instance
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = getattr(dggal, dggs_class_name)()

    dggal_ids_compact = dggal_compact(
        dggs_type, list(bags.keys()), depth=depth, bags=bags, verbose=verbose
    )

    if not dggal_ids_compact:
        print("Warning: Compaction returned no results, returning original data")
        gdf = process_input_data_compact(input_data, zone_id)
        return convert_to_output_format(gdf, output_format, f"{dggs_type}_original")

    rows = []
    for dggal_id_compact in tqdm(
        dggal_ids_compact,
        desc="Building DGGAL compact",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            # Get zone object to get resolution directly
            zone = dggrs.getZoneFromTextID(dggal_id_compact)
            cell_resolution = dggrs.getZoneLevel(zone)
            cell_polygon = dggal2geo(
                dggs_type, dggal_id_compact, split_antimeridian=split_antimeridian
            )
            num_edges = dggrs.countZoneEdges(zone)
            row = geodesic_dggs_to_geoseries(
                f"dggal_{dggs_type}",
                dggal_id_compact,
                cell_resolution,
                cell_polygon,
                num_edges,
            )
            row[agg_col] = aggregate_values(bags.get(dggal_id_compact, []), agg)
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_dggal_compacted"
        else:
            output_name = "dggal_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def dggalcompact_cli():
    parser = argparse.ArgumentParser(description="DGGAL Compact")
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGAL_TYPES.keys(),
        help="DGGAL DGGS type",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input DGGAL (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-zoneid", "--zoneid", type=str, help="DGGAL ID field")
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
    dggs_type = args.dggs_type
    zoneid = args.zoneid
    output_format = args.output_format

    result = dggalcompact(
        dggs_type,
        input_data,
        zone_id=zoneid,
        output_format=output_format,
        depth=args.depth,
        agg=args.agg,
        numeric_col=args.numeric_col,
        verbose=args.verbose,
    )
    if output_format in STRUCTURED_FORMATS:
        print(result)


def dggal_expand(dggs_type, zone_ids, resolution=None, depth=None, verbose=True):
    """
    Expand DGGAL zone IDs to a target resolution, or by a relative child depth.

    When ``resolution`` is set, ``depth`` is ignored and all cells are expanded
    to that absolute resolution. When only ``depth`` is set, ``resolution`` is
    ignored and each cell is expanded ``depth`` levels down (``1`` = direct
    children, ``2`` = grandchildren, and so on).
    """
    dggs_type = validate_dggal_type(dggs_type)
    max_res = int(DGGAL_TYPES[dggs_type]["max_res"])
    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = getattr(dggal, dggs_class_name)()

    if resolution is not None:
        resolution = validate_dggs_expand_resolution(
            dggs_type, resolution, max_res=max_res
        )
        expanded_cells = []
        for zid in tqdm(zone_ids, desc="Expanding DGGAL", unit=" cells", disable=not verbose):
            try:
                zone = dggrs.getZoneFromTextID(zid)
                current_res = dggrs.getZoneLevel(zone)

                if resolution < current_res:
                    print(
                        f"Warning: Target resolution {resolution} is lower than "
                        f"current resolution {current_res} for zone {zid}"
                    )
                    continue

                if resolution == current_res:
                    expanded_cells.append(zid)
                else:
                    sub_zones = dggrs.getSubZones(zone, resolution - current_res)
                    for sub_zone in sub_zones:
                        expanded_cells.append(dggrs.getZoneTextID(sub_zone))
            except Exception as e:
                print(f"Warning: Could not expand zone {zid}: {e}")
                continue
        return expanded_cells

    if depth is None:
        raise ValueError("Either resolution or depth must be specified.")
    depth = validate_dggs_expand_depth(dggs_type, depth, max_res=max_res)
    expanded_cells = []
    for zid in tqdm(zone_ids, desc="Expanding DGGAL", unit=" cells", disable=not verbose):
        try:
            zone = dggrs.getZoneFromTextID(zid)
            for sub_zone in dggrs.getSubZones(zone, depth):
                expanded_cells.append(dggrs.getZoneTextID(sub_zone))
        except Exception:
            continue
    return expanded_cells


def dggalexpand(
    dggs_type,
    input_data,
    resolution=None,
    zone_id=None,
    output_format="gpd",
    split_antimeridian=False,
    verbose=True,
    depth=None,
):
    """
    Expand (uncompact) DGGAL cells to a target resolution or by a relative depth.

    When ``resolution`` is set, ``depth`` is ignored and cells are expanded to
    that absolute resolution (must be >= the maximum input resolution). When
    only ``depth`` is set, ``resolution`` is ignored: mixed-resolution input is
    allowed and each cell is expanded to its descendants ``depth`` levels down.
    """
    dggs_type = validate_dggal_type(dggs_type)
    max_res = int(DGGAL_TYPES[dggs_type]["max_res"])
    if zone_id is None:
        zone_id = f"dggal_{dggs_type}"
    if resolution is not None:
        resolution = validate_dggs_expand_resolution(
            dggs_type, resolution, max_res=max_res
        )
    elif depth is not None:
        depth = validate_dggs_expand_depth(dggs_type, depth, max_res=max_res)
    else:
        raise ValueError("Either resolution or depth must be specified.")

    gdf = process_input_data_compact(input_data, zone_id)
    zone_ids = gdf[zone_id].drop_duplicates().tolist()

    if not zone_ids:
        print(f"No Zone IDs found in <{zone_id}> field.")
        return

    dggs_class_name = DGGAL_TYPES[dggs_type]["class_name"]
    dggrs = getattr(dggal, dggs_class_name)()

    try:
        if resolution is not None:
            max_input_res = 0
            for zid in zone_ids:
                try:
                    zone = dggrs.getZoneFromTextID(zid)
                    max_input_res = max(max_input_res, dggrs.getZoneLevel(zone))
                except Exception:
                    continue

            if resolution < max_input_res:
                print(f"Target expand resolution ({resolution}) must >= {max_input_res}.")
                return None
            zone_ids_expand = dggal_expand(
                dggs_type, zone_ids, resolution=resolution, verbose=verbose
            )
        else:
            zone_ids_expand = dggal_expand(dggs_type, zone_ids, depth=depth, verbose=verbose)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your Zone ID field, resolution, or depth."
        )
    if not zone_ids_expand:
        return None

    rows = []
    for zone_id_expand in tqdm(
        zone_ids_expand,
        desc="Building DGGAL expand",
        unit=" cells",
        disable=not verbose,
    ):
        try:
            zone = dggrs.getZoneFromTextID(zone_id_expand)
            cell_resolution = dggrs.getZoneLevel(zone)
            cell_polygon = dggal2geo(
                dggs_type, zone_id_expand, split_antimeridian=split_antimeridian
            )
            num_edges = dggrs.countZoneEdges(zone)
            row = geodesic_dggs_to_geoseries(
                f"dggal_{dggs_type}",
                zone_id_expand,
                cell_resolution,
                cell_polygon,
                num_edges,
            )
            rows.append(row)
        except Exception:
            continue
    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    ouput_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            ouput_name = f"{base}_dggal_expanded"
        else:
            ouput_name = "dggal_expanded"

    return convert_to_output_format(out_gdf, output_format, ouput_name)


def dggalexpand_cli():
    parser = argparse.ArgumentParser(description="DGGAL Expand (Uncompact)")
    parser.add_argument(
        "-dggs",
        "--dggs_type",
        type=str,
        required=True,
        choices=DGGAL_TYPES.keys(),
        help="DGGAL DGGS type",
    )
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input DGGAL (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument(
        "-r",
        "--resolution",
        type=int,
        help="Target DGGAL resolution to expand to (must be >= maximum input resolution). "
        "Ignores --depth.",
    )
    mode.add_argument(
        "-d",
        "--depth",
        type=int,
        help="Expand each cell by this many child levels (1 = direct children, "
        "2 = grandchildren, ...; 1 <= depth <= DGGAL max_res). "
        "Mixed input resolutions are allowed. Ignores --resolution.",
    )

    parser.add_argument("-zoneid", "--zoneid", type=str, help="DGGAL ID field")
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
    result = dggalexpand(
        args.dggs_type,
        args.input,
        resolution=args.resolution,
        zone_id=args.zoneid,
        output_format=args.output_format,
        depth=args.depth,
        verbose=args.verbose,
    )
    if args.output_format in STRUCTURED_FORMATS:
        print(result)
