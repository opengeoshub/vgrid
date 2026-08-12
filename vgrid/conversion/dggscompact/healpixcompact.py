"""
HEALPix Compact Module

This module provides functionality to compact and expand HEALPix UNIQ cells with
flexible input and output formats.

Key Functions:
    healpixcompact: Compact a set of HEALPix cells to their minimal covering set
    healpixexpand: Expand (uncompact) HEALPix cells to a target resolution
    healpixcompact_cli: Command-line interface for compaction
    healpixexpand_cli: Command-line interface for expansion
"""

import os
import argparse
from collections import defaultdict

import geopandas as gpd

from vgrid.conversion.dggs2geo.healpix2geo import healpix2geo
from vgrid.dggs.healpix import (
    nestDescendants,
    orderpix2uniq,
    uniq2orderpix,
    uniqChildren,
)
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    convert_to_output_format,
    process_input_data_compact,
    validate_healpix_resolution,
)
from vgrid.utils.constants import DGGS_TYPES, OUTPUT_FORMATS, STRUCTURED_FORMATS

min_res = DGGS_TYPES["healpix"]["min_res"]
max_res = DGGS_TYPES["healpix"]["max_res"]


def healpix_compact(healpix_ids):
    """
    Compact HEALPix UNIQ IDs by replacing complete sibling sets with their parent.

    Parameters
    ----------
    healpix_ids : list of int or str
        HEALPix UNIQ cell IDs.

    Returns
    -------
    list of int
        Compacted HEALPix UNIQ cell IDs.
    """
    ids = {int(uid) for uid in healpix_ids}
    changed = True
    while changed:
        changed = False
        children_of = defaultdict(set)
        for uid in ids:
            order, ipix = uniq2orderpix(uid)
            if order == 0:
                continue
            parent = orderpix2uniq(order - 1, ipix >> 2)
            children_of[parent].add(uid)

        for parent, kids in children_of.items():
            expected = set(uniqChildren(parent))
            if expected <= ids:
                ids -= expected
                ids.add(parent)
                changed = True
    return list(ids)


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
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Compact HEALPix cells to their minimal covering set.

    Parameters
    ----------
    input_data : str, dict, GeoDataFrame, or list
        Input containing HEALPix UNIQ cell IDs.
    healpix_id : str, optional
        Column name for HEALPix IDs. Defaults to ``"healpix"``.
    output_format : str, default "gpd"
        Output format.
    fix_antimeridian : str, optional
        Antimeridian fixing method.

    Returns
    -------
    GeoDataFrame, str, dict, or None
        Compacted HEALPix cells in the requested format.
    """
    if healpix_id is None:
        healpix_id = "healpix"
    gdf = process_input_data_compact(input_data, healpix_id)
    healpix_ids = gdf[healpix_id].drop_duplicates().tolist()
    if not healpix_ids:
        print(f"No HEALPix IDs found in <{healpix_id}> field.")
        return
    try:
        healpix_ids_compact = healpix_compact(healpix_ids)
    except Exception:
        raise Exception("Compact cells failed. Please check your HEALPix ID field.")
    if not healpix_ids_compact:
        return None

    rows = []
    for uniq_id in healpix_ids_compact:
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
    args = parser.parse_args()
    result = healpixcompact(
        args.input,
        healpix_id=args.cellid,
        output_format=args.output_format,
        fix_antimeridian=args.fix_antimeridian,
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
