"""
RHEALPix Compact Module

This module provides functionality to compact and expand RHEALPix cells with flexible input and output formats.

Key Functions:
    rhealpixcompact: Compact a set of RHEALPix cells to their minimal covering set
    rhealpixexpand: Expand (uncompact) a set of RHEALPix cells to a target resolution
    rhealpixcompact_cli: Command-line interface for compaction
    rhealpixexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from vgrid.dggs.rhealpixdggs.dggs import WGS84_003 as rhealpix_dggs
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    process_input_data_compact,
    convert_to_output_format,
    validate_rhealpix_resolution,
)
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.conversion.dggs2geo.rhealpix2geo import rhealpix2geo
from collections import defaultdict


def rhealpix_compact(rhealpix_ids):
    """
    Compact a list of RHEALPix cell IDs to their minimal covering set.

    Groups RHEALPix cells by their parents and replaces complete sets of children
    with their parent cells, repeating until no more compaction is possible.

    Parameters
    ----------
    rhealpix_ids : list of str
        List of RHEALPix cell IDs to compact.

    Returns
    -------
    list of str
        Sorted list of compacted RHEALPix cell IDs representing the minimal covering set.

    Examples
    --------
    >>> rhealpix_ids = ["A0", "A1", "A2", "A3"]
    >>> compacted = rhealpix_compact(rhealpix_ids)
    >>> print(f"Compacted {len(rhealpix_ids)} cells to {len(compacted)} cells")
    """
    rhealpix_ids = sorted(set(rhealpix_ids))  # dedupe, stable order

    # Main loop for compaction
    while True:
        grouped_rhealpix_ids = defaultdict(set)

        # Group cells by their parent
        for rhealpix_id in rhealpix_ids:
            if len(rhealpix_id) > 1:  # Ensure there's a valid parent
                parent = rhealpix_id[:-1]
                grouped_rhealpix_ids[parent].add(rhealpix_id)

        new_rhealpix_ids = set(rhealpix_ids)
        changed = False

        # Check if we can replace children with parent
        for parent, children in grouped_rhealpix_ids.items():
            parent_uids = (parent[0],) + tuple(
                map(int, parent[1:])
            )  # Assuming parent is a string like 'A0'
            parent_cell = rhealpix_dggs.cell(
                parent_uids
            )  # Retrieve the parent cell object

            # Generate the subcells for the parent at the next resolution
            subcells_at_next_res = set(
                str(subcell) for subcell in parent_cell.subcells()
            )  # Collect subcells as strings

            # Check if the current children match the subcells at the next resolution
            if children == subcells_at_next_res:
                new_rhealpix_ids.difference_update(children)  # Remove children
                new_rhealpix_ids.add(parent)  # Add the parent
                changed = True  # A change occurred

        if not changed:
            break  # Stop if no more compaction is possible
        rhealpix_ids = new_rhealpix_ids  # Continue compacting

    return sorted(rhealpix_ids)  # Sorted for consistency


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
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Compact RHEALPix cells to their minimal covering set.

    Compacts a set of RHEALPix cells by replacing complete sets of children with their parent cells,
    repeating until no more compaction is possible. Supports flexible input and output formats.

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
        The compacted RHEALPix cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = rhealpixcompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = rhealpixcompact(["A0", "A1", "A2", "A3"])

    >>> # Compact to GeoJSON file
    >>> result = rhealpixcompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    gdf = process_input_data_compact(input_data, rhealpix_id)
    rhealpix_ids = sorted(gdf[rhealpix_id].drop_duplicates().tolist())
    if not rhealpix_ids:
        print(f"No rHEALPix tokens found in <{rhealpix_id}> field.")
        return
    try:
        rhealpix_tokens_compact = rhealpix_compact(rhealpix_ids)
    except Exception:
        raise Exception("Compact cells failed. Please check your rHEALPix ID field.")
    if not rhealpix_tokens_compact:
        return None
    rows = []
    for rhealpix_token_compact in rhealpix_tokens_compact:
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
