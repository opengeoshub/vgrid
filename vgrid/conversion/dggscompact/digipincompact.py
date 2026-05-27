"""
Digipin Compact Module

This module provides functionality to compact and expand DIGIPIN cells with flexible input and output formats.

Key Functions:
    digipincompact: Compact a set of DIGIPIN cells to their minimal covering set
    digipinexpand: Expand (uncompact) a set of DIGIPIN cells to a target resolution
    digipincompact_cli: Command-line interface for compaction
    digipinexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from collections import defaultdict
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import process_input_data_compact, convert_to_output_format
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.dggs.digipin import digipin_parent, digipin_children, digipin_resolution
from vgrid.conversion.dggs2geo.digipin2geo import digipin2geo


def digipin_compact(digipin_ids):
    """
    Compact a list of DIGIPIN cell IDs to their minimal covering set.

    Groups DIGIPIN cells by their parents and replaces complete sets of children
    with their parent cells, repeating until no more compaction is possible.

    Parameters
    ----------
    digipin_ids : list of str
        List of DIGIPIN cell IDs to compact.

    Returns
    -------
    list of str
        Sorted list of compacted DIGIPIN cell IDs representing the minimal covering set.

    Examples
    --------
    >>> digipin_ids = ["F3K-F", "F3K-C", "F3K-9", "F3K-8", "F3K-J", "F3K-3", "F3K-2", "F3K-7",
    ...                "F3K-K", "F3K-4", "F3K-5", "F3K-6", "F3K-L", "F3K-M", "F3K-P", "F3K-T"]
    >>> compacted = digipin_compact(digipin_ids)
    >>> print(f"Compacted {len(digipin_ids)} cells to {len(compacted)} cells")
    """
    digipin_ids = set(digipin_ids)  # Remove duplicates

    # Main loop for compaction
    while True:
        grouped_digipin_ids = defaultdict(set)

        # Group cells by their parent
        for digipin_id in digipin_ids:
            parent = digipin_parent(digipin_id)
            if parent != "Invalid DIGIPIN":  # Ensure there's a valid parent
                grouped_digipin_ids[parent].add(digipin_id)

        new_digipin_ids = set(digipin_ids)
        changed = False

        # Check if we can replace children with parent
        for parent, children in grouped_digipin_ids.items():
            # Generate the subcells for the parent at the next resolution
            parent_resolution = digipin_resolution(parent)
            if isinstance(parent_resolution, str):
                continue  # Skip invalid resolutions

            childcells_at_next_res = set(
                childcell
                for childcell in digipin_children(parent, parent_resolution + 1)
            )  # Collect subcells as strings

            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_digipin_ids.difference_update(children)  # Remove children
                new_digipin_ids.add(parent)  # Add the parent
                changed = True  # A change occurred

        if not changed:
            break  # Stop if no more compaction is possible
        digipin_ids = new_digipin_ids  # Continue compacting

    return sorted(digipin_ids)  # Sorted for consistency


def digipincompact(
    input_data,
    digipin_id="digipin",
    output_format="gpd",
):
    """
    Compact DIGIPIN cells to their minimal covering set.

    Compacts a set of DIGIPIN cells by replacing complete sets of children with their parent cells,
    repeating until no more compaction is possible. Supports flexible input and output formats.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing DIGIPIN cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of DIGIPIN cell IDs
    digipin_id : str, default "digipin"
        Name of the column containing DIGIPIN cell IDs.
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
        The compacted DIGIPIN cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = digipincompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = digipincompact(["F3K-F", "F3K-C", "F3K-9", "F3K-8"])

    >>> # Compact to GeoJSON file
    >>> result = digipincompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """

    gdf = process_input_data_compact(input_data, digipin_id)
    digipin_ids = gdf[digipin_id].drop_duplicates().tolist()

    if not digipin_ids:
        print(f"No DIGIPIN IDs found in <{digipin_id}> field.")
        return

    try:
        digipin_ids_compact = digipin_compact(digipin_ids)
    except Exception:
        raise Exception("Compact cells failed. Please check your DIGIPIN ID field.")

    if not digipin_ids_compact:
        return None

    rows = []
    for digipin_id_compact in digipin_ids_compact:
        try:
            cell_polygon = digipin2geo(digipin_id_compact)
            cell_resolution = digipin_resolution(digipin_id_compact)
            if isinstance(cell_resolution, str):
                continue  # Skip invalid resolutions
            row = graticule_dggs_to_geoseries(
                "digipin", digipin_id_compact, cell_resolution, cell_polygon
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_digipin_compacted"
        else:
            output_name = "digipin_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def digipincompact_cli():
    """Command-line interface for DIGIPIN compaction."""
    parser = argparse.ArgumentParser(description="DIGIPIN Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input DIGIPIN (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="DIGIPIN ID field")
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
    cellid = args.cellid
    output_format = args.output_format

    result = digipincompact(
        input_data,
        digipin_id=cellid,
        output_format=output_format,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def digipin_expand(digipin_ids, resolution):
    """
    Expand a list of DIGIPIN cells to the target resolution.

    Takes DIGIPIN cells and expands them to their children at the specified resolution.

    Parameters
    ----------
    digipin_ids : list of str
        List of DIGIPIN cell IDs to expand.
    resolution : int
        Target resolution to expand the cells to.

    Returns
    -------
    list of str
        List of expanded DIGIPIN cell IDs at the target resolution.

    Examples
    --------
    >>> digipin_ids = ["F3K"]
    >>> expanded = digipin_expand(digipin_ids, 5)
    >>> print(f"Expanded to {len(expanded)} cells at resolution 5")
    """
    expand_cells = []
    for digipin_id in digipin_ids:
        current_resolution = digipin_resolution(digipin_id)
        if isinstance(current_resolution, str):
            raise ValueError("Invalid DIGIPIN format.")

        if current_resolution >= resolution:
            expand_cells.append(digipin_id)
        else:
            expand_cells.extend(
                digipin_children(digipin_id, resolution)
            )  # Expand to the target level
    return expand_cells


def digipinexpand(
    input_data,
    resolution,
    digipin_id="digipin",
    output_format="gpd",
):
    """
    Expand (uncompact) DIGIPIN cells to a target resolution.

    Expands DIGIPIN cells to their children at the specified resolution. The target resolution
    must be greater than or equal to the maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing DIGIPIN cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of DIGIPIN cell IDs
    resolution : int
        Target DIGIPIN resolution to expand the cells to. Must be >= maximum input resolution.
    digipin_id : str, default "digipin"
        Name of the column containing DIGIPIN cell IDs.
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
        The expanded DIGIPIN cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = digipinexpand("cells.geojson", resolution=5)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = digipinexpand(["F3K"], resolution=5)

    >>> # Expand to GeoJSON file
    >>> result = digipinexpand("cells.geojson", resolution=5, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """

    gdf = process_input_data_compact(input_data, digipin_id)
    digipin_ids = gdf[digipin_id].drop_duplicates().tolist()

    if not digipin_ids:
        print(f"No DIGIPIN IDs found in <{digipin_id}> field.")
        return

    try:
        max_res = max(digipin_resolution(tid) for tid in digipin_ids)
        if isinstance(max_res, str):
            raise ValueError("Invalid DIGIPIN format.")
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None

        digipin_ids_expand = digipin_expand(digipin_ids, resolution)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your DIGIPIN ID field and resolution."
        )

    if not digipin_ids_expand:
        return None

    rows = []
    for digipin_id_expand in digipin_ids_expand:
        try:
            cell_polygon = digipin2geo(digipin_id_expand)
            cell_resolution = resolution
            row = graticule_dggs_to_geoseries(
                "digipin", digipin_id_expand, cell_resolution, cell_polygon
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_digipin_expanded"
        else:
            output_name = "digipin_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def digipinexpand_cli():
    """Command-line interface for DIGIPIN expansion."""
    parser = argparse.ArgumentParser(description="DIGIPIN Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input DIGIPIN (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target DIGIPIN resolution to expand to (must be greater than input cells)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="DIGIPIN ID field")
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

    result = digipinexpand(
        input_data,
        resolution,
        digipin_id=cellid,
        output_format=output_format,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)
