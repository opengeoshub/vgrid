"""
OLC Compact Module

This module provides functionality to compact and expand OLC cells with flexible input and output formats.

Key Functions:
    olccompact: Compact a set of OLC cells to their minimal covering set
    olcexpand: Expand (uncompact) a set of OLC cells to a target resolution
    olccompact_cli: Command-line interface for compaction
    olcexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from collections import defaultdict

from vgrid.conversion.dggs2geo.olc2geo import olc2geo
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import process_input_data_compact, convert_to_output_format
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.dggs import olc


# --- OLC Compaction/Expansion Logic ---
def get_olc_resolution(olc_id):
    """Get the resolution of an OLC cell ID."""
    try:
        coord = olc.decode(olc_id)
        return coord.codeLength
    except Exception as e:
        raise ValueError(f"Invalid OLC ID <{olc_id}> : {e}")


def olc_compact(olc_ids):
    """
    Compact a list of OLC cell IDs to their minimal covering set.

    Groups OLC cells by their parents and replaces complete sets of children
    with their parent cells, repeating until no more compaction is possible.

    Parameters
    ----------
    olc_ids : list of str
        List of OLC cell IDs to compact.

    Returns
    -------
    list of str
        Sorted list of compacted OLC cell IDs representing the minimal covering set.

    Examples
    --------
    >>> olc_ids = ["7P28QPG4+4P7", "7P28QPG4+4P8", "7P28QPG4+4P9"]
    >>> compacted = olc_compact(olc_ids)
    >>> print(f"Compacted {len(olc_ids)} cells to {len(compacted)} cells")
    """
    olc_ids = set(olc_ids)  # Remove duplicates

    # Main loop for compaction
    while True:
        grouped_olc_ids = defaultdict(set)

        # Group cells by their parent
        for olc_id in olc_ids:
            parent = olc.olc_parent(olc_id)
            grouped_olc_ids[parent].add(olc_id)

        new_olc_ids = set(olc_ids)
        changed = False

        # Check if we can replace children with parent
        for parent, children in grouped_olc_ids.items():
            coord = olc.decode(parent)
            coord_len = coord.codeLength
            if coord_len <= 10:
                next_resolution = coord_len + 2
            else:
                next_resolution = coord_len + 1

            # Generate the subcells for the parent at the next resolution
            childcells_at_next_res = set(
                childcell for childcell in olc.olc_children(parent, next_resolution)
            )

            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_olc_ids.difference_update(children)  # Remove children
                new_olc_ids.add(parent)  # Add the parent
                changed = True  # A change occurred

        if not changed:
            break  # Stop if no more compaction is possible
        olc_ids = new_olc_ids  # Continue compacting

    return sorted(olc_ids)  # Sorted for consistency


def olccompact(
    input_data,
    olc_id=None,
    output_format="gpd",
):
    """
    Compact OLC cells to their minimal covering set.

    Compacts a set of OLC cells by replacing complete sets of children with their parent cells,
    repeating until no more compaction is possible. Supports flexible input and output formats.

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
        The compacted OLC cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = olccompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = olccompact(["7P28QPG4+4P7", "7P28QPG4+4P8"])

    >>> # Compact to GeoJSON file
    >>> result = olccompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not olc_id:
        olc_id = "olc"

    gdf = process_input_data_compact(input_data, olc_id)
    olc_ids = gdf[olc_id].drop_duplicates().tolist()

    if not olc_ids:
        print(f"No OLC IDs found in <{olc_id}> field.")
        return

    try:
        olc_ids_compact = olc_compact(olc_ids)
    except Exception:
        raise Exception("Compact cells failed. Please check your OLC ID field.")

    if not olc_ids_compact:
        return None

    rows = []
    for olc_id_compact in olc_ids_compact:
        try:
            cell_polygon = olc2geo(olc_id_compact)
            cell_resolution = get_olc_resolution(olc_id_compact)
            row = graticule_dggs_to_geoseries(
                "olc", olc_id_compact, cell_resolution, cell_polygon
            )
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

    args = parser.parse_args()
    input_data = args.input
    cellid = args.cellid
    output_format = args.output_format

    result = olccompact(
        input_data,
        olc_id=cellid,
        output_format=output_format,
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
