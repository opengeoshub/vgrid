"""
Geohash Compact Module

This module provides functionality to compact and expand Geohash cells with flexible input and output formats.

Key Functions:
    geohashcompact: Compact a set of Geohash cells to their minimal covering set
    geohashexpand: Expand (uncompact) a set of Geohash cells to a target resolution
    geohashcompact_cli: Command-line interface for compaction
    geohashexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from collections import defaultdict

from vgrid.conversion.dggs2geo.geohash2geo import geohash2geo
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import process_input_data_compact, convert_to_output_format
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.dggs import geohash


# --- Geohash Compaction/Expansion Logic ---
def get_geohash_resolution(geohash_id):
    """Get the resolution of a Geohash cell ID."""
    return len(geohash_id)


def geohash_compact(geohash_ids):
    """
    Compact a list of Geohash cell IDs to their minimal covering set.

    Groups Geohash cells by their parents and replaces complete sets of children
    with their parent cells, repeating until no more compaction is possible.

    Parameters
    ----------
    geohash_ids : list of str
        List of Geohash cell IDs to compact.

    Returns
    -------
    list of str
        Sorted list of compacted Geohash cell IDs representing the minimal covering set.

    Examples
    --------
    >>> geohash_ids = ["w3gvk1td8", "w3gvk1td9", "w3gvk1tdb"]
    >>> compacted = geohash_compact(geohash_ids)
    >>> print(f"Compacted {len(geohash_ids)} cells to {len(compacted)} cells")
    """
    geohash_ids = set(geohash_ids)  # Remove duplicates

    # Main loop for compaction
    while True:
        grouped_geohash_ids = defaultdict(set)

        # Group cells by their parent
        for geohash_id in geohash_ids:
            parent = geohash.geohash_parent(geohash_id)
            grouped_geohash_ids[parent].add(geohash_id)

        new_geohash_ids = set(geohash_ids)
        changed = False

        # Check if we can replace children with parent
        for parent, children in grouped_geohash_ids.items():
            parent_resolution = len(parent)
            # Generate the subcells for the parent at the next resolution
            childcells_at_next_res = set(
                childcell
                for childcell in geohash.geohash_children(parent, parent_resolution + 1)
            )

            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_geohash_ids.difference_update(children)  # Remove children
                new_geohash_ids.add(parent)  # Add the parent
                changed = True  # A change occurred

        if not changed:
            break  # Stop if no more compaction is possible
        geohash_ids = new_geohash_ids  # Continue compacting

    return sorted(geohash_ids)  # Sorted for consistency


def geohashcompact(
    input_data,
    geohash_id=None,
    output_format="gpd",
):
    """
    Compact Geohash cells to their minimal covering set.

    Compacts a set of Geohash cells by replacing complete sets of children with their parent cells,
    repeating until no more compaction is possible. Supports flexible input and output formats.

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
        The compacted Geohash cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = geohashcompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = geohashcompact(["w3gvk1td8", "w3gvk1td9"])

    >>> # Compact to GeoJSON file
    >>> result = geohashcompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not geohash_id:
        geohash_id = "geohash"

    gdf = process_input_data_compact(input_data, geohash_id)
    geohash_ids = gdf[geohash_id].drop_duplicates().tolist()

    if not geohash_ids:
        print(f"No Geohash IDs found in <{geohash_id}> field.")
        return

    try:
        geohash_ids_compact = geohash_compact(geohash_ids)
    except Exception:
        raise Exception("Compact cells failed. Please check your Geohash ID field.")

    if not geohash_ids_compact:
        return None

    rows = []
    for geohash_id_compact in geohash_ids_compact:
        try:
            cell_polygon = geohash2geo(geohash_id_compact)
            cell_resolution = get_geohash_resolution(geohash_id_compact)
            row = graticule_dggs_to_geoseries(
                "geohash", geohash_id_compact, cell_resolution, cell_polygon
            )
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

    args = parser.parse_args()
    input_data = args.input
    cellid = args.cellid
    output_format = args.output_format

    result = geohashcompact(
        input_data,
        geohash_id=cellid,
        output_format=output_format,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def geohash_expand(geohash_ids, resolution):
    """
    Expand a list of Geohash cells to the target resolution.

    Takes Geohash cells and expands them to their children at the specified resolution.

    Parameters
    ----------
    geohash_ids : list of str
        List of Geohash cell IDs to expand.
    resolution : int
        Target resolution to expand the cells to.

    Returns
    -------
    list of str
        List of expanded Geohash cell IDs at the target resolution.

    Examples
    --------
    >>> geohash_ids = ["w3gvk1td8"]
    >>> expanded = geohash_expand(geohash_ids, 5)
    >>> print(f"Expanded to {len(expanded)} cells at resolution 5")
    """
    expand_cells = []
    for geohash_id in geohash_ids:
        cell_resolution = len(geohash_id)
        if cell_resolution >= resolution:
            expand_cells.append(geohash_id)
        else:
            expand_cells.extend(
                geohash.geohash_children(geohash_id, resolution)
            )  # Expand to the target level
    return expand_cells


def geohashexpand(
    input_data,
    resolution,
    geohash_id=None,
    output_format="gpd",
):
    """
    Expand (uncompact) Geohash cells to a target resolution.

    Expands Geohash cells to their children at the specified resolution. The target resolution
    must be greater than or equal to the maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing Geohash cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of Geohash cell IDs
    resolution : int
        Target Geohash resolution to expand the cells to. Must be >= maximum input resolution.
    geohash_id : str, optional
        Name of the column containing Geohash cell IDs. Defaults to "geohash".
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
        The expanded Geohash cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = geohashexpand("cells.geojson", resolution=5)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = geohashexpand(["w3gvk1td8"], resolution=5)

    >>> # Expand to GeoJSON file
    >>> result = geohashexpand("cells.geojson", resolution=5, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if geohash_id is None:
        geohash_id = "geohash"

    gdf = process_input_data_compact(input_data, geohash_id)
    geohash_ids = gdf[geohash_id].drop_duplicates().tolist()

    if not geohash_ids:
        print(f"No Geohash IDs found in <{geohash_id}> field.")
        return

    try:
        max_res = max(len(geohash_id) for geohash_id in geohash_ids)
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None

        geohash_ids_expand = geohash_expand(geohash_ids, resolution)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your Geohash ID field and resolution."
        )

    if not geohash_ids_expand:
        return None

    rows = []
    for geohash_id_expand in geohash_ids_expand:
        try:
            cell_polygon = geohash2geo(geohash_id_expand)
            cell_resolution = resolution
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
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target Geohash resolution to expand to (must be greater than input cells)",
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

    args = parser.parse_args()
    input_data = args.input
    resolution = args.resolution
    cellid = args.cellid
    output_format = args.output_format

    result = geohashexpand(
        input_data,
        resolution,
        geohash_id=cellid,
        output_format=output_format,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)
