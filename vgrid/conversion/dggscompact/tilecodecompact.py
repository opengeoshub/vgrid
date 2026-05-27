"""
Tilecode Compact Module

This module provides functionality to compact and expand Tilecode cells with flexible input and output formats.

Key Functions:
    tilecodecompact: Compact a set of Tilecode cells to their minimal covering set
    tilecodeexpand: Expand (uncompact) a set of Tilecode cells to a target resolution
    tilecodecompact_cli: Command-line interface for compaction
    tilecodeexpand_cli: Command-line interface for expansion
"""

import os
import re
import argparse
import geopandas as gpd
from collections import defaultdict
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import process_input_data_compact, convert_to_output_format
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.dggs import tilecode
from vgrid.dggs.tilecode import tilecode_resolution
from vgrid.conversion.dggs2geo.tilecode2geo import tilecode2geo


def tilecode_compact(tilecode_ids):
    """
    Compact a list of Tilecode cell IDs to their minimal covering set.

    Groups Tilecode cells by their parents and replaces complete sets of children
    with their parent cells, repeating until no more compaction is possible.

    Parameters
    ----------
    tilecode_ids : list of str
        List of Tilecode cell IDs to compact.

    Returns
    -------
    list of str
        Sorted list of compacted Tilecode cell IDs representing the minimal covering set.

    Examples
    --------
    >>> tilecode_ids = ["z3x1y1", "z3x1y2", "z3x2y1", "z3x2y2"]
    >>> compacted = tilecode_compact(tilecode_ids)
    >>> print(f"Compacted {len(tilecode_ids)} cells to {len(compacted)} cells")
    """
    tilecode_ids = set(tilecode_ids)  # Remove duplicates

    # Main loop for compaction
    while True:
        grouped_tilecode_ids = defaultdict(set)

        # Group cells by their parent
        for tilecode_id in tilecode_ids:
            match = re.match(r"z(\d+)x(\d+)y(\d+)", tilecode_id)
            if match:  # Ensure there's a valid parent
                parent = tilecode.tilecode_parent(tilecode_id)
                grouped_tilecode_ids[parent].add(tilecode_id)

        new_tilecode_ids = set(tilecode_ids)
        changed = False

        # Check if we can replace children with parent
        for parent, children in grouped_tilecode_ids.items():
            # Generate the subcells for the parent at the next resolution
            match = re.match(r"z(\d+)x(\d+)y(\d+)", parent)
            parent_resolution = int(match.group(1))

            childcells_at_next_res = set(
                childcell
                for childcell in tilecode.tilecode_children(
                    parent, parent_resolution + 1
                )
            )  # Collect subcells as strings

            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_tilecode_ids.difference_update(children)  # Remove children
                new_tilecode_ids.add(parent)  # Add the parent
                changed = True  # A change occurred

        if not changed:
            break  # Stop if no more compaction is possible
        tilecode_ids = new_tilecode_ids  # Continue compacting

    return sorted(tilecode_ids)  # Sorted for consistency


def tilecodecompact(
    input_data,
    tilecode_id="tilecode",
    output_format="gpd",
):
    """
    Compact Tilecode cells to their minimal covering set.

    Compacts a set of Tilecode cells by replacing complete sets of children with their parent cells,
    repeating until no more compaction is possible. Supports flexible input and output formats.

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
        The compacted Tilecode cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = tilecodecompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = tilecodecompact(["z3x1y1", "z3x1y2", "z3x2y1", "z3x2y2"])

    >>> # Compact to GeoJSON file
    >>> result = tilecodecompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """

    gdf = process_input_data_compact(input_data, tilecode_id)
    tilecode_ids = gdf[tilecode_id].drop_duplicates().tolist()

    if not tilecode_ids:
        print(f"No Tilecode IDs found in <{tilecode_id}> field.")
        return

    try:
        tilecode_ids_compact = tilecode_compact(tilecode_ids)
    except Exception:
        raise Exception("Compact cells failed. Please check your Tilecode ID field.")

    if not tilecode_ids_compact:
        return None

    rows = []
    for tilecode_id_compact in tilecode_ids_compact:
        try:
            cell_polygon = tilecode2geo(tilecode_id_compact)
            cell_resolution = tilecode_resolution(tilecode_id_compact)
            row = graticule_dggs_to_geoseries(
                "tilecode", tilecode_id_compact, cell_resolution, cell_polygon
            )
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

    args = parser.parse_args()
    input_data = args.input
    cellid = args.cellid
    output_format = args.output_format

    result = tilecodecompact(
        input_data,
        tilecode_id=cellid,
        output_format=output_format,
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
