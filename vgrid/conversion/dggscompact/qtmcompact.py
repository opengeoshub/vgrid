"""
QTM Compact Module

This module provides functionality to compact and expand QTM cells with flexible input and output formats.

Key Functions:
    qtmcompact: Compact a set of QTM cells to their minimal covering set
    qtmexpand: Expand (uncompact) a set of QTM cells to a target resolution
    qtmcompact_cli: Command-line interface for compaction
    qtmexpand_cli: Command-line interface for expansion
"""

import os
import argparse
import geopandas as gpd
from collections import defaultdict

from vgrid.conversion.dggs2geo.qtm2geo import qtm2geo
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import process_input_data_compact, convert_to_output_format
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.dggs import qtm


# --- QTM Compaction/Expansion Logic ---
def get_qtm_resolution(qtm_id):
    """Get the resolution of a QTM cell ID."""
    try:
        return len(qtm_id)
    except Exception as e:
        raise ValueError(f"Invalid QTM ID <{qtm_id}> : {e}")


def qtm_compact(qtm_ids):
    """
    Compact a list of QTM cell IDs to their minimal covering set.

    Groups QTM cells by their parents and replaces complete sets of children
    with their parent cells, repeating until no more compaction is possible.

    Parameters
    ----------
    qtm_ids : list of str
        List of QTM cell IDs to compact.

    Returns
    -------
    list of str
        Sorted list of compacted QTM cell IDs representing the minimal covering set.

    Examples
    --------
    >>> qtm_ids = ["A0", "A1", "A2", "A3"]
    >>> compacted = qtm_compact(qtm_ids)
    >>> print(f"Compacted {len(qtm_ids)} cells to {len(compacted)} cells")
    """
    qtm_ids = set(qtm_ids)  # Remove duplicates

    # Main loop for compaction
    while True:
        grouped_qtm_ids = defaultdict(set)

        # Group cells by their parent
        for qtm_id in qtm_ids:
            parent = qtm.qtm_parent(qtm_id)
            grouped_qtm_ids[parent].add(qtm_id)

        new_qtm_ids = set(qtm_ids)
        changed = False

        # Check if we can replace children with parent
        for parent, children in grouped_qtm_ids.items():
            next_resolution = len(parent) + 1
            # Generate the subcells for the parent at the next resolution
            childcells_at_next_res = set(
                childcell for childcell in qtm.qtm_children(parent, next_resolution)
            )

            # Check if the current children match the subcells at the next resolution
            if children == childcells_at_next_res:
                new_qtm_ids.difference_update(children)  # Remove children
                new_qtm_ids.add(parent)  # Add the parent
                changed = True  # A change occurred

        if not changed:
            break  # Stop if no more compaction is possible
        qtm_ids = new_qtm_ids  # Continue compacting

    return sorted(qtm_ids)  # Sorted for consistency


def qtmcompact(
    input_data,
    qtm_id="qtm",
    output_format="gpd",
):
    """
    Compact QTM cells to their minimal covering set.

    Compacts a set of QTM cells by replacing complete sets of children with their parent cells,
    repeating until no more compaction is possible. Supports flexible input and output formats.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing QTM cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of QTM cell IDs
    qtm_id : str, default "qtm"
        Name of the column containing QTM cell IDs.
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
        The compacted QTM cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = qtmcompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = qtmcompact(["A0", "A1", "A2", "A3"])

    >>> # Compact to GeoJSON file
    >>> result = qtmcompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """

    gdf = process_input_data_compact(input_data, qtm_id)
    qtm_ids = gdf[qtm_id].drop_duplicates().tolist()

    if not qtm_ids:
        print(f"No QTM IDs found in <{qtm_id}> field.")
        return

    try:
        qtm_ids_compact = qtm_compact(qtm_ids)
    except Exception:
        raise Exception("Compact cells failed. Please check your QTM ID field.")

    if not qtm_ids_compact:
        return None

    rows = []
    for qtm_id_compact in qtm_ids_compact:
        try:
            cell_polygon = qtm2geo(qtm_id_compact)
            cell_resolution = get_qtm_resolution(qtm_id_compact)
            num_edges = 3  # QTM cells are triangular
            row = geodesic_dggs_to_geoseries(
                "qtm", qtm_id_compact, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_qtm_compacted"
        else:
            output_name = "qtm_compacted"

    return convert_to_output_format(out_gdf, output_format, output_name)


def qtmcompact_cli():
    """Command-line interface for QTM compaction."""
    parser = argparse.ArgumentParser(description="QTM Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input QTM (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="QTM ID field")
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

    result = qtmcompact(
        input_data,
        qtm_id=cellid,
        output_format=output_format,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)


def qtm_expand(qtm_ids, resolution):
    """
    Expand a list of QTM cells to the target resolution.

    Takes QTM cells and expands them to their children at the specified resolution.

    Parameters
    ----------
    qtm_ids : list of str
        List of QTM cell IDs to expand.
    resolution : int
        Target resolution to expand the cells to.

    Returns
    -------
    list of str
        List of expanded QTM cell IDs at the target resolution.

    Examples
    --------
    >>> qtm_ids = ["A0"]
    >>> expanded = qtm_expand(qtm_ids, 3)
    >>> print(f"Expanded to {len(expanded)} cells at resolution 3")
    """
    expand_cells = []
    for qtm_id in qtm_ids:
        cell_resolution = len(qtm_id)
        if cell_resolution >= resolution:
            expand_cells.append(qtm_id)
        else:
            expand_cells.extend(
                qtm.qtm_children(qtm_id, resolution)
            )  # Expand to the target level
    return expand_cells


def qtmexpand(
    input_data,
    resolution,
    qtm_id="qtm",
    output_format="gpd",
):
    """
    Expand (uncompact) QTM cells to a target resolution.

    Expands QTM cells to their children at the specified resolution. The target resolution
    must be greater than or equal to the maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing QTM cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of QTM cell IDs
    resolution : int
        Target QTM resolution to expand the cells to. Must be >= maximum input resolution.
    qtm_id : str, default "qtm"
        Name of the column containing QTM cell IDs.
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
        The expanded QTM cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = qtmexpand("cells.geojson", resolution=3)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = qtmexpand(["A0"], resolution=3)

    >>> # Expand to GeoJSON file
    >>> result = qtmexpand("cells.geojson", resolution=3, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if qtm_id is None:
        qtm_id = "qtm"

    gdf = process_input_data_compact(input_data, qtm_id)
    qtm_ids = gdf[qtm_id].drop_duplicates().tolist()

    if not qtm_ids:
        print(f"No QTM IDs found in <{qtm_id}> field.")
        return

    try:
        max_res = max(len(qtm_id) for qtm_id in qtm_ids)
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None

        qtm_ids_expand = qtm_expand(qtm_ids, resolution)
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your QTM ID field and resolution."
        )

    if not qtm_ids_expand:
        return None

    rows = []
    for qtm_id_expand in qtm_ids_expand:
        try:
            cell_polygon = qtm2geo(qtm_id_expand)
            cell_resolution = resolution
            num_edges = 3  # QTM cells are triangular
            row = geodesic_dggs_to_geoseries(
                "qtm", qtm_id_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    output_name = None
    if output_format in OUTPUT_FORMATS:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_qtm_expanded"
        else:
            output_name = "qtm_expanded"

    return convert_to_output_format(out_gdf, output_format, output_name)


def qtmexpand_cli():
    """Command-line interface for QTM expansion."""
    parser = argparse.ArgumentParser(description="QTM Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input QTM (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target QTM resolution to expand to (must be greater than input cells)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="QTM ID field")
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

    result = qtmexpand(
        input_data,
        resolution,
        qtm_id=cellid,
        output_format=output_format,
    )

    if output_format in STRUCTURED_FORMATS:
        print(result)
