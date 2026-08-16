"""
ISEA3H Compact Module

This module provides functionality to compact and expand ISEA3H cells with flexible input and output formats.

Key Functions:
    isea3hcompact: Compact a set of ISEA3H cells to their minimal covering set
    isea3hexpand: Expand (uncompact) a set of ISEA3H cells to a target resolution
    isea3hcompact_cli: Command-line interface for compaction
    isea3hexpand_cli: Command-line interface for expansion

Note: This module is only supported on Windows systems due to OpenEaggr dependency.
"""

import os
import argparse
import geopandas as gpd
from collections import defaultdict
import platform

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.dggs.eaggr.enums.model import Model
    from vgrid.utils.constants import ISEA3H_ACCURACY_RES_DICT

    isea3h_dggs = Eaggr(Model.ISEA3H)

from vgrid.conversion.dggs2geo.isea3h2geo import isea3h2geo
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import process_input_data_compact, convert_to_output_format
from pyproj import Geod

geod = Geod(ellps="WGS84")


# --- ISEA3H Compaction/Expansion Logic ---
def get_isea3h_cell_children(isea3h_cell, resolution):
    """Recursively expands a DGGS cell until all children reach the desired resolution."""
    isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)
    cell_accuracy = isea3h2point._accuracy
    cell_resolution = ISEA3H_ACCURACY_RES_DICT.get(cell_accuracy)

    if cell_resolution >= resolution:
        return [
            isea3h_cell
        ]  # Base case: return the cell if it meets/exceeds resolution

    expanded_cells = []
    children = isea3h_dggs.get_dggs_cell_children(isea3h_cell)

    for child in children:
        expanded_cells.extend(get_isea3h_cell_children(child, resolution))

    return expanded_cells


def get_isea3h_resolution(isea3h_id):
    """Get the resolution of an ISEA3H cell ID."""
    try:
        isea3h_cell = DggsCell(isea3h_id)
        cell_polygon = isea3h2geo(isea3h_id)
        cell_area_perimeter = geod.geometry_area_perimeter(cell_polygon)
        cell_perimeter = abs(cell_area_perimeter[1])

        isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)
        cell_accuracy = isea3h2point._accuracy

        avg_edge_len = cell_perimeter / 6
        cell_resolution = ISEA3H_ACCURACY_RES_DICT.get(cell_accuracy)

        if cell_resolution == 0:  # icosahedron faces at resolution = 0
            avg_edge_len = cell_perimeter / 3

        if cell_accuracy == 0.0:
            if round(avg_edge_len, 2) == 0.06:
                cell_resolution = 33
            elif round(avg_edge_len, 2) == 0.03:
                cell_resolution = 34
            elif round(avg_edge_len, 2) == 0.02:
                cell_resolution = 35
            elif round(avg_edge_len, 2) == 0.01:
                cell_resolution = 36
            elif round(avg_edge_len, 3) == 0.007:
                cell_resolution = 37
            elif round(avg_edge_len, 3) == 0.004:
                cell_resolution = 38
            elif round(avg_edge_len, 3) == 0.002:
                cell_resolution = 39
            elif round(avg_edge_len, 3) <= 0.001:
                cell_resolution = 40

        return cell_resolution
    except Exception as e:
        raise ValueError(f"Invalid cell ID <{isea3h_id}> : {e}")


def isea3h_compact(isea3h_ids):
    """
    Compact a list of ISEA3H cell IDs to their minimal covering set.

    Groups ISEA3H cells by their parents and replaces complete sets of children
    with their parent cells, repeating until no more compaction is possible.

    Parameters
    ----------
    isea3h_ids : list of str
        List of ISEA3H cell IDs to compact.

    Returns
    -------
    list of str
        Sorted list of compacted ISEA3H cell IDs representing the minimal covering set.

    Examples
    --------
    >>> isea3h_ids = ["A0", "A1", "A2", "A3", "A4", "A5"]
    >>> compacted = isea3h_compact(isea3h_ids)
    >>> print(f"Compacted {len(isea3h_ids)} cells to {len(compacted)} cells")
    """

    isea3h_ids = set(isea3h_ids)  # Remove duplicates
    cell_cache = {cell_id: DggsCell(cell_id) for cell_id in isea3h_ids}

    while True:
        grouped_by_parent = defaultdict(set)

        # Group cells by *all* their parents
        for cell_id in isea3h_ids:
            cell = cell_cache[cell_id]
            try:
                parents = isea3h_dggs.get_dggs_cell_parents(cell)
            except Exception as e:
                print(f"Error getting parents for {cell_id}: {e}")
                continue

            for parent in parents:
                parent_id = parent.get_cell_id()
                grouped_by_parent[parent_id].add(cell_id)

        new_isea3h_ids = set(isea3h_ids)
        changed = False

        for parent_id, children_ids in grouped_by_parent.items():
            parent_cell = DggsCell(parent_id)
            try:
                expected_children = set(
                    child.get_cell_id()
                    for child in isea3h_dggs.get_dggs_cell_children(parent_cell)
                )
            except Exception as e:
                print(f"Error getting children for parent {parent_id}: {e}")
                continue

            # Check for full match: only then compact
            if children_ids == expected_children:
                new_isea3h_ids.difference_update(children_ids)
                new_isea3h_ids.add(parent_id)
                cell_cache[parent_id] = parent_cell
                changed = True
            else:
                # Keep original children if they don't fully match expected subcells
                new_isea3h_ids.update(children_ids)

        if not changed:
            break  # Fully compacted

        isea3h_ids = new_isea3h_ids

    return sorted(isea3h_ids)


def isea3hcompact(
    input_data,
    isea3h_id=None,
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Compact ISEA3H cells to their minimal covering set.

    Compacts a set of ISEA3H cells by replacing complete sets of children with their parent cells,
    repeating until no more compaction is possible. Supports flexible input and output formats.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing ISEA3H cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of ISEA3H cell IDs
    isea3h_id : str, optional
        Name of the column containing ISEA3H cell IDs. Defaults to "isea3h".
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path
    fix_antimeridian : str, optional
        Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        Defaults to None when omitted.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The compacted ISEA3H cells in the specified format, or None if no valid cells found.

    Examples
    --------
    >>> # Compact from file
    >>> result = isea3hcompact("cells.geojson")
    >>> print(f"Compacted to {len(result)} cells")

    >>> # Compact from list
    >>> result = isea3hcompact(["A0", "A1", "A2", "A3", "A4", "A5"])

    >>> # Compact to GeoJSON file
    >>> result = isea3hcompact("cells.geojson", output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if not isea3h_id:
        isea3h_id = "isea3h"

    gdf = process_input_data_compact(input_data, isea3h_id)
    isea3h_ids = gdf[isea3h_id].drop_duplicates().tolist()

    if not isea3h_ids:
        print(f"No ISEA3H IDs found in <{isea3h_id}> field.")
        return

    try:
        isea3h_ids_compact = isea3h_compact(isea3h_ids)
    except Exception:
        raise Exception("Compact cells failed. Please check your ISEA3H ID field.")

    if not isea3h_ids_compact:
        return None

    rows = []
    for isea3h_id_compact in isea3h_ids_compact:
        try:
            cell_polygon = isea3h2geo(
                isea3h_id_compact, fix_antimeridian=fix_antimeridian
            )
            cell_resolution = get_isea3h_resolution(isea3h_id_compact)
            num_edges = 6  # ISEA3H cells are hexagonal
            row = geodesic_dggs_to_geoseries(
                "isea3h", isea3h_id_compact, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    file_formats = ["csv", "geojson", "shapefile", "gpkg", "parquet", "geoparquet"]
    output_name = None
    if output_format in file_formats:
        ext_map = {
            "csv": ".csv",
            "geojson": ".geojson",
            "shapefile": ".shp",
            "gpkg": ".gpkg",
            "parquet": ".parquet",
            "geoparquet": ".parquet",
        }
        ext = ext_map.get(output_format, "")
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_isea3h_compacted{ext}"
        else:
            output_name = f"isea3h_compacted{ext}"

    return convert_to_output_format(out_gdf, output_format, output_name)


def isea3hcompact_cli():
    """Command-line interface for ISEA3H compaction."""
    parser = argparse.ArgumentParser(description="ISEA3H Compact")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input ISEA3H (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="ISEA3H ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default=None,
        help="Output format (None, csv, geojson, shapefile, gpd, geojson_dict, gpkg, geoparquet)",
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

    result = isea3hcompact(
        input_data,
        isea3h_id=cellid,
        output_format=output_format,
        fix_antimeridian=args.fix_antimeridian,
    )

    if output_format is None:
        print(result)
    elif output_format in [
        "csv",
        "geojson",
        "geojson_dict",
        "shapefile",
        "gpkg",
        "geoparquet",
        "parquet",
    ]:
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            ext_map = {
                "csv": ".csv",
                "geojson": ".geojson",
                "geojson_dict": ".geojson",
                "shapefile": ".shp",
                "gpkg": ".gpkg",
                "parquet": ".parquet",
                "geoparquet": ".parquet",
            }
            ext = ext_map.get(output_format, "")
            output = f"{base}_isea3h_compacted{ext}"
        else:
            output = f"isea3h_compacted{ext_map.get(output_format, '')}"
        print(f"Output written to {output}")
    elif output_format in ["gpd", "geopandas"]:
        print(result)
    else:
        print("ISEA3H compact completed.")


def isea3h_expand(isea3h_ids, resolution):
    """
    Expand a list of ISEA3H cells to the target resolution.

    Takes ISEA3H cells and expands them to their children at the specified resolution.

    Parameters
    ----------
    isea3h_ids : list of str
        List of ISEA3H cell IDs to expand.
    resolution : int
        Target resolution to expand the cells to.

    Returns
    -------
    list of str
        List of expanded ISEA3H cell IDs at the target resolution.

    Examples
    --------
    >>> isea3h_ids = ["A0"]
    >>> expanded = isea3h_expand(isea3h_ids, 5)
    >>> print(f"Expanded to {len(expanded)} cells at resolution 5")
    """
    expand_cells = []
    for isea3h_id in isea3h_ids:
        isea3h_cell = DggsCell(isea3h_id)
        expand_cells.extend(get_isea3h_cell_children(isea3h_cell, resolution))
    return expand_cells


def isea3hexpand(
    input_data,
    resolution,
    isea3h_id=None,
    output_format="gpd",
    fix_antimeridian=None,
):
    """
    Expand (uncompact) ISEA3H cells to a target resolution.

    Expands ISEA3H cells to their children at the specified resolution. The target resolution
    must be greater than or equal to the maximum resolution of the input cells.

    Parameters
    ----------
    input_data : str, dict, geopandas.GeoDataFrame, or list
        Input data containing ISEA3H cell IDs. Can be:
        - File path (GeoJSON, Shapefile, CSV, Parquet)
        - URL to a file
        - GeoJSON dictionary
        - GeoDataFrame
        - List of ISEA3H cell IDs
    resolution : int
        Target ISEA3H resolution to expand the cells to. Must be >= maximum input resolution.
    isea3h_id : str, optional
        Name of the column containing ISEA3H cell IDs. Defaults to "isea3h".
    output_format : str, default "gpd"
        Output format. Options:
        - "gpd": Returns GeoPandas GeoDataFrame (default)
        - "csv": Returns CSV file path
        - "geojson": Returns GeoJSON file path
        - "geojson_dict": Returns GeoJSON FeatureCollection as Python dict
        - "parquet": Returns Parquet file path
        - "shapefile"/"shp": Returns Shapefile file path
        - "gpkg"/"geopackage": Returns GeoPackage file path
    fix_antimeridian : str, optional
        Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        Defaults to None when omitted.

    Returns
    -------
    geopandas.GeoDataFrame or str or dict or None
        The expanded ISEA3H cells in the specified format, or None if expansion fails.

    Examples
    --------
    >>> # Expand from file
    >>> result = isea3hexpand("cells.geojson", resolution=5)
    >>> print(f"Expanded to {len(result)} cells")

    >>> # Expand from list
    >>> result = isea3hexpand(["A0"], resolution=5)

    >>> # Expand to GeoJSON file
    >>> result = isea3hexpand("cells.geojson", resolution=5, output_format="geojson")
    >>> print(f"Saved to: {result}")
    """
    if isea3h_id is None:
        isea3h_id = "isea3h"

    gdf = process_input_data_compact(input_data, isea3h_id)
    isea3h_ids = gdf[isea3h_id].drop_duplicates().tolist()

    if not isea3h_ids:
        print(f"No ISEA3H IDs found in <{isea3h_id}> field.")
        return

    try:
        max_res = max(get_isea3h_resolution(isea3h_id) for isea3h_id in isea3h_ids)
        if resolution < max_res:
            print(f"Target expand resolution ({resolution}) must >= {max_res}.")
            return None

        isea3h_cells_expand = isea3h_expand(isea3h_ids, resolution)
        isea3h_ids_expand = [cell.get_cell_id() for cell in isea3h_cells_expand]
    except Exception:
        raise Exception(
            "Expand cells failed. Please check your ISEA3H ID field and resolution."
        )

    if not isea3h_ids_expand:
        return None

    rows = []
    for isea3h_id_expand in isea3h_ids_expand:
        try:
            cell_polygon = isea3h2geo(
                isea3h_id_expand, fix_antimeridian=fix_antimeridian
            )
            cell_resolution = resolution
            num_edges = 6  # ISEA3H cells are hexagonal
            row = geodesic_dggs_to_geoseries(
                "isea3h", isea3h_id_expand, cell_resolution, cell_polygon, num_edges
            )
            rows.append(row)
        except Exception:
            continue

    out_gdf = gpd.GeoDataFrame(rows, geometry="geometry", crs="EPSG:4326")

    file_formats = ["csv", "geojson", "shapefile", "gpkg", "parquet", "geoparquet"]
    output_name = None
    if output_format in file_formats:
        ext_map = {
            "csv": ".csv",
            "geojson": ".geojson",
            "shapefile": ".shp",
            "gpkg": ".gpkg",
            "parquet": ".parquet",
            "geoparquet": ".parquet",
        }
        ext = ext_map.get(output_format, "")
        if isinstance(input_data, str):
            base = os.path.splitext(os.path.basename(input_data))[0]
            output_name = f"{base}_isea3h_expanded{ext}"
        else:
            output_name = f"isea3h_expanded{ext}"

    return convert_to_output_format(out_gdf, output_format, output_name)


def isea3hexpand_cli():
    """Command-line interface for ISEA3H expansion."""
    parser = argparse.ArgumentParser(description="ISEA3H Expand (Uncompact)")
    parser.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Input ISEA3H (GeoJSON, Shapefile, CSV, Parquet, or pickled GeoDataFrame .gpd/.geopandas)",
    )
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=True,
        help="Target ISEA3H resolution to expand to (must be greater than input cells)",
    )
    parser.add_argument("-cellid", "--cellid", type=str, help="ISEA3H ID field")
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        default=None,
        help="Output format (None, csv, geojson, shapefile, gpd, geojson_dict, gpkg, geoparquet)",
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
    if platform.system() == "Windows":
        result = isea3hexpand(
            input_data,
            resolution,
            isea3h_id=cellid,
            output_format=output_format,
            fix_antimeridian=fix_antimeridian,
        )

        if output_format is None:
            print(result)
        elif output_format in [
            "csv",
            "geojson",
            "geojson_dict",
            "shapefile",
            "gpkg",
            "geoparquet",
            "parquet",
        ]:
            if isinstance(input_data, str):
                base = os.path.splitext(os.path.basename(input_data))[0]
                ext_map = {
                    "csv": ".csv",
                    "geojson": ".geojson",
                    "geojson_dict": ".geojson",
                    "shapefile": ".shp",
                    "gpkg": ".gpkg",
                    "parquet": ".parquet",
                    "geoparquet": ".parquet",
                }
                ext = ext_map.get(output_format, "")
                output = f"{base}_isea3h_expanded{ext}"
            else:
                output = f"isea3h_expanded{ext_map.get(output_format, '')}"
            print(f"Output written to {output}")
        elif output_format in ["gpd", "geopandas"]:
            print(result)
        else:
            print("ISEA3H expand completed.")
    else:
        print("ISEA3H is only supported on Windows systems")
