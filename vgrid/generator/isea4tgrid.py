"""
ISEA4T Grid Generator Module

Generates ISEA4T DGGS grids for specified resolutions with automatic cell generation and validation using hierarchical triangular grid system.

Key Functions:
- isea4t_grid(): Main grid generation function for whole world
- isea4t_grid_within_bbox(): Grid generation within bounding box
- isea4tgrid(): User-facing function with multiple output formats
- isea4tgrid_cli(): Command-line interface for grid generation
"""

import argparse
from tqdm import tqdm
from shapely.geometry import box
import geopandas as gpd
import platform

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.dggs.eaggr.enums.model import Model
    from vgrid.dggs.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.constants import ISEA4T_RES_ACCURACY_DICT

    isea4t_dggs = Eaggr(Model.ISEA4T)

from vgrid.utils.constants import (
    MAX_CELLS,
    ISEA4T_BASE_CELLS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
)
from vgrid.conversion.dggs2geo.isea4t2geo import isea4t2geo
from vgrid.utils.geometry import geodesic_dggs_to_geoseries

if platform.system() == "Windows":
    from vgrid.conversion.dggscompact.isea4tcompact import (
        isea4t_compact,
        get_isea4t_resolution,
    )
from vgrid.utils.io import (
    is_full_world_bbox,
    validate_bbox,
    validate_isea4t_resolution,
    convert_to_output_format,
    add_verbose_argument,
)


def get_isea4t_children_cells(base_cells, target_resolution):
    """
    Recursively generate DGGS cells for the desired resolution.
    """
    current_cells = base_cells
    for res in range(target_resolution):
        next_cells = []
        for cell in current_cells:
            children = isea4t_dggs.get_dggs_cell_children(DggsCell(cell))
            next_cells.extend([child._cell_id for child in children])
        current_cells = next_cells
    return current_cells


def get_isea4t_children_cells_within_bbox(bounding_cell, bbox, target_resolution):
    current_cells = [
        bounding_cell
    ]  # Start with a list containing the single bounding cell
    bounding_resolution = len(bounding_cell) - 2

    for res in range(bounding_resolution, target_resolution):
        next_cells = []
        for cell in current_cells:
            # Get the child cells for the current cell
            children = isea4t_dggs.get_dggs_cell_children(DggsCell(cell))
            for child in children:
                # Convert child cell to geometry
                child_shape = isea4t2geo(child._cell_id)
                if child_shape.intersects(bbox):
                    # Add the child cell ID to the next_cells list
                    next_cells.append(child._cell_id)
        if not next_cells:  # Break early if no cells remain
            break
        current_cells = (
            next_cells  # Update current_cells to process the next level of children
        )

    return current_cells


def _isea4t_row_from_id(isea4t_id, fix_antimeridian=None):
    cell_polygon = isea4t2geo(isea4t_id, fix_antimeridian=fix_antimeridian)
    cell_resolution = get_isea4t_resolution(isea4t_id)
    return geodesic_dggs_to_geoseries(
        "isea4t", isea4t_id, cell_resolution, cell_polygon, 3
    )


def isea4t_grid(resolution, fix_antimeridian=None, compact=False, verbose=True):
    resolution = validate_isea4t_resolution(resolution)
    cell_ids = get_isea4t_children_cells(ISEA4T_BASE_CELLS, resolution)
    if compact:
        cell_ids = isea4t_compact(cell_ids)
    isea4t_rows = []
    for cell_id in tqdm(cell_ids, desc="Generating ISEA4T DGGS", unit=" cells", disable=not verbose):
        isea4t_rows.append(_isea4t_row_from_id(cell_id, fix_antimeridian))
    return gpd.GeoDataFrame(isea4t_rows, geometry="geometry", crs="EPSG:4326")


def isea4t_grid_within_bbox(resolution, bbox, fix_antimeridian=None, compact=False, verbose=True):
    resolution = validate_isea4t_resolution(resolution)
    bbox = validate_bbox(bbox)
    if is_full_world_bbox(bbox):
        return isea4t_grid(
            resolution, fix_antimeridian=fix_antimeridian, compact=compact, verbose=verbose
        )

    accuracy = ISEA4T_RES_ACCURACY_DICT.get(resolution)
    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
    isea4t_shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(
        bounding_box_wkt, ShapeStringFormat.WKT, accuracy
    )
    isea4t_shape = isea4t_shapes[0]
    bbox_cells = isea4t_shape.get_shape().get_outer_ring().get_cells()
    bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
    bounding_children = get_isea4t_children_cells_within_bbox(
        bounding_cell.get_cell_id(), bounding_box, resolution
    )
    if compact:
        bounding_children = isea4t_compact(bounding_children)
    isea4t_rows = []
    for cell_id in tqdm(
        bounding_children, desc="Generating ISEA4T DGGS", unit=" cells", disable=not verbose
    ):
        isea4t_rows.append(_isea4t_row_from_id(cell_id, fix_antimeridian))
    return gpd.GeoDataFrame(isea4t_rows, geometry="geometry", crs="EPSG:4326")


def isea4t_grid_ids(resolution, compact=False):
    """
    Return a list of ISEA4T cell IDs for the whole world at a given resolution.
    """
    resolution = validate_isea4t_resolution(resolution)
    cell_ids = get_isea4t_children_cells(ISEA4T_BASE_CELLS, resolution)
    if compact:
        cell_ids = isea4t_compact(cell_ids)
    return [str(cid) for cid in cell_ids]


def isea4t_grid_within_bbox_ids(resolution, bbox, compact=False):
    """
    Return a list of ISEA4T cell IDs intersecting the given bounding box at a given resolution.
    """
    resolution = validate_isea4t_resolution(resolution)
    bbox = validate_bbox(bbox)
    if is_full_world_bbox(bbox):
        return isea4t_grid_ids(resolution, compact=compact)

    accuracy = ISEA4T_RES_ACCURACY_DICT.get(resolution)
    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt
    isea4t_shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(
        bounding_box_wkt, ShapeStringFormat.WKT, accuracy
    )
    isea4t_shape = isea4t_shapes[0]
    bbox_cells = isea4t_shape.get_shape().get_outer_ring().get_cells()
    bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
    bounding_children = get_isea4t_children_cells_within_bbox(
        bounding_cell.get_cell_id(), bounding_box, resolution
    )
    cell_ids = list(bounding_children or [])
    if compact:
        cell_ids = isea4t_compact(cell_ids)
    return cell_ids


def isea4tgrid(
    resolution, bbox=None, output_format="gpd", fix_antimeridian=None, compact=False, verbose=True
):
    """
    Generate ISEA4T DGGS grid for pure Python usage.
    Args:
        resolution (int): ISEA4T resolution [0..39]
        bbox (list[float]): [min_lon, min_lat, max_lon, max_lat]
        output_format (str): Output output_format ('geojson', 'csv', etc.)
        fix_antimeridian (str): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
        compact (bool, optional): Enable ISEA4T compact mode to reduce cell count.
    Returns:
        dict or list: GeoJSON FeatureCollection, list of ISEA4T cell IDs, or file path depending on output_format
    """
    if bbox is None:
        bbox = [-180, -90, 180, 90]
        total_cells = 20 * (4**resolution)
        if total_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {total_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
        gdf = isea4t_grid(
            resolution, fix_antimeridian=fix_antimeridian, compact=compact, verbose=verbose
        )
    else:
        gdf = isea4t_grid_within_bbox(
            resolution, bbox, fix_antimeridian=fix_antimeridian, compact=compact, verbose=verbose
        )
    output_name = f"isea4t_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def isea4tgrid_cli():
    parser = argparse.ArgumentParser(description="Generate OpenEaggr ISEA4T DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..39]"
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box in the output_format: min_lon min_lat max_lon max_lat (default is the whole world)",
    )
    parser.add_argument(
        "-f",
        "--output_format",
        type=str,
        choices=OUTPUT_FORMATS,
        default="gpd",
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
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable ISEA4T compact mode to reduce cell count",
    )
    add_verbose_argument(parser)
    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    fix_antimeridian = args.fix_antimeridian
    if platform.system() == "Windows":
        try:
            result = isea4tgrid(
                resolution,
                bbox,
                args.output_format,
                fix_antimeridian=fix_antimeridian,
                compact=args.compact,
                verbose=args.verbose,
            )
            if args.output_format in STRUCTURED_FORMATS:
                print(result)
        except ValueError as e:
            print(f"Error: {str(e)}")
            return
    else:
        print("ISEA4T is only supported on Windows systems")


if __name__ == "__main__":
    isea4tgrid_cli()
