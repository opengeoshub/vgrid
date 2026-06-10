"""
ISEA3H Grid Generator Module

Generates ISEA3H DGGS grids for specified resolutions with automatic cell generation and validation using hierarchical hexagonal grid system.

Key Functions:
- isea3h_grid(): Main grid generation function for whole world
- isea3h_grid_within_bbox(): Grid generation within bounding box
- isea3hgrid(): User-facing function with multiple output formats
- isea3hgrid_cli(): Command-line interface for grid generation
"""

import argparse
import geopandas as gpd
from shapely.geometry import box
from tqdm import tqdm
import platform

if platform.system() == "Windows":
    from vgrid.dggs.eaggr.eaggr import Eaggr
    from vgrid.dggs.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.dggs.eaggr.enums.model import Model
    from vgrid.dggs.eaggr.enums.shape_string_format import ShapeStringFormat

    isea3h_dggs = Eaggr(Model.ISEA3H)

from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import validate_bbox, validate_isea3h_resolution, convert_to_output_format
from vgrid.conversion.dggs2geo.isea3h2geo import isea3h2geo

if platform.system() == "Windows":
    from vgrid.conversion.dggscompact.isea3hcompact import (
        isea3h_compact,
        get_isea3h_resolution,
    )

from vgrid.utils.constants import (
    ISEA3H_ACCURACY_RES_DICT,
    ISEA3H_RES_ACCURACY_DICT,
    MAX_CELLS,
    ISEA3H_BASE_CELLS,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
)

from pyproj import Geod

geod = Geod(ellps="WGS84")


def get_isea3h_children_cells(base_cells, target_resolution):
    """
    Recursively generate DGGS cells for the desired resolution, returning only the cells at the target resolution.
    """
    current_cells = base_cells
    for res in range(target_resolution):
        next_cells = []
        seen_cells = set()
        for cell in current_cells:
            children = isea3h_dggs.get_dggs_cell_children(DggsCell(cell))
            for child in children:
                if child._cell_id not in seen_cells:
                    seen_cells.add(child._cell_id)
                    next_cells.append(child._cell_id)
        current_cells = next_cells
    return current_cells


def get_isea3h_children_cells_within_bbox(bounding_cell, bbox, target_resolution):
    """
    Recursively generate DGGS cells within a bounding box, returning only the cells at the target resolution.
    """
    current_cells = [
        bounding_cell
    ]  # Start with a list containing the single bounding cell
    bounding_cell2point = isea3h_dggs.convert_dggs_cell_to_point(
        DggsCell(bounding_cell)
    )
    accuracy = bounding_cell2point._accuracy
    bounding_resolution = ISEA3H_ACCURACY_RES_DICT.get(accuracy)

    if bounding_resolution <= target_resolution:
        for res in range(bounding_resolution, target_resolution):
            next_cells = []
            seen_cells = set()
            for cell in current_cells:
                # Get the child cells for the current cell
                children = isea3h_dggs.get_dggs_cell_children(DggsCell(cell))
                for child in children:
                    if child._cell_id not in seen_cells:
                        child_shape = isea3h2geo(child._cell_id)
                        if child_shape.intersects(bbox):
                            seen_cells.add(child._cell_id)
                            next_cells.append(child._cell_id)
            if not next_cells:  # Break early if no cells remain
                break
            current_cells = (
                next_cells  # Update current_cells to process the next level of children
            )

        return current_cells
    else:
        # print('Bounding box area is < 0.028 square meters. Please select a bigger bounding box')
        return None


def _isea3h_row_from_id(isea3h_id, fix_antimeridian=None):
    cell_polygon = isea3h2geo(isea3h_id, fix_antimeridian=fix_antimeridian)
    cell_resolution = get_isea3h_resolution(isea3h_id)
    num_edges = 6 if cell_resolution > 0 else 3
    return geodesic_dggs_to_geoseries(
        "isea3h", isea3h_id, cell_resolution, cell_polygon, num_edges
    )


def isea3h_grid(resolution, fix_antimeridian=None, compact=False):
    """
    Generate DGGS cells and convert them to GeoJSON features.
    """
    resolution = validate_isea3h_resolution(resolution)
    cell_ids = get_isea3h_children_cells(ISEA3H_BASE_CELLS, resolution)
    if compact:
        cell_ids = isea3h_compact(cell_ids)
    records = []
    for cell_id in tqdm(cell_ids, desc="Generating ISEA3H DGGS", unit=" cells"):
        try:
            records.append(_isea3h_row_from_id(cell_id, fix_antimeridian))
        except Exception as e:
            print(f"Error generating ISEA3H DGGS cell {cell_id}: {e}")
            continue
    return gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")


def isea3h_grid_within_bbox(resolution, bbox, fix_antimeridian=None, compact=False):
    resolution = validate_isea3h_resolution(resolution)
    accuracy = ISEA3H_RES_ACCURACY_DICT.get(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    bounding_box = box(min_lon, min_lat, max_lon, max_lat)
    bounding_box_wkt = bounding_box.wkt
    shapes = isea3h_dggs.convert_shape_string_to_dggs_shapes(
        bounding_box_wkt, ShapeStringFormat.WKT, accuracy
    )
    shape = shapes[0]
    bbox_cells = shape.get_shape().get_outer_ring().get_cells()
    bounding_cell = isea3h_dggs.get_bounding_dggs_cell(bbox_cells)
    bounding_children_cells = get_isea3h_children_cells_within_bbox(
        bounding_cell.get_cell_id(), bounding_box, resolution
    )
    if bounding_children_cells:
        if compact:
            bounding_children_cells = isea3h_compact(bounding_children_cells)
        records = []
        for cell_id in bounding_children_cells:
            records.append(_isea3h_row_from_id(cell_id, fix_antimeridian))
        return gpd.GeoDataFrame(records, geometry="geometry", crs="EPSG:4326")


def isea3h_grid_ids(resolution, compact=False):
    """
    Return a list of ISEA3H cell IDs for the whole world at a given resolution.
    """
    resolution = validate_isea3h_resolution(resolution)
    cell_ids = get_isea3h_children_cells(ISEA3H_BASE_CELLS, resolution)
    if compact:
        cell_ids = isea3h_compact(cell_ids)
    return [str(cid) for cid in cell_ids]


def isea3h_grid_within_bbox_ids(resolution, bbox, compact=False):
    """
    Return a list of ISEA3H cell IDs intersecting the given bounding box at a given resolution.
    """
    resolution = validate_isea3h_resolution(resolution)
    accuracy = ISEA3H_RES_ACCURACY_DICT.get(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    bounding_box = box(min_lon, min_lat, max_lon, max_lat)
    bounding_box_wkt = bounding_box.wkt
    shapes = isea3h_dggs.convert_shape_string_to_dggs_shapes(
        bounding_box_wkt, ShapeStringFormat.WKT, accuracy
    )
    shape = shapes[0]
    bbox_cells = shape.get_shape().get_outer_ring().get_cells()
    bounding_cell = isea3h_dggs.get_bounding_dggs_cell(bbox_cells)
    bounding_children_cells = get_isea3h_children_cells_within_bbox(
        bounding_cell.get_cell_id(), bounding_box, resolution
    )
    cell_ids = list(bounding_children_cells or [])
    if compact:
        cell_ids = isea3h_compact(cell_ids)
    return cell_ids


def isea3hgrid(
    resolution, bbox=None, output_format="gpd", fix_antimeridian=None, compact=False
):
    """
    Generate ISEA3H grid for pure Python usage.

    Args:
        resolution (int): ISEA3H resolution [0..40]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output output_format ('geojson', 'csv', etc). Defaults to None (list of IDs).
        fix_antimeridian (str, optional): Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none
            Defaults to False when None or omitted.
        compact (bool, optional): Enable ISEA3H compact mode to reduce cell count.

    Returns:
        dict or list: GeoJSON FeatureCollection, file path, or list of IDs depending on output_format
    """
    # Allow running on all platforms

    if bbox is None:
        bbox = [-180, -90, 180, 90]
        total_cells = 20 * (7**resolution)
        if total_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {total_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
        gdf = isea3h_grid(
            resolution, fix_antimeridian=fix_antimeridian, compact=compact
        )
    else:
        gdf = isea3h_grid_within_bbox(
            resolution, bbox, fix_antimeridian=fix_antimeridian, compact=compact
        )

    output_name = f"isea3h_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def isea3hgrid_cli():
    parser = argparse.ArgumentParser(description="Generate OpenEaggr ISEA3H DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..40]"
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
        help="Enable ISEA3H compact mode to reduce cell count",
    )
    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    fix_antimeridian = args.fix_antimeridian
    if platform.system() == "Windows":
        try:
            result = isea3hgrid(
                resolution,
                bbox,
                args.output_format,
                fix_antimeridian=fix_antimeridian,
                compact=args.compact,
            )
            if args.output_format in STRUCTURED_FORMATS:
                print(result)
        except ValueError as e:
            print(f"Error: {str(e)}")
            return
    else:
        print("ISEA3H is only supported on Windows systems")


if __name__ == "__main__":
    isea3hgrid_cli()
