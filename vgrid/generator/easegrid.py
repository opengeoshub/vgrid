"""
EASE Grid Generator Module

Generates EASE (Equal-Area Scalable Earth) DGGS grids for specified resolutions with automatic cell generation and validation using equal-area projection system.

Key Functions:
- ease_grid(): Main grid generation function for whole world
- ease_grid_within_bbox(): Grid generation within bounding box
- easegrid(): User-facing function with multiple output formats
- easegrid_cli(): Command-line interface for grid generation
"""

import argparse
import geopandas as gpd
from shapely.geometry import box
from tqdm import tqdm
from ease_dggs.constants import grid_spec, ease_crs, geo_crs, levels_specs
from ease_dggs.dggs.grid_addressing import geo_polygon_to_grid_ids
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries, get_ease_resolution
from vgrid.utils.io import (
    validate_bbox,
    validate_ease_resolution,
    convert_to_output_format,
    add_verbose_argument,
)
from vgrid.conversion.dggscompact.easecompact import ease_compact
from vgrid.conversion.dggs2geo.ease2geo import ease2geo
# Initialize the geodetic model

geo_bounds = grid_spec["geo"]
min_longitude = geo_bounds["min_x"]
min_lattitude = geo_bounds["min_y"]
max_longitude = geo_bounds["max_x"]
max_latitude = geo_bounds["max_y"]


def get_ease_cells(resolution):
    """
    Generate a list of cell IDs based on the resolution, row, and column.
    """
    n_row = levels_specs[resolution]["n_row"]
    n_col = levels_specs[resolution]["n_col"]

    # Generate list of cell IDs
    cell_ids = []

    # Loop through all rows and columns at the specified resolution
    for row in range(n_row):
        for col in range(n_col):
            # Generate base ID (e.g., L0.RRRCCC for res=0)
            base_id = f"L{resolution}.{row:03d}{col:03d}"

            # Add additional ".RC" for each higher resolution
            cell_id = base_id
            for i in range(1, resolution + 1):
                cell_id += f".{row:1d}{col:1d}"  # For res=1: L0.RRRCCC.RC, res=2: L0.RRRCCC.RC.RC, etc.

            # Append the generated cell ID to the list
            cell_ids.append(cell_id)

    return cell_ids


def get_ease_cells_bbox(resolution, bbox):
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    bounding_box = box(min_lon, min_lat, max_lon, max_lat)
    bounding_box_wkt = bounding_box.wkt
    cells_bbox = geo_polygon_to_grid_ids(
        bounding_box_wkt,
        level=resolution,
        source_crs=geo_crs,
        target_crs=ease_crs,
        levels_specs=levels_specs,
        return_centroids=True,
        wkt_geom=True,
    )
    return cells_bbox


def _ease_row_from_id(cell_id):
    cell_polygon = ease2geo(cell_id)
    cell_resolution = get_ease_resolution(cell_id)
    return geodesic_dggs_to_geoseries(
        "ease", str(cell_id), cell_resolution, cell_polygon, 4
    )


def ease_grid_ids(resolution, compact=False):
    """
    Return a list of EASE-DGGS cell IDs for the whole world at a given resolution.

    Args:
        resolution (int): EASE resolution [0..6]

    Returns:
        list[str]: List of EASE cell IDs
    """
    resolution = validate_ease_resolution(resolution)
    cell_ids = get_ease_cells(resolution)
    if compact:
        cell_ids = ease_compact(cell_ids)
    return cell_ids


def ease_grid_within_bbox_ids(resolution, bbox, compact=False):
    """
    Return a list of EASE-DGGS cell IDs that intersect a bounding box.

    Args:
        resolution (int): EASE resolution [0..6]
        bbox (list[float]): [min_lon, min_lat, max_lon, max_lat]

    Returns:
        list[str]: List of EASE cell IDs intersecting the bbox
    """
    resolution = validate_ease_resolution(resolution)
    cells_result = get_ease_cells_bbox(resolution, validate_bbox(bbox))
    cells = (cells_result or {}).get("result", {}).get("data", [])
    if compact:
        cells = ease_compact(cells)
    return cells


def ease_grid(resolution, compact=False, verbose=True):
    resolution = validate_ease_resolution(resolution)
    cell_ids = ease_grid_ids(resolution, compact=compact)
    ease_rows = []
    for cell_id in tqdm(
        cell_ids, total=len(cell_ids), desc="Generating EASE DGGS", unit=" cells", disable=not verbose
    ):
        ease_rows.append(_ease_row_from_id(cell_id))
    return gpd.GeoDataFrame(ease_rows, geometry="geometry", crs="EPSG:4326")


def ease_grid_within_bbox(resolution, bbox, compact=False, verbose=True):
    resolution = validate_ease_resolution(resolution)
    cell_ids = ease_grid_within_bbox_ids(resolution, bbox, compact=compact)
    ease_rows = []
    for cell_id in tqdm(cell_ids, desc="Generating EASE DGGS", unit=" cells", disable=not verbose):
        ease_rows.append(_ease_row_from_id(cell_id))
    return gpd.GeoDataFrame(ease_rows, geometry="geometry", crs="EPSG:4326")


def easegrid(resolution, bbox=None, output_format="gpd", compact=False, verbose=True):
    if bbox is None:
        bbox = [min_longitude, min_lattitude, max_longitude, max_latitude]
        level_spec = levels_specs[resolution]
        n_row = level_spec["n_row"]
        n_col = level_spec["n_col"]
        total_cells = n_row * n_col
        if total_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {total_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
        gdf = ease_grid(resolution, compact=compact, verbose=verbose)
    else:
        gdf = ease_grid_within_bbox(resolution, bbox, compact=compact, verbose=verbose)
    output_name = f"ease_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def easegrid_cli():
    parser = argparse.ArgumentParser(description="Generate EASE-DGGS DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="resolution [0..6]"
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
        "-c",
        "--compact",
        action="store_true",
        help="Enable EASE compact mode to reduce cell count",
    )
    add_verbose_argument(parser)
    args = parser.parse_args()
    resolution = args.resolution
    bbox = (
        args.bbox
        if args.bbox
        else [min_longitude, min_lattitude, max_longitude, max_latitude]
    )
    try:
        result = easegrid(resolution, bbox, args.output_format, compact=args.compact, verbose=args.verbose)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    easegrid_cli()
