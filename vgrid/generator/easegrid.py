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
from shapely.geometry import Polygon, box
from tqdm import tqdm
from ease_dggs.constants import grid_spec, ease_crs, geo_crs, levels_specs
from ease_dggs.dggs.grid_addressing import (
    grid_ids_to_geos,
    geo_polygon_to_grid_ids,
)
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import validate_bbox, validate_ease_resolution, convert_to_output_format
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


def ease_grid_ids(resolution):
    """
    Return a list of EASE-DGGS cell IDs for the whole world at a given resolution.

    Args:
        resolution (int): EASE resolution [0..6]

    Returns:
        list[str]: List of EASE cell IDs
    """
    resolution = validate_ease_resolution(resolution)
    return get_ease_cells(resolution)


def ease_grid_within_bbox_ids(resolution, bbox):
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
    return cells


def ease_grid(resolution):
    resolution = validate_ease_resolution(resolution)
    ease_rows = []
    level_spec = levels_specs[resolution]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]
    cells = get_ease_cells(resolution)
    for cell in tqdm(
        cells, total=len(cells), desc="Generating EASE DGGS", unit=" cells"
    ):
        geo = grid_ids_to_geos([cell])
        center_lon, center_lat = geo["result"]["data"][0]
        cell_min_lat = center_lat - (180 / (2 * n_row))
        cell_max_lat = center_lat + (180 / (2 * n_row))
        cell_min_lon = center_lon - (360 / (2 * n_col))
        cell_max_lon = center_lon + (360 / (2 * n_col))
        cell_polygon = Polygon(
            [
                [cell_min_lon, cell_min_lat],
                [cell_max_lon, cell_min_lat],
                [cell_max_lon, cell_max_lat],
                [cell_min_lon, cell_max_lat],
                [cell_min_lon, cell_min_lat],
            ]
        )
        if cell_polygon:
            num_edges = 4
            row = geodesic_dggs_to_geoseries(
                "ease", str(cell), resolution, cell_polygon, num_edges
            )
            ease_rows.append(row)
    return gpd.GeoDataFrame(ease_rows, geometry="geometry", crs="EPSG:4326")


def ease_grid_within_bbox(resolution, bbox):
    resolution = validate_ease_resolution(resolution)
    ease_rows = []
    level_spec = levels_specs[resolution]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]
    cells = get_ease_cells_bbox(resolution, bbox)["result"]["data"]
    if cells:
        for cell in tqdm(cells, desc="Generating EASE DGGS", unit=" cells"):
            geo = grid_ids_to_geos([cell])
            if geo:
                center_lon, center_lat = geo["result"]["data"][0]
                cell_min_lat = center_lat - (180 / (2 * n_row))
                cell_max_lat = center_lat + (180 / (2 * n_row))
                cell_min_lon = center_lon - (360 / (2 * n_col))
                cell_max_lon = center_lon + (360 / (2 * n_col))
                cell_polygon = Polygon(
                    [
                        [cell_min_lon, cell_min_lat],
                        [cell_max_lon, cell_min_lat],
                        [cell_max_lon, cell_max_lat],
                        [cell_min_lon, cell_max_lat],
                        [cell_min_lon, cell_min_lat],
                    ]
                )
                num_edges = 4
                row = geodesic_dggs_to_geoseries(
                    "ease", str(cell), resolution, cell_polygon, num_edges
                )
                ease_rows.append(row)
    return gpd.GeoDataFrame(ease_rows, geometry="geometry", crs="EPSG:4326")


def easegrid(resolution, bbox=None, output_format="gpd"):
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
        gdf = ease_grid(resolution)
    else:
        gdf = ease_grid_within_bbox(resolution, bbox)
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
    args = parser.parse_args()
    resolution = args.resolution
    bbox = (
        args.bbox
        if args.bbox
        else [min_longitude, min_lattitude, max_longitude, max_latitude]
    )
    try:
        result = easegrid(resolution, bbox, args.output_format)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    easegrid_cli()
