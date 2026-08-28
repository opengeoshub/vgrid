"""
rHEALPix Grid Generator Module

Generates rHEALPix DGGS grids for specified resolutions with automatic cell generation and validation using hierarchical equal-area grid system.

Key Functions:
- rhealpix_grid(): Main grid generation function for whole world
- rhealpix_grid_within_bbox(): Grid generation within bounding box
- rhealpixgrid(): User-facing function with multiple output formats
- rhealpixgrid_cli(): Command-line interface for grid generation
"""

import argparse
from vgrid.dggs.rhealpixdggs.dggs import RHEALPixDGGS
from shapely.geometry import box
from tqdm import tqdm
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    validate_bbox,
    validate_rhealpix_resolution,
    convert_to_output_format,
    add_verbose_argument,
)
from vgrid.conversion.dggs2geo.rhealpix2geo import rhealpix2geo
from vgrid.conversion.dggscompact.rhealpixcompact import rhealpix_compact
from collections import deque

from pyproj import Geod
import geopandas as gpd

geod = Geod(ellps="WGS84")
rhealpix_dggs = RHEALPixDGGS()


def _rhealpix_row_from_cell_id(cell_id, fix_antimeridian=None):
    cell_polygon = rhealpix2geo(cell_id, fix_antimeridian=fix_antimeridian)
    rhealpix_uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
    rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
    cell_resolution = rhealpix_cell.resolution
    num_edges = 4
    if rhealpix_cell.ellipsoidal_shape() == "dart":
        num_edges = 3
    return geodesic_dggs_to_geoseries(
        "rhealpix", cell_id, cell_resolution, cell_polygon, num_edges
    )


def rhealpix_grid(resolution, fix_antimeridian=None, compact=False, verbose=True):
    resolution = validate_rhealpix_resolution(resolution)
    cell_ids = [str(rhealpix_cell) for rhealpix_cell in rhealpix_dggs.grid(resolution)]
    if compact:
        cell_ids = rhealpix_compact(cell_ids)

    rhealpix_rows = []
    for cell_id in tqdm(cell_ids, desc="Generating rHEALPix DGGS", unit=" cells", disable=not verbose):
        rhealpix_rows.append(_rhealpix_row_from_cell_id(cell_id, fix_antimeridian))
    return gpd.GeoDataFrame(rhealpix_rows, geometry="geometry", crs="EPSG:4326")


def rhealpix_grid_within_bbox(resolution, bbox, fix_antimeridian=None, compact=False, verbose=True):
    resolution = validate_rhealpix_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    bbox_polygon = box(min_lon, min_lat, max_lon, max_lat)
    bbox_center_lon = bbox_polygon.centroid.x
    bbox_center_lat = bbox_polygon.centroid.y
    seed_point = (bbox_center_lon, bbox_center_lat)
    rhealpix_rows = []
    seed_cell = rhealpix_dggs.cell_from_point(resolution, seed_point, plane=False)
    seed_cell_id = str(seed_cell)
    seed_cell_polygon = rhealpix2geo(seed_cell_id, fix_antimeridian=fix_antimeridian)
    if seed_cell_polygon.contains(bbox_polygon):
        num_edges = 4
        if seed_cell.ellipsoidal_shape() == "dart":
            num_edges = 3
        row = geodesic_dggs_to_geoseries(
            "rhealpix", seed_cell_id, resolution, seed_cell_polygon, num_edges
        )
        rhealpix_rows.append(row)
        return gpd.GeoDataFrame(rhealpix_rows, geometry="geometry", crs="EPSG:4326")

    # Store intersecting cells with their polygons
    intersecting_cells = {}  # {cell_id: (cell, polygon)}
    covered_cells = set()
    queue = deque([seed_cell])  # Use deque for BFS

    while queue:
        current_cell = queue.popleft()  # BFS: FIFO
        current_cell_id = str(current_cell)
        if current_cell_id in covered_cells:
            continue
        covered_cells.add(current_cell_id)

        # Convert polygon once
        cell_polygon = rhealpix2geo(current_cell_id, fix_antimeridian=fix_antimeridian)

        # Only process if intersects
        if cell_polygon.intersects(bbox_polygon):
            # Store for later processing
            intersecting_cells[current_cell_id] = (current_cell, cell_polygon)

            # Add neighbors to queue
            neighbors = current_cell.neighbors(plane=False)
            for _, neighbor in neighbors.items():
                neighbor_id = str(neighbor)
                if neighbor_id not in covered_cells:
                    queue.append(neighbor)

    cell_ids = list(intersecting_cells.keys())
    if compact:
        cell_ids = rhealpix_compact(cell_ids)

    for cell_id in tqdm(cell_ids, desc="Generating rHEALPix DGGS", unit=" cells", disable=not verbose):
        rhealpix_rows.append(_rhealpix_row_from_cell_id(cell_id, fix_antimeridian))

    return gpd.GeoDataFrame(rhealpix_rows, geometry="geometry", crs="EPSG:4326")


def rhealpix_grid_ids(resolution, compact=False):
    """
    Return a list of rHEALPix cell IDs for the whole world at a given resolution.
    """
    resolution = validate_rhealpix_resolution(resolution)
    cell_ids = [str(rhealpix_cell) for rhealpix_cell in rhealpix_dggs.grid(resolution)]
    if compact:
        cell_ids = rhealpix_compact(cell_ids)
    return cell_ids


def rhealpix_grid_within_bbox_ids(resolution, bbox, compact=False, verbose=True):
    """
    Return a list of rHEALPix cell IDs intersecting the given bounding box at a given resolution.
    """
    gdf = rhealpix_grid_within_bbox(resolution, bbox, compact=compact, verbose=verbose)
    if gdf.empty:
        return []
    return gdf["rhealpix"].tolist()


# Remove convert_rhealpixgrid_output_format and handle output logic in rhealpixgrid


def rhealpixgrid(
    resolution,
    bbox=None,
    output_format="gpd",
    fix_antimeridian=None,
    compact=False,
    verbose=True,
):
    """
    Generate rHEALPix grid for pure Python usage.

    Args:
        resolution (int): rHEALPix resolution [0..15]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output output_format ('geojson', 'csv', 'geo', 'gpd', 'shapefile', 'gpkg', 'parquet', or None for list of rHEALPix IDs). Defaults to None.
        fix_antimeridian (Antimeridian fixing method: shift, shift_balanced, shift_west, shift_east, split, none, optional): When True, apply antimeridian fixing to the resulting polygons.
            Defaults to False when None or omitted.
        compact (bool, optional): Enable rHEALPix compact mode to reduce cell count.

    Returns:
        dict, list, or str: Output in the requested output_format (GeoJSON FeatureCollection, list of IDs, file path, etc.)
    """
    if bbox is None:
        bbox = [-180, -90, 180, 90]
        num_cells = rhealpix_dggs.num_cells(resolution)
        if num_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {num_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
        gdf = rhealpix_grid(
            resolution, fix_antimeridian=fix_antimeridian, compact=compact, verbose=verbose
        )
    else:
        gdf = rhealpix_grid_within_bbox(
            resolution,
            bbox,
            fix_antimeridian=fix_antimeridian,
            compact=compact,
            verbose=verbose,
        )
    output_name = f"rhealpix_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def rhealpixgrid_cli():
    """CLI interface for generating rHEALPix grid."""
    parser = argparse.ArgumentParser(description="Generate rHEALPix DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..15]"
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
        help="Enable rHEALPix compact mode to reduce cell count",
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

    add_verbose_argument(parser)
    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    output_format = args.output_format
    fix_antimeridian = args.fix_antimeridian
    try:
        result = rhealpixgrid(
            resolution,
            bbox,
            output_format,
            fix_antimeridian=fix_antimeridian,
            compact=args.compact,
            verbose=args.verbose,
        )
        if output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    rhealpixgrid_cli()
