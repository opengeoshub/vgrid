"""
A5 Grid Generator Module

Generates A5 (Adaptive 5) DGGS grids for specified resolutions and bounding boxes with automatic cell generation and validation.

Key Functions:
- a5_grid(): Main grid generation function with bounding box support
- a5_grid_ids(): Same coverage as a5_grid, returns cell IDs only (hex strings)
- a5grid(): User-facing function with multiple output formats
- a5grid_cli(): Command-line interface for grid generation
"""

import argparse
import json
from typing import Optional

import geopandas as gpd
from tqdm import tqdm
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import validate_a5_resolution, validate_bbox, convert_to_output_format
from a5.core.cell_info import get_num_cells
from vgrid.conversion.latlon2dggs import latlon2a5
from vgrid.conversion.dggs2geo.a52geo import a52geo, a52geo_u64
import a5
from collections import deque
from pyproj import Geod
geod = Geod(ellps="WGS84")
from shapely.geometry import box
from shapely.ops import unary_union

# def a5_grid_old(resolution, bbox, options=None, split_antimeridian=False):
#     resolution = validate_a5_resolution(resolution)
#     """
#     Generate an A5 DGGS grid for a given resolution and bounding box.
#     Based on JavaScript logic that creates a regular grid and converts centroids to A5 cells.
    
#     Args:
#         resolution (int): A5 resolution [0..29]
#         bbox (list): Bounding box [min_lon, min_lat, max_lon, max_lat]
#         compact (bool, optional): When True, compact the cells.
#         options (dict, optional): Options for a52geo.
#         split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
#     Returns:
#         GeoDataFrame: A5 grid cells within the bounding box
#     """
#     min_lng, min_lat, max_lng, max_lat = bbox

#     # Calculate longitude and latitude width based on resolution
#     if resolution == 0:
#         lon_width = 35
#         lat_width = 35
#     elif resolution == 1:
#         lon_width = 18
#         lat_width = 18
#     elif resolution == 2:
#         lon_width = 10
#         lat_width = 10
#     elif resolution == 3:
#         lon_width = 5
#         lat_width = 5
#     elif resolution > 3:
#         base_width = 5  # at resolution 3
#         factor = 0.5 ** (resolution - 3)
#         lon_width = base_width * factor
#         lat_width = base_width * factor

#     # Generate longitude and latitude arrays
#     longitudes = []
#     latitudes = []

#     lon = min_lng
#     while lon < max_lng:
#         longitudes.append(lon)
#         lon += lon_width

#     lat = min_lat
#     while lat < max_lat:
#         latitudes.append(lat)
#         lat += lat_width

#     a5_rows = []
#     seen_a5_hex = set()  # Track unique A5 hex codes

#     # Generate features for each grid cell
#     total_cells = len(longitudes) * len(latitudes)
#     with tqdm(total=total_cells, desc="Generating A5 DGGS", unit=" cells") as pbar:
#         for lon in longitudes:
#             for lat in latitudes:
#                 min_lon = lon
#                 min_lat = lat
#                 max_lon = lon + lon_width
#                 max_lat = lat + lat_width

#                 # Calculate centroid
#                 centroid_lat = (min_lat + max_lat) / 2
#                 centroid_lon = (min_lon + max_lon) / 2

#                 try:
#                     # Convert centroid to A5 cell ID using direct A5 functions
#                     a5_hex = latlon2a5(centroid_lat, centroid_lon, resolution)
#                     cell_polygon = a52geo(
#                         a5_hex, options, split_antimeridian=split_antimeridian
#                     )

#                     if cell_polygon is not None:
#                         # Only add if this A5 hex code hasn't been seen before
#                         if a5_hex not in seen_a5_hex:
#                             seen_a5_hex.add(a5_hex)
#                             num_edges = 5
#                             if resolution == 1:
#                                 num_edges = 3   
#                             # Create row data
#                             row = geodesic_dggs_to_geoseries(
#                                 "a5", a5_hex, resolution, cell_polygon, num_edges
#                             )
#                             a5_rows.append(row)

#                 except Exception as e:
#                     # Skip cells that can't be processed
#                     print(
#                         f"Error processing cell at ({centroid_lon}, {centroid_lat}): {e}"
#                     )
#                 finally:
#                     pbar.update(1)

#     if not a5_rows:
#         raise ValueError(
#             "No A5 cells were generated. Check the input parameters and A5 library functions."
#         )
  
#     return gpd.GeoDataFrame(a5_rows, geometry="geometry", crs="EPSG:4326")

def a5_grid(resolution, bbox, options=None, split_antimeridian=False):
    resolution = validate_a5_resolution(resolution)
    """
    Generate an A5 DGGS grid for a given resolution and bounding box.
    Based on JavaScript logic that creates a regular grid and converts centroids to A5 cells.
    
    Args:
        resolution (int): A5 resolution [0..30]
        bbox (list): Bounding box [min_lon, min_lat, max_lon, max_lat]
        options (dict, optional): Options for a52geo.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
    Returns:
        GeoDataFrame: A5 grid cells within the bounding box
    """
    a5_rows = []    
    if bbox is None:
        min_lon, min_lat, max_lon, max_lat = -180, -90, 180, 90
    else:
        min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)

    bbox_polygon = box(min_lon, min_lat, max_lon, max_lat)
    bbox_center_lon = bbox_polygon.centroid.x
    bbox_center_lat = bbox_polygon.centroid.y
    a5_rows = []
    seed_cell_id = a5.lonlat_to_cell((bbox_center_lon, bbox_center_lat), resolution)
    seed_cell_polygon = a52geo_u64(seed_cell_id, options=options, split_antimeridian=split_antimeridian)
    cell_resolution = a5.get_resolution(seed_cell_id)      
    if seed_cell_polygon.contains(bbox_polygon):
        num_edges = 5  
        if cell_resolution == 1:
            num_edges = 3
        row = geodesic_dggs_to_geoseries(
            "a5", seed_cell_id, cell_resolution, seed_cell_polygon, num_edges
        )
        a5_rows.append(row)
        return gpd.GeoDataFrame(a5_rows, geometry="geometry", crs="EPSG:4326")

    # Store intersecting cells with their polygons
    intersecting_cells = {}  # {cell_id: (cell, polygon)}
    covered_cells = set()
    queue = deque([seed_cell_id])  # Use deque for BFS

    while queue:
        current_cell_id = queue.popleft()  # BFS: FIFO
        if current_cell_id in covered_cells:
            continue
        covered_cells.add(current_cell_id)

        # Convert polygon once
        cell_polygon = a52geo_u64(current_cell_id, options=options, split_antimeridian=split_antimeridian)

        # Only process if intersects
        if cell_polygon.intersects(bbox_polygon):
            # Store for later processing
            intersecting_cells[current_cell_id] = (cell_polygon)
            # Add neighbors to queue
            neighbors = a5.uncompact(a5.grid_disk_vertex(current_cell_id, 1), resolution)
            for neighbor_id in neighbors:
                if neighbor_id not in covered_cells:
                    queue.append(neighbor_id)

    # Process only intersecting cells (no double conversion)
    # Note: fix_antimeridian already applied when creating polygon in BFS loop    
    for cell_id, cell_polygon in tqdm(
        intersecting_cells.items(), desc="Generating A5 cells", unit=" cells"
    ):
        cell_resolution = a5.get_resolution(cell_id)
        num_edges = 5
        if cell_resolution == 1:
            num_edges = 3       
        row = geodesic_dggs_to_geoseries(
            "a5", a5.u64_to_hex(cell_id), resolution, cell_polygon, num_edges
        )
        a5_rows.append(row)


    return gpd.GeoDataFrame(a5_rows, geometry="geometry", crs="EPSG:4326")


def a5_grid_ids(resolution, bbox=None, options=None, split_antimeridian=False):
    """
    Return A5 cell IDs (hex strings) for the same cells as `a5_grid`.

    Uses the centroid-seed BFS from `a5_grid`: expand neighbors until no cell
    intersecting the bbox remains to discover.

    Args:
        resolution (int): A5 resolution [0..30]
        bbox (list, optional): [min_lon, min_lat, max_lon, max_lat]. Defaults to world.
        options (dict, optional): Options for a52geo_u64.
        split_antimeridian (bool, optional): Passed to a52geo_u64.

    Returns:
        list[str]: A5 cell IDs in hex form, in discovery order. Does not enforce MAX_CELLS.
    """
    resolution = validate_a5_resolution(resolution)
    if bbox is None:
        min_lon, min_lat, max_lon, max_lat = -180, -90, 180, 90
    else:
        min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)

    bbox_polygon = box(min_lon, min_lat, max_lon, max_lat)
    bbox_center_lon = bbox_polygon.centroid.x
    bbox_center_lat = bbox_polygon.centroid.y

    seed_cell_id = a5.lonlat_to_cell((bbox_center_lon, bbox_center_lat), resolution)
    seed_cell_polygon = a52geo_u64(
        seed_cell_id, options=options, split_antimeridian=split_antimeridian
    )
    if seed_cell_polygon.contains(bbox_polygon):
        return [a5.u64_to_hex(seed_cell_id)]

    intersecting_cells = {}
    covered_cells = set()
    queue = deque([seed_cell_id])

    while queue:
        current_cell_id = queue.popleft()
        if current_cell_id in covered_cells:
            continue
        covered_cells.add(current_cell_id)

        cell_polygon = a52geo_u64(
            current_cell_id, options=options, split_antimeridian=split_antimeridian
        )

        if cell_polygon.intersects(bbox_polygon):
            intersecting_cells[current_cell_id] = cell_polygon
            neighbors = a5.uncompact(
                a5.grid_disk_vertex(current_cell_id, 1), resolution
            )
            for neighbor_id in neighbors:
                if neighbor_id not in covered_cells:
                    queue.append(neighbor_id)

    return [a5.u64_to_hex(cid) for cid in intersecting_cells.keys()]


def a5_grid_ids_old(resolution, bbox):
    """
    Generate a list of unique A5 cell IDs intersecting the given bounding box.

    Note: Intentionally does not enforce MAX_CELLS limit for ID generation.

    Args:
        resolution (int): A5 resolution [0..29]
        bbox (list): [min_lon, min_lat, max_lon, max_lat]

    Returns:
        list[str]: List of A5 cell IDs
    """
    resolution = validate_a5_resolution(resolution)

    min_lng, min_lat, max_lng, max_lat = bbox

    if resolution == 0:
        lon_width = 35
        lat_width = 35
    elif resolution == 1:
        lon_width = 18
        lat_width = 18
    elif resolution == 2:
        lon_width = 10
        lat_width = 10
    elif resolution == 3:
        lon_width = 5
        lat_width = 5
    elif resolution > 3:
        base_width = 5
        factor = 0.5 ** (resolution - 3)
        lon_width = base_width * factor
        lat_width = base_width * factor

    longitudes = []
    latitudes = []

    lon = min_lng
    while lon < max_lng:
        longitudes.append(lon)
        lon += lon_width

    lat = min_lat
    while lat < max_lat:
        latitudes.append(lat)
        lat += lat_width

    seen_ids = set()
    ids = []
    total_cells = len(longitudes) * len(latitudes)
    with tqdm(total=total_cells, desc="Generating A5 IDs", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                min_lon = lon
                min_lat = lat
                max_lon = lon + lon_width
                max_lat = lat + lat_width

                centroid_lat = (min_lat + max_lat) / 2
                centroid_lon = (min_lon + max_lon) / 2

                try:
                    a5_hex = latlon2a5(centroid_lat, centroid_lon, resolution)
                    if a5_hex and a5_hex not in seen_ids:
                        seen_ids.add(a5_hex)
                        ids.append(a5_hex)
                except Exception:
                    pass
                finally:
                    pbar.update(1)

    return ids


def a5grid(
    resolution, bbox=None, output_format="gpd", options=None, split_antimeridian=False
):
    """
    Generate A5 grid for pure Python usage.

    Args:
        resolution (int): A5 resolution [0..30]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format (gpd, gdf, geojson_dict/json_dict, geojson/json, csv, shp/shapefile, gpkg/geopackage, parquet/geoparquet, or None)     
        options (dict, optional): Options for a52geo.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
    Returns:
        Depends on output_format. If None, returns a GeoDataFrame (gpd)
    """
    if bbox is None:
        bbox = [-180, -90, 180, 90]
        num_cells = get_num_cells(resolution)
        if num_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {num_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
    
    a5_gdf = a5_grid(resolution, bbox, options=options, split_antimeridian=split_antimeridian)


    output_name = f"a5_grid_{resolution}"
    return convert_to_output_format(a5_gdf, output_format, output_name)


def a5grid_cli():
    """CLI interface for generating A5 DGGS."""
    parser = argparse.ArgumentParser(description="Generate A5 DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..29]"
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
        "-split",
        "--split_antimeridian",
        action="store_true",
        default=False,
        help="Apply antimeridian fixing to the resulting polygons",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to a52geo. "
             "Example: '{\"segments\": 1000}'",
    )
    args = parser.parse_args()
    
    # Parse options JSON if provided
    options = None
    if args.options:
        try:
            options = json.loads(args.options)
        except json.JSONDecodeError as e:
            print(f"Error: Invalid JSON in options: {str(e)}")
            return
    
    try:
        result = a5grid(
            resolution=args.resolution, 
            bbox=args.bbox, 
            output_format=args.output_format, 
            options=options,
            split_antimeridian=args.split_antimeridian
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    a5grid_cli()
