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
from collections import deque

import geopandas as gpd
from tqdm import tqdm
from shapely.geometry import box
import a5
from a5.core.cell_info import get_num_cells
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import geodesic_dggs_to_geoseries
from vgrid.utils.io import (
    validate_a5_resolution,
    validate_bbox,
    convert_to_output_format,
    add_verbose_argument,
)
from vgrid.conversion.dggs2geo.a52geo import a52geo_u64
from vgrid.conversion.dggscompact.a5compact import a5_compact


def _a5_compact_cell_ids(cell_ids, verbose=True):
    """Compact A5 cell IDs via ``a5_compact``."""
    if not cell_ids:
        return []
    as_hex = isinstance(cell_ids[0], str)
    hexes = list(cell_ids) if as_hex else [a5.u64_to_hex(cell_id) for cell_id in cell_ids]
    compacted = a5_compact(hexes, verbose=verbose)
    if as_hex:
        return compacted
    return [a5.hex_to_u64(cell_id) for cell_id in compacted]


def a5_grid(resolution, bbox, options=None, split_antimeridian=False, compact=False, verbose=True):
    resolution = validate_a5_resolution(resolution)
    """
    Generate an A5 DGGS grid for a given resolution and bounding box.

    Args:
        resolution (int): A5 resolution [0..30]
        bbox (list): Bounding box [min_lon, min_lat, max_lon, max_lat]
        options (dict, optional): Options for a52geo.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
        compact (bool, optional): Enable A5 compact mode to reduce cell count.
    Returns:
        GeoDataFrame: A5 grid cells within the bounding box
    """
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
        cell_resolution = a5.get_resolution(seed_cell_id)
        num_edges = 5
        if cell_resolution == 1:
            num_edges = 3
        row = geodesic_dggs_to_geoseries(
            "a5",
            a5.u64_to_hex(seed_cell_id),
            cell_resolution,
            seed_cell_polygon,
            num_edges,
        )
        return gpd.GeoDataFrame([row], geometry="geometry", crs="EPSG:4326")

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

    cell_ids = list(intersecting_cells.keys())
    if compact:
        cell_ids = _a5_compact_cell_ids(cell_ids, verbose=verbose)

    a5_rows = []
    for cell_id in tqdm(cell_ids, desc="Generating A5 cells", unit=" cells", disable=not verbose):
        cell_polygon = intersecting_cells.get(cell_id)
        if cell_polygon is None or cell_polygon.is_empty:
            cell_polygon = a52geo_u64(
                cell_id, options=options, split_antimeridian=split_antimeridian
            )
        if cell_polygon is None or cell_polygon.is_empty:
            continue
        cell_resolution = a5.get_resolution(cell_id)
        if cell_resolution < 0:
            continue
        num_edges = 5
        if cell_resolution == 1:
            num_edges = 3
        row = geodesic_dggs_to_geoseries(
            "a5", a5.u64_to_hex(cell_id), cell_resolution, cell_polygon, num_edges
        )
        a5_rows.append(row)

    return gpd.GeoDataFrame(a5_rows, geometry="geometry", crs="EPSG:4326")


def a5_grid_ids(
    resolution, bbox=None, options=None, split_antimeridian=False, compact=False, verbose=True
):
    """
    Return A5 cell IDs (hex strings) for the same cells as `a5_grid`.

    Uses the centroid-seed BFS from `a5_grid`: expand neighbors until no cell
    intersecting the bbox remains to discover.

    Args:
        resolution (int): A5 resolution [0..30]
        bbox (list, optional): [min_lon, min_lat, max_lon, max_lat]. Defaults to world.
        options (dict, optional): Options for a52geo_u64.
        split_antimeridian (bool, optional): Passed to a52geo_u64.
        compact (bool, optional): Enable A5 compact mode to reduce cell count.

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

    cell_ids = list(intersecting_cells.keys())
    if compact:
        cell_ids = _a5_compact_cell_ids(cell_ids, verbose=verbose)
    return [
        a5.u64_to_hex(cell_id)
        for cell_id in cell_ids
        if a5.get_resolution(cell_id) >= 0
    ]


def a5grid(
    resolution,
    bbox=None,
    output_format="gpd",
    options=None,
    split_antimeridian=False,
    compact=False,
    verbose=True,
):
    """
    Generate A5 grid for pure Python usage.

    Args:
        resolution (int): A5 resolution [0..30]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format (gpd, gdf, geojson_dict/json_dict, geojson/json, csv, shp/shapefile, gpkg/geopackage, parquet/geoparquet, or None)
        options (dict, optional): Options for a52geo.
        split_antimeridian (bool, optional): When True, apply antimeridian fixing to the resulting polygons.
        compact (bool, optional): Enable A5 compact mode to reduce cell count.
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

    a5_gdf = a5_grid(
        resolution,
        bbox,
        options=options,
        split_antimeridian=split_antimeridian,
        compact=compact,
        verbose=verbose,
    )

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
        "-c",
        "--compact",
        action="store_true",
        help="Enable A5 compact mode to reduce cell count",
    )
    parser.add_argument(
        "-options",
        "--options",
        type=str,
        default=None,
        help="JSON string of options to pass to a52geo. "
        "Example: '{\"segments\": 1000}'",
    )
    add_verbose_argument(parser)
    args = parser.parse_args()

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
            split_antimeridian=args.split_antimeridian,
            compact=args.compact,
            verbose=args.verbose,
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    a5grid_cli()
