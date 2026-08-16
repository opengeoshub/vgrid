"""
Geohash Grid Generator Module

Generates Geohash DGGS grids for specified resolutions with automatic cell generation and validation using hierarchical geocoding system.

Key Functions:
- geohash_grid(): Main grid generation function for whole world
- geohash_grid_within_bbox(): Grid generation within bounding box
- geohashgrid(): User-facing function with multiple output formats
- geohashgrid_cli(): Command-line interface for grid generation

Reference: https://geohash.softeng.co/uekkn, https://github.com/vinsci/geohash, https://www.movable-type.co.uk/scripts/geohash.html?geohash=dp3
"""

import argparse
from shapely.geometry import Polygon
from tqdm import tqdm
from vgrid.utils.constants import (
    MAX_CELLS,
    INITIAL_GEOHASHES,
    OUTPUT_FORMATS,
    STRUCTURED_FORMATS,
)
from vgrid.utils.geometry import graticule_dggs_to_geoseries
import geopandas as gpd
from vgrid.conversion.dggs2geo.geohash2geo import geohash2geo
from vgrid.utils.io import (
    validate_bbox,
    validate_geohash_resolution,
    convert_to_output_format,
)
from vgrid.conversion.dggscompact.geohashcompact import (
    geohash_compact,
    get_geohash_resolution,
)


def _geohash_row_from_id(geohash_id):
    cell_polygon = geohash2geo(geohash_id)
    cell_resolution = get_geohash_resolution(geohash_id)
    return graticule_dggs_to_geoseries(
        "geohash", geohash_id, cell_resolution, cell_polygon
    )


def expand_geohash(gh, target_length, geohashes):
    if len(gh) == target_length:
        geohashes.add(gh)
        return
    for char in "0123456789bcdefghjkmnpqrstuvwxyz":
        expand_geohash(gh + char, target_length, geohashes)


def geohash_grid(resolution, compact=False):
    """Generate GeoJSON for the entire world at the given geohash resolution."""
    resolution = validate_geohash_resolution(resolution)
    geohashes = set()
    for gh in INITIAL_GEOHASHES:
        expand_geohash(gh, resolution, geohashes)

    geohash_ids = list(geohashes)
    if compact:
        geohash_ids = geohash_compact(geohash_ids)

    geohash_records = []
    for gh in tqdm(geohash_ids, desc="Generating Geohash DGGS", unit=" cells"):
        geohash_records.append(_geohash_row_from_id(gh))
    return gpd.GeoDataFrame(geohash_records, geometry="geometry", crs="EPSG:4326")


def expand_geohash_bbox(gh, target_length, geohashes, bbox_polygon):
    """Expand geohash only if it intersects the bounding box."""
    polygon = geohash2geo(gh)
    if not polygon.intersects(bbox_polygon):
        return

    if len(gh) == target_length:
        geohashes.add(gh)  # Add to the set if it reaches the target resolution
        return

    for char in "0123456789bcdefghjkmnpqrstuvwxyz":
        expand_geohash_bbox(gh + char, target_length, geohashes, bbox_polygon)


def geohash_grid_within_bbox(resolution, bbox, compact=False):
    """Generate GeoJSON for geohashes within a bounding box at the given resolution."""
    resolution = validate_geohash_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    geohash_records = []
    bbox_polygon = Polygon.from_bounds(min_lon, min_lat, max_lon, max_lat)
    intersected_geohashes = {
        gh for gh in INITIAL_GEOHASHES if geohash2geo(gh).intersects(bbox_polygon)
    }
    geohashes_bbox = set()
    for gh in intersected_geohashes:
        expand_geohash_bbox(gh, resolution, geohashes_bbox, bbox_polygon)
    geohash_ids = list(geohashes_bbox)
    if compact:
        geohash_ids = geohash_compact(geohash_ids)
    for gh in tqdm(geohash_ids, desc="Generating Geohash DGGS", unit=" cells"):
        geohash_records.append(_geohash_row_from_id(gh))
    return gpd.GeoDataFrame(geohash_records, geometry="geometry", crs="EPSG:4326")


def geohash_grid_ids(resolution, compact=False):
    """
    Return a list of Geohash IDs for the whole world at the given resolution.
    """
    resolution = validate_geohash_resolution(resolution)
    geohashes = set()
    for gh in INITIAL_GEOHASHES:
        expand_geohash(gh, resolution, geohashes)
    geohash_ids = list(geohashes)
    if compact:
        geohash_ids = geohash_compact(geohash_ids)
    return geohash_ids


def geohash_grid_within_bbox_ids(resolution, bbox, compact=False):
    """
    Return a list of Geohash IDs intersecting the given bounding box at the given resolution.
    """
    resolution = validate_geohash_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    bbox_polygon = Polygon.from_bounds(min_lon, min_lat, max_lon, max_lat)
    intersected_geohashes = {
        gh for gh in INITIAL_GEOHASHES if geohash2geo(gh).intersects(bbox_polygon)
    }
    geohashes_bbox = set()
    for gh in intersected_geohashes:
        expand_geohash_bbox(gh, resolution, geohashes_bbox, bbox_polygon)
    geohash_ids = list(geohashes_bbox)
    if compact:
        geohash_ids = geohash_compact(geohash_ids)
    return geohash_ids


def geohashgrid(resolution, bbox=None, output_format="gpd", compact=False):
    """
    Generate Geohash grid for pure Python usage.

    Args:
        resolution (int): Geohash resolution [1..10]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format ('geojson', 'csv', 'gpd', 'shapefile', 'gpkg', 'parquet', or None for list of Geohash IDs).
        compact (bool, optional): Enable Geohash compact mode to reduce cell count.

    Returns:
        dict, list, or str: Output in the requested format or file path.
    """
    if bbox is None:
        bbox = [-180, -90, 180, 90]
        total_cells = 32**resolution
        if total_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {total_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
        gdf = geohash_grid(resolution, compact=compact)
    else:
        gdf = geohash_grid_within_bbox(resolution, bbox, compact=compact)
    output_name = f"geohash_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def geohashgrid_cli():
    parser = argparse.ArgumentParser(description="Generate Geohash DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [1..10]"
    )
    parser.add_argument(
        "-b",
        "--bbox",
        type=float,
        nargs=4,
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)",
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
        help="Enable Geohash compact mode to reduce cell count",
    )
    args = parser.parse_args()
    try:
        result = geohashgrid(
            args.resolution, args.bbox, args.output_format, compact=args.compact
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    geohashgrid_cli()
