"""
Tilecode Grid Generator Module

Generates Tilecode DGGS grids for specified resolutions with automatic cell generation and validation using hierarchical geospatial indexing system.

Key Functions:
- tilecode_grid(): Main grid generation function with bounding box support
- tilecode_grid_resample(): Grid generation within GeoJSON features
- tilecodegrid(): User-facing function with multiple output formats
- tilecodegrid_cli(): Command-line interface for grid generation
"""

import argparse
from shapely.geometry import shape, Polygon
import geopandas as gpd
from tqdm import tqdm
from shapely.ops import unary_union
from vgrid.dggs import mercantile
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import validate_bbox, validate_tilecode_resolution, convert_to_output_format


def tilecode_grid(resolution, bbox):
    resolution = validate_tilecode_resolution(resolution)
    tilecode_records = []
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
    for tile in tqdm(tiles, desc="Generating Tilecode DGGS", unit=" cells"):
        z, x, y = tile.z, tile.x, tile.y
        tilecode_id = f"z{tile.z}x{tile.x}y{tile.y}"
        bounds = mercantile.bounds(x, y, z)
        if bounds:
            # Create the bounding box coordinates for the polygon
            min_lat, min_lon = bounds.south, bounds.west
            max_lat, max_lon = bounds.north, bounds.east

            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],  # Bottom-left corner
                    [max_lon, min_lat],  # Bottom-right corner
                    [max_lon, max_lat],  # Top-right corner
                    [min_lon, max_lat],  # Top-left corner
                    [min_lon, min_lat],  # Closing the polygon (same as the first point)
                ]
            )
            tilecode_record = graticule_dggs_to_geoseries(
                "tilecode", tilecode_id, resolution, cell_polygon
            )
            tilecode_records.append(tilecode_record)

    return gpd.GeoDataFrame(tilecode_records, geometry="geometry", crs="EPSG:4326")


def tilecode_grid_ids(resolution):
    """
    Return a list of Tilecode IDs for the whole world at the given resolution.
    """
    resolution = validate_tilecode_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = [-180.0, -85.05112878, 180.0, 85.05112878]
    tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
    ids = []
    for tile in tqdm(tiles, desc="Generating Tilecode IDs", unit=" cells"):
        z, x, y = tile.z, tile.x, tile.y
        ids.append(f"z{z}x{x}y{y}")
    return ids


def tilecode_grid_within_bbox_ids(resolution, bbox):
    """
    Return a list of Tilecode IDs intersecting the given bounding box at the given resolution.
    """
    resolution = validate_tilecode_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
    ids = []
    for tile in tqdm(tiles, desc="Generating Tilecode IDs", unit=" cells"):
        z, x, y = tile.z, tile.x, tile.y
        ids.append(f"z{z}x{x}y{y}")
    return ids


def tilecodegrid(resolution, bbox=None, output_format="gpd"):
    """
    Generate Tilecode grid for pure Python usage.

    Args:
        resolution (int): Tilecode resolution [0..26]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format ('geojson', 'csv', etc.). Defaults to None (list of Tilecode IDs).

    Returns:
        dict, list, or str: Output depending on output_format
    """
    if bbox is None:
        bbox = [-180.0, -85.05112878, 180.0, 85.05112878]
        num_cells = 4**resolution
        if num_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution {resolution} will generate {num_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
        gdf = tilecode_grid(resolution, bbox)
    else:
        gdf = tilecode_grid(resolution, bbox)

    output_name = f"tilecode_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def tilecodegrid_cli():
    parser = argparse.ArgumentParser(description="Generate Tilecode DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="resolution [0..29]"
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
    bbox = args.bbox if args.bbox else [-180.0, -85.05112878, 180.0, 85.05112878]

    if bbox == [-180.0, -85.05112878, 180.0, 85.05112878]:
        num_cells = 4**resolution
        if num_cells > MAX_CELLS:
            print(
                f"Resolution {resolution} will generate {num_cells} cells "
                f"which exceeds the limit of {MAX_CELLS}."
            )
            print("Please select a smaller resolution and try again.")
            return
    try:
        result = tilecodegrid(resolution, bbox, args.output_format)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    tilecodegrid_cli()
