"""
Quadkey Grid Generator Module

Generates Quadkey DGGS grids for specified resolutions with automatic cell generation and validation using hierarchical geospatial indexing system.

Key Functions:
- quadkey_grid(): Main grid generation function with bounding box support
- quadkey_grid_resample(): Grid generation within GeoJSON features
- quadkeygrid(): User-facing function with multiple output formats
- quadkeygrid_cli(): Command-line interface for grid generation
"""

import argparse
import geopandas as gpd
from shapely.geometry import shape, Polygon
from shapely.ops import unary_union
from tqdm import tqdm
from vgrid.dggs import mercantile
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import validate_bbox, validate_quadkey_resolution, convert_to_output_format
from pyproj import Geod

geod = Geod(ellps="WGS84")


def quadkey_grid(resolution, bbox):
    resolution = validate_quadkey_resolution(resolution)
    quadkey_records = []
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)

    for tile in tqdm(tiles, desc="Generating Quadkey DGGS", unit=" cells"):
        z, x, y = tile.z, tile.x, tile.y
        bounds = mercantile.bounds(x, y, z)
        if bounds:
            min_lat, min_lon = bounds.south, bounds.west
            max_lat, max_lon = bounds.north, bounds.east
            quadkey_id = mercantile.quadkey(tile)
            cell_polygon = Polygon(
                [
                    [min_lon, min_lat],
                    [max_lon, min_lat],
                    [max_lon, max_lat],
                    [min_lon, max_lat],
                    [min_lon, min_lat],
                ]
            )
            quadkey_record = graticule_dggs_to_geoseries(
                "quadkey", quadkey_id, resolution, cell_polygon
            )
            quadkey_records.append(quadkey_record)
    return gpd.GeoDataFrame(quadkey_records, geometry="geometry", crs="EPSG:4326")


def quadkey_grid_ids(resolution):
    """
    Return a list of Quadkey IDs for the whole world at the given resolution.
    """
    resolution = validate_quadkey_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = [-180.0, -85.05112878, 180.0, 85.05112878]
    tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
    ids = []
    for tile in tqdm(tiles, desc="Generating Quadkey IDs", unit=" cells"):
        ids.append(mercantile.quadkey(tile))
    return ids


def quadkey_grid_within_bbox_ids(resolution, bbox):
    """
    Return a list of Quadkey IDs intersecting the given bounding box at the given resolution.
    """
    resolution = validate_quadkey_resolution(resolution)
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
    ids = []
    for tile in tqdm(tiles, desc="Generating Quadkey IDs", unit=" cells"):
        ids.append(mercantile.quadkey(tile))
    return ids


def quadkeygrid(resolution, bbox=None, output_format="gpd"):
    """
    Generate Quadkey grid for pure Python usage.

    Args:
        resolution (int): Quadkey resolution [0..26]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format ('geojson', 'csv', etc.). Defaults to None (list of Quadkey IDs).

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
        gdf = quadkey_grid(resolution, bbox)
    else:
        gdf = quadkey_grid(resolution, bbox)

    output_name = f"quadkey_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def quadkeygrid_cli():
    parser = argparse.ArgumentParser(description="Generate Quadkey DGGS.")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="resolution [0..26]"
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
    gdf = quadkey_grid(resolution, bbox)
    try:
        result = convert_to_output_format(
            gdf, args.output_format, f"quadkey_grid_{resolution}"
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    quadkeygrid_cli()
