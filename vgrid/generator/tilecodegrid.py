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
import geopandas as gpd
from tqdm import tqdm
from vgrid.dggs import mercantile
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import (
    validate_bbox,
    validate_tilecode_resolution,
    convert_to_output_format,
    add_verbose_argument,
)
from vgrid.conversion.dggscompact.tilecodecompact import tilecode_compact
from vgrid.dggs.tilecode import tilecode_resolution
from vgrid.conversion.dggs2geo.tilecode2geo import tilecode2geo


def _tilecode_ids_for_bbox(resolution, bbox):
    min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    return [
        f"z{tile.z}x{tile.x}y{tile.y}"
        for tile in mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
    ]


def tilecode_grid(resolution, bbox, compact=False, verbose=True):
    resolution = validate_tilecode_resolution(resolution)
    tilecode_ids = _tilecode_ids_for_bbox(resolution, bbox)
    if compact:
        tilecode_ids = tilecode_compact(tilecode_ids, verbose=verbose)

    tilecode_records = []
    for tilecode_id in tqdm(
        tilecode_ids, desc="Generating Tilecode DGGS", unit=" cells", disable=not verbose
    ):
        cell_polygon = tilecode2geo(tilecode_id)
        if cell_polygon is None or cell_polygon.is_empty:
            continue
        cell_resolution = tilecode_resolution(tilecode_id)
        tilecode_record = graticule_dggs_to_geoseries(
            "tilecode", tilecode_id, cell_resolution, cell_polygon
        )
        tilecode_records.append(tilecode_record)

    return gpd.GeoDataFrame(tilecode_records, geometry="geometry", crs="EPSG:4326")


def tilecode_grid_ids(resolution, compact=False, verbose=True):
    """
    Return a list of Tilecode IDs for the whole world at the given resolution.
    """
    resolution = validate_tilecode_resolution(resolution)
    bbox = [-180.0, -85.05112878, 180.0, 85.05112878]
    tilecode_ids = _tilecode_ids_for_bbox(resolution, bbox)
    if compact:
        tilecode_ids = tilecode_compact(tilecode_ids, verbose=verbose)
    return tilecode_ids


def tilecode_grid_within_bbox_ids(resolution, bbox, compact=False, verbose=True):
    """
    Return a list of Tilecode IDs intersecting the given bounding box at the given resolution.
    """
    resolution = validate_tilecode_resolution(resolution)
    tilecode_ids = _tilecode_ids_for_bbox(resolution, bbox)
    if compact:
        tilecode_ids = tilecode_compact(tilecode_ids, verbose=verbose)
    return tilecode_ids


def tilecodegrid(resolution, bbox=None, output_format="gpd", compact=False, verbose=True):
    """
    Generate Tilecode grid for pure Python usage.

    Args:
        resolution (int): Tilecode resolution [0..26]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format ('geojson', 'csv', etc.). Defaults to None (list of Tilecode IDs).
        compact (bool, optional): Enable Tilecode compact mode to reduce cell count.

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

    gdf = tilecode_grid(resolution, bbox, compact=compact, verbose=verbose)

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
    parser.add_argument(
        "-c",
        "--compact",
        action="store_true",
        help="Enable Tilecode compact mode to reduce cell count",
    )
    add_verbose_argument(parser)
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
        result = tilecodegrid(
            resolution, bbox, args.output_format, compact=args.compact, verbose=args.verbose
        )
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    tilecodegrid_cli()
