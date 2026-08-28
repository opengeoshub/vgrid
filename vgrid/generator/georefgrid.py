"""
GEOREF Grid Generator Module

Generates GEOREF DGGS grids for specified resolutions with automatic cell generation and validation using World Geographic Reference System.

Key Functions:
- georef_grid(): Main grid generation function with bounding box support
- georef_grid_ids(): Returns list of GEOREF IDs for given bbox and resolution
- georefgrid(): User-facing function with multiple output formats
- georefgrid_cli(): Command-line interface for grid generation
"""

import argparse
from tqdm import tqdm
import numpy as np
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import graticule_dggs_to_geoseries
import geopandas as gpd
from vgrid.utils.io import (
    validate_bbox,
    validate_georef_resolution,
    convert_to_output_format,
    add_verbose_argument,
)
from vgrid.utils.constants import GEOREF_RESOLUTION_DEGREES
from vgrid.conversion.latlon2dggs import latlon2georef
from vgrid.conversion.dggs2geo.georef2geo import georef2geo


def georef_grid(resolution, bbox=None, verbose=True):
    resolution = validate_georef_resolution(resolution)
    if bbox is None:
        min_lon, min_lat, max_lon, max_lat = -180, -90, 180, 90
    else:
        min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    resolution_degrees = GEOREF_RESOLUTION_DEGREES.get(resolution)
    longitudes = np.arange(min_lon, max_lon, resolution_degrees)
    latitudes = np.arange(min_lat, max_lat, resolution_degrees)
    num_cells = len(longitudes) * len(latitudes)

    georef_records = []

    with tqdm(total=num_cells, desc="Generating GEOREF DGGS", unit=" cells", disable=not verbose) as pbar:
        for lon in longitudes:
            for lat in latitudes:
                georef_id = latlon2georef(lat, lon, resolution)
                cell_polygon = georef2geo(georef_id)
                georef_record = graticule_dggs_to_geoseries(
                    "georef", georef_id, resolution, cell_polygon
                )
                georef_records.append(georef_record)
                pbar.update(1)

    return gpd.GeoDataFrame(georef_records, geometry="geometry", crs="EPSG:4326")


def georef_grid_ids(resolution, bbox=None, verbose=True):
    """
    Return a list of GEOREF IDs for a given bounding box at the specified resolution.

    Args:
        bbox (list[float]): [min_lon, min_lat, max_lon, max_lat]
        resolution (int): GEOREF resolution [0..4]

    Returns:
        list[str]: List of GEOREF IDs
    """
    resolution = validate_georef_resolution(resolution)
    if bbox is None:
        min_lon, min_lat, max_lon, max_lat = -180, -90, 180, 90
    else:
        min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    resolution_degrees = GEOREF_RESOLUTION_DEGREES.get(resolution)
    longitudes = np.arange(min_lon, max_lon, resolution_degrees)
    latitudes = np.arange(min_lat, max_lat, resolution_degrees)

    num_cells = len(longitudes) * len(latitudes)
    ids = []
    with tqdm(total=num_cells, desc="Generating GEOREF IDs", unit=" cells", disable=not verbose) as pbar:
        for lon in longitudes:
            for lat in latitudes:
                georef_id = latlon2georef(lat, lon, resolution)
                ids.append(georef_id)
                pbar.update(1)

    return ids


def georefgrid(resolution, bbox=None, output_format="gpd", verbose=True):
    """
    Generate GEOREF grid for pure Python usage.

    Args:
        resolution (int): GEOREF resolution [0..4]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format ('geojson', 'csv', 'geo', 'gpd', 'shapefile', 'gpkg', 'parquet', or None for list of GEOREF IDs).

    Returns:
        dict, list, or str: Output in the requested format or file path.
    """
    if bbox is None:
        bbox = [-180, -90, 180, 90]
    gdf = georef_grid(resolution, bbox, verbose=verbose)
    output_name = f"georef_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def georefgrid_cli():
    parser = argparse.ArgumentParser(description="Generate GEOREF DGGS")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True, help="Resolution [0..4]"
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
    add_verbose_argument(parser)
    args = parser.parse_args()
    try:
        result = georefgrid(args.resolution, args.bbox, args.output_format, verbose=args.verbose)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    georefgrid_cli()
