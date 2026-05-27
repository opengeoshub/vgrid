"""
GARS Grid Generator Module

Generates GARS (Global Area Reference System) DGGS grids for specified resolutions with automatic cell generation and validation using military grid reference system.

Key Functions:
- gars_grid(): Main grid generation function for whole world
- gars_grid_within_bbox(): Grid generation within bounding box
- garsgrid(): User-facing function with multiple output formats
- garsgrid_cli(): Command-line interface for grid generation
"""

import argparse
from tqdm import tqdm
from shapely.geometry import Polygon
import numpy as np
from gars_field.garsgrid import GARSGrid  # Ensure the correct import path
from vgrid.utils.constants import MAX_CELLS, OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import (
    validate_gars_resolution,
    convert_to_output_format,
    gars_num_cells,
)
import geopandas as gpd
from vgrid.utils.constants import GARS_RESOLUTION_MINUTES
from vgrid.utils.io import validate_bbox


def gars_grid(resolution, bbox=None):
    resolution = validate_gars_resolution(resolution)
    # Default to the whole world if no bounding box is provided
    if bbox is None:
        min_lon, min_lat, max_lon, max_lat = -180, -90, 180, 90
    else:
        min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)

    resolution_minutes = GARS_RESOLUTION_MINUTES.get(resolution)
    resolution_degrees = resolution_minutes / 60.0

    # Generate ranges for longitudes and latitudes
    longitudes = np.arange(min_lon, max_lon, resolution_degrees)
    latitudes = np.arange(min_lat, max_lat, resolution_degrees)

    total_cells = len(longitudes) * len(latitudes)

    gars_records = []
    # Loop over longitudes and latitudes with tqdm progress bar
    with tqdm(total=total_cells, desc="Generating GARS DGGS", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                # Create the GARS grid code
                gars_cell = GARSGrid.from_latlon(lat, lon, resolution_minutes)
                wkt_polygon = gars_cell.polygon

                if wkt_polygon:
                    cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
                    gars_id = gars_cell.gars_id
                    gars_record = graticule_dggs_to_geoseries(
                        "gars", gars_id, resolution, cell_polygon
                    )
                    gars_records.append(gars_record)
                    pbar.update(1)

    # Create a FeatureCollection
    return gpd.GeoDataFrame(gars_records, geometry="geometry", crs="EPSG:4326")


def gars_grid_ids(resolution, bbox=None):
    """
    Return a list of GARS IDs for the whole world at the given resolution.
    """
    resolution = validate_gars_resolution(resolution)
    if bbox is None:
        min_lon, min_lat, max_lon, max_lat = -180, -90, 180, 90
    else:
        min_lon, min_lat, max_lon, max_lat = validate_bbox(bbox)
    resolution_minutes = GARS_RESOLUTION_MINUTES.get(resolution)
    resolution_degrees = resolution_minutes / 60.0

    longitudes = np.arange(min_lon, max_lon, resolution_degrees)
    latitudes = np.arange(min_lat, max_lat, resolution_degrees)

    total_cells = len(longitudes) * len(latitudes)
    ids = []
    with tqdm(total=total_cells, desc="Generating GARS IDs", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                cell = GARSGrid.from_latlon(lat, lon, resolution_minutes)
                ids.append(cell.gars_id)
                pbar.update(1)

    return ids


def garsgrid(resolution, bbox=None, output_format="gpd"):
    """
    Generate GARS grid for pure Python usage.

    Args:
        resolution (int): GARS resolution [1..4]
        bbox (list, optional): Bounding box [min_lon, min_lat, max_lon, max_lat]. Defaults to None (whole world).
        output_format (str, optional): Output format ('geojson', 'csv', etc.). Defaults to None (list of GARS IDs).

    Returns:
        dict, list, or str: Output depending on output_format
    """
    if bbox is None:
        resolution_minutes = GARS_RESOLUTION_MINUTES.get(resolution)
        total_cells = gars_num_cells(resolution)
        if total_cells > MAX_CELLS:
            raise ValueError(
                f"Resolution level {resolution} ({resolution_minutes} minutes) will generate {total_cells} cells which exceeds the limit of {MAX_CELLS}"
            )
        gdf = gars_grid(resolution)
    else:
        gdf = gars_grid(resolution, bbox)
    output_name = f"gars_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def garsgrid_cli():
    parser = argparse.ArgumentParser(description="Generate GARS DGGS")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        choices=[1, 2, 3, 4],
        required=True,
        help="Resolution level (1=30min, 2=15min, 3=5min, 4=1min)",
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
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    try:
        result = garsgrid(resolution, bbox, args.output_format)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    garsgrid_cli()
