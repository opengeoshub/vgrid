"""
MGRS Grid Generator Module

Generates MGRS (Military Grid Reference System) DGGS grids for specified resolutions with automatic cell generation and validation using NATO military coordinate system.

Key Functions:
- mgrs_grid(): Main grid generation function with GZD support
- mgrsgrid(): User-facing function with multiple output formats
- mgrsgrid_cli(): Command-line interface for grid generation
"""

import argparse
import json
import re
import os
import numpy as np
import geopandas as gpd
from shapely.geometry import shape, Polygon
from shapely.ops import transform
from pyproj import CRS, Transformer
from tqdm import tqdm
from vgrid.utils.geometry import graticule_dggs_to_geoseries
from vgrid.utils.io import validate_mgrs_resolution, convert_to_output_format, add_verbose_argument
from vgrid.utils.constants import OUTPUT_FORMATS, STRUCTURED_FORMATS
from vgrid.dggs import mgrs


def is_valid_gzd(gzd):
    """Check if a Grid Zone Designator (GZD) is valid."""
    pattern = r"^(?:0[1-9]|[1-5][0-9]|60)[C-HJ-NP-X]$"
    return bool(re.match(pattern, gzd))


def mgrs_grid(gzd, resolution, verbose=True):
    resolution = validate_mgrs_resolution(resolution)
    # Reference: https://www.maptools.com/tutorials/utm/details
    cell_size = 100_000 // (10**resolution)
    north_bands = "NPQRSTUVWX"
    south_bands = "MLKJHGFEDC"
    band_distance = 111_132 * 8
    gzd_band = gzd[2]

    if gzd_band >= "N":  # North Hemesphere
        epsg_code = int("326" + gzd[:2])
        min_x, min_y, max_x, max_y = 100000, 0, 900000, 9500000  # for the North
        north_band_idx = north_bands.index(gzd_band)
        max_y = band_distance * (north_band_idx + 1)
        if gzd_band == "X":
            max_y += band_distance  # band X = 12 deggrees instead of 8 degrees

    else:  # South Hemesphere
        epsg_code = int("327" + gzd[:2])
        min_x, min_y, max_x, max_y = 100000, 0, 900000, 10000000  # for the South
        south_band_idx = south_bands.index(gzd_band)
        max_y = band_distance * (south_band_idx + 1)

    utm_crs = CRS.from_epsg(epsg_code)
    wgs84_crs = CRS.from_epsg(4326)
    transformer = Transformer.from_crs(utm_crs, wgs84_crs, always_xy=True).transform

    gzd_json_path = os.path.join(os.path.dirname(__file__), "gzd.geojson")
    with open(gzd_json_path, encoding="utf-8") as f:
        gzd_data = json.load(f)

    gzd_features = gzd_data["features"]
    gzd_feature = [
        feature for feature in gzd_features if feature["properties"].get("gzd") == gzd
    ][0]
    gzd_geom = shape(gzd_feature["geometry"])

    # Create grid polygons
    mgrs_records = []
    x_coords = np.arange(min_x, max_x, cell_size)
    y_coords = np.arange(min_y, max_y, cell_size)
    num_cells = len(x_coords) * len(y_coords)
    with tqdm(total=num_cells, desc="Generating MGRS DGGS", unit=" cells", disable=not verbose) as pbar:
        for x in x_coords:
            for y in y_coords:
                cell_polygon_utm = Polygon(
                    [
                        (x, y),
                        (x + cell_size, y),
                        (x + cell_size, y + cell_size),
                        (x, y + cell_size),
                        (x, y),  # Close the polygon
                    ]
                )
                cell_polygon = transform(transformer, cell_polygon_utm)

                if cell_polygon.intersects(gzd_geom):
                    centroid_lat, centroid_lon = (
                        cell_polygon.centroid.y,
                        cell_polygon.centroid.x,
                    )
                    mgrs_id = mgrs.toMgrs(centroid_lat, centroid_lon, resolution)
                    mgrs_record = graticule_dggs_to_geoseries(
                        "mgrs", mgrs_id, resolution, cell_polygon
                    )
                    # clip inside GZD:
                    if not gzd_geom.contains(cell_polygon):
                        intersected_polygon = cell_polygon.intersection(gzd_geom)
                        if intersected_polygon:
                            intersected_centroid_lat, intersected_centroid_lon = (
                                intersected_polygon.centroid.y,
                                intersected_polygon.centroid.x,
                            )
                            interescted_mgrs_id = mgrs.toMgrs(
                                intersected_centroid_lat,
                                intersected_centroid_lon,
                                resolution,
                            )
                            mgrs_record = graticule_dggs_to_geoseries(
                                "mgrs",
                                interescted_mgrs_id,
                                resolution,
                                intersected_polygon,
                            )
                    mgrs_records.append(mgrs_record)
                pbar.update(1)
    return gpd.GeoDataFrame(mgrs_records, geometry="geometry", crs="EPSG:4326")


def mgrs_grid_ids(gzd, resolution, verbose=True):
    """
    Return a list of MGRS IDs for a given GZD and resolution.
    """
    if not is_valid_gzd(gzd):
        raise ValueError("Invalid GZD. Please input a valid GZD.")

    resolution = validate_mgrs_resolution(resolution)
    cell_size = 100_000 // (10**resolution)
    north_bands = "NPQRSTUVWX"
    south_bands = "MLKJHGFEDC"
    band_distance = 111_132 * 8
    gzd_band = gzd[2]

    if gzd_band >= "N":  # North Hemesphere
        epsg_code = int("326" + gzd[:2])
        min_x, min_y, max_x, max_y = 100000, 0, 900000, 9500000
        north_band_idx = north_bands.index(gzd_band)
        max_y = band_distance * (north_band_idx + 1)
        if gzd_band == "X":
            max_y += band_distance  # band X = 12 degrees instead of 8 degrees
    else:  # South Hemesphere
        epsg_code = int("327" + gzd[:2])
        min_x, min_y, max_x, max_y = 100000, 0, 900000, 10000000
        south_band_idx = south_bands.index(gzd_band)
        max_y = band_distance * (south_band_idx + 1)

    utm_crs = CRS.from_epsg(epsg_code)
    wgs84_crs = CRS.from_epsg(4326)
    transformer = Transformer.from_crs(utm_crs, wgs84_crs, always_xy=True).transform

    gzd_json_path = os.path.join(os.path.dirname(__file__), "gzd.geojson")
    with open(gzd_json_path, encoding="utf-8") as f:
        gzd_data = json.load(f)

    gzd_features = gzd_data["features"]
    gzd_feature = [
        feature for feature in gzd_features if feature["properties"].get("gzd") == gzd
    ][0]
    gzd_geom = shape(gzd_feature["geometry"])

    ids = []
    x_coords = np.arange(min_x, max_x, cell_size)
    y_coords = np.arange(min_y, max_y, cell_size)
    num_cells = len(x_coords) * len(y_coords)
    with tqdm(total=num_cells, desc="Generating MGRS IDs", unit=" cells", disable=not verbose) as pbar:
        for x in x_coords:
            for y in y_coords:
                cell_polygon_utm = Polygon(
                    [
                        (x, y),
                        (x + cell_size, y),
                        (x + cell_size, y + cell_size),
                        (x, y + cell_size),
                        (x, y),
                    ]
                )
                cell_polygon = transform(transformer, cell_polygon_utm)

                if cell_polygon.intersects(gzd_geom):
                    centroid_lat, centroid_lon = (
                        cell_polygon.centroid.y,
                        cell_polygon.centroid.x,
                    )
                    final_id = mgrs.toMgrs(centroid_lat, centroid_lon, resolution)
                    if not gzd_geom.contains(cell_polygon):
                        intersected_polygon = cell_polygon.intersection(gzd_geom)
                        if intersected_polygon:
                            intersected_centroid_lat, intersected_centroid_lon = (
                                intersected_polygon.centroid.y,
                                intersected_polygon.centroid.x,
                            )
                            final_id = mgrs.toMgrs(
                                intersected_centroid_lat,
                                intersected_centroid_lon,
                                resolution,
                            )
                    ids.append(final_id)
                pbar.update(1)

    return ids


def mgrsgrid(gzd, resolution, output_format="gpd", verbose=True):
    """
    Generate MGRS grid for pure Python usage.

    Args:
        gzd (str): Grid Zone Designator, e.g. '48P'.
        resolution (int): MGRS resolution [0..5].
        output_format (str, optional): Output format ('geojson', 'csv', 'geo', 'gpd', 'shapefile', 'gpkg', 'parquet', or None for list of MGRS IDs). Defaults to None.

    Returns:
        Depends on output_format: list, GeoDataFrame, file path, or GeoJSON FeatureCollection.
    """
    if not is_valid_gzd(gzd):
        raise ValueError("Invalid GZD. Please input a valid GZD.")
    gdf = mgrs_grid(gzd, resolution, verbose=verbose)

    output_name = f"mgrs_grid_{resolution}"
    return convert_to_output_format(gdf, output_format, output_name)


def mgrsgrid_cli():
    parser = argparse.ArgumentParser(description="Generate MGRS DGGS.")
    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        default=0,
        required=True,
        help="Resolution [0..5]",
    )
    parser.add_argument(
        "-gzd",
        type=str,
        default="48P",
        required=True,
        help="GZD - Grid Zone Designator, e.g. 48P",
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

    gzd = args.gzd
    if not is_valid_gzd(gzd):
        print("Invalid GZD. Please input a valid GZD and try again.")
        return
    resolution = args.resolution

    try:
        result = mgrsgrid(gzd, resolution, args.output_format, verbose=args.verbose)
        if args.output_format in STRUCTURED_FORMATS:
            print(result)
    except ValueError as e:
        print(f"Error: {str(e)}")
        return


if __name__ == "__main__":
    mgrsgrid_cli()
