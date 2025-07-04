import os
import argparse
import json
from tqdm import tqdm
import rasterio
from vgrid.utils import tilecode, mercantile
import numpy as np
from shapely.geometry import Polygon
from vgrid.stats.quadkeystats import quadkey_metrics
from vgrid.generator.settings import graticule_dggs_to_feature
from math import cos, radians
import csv


def get_nearest_quadkey_resolution(raster_path):
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        pixel_width = transform.a
        pixel_height = -transform.e
        cell_size = pixel_width * pixel_height

        if crs.is_geographic:
            # Latitude of the raster center
            center_latitude = (src.bounds.top + src.bounds.bottom) / 2
            # Convert degrees to meters
            meter_per_degree_lat = 111_320  # Roughly 1 degree latitude in meters
            meter_per_degree_lon = meter_per_degree_lat * cos(radians(center_latitude))

            pixel_width_m = pixel_width * meter_per_degree_lon
            pixel_height_m = pixel_height * meter_per_degree_lat
            cell_size = pixel_width_m * pixel_height_m

    nearest_resolution = None
    min_diff = float("inf")

    # Check resolutions from 0 to 29
    for res in range(30):
        _, _, avg_area = quadkey_metrics(res)
        diff = abs(avg_area - cell_size)
        # If the difference is smaller than the current minimum, update the nearest resolution
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res

    return nearest_resolution


def convert_numpy_types(obj):
    """Recursively convert NumPy types to native Python types"""
    if isinstance(obj, np.generic):
        return obj.item()  # Convert numpy types like np.uint8 to native Python int
    elif isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(v) for v in obj]
    else:
        return obj


def raster2quadkey(raster_path, resolution=None, format="geojson"):
    """Convert raster to quadkey format

    Args:
        raster_path (str): Path to input raster file
        resolution (int, optional): Quadkey resolution level [0-29]. If None, will be determined automatically
        format (str, optional): Output format, either 'geojson' or 'csv'. Defaults to 'geojson'

    Returns:
        dict: For geojson format, returns a GeoJSON FeatureCollection
        str: For csv format, returns a CSV string
    """
    # Step 1: Determine the nearest quadkey resolution if none is provided
    if resolution is None:
        resolution = get_nearest_quadkey_resolution(raster_path)
        print(f"Nearest quadkey resolution determined: {resolution}")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    quadkey_ids = set()

    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            quadkey_id = tilecode.latlon2quadkey(lat, lon, resolution)
            quadkey_ids.add(quadkey_id)

    # Sample the raster values at the centroids of the quadkey hexagons
    quadkey_data = []

    for quadkey_id in tqdm(quadkey_ids, desc="Resampling", unit=" cells"):
        # Get the centroid of the quadkey cell
        centroid_lat, centroid_lon = tilecode.quadkey2latlon(quadkey_id)

        # Sample the raster values at the centroid (lat, lon)
        col, row = ~transform * (centroid_lon, centroid_lat)

        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            quadkey_data.append(
                {
                    "quadkey": quadkey_id,
                    **{
                        f"band_{i + 1}": values[i] for i in range(band_count)
                    },  # Create separate columns for each band
                }
            )

    if format == "csv":
        import io

        output = io.StringIO()
        if quadkey_data:
            writer = csv.DictWriter(output, fieldnames=quadkey_data[0].keys())
            writer.writeheader()
            writer.writerows(quadkey_data)
        return output.getvalue()

    # Create the GeoJSON-like structure
    quadkey_features = []
    for data in tqdm(quadkey_data, desc="Converting to GeoJSON", unit=" cells"):
        quadkey_id = data["quadkey"]
        tile = mercantile.quadkey_to_tile(quadkey_id)
        # Format as tilecode_id
        z = tile.z
        x = tile.x
        y = tile.y
        # Get the bounds of the tile in (west, south, east, north)
        bounds = mercantile.bounds(x, y, z)
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

        cell_resolution = z
        quadkey_feature = graticule_dggs_to_feature(
            "quadkey", quadkey_id, cell_resolution, cell_polygon
        )
        band_properties = {
            f"band_{i + 1}": data[f"band_{i + 1}"] for i in range(band_count)
        }
        quadkey_feature["properties"].update(convert_numpy_types(band_properties))
        quadkey_features.append(quadkey_feature)

    return {
        "type": "FeatureCollection",
        "features": quadkey_features,
    }


def raster2quadkey_cli():
    parser = argparse.ArgumentParser(
        description="Convert Raster in Geographic CRS to Quadkey DGGS"
    )
    parser.add_argument("-raster", type=str, required=True, help="Raster file path")

    parser.add_argument(
        "-r",
        "--resolution",
        type=int,
        required=False,
        default=None,
        help="Resolution [0..29]",
    )

    parser.add_argument(
        "-f",
        "--format",
        type=str,
        required=False,
        default="geojson",
        choices=["geojson", "csv"],
        help="Output format: 'geojson' for GeoJSON format with geometry, 'csv' for tabular data with quadkey and band values",
    )

    args = parser.parse_args()
    raster = args.raster
    resolution = args.resolution
    format = args.format.lower()

    if not os.path.exists(raster):
        print(f"Error: The file {raster} does not exist.")
        return

    if resolution is not None:
        if resolution < 0 or resolution > 29:
            print("Please select a resolution in [0..29] range and try again ")
            return

    result = raster2quadkey(raster, resolution, format)

    output_name = os.path.splitext(os.path.basename(raster))[0]
    output_path = f"{output_name}2quadkey.{format}"

    if format.lower() == "csv":
        with open(output_path, "w", newline="") as f:
            f.write(result)
    else:  # geojson
        with open(output_path, "w") as f:
            json.dump(result, f)

    print(f"Output saved as {output_path}")
