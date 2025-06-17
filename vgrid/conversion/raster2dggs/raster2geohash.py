import os, argparse, json
from tqdm import tqdm
import rasterio
import numpy as np
from shapely.geometry import Polygon
import json
from vgrid.stats.geohashstats import geohash_metrics
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.conversion.latlon2dggs import latlon2geohash
from math import cos, radians
import re
from vgrid.utils import geohash
import csv

def get_nearest_geohash_resolution(raster_path):
    with rasterio.open(raster_path) as src:
        transform = src.transform
        crs = src.crs
        pixel_width = transform.a
        pixel_height = -transform.e
        cell_size = pixel_width*pixel_height
        
        if crs.is_geographic: 
            # Latitude of the raster center
            center_latitude = (src.bounds.top + src.bounds.bottom) / 2
            # Convert degrees to meters
            meter_per_degree_lat = 111_320  # Roughly 1 degree latitude in meters
            meter_per_degree_lon = meter_per_degree_lat * cos(radians(center_latitude))

            pixel_width_m = pixel_width * meter_per_degree_lon
            pixel_height_m = pixel_height * meter_per_degree_lat
            cell_size = pixel_width_m*pixel_height_m    
       
    nearest_resolution = None
    min_diff = float('inf')
        
    # Check resolutions from 1 to 10
    for res in range(1,11):
        _, _, avg_area = geohash_metrics(res)
        diff = abs(avg_area - cell_size)        
        # If the difference is smaller than the current minimum, update the nearest resolution
        if diff < min_diff:
            min_diff = diff
            nearest_resolution = res
    
    return nearest_resolution

def convert_numpy_types(obj):
    """ Recursively convert NumPy types to native Python types """
    if isinstance(obj, np.generic):
        return obj.item()  # Convert numpy types like np.uint8 to native Python int
    elif isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(v) for v in obj]
    else:
        return obj

def raster2geohash(raster_path, resolution=None, format='geojson'):
    # Step 1: Determine the nearest geohash resolution if none is provided
    if resolution is None:
        resolution = get_nearest_geohash_resolution(raster_path)
        print(f"Nearest geohash resolution determined: {resolution}")
    
    # Validate resolution is in valid range
    if resolution < 1 or resolution > 10:
        raise ValueError("Resolution must be in range [1..10]")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    geohash_ids = set()
    
    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            geohash_id = latlon2geohash(lat,lon,resolution)            
            geohash_ids.add(geohash_id)

    # Sample the raster values at the centroids of the geohash hexagons
    geohash_data = []
    
    for geohash_id in tqdm(geohash_ids, desc="Resampling", unit=" cells"):
        # Get the centroid of the geohash cell
        centroid_lat, centroid_lon = geohash.decode(geohash_id)        
        # Sample the raster values at the centroid (lat, lon)
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            geohash_data.append({
                "geohash": geohash_id,
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    if format.lower() == 'csv':
        import io
        output = io.StringIO()
        if geohash_data:
            writer = csv.DictWriter(output, fieldnames=geohash_data[0].keys())
            writer.writeheader()
            writer.writerows(geohash_data)
        return output.getvalue()
    
    # Create the GeoJSON-like structure
    geohash_features = []
    for data in tqdm(geohash_data, desc="Converting to GeoJSON", unit=" cells"):
        geohash_id = data["geohash"]
        geohash_bbox =  geohash.bbox(geohash_id)
        if geohash_bbox:
            min_lat, min_lon = geohash_bbox['s'], geohash_bbox['w']  # Southwest corner
            max_lat, max_lon = geohash_bbox['n'], geohash_bbox['e']  # Northeast corner        
            cell_resolution =  len(geohash_id)
            # Define the polygon based on the bounding box
            cell_polygon = Polygon([
                [min_lon, min_lat],  # Bottom-left corner
                [max_lon, min_lat],  # Bottom-right corner
                [max_lon, max_lat],  # Top-right corner
                [min_lon, max_lat],  # Top-left corner
                [min_lon, min_lat]   # Closing the polygon (same as the first point)
            ])
            geohash_feature = graticule_dggs_to_feature("geohash",geohash_id,cell_resolution,cell_polygon)   
            band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
            geohash_feature["properties"].update(convert_numpy_types(band_properties) )
            geohash_features.append(geohash_feature)               
            
    return {
        "type": "FeatureCollection",
        "features": geohash_features,
    }

def raster2geohash_cli():
    parser = argparse.ArgumentParser(description="Convert Raster in Geographic CRS to Geohash DGGS")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default=None, 
        help="Resolution [1..10]", min=1, max=10
    )

    parser.add_argument(
        '-f', '--format', type=str, required=False, default='geojson', 
        choices=['geojson', 'csv'], help="Output format (geojson or csv)"
    )

    args = parser.parse_args()
    raster = args.raster
    resolution = args.resolution
    output_format = args.format
    
    if not os.path.exists(raster):
        print(f"Error: The file {raster} does not exist.")
        return

    result = raster2geohash(raster, resolution, output_format)
    base_name = os.path.splitext(os.path.basename(raster))[0]
    
    if output_format.lower() == 'csv':
        output_path = f"{base_name}2geohash.csv"
        # Get all possible field names from the first row
        fieldnames = list(result[0].keys())
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(result)
    else:
        output_path = f"{base_name}2geohash.geojson"
        with open(output_path, 'w') as f:
            json.dump(result, f)
    
    print(f"Output saved as {output_path}")