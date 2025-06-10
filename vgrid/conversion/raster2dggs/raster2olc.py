import os, argparse, json
from tqdm import tqdm
import rasterio
from vgrid.utils import olc, mercantile
import numpy as np
from shapely.geometry import Polygon, Point, mapping
import json
from vgrid.stats.olcstats import olc_metrics
from vgrid.generator.settings import graticule_dggs_to_feature
from math import cos, radians
import re
from vgrid.conversion.latlon2dggs import latlon2olc
import csv

def get_nearest_olc_resolution(raster_path):
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
       
    # Find the nearest s2 resolution by comparing the pixel size to the s2 edge lengths
    nearest_resolution = None
    min_diff = float('inf')
        
    # Check resolutions from 0 to 29
    for res in range(10,13):
        _, _, avg_area = olc_metrics(res)
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

def raster2olc(raster_path, resolution=None, format='geojson'):
    # Step 1: Determine the nearest olc resolution if none is provided
    if resolution is None:
        resolution = get_nearest_olc_resolution(raster_path)
        print(f"Nearest olc resolution determined: {resolution}")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    olc_ids = set()
    
    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            olc_id = latlon2olc(lat,lon,resolution)            
            olc_ids.add(olc_id)

    # Sample the raster values at the centroids of the olc hexagons
    olc_data = []
    
    for olc_id in tqdm(olc_ids, desc="Resampling", unit=" cells"): 
        # Get the centroid of the olc cell
        coord = olc.decode(olc_id)
        centroid_lat, centroid_lon = coord.latitudeCenter, coord.longitudeCenter
        
        # Sample the raster values at the centroid (lat, lon)
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            olc_data.append({
                "olc": olc_id,
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    if format.lower() == 'csv':
        import io
        output = io.StringIO()
        if olc_data:
            writer = csv.DictWriter(output, fieldnames=olc_data[0].keys())
            writer.writeheader()
            writer.writerows(olc_data)
        return output.getvalue()
    
    
    # Create the GeoJSON-like structure
    olc_features = []
    for data in tqdm(olc_data, desc="Converting to GeoJSON", unit=" cells"):
        olc_id = data["olc"]
        coord = olc.decode(olc_id)    
        if coord:
            # Create the bounding box coordinates for the polygon
            min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
            max_lat, max_lon = coord.latitudeHi, coord.longitudeHi        
            cell_resolution = coord.codeLength 

            # Define the polygon based on the bounding box
            cell_polygon = Polygon([
                [min_lon, min_lat],  # Bottom-left corner
                [max_lon, min_lat],  # Bottom-right corner
                [max_lon, max_lat],  # Top-right corner
                [min_lon, max_lat],  # Top-left corner
                [min_lon, min_lat]   # Closing the polygon (same as the first point)
            ])
            
            olc_feature = graticule_dggs_to_feature("olc",olc_id,cell_resolution,cell_polygon)   
            band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
            olc_feature["properties"].update(convert_numpy_types(band_properties))
            olc_features.append(olc_feature)               
            
    return {
        "type": "FeatureCollection",
        "features": olc_features,
    }

def raster2olc_cli():
    """Command line interface for raster to OLC conversion"""
    parser = argparse.ArgumentParser(description="Convert Raster in Geographic CRS to OLC/ Google Plus Code DGGS")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default=None, help="Resolution [10..12]"
    )

    parser.add_argument(
        '-f', '--format', type=str, required=False, default='geojson',
        choices=['geojson', 'csv'], help="Output format (geojson or csv)"
    )

    args = parser.parse_args()
    raster = args.raster
    resolution = args.resolution
    format = args.format
    
    if not os.path.exists(raster):
        print(f"Error: The file {raster} does not exist.")
        return
        
    if resolution is not None:
        if resolution < 10 or resolution > 12:
            print(f"Please select a resolution in [10..12] range and try again ")
            return

    output_data = raster2olc(raster, resolution, format)
    output_name = os.path.splitext(os.path.basename(raster))[0]
    output_path = f"{output_name}2olc.{format}"
   
    if format.lower() == 'csv':
        with open(output_path, 'w', newline='') as f:
            if output_data:
                writer = csv.DictWriter(f, fieldnames=output_data[0].keys())
                writer.writeheader()
                writer.writerows(output_data)
    else:  # geojson
        with open(output_path, 'w') as f:
            json.dump(output_data, f)
            
    print(f"Output saved as {output_path}")

if __name__ == "__main__":
    raster2olc_cli()
