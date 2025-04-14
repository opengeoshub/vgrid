import os, argparse, json
from tqdm import tqdm
import rasterio
from vgrid.utils import tilecode, mercantile
import numpy as np
from shapely.geometry import Polygon
import json
from vgrid.stats.quadkeystats import quadkey_metrics
from vgrid.generator.settings import graticule_dggs_to_feature
from math import cos, radians
import re

def get_nearest_quadkey_resolution(raster_path):
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
    """ Recursively convert NumPy types to native Python types """
    if isinstance(obj, np.generic):
        return obj.item()  # Convert numpy types like np.uint8 to native Python int
    elif isinstance(obj, dict):
        return {k: convert_numpy_types(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [convert_numpy_types(v) for v in obj]
    else:
        return obj

def raster_to_quadkey(raster_path, resolution=None):
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
            quadkey_id = tilecode.latlon2quadkey(lat,lon,resolution)            
            quadkey_ids.add(quadkey_id)

    # Sample the raster values at the centroids of the quadkey hexagons
    quadkey_data = []
    
    for quadkey_id in quadkey_ids:
        # Get the centroid of the quadkey cell
        centroid_lat, centroid_lon = tilecode.quadkey2latlon(quadkey_id)
        
        # Sample the raster values at the centroid (lat, lon)
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            quadkey_data.append({
                "quadkey": quadkey_id,
                # "centroid": Point(centroid_lon, centroid_lat),
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    # Create the GeoJSON-like structure
    quadkey_features = []
    for data in tqdm(quadkey_data, desc="Resampling", unit=" cells"):
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
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        
        cell_resolution = z
        quadkey_feature = graticule_dggs_to_feature("quadkey",quadkey_id,cell_resolution,cell_polygon)   
        band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
        quadkey_feature["properties"].update(convert_numpy_types(band_properties) )
        quadkey_features.append(quadkey_feature)               
            
    return {
        "type": "FeatureCollection",
        "features": quadkey_features,
    }

       
# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Convert Raster in Geographic CRS to Quadkey Grid")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default= None, help="Resolution of Quadkey [0..29] to be generated"
    )


    args = parser.parse_args()
    raster = args.raster
    resolution = args.resolution
    
    if not os.path.exists(raster):
        print(f"Error: The file {raster} does not exist.")
        return
    if resolution is not None:
        if resolution < 0 or resolution > 29:
            print(f"Please select a resolution in [0..29] range and try again ")
            return


    quadkey_geojson = raster_to_quadkey(raster, resolution)
    geojson_name = os.path.splitext(os.path.basename(raster))[0]
    geojson_path = f"{geojson_name}2quadkey.geojson"
   
    with open(geojson_path, 'w') as f:
        json.dump(quadkey_geojson, f)
    
    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
