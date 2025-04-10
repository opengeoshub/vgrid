import os, argparse, json
from tqdm import tqdm
import rasterio
from vgrid.utils.qtm import constructGeometry, qtm_id_to_facet

import numpy as np
from shapely.geometry import Polygon, Point, mapping
import json
from vgrid.stats.qtmstats import qtm_metrics
from vgrid.generator.settings import geodesic_dggs_to_feature
from math import cos, radians
import re
from vgrid.conversion.latlon2dggs import latlon2qtm

def get_nearest_qtm_resolution(raster_path):
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
    for res in range(1, 25):
        _, _, avg_area = qtm_metrics(res)
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

def raster_to_qtm(raster_path, resolution=None):
    # Step 1: Determine the nearest qtm resolution if none is provided
    if resolution is None:
        resolution = get_nearest_qtm_resolution(raster_path)
        print(f"Nearest qtm resolution determined: {resolution}")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    qtm_ids = set()
    
    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            qtm_id = latlon2qtm(lat,lon,resolution)            
            qtm_ids.add(qtm_id)

    # Sample the raster values at the centroids of the qtm hexagons
    qtm_data = []
    
    for qtm_id in qtm_ids:
        # Get the centroid of the qtm cell
        facet = qtm_id_to_facet(qtm_id)
        cell_polygon = constructGeometry(facet)    
        centroid = cell_polygon.centroid
        centroid_lon, centroid_lat = centroid.x, centroid.y
                
        # Sample the raster values at the centroid (lat, lon)
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            qtm_data.append({
                "qtm": qtm_id,
                # "centroid": Point(centroid_lon, centroid_lat),
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    # Create the GeoJSON-like structure
    qtm_features = []
    for data in tqdm(qtm_data, desc="Resampling", unit=" cells"):
        qtm_id = data["qtm"]
        facet = qtm_id_to_facet(qtm_id)
        cell_polygon = constructGeometry(facet)    
        cell_resolution = len(qtm_id)
        num_edges = 3
        qtm_feature = geodesic_dggs_to_feature("qtm",qtm_id,cell_resolution,cell_polygon,num_edges)   
        band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
        qtm_feature["properties"].update(convert_numpy_types(band_properties) )
        qtm_features.append(qtm_feature)               
            
    return {
        "type": "FeatureCollection",
        "features": qtm_features,
    }
 
       
# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Convert Raster to qtm Grid")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default= None, help="Resolution of qtm [1..24] to be generated"
    )


    args = parser.parse_args()
    raster = args.raster
    resolution = args.resolution
    
    if not os.path.exists(raster):
        print(f"Error: The file {raster} does not exist.")
        return
    if resolution is not None:
        if resolution < 1 or resolution > 24:
            print(f"Please select a resolution in [1..24] range and try again ")
            return


    qtm_geojson = raster_to_qtm(raster, resolution)
    geojson_name = os.path.splitext(os.path.basename(raster))[0]
    geojson_path = f"{geojson_name}2qtm.geojson"
   
    with open(geojson_path, 'w') as f:
        json.dump(qtm_geojson, f)
    
    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
