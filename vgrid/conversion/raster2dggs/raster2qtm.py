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
import csv


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

def raster2qtm(raster_path, resolution=None, format='geojson'):
    """
    Convert raster data to QTM DGGS format.
    
    Args:
        raster_path (str): Path to the raster file
        resolution (int, optional): QTM resolution level. If None, will be determined automatically.
        format (str): Output format, either 'geojson' or 'csv'
        
    Returns:
        dict or str: A dictionary containing the QTM data in GeoJSON format, or a CSV string if format is 'csv'
    """
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
    
    for qtm_id in tqdm(qtm_ids, desc="Resampling", unit=" cells"):
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
                **{f"band_{i+1}": values[i] for i in range(band_count)}
            })
    
    if format.lower() == 'csv':
        import io
        output = io.StringIO()
        if qtm_data:
            writer = csv.DictWriter(output, fieldnames=qtm_data[0].keys())
            writer.writeheader()
            writer.writerows(qtm_data)
        return output.getvalue()
    
    # Create the GeoJSON-like structure
    qtm_features = []
    for data in tqdm(qtm_data, desc="Converting to GeoJSON", unit=" cells"):
        qtm_id = data["qtm"]
        facet = qtm_id_to_facet(qtm_id)
        cell_polygon = constructGeometry(facet)    
        cell_resolution = len(qtm_id)
        num_edges = 3
        qtm_feature = geodesic_dggs_to_feature("qtm",qtm_id,cell_resolution,cell_polygon,num_edges)   
        band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
        qtm_feature["properties"].update(convert_numpy_types(band_properties))
        qtm_features.append(qtm_feature)               
            
    return {
        "type": "FeatureCollection",
        "features": qtm_features
    }

def raster2qtm_cli():
    """Command line interface for raster2qtm conversion"""
    parser = argparse.ArgumentParser(description="Convert Raster in Geographic CRS to QTM DGGS")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default=None, help="Resolution [1..24]"
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
        if resolution < 1 or resolution > 24:
            print(f"Please select a resolution in [1..24] range and try again ")
            return

    # Process the raster
    result = raster2qtm(raster, resolution, format)
    
    # Generate output filename
    base_name = os.path.splitext(os.path.basename(raster))[0]
    output_path = f"{base_name}2qtm.{format}"
    
    # Save the output
    if format.lower() == 'geojson':
        with open(output_path, 'w') as f:
            json.dump(result, f)
    else:  # csv
        with open(output_path, 'w', newline='') as f:
            f.write(result)
    
    print(f"Output saved as {output_path}")

if __name__ == "__main__":
    raster2qtm_cli()
