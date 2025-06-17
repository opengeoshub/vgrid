import os, argparse, json
from tqdm import tqdm
import rasterio
import h3
import numpy as np
from shapely.geometry import Polygon, Point, mapping
import json
import csv
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import geodesic_dggs_to_feature
from math import cos, radians

def get_nearest_h3_resolution(raster_path):
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
        
    # Check resolutions from 0 to 15
    for res in range(16):
        avg_area = h3.average_hexagon_area(res, unit='m^2')
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

def raster2h3(raster_path, resolution=None, format='geojson'):
    # Step 1: Determine the nearest H3 resolution if none is provided
    if resolution is None:
        resolution = get_nearest_h3_resolution(raster_path)
        print(f"Nearest H3 resolution determined: {resolution}")
    
    # Validate resolution is in valid range
    if resolution < 0 or resolution > 15:
        raise ValueError("Resolution must be in range [0..15]")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    h3_ids = set()
    
    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            h3_id = h3.latlng_to_cell(lat, lon,resolution)
            h3_ids.add(h3_id)

    # Sample the raster values at the centroids of the H3 hexagons
    h3_data = []
    
    for h3_id in tqdm(h3_ids, desc="Resampling", unit=" cells"):
        # Get the centroid of the H3 cell
        centroid_lat, centroid_lon = h3.cell_to_latlng(h3_id)
        # Sample the raster values at the centroid (lat, lon)
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            h3_data.append({
                "h3": h3_id,
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    if format.lower() == 'csv':
        import io
        output = io.StringIO()
        if h3_data:
            writer = csv.DictWriter(output, fieldnames=h3_data[0].keys())
            writer.writeheader()
            writer.writerows(h3_data)
        return output.getvalue()
    
    # Create the GeoJSON-like structure
    h3_features = []
    for data in tqdm(h3_data, desc="Converting to GeoJSON", unit=" cells"):
        cell_boundary = h3.cell_to_boundary(data["h3"])   
        if cell_boundary:
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            # Reverse lat/lon to lon/lat for GeoJSON compatibility
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)            
            num_edges = 6
            if (h3.is_pentagon(data["h3"])):
                num_edges = 5
            
            h3_feature = geodesic_dggs_to_feature("h3",data["h3"],resolution,cell_polygon,num_edges)   
        
            band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
            h3_feature["properties"].update(convert_numpy_types(band_properties) )
            h3_features.append(h3_feature)               
          
    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }

def raster2h3_cli():
    parser = argparse.ArgumentParser(description="Convert Raster in Geographic CRS to H3 DGGS")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default=None, 
        help="Resolution [0..15]", min=0, max=15
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
        if resolution < 0 or resolution > 15:
            print(f"Please select a resolution in [0..15] range and try again ")
            return

    result = raster2h3(raster, resolution, format)
    base_name = os.path.splitext(os.path.basename(raster))[0]
    output_path = f"{base_name}2h3.{format}"
    
    
    if format.lower() == 'geojson':
        with open(output_path, 'w') as f:
            json.dump(result, f)
    
    elif format.lower() == 'csv':
        # Get all possible field names from the first row
        fieldnames = list(result[0].keys())
        with open(output_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(result)
    
    print(f"Output saved as {output_path}")