import os, argparse, json
from tqdm import tqdm
import rasterio
import h3
import numpy as np
from shapely.geometry import Polygon, Point, mapping
import json
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

def raster_to_h3(raster_path, resolution=None):
    # Step 1: Determine the nearest H3 resolution if none is provided
    if resolution is None:
        resolution = get_nearest_h3_resolution(raster_path)
        print(f"Nearest H3 resolution determined: {resolution}")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    h3_cells = set()
    
    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            h3_index = h3.latlng_to_cell(lat, lon,resolution)
            h3_cells.add(h3_index)

    # Sample the raster values at the centroids of the H3 hexagons
    h3_data = []
    
    for h3_index in h3_cells:
        # Get the centroid of the H3 cell
        centroid_lat, centroid_lon = h3.cell_to_latlng(h3_index)
        
        # Sample the raster values at the centroid (lat, lon)
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            h3_data.append({
                "h3": h3_index,
                # "centroid": Point(centroid_lon, centroid_lat),
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    # Create the GeoJSON-like structure
    h3_features = []
    for data in tqdm(h3_data, desc="Resampling", unit=" cells"):
        cell_boundary = h3.cell_to_boundary(data["h3"])   
        if cell_boundary:
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            # Reverse lat/lon to lon/lat for GeoJSON compatibility
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)
            resolution = h3.get_resolution(data["h3"]) 
            h3_feature = {
                "type": "Feature",
                "geometry":mapping(cell_polygon),
                "properties": {
                    "h3": data["h3"],
                    "resolution": resolution                    
                }
            }
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

       
# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Convert Raster in Geographic CRS to H3 Grid")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default= None, help="Resolution of H3 to be generated"
    )


    args = parser.parse_args()
    raster = args.raster
    resolution = args.resolution
    
    if not os.path.exists(raster):
        print(f"Error: The file {raster} does not exist.")
        return
    if resolution is not None:
        if resolution < 0 or resolution > 15:
            print(f"Please select a resolution in [0..15] range and try again ")
            return


    h3_geojson = raster_to_h3(raster, resolution)
    geojson_name = os.path.splitext(os.path.basename(raster))[0]
    geojson_path = f"{geojson_name}2h3.geojson"
   
    with open(geojson_path, 'w') as f:
        json.dump(h3_geojson, f)
    
    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
