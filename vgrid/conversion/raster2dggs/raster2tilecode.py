import os, argparse, json
from tqdm import tqdm
import rasterio
from vgrid.utils import tilecode
import numpy as np
from shapely.geometry import Polygon, Point, mapping
import json
from vgrid.generator.tilecodegrid import fix_tilecode_antimeridian_cells
from vgrid.generator.settings import geodesic_dggs_to_feature
from math import cos, radians

def get_nearest_tilecode_resolution(raster_path):
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
        
    # Check resolutions from 0 to 15
    for res in range(16):
        avg_area = tilecode.average_hexagon_area(res, unit='m^2')
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

def resample_raster_to_tilecode(raster_path, resolution=None):
    # Step 1: Determine the nearest tilecode resolution if none is provided
    if resolution is None:
        resolution = get_nearest_tilecode_resolution(raster_path)
        print(f"Nearest tilecode resolution determined: {resolution}")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    tilecode_cells = set()
    
    for row in range(height):
        for col in range(width):
            lon, lat = transform * (col, row)
            tilecode_index = tilecode.latlng_to_cell(lat, lon,resolution)
            tilecode_cells.add(tilecode_index)

    # Sample the raster values at the centroids of the tilecode hexagons
    tilecode_data = []
    
    for tilecode_index in tilecode_cells:
        # Get the centroid of the tilecode cell
        centroid_lat, centroid_lon = tilecode.cell_to_latlng(tilecode_index)
        
        # Sample the raster values at the centroid (lat, lon)
        col, row = ~transform * (centroid_lon, centroid_lat)
        
        if 0 <= col < width and 0 <= row < height:
            # Get the values for all bands at this centroid
            values = raster_data[:, int(row), int(col)]
            tilecode_data.append({
                "tilecode": tilecode_index,
                # "centroid": Point(centroid_lon, centroid_lat),
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    # Create the GeoJSON-like structure
    tilecode_features = []
    for data in tqdm(tilecode_data, desc="Resampling", unit=" cells"):
        cell_boundary = tilecode.cell_to_boundary(data["tilecode"])   
        if cell_boundary:
            filtered_boundary = fix_tilecode_antimeridian_cells(cell_boundary)
            # Reverse lat/lon to lon/lat for GeoJSON compatibility
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)
            resolution = tilecode.get_resolution(data["tilecode"]) 
            tilecode_feature = {
                "type": "Feature",
                "geometry":mapping(cell_polygon),
                "properties": {
                    "tilecode": data["tilecode"],
                    "resolution": resolution                    
                }
            }
            num_edges = 6
            if (tilecode.is_pentagon(data["tilecode"])):
                num_edges = 5
            
            tilecode_feature = geodesic_dggs_to_feature("tilecode",data["tilecode"],resolution,cell_polygon,num_edges)   
        
            band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
            tilecode_feature["properties"].update(convert_numpy_types(band_properties) )
            tilecode_features.append(tilecode_feature)               
          
    return {
        "type": "FeatureCollection",
        "features": tilecode_features,
    }

       
# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Convert Raster to tilecode Grid")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path"
    )
    
    parser.add_argument(
        '-r', '--resolution', type=int, required=False, default= None, help="Resolution of tilecode to be generated"
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


    tilecode_geojson = resample_raster_to_tilecode(raster, resolution)
    geojson_name = os.path.splitext(os.path.basename(raster))[0]
    geojson_path = f"{geojson_name}2tilecode.geojson"
   
    with open(geojson_path, 'w') as f:
        json.dump(tilecode_geojson, f)
    
    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
