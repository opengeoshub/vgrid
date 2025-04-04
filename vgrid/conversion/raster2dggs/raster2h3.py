import os, argparse, json
from tqdm import tqdm
import rasterio
import h3
import numpy as np
import geopandas as gpd
from shapely.geometry import Polygon, Point, mapping
from collections import defaultdict
import json
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.rhealpixgrid import fix_rhealpix_antimeridian_cells


def get_nearest_h3_resolution(raster_path):
    # Open the raster to get its resolution
    with rasterio.open(raster_path) as src:
        # Get the pixel size (resolution) in the raster (in meters)
        pixel_size_x, pixel_size_y = src.res[0], src.res[1]
        
    # Calculate the diagonal pixel size (Pythagorean theorem)
    pixel_size = (pixel_size_x**2 + pixel_size_y**2)**0.5
    print (pixel_size)
    # Find the nearest H3 resolution by comparing the pixel size to the H3 edge lengths
    nearest_resolution = None
    min_diff = float('inf')
    
    # Check resolutions from 0 to 15 (since those are the typical H3 resolutions)
    for res in range(16):
        # Get the average edge length for the H3 resolution in meters
        cell_width = h3.average_hexagon_edge_length(res, unit='m')
        
        # Calculate the difference between the raster's pixel size and the H3 edge length
        diff = abs(cell_width - pixel_size)
        
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

def resample_raster_to_h3(raster_path, resolution=None):
    # Step 1: Determine the nearest H3 resolution if none is provided
    if resolution is None:
        resolution = get_nearest_h3_resolution(raster_path)
        print(f"Nearest H3 resolution determined: {resolution}")

    # Open the raster file to get metadata and data
    with rasterio.open(raster_path) as src:
        raster_data = src.read()  # Read all bands
        transform = src.transform
        crs = src.crs
        bounds = src.bounds  # (min_x, min_y, max_x, max_y)
        width, height = src.width, src.height
        band_count = src.count  # Number of bands in the raster

    # Generate H3 hexagons that intersect with the raster's bounding box
    min_lon, min_lat = transform * (0, 0)  # Top-left corner (pixel (0, 0))
    max_lon, max_lat = transform * (width, height)  # Bottom-right corner (pixel (width, height))
    
    # Create a grid of points covering the raster's bounds
    lon_step = (max_lon - min_lon) / width
    lat_step = (max_lat - min_lat) / height
    
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
                "centroid": Point(centroid_lon, centroid_lat),
                **{f"band_{i+1}": values[i] for i in range(band_count)}  # Create separate columns for each band
            })
    
    # Create the GeoJSON-like structure
    h3_features = []
    for data in h3_data:
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
            band_properties = {f"band_{i+1}": data[f"band_{i+1}"] for i in range(band_count)}
            h3_feature["properties"].update(convert_numpy_types(band_properties) )
            h3_features.append(h3_feature)               
          
    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }


       
# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Convert Raster to H3 Grid")
    parser.add_argument(
        '-raster', type=str, required=True, help="Raster file path (Point, Polyline or Polygon)"
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


    h3_geojson = resample_raster_to_h3(raster, resolution)
    print(h3_geojson)
    geojson_name = os.path.splitext(os.path.basename(raster))[0]
    geojson_path = f"{geojson_name}2h3.geojson"
   
    with open(geojson_path, 'w') as f:
        json.dump(h3_geojson, f)
    
    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
