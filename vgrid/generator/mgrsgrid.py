import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon
from pyproj import CRS
import argparse
import re
from tqdm import tqdm 
import fiona
import os

def clip_grid_to_gzd(grid_gdf, gzd_geojson, input_gzd):
    # Load GZD polygons
    # gzd_gdf = gpd.read_file(gzd_geojson) # prone to fiona and geopandas version conflicts
    with fiona.open(gzd_geojson, "r") as src:
        gzd_gdf = gpd.GeoDataFrame.from_features(src, crs=src.crs)
        

    # Filter GZD based on input_gzd
    selected_gzd = gzd_gdf[gzd_gdf['gzd'] == input_gzd]

    if selected_gzd.empty:
        print(f"No matching GZD found for {input_gzd}.")
        return None

    # Clip the grid using the selected GZD polygon
    clipped_grid = gpd.clip(grid_gdf, selected_gzd)

    return clipped_grid


def generate_grid(gzd, cell_size):    # Define the UTM CRS  
     
    min_x, min_y, max_x, max_y = 100000, 0, 900000, 9500000 # for the North      

    if gzd[2] >='N': # North Hemesphere
        epsg_code = int('326' + gzd[:2])       
        utm_crs = CRS.from_epsg(epsg_code)
    else:  # South Hemesphere
        epsg_code = int('327' + gzd[:2])     
        utm_crs = CRS.from_epsg(epsg_code)        
        min_x, min_y, max_x, max_y = 100000, 100000, 900000, 10000000 # for the South 
        
    # Generate x and y coordinates
    x_coords = np.arange(min_x, max_x, cell_size)
    y_coords = np.arange(min_y, max_y, cell_size)

    # Create grid polygons
    polygons = []
    for x in x_coords:
        for y in tqdm(y_coords, desc="Processing cells ", unit="cells"):
            polygons.append(Polygon([
                (x, y),
                (x + cell_size, y),
                (x + cell_size, y + cell_size),
                (x, y + cell_size),
                (x, y)  # Close the polygon
            ]))

    # Create a GeoDataFrame
    grid_gdf = gpd.GeoDataFrame(geometry=polygons, crs=utm_crs)
    wg284_grid_gdf = grid_gdf.to_crs(epsg=4326)
    
    gzd_geojson =  os.path.join(os.path.dirname(__file__), 'gzd.geojson')
    clipped_grid = clip_grid_to_gzd(wg284_grid_gdf, gzd_geojson, gzd)

    return clipped_grid

def is_valid_gzd(gzd):
    """Check if a Grid Zone Designator (GZD) is valid."""
    pattern = r"^(?:0[1-9]|[1-5][0-9]|60)[C-HJ-NP-X]$"
    return bool(re.match(pattern, gzd))


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate MGRS grid.")
    parser.add_argument("-cellsize", type=int, default='100_000', required=True, help="Cell size in meters, e.g. 100_0000")
    parser.add_argument("-gzd", type = str, default='48P', required=True, help="GZD - Grid Zone Designator, e.g. 48P")
    # Parse the arguments
    args = parser.parse_args()
    gzd = args.gzd   
    if not is_valid_gzd(gzd):
        print(f"Invalid GZD. Please input a valid GZD and try again.")
        return
    
    cellsize = args.cellsize

    # Create the grid with the specified cell size
    mgrs_grid = generate_grid(gzd,cellsize)
  
    geojson_path = f"mgrs_grid_{gzd}_{cellsize}.geojson"
    # Save as GeoJSON
    mgrs_grid.to_file(geojson_path, driver="GeoJSON")

    print(f"Grid saved to {geojson_path}")


if __name__ == "__main__":
    main()