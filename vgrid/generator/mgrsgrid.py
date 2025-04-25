import numpy as np
from shapely.geometry import shape, Polygon
from shapely.ops import transform
from pyproj import CRS, Transformer
import argparse
import re
from tqdm import tqdm 
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.utils import mgrs
import json
import os

def is_valid_gzd(gzd):
    """Check if a Grid Zone Designator (GZD) is valid."""
    pattern = r"^(?:0[1-9]|[1-5][0-9]|60)[C-HJ-NP-X]$"
    return bool(re.match(pattern, gzd))

def generate_grid(gzd, resolution):    # Define the UTM CRS  
    # Reference: https://www.maptools.com/tutorials/utm/details
    cell_size = 100_000 // (10 ** resolution)   
    north_bands = 'NPQRSTUVWX'
    south_bands = 'MLKJHGFEDC'
    band_distance = 111_132 * 8 
    gzd_band = gzd[2]

    if gzd_band >='N': # North Hemesphere
        epsg_code = int('326' + gzd[:2])      
        min_x, min_y, max_x, max_y = 100000, 0, 900000, 9500000 # for the North   
        north_band_idx = north_bands.index(gzd_band)
        max_y = band_distance * (north_band_idx + 1)
        if gzd_band == 'X':
            max_y += band_distance # band X = 12 deggrees instead of 8 degrees
    
    else:  # South Hemesphere
        epsg_code = int('327' + gzd[:2])     
        min_x, min_y, max_x, max_y = 100000, 0, 900000, 10000000 # for the South         
        south_band_idx = south_bands.index(gzd_band)
        max_y = band_distance * (south_band_idx + 1 )
   
    utm_crs = CRS.from_epsg(epsg_code)
    wgs84_crs = CRS.from_epsg(4326)
    transformer = Transformer.from_crs(utm_crs, wgs84_crs, always_xy=True).transform

    gzd_json_path = os.path.join(os.path.dirname(__file__), 'gzd.geojson')            
    with open(gzd_json_path, 'r') as f:
        gzd_data = json.load(f)
    
    gzd_features = gzd_data["features"]
    gzd_feature = [feature for feature in gzd_features if feature["properties"].get("gzd") == gzd][0]
    gzd_geom = shape(gzd_feature["geometry"])

    # Create grid polygons
    mgrs_features = []
    x_coords = np.arange(min_x, max_x, cell_size)
    y_coords = np.arange(min_y, max_y, cell_size)
    num_cells = len(x_coords) * len(y_coords)
    with tqdm(total=num_cells, desc="Generating MGRS DGGS", unit=" cells") as pbar:
        for x in x_coords:
            for y in y_coords:
                cell_polygon_utm = Polygon([
                    (x, y),
                    (x + cell_size, y),
                    (x + cell_size, y + cell_size),
                    (x, y + cell_size),
                    (x, y)  # Close the polygon
                ])            
                cell_polygon = transform(transformer, cell_polygon_utm)
                                
                if cell_polygon.intersects(gzd_geom):
                    centroid_lat, centroid_lon  =  cell_polygon.centroid.y, cell_polygon.centroid.x,
                    mgrs_id = mgrs.toMgrs(centroid_lat, centroid_lon, resolution)
                    mgrs_feature = graticule_dggs_to_feature("mgrs",mgrs_id,resolution,cell_polygon)             
                    # clip inside GZD:
                    if not gzd_geom.contains(cell_polygon):
                        intersected_polygon = cell_polygon.intersection(gzd_geom)  
                        if intersected_polygon:
                            intersected_centroid_lat, intersected_centroid_lon  =  intersected_polygon.centroid.y, intersected_polygon.centroid.x,
                            interescted_mgrs_id = mgrs.toMgrs(intersected_centroid_lat, intersected_centroid_lon, resolution)            
                            mgrs_feature = graticule_dggs_to_feature("mgrs",interescted_mgrs_id,resolution,intersected_polygon)
                    mgrs_features.append(mgrs_feature)
                pbar.update(1)
    return {
        "type": "FeatureCollection",
        "features": mgrs_features
    }
    

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate MGRS DGGS.")
    parser.add_argument("-r", "--resolution", type=int, default=0, required=True, help="Resolution [0..5]")
    parser.add_argument("-gzd", type = str, default='48P', required=True, help="GZD - Grid Zone Designator, e.g. 48P")
    # Parse the arguments
    args = parser.parse_args()
   
    gzd = args.gzd   
    if not is_valid_gzd(gzd):
        print(f"Invalid GZD. Please input a valid GZD and try again.")
        return
    
    resolution = args.resolution
    if resolution < 0 or resolution > 5:
        print(f"Please select a resolution in [0..5] range and try again ")
        return

    geojson_features = generate_grid(gzd, resolution)
    
    if geojson_features:
        geojson_path = f"mgrs_grid_{gzd}_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
    
if __name__ == "__main__":
    main()