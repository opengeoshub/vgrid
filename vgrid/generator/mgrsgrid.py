import numpy as np
from shapely.geometry import shape, Polygon,mapping
from pyproj import CRS, Transformer
import argparse
import re
from tqdm import tqdm 
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.utils import mgrs
import json
import os

def mgrs_is_fully_within(mgrs_feature, gzd_feature):
    mgrs_geom = shape(mgrs_feature["geometry"])
    gzd_geom = shape(gzd_feature["geometry"])
    if gzd_geom.contains(mgrs_geom):  # Check if gzd_geom fully contains mgrs_geom
        return True  # At least one GZD feature fully contains the MGRS feature
    return False  # No GZD feature fully contains the MGRS feature

def mgrs_get_intersection(mgrs_feature, gzd_feature):
    mgrs_geom = shape(mgrs_feature["geometry"])    
    gzd_geom = shape(gzd_feature["geometry"])        
    intersection_geom = mgrs_geom.intersection(gzd_geom)  # Get intersection geometry
    if not intersection_geom.is_empty:
        return intersection_geom 
    return None

def utm_to_wgs84(polygon, transformer):
    """Transform a polygon from UTM to WGS84."""
    return Polygon([transformer.transform(x, y) for x, y in polygon.exterior.coords])

def generate_grid(gzd, resolution):    # Define the UTM CRS  
    cell_size = 100_000 // (10 ** resolution)
    
     # Reference: https://www.maptools.com/tutorials/utm/details
    min_x, min_y, max_x, max_y = 100000, 0, 900000, 9500000 # for the North   
    # min_x, min_y, max_x, max_y = 160000, 0, 834000, 9500000 # for the North    
        
    if gzd[2] >='N': # North Hemesphere
        epsg_code = int('326' + gzd[:2])      
        utm_crs = CRS.from_epsg(epsg_code)
    else:  # South Hemesphere
        epsg_code = int('327' + gzd[:2])     
        utm_crs = CRS.from_epsg(epsg_code)        
        min_x, min_y, max_x, max_y = 100000, 0, 900000, 10000000 # for the South 
        # min_x, min_y, max_x, max_y = 160000, 0, 834000, 10000000 # for the South 
   
    utm_crs = CRS.from_epsg(epsg_code)
    wgs84_crs = CRS.from_epsg(4326)

    # Create transformer from UTM to WGS84
    transformer = Transformer.from_crs(utm_crs, wgs84_crs, always_xy=True)
   
    # Generate x and y coordinates       
    x_coords = np.arange(min_x, max_x, cell_size)
    y_coords = np.arange(min_y, max_y, cell_size)

    # Create grid polygons
    mgrs_features = []
    for x in x_coords:
        for y in tqdm(y_coords, desc="Processing cells", unit="cells"):
            cell_polygon_utm = Polygon([
                (x, y),
                (x + cell_size, y),
                (x + cell_size, y + cell_size),
                (x, y + cell_size),
                (x, y)  # Close the polygon
            ])
            cell_polygon = utm_to_wgs84(cell_polygon_utm, transformer)
            centroid_lat, centroid_lon  =  cell_polygon.centroid.y, cell_polygon.centroid.x,
            mgrs_id = mgrs.toMgrs(centroid_lat, centroid_lon, resolution)            
            mgrs_feature = graticule_dggs_to_feature("mgrs",mgrs_id,resolution,cell_polygon)   

            # Load the GZD GeoJSON file
            gzd_json_path = os.path.join(os.path.dirname(__file__), 'gzd.geojson')            
            with open(gzd_json_path, 'r') as f:
                gzd_data = json.load(f)
            
            gzd_features = gzd_data["features"]
            gzd_feature = [feature for feature in gzd_features if feature["properties"].get("gzd") == gzd][0]
            
            intersected_polygon = mgrs_get_intersection(mgrs_feature, gzd_feature)   
            if intersected_polygon:
                intersected_centroid_lat, intersected_centroid_lon  =  intersected_polygon.centroid.y, intersected_polygon.centroid.x,
                interescted_mgrs_id = mgrs.toMgrs(intersected_centroid_lat, intersected_centroid_lon, resolution)            
                mgrs_feature = graticule_dggs_to_feature("mgrs",interescted_mgrs_id,resolution,intersected_polygon)
                mgrs_features.append(mgrs_feature)
           

    return {
        "type": "FeatureCollection",
        "features": mgrs_features
    }
    

def is_valid_gzd(gzd):
    """Check if a Grid Zone Designator (GZD) is valid."""
    pattern = r"^(?:0[1-9]|[1-5][0-9]|60)[C-HJ-NP-X]$"
    return bool(re.match(pattern, gzd))


def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Generate MGRS grid.")
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