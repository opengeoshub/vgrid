from vgrid.utils import mgrs
from shapely.geometry import Point, LineString, Polygon, mapping, shape
import argparse
import json
from tqdm import tqdm
import os
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.generator.geohashgrid import initial_geohashes, geohash_to_polygon, expand_geohash_bbox

from vgrid.conversion.dggs2geojson import mgrs_is_fully_within, mgrs_get_intersection
from vgrid.conversion.latlon2dggs import latlon2mgrs

def mgrs2geojson(mgrs_id):
    # Assuming mgrs.mgrscell returns cell bounds and origin
    min_lat, min_lon, max_lat, max_lon, resolution = mgrscell(mgrs_id)
    
    # Define the polygon coordinates for the MGRS cell
    cell_polygon = Polygon([
        (min_lon, min_lat),  # Bottom-left corner
        (max_lon, min_lat),  # Bottom-right corner
        (max_lon, max_lat),  # Top-right corner
        (min_lon, max_lat),  # Top-left corner
        (min_lon, min_lat)   # Closing the polygon
    ])   
    if cell_polygon.is_valid:
        mgrs_feature = graticule_dggs_to_feature("mgrs",mgrs_id,resolution,cell_polygon)
        
        try:
            # Load the GZD GeoJSON file
            gzd_json_path = os.path.join(os.path.dirname(__file__), './generator/gzd.geojson')
            
            with open(gzd_json_path, 'r') as f:
                gzd_data = json.load(f)
            
            gzd_features = gzd_data["features"]
            
            if mgrs_id[2] not in {"A", "B", "Y", "Z"}:
                if not mgrs_is_fully_within(mgrs_feature, gzd_features):
                    mgrs_feature = mgrs_get_intersection(mgrs_feature, gzd_features)
        except:
            pass    
        return mgrs_feature

import pyproj

def mgrscell(mgrs_code):
    """Get bounding box (min_lat, min_lon, max_lat, max_lon, precision) for an MGRS grid cell."""
    geod = pyproj.Geod(ellps="WGS84")

    min_lat, min_lon = mgrs.toWgs(mgrs_code)
    precision, grid_size = mgrs.get_precision_and_grid_size(mgrs_code)

    # Move north by grid_size meters
    max_lat, _, _ = geod.fwd(min_lon, min_lat, 0, grid_size)

    # Move east by grid_size meters
    _, max_lon, _ = geod.fwd(min_lon, min_lat, 90, grid_size)

    return min_lat, min_lon, max_lat, max_lon, precision


# Function to generate grid for Point
def point_to_grid(resolution, point, feature_properties):  
    geohash_features = []
    longitude = point.x
    latitude = point.y    
    mgrs_id = latlon2mgrs(latitude, longitude, resolution) 
    _,_,min_lat, min_lon, max_lat, max_lon, resolution = mgrs.mgrscell(mgrs_id)
    
    # Define the polygon coordinates for the MGRS cell
    cell_polygon = Polygon([
        (min_lon, min_lat),  # Bottom-left corner
        (max_lon, min_lat),  # Bottom-right corner
        (max_lon, max_lat),  # Top-right corner
        (min_lon, max_lat),  # Top-left corner
        (min_lon, min_lat)   # Closing the polygon
    ])

    mgrs_feature = graticule_dggs_to_feature("mgrs",mgrs_id,resolution,cell_polygon)
    
    try:
        # Load the GZD GeoJSON file
        gzd_json_path = os.path.join(os.path.dirname(__file__), './generator/gzd.geojson')
        
        with open(gzd_json_path, 'r') as f:
            gzd_data = json.load(f)
        
        gzd_features = gzd_data["features"]
        
        if mgrs_id[2] not in {"A", "B", "Y", "Z"}:
            if not mgrs_is_fully_within(mgrs_feature, gzd_features):
                mgrs_feature = mgrs_get_intersection(mgrs_feature, gzd_features)
    except:
        pass    
    
    mgrs_feature["properties"].update(feature_properties)
    geohash_features.append(mgrs_feature)

    return {
        "type": "FeatureCollection",
        "features": geohash_features
    }
       
# Function to generate grid for Polyline
def poly_to_grid(resolution, geometry,feature_properties):
    mgrs_features = []
    return {
        "type": "FeatureCollection",
        "features": mgrs_features
    }

def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to mgrs Grid")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..5]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Point, Polyline or Polygon in GeoJSON format"
    )
    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    
    if resolution < 0 or resolution > 5:
        print(f"Please select a resolution in [0..5] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, "r", encoding="utf-8") as f:
        try:
            geojson_data = json.load(f)  # Attempt to parse the JSON
        except json.JSONDecodeError as e:
            print(f"Invalid GeoJSON file: {e}")
            return

    
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing GeoJSON features"):
        feature_properties = feature['properties']
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)                
                point_features = point_to_grid(resolution, point,feature_properties)
                geojson_features.extend(point_features['features'])   


            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)  # Create Point for each coordinate set
                    point_features = point_to_grid(resolution, point)
                    geojson_features.extend(point_features['features'])
        
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                # Directly process LineString geometry
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(resolution, polyline,feature_properties)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                # Iterate through each line in MultiLineString geometry
                for line_coords in coordinates:
                    polyline = LineString(line_coords)  # Use each part's coordinates
                    polyline_features = poly_to_grid(resolution, polyline,feature_properties)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                # Create Polygon with exterior and interior rings
                exterior_ring = coordinates[0]  # The first coordinate set is the exterior ring
                interior_rings = coordinates[1:]  # Remaining coordinate sets are interior rings (holes)
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(resolution, polygon,feature_properties)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(resolution, polygon,feature_properties)
                    geojson_features.extend(polygon_features['features'])

    # Save the results to GeoJSON
    geojson_name = os.path.splitext(os.path.basename(geojson))[0]
    geojson_path = f"{geojson_name}2mgrs_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
