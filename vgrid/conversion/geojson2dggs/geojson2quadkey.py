from vgrid.utils import tilecode
from shapely.geometry import Point, LineString, Polygon
import argparse
import json
from tqdm import tqdm
import os
from vgrid.utils import mercantile
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.conversion.dggscompact import quadkeycompact

# Function to generate grid for Point
def point_to_grid(resolution, point, feature_properties):  
    quadkey_features = []
     # res: [0..29]        
    quadkey_id = tilecode.latlon2quadkey(point.y, point.x,resolution)
    quadkey_cell = mercantile.tile(point.x, point.y, resolution)
    bounds = mercantile.bounds(quadkey_cell)
    if bounds:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east
      
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        
        quadkey_feature = graticule_dggs_to_feature("quadkey",quadkey_id,resolution,cell_polygon)   
        quadkey_feature["properties"].update(feature_properties)

        quadkey_features.append(quadkey_feature)

    return {
        "type": "FeatureCollection",
        "features": quadkey_features
    }

# Function to generate grid for Polyline
def polyline_to_grid(resolution, geometry,feature_properties):
    quadkey_features = []
    # Extract points from polyline
    if geometry.geom_type == 'LineString':
        # Handle single Polygon as before
        polygons = [geometry]
    elif geometry.geom_type == 'MultiLineString':
        # Handle MultiPolygon: process each polygon separately
        polygons = list(geometry)

    for polygon in polygons:    
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        for tile in tiles:
            z, x, y = tile.z, tile.x, tile.y
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                # Create the bounding box coordinates for the polygon
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                quadkey_id = mercantile.quadkey(tile)
                
                cell_polygon = Polygon([
                    [min_lon, min_lat],  # Bottom-left corner
                    [max_lon, min_lat],  # Bottom-right corner
                    [max_lon, max_lat],  # Top-right corner
                    [min_lon, max_lat],  # Top-left corner
                    [min_lon, min_lat]   # Closing the polygon (same as the first point)
                ])
                if cell_polygon.intersects(polygon):
                    quadkey_feature = graticule_dggs_to_feature("quadkey",quadkey_id,resolution,cell_polygon) 
                    quadkey_feature["properties"].update(feature_properties)
                    quadkey_features.append(quadkey_feature)

    return {
        "type": "FeatureCollection",
        "features": quadkey_features
    }


def polygon_to_grid(resolution, geometry,feature_properties,compact):
    quadkey_features = []
    if geometry.geom_type == 'Polygon':
        # Handle single Polygon as before
        polygons = [geometry]
    elif geometry.geom_type == 'MultiPolygon':
        # Handle MultiPolygon: process each polygon separately
        polygons = list(geometry)

    for polygon in polygons:    
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        for tile in tiles:
            z, x, y = tile.z, tile.x, tile.y
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                # Create the bounding box coordinates for the polygon
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                quadkey_id = mercantile.quadkey(tile)
                
                cell_polygon = Polygon([
                    [min_lon, min_lat],  # Bottom-left corner
                    [max_lon, min_lat],  # Bottom-right corner
                    [max_lon, max_lat],  # Top-right corner
                    [min_lon, max_lat],  # Top-left corner
                    [min_lon, min_lat]   # Closing the polygon (same as the first point)
                ])
                if cell_polygon.intersects(polygon):
                    quadkey_feature = graticule_dggs_to_feature("quadkey",quadkey_id,resolution,cell_polygon) 
                    quadkey_feature["properties"].update(feature_properties)
                    quadkey_features.append(quadkey_feature)

    quadkey_geosjon = {
            "type": "FeatureCollection",
            "features": quadkey_features
    }
    if compact:
        return quadkeycompact(quadkey_geosjon)

    else: return quadkey_geosjon

def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to Quadkey Grid")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..29]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON file path (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable Tilecode compact mode")

    args = parser.parse_args()
    geojson = args.geojson
     # Initialize h3 DGGS
    resolution = args.resolution
    compact = args.compact  

    if resolution < 0 or resolution > 29:
        print(f"Please select a resolution in [0..29] range and try again ")
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
                polyline_features = polyline_to_grid(resolution, polyline,feature_properties)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                # Iterate through each line in MultiLineString geometry
                for line_coords in coordinates:
                    polyline = LineString(line_coords)  # Use each part's coordinates
                    polyline_features = polyline_to_grid(resolution, polyline,feature_properties)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                # Create Polygon with exterior and interior rings
                exterior_ring = coordinates[0]  # The first coordinate set is the exterior ring
                interior_rings = coordinates[1:]  # Remaining coordinate sets are interior rings (holes)
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = polygon_to_grid(resolution, polygon,feature_properties,compact)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = polygon_to_grid(resolution, polygon,feature_properties,compact)
                    geojson_features.extend(polygon_features['features'])

    # Save the results to GeoJSON
    geojson_name = os.path.splitext(os.path.basename(geojson))[0]
    geojson_path = f"{geojson_name}2quadkey_{resolution}.geojson"
    if compact:
        geojson_path = f"{geojson_name}2quadkey_{resolution}_compacted.geojson"
    
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
