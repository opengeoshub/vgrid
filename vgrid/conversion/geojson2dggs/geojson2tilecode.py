from vgrid.utils import tilecode
from shapely.geometry import Point, LineString, Polygon
import argparse
import json
from tqdm import tqdm
import os
from vgrid.utils import mercantile
from vgrid.generator.settings import graticule_dggs_to_feature
from vgrid.conversion.dggscompact import tilecodecompact
import re

# Function to generate grid for Point
def point_to_grid(resolution, point, feature_properties):  
    tilecode_features = []
     # res: [0..29]        
    tilecode_id = tilecode.latlon2tilecode(point.y, point.x,resolution)
    tilecode_cell = mercantile.tile(point.x, point.y, resolution)
    bounds = mercantile.bounds(tilecode_cell)
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
        
        tilecode_feature = graticule_dggs_to_feature("tilecode",tilecode_id,resolution,cell_polygon)   
        tilecode_feature["properties"].update(feature_properties)

        tilecode_features.append(tilecode_feature)

    return {
        "type": "FeatureCollection",
        "features": tilecode_features
    }

def polyline_to_grid(resolution, geometry,feature_properties):
    tilecode_features = []
    # Extract points from polyline
    if geometry.geom_type == 'LineString':
        # Handle single Polygon as before
        polylines = [geometry]
    elif geometry.geom_type == 'MultiLineString' :
        # Handle MultiPolygon: process each polygon separately
        polylines = list(geometry)

    for polyline in polylines:    
        min_lon, min_lat, max_lon, max_lat = polyline.bounds
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        for tile in tqdm(tiles, desc="Processing cells",unit=" cells"):
            z, x, y = tile.z, tile.x, tile.y
            tilecode_id = f"z{tile.z}x{tile.x}y{tile.y}"
            bounds = mercantile.bounds(x, y, z)
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
                if cell_polygon.intersects(polyline):
                    tilecode_feature = graticule_dggs_to_feature("tilecode",tilecode_id,resolution,cell_polygon) 
                    tilecode_feature["properties"].update(feature_properties)
                    tilecode_features.append(tilecode_feature)

    return {
        "type": "FeatureCollection",
        "features": tilecode_features
    }

def polygon_to_grid(resolution, geometry,feature_properties,compact):
    tilecode_features = []
    # Extract points from polyline
    if geometry.geom_type == 'Polygon':
        # Handle single Polygon as before
        polygons = [geometry]
    elif geometry.geom_type == 'MultiPolygon':
        # Handle MultiPolygon: process each polygon separately
        polygons = list(geometry)

    for polygon in polygons:    
        min_lon, min_lat, max_lon, max_lat = polygon.bounds
        tilecodes = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        tilecode_ids = []
        for tilecode in tilecodes:
            tilecode_id = f"z{tilecode.z}x{tilecode.x}y{tilecode.y}"
            tilecode_ids.append(tilecode_id)
       
        for tilecode_id in tqdm(tilecode_ids, desc="Processing cells",unit=" cells"):
            match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id)
            if not match:
                raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")
            cell_resolution= int(match.group(1))
              # Convert matched groups to integers
            z = int(match.group(1))
            x = int(match.group(2))
            y = int(match.group(3))

            # Get the bounds of the tile in (west, south, east, north)
            bounds = mercantile.bounds(x, y, z)
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
                if cell_polygon.intersects(polygon):
                    tilecode_feature = graticule_dggs_to_feature("tilecode",tilecode_id,cell_resolution,cell_polygon) 
                    tilecode_feature["properties"].update(feature_properties)
                    tilecode_features.append(tilecode_feature)

        tilecode_geosjon = {
            "type": "FeatureCollection",
            "features": tilecode_features
        }

        if compact:
            return tilecodecompact(tilecode_geosjon)


    return tilecode_geosjon


def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to Tilecode Grid")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..29]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON file path (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable Tilecode compact mode")

    args = parser.parse_args()
    geojson = args.geojson
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
    geojson_path = f"{geojson_name}2tilecode_{resolution}.geojson"
    if compact:
        geojson_path = f"{geojson_name}2s2_{resolution}_compacted.geojson"
    
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
