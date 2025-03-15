from vgrid.utils import tilecode
from shapely.geometry import Point, LineString, Polygon, mapping, box
import argparse
import json
from tqdm import tqdm
import os
from vgrid.utils import mercantile
from pyproj import Geod
geod = Geod(ellps="WGS84")


# Function to generate grid for Point
def point_to_grid(resolution, point):    
    features = []
     # res: [0..29]        
    tilecode_id = tilecode.latlon2tilecode(point.y, point.x,resolution)
    tilecode_cell = mercantile.tile(point.x, point.y, resolution)
    bounds = mercantile.bounds(tilecode_cell)
    if bounds:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east
        
        quadkey = mercantile.quadkey(tilecode_cell)

        center_lat = round((min_lat + max_lat) / 2,7)
        center_lon = round((min_lon + max_lon) / 2,7)
        
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)  # Area in square meters     
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        avg_edge_len = round(cell_perimeter/6,2)
        resolution = tilecode_cell.z 
              
        features.append({
                "type": "Feature",
                "geometry": mapping(cell_polygon),          
                "properties": {
                    "tilecode": tilecode_id,  
                    "quadkey": quadkey,
                    "resolution": resolution,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "avg_edge_len": avg_edge_len,
                    "cell_area": cell_area                
                }
            })

    geojson_features = {
        'type': 'FeatureCollection',
        'features': features
    }

    return geojson_features

# Function to generate grid for Polyline
def polyline_to_grid(resolution, geometry):
    features = []
    # Extract points from polyline
    if geometry.geom_type == 'LineString':
        # Handle single Polygon as before
        polylines = [geometry]
    elif geometry.geom_type == 'MultiLineString':
        # Handle MultiPolygon: process each polygon separately
        polylines = list(geometry)

    for polyline in polylines:    
        min_lon, min_lat, max_lon, max_lat = polyline.bounds
        tiles = mercantile.tiles(min_lon, min_lat, max_lon, max_lat, resolution)
        for tile in tiles:
            z, x, y = tile.z, tile.x, tile.y
            tilecode_id = f"z{tile.z}x{tile.x}y{tile.y}"
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                # Create the bounding box coordinates for the polygon
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                
                quadkey = mercantile.quadkey(tile)

                center_lat = round((min_lat + max_lat) / 2,7)
                center_lon = round((min_lon + max_lon) / 2,7)
                
                cell_polygon = Polygon([
                    [min_lon, min_lat],  # Bottom-left corner
                    [max_lon, min_lat],  # Bottom-right corner
                    [max_lon, max_lat],  # Top-right corner
                    [min_lon, max_lat],  # Top-left corner
                    [min_lon, min_lat]   # Closing the polygon (same as the first point)
                ])
                if cell_polygon.intersects(polyline):
                    cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)  # Area in square meters     
                    cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
                    avg_edge_len = round(cell_perimeter/6,2)        
                    features.append({
                        "type": "Feature",
                        "geometry": mapping(cell_polygon),          
                        "properties": {
                            "tilecode": tilecode_id,  
                            "quadkey": quadkey,
                            "resolution": z, 
                            "center_lat": center_lat,
                            "center_lon": center_lon,
                            "avg_edge_len": avg_edge_len,
                            "cell_area": cell_area
                        },
                    })     

        geojson_features = {
            'type': 'FeatureCollection',
            'features': features
        }

    return geojson_features

# # Function to generate grid for Polygon
def polygon_to_grid(resolution, geometry):
    features = []
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
            tilecode_id = f"z{tile.z}x{tile.x}y{tile.y}"
            bounds = mercantile.bounds(x, y, z)
            if bounds:
                # Create the bounding box coordinates for the polygon
                min_lat, min_lon = bounds.south, bounds.west
                max_lat, max_lon = bounds.north, bounds.east
                
                quadkey = mercantile.quadkey(tile)

                center_lat = round((min_lat + max_lat) / 2,7)
                center_lon = round((min_lon + max_lon) / 2,7)
                
                cell_polygon = Polygon([
                    [min_lon, min_lat],  # Bottom-left corner
                    [max_lon, min_lat],  # Bottom-right corner
                    [max_lon, max_lat],  # Top-right corner
                    [min_lon, max_lat],  # Top-left corner
                    [min_lon, min_lat]   # Closing the polygon (same as the first point)
                ])
                if cell_polygon.intersects(polygon):
                    cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)  # Area in square meters     
                    # Calculate width (longitude difference at a constant latitude)
                    cell_width = round(geod.line_length([min_lon, max_lon], [min_lat, min_lat]),2)
                    # Calculate height (latitude difference at a constant longitude)
                    cell_height = round(geod.line_length([min_lon, min_lon], [min_lat, max_lat]),2)


                    features.append({
                        "type": "Feature",
                        "geometry": mapping(cell_polygon),          
                        "properties": {
                            "tilecode": tilecode_id,  
                            "quadkey": quadkey,
                            "center_lat": center_lat,
                            "center_lon": center_lon,
                            "cell_area": cell_area,
                            "cell_width": cell_width,
                            "cell_height": cell_height,
                            "resolution": z  
                        },
                    })     

        geojson_features = {
            'type': 'FeatureCollection',
            'features': features
        }

    return geojson_features


# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to Tilecode Grid")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..29]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Point, Polyline or Polygon in GeoJSOn format"
    )
    args = parser.parse_args()
    geojson = args.geojson
     # Initialize h3 DGGS
    resolution = args.resolution
    
    if resolution < 0 or resolution > 29:
        print(f"Please select a resolution in [0..29] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing GeoJSON features"):
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(resolution, point)
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
                polyline_features = polyline_to_grid(resolution, polyline)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                # Iterate through each line in MultiLineString geometry
                for line_coords in coordinates:
                    polyline = LineString(line_coords)  # Use each part's coordinates
                    polyline_features = polyline_to_grid(resolution, polyline)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                # Create Polygon with exterior and interior rings
                exterior_ring = coordinates[0]  # The first coordinate set is the exterior ring
                interior_rings = coordinates[1:]  # Remaining coordinate sets are interior rings (holes)
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = polygon_to_grid(resolution, polygon)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = polygon_to_grid(resolution, polygon)
                    geojson_features.extend(polygon_features['features'])

    # Save the results to GeoJSON
    geojson_path = f"geojson2tilecode_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
