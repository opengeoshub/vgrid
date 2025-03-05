from vgrid.utils import s2
from shapely.geometry import Point, LineString, Polygon, mapping, box
import argparse
import json
from tqdm import tqdm
import os
from vgrid.generator.s2grid import cell_to_polygon
from pyproj import Geod
geod = Geod(ellps="WGS84")
 
# Function to generate grid for Point
def point_to_grid(resolution, point):    
    features = []
    # Convert point to the seed cell
    latitude = point.y
    longitude = point.x
    lat_lng = s2.LatLng.from_degrees(latitude, longitude)
    cell_id_max_res = s2.CellId.from_lat_lng(lat_lng)
    cell_id = cell_id_max_res.parent(resolution)
    s2_cell = s2.Cell(cell_id)
    cell_token = s2.CellId.to_token(s2_cell.id())
    
    if s2_cell:
        cell_polygon = cell_to_polygon(cell_id) # Fix antimeridian
        lat_lng = cell_id.to_lat_lng()            
        # Extract latitude and longitude in degrees
        center_lat = round(lat_lng.lat().degrees,7)
        center_lon = round(lat_lng.lng().degrees,7)

        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2) # Area in square meters     
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])  # Perimeter in meters  
        avg_edge_len = round(cell_perimeter/4,2)

        features.append({
                "type": "Feature",
                "geometry": mapping(cell_polygon),
                "properties": {
                    "s2_token": cell_token,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "area": cell_area,
                    "avg_edge_len": avg_edge_len,
                    "resolution": cell_id.level()
                        },
             })
            
        return {
            "type": "FeatureCollection",
            "features": features,
        }


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
        min_lng, min_lat, max_lng, max_lat = polyline.bounds
        # Define the cell level (S2 uses a level system for zoom, where level 30 is the highest resolution)
        level = resolution
        # Create a list to store the S2 cell IDs
        cell_ids = []
        # Define the cell covering
        coverer = s2.RegionCoverer()
        coverer.min_level = level
        coverer.max_level = level
        # coverer.max_cells = 1000_000  # Adjust as needed
        # coverer.max_cells = 0  # Adjust as needed

        # Define the region to cover (in this example, we'll use the entire world)
        region = s2.LatLngRect(
            s2.LatLng.from_degrees(min_lat, min_lng),
            s2.LatLng.from_degrees(max_lat, max_lng)
        )

        # Get the covering cells
        covering = coverer.get_covering(region)

        # Convert the covering cells to S2 cell IDs
        for cell_id in covering:
            cell_ids.append(cell_id)

        for cell_id in tqdm(cell_ids, desc="processing cells"):
            cell_polygon = cell_to_polygon(cell_id)
            lat_lng = cell_id.to_lat_lng()      
            cell_token = s2.CellId.to_token(cell_id)      
            # Extract latitude and longitude in degrees
            center_lat = round(lat_lng.lat().degrees,7)
            center_lon = round(lat_lng.lng().degrees,7)

            cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2) # Area in square meters     
            cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])  # Perimeter in meters  
            avg_edge_len = round(cell_perimeter/4,2)
                        
            if cell_polygon.intersects(polyline):
                features.append({
                    "type": "Feature",
                    "geometry": mapping(cell_polygon),
                    "properties": {
                            "s2_token": cell_token,
                            "center_lat": center_lat,
                            "center_lon": center_lon,
                            "area": cell_area,
                            "avg_edge_len": avg_edge_len,
                            "resolution": cell_id.level()
                            },
                })
            
        # Create a FeatureCollection
        return {
                "type": "FeatureCollection",
                "features": features,
            }

# Function to generate grid for Polygon
def polygon_to_grid(resolution, geometry):
    features = []
    if geometry.geom_type == 'Polygon':
        # Handle single Polygon as before
        polygons = [geometry]
    elif geometry.geom_type == 'MultiPolygon':
        # Handle MultiPolygon: process each polygon separately
        polygons = list(geometry)

    for polygon in polygons:
        min_lng, min_lat, max_lng, max_lat = polygon.bounds
        # Define the cell level (S2 uses a level system for zoom, where level 30 is the highest resolution)
        level = resolution
        # Create a list to store the S2 cell IDs
        cell_ids = []
        # Define the cell covering
        coverer = s2.RegionCoverer()
        coverer.min_level = level
        coverer.max_level = level
        # coverer.max_cells = 1000_000  # Adjust as needed
        # coverer.max_cells = 0  # Adjust as needed

        # Define the region to cover (in this example, we'll use the entire world)
        region = s2.LatLngRect(
            s2.LatLng.from_degrees(min_lat, min_lng),
            s2.LatLng.from_degrees(max_lat, max_lng)
        )


        # Get the covering cells
        covering = coverer.get_covering(region)

        # Convert the covering cells to S2 cell IDs
        for cell_id in covering:
            cell_ids.append(cell_id)

        for cell_id in tqdm(cell_ids, desc="processing cells"):
            cell_polygon = cell_to_polygon(cell_id)
            lat_lng = cell_id.to_lat_lng()      
            cell_token = s2.CellId.to_token(cell_id)      
            # Extract latitude and longitude in degrees
            center_lat = round(lat_lng.lat().degrees,7)
            center_lon = round(lat_lng.lng().degrees,7)

            cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2) # Area in square meters     
            cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])  # Perimeter in meters  
            avg_edge_len = round(cell_perimeter/4,2)
                        
            if cell_polygon.intersects(polygon):
                features.append({
                    "type": "Feature",
                    "geometry": mapping(cell_polygon),
                    "properties": {
                            "s2_token": cell_token,
                            "center_lat": center_lat,
                            "center_lon": center_lon,
                            "area": cell_area,
                            "avg_edge_len": avg_edge_len,
                            "resolution": cell_id.level()
                            },
                })
            
        # Create a FeatureCollection
        return {
                "type": "FeatureCollection",
                "features": features,
            }


# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Generate S2 grid for shapes in GeoJSON format")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..30]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON string with Point, Polyline or Polygon"
    )
    args = parser.parse_args()
    geojson = args.geojson
     # Initialize h3 DGGS
    resolution = args.resolution
    
    if resolution < 0 or resolution > 30:
        print(f"Please select a resolution in [0..30] range and try again ")
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
    geojson_path = f"geojson2s2_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
