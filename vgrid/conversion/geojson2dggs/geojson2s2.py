from vgrid.utils import s2
from shapely.geometry import Point, LineString, Polygon
import argparse
import json
from tqdm import tqdm
import os
from vgrid.generator.s2grid import s2_cell_to_polygon
from vgrid.generator.settings import geodesic_dggs_to_feature

# Function to generate grid for Point
def point_to_grid(resolution, point,feature_properties):    
    s2_features = []
    # Convert point to the seed cell
    latitude = point.y
    longitude = point.x
    lat_lng = s2.LatLng.from_degrees(latitude, longitude)
    cell_id_max_res = s2.CellId.from_lat_lng(lat_lng)
    cell_id = cell_id_max_res.parent(resolution)
    s2_cell = s2.Cell(cell_id)
    cell_token = s2.CellId.to_token(s2_cell.id())    
    if s2_cell:
        cell_polygon = s2_cell_to_polygon(cell_id) # Fix antimeridian
        resolution = cell_id.level()
        num_edges = 4
        s2_feature = geodesic_dggs_to_feature("s2",cell_token,resolution,cell_polygon,num_edges)   
        s2_feature["properties"].update(feature_properties)
        s2_features.append(s2_feature)

    return {
        "type": "FeatureCollection",
        "features": s2_features
    }
        
def poly_to_grid(resolution, geometry,feature_properties,compact = None):
    s2_features = []
    
    if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon' :
        polys = [geometry]
    elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'MultiPolygon' :
        polys = list(geometry)

    for poly in polys:    
        min_lng, min_lat, max_lng, max_lat = poly.bounds
        # Define the cell level (S2 uses a level system for zoom, where level 30 is the highest resolution)
        level = resolution
        # Create a list to store the S2 cell IDs
        cell_ids = []
        # Define the cell covering
        coverer = s2.RegionCoverer()
        coverer.min_level = level
        coverer.max_level = level
        # coverer.max_cells = 1000_000  # Adjust as needed

        # Define the region to cover (in this example, we'll use the entire world)
        region = s2.LatLngRect(
            s2.LatLng.from_degrees(min_lat, min_lng),
            s2.LatLng.from_degrees(max_lat, max_lng)
        )

        # Get the covering cells
        covering = coverer.get_covering(region)
        cell_ids = covering  
        if compact:
            covering = s2.CellUnion(covering)
            covering.normalize()
            cell_ids = covering.cell_ids()  
            
        for cell_id in tqdm(cell_ids, desc="processing cells"):
            cell_polygon = s2_cell_to_polygon(cell_id)
          
            if cell_polygon.intersects(poly):
                cell_token = s2.CellId.to_token(cell_id)  
                resolution = cell_id.level()
                num_edges = 4
                s2_feature = geodesic_dggs_to_feature("s2",cell_token,resolution,cell_polygon,num_edges)   
                s2_feature["properties"].update(feature_properties)    
                s2_features.append(s2_feature)
                            
        return {
                "type": "FeatureCollection",
                "features": s2_features,
            }


def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to S2 Grid")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..30]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON file path (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable S2 compact mode - for polygon only")

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    compact = args.compact  
    
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
                    point_features = point_to_grid(resolution, point,feature_properties)
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
                polygon_features = poly_to_grid(resolution, polygon,feature_properties,compact)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(resolution, polygon,feature_properties,compact)
                    geojson_features.extend(polygon_features['features'])

    # Save the results to GeoJSON
    geojson_name = os.path.splitext(os.path.basename(geojson))[0]
    geojson_path = f"{geojson_name}2s2_{resolution}.geojson"
    if compact:
        geojson_path = f"{geojson_name}2s2_{resolution}_compacted.geojson"
    
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
