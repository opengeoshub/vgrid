import os, argparse, json
from tqdm import tqdm
from shapely.geometry import Point, LineString, Polygon, mapping, box
import h3
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.settings import geodesic_dggs_to_feature
from pyproj import Geod
geod = Geod(ellps="WGS84")

# Function to generate grid for Point
def point_to_grid(resolution, point,feature_properties):
    h3_features = []
    # Convert point to the seed cell
    latitude = point.y
    longitude = point.x
    h3_id = h3.latlng_to_cell(latitude, longitude, resolution)
    
    cell_boundary = h3.cell_to_boundary(h3_id)     
    # Wrap and filter the boundary
    filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
    # Reverse lat/lon to lon/lat for GeoJSON compatibility
    reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
    cell_polygon = Polygon(reversed_boundary)
    if cell_polygon:
        num_edges = 6       
        if (h3.is_pentagon(h3_id)):
            num_edges = 5           
        h3_feature = geodesic_dggs_to_feature("h3",h3_id,resolution,cell_polygon,num_edges)   
        h3_feature["properties"].update(feature_properties)
        h3_features.append(h3_feature)    
    
    return {
        "type": "FeatureCollection",
        "features": h3_features,
    } 
        
def geodesic_buffer(polygon, distance):
    buffered_coords = []
    for lon, lat in polygon.exterior.coords:
        # Generate points around the current vertex to approximate a circle
        circle_coords = [
            geod.fwd(lon, lat, azimuth, distance)[:2]  # Forward calculation: returns (lon, lat, back_azimuth)
            for azimuth in range(0, 360, 10)  # Generate points every 10 degrees
        ]
        buffered_coords.append(circle_coords)
    
    # Flatten the list of buffered points and form a Polygon
    all_coords = [coord for circle in buffered_coords for coord in circle]
    return Polygon(all_coords).convex_hull

def poly_to_grid(resolution, geometry,feature_properties, compact= None):
    h3_features = []

    if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon' :
        polys = [geometry]
    elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'MultiPolygon' :
        polys = list(geometry)

    for poly in polys:
        bbox = box(*poly.bounds)
        distance = h3.average_hexagon_edge_length(resolution, unit='m') * 2
        bbox_buffer = geodesic_buffer(bbox, distance)
        bbox_buffer_cells = h3.geo_to_cells(bbox_buffer, resolution)
        if compact:
            bbox_buffer_cells = h3.compact_cells(bbox_buffer_cells)
            
        for bbox_buffer_cell in bbox_buffer_cells:
            cell_boundary = h3.cell_to_boundary(bbox_buffer_cell)
            filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            cell_polygon = Polygon(reversed_boundary)
            if cell_polygon.intersects(poly):
                cell_resolution = h3.get_resolution(bbox_buffer_cell)   
                num_edges = 6       
                if (h3.is_pentagon(bbox_buffer_cell)):
                    num_edges = 5           
                h3_feature = geodesic_dggs_to_feature("h3",bbox_buffer_cell,cell_resolution,cell_polygon,num_edges)   
                h3_feature["properties"].update(feature_properties)
                h3_features.append(h3_feature)               
          
    return {
        "type": "FeatureCollection",
        "features": h3_features,
    }

def geojson2h3(geojson_data, resolution, compact=False):
    """
    Convert GeoJSON data to H3 grid cells.
    
    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): H3 resolution [0..15]
        compact (bool): Enable H3 compact mode - for polygon only
        
    Returns:
        dict: GeoJSON FeatureCollection containing H3 grid cells
    """
    if resolution < 0 or resolution > 15:
        raise ValueError("Resolution must be in range [0..15]")
    
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing GeoJSON features"):
        feature_properties = feature['properties'] 
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(resolution, point, feature_properties)
                geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)
                    point_features = point_to_grid(resolution, point, feature_properties)
                    geojson_features.extend(point_features['features'])

        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(resolution, polyline, feature_properties)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                for line_coords in coordinates:
                    polyline = LineString(line_coords)
                    polyline_features = poly_to_grid(resolution, polyline, feature_properties)
                    geojson_features.extend(polyline_features['features'])

        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                exterior_ring = coordinates[0]
                interior_rings = coordinates[1:]
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(resolution, polygon, feature_properties, compact)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]
                    interior_rings = sub_polygon_coords[1:]
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(resolution, polygon, feature_properties, compact)
                    geojson_features.extend(polygon_features['features'])

    return {
        "type": "FeatureCollection",
        "features": geojson_features,
    }

def geojson2h3_cli():
    """
    Command-line interface for converting GeoJSON to H3 grid cells.
    """
    parser = argparse.ArgumentParser(description="Convert GeoJSON to H3 DGGS")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..15]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON file path (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable H3 compact mode - for polygon only")

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    compact = args.compact  

    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    try:
        result = geojson2h3(geojson_data, resolution, compact)
        
        geojson_name = os.path.splitext(os.path.basename(geojson))[0]
        geojson_path = f"{geojson_name}2h3_{resolution}.geojson"
        if compact:
            geojson_path = f"{geojson_name}2h3_{resolution}_compacted.geojson"
        
        with open(geojson_path, 'w') as f:
            json.dump(result, f)

        print(f"GeoJSON saved as {geojson_path}")
    except ValueError as e:
        print(f"Error: {str(e)}")
        
if __name__ == "__main__":
    geojson2h3_cli()