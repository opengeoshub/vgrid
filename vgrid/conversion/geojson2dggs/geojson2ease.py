import argparse, json, os
from tqdm import tqdm
from shapely.geometry import Polygon, box, Point, LineString
import requests
from urllib.parse import urlparse
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.utils.easedggs.constants import grid_spec, ease_crs, geo_crs, levels_specs
from vgrid.utils.easedggs.dggs.grid_addressing import grid_ids_to_geos, geos_to_grid_ids, geo_polygon_to_grid_ids
from vgrid.conversion.dggscompact import ease_compact

def point_to_grid(resolution, point, feature_properties):
    """
    Convert a Point geometry to EASE-DGGS grid cells.
    
    Args:
        resolution (int): Resolution level [0..6]
        point (Point): Shapely Point geometry
        feature_properties (dict): Properties to be added to the output features
        
    Returns:
        dict: GeoJSON FeatureCollection containing the grid cells
    """
    ease_features = []
    latitude = point.y
    longitude = point.x
    ease_cell = geos_to_grid_ids([(longitude,latitude)],level = resolution)
    ease_id = ease_cell['result']['data'][0]

    level = int(ease_id[1])  # Get the level (e.g., 'L0' -> 0)
    # Get level specs
    level_spec = levels_specs[level]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]
        
    geo = grid_ids_to_geos([ease_id])
    center_lon, center_lat = geo['result']['data'][0] 

    cell_min_lat = center_lat - (180 / (2 * n_row))
    cell_max_lat = center_lat + (180 / (2 * n_row))
    cell_min_lon = center_lon - (360 / (2 * n_col))
    cell_max_lon = center_lon + (360 / (2 * n_col))

    cell_polygon = Polygon([
        [cell_min_lon, cell_min_lat],
        [cell_max_lon, cell_min_lat],
        [cell_max_lon, cell_max_lat],
        [cell_min_lon, cell_max_lat],
        [cell_min_lon, cell_min_lat]
     ])

    if cell_polygon:   
        resolution = level
        num_edges = 4
        ease_feature = geodesic_dggs_to_feature("ease",ease_id,resolution,cell_polygon,num_edges)  
        ease_feature["properties"].update(feature_properties) 
        ease_features.append(ease_feature)

    return {
        "type": "FeatureCollection",
        "features": ease_features
    }

def poly_to_grid(resolution, geometry, feature_properties, compact=None):
    """
    Convert a Polygon or LineString geometry to EASE-DGGS grid cells.
    
    Args:
        resolution (int): Resolution level [0..6]
        geometry (Polygon/LineString): Shapely geometry
        feature_properties (dict): Properties to be added to the output features
        compact (bool, optional): Whether to enable EASE compact mode
        
    Returns:
        dict: GeoJSON FeatureCollection containing the grid cells
    """
    ease_features = []
    if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon':
        polys = [geometry]
    elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'MultiPolygon':
        polys = list(geometry)

    for poly in polys:
        poly_bbox = box(*poly.bounds)
        # Get all grid cells within the bounding box
        polygon_bbox_wkt = poly_bbox.wkt
        cells_bbox = geo_polygon_to_grid_ids(polygon_bbox_wkt, resolution, geo_crs, ease_crs, levels_specs, return_centroids = True, wkt_geom=True)
        ease_cells = cells_bbox['result']['data']   
        
        if compact:
            ease_cells = ease_compact(ease_cells)
             
        for ease_cell in ease_cells:
            cell_resolution =  int(ease_cell[1])  # Get the level (e.g., 'L0' -> 0)
            level_spec = levels_specs[cell_resolution]
            n_row = level_spec["n_row"]
            n_col = level_spec["n_col"]

            geo = grid_ids_to_geos([ease_cell])
            center_lon, center_lat = geo['result']['data'][0] 

            cell_min_lat = center_lat - (180 / (2 * n_row))
            cell_max_lat = center_lat + (180 / (2 * n_row))
            cell_min_lon = center_lon - (360 / (2 * n_col))
            cell_max_lon = center_lon + (360 / (2 * n_col))

            cell_polygon = Polygon([
                [cell_min_lon, cell_min_lat],
                [cell_max_lon, cell_min_lat],
                [cell_max_lon, cell_max_lat],
                [cell_min_lon, cell_max_lat],
                [cell_min_lon, cell_min_lat]
            ])
            if cell_polygon.intersects(poly):
                num_edges = 4
                ease_feature = geodesic_dggs_to_feature('ease', str(ease_cell), cell_resolution, cell_polygon, num_edges)
                ease_feature["properties"].update(feature_properties)
                ease_features.append(ease_feature)            
   
    return {
        "type": "FeatureCollection",
        "features": ease_features
    }

def geojson2ease(geojson_data, resolution, compact=False):
    """
    Convert GeoJSON data to EASE-DGGS format.
    
    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): Resolution level [0..6]
        compact (bool): Whether to enable EASE compact mode
        
    Returns:
        dict: GeoJSON data in EASE-DGGS format
        
    Raises:
        ValueError: If resolution is not in range [0..6]
    """
    if resolution < 0 or resolution > 6:
        raise ValueError("Resolution must be in range [0..6]")
        
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing features"):
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
        "features": geojson_features
    }

def is_url(path):
    """Check if the given path is a URL."""
    try:
        result = urlparse(path)
        return all([result.scheme, result.netloc])
    except:
        return False

def read_geojson_file(geojson_path):
    """Read GeoJSON from either a local file or URL."""
    if is_url(geojson_path):
        try:
            response = requests.get(geojson_path)
            response.raise_for_status()
            return json.loads(response.text)
        except requests.RequestException as e:
            print(f"Error: Failed to download GeoJSON from URL {geojson_path}: {str(e)}")
            return None
    else:
        if not os.path.exists(geojson_path):
            print(f"Error: The file {geojson_path} does not exist.")
            return None
        try:
            with open(geojson_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error reading GeoJSON file: {e}")
            return None

def geojson2ease_cli():
    """
    Command-line interface for converting GeoJSON to EASE-DGGS format.
    """
    parser = argparse.ArgumentParser(description="Convert GeoJSON to EASE-DGGS")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..6]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, 
        help="GeoJSON file path or URL (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable EASE compact mode")

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    compact = args.compact  

    try:
        # Load GeoJSON data from file or URL
        geojson_data = read_geojson_file(geojson)
        
        # Convert the GeoJSON data
        result = geojson2ease(geojson_data, resolution, compact)
        
        # Save the result
        # Extract base name from URL or file path
        geojson_name = os.path.splitext(os.path.basename(args.geojson))[0]
        geojson_path = f"{geojson_name}2ease_{resolution}.geojson"
        if compact:        
            geojson_path = f"{geojson_name}2ease_{resolution}_compacted.geojson"
            
        with open(geojson_path, 'w') as f:
            json.dump(result, f)

        print(f"GeoJSON saved as {geojson_path}")
        
    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
