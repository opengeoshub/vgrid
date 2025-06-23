import argparse, json, os
from tqdm import tqdm
from shapely.geometry import box, Polygon, Point, LineString
from vgrid.generator.settings import geodesic_dggs_to_feature
import platform
import requests
from urllib.parse import urlparse

if (platform.system() == 'Windows'):
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
    from vgrid.generator.isea4tgrid import isea4t_cell_to_polygon, isea4t_res_accuracy_dict,\
                                            fix_isea4t_antimeridian_cells, get_isea4t_children_cells_within_bbox
    from vgrid.conversion.dggscompact import isea4t_compact
                                          
# Function to generate grid for Point
def point_to_grid(isea4t_dggs, resolution, point,feature_properties):
    isea4t_features = []   
    accuracy = isea4t_res_accuracy_dict.get(resolution)
    lat_long_point = LatLongPoint(point.y, point.x,accuracy)
    isea4t_cell = isea4t_dggs.convert_point_to_dggs_cell(lat_long_point)
    isea4t_id = isea4t_cell.get_cell_id() # Unique identifier for the current cell
    cell_polygon = isea4t_cell_to_polygon(isea4t_dggs,isea4t_cell)
    
    if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14') or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
            cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
    
    if cell_polygon:    
        num_edges = 3
        isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id,resolution,cell_polygon,num_edges)   
        isea4t_feature["properties"].update(feature_properties)
        isea4t_features.append(isea4t_feature)    
    
    return {
        "type": "FeatureCollection",
        "features": isea4t_features,
    }


def poly_to_grid(isea4t_dggs, resolution, geometry,feature_properties,compact=None):    
    isea4t_features = []

    if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon':
        polys = [geometry]
    elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'MultiPolygon':
        polys = list(geometry)

    for poly in polys:
        accuracy = isea4t_res_accuracy_dict.get(resolution)
        bounding_box = box(*poly.bounds)
        bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
        shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
        shape =  shapes[0]
        # for shape in shapes:
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_child_cells = get_isea4t_children_cells_within_bbox(isea4t_dggs,bounding_cell.get_cell_id(), bounding_box,resolution)
       
        if compact:
            bounding_child_cells = isea4t_compact(isea4t_dggs,bounding_child_cells)

        for child in bounding_child_cells:
            isea4t_cell = DggsCell(child)
            cell_polygon = isea4t_cell_to_polygon(isea4t_dggs,isea4t_cell)
            isea4t_id = isea4t_cell.get_cell_id()

            if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14') or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
                cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
            
            if cell_polygon.intersects(poly):
                num_edges = 3
                cell_resolution = len(isea4t_id)-2
                isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id,cell_resolution,cell_polygon,num_edges)   
                isea4t_feature["properties"].update(feature_properties)
                isea4t_features.append(isea4t_feature)          
               
    return {
        "type": "FeatureCollection",
        "features": isea4t_features,
    }

def geojson2isea4t(geojson_data, resolution, compact=False):
    """
    Convert GeoJSON data to ISEA4T DGGS format.
    
    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): Resolution level [0..25]
        compact (bool, optional): Enable ISEA4T compact mode for polygons. Defaults to False.
    
    Returns:
        dict: GeoJSON data in ISEA4T DGGS format
    """
    if platform.system() != 'Windows':
        raise NotImplementedError("ISEA4T DGGS conversion is only supported on Windows")
        
    if resolution < 0 or resolution > 25:
        raise ValueError("Resolution must be in range [0..25]")
        
    isea4t_dggs = Eaggr(Model.ISEA4T)
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing GeoJSON features"):
        feature_properties = feature['properties']    
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(isea4t_dggs, resolution, point, feature_properties)
                geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)
                    point_features = point_to_grid(isea4t_dggs, resolution, point, feature_properties)
                    geojson_features.extend(point_features['features'])
        
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(isea4t_dggs, resolution, polyline, feature_properties)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                for line_coords in coordinates:
                    polyline = LineString(line_coords)
                    polyline_features = poly_to_grid(isea4t_dggs, resolution, polyline, feature_properties)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                exterior_ring = coordinates[0]
                interior_rings = coordinates[1:]
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(isea4t_dggs, resolution, polygon, feature_properties, compact)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]
                    interior_rings = sub_polygon_coords[1:]
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(isea4t_dggs, resolution, polygon, feature_properties, compact)
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

def geojson2isea4t_cli():
    """
    Command-line interface for converting GeoJSON to ISEA4T DGGS format.
    Supports both local files and remote URLs.
    """
    parser = argparse.ArgumentParser(description="Convert GeoJSON to Open-Eaggr ISEA4T DGGS")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..25]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, 
        help="GeoJSON file path or URL (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable ISEA4T compact mode - for polygon only")

    if platform.system() != 'Windows':
        print("Error: ISEA4T DGGS conversion is only supported on Windows")
        return

    args = parser.parse_args()
    
    # Read GeoJSON data from file or URL
    geojson_data = read_geojson_file(args.geojson)
    if geojson_data is None:
        return

    try:
        result = geojson2isea4t(geojson_data, args.resolution, args.compact)
        
        # Save the results to GeoJSON
        geojson_name = os.path.splitext(os.path.basename(args.geojson))[0]
        geojson_path = f"{geojson_name}2isea4t_{args.resolution}.geojson"
        if args.compact:
            geojson_path = f"{geojson_name}2isea4t_{args.resolution}_compacted.geojson"
    
        with open(geojson_path, 'w') as f:
            json.dump(result, f)

        print(f"GeoJSON saved as {geojson_path}")
        
    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
