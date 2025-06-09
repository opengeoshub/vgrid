import argparse, json, os
from tqdm import tqdm
from shapely.geometry import Polygon, box, Point, LineString
from vgrid.generator.settings import geodesic_dggs_to_feature
import platform

if (platform.system() == 'Windows'): 
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
    from vgrid.generator.isea3hgrid import isea3h_cell_to_polygon, isea3h_res_accuracy_dict, isea3h_accuracy_res_dict, get_isea3h_children_cells_within_bbox
    from vgrid.conversion.dggscompact import isea3h_compact
    
from pyproj import Geod
geod = Geod(ellps="WGS84")
from shapely.geometry import Polygon,mapping


# Function to generate grid for Point
def point_to_grid(isea3h_dggs,resolution, point,feature_properties):
    if (platform.system() == 'Windows'):
        isea3h_features = []   
        accuracy = isea3h_res_accuracy_dict.get(resolution)
        lat_long_point = LatLongPoint(point.y, point.x, accuracy)
        isea3h_cell = isea3h_dggs.convert_point_to_dggs_cell(lat_long_point)
        cell_polygon = isea3h_cell_to_polygon(isea3h_dggs,isea3h_cell)
        
        if cell_polygon:
            isea3h_id = isea3h_cell.get_cell_id() 
            cell_resolution = resolution
            num_edges = 3 if cell_resolution == 0 else 6       
            isea4t_feature = geodesic_dggs_to_feature("isea3h",isea3h_id,cell_resolution,cell_polygon,num_edges)   
            isea4t_feature["properties"].update(feature_properties)
            isea3h_features.append(isea4t_feature)
       
        return {
            "type": "FeatureCollection",
            "features": isea3h_features,
        }

def poly_to_grid(isea3h_dggs,resolution, geometry, feature_properties,compact):    
    if (platform.system() == 'Windows'):
        isea3h_features = []
        if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon':
            polys = [geometry]
        elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'MultiPolygon':
            polys = list(geometry)

        for poly in polys:
            accuracy = isea3h_res_accuracy_dict.get(resolution)
            bounding_box = box(*poly.bounds)
            bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
            shapes = isea3h_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
            shape =  shapes[0]
            # for shape in shapes:
            bbox_cells = shape.get_shape().get_outer_ring().get_cells()
            bounding_cell = isea3h_dggs.get_bounding_dggs_cell(bbox_cells)
            bounding_child_cells = get_isea3h_children_cells_within_bbox(isea3h_dggs,bounding_cell.get_cell_id(), bounding_box,resolution)
            if compact:
                bounding_child_cells = isea3h_compact(isea3h_dggs,bounding_child_cells)

            for child in bounding_child_cells:
                isea3h_cell = DggsCell(child)
                cell_polygon = isea3h_cell_to_polygon(isea3h_dggs,isea3h_cell)
                if cell_polygon.intersects(poly):
                    isea3h_id = isea3h_cell.get_cell_id()
                    isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(isea3h_cell)      
                    cell_accuracy = isea3h2point._accuracy        
                    cell_resolution  = isea3h_accuracy_res_dict.get(cell_accuracy)                    
                    num_edges = 3 if cell_resolution == 0 else 6  
                         
                    isea4t_feature = geodesic_dggs_to_feature("isea3h",isea3h_id,cell_resolution,cell_polygon,num_edges)   
                    isea4t_feature["properties"].update(feature_properties)
                    isea3h_features.append(isea4t_feature)
       
        return {
            "type": "FeatureCollection",
            "features": isea3h_features,
        }
            
def geojson2isea3h(geojson_data, resolution, compact=False):
    """
    Convert GeoJSON data to ISEA3H DGGS format.
    
    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): Resolution level [0..32]
        compact (bool): Whether to enable ISEA3H compact mode
        
    Returns:
        dict: GeoJSON data in ISEA3H DGGS format
    """
    if platform.system() != 'Windows':
        raise NotImplementedError("ISEA3H DGGS conversion is only supported on Windows")
        
    if resolution < 0 or resolution > 32:
        raise ValueError("Resolution must be in range [0..32]")
        
    isea3h_dggs = Eaggr(Model.ISEA3H)
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing GeoJSON features"):   
        feature_properties = feature['properties'] 
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(isea3h_dggs, resolution, point, feature_properties)
                geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)
                    point_features = point_to_grid(isea3h_dggs, resolution, point, feature_properties)
                    geojson_features.extend(point_features['features'])
        
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(isea3h_dggs, resolution, polyline, feature_properties, compact)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                for line_coords in coordinates:
                    polyline = LineString(line_coords)
                    polyline_features = poly_to_grid(isea3h_dggs, resolution, polyline, feature_properties, compact)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                exterior_ring = coordinates[0]
                interior_rings = coordinates[1:]
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(isea3h_dggs, resolution, polygon, feature_properties, compact)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]
                    interior_rings = sub_polygon_coords[1:]
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(isea3h_dggs, resolution, polygon, feature_properties, compact)
                    geojson_features.extend(polygon_features['features'])

    return {
        "type": "FeatureCollection",
        "features": geojson_features
    }

def geojson2isea3h_cli():
    """Command line interface for converting GeoJSON to ISEA3H DGGS format."""
    parser = argparse.ArgumentParser(description="Convert GeoJSON to Open-Eaggr ISEA3H DGGS")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..32]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON file path (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable ISEA3H compact mode")

    args = parser.parse_args()
    
    if not os.path.exists(args.geojson):
        print(f"Error: The file {args.geojson} does not exist.")
        return

    try:
        with open(args.geojson, 'r', encoding='utf-8') as f:
            geojson_data = json.load(f)
            
        result = geojson2isea3h(geojson_data, args.resolution, args.compact)
        
        geojson_name = os.path.splitext(os.path.basename(args.geojson))[0]
        geojson_path = f"{geojson_name}2isea3h_{args.resolution}.geojson"
        if args.compact:
            geojson_path = f"{geojson_name}2isea3h_{args.resolution}_compacted.geojson"

        with open(geojson_path, 'w') as f:
            json.dump(result, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
        
    except Exception as e:
        print(f"Error: {str(e)}")
        
if __name__ == "__main__":
    geojson2isea3h_cli()