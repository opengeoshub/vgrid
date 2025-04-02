import argparse, json, os
from tqdm import tqdm
from shapely.geometry import box, Polygon, Point, LineString
from vgrid.generator.settings import geodesic_dggs_to_feature
import platform

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


# Function to generate grid for Polyline
def polyline_to_grid(isea4t_dggs, resolution, geometry,feature_properties):    
    isea4t_features = []
    if geometry.geom_type == 'LineString':
        polylines = [geometry]
    elif geometry.geom_type == 'MultiLineString':
        polylines = list(geometry)

    for polyline in polylines:
        accuracy = isea4t_res_accuracy_dict.get(resolution)
        bounding_box = box(*polyline.bounds)
        bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
        shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
        shape =  shapes[0]
        # for shape in shapes:
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_children_cells = get_isea4t_children_cells_within_bbox(isea4t_dggs,bounding_cell.get_cell_id(), bounding_box,resolution)
        for child in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
            isea4t_cell = DggsCell(child)
            cell_polygon = isea4t_cell_to_polygon(isea4t_dggs,isea4t_cell)
            isea4t_id = isea4t_cell.get_cell_id()

            if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14') or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
                cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
            
            if cell_polygon.intersects(polyline):
                num_edges = 3
                isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id,resolution,cell_polygon,num_edges)   
                isea4t_feature["properties"].update(feature_properties)
                isea4t_features.append(isea4t_feature)          
               
    return {
        "type": "FeatureCollection",
        "features": isea4t_features,
    }

# Function to generate grid for Polygon
def polygon_to_grid(isea4t_dggs, resolution, geometry,feature_properties,compact):    
    isea4t_features = []
    if geometry.geom_type == 'Polygon':
        polygons = [geometry]
    elif geometry.geom_type == 'MultiPolygon':
        polygons = list(geometry)

    for polygon in polygons:
        accuracy = isea4t_res_accuracy_dict.get(resolution)
        bounding_box = box(*polygon.bounds)
        bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
        shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
        shape =  shapes[0]
        # for shape in shapes:
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_children_cells = get_isea4t_children_cells_within_bbox(isea4t_dggs,bounding_cell.get_cell_id(), bounding_box,resolution)
       
        if compact:
            bounding_children_cells = isea4t_compact(isea4t_dggs,bounding_children_cells)

        for child in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
            isea4t_cell = DggsCell(child)
            cell_polygon = isea4t_cell_to_polygon(isea4t_dggs,isea4t_cell)
            isea4t_id = isea4t_cell.get_cell_id()

            if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14') or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
                cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
            
            if cell_polygon.intersects(polygon):
                num_edges = 3
                isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id,resolution,cell_polygon,num_edges)   
                isea4t_feature["properties"].update(feature_properties)
                isea4t_features.append(isea4t_feature)          
               
    return {
        "type": "FeatureCollection",
        "features": isea4t_features,
    }


def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to Open-Eaggr ISEA4T Grid")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..25]")
    # actual resolution range: [0..39]
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON string with Point, Polyline or Polygon"
    )
    parser.add_argument('-compact', action='store_true', help="Enable ISEA4T compact mode")

    if (platform.system() == 'Windows'):
        isea4t_dggs = Eaggr(Model.ISEA4T)
        args = parser.parse_args()
        geojson = args.geojson
        resolution = args.resolution
        compact = args.compact  

        if resolution < 0 or resolution > 25:
            print(f"Please select a resolution in [0..25] range and try again ")
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
                    point_features = point_to_grid(isea4t_dggs,resolution, point,feature_properties)
                    geojson_features.extend(point_features['features'])

                elif feature['geometry']['type'] == 'MultiPoint':
                    for point_coords in coordinates:
                        point = Point(point_coords)  # Create Point for each coordinate set
                        point_features = point_to_grid(isea4t_dggs,resolution, point,feature_properties)
                        geojson_features.extend(point_features['features'])
            
            elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
                coordinates = feature['geometry']['coordinates']
                if feature['geometry']['type'] == 'LineString':
                    # Directly process LineString geometry
                    polyline = LineString(coordinates)
                    polyline_features = polyline_to_grid(isea4t_dggs,resolution, polyline,feature_properties)
                    geojson_features.extend(polyline_features['features'])

                elif feature['geometry']['type'] == 'MultiLineString':
                    # Iterate through each line in MultiLineString geometry
                    for line_coords in coordinates:
                        polyline = LineString(line_coords)  # Use each part's coordinates
                        polyline_features = polyline_to_grid(isea4t_dggs,resolution, polyline,feature_properties)
                        geojson_features.extend(polyline_features['features'])
                
            elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
                coordinates = feature['geometry']['coordinates']

                if feature['geometry']['type'] == 'Polygon':
                    # Create Polygon with exterior and interior rings
                    exterior_ring = coordinates[0]  # The first coordinate set is the exterior ring
                    interior_rings = coordinates[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = polygon_to_grid(isea4t_dggs,resolution, polygon,feature_properties,compact)
                    geojson_features.extend(polygon_features['features'])

                elif feature['geometry']['type'] == 'MultiPolygon':
                    # Handle each sub-polygon in MultiPolygon geometry
                    for sub_polygon_coords in coordinates:
                        exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                        interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                        polygon = Polygon(exterior_ring, interior_rings)
                        polygon_features = polygon_to_grid(isea4t_dggs,resolution, polygon,feature_properties,compact)
                        geojson_features.extend(polygon_features['features'])

        geojson_name = os.path.splitext(os.path.basename(geojson))[0]
        geojson_path = f"{geojson_name}2isea4t_{resolution}.geojson"
        if compact:
            geojson_path = f"{geojson_name}2isea4t_{resolution}_compacted.geojson"
    

        with open(geojson_path, 'w') as f:
            json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")

if __name__ == "__main__":
    main()
