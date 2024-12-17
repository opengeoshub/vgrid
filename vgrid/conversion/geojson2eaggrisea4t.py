import argparse
import json
from shapely.geometry import Polygon
from shapely.wkt import loads
from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.dggs_shape import DggsShape
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.dggs_shape_location import DggsShapeLocation
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.conversion.cell2geojson import fix_eaggr_wkt
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from vgrid.generator.eaggrisea4tgrid import cell_to_feature, cell_to_polygon, length_accuracy_dict,\
                                                get_children_cells_within_bbox
from tqdm import tqdm
from shapely.geometry import shape, Polygon, box, Point, LineString, mapping
from pyproj import Geod
import os

# Function to generate grid for Point
def point_to_grid(eaggr_dggs, resolution, point):
    features = []
   
    cell_id_len = resolution +2
    accuracy = length_accuracy_dict.get(cell_id_len)

    lat_long_point = LatLongPoint(point.y, point.x,accuracy)

    eaggr_cell = eaggr_dggs.convert_point_to_dggs_cell(lat_long_point)

    # Convert point to the seed cell
    cell_id = eaggr_cell.get_cell_id() # Unique identifier for the current cell
    cell_polygon = cell_to_polygon(eaggr_dggs,eaggr_cell)
    geod = Geod(ellps="WGS84")
    
    # Get the bounds and area of the cell
    min_x, min_y, max_x, max_y = cell_polygon.bounds
    center_lon = cell_polygon.centroid.x
    center_lat = cell_polygon.centroid.y
    cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
    _, _, cell_width = geod.inv(min_x, min_y, max_x, min_y)
    _, _, cell_height = geod.inv(min_x, min_y, min_x, max_y)
    
    features.append({
        "type": "Feature",
        "geometry": mapping(cell_polygon),
        "properties": {
            "isea4t": cell_id,
            "center_lat": center_lat,
            "center_lon": center_lon,
            "cell_width": cell_width,
            "cell_height": cell_height,
            "cell_area": cell_area
        },
    })
    
    return {
        "type": "FeatureCollection",
        "features": features,
    }


# Function to generate grid for Polyline
def polyline_to_grid(eaggr_dggs, resolution, geometry):
    features = []
    geod = Geod(ellps="WGS84")
    # Extract points from polyline
    if geometry.geom_type == 'LineString':
        # Handle single Polygon as before
        polylines = [geometry]
    elif geometry.geom_type == 'MultiLineString':
        # Handle MultiPolygon: process each polygon separately
        polylines = list(geometry)

    for polyline in polylines:
        cell_id_len = resolution +2
        accuracy = length_accuracy_dict.get(cell_id_len)
        bounding_box = box(*polyline.bounds)
        bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
        shapes = eaggr_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
        for shape in shapes:
            bbox_cells = shape.get_shape().get_outer_ring().get_cells()
            bounding_cell = eaggr_dggs.get_bounding_dggs_cell(bbox_cells)
            bounding_children_cells = get_children_cells_within_bbox(bounding_cell.get_cell_id(), bounding_box,resolution)
            features = []
            for cell in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
                cell_feature = cell_to_feature(DggsCell(cell))
                cell_shape = cell_to_polygon(DggsCell(cell))
                if cell_shape.intersects(geometry):
                    features.append(cell_feature)
            
            return {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": json.loads(json.dumps(feature["geometry"].__geo_interface__)),
                        "properties": feature["properties"],
                    }
                    for feature in features
                ]
            }

    return {
        "type": "FeatureCollection",
        "features": features,
    }
        
# Function to generate grid for Polygon
def polygon_to_grid(eaggr_dggs, resolution, geometry):
    features = []
    geod = Geod(ellps="WGS84")
    
    if geometry.geom_type == 'Polygon':
        # Handle single Polygon as before
        polygons = [geometry]
    elif geometry.geom_type == 'MultiPolygon':
        # Handle MultiPolygon: process each polygon separately
        polygons = list(geometry)

    for polygon in polygons:
        cell_id_len = resolution +2
        accuracy = length_accuracy_dict.get(cell_id_len)
        bounding_box = box(*polygon.bounds)
        bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
        shapes = eaggr_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
        for shape in shapes:
            bbox_cells = shape.get_shape().get_outer_ring().get_cells()
            bounding_cell = eaggr_dggs.get_bounding_dggs_cell(bbox_cells)
            bounding_children_cells = get_children_cells_within_bbox(bounding_cell.get_cell_id(), bounding_box,resolution)
            features = []
            for cell in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
                cell_feature = cell_to_feature(DggsCell(cell))
                cell_shape = cell_to_polygon(DggsCell(cell))
                if cell_shape.intersects(geometry):
                    features.append(cell_feature)
            
            return {
                "type": "FeatureCollection",
                "features": [
                    {
                        "type": "Feature",
                        "geometry": json.loads(json.dumps(feature["geometry"].__geo_interface__)),
                        "properties": feature["properties"],
                    }
                    for feature in features
                ]
            }

    return {
        "type": "FeatureCollection",
        "features": features,
    }
    
def get_bounding_box(geojson_file):
    # Load GeoJSON data
    with open(geojson_file, 'r') as f:
        data = json.load(f)
    
    # Initialize variables to store min/max bounds
    min_x, min_y, max_x, max_y = float('inf'), float('inf'), float('-inf'), float('-inf')
    
    # Iterate through each feature
    for feature in data['features']:
        geom = shape(feature['geometry'])  # Convert to Shapely geometry
        
        # Get bounds for the geometry
        geom_min_x, geom_min_y, geom_max_x, geom_max_y = geom.bounds
        
        # Update overall bounding box
        min_x = min(min_x, geom_min_x)
        min_y = min(min_y, geom_min_y)
        max_x = max(max_x, geom_max_x)
        max_y = max(max_y, geom_max_y)
    
    # Return the bounding box as a tuple
    return (min_x, min_y, max_x, max_y)

# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Generate EaggrISEA4T grid for shapes in GeoJSON format")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON string with Point, Polyline or Polygon"
    )
    args = parser.parse_args()
    geojson = args.geojson
    eaggr_dggs = Eaggr(Model.ISEA4T)

    resolution = args.resolution
    
    if resolution < 1 or resolution > 22:
        print(f"Please select a resolution in [1..22] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = []

    for feature in geojson_data['features']:      
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(eaggr_dggs, resolution, point)
                geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)  # Create Point for each coordinate set
                    point_features = point_to_grid(eaggr_dggs, resolution, point)
                    geojson_features.extend(point_features['features'])
        
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                # Directly process LineString geometry
                polyline = LineString(coordinates)
                polyline_features = polyline_to_grid(eaggr_dggs, resolution, polyline)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                # Iterate through each line in MultiLineString geometry
                for line_coords in coordinates:
                    polyline = LineString(line_coords)  # Use each part's coordinates
                    polyline_features = polyline_to_grid(eaggr_dggs, resolution, polyline)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                # Create Polygon with exterior and interior rings
                exterior_ring = coordinates[0]  # The first coordinate set is the exterior ring
                interior_rings = coordinates[1:]  # Remaining coordinate sets are interior rings (holes)
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = polygon_to_grid(eaggr_dggs, resolution, polygon)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = polygon_to_grid(eaggr_dggs, resolution, polygon)
                    geojson_features.extend(polygon_features['features'])

                    
    # Save the results to GeoJSON
    geojson_path = f"geojson2eaggrisea4t_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
