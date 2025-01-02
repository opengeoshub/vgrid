from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from shapely.geometry import Polygon,mapping
from pyproj import Geod
from vgrid.utils.antimeridian import fix_polygon
import json
import argparse

isea3h_dggs = Eaggr(Model.ISEA3H)
geod = Geod(ellps="WGS84")

def cell_to_polygon(isea3h_cell):
    cell_to_shape =  isea3h_dggs.convert_dggs_cell_outline_to_shape_string(isea3h_cell, ShapeStringFormat.WKT)
    if cell_to_shape:
        coordinates_part = cell_to_shape.replace("POLYGON ((", "").replace("))", "")
        coordinates = []
        for coord_pair in coordinates_part.split(","):
            lon, lat = map(float, coord_pair.strip().split())
            coordinates.append([lon, lat])

        # Ensure the polygon is closed (first and last point must be the same)
        if coordinates[0] != coordinates[-1]:
            coordinates.append(coordinates[0])

    cell_polygon = Polygon(coordinates)
    fixed_polygon = fix_polygon(cell_polygon)    
    return fixed_polygon

accuracy_res_dict = {
        25503281086204.43: 0,
        629710644103.8047: 1,
        69967849344.8546: 2,
        7774205482.77106: 3,
        863800609.1842003: 4,
        95977845.45861907: 5,
        31992615.152873024: 6,
        131656.84875232293: 7,
        43885.62568888426: 8,
        14628.541896294753: 9,
        541.7947019742651: 10,
        60.196265293822194: 11,
        6.6821818482323785: 12,
        0.7361725765001773: 13,
        0.0849429895961743: 14,
        0.0:                15
        # 0.0:                16,                             
        # 0.0:                17,                             
        # 0.0:                18                   
        }

def children2geojson(isea3h):
    # Convert the input cell ID to a DggsCell object
    parent_cell = DggsCell(isea3h)
    # Get all parent cells for the input cell
    children = isea3h_dggs.get_dggs_cell_children(parent_cell)
    # List to store GeoJSON features
    features = []

    # Process each parent cell
    for child in children:
        # Convert parent cell to WKT shape string
        cell_polygon = cell_to_polygon(child)
        child_id = child.get_cell_id()
        # if child_id.startswith('00') or child_id.startswith('03') or child_id.startswith('09')\
        # or child_id.startswith('14') or child_id.startswith('19') or child_id.startswith('04') :
        #     cell_polygon = fix_isea3h_antimeridian_cells(cell_polygon)            
   
        cell_centroid = cell_polygon.centroid
        center_lat =  round(cell_centroid.y, 7)
        center_lon = round(cell_centroid.x, 7)
        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        
        isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(DggsCell(child_id))
        accuracy = isea3h2point._accuracy
        print (accuracy)
        avg_edge_len = round(cell_perimeter / 6,3)
        if (accuracy== 25503281086204.43): # icosahedron faces at resolution = 0
            avg_edge_len = round(cell_perimeter / 3,3)
        
        resolution  = accuracy_res_dict.get(accuracy)
        if accuracy == 0.0:
            # if avg_edge_len ==0.06:
            #     resolution = 15
            if avg_edge_len ==0.02:
                resolution = 16
            elif avg_edge_len ==0.01:
                resolution = 17
            elif avg_edge_len ==0.0:
                resolution = 18
        
        # Step 3: Construct the GeoJSON feature
        feature = {
            "type": "Feature",
            "geometry": mapping(cell_polygon),
            "properties": {
                    "isea3h": child_id,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "cell_area": cell_area,
                    "avg_edge_len": avg_edge_len,
                    "accuracy": accuracy,
                    "resolution": resolution
                    }
        }

        features.append(feature)

    # Construct and return the GeoJSON FeatureCollection
    feature_collection = {
        "type": "FeatureCollection",
        "features": features,
    }
    return feature_collection

def children2geojson_cli():
    """
    Command-line interface for isea3h2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert OpenEAGGR ISEA3H code to GeoJSON")
    parser.add_argument("isea3h", help="Input ISEA3H code, e.g., isea3h2geojson '07024,0'")
    args = parser.parse_args()
    geojson_data = json.dumps(children2geojson(args.isea3h))
    print(geojson_data)