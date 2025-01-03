from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from shapely.geometry import Polygon, mapping
from pyproj import Geod
import json, argparse
from vgrid.utils.antimeridian import fix_polygon

isea3h_dggs = Eaggr(Model.ISEA3H)
from pyproj import Geod
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
        # 10664205.060395785: 6,
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

def siblings2geojson(isea3h):
    input_cell = DggsCell(isea3h)
    siblings = isea3h_dggs.get_dggs_cell_siblings(input_cell)
    features = []

    # Process each parent cell
    for sibling in siblings:
        # Convert parent cell to WKT shape string
        cell_polygon = cell_to_polygon(sibling)
        sibling_id = sibling.get_cell_id()
        # if sibling_id.startswith('00') or sibling_id.startswith('03') or sibling_id.startswith('09')\
        # or sibling_id.startswith('14') or sibling_id.startswith('19') or sibling_id.startswith('04') :
        #     cell_polygon = fix_isea3h_antimeridian_cells(cell_polygon)            
   
        cell_centroid = cell_polygon.centroid
        center_lat =  round(cell_centroid.y, 7)
        center_lon = round(cell_centroid.x, 7)
        cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(DggsCell(isea3h))      
        
        accuracy = isea3h2point._accuracy
            
        avg_edge_len = cell_perimeter / 6
        if (accuracy== 25_503_281_086_204.43): # icosahedron faces at resolution = 0
            avg_edge_len = cell_perimeter / 3
        
        resolution  = accuracy_res_dict.get(accuracy)
        if accuracy == 0.0:
            if round(avg_edge_len,2) == 0.06:
                resolution = 33
            elif round(avg_edge_len,2) == 0.03:
                resolution = 34
            elif round(avg_edge_len,2) == 0.02:
                resolution = 35
            elif round(avg_edge_len,2) == 0.01:
                resolution = 36
            
            elif round(avg_edge_len,3) == 0.007:
                resolution = 37
            elif round(avg_edge_len,3) == 0.004:
                resolution = 38
            elif round(avg_edge_len,3) == 0.002:
                resolution = 39
            elif round(avg_edge_len,3) <= 0.001:
                resolution = 40
                
        # Step 3: Construct the GeoJSON feature
        feature = {
            "type": "Feature",
            "geometry": mapping(cell_polygon),
            "properties": {
                    "isea3h": sibling_id,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "cell_area": round(cell_area,3),
                    "avg_edge_len": round(avg_edge_len,3),
                    "resolution": resolution,
                    "accuracy": accuracy
                    }
        }

        features.append(feature)
        
    # Step 4: Construct the FeatureCollection
    feature_collection = {
        "type": "FeatureCollection",
        "features": features
    }
    return  feature_collection

def siblings2geojson_cli():
    parser = argparse.ArgumentParser(description="Convert OpenEAGGR ISEA3H code to GeoJSON")
    parser.add_argument("isea3h", help="Input ISEA3H code, e.g., isea3h2geojson '07024,0'")
    args = parser.parse_args()
    geojson_data = json.dumps(siblings2geojson(args.isea3h))
    print(geojson_data)