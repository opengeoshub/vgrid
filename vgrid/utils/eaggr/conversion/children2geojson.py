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
            25_503_281_086_204.43: 0,
            17_002_187_390_802.953: 1,
            5_667_395_796_934.327: 2,
            1_889_131_932_311.4424: 3,
            629_710_644_103.8047: 4,
            209_903_548_034.5921: 5,
            69_967_849_344.8546: 6,
            23_322_616_448.284866: 7,
            7_774_205_482.77106: 8,
            2_591_401_827.5809155: 9,
            863_800_609.1842003: 10,
            287_933_536.4041716: 11,
            95_977_845.45861907: 12,
            31_992_615.152873024: 13,
            10_664_205.060395785: 14,
            3_554_735.0295700384: 15,
            1_184_911.6670852362: 16,
            394_970.54625696875: 17,
            131_656.84875232293: 18,
            43_885.62568888426: 19,
            14628.541896294753: 20,
            4_876.180632098251: 21,
            1_625.3841059227952: 22,
            541.7947019742651: 23,
            180.58879588146658: 24,
            60.196265293822194: 25,
            20.074859874562527: 26,
            6.6821818482323785: 27,
            
            2.2368320593659234: 28,
            0.7361725765001773: 29,
            0.2548289687885229: 30,
            0.0849429895961743: 31,
            0.028314329865391435: 32,
            
            0.0: 33, # isea3h2point._accuracy always returns 0.0 from res 33
            0.0: 34,
            0.0: 35,
            0.0: 36,
            0.0: 37,
            0.0: 38,
            0.0: 39,
            0.0: 40
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
                    "isea3h": child_id,
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

def children2geojson_cli():
    """
    Command-line interface for isea3h2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert OpenEAGGR ISEA3H code to GeoJSON")
    parser.add_argument("isea3h", help="Input ISEA3H code, e.g., isea3h2geojson '07024,0'")
    args = parser.parse_args()
    geojson_data = json.dumps(children2geojson(args.isea3h))
    print(geojson_data)