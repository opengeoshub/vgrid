import argparse
import json
from shapely.geometry import Polygon
from shapely.wkt import loads
from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from pyproj import Geod
from tqdm import tqdm
from shapely.geometry import Polygon, box, mapping
from vgrid.utils.antimeridian import fix_polygon
import locale
current_locale = locale.getlocale()  # Get the current locale setting
locale.setlocale(locale.LC_ALL,current_locale)  # Use the system's default locale
geod = Geod(ellps="WGS84")

# Initialize the DGGS system
base_cells = [
    '00000,0', '01000,0', '02000,0', '03000,0', '04000,0', '05000,0', '06000,0', '07000,0', '08000,0', '09000,0',
    '10000,0', '11000,0', '12000,0', '13000,0', '14000,0', '15000,0', '16000,0', '17000,0', '18000,0', '19000,0'
]

max_cells = 1_000_000
isea3h_dggs = Eaggr(Model.ISEA3H)

res_accuracy_dict = {
        0: 25_503_281_086_204.43,
        1: 17_002_187_390_802.953,
        2: 5_667_395_796_934.327,
        3: 1_889_131_932_311.4424,
        4: 629_710_644_103.8047,
        5: 209_903_548_034.5921,
        6: 69_967_849_344.8546,
        7: 23_322_616_448.284866,
        8: 7_774_205_482.77106,
        9: 2_591_401_827.5809155,
        10: 863_800_609.1842003,
        11: 287_933_536.4041716,
        12: 95_977_845.45861907,
        13: 31_992_615.152873024,
        14: 10_664_205.060395785,
        15: 3_554_735.0295700384,
        16: 1_184_911.6670852362,
        17: 394_970.54625696875,
        18: 131_656.84875232293,
        19: 43_885.62568888426, 
        20: 14628.541896294753,
        21: 4_876.180632098251,
        22: 1_625.3841059227952,
        23: 541.7947019742651,
        24: 180.58879588146658,
        25: 60.196265293822194,
        26: 20.074859874562527,
        27: 6.6821818482323785,
        
        28: 2.2368320593659234,
        29: 0.7361725765001773,
        30: 0.2548289687885229,
        31: 0.0849429895961743,
        32: 0.028314329865391435,
       
        33: 0.009438109955130478,  
        34: 0.0031460366517101594,  
        35: 0.0010486788839033865,      
        36: 0.0003495596279677955, 
        37: 0.0001165198769892652,   
        38: 0.0000388399589964217,
        39: 0.0000129466529988072,      
        40: 0.0000043155509996024
    }    
    
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

def cell_to_polygon(eaggr_cell):
    cell_to_shape =  isea3h_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell, ShapeStringFormat.WKT)
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


def get_children_cells(base_cells, target_resolution):
    """
    Recursively generate DGGS cells for the desired resolution, avoiding duplicates.
    """
    current_cells = base_cells
    seen_cells = set(base_cells)  # Track already processed cells

    for res in range(target_resolution):
        next_cells = []
        for cell in tqdm(current_cells, desc=f"Generating child cells at resolution {res}", unit=" cells"):
            children = isea3h_dggs.get_dggs_cell_children(DggsCell(cell))
            for child in children:
                if child._cell_id not in seen_cells:
                    seen_cells.add(child._cell_id)  # Mark as seen
                    next_cells.append(child._cell_id)
        current_cells = next_cells
    return current_cells


def get_children_cells_within_bbox(bounding_cell, bbox, target_resolution):
    """
    Recursively generate DGGS cells within a bounding box, avoiding duplicates.
    """
    current_cells = [bounding_cell]  # Start with a list containing the single bounding cell
    seen_cells = set(current_cells)  # Track already processed cells
    bounding_cell2point = isea3h_dggs.convert_dggs_cell_to_point(DggsCell(bounding_cell))
    accuracy = bounding_cell2point._accuracy
    bounding_resolution = accuracy_res_dict.get(accuracy)

    if bounding_resolution <= target_resolution:
        for res in range(bounding_resolution, target_resolution):
            next_cells = []
            for cell in tqdm(current_cells, desc=f"Generating child cells at resolution {res}", unit=" cells"):
                # Get the child cells for the current cell
                children = isea3h_dggs.get_dggs_cell_children(DggsCell(cell))
                for child in children:
                    if child._cell_id not in seen_cells:  # Check if the child is already processed
                        child_shape = cell_to_polygon(child)
                        if child_shape.intersects(bbox):
                            seen_cells.add(child._cell_id)  # Mark as seen
                            next_cells.append(child._cell_id)
            if not next_cells:  # Break early if no cells remain
                break
            current_cells = next_cells  # Update current_cells to process the next level of children

        return current_cells
    else:
        print('Bounding box area is < 0.028 square meters. Please select a bigger bounding box')
        return None

def generate_grid(resolution):
    """
    Generate DGGS cells and convert them to GeoJSON features.
    """
    children = get_children_cells(base_cells, resolution)
    features = []
    for child in tqdm(children, desc="Processing cells", unit=" cells"):
        eaggr_cell = DggsCell(child)
        cell_polygon = cell_to_polygon(eaggr_cell)
        isea3h_id = eaggr_cell.get_cell_id()

        cell_centroid = cell_polygon.centroid
        center_lat =  round(cell_centroid.y, 7)
        center_lon = round(cell_centroid.x, 7)
        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        avg_edge_len = round(cell_perimeter / 6,3)
        if resolution == 0:
            avg_edge_len = round(cell_perimeter / 3,3) # icosahedron faces
            
        features.append({
            "type": "Feature",
            "geometry": mapping(cell_polygon),
            "properties": {
                    "isea3h": isea3h_id,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "cell_area": cell_area,
                    "avg_edge_len": avg_edge_len,
                    "resolution": resolution
                    },
        })
    
    
    return {
            "type": "FeatureCollection",
            "features": features
        }

def generate_grid_within_bbox(resolution,bbox):
    accuracy = res_accuracy_dict.get(resolution)
    # print(accuracy)
    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
    # print (bounding_box_wkt)
    shapes = isea3h_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
    
    for shape in shapes:
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea3h_dggs.get_bounding_dggs_cell(bbox_cells)
        # print("boudingcell: ", bounding_cell.get_cell_id())
        bounding_children_cells = get_children_cells_within_bbox(bounding_cell.get_cell_id(), bounding_box,resolution)
        # print (bounding_children_cells)
        if bounding_children_cells:
            features = []
            for child in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
                isea3h_cell = DggsCell(child)
                cell_polygon = cell_to_polygon(isea3h_cell)
                isea3h_id = isea3h_cell.get_cell_id()

                cell_centroid = cell_polygon.centroid
                center_lat =  round(cell_centroid.y, 7)
                center_lon = round(cell_centroid.x, 7)
                cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)
                cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
                avg_edge_len = round(cell_perimeter / 6,3)
                if resolution == 0:
                    avg_edge_len = round(cell_perimeter / 3,3) # icosahedron faces
            
                if cell_polygon.intersects(bounding_box):
                    features.append({
                        "type": "Feature",
                        "geometry": mapping(cell_polygon),
                        "properties": {
                                "isea3h": isea3h_id,
                                "center_lat": center_lat,
                                "center_lon": center_lon,
                                "cell_area": cell_area,
                                "avg_edge_len": avg_edge_len,
                                "resolution": resolution
                                },
                    })
                    
            return {
                "type": "FeatureCollection",
                "features": features
            }

def main():
    """
    Main function to parse arguments and generate the DGGS grid.
    """
    parser = argparse.ArgumentParser(description="Generate full DGGS grid at a specified resolution.")
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution [0..32] of the grid")
    # Resolution max range: [0..40]
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )

    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    if resolution < 0 or resolution > 32:
        print(f"Please select a resolution in [0..32] range and try again ")
        return
    
    if bbox == [-180, -90, 180, 90]:        
        num_cells =  20*(7**resolution)
        if num_cells > max_cells:
            print(
                f"The selected resolution will generate "
                f"{locale.format_string('%d', num_cells, grouping=True)} cells, "
                f"which exceeds the limit of {locale.format_string('%d', max_cells, grouping=True)}."
            )
            print("Please select a smaller resolution and try again.")
            return
        
        geojson = generate_grid(resolution)
        geojson_path = f"isea3h_grid_{resolution}.geojson"

        with open(geojson_path, 'w', encoding='utf-8') as f:
            json.dump(geojson, f, ensure_ascii=False, indent=4)

        print(f"GeoJSON saved as {geojson_path}")
    else:       
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        # Define the GeoJSON file path
        geojson_path = f"isea3h_grid_{resolution}_bbox.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print (f"GeoJSON saved as {geojson_path}")
        
if __name__ == "__main__":
    main()
