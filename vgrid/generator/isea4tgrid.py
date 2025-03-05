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
import platform

geod = Geod(ellps="WGS84")

# Initialize the DGGS system
base_cells = [
    '00', '01', '02', '03', '04', '05', '06', '07', '08', '09',
    '10', '11', '12', '13', '14', '15', '16', '17', '18', '19'
]
max_cells = 1_000_000

def fix_isea4t_wkt(isea4t_wkt):
    # Extract the coordinate section
    coords_section = isea4t_wkt[isea4t_wkt.index("((") + 2 : isea4t_wkt.index("))")]
    coords = coords_section.split(",")
    # Append the first point to the end if not already closed
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    fixed_coords = ", ".join(coords)
    return f"POLYGON (({fixed_coords}))"

def fix_isea4t_antimeridian_cells(isea4t_boundary, threshold=-100):
    """
    Adjusts polygon coordinates to handle antimeridian crossings.
    """
    lon_lat = [(float(lon), float(lat)) for lon, lat in isea4t_boundary.exterior.coords]

    if any(lon < threshold for lon, _ in lon_lat):
        adjusted_coords = [(lon - 360 if lon > 0 else lon, lat) for lon, lat in lon_lat]
    else:
        adjusted_coords = lon_lat

    return Polygon(adjusted_coords)

def isea4t_cell_to_polygon(isea4t_dggs,isea4t_cell):
    cell_to_shp =  isea4t_dggs.convert_dggs_cell_outline_to_shape_string(isea4t_cell, ShapeStringFormat.WKT)
    cell_to_shp_fixed = fix_isea4t_wkt(cell_to_shp)
    cell_polygon = loads(cell_to_shp_fixed)
    return cell_polygon


def get_isea4t_children_cells(isea4t_dggs,base_cells, target_resolution):
    """
    Recursively generate DGGS cells for the desired resolution.
    """
    current_cells = base_cells
    for res in range(target_resolution):
        next_cells = []
        for cell in tqdm(current_cells, desc= f"Generating child cells at resolution {res}", unit=" cells"):
            children = isea4t_dggs.get_dggs_cell_children(DggsCell(cell))
            next_cells.extend([child._cell_id for child in children])
        current_cells = next_cells
    return current_cells

def get_isea4t_children_cells_within_bbox(isea4t_dggs,bounding_cell, bbox, target_resolution):
    """
    Recursively generate child cells for the given bounding cell up to the target resolution,
    considering only cells that intersect with the given bounding box.

    Parameters:
        bounding_cell (str): The starting cell ID.
        bbox (Polygon): The bounding box as a Shapely Polygon.
        target_resolution (int): The target resolution for cell generation.

    Returns:
        list: List of cell IDs that intersect with the bounding box.
    """
    current_cells = [bounding_cell]  # Start with a list containing the single bounding cell
    bounding_resolution = len(bounding_cell) - 2

    for res in range(bounding_resolution, target_resolution):
        next_cells = []
        for cell in tqdm(current_cells, desc=f"Generating child cells at resolution {res}", unit=" cells"):
            # Get the child cells for the current cell
            children = isea4t_dggs.get_dggs_cell_children(DggsCell(cell))
            for child in children:
                # Convert child cell to geometry
                child_shape = isea4t_cell_to_polygon(isea4t_dggs, child)
                if child_shape.intersects(bbox):              
                    # Add the child cell ID to the next_cells list
                    next_cells.append(child._cell_id)  
        if not next_cells:  # Break early if no cells remain
            break
        current_cells = next_cells  # Update current_cells to process the next level of children
    
    return current_cells

# length_accuracy_dict = {
#     41: 10**-10,
#     40: 5*10**-10,
#     39: 10**-9,
#     38: 10**-8,
#     37: 5*10**-8,
#     36: 10**-7,
#     35: 5*10**-7,
#     34: 10**-6,
#     33: 5*10**-6,
#     32: 5*10**-5,
#     31: 10**-4,
#     30: 5*10**-4,
#     29: 9*10**-4,
#     28: 5*10**-3,
#     27: 2*10**-2,
#     26: 5*10**-2,
#     25: 5*10**-1,
#     24: 1,
#     23: 10,
#     22: 5*10,
#     21: 10**2,
#     20: 5*10**2,
#     19: 10**3,
#     18: 5*10**3,
#     17: 5*10**4,
#     16: 10**5,
#     15: 5*10**5,
#     14: 10**6,
#     13: 5*10**6,
#     12: 5*10**7,
#     11: 10**8,
#     10: 5*10**8,
#      9: 10**9,
#      8: 10**10,
#      7: 5*10**10,
#      6: 10**11,
#      5: 5*10**11,
#      4: 10**12,
#      3: 5*10**12,
#      2: 5*10**13
# }

isea4t_res_accuracy_dict = {
    0: 25_503_281_086_204.43,
    1: 6_375_820_271_551.114,    
    2: 1_593_955_067_887.7715,
    3: 398_488_766_971.94995,
    4: 99_622_191_742.98041,
    5: 24905_547_935.752182,
    6: 6_226_386_983.930966,
    7: 1_556_596_745.9898202,
    8: 389_149_186.4903765,
    9: 97_287_296.6296727,
    10: 24_321_824.150339592,
    11: 6_080_456.0446634805,
    12: 1_520_114.0040872877,
    13: 380_028.5081004044,
    14: 95_007.11994651864,
    15: 23_751.787065212124,
    16: 5_937.9396877205645,
    17: 1_484.492000512607,
    18: 371.1159215456855,
    19: 92.78605896888773,    
    20: 23.189436159755584,
    21: 5.804437622405244,
    22: 1.4440308231349632,
    23: 0.36808628825008866,
    24: 0.0849429895961743,
    25: 0.028314329865391435,
    
    26: 7.08*10**-3, # accuracy returns 0.0, avg_edge_len =  0.11562
    27: 1.77*10**-3, # accuracy returns 0.0, avg_edge_len =  0.05781
    28: 4.42*10**-4, # accuracy returns 0.0, avg_edge_len =  0.0289
    29: 1.11*10**-4, # accuracy returns 0.0, avg_edge_len =  0.01445
    30: 2.77*10**-5, # accuracy returns 0.0, avg_edge_len = 0.00723
    31: 6.91*10**-6, # accuracy returns 0.0, avg_edge_len =  0.00361
    32: 1.73*10**-6, # accuracy returns 0.0, avg_edge_len =  0.00181
    33: 5.76*10**-7, # accuracy returns 0.0, avg_edge_len = 0.0009
    34: 1.92*10**-7, # accuracy returns 0.0, avg_edge_len = 0.00045
    35: 6.40*10**-8, # accuracy returns 0.0, avg_edge_len = 0.00023
    36: 2.13*10**-8, # accuracy returns 0, avg_edge_len = 0.00011
    37: 7.11*10**-9, # accuracy returns 0.0, avg_edge_len = 6*10**(-5)
    38: 2.37*10**-9, # accuracy returns 0.0, avg_edge_len = 3*10**(-5)
    39: 7.90*10**-10 # accuracy returns 0.0, avg_edge_len = 10**(-5)
    }       

def generate_grid(isea4t_dggs, resolution):
    """
    Generate DGGS cells and convert them to GeoJSON features.
    """
    accuracy = isea4t_res_accuracy_dict.get(resolution)
    children = get_isea4t_children_cells(isea4t_dggs, base_cells, resolution)
    features = []
    for child in tqdm(children, desc="Processing cells", unit=" cells"):
        isea4t_cell = DggsCell(child)
        cell_polygon = isea4t_cell_to_polygon(isea4t_dggs, isea4t_cell)
        isea4t_cell_id = isea4t_cell.get_cell_id()

        if resolution == 0:
            cell_polygon = fix_polygon(cell_polygon)
        elif isea4t_cell_id.startswith('00') or isea4t_cell_id.startswith('09')\
            or isea4t_cell_id.startswith('14') or isea4t_cell_id.startswith('04') or isea4t_cell_id.startswith('19'):
            cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
        
        cell_centroid = cell_polygon.centroid
        center_lat =  round(cell_centroid.y, 7)
        center_lon = round(cell_centroid.x, 7)
        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),5)
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        avg_edge_len = round(cell_perimeter / 3,5)
        
        features.append({
            "type": "Feature",
            "geometry": mapping(cell_polygon),
            "properties": {
                    "eaggr_isea4t": isea4t_cell_id,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "cell_area": cell_area,
                    "avg_edge_len": avg_edge_len,
                    "resolution": resolution,
                    "accuracy": accuracy
                    },
        })
    
    
    return {
            "type": "FeatureCollection",
            "features": features
        }

   
def generate_grid_within_bbox(isea4t_dggs, resolution,bbox):
    accuracy = isea4t_res_accuracy_dict.get(resolution)

    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
    shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
    for shape in shapes:
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_children_cells = get_isea4t_children_cells_within_bbox(isea4t_dggs, bounding_cell.get_cell_id(), bounding_box,resolution)
        features = []
        for child in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
            isea4t_cell = DggsCell(child)
            cell_polygon = isea4t_cell_to_polygon(isea4t_dggs,isea4t_cell)
            isea4t_cell_id = isea4t_cell.get_cell_id()
            if resolution == 0:
                cell_polygon = fix_polygon(cell_polygon)
            
            elif isea4t_cell_id.startswith('00') or isea4t_cell_id.startswith('09') or isea4t_cell_id.startswith('14') or isea4t_cell_id.startswith('04') or isea4t_cell_id.startswith('19'):
                cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
            
            cell_centroid = cell_polygon.centroid
            center_lat =  round(cell_centroid.y, 7)
            center_lon = round(cell_centroid.x, 7)
            cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)
            cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
            avg_edge_len = round(cell_perimeter / 3,3)
            
            if cell_polygon.intersects(bounding_box):
                features.append({
                    "type": "Feature",
                    "geometry": mapping(cell_polygon),
                    "properties": {
                            "isea4t": isea4t_cell_id,
                            "center_lat": center_lat,
                            "center_lon": center_lon,
                            "cell_area": cell_area,
                            "avg_edge_len": avg_edge_len,
                            "resolution": resolution,
                            "accuracy": accuracy
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
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution [0..39] of the grid")
    # Resolution max range: [0..39]
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )
    
    if (platform.system() == 'Windows'):
        isea4t_dggs = Eaggr(Model.ISEA4T)
        args = parser.parse_args()
        resolution = args.resolution
        bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
        if bbox == [-180, -90, 180, 90]:        
            total_cells = 20*(4**resolution)
            print(f"The selected resolution will generate {total_cells} cells ")                    
            if total_cells > max_cells:
                print(f"which exceeds the limit of {max_cells}.")
                print("Please select a smaller resolution and try again.")
                return
            
            geojson = generate_grid(isea4t_dggs,resolution)
            geojson_path = f"isea4t_grid_{resolution}.geojson"

            with open(geojson_path, 'w', encoding='utf-8') as f:
                json.dump(geojson, f, ensure_ascii=False, indent=4)

            print(f"GeoJSON saved as {geojson_path}")
        else:
            if resolution < 1 or resolution > 25:
                print(f"Please select a resolution in [1..25] range and try again ")
                return
            # Generate grid within the bounding box
            geojson_features = generate_grid_within_bbox(isea4t_dggs,resolution, bbox)
            # Define the GeoJSON file path
            geojson_path = f"isea4t_grid_{resolution}_bbox.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print (f"GeoJSON saved as {geojson_path}")

if __name__ == "__main__":
    main()
