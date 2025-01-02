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
    0: 25503281086204.43,
    1: 629710644103.8047,
    2: 69967849344.8546,
    3: 7774205482.77106,
    4: 863800609.1842003,
    5: 95977845.45861907,
    # 6: 10664205.060395785,
    6: 31992615.152873024,
    7: 131656.84875232293,
    8: 43885.62568888426,
    9: 14628.541896294753,
    10: 541.7947019742651,
    11: 60.196265293822194,
    12: 6.6821818482323785,
    13: 0.7361725765001773,
    14: 0.0849429895961743,
    15: 0.0
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
    Recursively generate DGGS cells for the desired resolution.
    """
    current_cells = base_cells
    for res in range(target_resolution):
        next_cells = []
        for cell in tqdm(current_cells, desc= f"Generating child cells at resolution {res}", unit=" cells"):
            children = isea3h_dggs.get_dggs_cell_children(DggsCell(cell))
            next_cells.extend([child._cell_id for child in children])
        current_cells = next_cells
    return current_cells

def get_resolution_from_accuracy(accuracy, dictionary):
    """Get resolution from accuracy level."""
    # Invert the dictionary temporarily for value lookup
    inverted_dict = {value: key for key, value in dictionary.items()}
    return inverted_dict.get(accuracy, "Resolution not found")


def get_children_cells_within_bbox(bounding_cell, bbox, target_resolution):
    current_cells = [bounding_cell]  # Start with a list containing the single bounding cell
    bounding_cell2point = isea3h_dggs.convert_dggs_cell_to_point(DggsCell(bounding_cell))    
    accuracy = bounding_cell2point._accuracy
    bounding_resolution = get_resolution_from_accuracy(accuracy,res_accuracy_dict)

    for res in range(bounding_resolution, target_resolution):
        next_cells = []
        for cell in tqdm(current_cells, desc=f"Generating child cells at resolution {res}", unit=" cells"):
            # Get the child cells for the current cell
            print(cell)
            children = isea3h_dggs.get_dggs_cell_children(DggsCell(cell))
            for child in children:
                child_shape = cell_to_polygon(child)
                if child_shape.intersects(bbox):              
                    next_cells.append(child._cell_id)  # Use append instead of extend
        if not next_cells:  # Break early if no cells remain
            break
        current_cells = next_cells  # Update current_cells to process the next level of children
    
    return current_cells


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
        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        avg_edge_len = round(cell_perimeter / 6,2)
        
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
    print(accuracy)
    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
    print (bounding_box_wkt)
    shapes = isea3h_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
    
    for shape in shapes:
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea3h_dggs.get_bounding_dggs_cell(bbox_cells)
        print("boudingcell: ", bounding_cell.get_cell_id())
        bounding_children_cells = get_children_cells_within_bbox(bounding_cell.get_cell_id(), bounding_box,resolution)
        # print (bounding_children_cells)
        features = []
        for child in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
            isea3h_cell = DggsCell(child)
            cell_polygon = cell_to_polygon(isea3h_cell)
            isea3h_id = isea3h_cell.get_cell_id()

            cell_centroid = cell_polygon.centroid
            center_lat =  round(cell_centroid.y, 7)
            center_lon = round(cell_centroid.x, 7)
            cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)
            cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
            avg_edge_len = round(cell_perimeter / 6,2)
            
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
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution [0..13] of the grid")
    # Resolution max range: [0..18]
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )

    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    
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
        if resolution < 0 or resolution > 13:
            print(f"Please select a resolution in [0..13] range and try again ")
            return
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        # Define the GeoJSON file path
        geojson_path = f"isea3h_grid_{resolution}_bbox.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print (f"GeoJSON saved as {geojson_path}")
        
if __name__ == "__main__":
    main()
