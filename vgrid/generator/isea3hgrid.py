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
from vgrid.generator.settings import max_cells, isea3h_accuracy_res_dict, isea3h_res_accuracy_dict

geod = Geod(ellps="WGS84")

# Initialize the DGGS system
base_cells = [
    '00000,0', '01000,0', '02000,0', '03000,0', '04000,0', '05000,0', '06000,0', '07000,0', '08000,0', '09000,0',
    '10000,0', '11000,0', '12000,0', '13000,0', '14000,0', '15000,0', '16000,0', '17000,0', '18000,0', '19000,0'
]


def isea3h_cell_to_polygon(isea3h_dggs,isea3h_cell):
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


def get_isea3h_children_cells(isea3h_dggs, base_cells, target_resolution):
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


def get_isea3h_children_cells_within_bbox(isea3h_dggs,bounding_cell, bbox, target_resolution):
    """
    Recursively generate DGGS cells within a bounding box, avoiding duplicates.
    """
    current_cells = [bounding_cell]  # Start with a list containing the single bounding cell
    seen_cells = set(current_cells)  # Track already processed cells
    bounding_cell2point = isea3h_dggs.convert_dggs_cell_to_point(DggsCell(bounding_cell))
    accuracy = bounding_cell2point._accuracy
    bounding_resolution = isea3h_accuracy_res_dict.get(accuracy)

    if bounding_resolution <= target_resolution:
        for res in range(bounding_resolution, target_resolution):
            next_cells = []
            for cell in tqdm(current_cells, desc=f"Generating child cells at resolution {res}", unit=" cells"):
                # Get the child cells for the current cell
                children = isea3h_dggs.get_dggs_cell_children(DggsCell(cell))
                for child in children:
                    if child._cell_id not in seen_cells:  # Check if the child is already processed
                        child_shape = isea3h_cell_to_polygon(isea3h_dggs,child)
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

def generate_grid(isea3h_dggs, resolution):
    """
    Generate DGGS cells and convert them to GeoJSON features.
    """
    children = get_isea3h_children_cells(isea3h_dggs,base_cells, resolution)
    features = []
    for child in tqdm(children, desc="Processing cells", unit=" cells"):
        isea3h_cell = DggsCell(child)
        cell_polygon = isea3h_cell_to_polygon(isea3h_dggs,isea3h_cell)
        isea3h_id = isea3h_cell.get_cell_id()

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

def generate_grid_within_bbox(isea3h_dggs, resolution,bbox):
    accuracy = isea3h_res_accuracy_dict.get(resolution)
    # print(accuracy)
    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
    # print (bounding_box_wkt)
    shapes = isea3h_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
    
    for shape in shapes:
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = isea3h_dggs.get_bounding_dggs_cell(bbox_cells)
        # print("boudingcell: ", bounding_cell.get_cell_id())
        bounding_children_cells = get_isea3h_children_cells_within_bbox(isea3h_dggs, bounding_cell.get_cell_id(), bounding_box,resolution)
        # print (bounding_children_cells)
        if bounding_children_cells:
            features = []
            for child in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
                isea3h_cell = DggsCell(child)
                cell_polygon = isea3h_cell_to_polygon(isea3h_dggs, isea3h_cell)
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
    parser = argparse.ArgumentParser(description="Generate Open-Eaggr ISEA3H grid.")
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution [0..32] of the grid")
    # Resolution max range: [0..40]
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )
    if (platform.system() == 'Windows'):
        isea3h_dggs = Eaggr(Model.ISEA3H)
        args = parser.parse_args()
        resolution = args.resolution
        bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
        if resolution < 0 or resolution > 32:
            print(f"Please select a resolution in [0..32] range and try again ")
            return
        
        if bbox == [-180, -90, 180, 90]:        
            total_cells =  20*(7**resolution)
            print(f"Resolution {resolution} within bounding box {bbox} will generate {total_cells} cells ")
            
            if total_cells > max_cells:
                print(f"which exceeds the limit of {max_cells}. ")
                print("Please select a smaller resolution and try again.")
                return   
            
            geojson = generate_grid(isea3h_dggs,resolution)
            geojson_path = f"isea3h_grid_{resolution}.geojson"

            with open(geojson_path, 'w', encoding='utf-8') as f:
                json.dump(geojson, f, ensure_ascii=False, indent=4)

            print(f"GeoJSON saved as {geojson_path}")
        else:       
            # Generate grid within the bounding box
            geojson_features = generate_grid_within_bbox(isea3h_dggs,resolution, bbox)
            # Define the GeoJSON file path
            geojson_path = f"isea3h_grid_{resolution}_bbox.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print (f"GeoJSON saved as {geojson_path}")

if __name__ == "__main__":
    main()
