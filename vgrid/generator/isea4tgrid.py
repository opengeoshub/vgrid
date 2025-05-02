import argparse
import json
from shapely.geometry import Polygon
from shapely.wkt import loads
import platform

if (platform.system() == 'Windows'):
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.generator.settings import isea4t_res_accuracy_dict

from shapely.ops import unary_union
from tqdm import tqdm
from shapely.geometry import Polygon, box
from vgrid.utils.antimeridian import fix_polygon
from vgrid.generator.settings import max_cells, isea4t_base_cells, geodesic_dggs_to_feature


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
        for cell in current_cells:
            children = isea4t_dggs.get_dggs_cell_children(DggsCell(cell))
            next_cells.extend([child._cell_id for child in children])
        current_cells = next_cells
    return current_cells

def get_isea4t_children_cells_within_bbox(isea4t_dggs,bounding_cell, bbox, target_resolution):    
    current_cells = [bounding_cell]  # Start with a list containing the single bounding cell
    bounding_resolution = len(bounding_cell) - 2

    for res in range(bounding_resolution, target_resolution):
        next_cells = []
        for cell in current_cells:
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

   
def generate_grid(isea4t_dggs, resolution):
    # accuracy = isea4t_res_accuracy_dict.get(resolution)
    children = get_isea4t_children_cells(isea4t_dggs, isea4t_base_cells, resolution)
    isea4t_features = []
    for child in tqdm(children, desc="Generating ISEA4T DGGS", unit=" cells"):
        isea4t_cell = DggsCell(child)
        cell_polygon = isea4t_cell_to_polygon(isea4t_dggs, isea4t_cell)
        isea4t_id = isea4t_cell.get_cell_id()
        num_edges = 3
        if resolution == 0:
            cell_polygon = fix_polygon(cell_polygon)
        elif isea4t_id.startswith('00') or isea4t_id.startswith('09')\
            or isea4t_id.startswith('14') or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
            cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
        
        isea4t_feature = geodesic_dggs_to_feature('isea4t', isea4t_id,resolution, cell_polygon, num_edges)
        isea4t_features.append(isea4t_feature)
    
    return {
            "type": "FeatureCollection",
            "features": isea4t_features
        }

   
def generate_grid_within_bbox(isea4t_dggs, resolution,bbox):
    accuracy = isea4t_res_accuracy_dict.get(resolution)

    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
    shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
    shape = shapes[0]
    # for shape in shapes:
    bbox_cells = shape.get_shape().get_outer_ring().get_cells()
    bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)
    bounding_children = get_isea4t_children_cells_within_bbox(isea4t_dggs, bounding_cell.get_cell_id(), bounding_box,resolution)
    isea4t_features = []
    for child in tqdm(bounding_children, desc="Generating ISEA4T DGGS", unit=" cells"):
        isea4t_cell = DggsCell(child)
        cell_polygon = isea4t_cell_to_polygon(isea4t_dggs,isea4t_cell)
        isea4t_id = isea4t_cell.get_cell_id()
        if resolution == 0:
            cell_polygon = fix_polygon(cell_polygon)
        
        elif isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14') or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
            cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
        num_edges = 3
                
        # if cell_polygon.intersects(bounding_box):
        isea4t_feature = geodesic_dggs_to_feature('isea4t', isea4t_id, resolution, cell_polygon, num_edges)
        isea4t_features.append(isea4t_feature)           
            
    return {
        "type": "FeatureCollection",
        "features": isea4t_features
    }

def generate_grid_resample(isea4t_dggs, resolution, geojson_features):
    accuracy = isea4t_res_accuracy_dict.get(resolution)
    # Step 1: Unify all geometries into a single shape
    geometries = [shape(feature["geometry"]) for feature in geojson_features["features"]]
    unified_geom = unary_union(geometries)
    unified_geom_wkt = unified_geom.wkt

    # Step 2: Generate DGGS shapes from WKT geometry
    shapes = isea4t_dggs.convert_shape_string_to_dggs_shapes(unified_geom_wkt, ShapeStringFormat.WKT, accuracy)
    shape = shapes[0]
    bbox_cells = shape.get_shape().get_outer_ring().get_cells()
    bounding_cell = isea4t_dggs.get_bounding_dggs_cell(bbox_cells)

    # Step 3: Generate children cells within geometry bounds
    bounding_children = get_isea4t_children_cells_within_bbox(
        isea4t_dggs, bounding_cell.get_cell_id(), unified_geom, resolution
    )

    isea4t_features = []
    for child in tqdm(bounding_children, desc="Generating ISEA4T DGGS", unit=" cells"):
        isea4t_cell = DggsCell(child)
        cell_polygon = isea4t_cell_to_polygon(isea4t_dggs, isea4t_cell)
        isea4t_id = isea4t_cell.get_cell_id()

        if resolution == 0:
            cell_polygon = fix_polygon(cell_polygon)
        elif isea4t_id.startswith(('00', '09', '14', '04', '19')):
            cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)

        num_edges = 3

        # Optional: only include cells intersecting original geometry
        if not cell_polygon.intersects(unified_geom):
            continue

        isea4t_feature = geodesic_dggs_to_feature('isea4t', isea4t_id, resolution, cell_polygon, num_edges)
        isea4t_features.append(isea4t_feature)

    return {
        "type": "FeatureCollection",
        "features": isea4t_features
    }

def main():
    parser = argparse.ArgumentParser(description="Generate Open-Eaggr ISEA4T DGGS.")
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution [0..25]")
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
        if resolution < 0 or resolution > 25:
            print(f"Please select a resolution in [0..25] range and try again ")
            return
    
        if bbox == [-180, -90, 180, 90]:        
            total_cells = 20*(4**resolution)
            print(f"Resolution {resolution} will generate {total_cells} cells ")                    
            if total_cells > max_cells:
                print(f"which exceeds the limit of {max_cells}.")
                print("Please select a smaller resolution and try again.")
                return
            
            geojson_features = generate_grid(isea4t_dggs,resolution)          
        
        else:           
            geojson_features = generate_grid_within_bbox(isea4t_dggs,resolution, bbox)
       
        # Define the GeoJSON file path
        geojson_path = f"isea4t_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print (f"GeoJSON saved as {geojson_path}")

if __name__ == "__main__":
    main()
