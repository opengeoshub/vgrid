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
import locale
current_locale = locale.getlocale()  # Get the current locale setting
locale.setlocale(locale.LC_ALL,current_locale)  # Use the system's default locale
geod = Geod(ellps="WGS84")

# Initialize the DGGS system
base_cells = [
    '00', '01', '02', '03', '04', '05', '06', '07', '08', '09',
    '10', '11', '12', '13', '14', '15', '16', '17', '18', '19'
]
eaggr_dggs = Eaggr(Model.ISEA4T)

def fix_isea4t_wkt(eaggr_wkt):
    # Extract the coordinate section
    coords_section = eaggr_wkt[eaggr_wkt.index("((") + 2 : eaggr_wkt.index("))")]
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

def cell_to_polygon(eaggr_cell):
    cell_to_shp =  eaggr_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell, ShapeStringFormat.WKT)
    cell_to_shp_fixed = fix_isea4t_wkt(cell_to_shp)
    cell_polygon = loads(cell_to_shp_fixed)
    return Polygon(cell_polygon)

def get_children_cells(base_cells, target_resolution):
    """
    Recursively generate DGGS cells for the desired resolution.
    """
    current_cells = base_cells
    for res in range(target_resolution):
        next_cells = []
        for cell in tqdm(current_cells, desc= f"Generating child cells at resolution {res}", unit="cell"):
            children = eaggr_dggs.get_dggs_cell_children(DggsCell(cell))
            next_cells.extend([child._cell_id for child in children])
        current_cells = next_cells
    return current_cells

def get_children_cells_within_bbox(bounding_cell, bbox, target_resolution):
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
        for cell in tqdm(current_cells, desc=f"Generating child cells at resolution {res}", unit="cell"):
            # Get the child cells for the current cell
            children = eaggr_dggs.get_dggs_cell_children(DggsCell(cell))
            for child in children:
                # Convert child cell to geometry
                child_shape = cell_to_polygon(child)
                if child_shape.intersects(bbox):              
                    # Add the child cell ID to the next_cells list
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
        eaggr_cell_id = eaggr_cell.get_cell_id()

        if eaggr_cell_id.startswith('00') or eaggr_cell_id.startswith('09') or eaggr_cell_id.startswith('14') or eaggr_cell_id.startswith('04') or eaggr_cell_id.startswith('19'):
            cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
        
        # cell_centroid = cell_polygon.centroid
        # center_lat =  round(cell_centroid.y, 7)
        # center_lon = round(cell_centroid.x, 7)
        # cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)
        # cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        # avg_edge_len = round(cell_perimeter / 3,2)
        
        features.append({
            "type": "Feature",
            "geometry": mapping(cell_polygon),
            "properties": {
                    "isea4t": eaggr_cell_id,
                    # "center_lat": center_lat,
                    # "center_lon": center_lon,
                    # "cell_area": cell_area,
                    # "avg_edge_len": avg_edge_len,
                    # "resolution": resolution
                    },
        })
    
    
    return {
            "type": "FeatureCollection",
            "features": features
        }

length_accuracy_dict = {
    41: 10**-10,
    40: 5*10**-10,
    39: 10**-9,
    38: 10**-8,
    37: 5*10**-8,
    36: 10**-7,
    35: 5*10**-7,
    34: 10**-6,
    33: 5*10**-6,
    32: 5*10**-5,
    31: 10**-4,
    30: 5*10**-4,
    29: 9*10**-4,
    28: 5*10**-3,
    27: 2*10**-2,
    26: 5*10**-2,
    25: 5*10**-1,
    24: 1,
    23: 10,
    22: 5*10,
    21: 10**2,
    20: 5*10**2,
    19: 10**3,
    18: 5*10**3,
    17: 5*10**4,
    16: 10**5,
    15: 5*10**5,
    14: 10**6,
    13: 5*10**6,
    12: 5*10**7,
    11: 10**8,
    10: 5*10**8,
     9: 10**9,
     8: 10**10,
     7: 5*10**10,
     6: 10**11,
     5: 5*10**11,
     4: 10**12,
     3: 5*10**12,
     2: 5*10**13
}

def generate_grid_within_bbox(resolution,bbox):
    cell_id_len = resolution +2
    accuracy = length_accuracy_dict.get(cell_id_len)

    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt  # Create a bounding box polygon
    shapes = eaggr_dggs.convert_shape_string_to_dggs_shapes(bounding_box_wkt, ShapeStringFormat.WKT, accuracy)
    for shape in shapes:
        bbox_cells = shape.get_shape().get_outer_ring().get_cells()
        bounding_cell = eaggr_dggs.get_bounding_dggs_cell(bbox_cells)
        bounding_children_cells = get_children_cells_within_bbox(bounding_cell.get_cell_id(), bounding_box,resolution)
        features = []
        for child in tqdm(bounding_children_cells, desc="Processing cells", unit=" cells"):
            eaggr_cell = DggsCell(child)
            cell_polygon = cell_to_polygon(eaggr_cell)
            eaggr_cell_id = eaggr_cell.get_cell_id()

            if eaggr_cell_id.startswith('00') or eaggr_cell_id.startswith('09') or eaggr_cell_id.startswith('14') or eaggr_cell_id.startswith('04') or eaggr_cell_id.startswith('19'):
                cell_polygon = fix_isea4t_antimeridian_cells(cell_polygon)
            
            # cell_centroid = cell_polygon.centroid
            # center_lat =  round(cell_centroid.y, 7)
            # center_lon = round(cell_centroid.x, 7)
            # cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)
            # cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
            # avg_edge_len = round(cell_perimeter / 3,2)
            
            if cell_polygon.intersects(bounding_box):
                features.append({
                    "type": "Feature",
                    "geometry": mapping(cell_polygon),
                    "properties": {
                            "isea4t": eaggr_cell_id,
                            # "center_lat": center_lat,
                            # "center_lon": center_lon,
                            # "cell_area": cell_area,
                            # "avg_edge_len": avg_edge_len,
                            # "resolution": resolution
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
    parser.add_argument("-r", "--resolution", type=int, required=True, help="Resolution [0..22] of the grid")
    # Resolution max range: [0..39]
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )

    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    
    if bbox == [-180, -90, 180, 90]:        
        num_cells = 20*(4**resolution)
        max_cells = 1_000_000
        if num_cells > max_cells:
            print(
                f"The selected resolution will generate "
                f"{locale.format_string('%d', num_cells, grouping=True)} cells, "
                f"which exceeds the limit of {locale.format_string('%d', max_cells, grouping=True)}."
            )
            print("Please select a smaller resolution and try again.")
            return
        
        geojson = generate_grid(resolution)
        geojson_path = f"isea4t_grid_{resolution}.geojson"

        with open(geojson_path, 'w', encoding='utf-8') as f:
            json.dump(geojson, f, ensure_ascii=False, indent=4)

        print(f"GeoJSON saved as {geojson_path}")
    else:
        if resolution < 1 or resolution > 22:
            print(f"Please select a resolution in [1..22] range and try again ")
            return
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        # Define the GeoJSON file path
        geojson_path = f"isea4t_grid_{resolution}_bbox.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print (f"GeoJSON saved as {geojson_path}")
        
if __name__ == "__main__":
    main()
