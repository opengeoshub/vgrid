import argparse 
import json
from shapely.geometry import mapping, Point, Polygon, box
from tqdm import tqdm
from pyproj import Geod
from vgrid.utils.easedggs.constants import grid_spec, ease_crs, geo_crs, levels_specs
from vgrid.utils.easedggs.dggs.grid_addressing import grid_ids_to_geos, geo_polygon_to_grid_ids

# Initialize the geodetic model
geod = Geod(ellps="WGS84")

max_cells = 1_000_000
chunk_size=10_000

geo_bounds = grid_spec['geo']
min_longitude = geo_bounds['min_x']
min_lattitude = geo_bounds['min_y']
max_longitude = geo_bounds['max_x']
max_latitude = geo_bounds['max_y']

def get_cells(resolution):
    """
    Generate a list of cell IDs based on the resolution, row, and column.
    """
    n_row = levels_specs[resolution]["n_row"]
    n_col = levels_specs[resolution]["n_col"]

    # Generate list of cell IDs
    cell_ids = []

    # Loop through all rows and columns at the specified resolution
    for row in range(n_row):
        for col in range(n_col):
            # Generate base ID (e.g., L0.RRRCCC for res=0)
            base_id = f"L{resolution}.{row:03d}{col:03d}"

            # Add additional ".RC" for each higher resolution
            cell_id = base_id
            for i in range(1, resolution + 1):
                cell_id += f".{row:1d}{col:1d}"  # For res=1: L0.RRRCCC.RC, res=2: L0.RRRCCC.RC.RC, etc.

            # Append the generated cell ID to the list
            cell_ids.append(cell_id)

    return cell_ids


def get_cells_bbox(resolution, bbox):
    bounding_box = box(*bbox)
    bounding_box_wkt = bounding_box.wkt
    cells_bbox = geo_polygon_to_grid_ids(bounding_box_wkt, level=resolution, source_crs = geo_crs, target_crs = ease_crs, levels_specs = levels_specs, return_centroids = True, wkt_geom=True)
    return cells_bbox

def generate_grid(resolution):
    features = []
    min_lon, min_lat, max_lon, max_lat = [min_longitude, min_lattitude, max_longitude, max_latitude]
    
    level_spec = levels_specs[resolution]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]
    x_length = level_spec["x_length"]
    y_length = level_spec["y_length"]

    cells = get_cells(resolution)

    # Process cells in chunks with tqdm progress bar
    for i in tqdm(range(0, len(cells), chunk_size), total=(len(cells) // chunk_size) + 1, desc="Processing cells in chunks", unit="chunk"):
        chunk = cells[i:i + chunk_size]
        for cell in chunk:
            geo = grid_ids_to_geos([cell])
            center_lon, center_lat = geo['result']['data'][0]
            cell_min_lat = center_lat - (180 / (2 * n_row))
            cell_max_lat = center_lat + (180 / (2 * n_row))
            cell_min_lon = center_lon - (360 / (2 * n_col))
            cell_max_lon = center_lon + (360 / (2 * n_col))

        #    # Calculate the cell indices
        #     row = int((max_lat - center_lat) / (180 / n_row))
        #     col = int((center_lon - min_lon) / (360 / n_col))

        #     # Validate row and col within bounds
        #     row = max(0, min(row, n_row - 1))
        #     col = max(0, min(col, n_col - 1))

        #     # Optional: Create a GeoJSON for the cell (bounding box)
        #     cell_min_lat = max_lat - (row + 1) * (180 / n_row)
        #     cell_max_lat = max_lat - row * (180 / n_row)
        #     cell_min_lon = min_lon + col * (360 / n_col)
        #     cell_max_lon = min_lon + (col + 1) * (360 / n_col)

            cell_polygon = Polygon([
                [cell_min_lon, cell_min_lat],
                [cell_max_lon, cell_min_lat],
                [cell_max_lon, cell_max_lat],
                [cell_min_lon, cell_max_lat],
                [cell_min_lon, cell_min_lat]
            ])

            cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
            # Calculate width (longitude difference at a constant latitude)
            cell_width = round(geod.line_length([cell_min_lon, cell_max_lon], [cell_min_lat, cell_max_lat]),3)
                
            # Calculate height (latitude difference at a constant longitude)
            cell_height = round(geod.line_length([cell_min_lon, cell_min_lon], [cell_min_lat, cell_max_lat]),3)

            features.append({
                "type": "Feature",
                "geometry": mapping(cell_polygon),
                "properties": {
                    "ease": cell,
                    "center_lat": round(center_lat, 7),
                    "center_lon": round(center_lon, 7),
                    "cell_area": round(cell_area,3),
                    "cell_width": cell_width,
                    "cell_height": cell_height,
                    "resolution": resolution,
                }
            })

    # Create GeoJSON FeatureCollection
    geojson_features = {
        "type": "FeatureCollection",
        "features": features
    }

    return geojson_features

def generate_grid_bbox_point(resolution, bbox):
    features = []  
    # Get all grid cells within the bounding box
    cells = get_cells_bbox(resolution, bbox)['result']['data']   
   
    if cells:
        # Use tqdm for progress bar, processing cells sequentially
        for cell in tqdm(cells, desc="Processing cells", unit=" cells"):
            geo = grid_ids_to_geos([cell])
            if geo:
                center_lon, center_lat = geo['result']['data'][0]
                cell_point = Point(center_lon, center_lat)           
                features.append({
                "type": "Feature",
                "geometry": mapping(cell_point),
                "properties": {
                    "ease": cell,
                    "center_lat": round(center_lat, 7),
                    "center_lon": round(center_lon, 7),               
                    "resolution": resolution,
                    }
                })

        # Create GeoJSON FeatureCollection
        geojson_features = {
            "type": "FeatureCollection",
            "features": features
        }

        return geojson_features

def generate_grid_bbox(resolution, bbox):
    features = []
    level_spec = levels_specs[resolution]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]

    # Get all grid cells within the bounding box
    cells = get_cells_bbox(resolution, bbox)['result']['data']   
   
    if cells:
        # Use tqdm for progress bar, processing cells sequentially
        for cell in tqdm(cells, desc="Processing cells", unit=" cells"):
            geo = grid_ids_to_geos([cell])
            if geo:
                center_lon, center_lat = geo['result']['data'][0]            
                cell_min_lat = center_lat - (180 / (2 * n_row))
                cell_max_lat = center_lat + (180 / (2 * n_row))
                cell_min_lon = center_lon - (360 / (2 * n_col))
                cell_max_lon = center_lon + (360 / (2 * n_col))

                cell_polygon = Polygon([
                    [cell_min_lon, cell_min_lat],
                    [cell_max_lon, cell_min_lat],
                    [cell_max_lon, cell_max_lat],
                    [cell_min_lon, cell_max_lat],
                    [cell_min_lon, cell_min_lat]
                ])

                cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
                # Calculate width (longitude difference at a constant latitude)
                cell_width = round(geod.line_length([cell_min_lon, cell_max_lon], [cell_min_lat, cell_max_lat]),3)
                
                # Calculate height (latitude difference at a constant longitude)
                cell_height = round(geod.line_length([cell_min_lon, cell_min_lon], [cell_min_lat, cell_max_lat]),3)


                # Append feature to list
                features.append({
                    "type": "Feature",
                    "geometry": mapping(cell_polygon),
                    "properties": {
                        "ease": cell,
                        "center_lat": round(center_lat, 7),
                        "center_lon": round(center_lon, 7),
                        "cell_area": round(cell_area,3),
                        "cell_width": cell_width,
                        "cell_height": cell_height,
                        "resolution": resolution,
                    }
                })

        # Create GeoJSON FeatureCollection
        geojson_features = {
            "type": "FeatureCollection",
            "features": features
        }

        return geojson_features


def main():
    parser = argparse.ArgumentParser(description='Generate EASE-DGGS grid.')
    parser.add_argument('-r', '--resolution', type=int, required=True, help='zoom level/ resolution= [0..6]')
    parser.add_argument('-b', '--bbox', type=float, nargs=4, help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)")

    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [min_longitude, min_lattitude, max_longitude, max_latitude]

    if resolution < 0 or resolution > 6:
        print(f"Please select a resolution in [0..6] range and try again")
        return

    if bbox == [min_longitude, min_lattitude, max_longitude, max_latitude]:
        level_spec = levels_specs[resolution]        
        n_row = level_spec["n_row"]
        n_col = level_spec["n_col"]
        total_cells = n_row*n_col
        
        print(f"Resolution {resolution} will generate {total_cells} cells ")
        if total_cells > max_cells:
            print(f"which exceeds the limit of {max_cells}.")
            print("Please select a smaller resolution and try again.")
            return
        # Start generating and saving the grid in chunks
        geojson_features = generate_grid(resolution)
        
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = f"ease_grid_{resolution}.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")
    else: 
        # total_cells = len(get_cells_bbox(resolution,bbox))        
        # print(f"Resolution {resolution} will generate {total_cells} cells ")
        # if total_cells > max_cells:
        #     print(f"which exceeds the limit of {max_cells}.")
        #     print("Please select a smaller resolution and try again.")
        #     return                       

        # Start generating and saving the grid in chunks
        geojson_features = generate_grid_bbox_point(resolution, bbox)
        
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = f"ease_grid_bbox_{resolution}.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")



if __name__ == '__main__':
    main()
