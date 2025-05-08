# https://github.com/ha8tks/Leaflet.Maidenhead
# https://ha8tks.github.io/Leaflet.Maidenhead/examples/
# https://www.sotamaps.org/

import json, math
import argparse
from vgrid.utils import maidenhead
from tqdm import tqdm  
from vgrid.generator.settings import max_cells, graticule_dggs_to_feature
from shapely.geometry import Polygon
from vgrid.stats.maidenheadstats import maidenhead_metrics

def generate_grid(resolution):
    if resolution == 1:
        lon_width, lat_width = 20, 10
    elif resolution == 2:
        lon_width, lat_width = 2, 1
    elif resolution == 3:
        lon_width, lat_width = 0.083333, 0.041666  # 5 minutes x 2.5 minutes
    elif resolution == 4:
        lon_width, lat_width = 0.008333, 0.004167  # 30 seconds x 15 seconds
    else:
        raise ValueError("Unsupported resolution")

    # Determine the bounding box
    min_lon, min_lat, max_lon, max_lat = [-180, -90, 180, 90]
    x_cells = int((max_lon - min_lon) / lon_width)
    y_cells = int((max_lat - min_lat) / lat_width)
    total_cells = x_cells * y_cells

    maidenhead_features = []
    with tqdm(total=total_cells, desc="Generating Maidenhead DGGS", unit=" cells") as pbar:
        for i in range(x_cells):
            for j in range(y_cells):
                cell_min_lon = min_lon + i * lon_width
                cell_max_lon = cell_min_lon + lon_width
                cell_min_lat = min_lat + j * lat_width
                cell_max_lat = cell_min_lat + lat_width

                cell_center_lat = (cell_min_lat + cell_max_lat) / 2
                cell_center_lat = (cell_min_lon + cell_max_lon) / 2
                maidenhead_id = maidenhead.toMaiden(cell_center_lat, cell_center_lat, resolution)
                _, _, min_lat_maiden, min_lon_maiden, max_lat_maiden, max_lon_maiden, _ = maidenhead.maidenGrid(maidenhead_id)
                # Define the polygon based on the bounding box
                cell_polygon = Polygon([
                    [min_lon_maiden, min_lat_maiden],  # Bottom-left corner
                    [max_lon_maiden, min_lat_maiden],  # Bottom-right corner
                    [max_lon_maiden, max_lat_maiden],  # Top-right corner
                    [min_lon_maiden, max_lat_maiden],  # Top-left corner
                    [min_lon_maiden, min_lat_maiden]   # Closing the polygon (same as the first point)
                ])
                maidenhead_feature = graticule_dggs_to_feature('maidenhead',maidenhead_id,resolution,cell_polygon)
                maidenhead_features.append(maidenhead_feature)
                pbar.update(1)

    return {
        "type": "FeatureCollection",
        "features": maidenhead_features
    }

def generate_grid_within_bbox(resolution, bbox):
    # Define the grid parameters based on the resolution
    if resolution == 1:
        lon_width, lat_width = 20, 10
    elif resolution == 2:
        lon_width, lat_width = 2, 1
    elif resolution == 3:
        lon_width, lat_width = 0.083333, 0.041666  # 5 minutes x 2.5 minutes
    elif resolution == 4:
        lon_width, lat_width = 0.008333, 0.004167  # 30 seconds x 15 seconds
    else:
        raise ValueError("Unsupported resolution")

    min_lon, min_lat, max_lon, max_lat = bbox

    # Calculate grid cell indices for the bounding box
    base_lat, base_lon = -90, -180
    start_x = math.floor((min_lon - base_lon) / lon_width)
    end_x = math.floor((max_lon - base_lon) / lon_width)
    start_y = math.floor((min_lat - base_lat) / lat_width)
    end_y = math.floor((max_lat - base_lat) / lat_width)

    maidenhead_features = []

    total_cells = (end_x - start_x + 1) * (end_y - start_y + 1)

    # Loop through all intersecting grid cells with tqdm progress bar
    with tqdm(total=total_cells, desc="Generating Maidenhead DGGS") as pbar:
        for x in range(start_x, end_x + 1):
            for y in range(start_y, end_y + 1):
                # Calculate the cell bounds
                cell_min_lon = base_lon + x * lon_width
                cell_max_lon = cell_min_lon + lon_width
                cell_min_lat = base_lat + y * lat_width
                cell_max_lat = cell_min_lat + lat_width

                # Ensure the cell intersects with the bounding box
                if not (cell_max_lon < min_lon or cell_min_lon > max_lon or
                        cell_max_lat < min_lat or cell_min_lat > max_lat):
                    # Center point for the Maidenhead code
                    cell_center_lat = (cell_min_lat + cell_max_lat) / 2
                    cell_center_lat = (cell_min_lon + cell_max_lon) / 2
                    
                    maidenhead_id = maidenhead.toMaiden(cell_center_lat, cell_center_lat, resolution)
                    _, _, min_lat_maiden, min_lon_maiden, max_lat_maiden, max_lon_maiden, _ = maidenhead.maidenGrid(maidenhead_id)
                    # Define the polygon based on the bounding box
                    cell_polygon = Polygon([
                        [min_lon_maiden, min_lat_maiden],  # Bottom-left corner
                        [max_lon_maiden, min_lat_maiden],  # Bottom-right corner
                        [max_lon_maiden, max_lat_maiden],  # Top-right corner
                        [min_lon_maiden, max_lat_maiden],  # Top-left corner
                        [min_lon_maiden, min_lat_maiden]   # Closing the polygon (same as the first point)
                    ])
                    
                    maidenhead_feature = graticule_dggs_to_feature('maidenhead',maidenhead_id,resolution,cell_polygon)            
                 
                    maidenhead_features.append(maidenhead_feature)

                pbar.update(1)

    return {
        "type": "FeatureCollection",
        "features": maidenhead_features
    }
   
def main():
    parser = argparse.ArgumentParser(description="Generate Maidenhead DGGS")
    parser.add_argument('-r', '--resolution', type=int, choices=[1, 2, 3, 4], default=1,
                        help="resolution [1..4]")
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    ) 
    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    
    if bbox == [-180, -90, 180, 90]:
        # Calculate the number of cells at the given resolution
        num_cells,_,_ = maidenhead_metrics(resolution)
        print(f"Resolution {resolution} will generate {num_cells} cells ")
        if num_cells > max_cells:
            print(f"which exceeds the limit of {max_cells}.")
            print("Please select a smaller resolution and try again.")
            return
        geojson_features = generate_grid(resolution)    
    else:
        geojson_features = generate_grid_within_bbox(resolution, bbox)    
   
    if geojson_features:
        # Define the GeoJSON file path
        geojson_path = f"maidenhead_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
        
if __name__ == "__main__":
    main()