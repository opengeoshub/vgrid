# https://github.com/ha8tks/Leaflet.Maidenhead
# https://ha8tks.github.io/Leaflet.Maidenhead/examples/
# https://www.sotamaps.org/

import json, math
import argparse
from vgrid.utils import maidenhead
from vgrid.stats.maidenheadstats import maidenhead_metrics
from tqdm import tqdm  
import locale

current_locale = locale.getlocale()  # Get the current locale setting
locale.setlocale(locale.LC_ALL,current_locale)  # Use the system's default locale
max_cells = 1_000_000

def generate_grid(resolution, bbox=None):
    if resolution == 1:
        lon_width, lat_width = 20, 10
    elif resolution == 2:
        lon_width, lat_width = 2, 1
    elif resolution == 3:
        lon_width, lat_width = 0.2, 0.1
    elif resolution == 4:
        lon_width, lat_width = 0.02, 0.01
    else:
        raise ValueError("Unsupported resolution")

    # Determine the bounding box
    min_lon, min_lat, max_lon, max_lat = bbox or [-180, -90, 180, 90]
    x_cells = int((max_lon - min_lon) / lon_width)
    y_cells = int((max_lat - min_lat) / lat_width)
    total_cells = x_cells * y_cells

    features = []
    with tqdm(total=total_cells, desc="Generating grid", unit=" cells") as pbar:
        for i in range(x_cells):
            for j in range(y_cells):
                cell_min_lon = min_lon + i * lon_width
                cell_max_lon = cell_min_lon + lon_width
                cell_min_lat = min_lat + j * lat_width
                cell_max_lat = cell_min_lat + lat_width

                center_lat = (cell_min_lat + cell_max_lat) / 2
                center_lon = (cell_min_lon + cell_max_lon) / 2
                maidenhead_code = maidenhead.toMaiden(center_lat, center_lon, resolution)

                polygon_coords = [
                    [cell_min_lon, cell_min_lat],
                    [cell_max_lon, cell_min_lat],
                    [cell_max_lon, cell_max_lat],
                    [cell_min_lon, cell_max_lat],
                    [cell_min_lon, cell_min_lat]
                ]

                features.append({
                    "type": "Feature",
                    "geometry": {
                        "type": "Polygon",
                        "coordinates": [polygon_coords]
                    },
                    "properties": {
                        "maidenhead": maidenhead_code
                    }
                })
                pbar.update(1)

    return {
        "type": "FeatureCollection",
        "features": features
    }


def generate_grid_within_bbox(resolution, bbox):
    # Define the grid parameters based on the resolution
    if resolution == 1:
        lon_width, lat_width = 20, 10
    elif resolution == 2:
        lon_width, lat_width = 2, 1
    elif resolution == 3:
        lon_width, lat_width = 0.2, 0.1
    elif resolution == 4:
        lon_width, lat_width = 0.02, 0.01
    else:
        raise ValueError("Unsupported resolution")

    min_lon, min_lat, max_lon, max_lat = bbox

    # Calculate grid cell indices for the bounding box
    base_lat, base_lon = -90, -180
    start_x = math.floor((min_lon - base_lon) / lon_width)
    end_x = math.floor((max_lon - base_lon) / lon_width)
    start_y = math.floor((min_lat - base_lat) / lat_width)
    end_y = math.floor((max_lat - base_lat) / lat_width)

    features = []

    total_cells = (end_x - start_x + 1) * (end_y - start_y + 1)

    # Loop through all intersecting grid cells with tqdm progress bar
    with tqdm(total=total_cells, desc="Generating cells") as pbar:
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
                    center_lat = (cell_min_lat + cell_max_lat) / 2
                    center_lon = (cell_min_lon + cell_max_lon) / 2

                    maidenhead_code = maidenhead.toMaiden(center_lat, center_lon, resolution)

                    # Create GeoJSON feature for the cell
                    polygon_coords = [
                        [cell_min_lon, cell_min_lat],  # Bottom-left
                        [cell_max_lon, cell_min_lat],  # Bottom-right
                        [cell_max_lon, cell_max_lat],  # Top-right
                        [cell_min_lon, cell_max_lat],  # Top-left
                        [cell_min_lon, cell_min_lat]   # Closing the polygon
                    ]

                    feature = {
                        "type": "Feature",
                        "geometry": {
                            "type": "Polygon",
                            "coordinates": [polygon_coords]
                        },
                        "properties": {
                            "maidenhead": maidenhead_code
                        }
                    }

                    features.append(feature)

                # Update the progress bar
                pbar.update(1)

    return {
        "type": "FeatureCollection",
        "features": features
    }
   
def main():
    parser = argparse.ArgumentParser(description="Generate Maidenhead grid cells and save as GeoJSON")
    parser.add_argument('-r', '--resolution', type=int, choices=[1, 2, 3, 4], default=1,
                        help="resolution for Maidenhead grid (1 to 4)")
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
        print(f"Resolution {resolution} will generate "
              f"{locale.format_string('%d', num_cells, grouping=True)} cells, ")
        if num_cells > max_cells:
            print(f"which exceeds the limit of {locale.format_string('%d', max_cells, grouping=True)}.")
            print("Please select a smaller resolution and try again.")
            return

        # Generate grid within the bounding box
        geojson_features = generate_grid(resolution)

        # Define the GeoJSON file path
        geojson_path = f"maidenhead_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
    
    else:
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = f"maidenhead_grid_{resolution}_bbox.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")
            
if __name__ == "__main__":
    main()