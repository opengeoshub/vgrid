
# https://github.com/corteva/georef-field
from vgrid.utils import georef
import json, argparse
from tqdm import tqdm
from shapely.geometry import Polygon, mapping, box
import numpy as np
from pyproj import Geod

max_cells = 1_000_000

geod = Geod(ellps="WGS84")  # Initialize a Geod object for calculations

def generate_grid(resolution):
    # Default to the whole world if no bounding box is provided
    lon_min, lat_min, lon_max, lat_max = -180, -90, 180, 90
    
    resolution_degrees = {
        -1: 15.0,    # 15° x 15°
        0: 1.0,     # 1° x 1°
        1: 1 / 60,  # 1-minute
        2: 1 / 600, # 0.1-minute
        3: 1 / 6000, # 0.01-minute
        4: 1 / 60_000, # 0.001-minute,
        5: 1 / 600_000, 
        6: 1 / 6_000_000, 
        7: 1 / 60_000_000,
        8: 1 / 600_000_000, 
        9: 1 / 6000_000_000
    }[resolution]

    longitudes = np.arange(lon_min, lon_max, resolution_degrees)
    latitudes = np.arange(lat_min, lat_max, resolution_degrees)
    # Use numpy to generate ranges with floating-point steps
    total_cells = len(longitudes) * len(latitudes)
    
    features = []
    # Loop over longitudes and latitudes in 30-minute intervals
    with tqdm(total=total_cells, desc="Generating GEOREF grid", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                # Create a polygon for each GARS cell
                cell_polygon = Polygon(box(lon, lat, lon + resolution_degrees, lat + resolution_degrees))
                georef_code = georef.encode(lat,lon,resolution) 
                # print(poly)
                features.append({
                        "type": "Feature",
                        "geometry": mapping(cell_polygon),
                        "properties": {
                            "georef": georef_code
                        },
                    })

                pbar.update(1)

    return {
            "type": "FeatureCollection",
            "features": features,
        }
 
def generate_grid_within_bbox(bbox, resolution_minutes):
    features = []
    return {
            "type": "FeatureCollection",
            "features": features,
        }


def main():
    parser = argparse.ArgumentParser(description="Generate GEOREF grid")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True,
        help="Resolution in range[0..10]"
    )
    parser.add_argument(
        "-b", "--bbox", type=float, nargs=4,
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )

    args = parser.parse_args()
    resolution =args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]

    if resolution < 0 or resolution > 10:
        print(f"Please select a resolution in [0..10] range and try again ")
        return

    if bbox == [-180, -90, 180, 90]:
        #  # Calculate the number of cells at the given resolution
        # lon_min, lat_min, lon_max, lat_max = bbox
        # resolution_degrees = resolution / 60.0
        # longitudes = np.arange(lon_min, lon_max, resolution_degrees)
        # latitudes = np.arange(lat_min, lat_max, resolution_degrees)

        # total_cells = len(longitudes) * len(latitudes)

        # print(f"Resolution {resolution} will generate {total_cells} cells "
        # if total_cells > max_cells:
        #     print(f"which exceed the limit of {max_cells}.")
        #     print("Please select a smaller resolution and try again.")
        #     return
        feature_collection = generate_grid(resolution)
        output_filename = f'georef_grid_{resolution}.geojson'
    else: 
        feature_collection = generate_grid_within_bbox(bbox, resolution)
        output_filename = f'georef_grid_{resolution}_bbox.geojson'
   
    with open(output_filename, 'w') as f:
        json.dump(feature_collection, f, indent=2)

    print(f"georef grid saved to {output_filename}")

if __name__ == "__main__":
    main()