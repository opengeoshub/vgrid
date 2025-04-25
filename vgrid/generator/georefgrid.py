from vgrid.utils import georef
import json, argparse
from tqdm import tqdm
from shapely.geometry import Polygon, box
import numpy as np

from vgrid.generator.settings import max_cells, graticule_dggs_to_feature

RESOLUTION_DEGREES = {
    -1: 15.0,       # 15° x 15°
    0: 1.0,        # 1° x 1°
    1: 1 / 60,     # 1-minute
    2: 1 / 600,    # 0.1-minute
    3: 1 / 6000,   # 0.01-minute
    4: 1 / 60_000,  # 0.001-minute
    5: 1 / 600_000  # 0.0001-minute
}

# RESOLUTION_DEGREES = {
#     0: 15.0,       # 15° x 15°
#     1: 1.0,        # 1° x 1°
#     2: 1 / 60,     # 1-minute
#     3: 1 / 600,    # 0.1-minute
#     5: 1 / 6000,   # 0.01-minute
#     5: 1 / 60_000,  # 0.001-minute
#     # 5: 1 / 600_000  # 0.0001-minute
# }


def generate_grid(bbox, resolution):
    lon_min, lat_min, lon_max, lat_max = bbox    
    resolution_degrees = RESOLUTION_DEGREES[resolution]
    longitudes = np.arange(lon_min, lon_max, resolution_degrees)
    latitudes = np.arange(lat_min, lat_max, resolution_degrees)
    num_cells = len(longitudes) * len(latitudes)
    
    print(f"Resolution {resolution} will generate {num_cells} cells ")
    if num_cells > max_cells:
        print(f"which exceeds the limit of {max_cells}.")
        print("Please select a smaller resolution and try again.")
        return
    georef_features = []
    
    with tqdm(total=num_cells, desc="Generating GEOREF DGGS", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                cell_polygon = Polygon(box(lon, lat, lon + resolution_degrees, lat + resolution_degrees))
                georef_id = georef.encode(lat, lon, resolution)
                georef_feature = graticule_dggs_to_feature('georef', georef_id, resolution,cell_polygon)
                georef_features.append(georef_feature)
                pbar.update(1)
    
    return {
        "type": "FeatureCollection",
        "features": georef_features,
    }


def main():
    parser = argparse.ArgumentParser(description="Generate GEOREF DGGS")
    parser.add_argument(
        "-r", "--resolution", type=int, required=True,
        help="Resolution [-1..5]"
    )
    parser.add_argument(
        "-b", "--bbox", type=float, nargs=4,
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )

    args = parser.parse_args()
    resolution =args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]

    if resolution < -1 or resolution > 5:
        print(f"Please select a resolution in [-1..5] range and try again ")
        return
    
    geojson_features = generate_grid(bbox, resolution)
    output_filename = f'georef_grid_{resolution}.geojson'
        
    with open(output_filename, 'w') as f:
        json.dump(geojson_features, f, indent=2)

    print(f"GEOREF grid saved to {output_filename}")

if __name__ == "__main__":
    main()