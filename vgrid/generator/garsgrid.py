import json, argparse
from tqdm import tqdm
from shapely.geometry import Polygon, mapping, box
import numpy as np
from vgrid.utils.gars.garsgrid import GARSGrid  # Ensure the correct import path

from vgrid.generator.settings import max_cells, graticule_dggs_to_feature


def generate_grid(resolution_minutes):
    # Default to the whole world if no bounding box is provided
    lon_min, lat_min, lon_max, lat_max = -180, -90, 180, 90

    resolution_degrees = resolution_minutes / 60.0

    # Initialize a list to store GARS grid features

    # Generate ranges for longitudes and latitudes
    longitudes = np.arange(lon_min, lon_max, resolution_degrees)
    latitudes = np.arange(lat_min, lat_max, resolution_degrees)

    total_cells = len(longitudes) * len(latitudes)
        
    gars_features = []
    # Loop over longitudes and latitudes with tqdm progress bar
    with tqdm(total=total_cells, desc="Generating GARS grid", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                # Create the GARS grid code
                gars_cell= GARSGrid.from_latlon(lat, lon, resolution_minutes)
                wkt_polygon = gars_cell.polygon
                
                if wkt_polygon:
                    cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
                    gars_id = gars_cell.gars_id
                    gars_feature  = graticule_dggs_to_feature('gars',gars_id, resolution_minutes, cell_polygon)
                    gars_features.append(gars_feature)
                    pbar.update(1)

    # Create a FeatureCollection
    return {
            "type": "FeatureCollection",
            "features": gars_features,
        }
 
def generate_grid_within_bbox(bbox, resolution_minutes):
    # Default to the whole world if no bounding box is provided
    bbox_polygon = box(*bbox)
    lon_min, lat_min, lon_max, lat_max = bbox

    resolution_degrees = resolution_minutes / 60.0
    
    longitudes = np.arange(lon_min-resolution_degrees, lon_max + resolution_degrees, resolution_degrees)
    latitudes = np.arange(lat_min-resolution_degrees, lat_max + resolution_degrees, resolution_degrees)

    # total_cells = len(longitudes) * len(latitudes)
    gars_features = []
    # Loop over longitudes and latitudes with tqdm progress bar
    with tqdm(desc="Generating GARS grid", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                # Create the GARS grid code
                gars_cell= GARSGrid.from_latlon(lat, lon, resolution_minutes)
                wkt_polygon = gars_cell.polygon
                
                if wkt_polygon:                   
                    cell_polygon = Polygon(list(wkt_polygon.exterior.coords))                
                   
                    if bbox_polygon.intersects(cell_polygon):
                        gars_id = gars_cell.gars_id                       
                        gars_feature  = graticule_dggs_to_feature('gars',gars_id, resolution_minutes, cell_polygon)
                        gars_features.append(gars_feature)
                        pbar.update(1)

    # Create a FeatureCollection
    return {
            "type": "FeatureCollection",
            "features": gars_features,
        }


def main():
    parser = argparse.ArgumentParser(description="Generate GARS grid")
    parser.add_argument(
        "-r", "--resolution", type=int, choices=[30, 15, 5, 1], required=True,
        help="Resolution in minutes (30, 15, 5, 1)"
    )
    parser.add_argument(
        "-b", "--bbox", type=float, nargs=4,
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )

    args = parser.parse_args()
    resolution_minutes=args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]

    # Write to a GeoJSON file
    if bbox == [-180, -90, 180, 90]:
         # Calculate the number of cells at the given resolution
        lon_min, lat_min, lon_max, lat_max = -180, -90, 180, 90

        resolution_degrees = resolution_minutes / 60.0

        # Initialize a list to store GARS grid features

        # Generate ranges for longitudes and latitudes
        longitudes = np.arange(lon_min, lon_max, resolution_degrees)
        latitudes = np.arange(lat_min, lat_max, resolution_degrees)

        total_cells = len(longitudes) * len(latitudes)
        print(f"Resolution {resolution_minutes} minutes will generate {total_cells} cells ")
        if total_cells > max_cells:
            print(f"which exceeds the limit of {max_cells}.")
            print("Please select a smaller resolution and try again.")
            return
    
        geojson_features = generate_grid(resolution_minutes)
    
    else: 
        geojson_features = generate_grid_within_bbox(bbox, resolution_minutes)
    
    output_filename = f'gars_grid_{resolution_minutes}_minutes.geojson'
   
    with open(output_filename, 'w') as f:
        json.dump(geojson_features, f, indent=2)

    print(f"GARS grid saved to {output_filename}")


if __name__ == "__main__":
    main()
