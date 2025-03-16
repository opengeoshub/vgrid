# from gars_field import EDGARSGrid

# # from latlon
# ggrid = EDGARSGrid.from_latlon(-89.55, -179.57, resolution=6)
# print(ggrid)
# # from GARS ID
# # ggrid = EDGARSGrid("D01AA23")

# # get bounding poly
# grid_poly = ggrid.polygon
# print(grid_poly)
# # get GARS ID
# gars_id = str(ggrid)
# print(gars_id)
# # UTM CRS EPSG Code
# epsg_code = ggrid.utm_epsg
# print(epsg_code)
# https://github.com/Moustikitos/gryd/blob/c79edde94f19d46e3b3532ae14eb351e91d55322/Gryd/geodesy.py

import json, argparse
from tqdm import tqdm
from shapely.geometry import Polygon, mapping, box
import numpy as np
from pyproj import Geod
from vgrid.utils.gars.garsgrid import GARSGrid  # Ensure the correct import path

from vgrid.generator.settings import max_cells


geod = Geod(ellps="WGS84")  # Initialize a Geod object for calculations

def generate_grid(resolution_minutes):
    # Default to the whole world if no bounding box is provided
    lon_min, lat_min, lon_max, lat_max = -180, -90, 180, 90

    resolution_degrees = resolution_minutes / 60.0

    # Initialize a list to store GARS grid features

    # Generate ranges for longitudes and latitudes
    longitudes = np.arange(lon_min, lon_max, resolution_degrees)
    latitudes = np.arange(lat_min, lat_max, resolution_degrees)

    total_cells = len(longitudes) * len(latitudes)
    features = []
    # Loop over longitudes and latitudes with tqdm progress bar
    with tqdm(total=total_cells, desc="Generating GARS grid", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                # Create the GARS grid code
                gars_cell= GARSGrid.from_latlon(lat, lon, resolution_minutes)
                wkt_polygon = gars_cell.polygon
                
                if wkt_polygon:
                    # Extract coordinates and create polygon
                    # x, y = wkt_polygon.exterior.xy
                    # min_lon = min(x)
                    # max_lon = max(x)
                    # min_lat = min(y)
                    # max_lat = max(y)

                    # Calculate center coordinates
                    # center_lon = round((min_lon + max_lon) / 2, 7)
                    # center_lat = round((min_lat + max_lat) / 2, 7)

                    # Create a shapely polygon
                    cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
                    
                    # Calculate area, width, and height
                    # cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]), 2)
                    # cell_width = round(geod.line_length([min_lon, max_lon], [min_lat, min_lat]), 2)
                    # cell_height = round(geod.line_length([min_lon, min_lon], [min_lat, max_lat]), 2)

                    # Add feature to the list
                    features.append({
                        "type": "Feature",
                        "geometry": mapping(cell_polygon),
                        "properties": {
                            "gars": gars_cell.gars_id,
                            # "center_lat": center_lat,
                            # "center_lon": center_lon,
                            # "cell_area": cell_area,
                            # "cell_width": cell_width,
                            # "cell_height": cell_height,
                            # "resolution_minute": gars_grid.resolution,
                        },
                    })

                pbar.update(1)

    # Create a FeatureCollection
    return {
            "type": "FeatureCollection",
            "features": features,
        }
 
def generate_grid_within_bbox(bbox, resolution_minutes):
    # Default to the whole world if no bounding box is provided
    bbox_polygon = box(*bbox)
    lon_min, lat_min, lon_max, lat_max = bbox

    resolution_degrees = resolution_minutes / 60.0

    # Generate ranges for longitudes and latitudes
    # longitudes = np.arange(lon_min, lon_max, resolution_degrees)
    # latitudes = np.arange(lat_min, lat_max, resolution_degrees)
    
    longitudes = np.arange(lon_min-resolution_degrees, lon_max + resolution_degrees, resolution_degrees)
    latitudes = np.arange(lat_min-resolution_degrees, lat_max + resolution_degrees, resolution_degrees)

    # total_cells = len(longitudes) * len(latitudes)
    features = []
    # Loop over longitudes and latitudes with tqdm progress bar
    with tqdm(desc="Generating GARS grid", unit=" cells") as pbar:
        for lon in longitudes:
            for lat in latitudes:
                # Create the GARS grid code
                gars_cell= GARSGrid.from_latlon(lat, lon, resolution_minutes)
                wkt_polygon = gars_cell.polygon
                
                if wkt_polygon:
                    # Extract coordinates and create polygon
                    # x, y = wkt_polygon.exterior.xy
                    # min_lon = min(x)
                    # max_lon = max(x)
                    # min_lat = min(y)
                    # max_lat = max(y)

                    # Calculate center coordinates
                    # center_lon = round((min_lon + max_lon) / 2, 7)
                    # center_lat = round((min_lat + max_lat) / 2, 7)

                    # Create a shapely polygon
                    cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
                    
                    # Calculate area, width, and height
                    # cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]), 2)
                    # cell_width = round(geod.line_length([min_lon, max_lon], [min_lat, min_lat]), 2)
                    # cell_height = round(geod.line_length([min_lon, min_lon], [min_lat, max_lat]), 2)

                    # Add feature to the list
                    if bbox_polygon.intersects(cell_polygon):
                        features.append({
                            "type": "Feature",
                            "geometry": mapping(cell_polygon),
                            "properties": {
                                "gars": gars_cell.gars_id,
                                # "center_lat": center_lat,
                                # "center_lon": center_lon,
                                # "cell_area": cell_area,
                                # "cell_width": cell_width,
                                # "cell_height": cell_height,
                                # "resolution_minute": gars_grid.resolution,
                            },
                        })

                        pbar.update(1)

    # Create a FeatureCollection
    return {
            "type": "FeatureCollection",
            "features": features,
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
        lon_min, lat_min, lon_max, lat_max = bbox
        resolution_degrees = resolution_minutes / 60.0
        longitudes = np.arange(lon_min, lon_max, resolution_degrees)
        latitudes = np.arange(lat_min, lat_max, resolution_degrees)

        total_cells = len(longitudes) * len(latitudes)

        print(f"Resolution {resolution_minutes} minutes will generate{total_cells} cells ")
        if total_cells > max_cells:
            print(f"which exceeds the limit of {max_cells}.")
            print("Please select a smaller resolution and try again.")
            return
   
        feature_collection = generate_grid(resolution_minutes)
        output_filename = f'gars_grid_{resolution_minutes}_minutes.geojson'
    else: 
        feature_collection = generate_grid_within_bbox(bbox, resolution_minutes)
        output_filename = f'gars_grid_{resolution_minutes}_minutes_bbox.geojson'
   
    with open(output_filename, 'w') as f:
        json.dump(feature_collection, f, indent=2)

    print(f"GARS grid saved to {output_filename}")


if __name__ == "__main__":
    main()
