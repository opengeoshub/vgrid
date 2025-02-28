#Reference: https://observablehq.com/@claude-ducharme/h3-map
# https://h3-snow.streamlit.app/

import h3
from shapely.geometry import Polygon, mapping, box
from shapely import buffer
import argparse
import json
from tqdm import tqdm
import locale
from pyproj import Geod
geod = Geod(ellps="WGS84")

current_locale = locale.getlocale()  # Get the current locale setting
locale.setlocale(locale.LC_ALL,current_locale)  # Use the system's default locale
max_cells = 10_000_000

def fix_h3_antimeridian_cells(hex_boundary, threshold=-128):
    if any(lon < threshold for _, lon in hex_boundary):
        # Adjust all longitudes accordingly
        return [(lat, lon - 360 if lon > 0 else lon) for lat, lon in hex_boundary]
    return hex_boundary

def generate_grid(resolution):
    base_cells = h3.get_res0_cells()
    num_base_cells = len(base_cells)
    features = []
    # Progress bar for base cells
    with tqdm(total=num_base_cells, desc="Processing base cells", unit=" cells") as pbar:
        for cell in base_cells:
            child_cells = h3.cell_to_children(cell, resolution)
            # Progress bar for child cells
            for child_cell in child_cells:
                # Get the boundary of the cell
                hex_boundary = h3.cell_to_boundary(child_cell)
                # Wrap and filter the boundary
                filtered_boundary = fix_h3_antimeridian_cells(hex_boundary)
                # Reverse lat/lon to lon/lat for GeoJSON compatibility
                reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
                polygon = Polygon(reversed_boundary)
                if polygon.is_valid:
                    features.append({
                        "type": "Feature",
                        "geometry": mapping(polygon),
                        "properties": {
                            "h3": child_cell
                        }
                    })
            pbar.update(1)

    return {
        "type": "FeatureCollection",
        "features": features,
    }
    
def geodesic_buffer(polygon, distance):
    buffered_coords = []
    for lon, lat in polygon.exterior.coords:
        # Generate points around the current vertex to approximate a circle
        circle_coords = [
            geod.fwd(lon, lat, azimuth, distance)[:2]  # Forward calculation: returns (lon, lat, back_azimuth)
            for azimuth in range(0, 360, 10)  # Generate points every 10 degrees
        ]
        buffered_coords.append(circle_coords)
    
    # Flatten the list of buffered points and form a Polygon
    all_coords = [coord for circle in buffered_coords for coord in circle]
    return Polygon(all_coords).convex_hull

def generate_grid_within_bbox(resolution,bbox):
    bbox_polygon = box(*bbox)  # Create a bounding box polygon
    distance = h3.average_hexagon_edge_length(resolution,unit='m')*2
    bbox_buffer = geodesic_buffer(bbox_polygon, distance)
    bbox_buffer_cells  = h3.geo_to_cells(bbox_buffer,resolution)
    num_cells = len(bbox_buffer_cells)
    print(f"Resolution {resolution} within bounding box {bbox} will generate "
        f"{locale.format_string('%d', num_cells, grouping=True)} cells, ")
    
    if num_cells > max_cells:
        print(f"which exceeds the limit of {locale.format_string('%d', max_cells, grouping=True)}.")
        print("Please select a smaller resolution and try again.")
        return
    else:    
        features = []

        # Progress bar for base cells
        for bbox_buffer_cell in tqdm(bbox_buffer_cells, desc="Processing cells"):
                # Get the boundary of the cell
                hex_boundary = h3.cell_to_boundary(bbox_buffer_cell)
                # Wrap and filter the boundary
                filtered_boundary = fix_h3_antimeridian_cells(hex_boundary)
                # Reverse lat/lon to lon/lat for GeoJSON compatibility
                reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
                cell_polygon = Polygon(reversed_boundary)
                if cell_polygon.intersects(bbox_polygon):
                    features.append({
                        "type": "Feature",
                        "geometry": mapping(cell_polygon),
                        "properties": {
                            "h3": bbox_buffer_cell
                        }
                    })

        return {
            "type": "FeatureCollection",
            "features": features,
        }


# Example Usage
def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate H3 grid within a bounding box and save as a GeoJSON.")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..15] of the grid")
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    ) 
    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    
    if resolution < 0 or resolution > 15:
        print(f"Please select a resolution in [0..15] range and try again ")
        return

    if bbox == [-180, -90, 180, 90]:
        # Calculate the number of cells at the given resolution
        num_cells = h3.get_num_cells(resolution)
        print(f"Resolution {resolution} will generate "
              f"{locale.format_string('%d', num_cells, grouping=True)} cells, ")
        if num_cells > max_cells:
            print(f"which exceeds the limit of {locale.format_string('%d', max_cells, grouping=True)}.")
            print("Please select a smaller resolution and try again.")
            return

        # Generate grid within the bounding box
        geojson_features = generate_grid(resolution)

        # Define the GeoJSON file path
        geojson_path = f"h3_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
    
    else:
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        if geojson_features:
            # Define the GeoJSON file path
            geojson_path = f"h3_grid_{resolution}_bbox.geojson"
            with open(geojson_path, 'w') as f:
                json.dump(geojson_features, f, indent=2)

            print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()