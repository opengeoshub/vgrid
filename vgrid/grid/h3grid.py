#Reference: https://observablehq.com/@claude-ducharme/h3-map
# https://h3-snow.streamlit.app/

import h3
from shapely.geometry import Polygon, mapping
import json

# def wrap_longitudes(boundary):
#     """
#     Adjust longitudes of a hexagon boundary to wrap correctly around the antimeridian.

#     Parameters:
#         boundary (list): List of (lat, lon) tuples representing the hexagon boundary.

#     Returns:
#         list: Adjusted boundary with longitudes wrapped as needed.
#     """
#     wrapped_boundary = []
#     for lat, lon in boundary:
#         # Wrap longitude if greater than 180 degrees
#         if lon > 180:
#             lon -= 360
#         wrapped_boundary.append((lat, lon))
#     return wrapped_boundary

def filter_antimeridian_cells(hex_boundary, threshold=-128):
    """
    Filters and adjusts hexagons crossing the antimeridian.

    Parameters:
        hex_boundary (list): List of (lat, lon) tuples representing the hexagon boundary.
        threshold (float): Longitude threshold to identify crossing cells.

    Returns:
        list: Adjusted boundary if it crosses the antimeridian.
    """
    # Check if any longitude in the boundary is below the threshold
    if any(lon < threshold for _, lon in hex_boundary):
        # Adjust all longitudes accordingly
        return [(lat, lon - 360 if lon > 0 else lon) for lat, lon in hex_boundary]
    return hex_boundary

def generate_h3_geojson_with_filter(precision, output_file):
    """
    Generate H3 cells at a specified precision, filter for antimeridian crossings, 
    and save them as a GeoJSON file.

    Parameters:
        precision (int): The H3 resolution (0-15).
        output_file (str): Path to save the GeoJSON file.

    Returns:
        None
    """
    if not (0 <= precision <= 15):
        raise ValueError("Precision must be between 0 and 15.")

    base_cells = h3.get_res0_cells()
    features = []

    for cell in base_cells:
        child_cells = h3.cell_to_children(cell, precision)
        for child_cell in child_cells:
            # Get the boundary of the cell
            hex_boundary = h3.cell_to_boundary(child_cell)
            # Wrap and filter the boundary
            filtered_boundary = filter_antimeridian_cells(hex_boundary)
            # Reverse lat/lon to lon/lat for GeoJSON compatibility
            reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
            polygon = Polygon(reversed_boundary)
            if polygon.is_valid:
                features.append({
                    "type": "Feature",
                    "geometry": mapping(polygon),
                    "properties": {
                        "h3_index": child_cell
                    }
                })

    geojson_data = {
        "type": "FeatureCollection",
        "features": features
    }

    with open(output_file, 'w') as f:
        json.dump(geojson_data, f)

    print(f"GeoJSON saved to {output_file}")

# Example Usage
precision = 1 # Choose a precision level
output_file = f"h3_filtered_cells_precision_{precision}.geojson"
generate_h3_geojson_with_filter(precision, output_file)
