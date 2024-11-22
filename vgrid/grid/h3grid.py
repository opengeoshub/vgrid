#Reference: https://observablehq.com/@claude-ducharme/h3-map
# https://h3-snow.streamlit.app/

import h3
import json
from shapely.geometry import Polygon, mapping
import argparse
import geopandas as gdp

def generate_h3_hexagons(bbox, resolution):
    """
    Generate H3 hexagons within a bounding box.

    Args:
        bbox (list): Bounding box [min_lat, min_lon, max_lat, max_lon].
        resolution (int): H3 resolution (0 to 15).

    Returns:
        list: List of GeoJSON-like feature dictionaries.
    """
    # Define the bounding box polygon in GeoJSON format
    # geojson_bbox = {
    #     "type": "Polygon",
    #     "coordinates": [[
    #         [bbox[1], bbox[0]],  # [min_lon, min_lat]
    #         [bbox[1], bbox[2]],  # [min_lon, max_lat]
    #         [bbox[3], bbox[2]],  # [max_lon, max_lat]
    #         [bbox[3], bbox[0]],  # [max_lon, min_lat]
    #         [bbox[1], bbox[0]]   # Closing the loop
    #     ]]
    # }
    geojson_bbox = [
        [bbox[1], bbox[0]],  # [min_lon, min_lat]
        [bbox[1], bbox[2]],  # [min_lon, max_lat]
        [bbox[3], bbox[2]],  # [max_lon, max_lat]
        [bbox[3], bbox[0]],  # [max_lon, min_lat]
        [bbox[1], bbox[0]]   # Closing the loop
    ]

    # Generate H3 hexagons using h3.polygon_to_cells
    h3_hexagons = h3.polygon_to_cells(geojson_bbox, resolution)

    # Create GeoJSON features
    features = []
    for hexagon in h3_hexagons:
        hex_boundary = h3.h3_to_geo_boundary(hexagon, geo_json=True)  # Correct method to get the boundary
        polygon = Polygon(hex_boundary)
        features.append({
            "type": "Feature",
            "geometry": mapping(polygon),
            "properties": {"h3_index": hexagon}
        })
    return features

def save_as_geojson(features, output_file):
    """
    Save a list of features to a GeoJSON file.

    Args:
        features (list): List of feature dictionaries.
        output_file (str): Path to the output GeoJSON file.
    """
    feature_collection = {
        "type": "FeatureCollection",
        "features": features
    }
    with open(output_file, 'w') as f:
        json.dump(feature_collection, f, indent=2)

def main():
    parser = argparse.ArgumentParser(description="Generate H3 hexagons and save as GeoJSON.")
    parser.add_argument("-r", "--resolution", type=int, required=True,
                        help="H3 resolution (0 to 15).")
    args = parser.parse_args()

    # Parameters
    bounding_box = [-90, -180, 90, 180]  # Global coverage (min_lat, min_lon, max_lat, max_lon)
    output_file = "h3_hexagons.geojson"  # Output file name

    # Generate and save hexagons
    hex_features = generate_h3_hexagons(bounding_box, args.resolution)
    save_as_geojson(hex_features, output_file)

    print(f"H3 hexagons saved to {output_file} at resolution {args.resolution}")

if __name__ == "__main__":
    main()
