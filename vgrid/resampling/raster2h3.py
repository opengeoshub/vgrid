from vgrid.utils import s2
from shapely.geometry import Point, LineString, Polygon, mapping, box
import argparse
import json
from tqdm import tqdm
import os
from vgrid.generator.s2grid import cell_to_polygon
from pyproj import Geod
geod = Geod(ellps="WGS84")
from shapely.ops import unary_union
from shapely.geometry import shape, mapping, MultiPolygon

def extract_bbox(geojson_data):
    """Extracts the bounding box from a GeoJSON file."""
    geometries = [shape(feature["geometry"]) for feature in geojson_data["features"]]
    unified_geom = unary_union(geometries)  # Merge multiple features if needed
    min_lng, min_lat, max_lng, max_lat = unified_geom.bounds  # Get BBOX
    return [min_lng, min_lat, max_lng, max_lat]

def generate_grid(resolution, geojson_data):
    """Generates S2 grid cells within the GeoJSON bounding box and keeps only fully contained cells."""
    bbox = extract_bbox(geojson_data)  # Extract bounding box from GeoJSON
    min_lng, min_lat, max_lng, max_lat = bbox
    level = resolution
    cell_ids = []

    coverer = s2.RegionCoverer()
    coverer.min_level = level
    coverer.max_level = level

    region = s2.LatLngRect(
        s2.LatLng.from_degrees(min_lat, min_lng),
        s2.LatLng.from_degrees(max_lat, max_lng)
    )

    covering = coverer.get_covering(region)

    for cell_id in covering:
        cell_ids.append(cell_id)

    # Merge all features in geojson_data into a single geometry
    geojson_shape = unary_union([shape(feature["geometry"]) for feature in geojson_data["features"]])

    features = []
    for cell_id in tqdm(cell_ids, desc="Processing cells"):
        polygon = cell_to_polygon(cell_id)  # Convert S2 cell to Polygon
        if polygon.within(geojson_shape):  # Keep only fully contained cells
            geometry = mapping(polygon)
            feature = {
                "type": "Feature",
                "geometry": geometry,
                "properties": {"s2": cell_id.to_token()}
            }
            features.append(feature)

    return {
        "type": "FeatureCollection",
        "features": features,
    }


def calculate_population(layer1, layer2):
    """Filters Layer 2 and assigns population based on area proportion from Layer 1."""
    
    layer1_features = [
        (shape(feature["geometry"]), feature["properties"].get("population", 0))
        for feature in layer1["features"]
    ]
    
    filtered_features = []

    for feature in layer2["features"]:
        layer2_shape = shape(feature["geometry"])
        total_population = 0
        intersected_parts = []
        
        for l1_shape, l1_population in layer1_features:
            if layer2_shape.intersects(l1_shape):
                intersection = layer2_shape.intersection(l1_shape)
                intersected_parts.append(intersection)

                # Area-weighted population calculation
                if not intersection.is_empty:
                    proportion = intersection.area / l1_shape.area
                    total_population += l1_population * proportion

        if not intersected_parts:
            continue  # Skip features that have no intersection

        # Add population to properties
        feature["properties"]["population"] = round(total_population, 2)
        filtered_features.append(feature)

    return {
        "type": "FeatureCollection",
        "features": filtered_features
    }


# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Generate H3 grid for shapes in GeoJSON format")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..30]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON string with Point, Polyline or Polygon"
    )
    args = parser.parse_args()
    geojson = args.geojson
     # Initialize h3 DGGS
    resolution = args.resolution
    
    if resolution < 0 or resolution > 30:
        print(f"Please select a resolution in [0..30] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson) as f:
        geojson_data = json.load(f)
    
    geojson_features = generate_grid(resolution,geojson_data)
    geojson_features = calculate_population(geojson_data, geojson_features)
    
    # Define the GeoJSON file path

    geojson_path = f"s2_resample_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump(geojson_features, f, indent=2)

    print(f"Resample GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
