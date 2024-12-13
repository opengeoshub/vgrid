#Reference: 
# https://github.com/aaliddell/s2cell, 
# https://medium.com/@claude.ducharme/selecting-a-geo-representation-81afeaf3bf01
# https://github.com/sidewalklabs/s2
# https://github.com/google/s2geometry/tree/master/src/python
# https://github.com/google/s2geometry
# https://gis.stackexchange.com/questions/293716/creating-shapefile-of-s2-cells-for-given-level
# https://s2.readthedocs.io/en/latest/quickstart.html
from vgrid.utils import s2
import json
import argparse
from tqdm import tqdm

def create_s2_grid(resolution):
    # Define the cell level (S2 uses a level system for zoom, where level 30 is the highest resolution)
    level = resolution

    # Create a list to store the S2 cell IDs
    cell_ids = []

    # Define the cell covering
    coverer = s2.RegionCoverer()
    coverer.min_level = level
    coverer.max_level = level
    coverer.max_cells = 1000000  # Adjust as needed
    # coverer.max_cells = 0  # Adjust as needed


    # Define the region to cover (in this example, we'll use the entire world)
    region = s2.LatLngRect(s2.LatLng.from_degrees(-90, -180),
                                 s2.LatLng.from_degrees(90, 180))

    # Get the covering cells
    covering = coverer.get_covering(region)

    # Convert the covering cells to S2 cell IDs
    for cell_id in covering:
        cell_ids.append(cell_id)

    return cell_ids

def cell_to_polygon(cell_id):
    cell = s2.Cell(cell_id)
    vertices = []
    for i in range(4):
        vertex = s2.LatLng.from_point(cell.get_vertex(i))
        vertices.append((vertex.lng().degrees, vertex.lat().degrees))
    
    vertices.append(vertices[0])  # Close the polygon
  # Return a dictionary representing the GeoJSON polygon
    polygon = {
        "type": "Polygon",
        "coordinates": [vertices]
    }
    
    return polygon

def save_s2_grid_as_geojson(cell_ids, output_filename):
    features = []
    for cell_id in tqdm(cell_ids, desc="Saving cells"):
        polygon = cell_to_polygon(cell_id)
        
        # Create a feature dictionary
        feature = {
            "type": "Feature",
            "geometry": polygon,
            "properties": {"cell_id": cell_id.to_token()}
        }
        
        features.append(feature)
    
    # Create a FeatureCollection
    feature_collection = {
        "type": "FeatureCollection",
        "features": features
    }
    
    # Save the feature collection to a JSON file
    with open(output_filename, 'w') as f:
        json.dump(feature_collection, f)

def main():
    parser = argparse.ArgumentParser(description="Generate S2 grid at a specific resolution and save as GeoJSON.")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="resolution for the S2 grid.")
    parser.add_argument('-o', '--output', type=str, required=True, help="Output GeoJSON file name.")
    args = parser.parse_args()
    
    res = args.resolution
    output_filename = args.output
    
    cell_ids = create_s2_grid(res)
    save_s2_grid_as_geojson(cell_ids, output_filename)
    print(f"S2 grid at resolution {res} saved to {output_filename}")

if __name__ == "__main__":
    main()

