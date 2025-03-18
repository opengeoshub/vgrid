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
from vgrid.utils.antimeridian import fix_polygon
from shapely.geometry import Polygon, mapping
from vgrid.generator.settings import max_cells, geodesic_dggs_to_feature

def s2_cell_to_polygon(cell_id):
    cell = s2.Cell(cell_id)
    vertices = []
    for i in range(4):
        vertex = s2.LatLng.from_point(cell.get_vertex(i))
        vertices.append((vertex.lng().degrees, vertex.lat().degrees))
    
    vertices.append(vertices[0])  # Close the polygon
    
    # Create a Shapely Polygon
    polygon = Polygon(vertices)
    #  Fix Antimerididan:
    fixed_polygon = fix_polygon(polygon)    
    return fixed_polygon

def generate_grid(resolution,bbox):
    min_lng, min_lat, max_lng, max_lat = bbox
     # Define the cell level (S2 uses a level system for zoom, where level 30 is the highest resolution)
    level = resolution
    # Create a list to store the S2 cell IDs
    cell_ids = []
    # Define the cell covering
    coverer = s2.RegionCoverer()
    coverer.min_level = level
    coverer.max_level = level
    # coverer.max_cells = 1000_000  # Adjust as needed
    # coverer.max_cells = 0  # Adjust as needed

    # Define the region to cover (in this example, we'll use the entire world)
    region = s2.LatLngRect(
        s2.LatLng.from_degrees(min_lat, min_lng),
        s2.LatLng.from_degrees(max_lat, max_lng)
    )

    # Get the covering cells
    covering = coverer.get_covering(region)

    # Convert the covering cells to S2 cell IDs
    for cell_id in covering:
        cell_ids.append(cell_id)

    s2_features = []
    s2_token = cell_id.to_token()
    num_edges = 4
    
    for cell_id in tqdm(cell_ids, desc="Processing cells"):
        # Generate a Shapely Polygon
        cell_polygon = s2_cell_to_polygon(cell_id)
        s2_feature = geodesic_dggs_to_feature("s2",s2_token,resolution,cell_polygon,num_edges)   
        s2_features.append(s2_feature)

    return {
        "type": "FeatureCollection",
        "features": s2_features
    }


def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate S2 grid.")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..30] of the grid")
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    ) 
    args = parser.parse_args()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    
    if resolution < 0 or resolution > 30:
        print(f"Please select a resolution in [0..30] range and try again ")
        return
    
    geojson_features = generate_grid(resolution,bbox)
    # Define the GeoJSON file path
    geojson_path = f"s2_grid_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump(geojson_features, f, indent=2)

if __name__ == "__main__":
    main()