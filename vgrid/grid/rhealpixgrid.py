import argparse
import json
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from shapely.geometry import Polygon, box, mapping
from tqdm import tqdm

# Function to filter cells crossing the antimeridian
def filter_antimeridian_cells(boundary, threshold=-128):
    if any(lon < threshold for lon, _ in boundary):
        return [(lon - 360 if lon > 0 else lon, lat) for lon, lat in boundary]
    return boundary

# Function to convert cell vertices to a Shapely Polygon
def cell_to_polygon(cell):
    vertices = [tuple(my_round(coord, 14) for coord in vertex) for vertex in cell.vertices(plane=False)]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    vertices = filter_antimeridian_cells(vertices)
    return Polygon(vertices)

# Function to process grid cells and return GeoJSON features within a bounding box
def generate_grid_within_bbox(rhealpix_dggs, resolution, bbox):
    features = []
    total_cells = rhealpix_dggs.num_cells(resolution)
    rhealpix_grid = rhealpix_dggs.grid(resolution)

    bbox_polygon = box(*bbox)  # Create a bounding box polygon
    with tqdm(total=total_cells, desc="Processing grid cells", unit="cell") as pbar:
        for cell in rhealpix_grid:
            polygon = cell_to_polygon(cell)
            if polygon.intersects(bbox_polygon):  # Only include cells that intersect the bbox
                features.append({
                    "type": "Feature",
                    "geometry": mapping(polygon),
                    "properties": {"rhealpix": str(cell)},
                })
            pbar.update(1)

    return {
        "type": "FeatureCollection",
        "features": features,
    }

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate RHEALPix grid within a bounding box and save as a GeoJSON.")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid")
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )
    args = parser.parse_args()

    # Initialize RHEALPix DGGS
    rhealpix_dggs = RHEALPixDGGS()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]

    # Calculate the number of cells at the given resolution
    num_cells = rhealpix_dggs.num_cells(resolution)
    max_cells = 1_000_000
    if num_cells > max_cells:
        print(f"Error: The selected resolution will generate {num_cells:,} cells, which exceeds the limit of {max_cells:,}.")
        print("Please select a smaller resolution and try again.")
        return

    # Generate grid within the bounding box
    geojson_features = generate_grid_within_bbox(rhealpix_dggs, resolution, bbox)

    # Define the GeoJSON file path
    geojson_path = f"rhealpix_grid_{resolution}_bbox.geojson"
    with open(geojson_path, 'w') as f:
        json.dump(geojson_features, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")

if __name__ == "__main__":
    main()
