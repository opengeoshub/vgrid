import argparse
import json
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from shapely.geometry import Polygon, box, mapping
from tqdm import tqdm
from pyproj import Geod

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

def generate_grid(rhealpix_dggs, resolution):
    features = []
    total_cells = rhealpix_dggs.num_cells(resolution)
    rhealpix_grid = rhealpix_dggs.grid(resolution)
    with tqdm(total=total_cells, desc="Processing grid cells", unit=" cells") as pbar:
        for cell in rhealpix_grid:
            polygon = cell_to_polygon(cell)
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


def generate_grid_within_bbox(rhealpix_dggs, resolution, bbox):    
    bbox_polygon = box(*bbox)  # Create a bounding box polygon
    bbox_center_lon = bbox_polygon.centroid.x
    bbox_center_lat = bbox_polygon.centroid.y
    seed_point = (bbox_center_lon, bbox_center_lat)

    features = []
    seed_cell = rhealpix_dggs.cell_from_point(resolution, seed_point, plane=False)
    seed_cell_id = str(seed_cell)  # Unique identifier for the current cell
    seed_cell_polygon = cell_to_polygon(seed_cell)
    geod = Geod(ellps="WGS84")

    if seed_cell_polygon.contains(bbox_polygon):
        center_lon = seed_cell_polygon.centroid.x
        center_lat = seed_cell_polygon.centroid.y
        min_x, min_y, max_x, max_y = seed_cell_polygon.bounds
        cell_area = abs(geod.geometry_area_perimeter(seed_cell_polygon)[0])  # Area in square meters                
        _, _, cell_width = geod.inv(min_x, min_y, max_x, min_y)
        _, _, cell_height = geod.inv(min_x, min_y, min_x, max_y)

        features.append({
            "type": "Feature",
            "geometry": mapping(seed_cell_polygon),
            "properties": {
                    "rhealpix": seed_cell_id,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "cell_width": cell_width,
                    "cell_height": cell_height,
                    "cell_area": cell_area
                    },
        })
        return {
            "type": "FeatureCollection",
            "features": features,
        }

    else:
        # Initialize sets and queue
        covered_cells = set()  # Cells that have been processed (by their unique ID)
        queue = [seed_cell]  # Queue for BFS exploration
        while queue:
            current_cell = queue.pop()
            current_cell_id = str(current_cell)  # Unique identifier for the current cell

            if current_cell_id in covered_cells:
                continue

            # Add current cell to the covered set
            covered_cells.add(current_cell_id)

            # Convert current cell to polygon
            cell_polygon = cell_to_polygon(current_cell)

            # Skip cells that do not intersect the bounding box
            if not cell_polygon.intersects(bbox_polygon):
                continue

            # Get neighbors and add to queue
            neighbors = current_cell.neighbors(plane=False)
            for _, neighbor in neighbors.items():
                neighbor_id = str(neighbor)  # Unique identifier for the neighbor
                if neighbor_id not in covered_cells:
                    queue.append(neighbor)

        for cell_id in covered_cells:
            rhealpix_uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
            cell = rhealpix_dggs.cell(rhealpix_uids)    
            cell_polygon = cell_to_polygon(cell)
            min_x, min_y, max_x, max_y = cell_polygon.bounds           
            center_lon = cell_polygon.centroid.x
            center_lat = cell_polygon.centroid.y
            cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters                
            _, _, cell_width = geod.inv(min_x, min_y, max_x, min_y)
            _, _, cell_height = geod.inv(min_x, min_y, min_x, max_y)
            if cell_polygon.intersects(bbox_polygon):
                features.append({
                    "type": "Feature",
                    "geometry": mapping(cell_polygon),
                    "properties": {
                                    "rhealpix": cell_id,
                                    "center_lat": center_lat,
                                    "center_lon": center_lon,
                                    "cell_width": cell_width,
                                    "cell_height": cell_height,
                                    "cell_area": cell_area
                                    },
                })

        return {
            "type": "FeatureCollection",
            "features": features,
        }

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate RHEALPix grid within a bounding box and save as a GeoJSON.")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..15] of the grid")
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )
    args = parser.parse_args()

    # Initialize RHEALPix DGGS
    rhealpix_dggs = RHEALPixDGGS()
    resolution = args.resolution
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    if bbox == [-180, -90, 180, 90]:
        # Calculate the number of cells at the given resolution
        num_cells = rhealpix_dggs.num_cells(resolution)
        max_cells = 1_000_000
        if num_cells > max_cells:
            print(f"The selected resolution will generate {num_cells:,} cells, which exceeds the limit of {max_cells:,}.")
            print("Please select a smaller resolution and try again.")
            return

        # Generate grid within the bounding box
        geojson_features = generate_grid(rhealpix_dggs, resolution)

        # Define the GeoJSON file path
        geojson_path = f"rhealpix_grid_{resolution}.geojson"
        with open(geojson_path, 'w') as f:
            json.dump(geojson_features, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
    
    else:
        if resolution < 1 or resolution > 15:
            print(f"Please select a resolution in [1..15] range and try again ")
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
