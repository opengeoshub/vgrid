import argparse
import json
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from shapely.geometry import Polygon, box, mapping
from tqdm import tqdm
from vgrid.generator.settings import max_cells, geodesic_dggs_to_feature

from pyproj import Geod
geod = Geod(ellps="WGS84")
rhealpix_dggs = RHEALPixDGGS()

# Function to filter cells crossing the antimeridian
def fix_rhealpix_antimeridian_cells(boundary, threshold=-128):
    if any(lon < threshold for lon, _ in boundary):
        return [(lon - 360 if lon > 0 else lon, lat) for lon, lat in boundary]
    return boundary

# Function to convert cell vertices to a Shapely Polygon
def rhealpix_cell_to_polygon(cell):
    vertices = [tuple(my_round(coord, 14) for coord in vertex) for vertex in cell.vertices(plane=False)]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    vertices = fix_rhealpix_antimeridian_cells(vertices)
    return Polygon(vertices)

def generate_grid(resolution):
    rhealpix_features = []
    total_cells = rhealpix_dggs.num_cells(resolution)
    rhealpix_grid = rhealpix_dggs.grid(resolution)
            
    with tqdm(total=total_cells, desc="Processing cells", unit=" cells") as pbar:
        for rhealpix_cell in rhealpix_grid:
            cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
            rhealpix_id = str(rhealpix_cell)
            num_edges = 4
            if rhealpix_cell.ellipsoidal_shape() == 'dart':
                num_edges = 3
            rhealpix_feature = geodesic_dggs_to_feature('rhealpix', rhealpix_id, resolution, cell_polygon, num_edges)
            rhealpix_features.append(rhealpix_feature)
            pbar.update(1)
        
        return {
            "type": "FeatureCollection",
            "features": rhealpix_features,
        }
      
def generate_grid_within_bbox(resolution, bbox):    
    bbox_polygon = box(*bbox)  # Create a bounding box polygon
    bbox_center_lon = bbox_polygon.centroid.x
    bbox_center_lat = bbox_polygon.centroid.y
    seed_point = (bbox_center_lon, bbox_center_lat)

    rhealpix_features = []
    seed_cell = rhealpix_dggs.cell_from_point(resolution, seed_point, plane=False)
    seed_cell_id = str(seed_cell)  # Unique identifier for the current cell
    seed_cell_polygon = rhealpix_cell_to_polygon(seed_cell)

    if seed_cell_polygon.contains(bbox_polygon):             
        num_edges = 4
        if seed_cell.ellipsoidal_shape() == 'dart':
            num_edges = 3

        rhealpix_feature = geodesic_dggs_to_feature('rhealpix', seed_cell_id,resolution, seed_cell_polygon, num_edges)
        rhealpix_features.append(rhealpix_feature)
        return {
            "type": "FeatureCollection",
            "features": rhealpix_features,
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
            cell_polygon = rhealpix_cell_to_polygon(current_cell)

            # Skip cells that do not intersect the bounding box
            if not cell_polygon.intersects(bbox_polygon):
                continue

            # Get neighbors and add to queue
            neighbors = current_cell.neighbors(plane=False)
            for _, neighbor in neighbors.items():
                neighbor_id = str(neighbor)  # Unique identifier for the neighbor
                if neighbor_id not in covered_cells:
                    queue.append(neighbor)

        for cell_id in tqdm(covered_cells, desc="Processing Cells"):
            rhealpix_uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
            cell = rhealpix_dggs.cell(rhealpix_uids)    
            cell_polygon = rhealpix_cell_to_polygon(cell)          
            if cell_polygon.intersects(bbox_polygon):             
                num_edges = 4
                if seed_cell.ellipsoidal_shape() == 'dart':
                    num_edges = 3
                rhealpix_feature = geodesic_dggs_to_feature('rhealpix', cell_id,resolution, cell_polygon, num_edges)
                rhealpix_features.append(rhealpix_feature)
                
        return {
            "type": "FeatureCollection",
            "features": rhealpix_features,
        }


def generate_grid_within_bbox_new(resolution, bbox):    
    bbox_polygon = box(*bbox)  # Create a bounding box polygon
    bbox_center_lon = bbox_polygon.centroid.x
    bbox_center_lat = bbox_polygon.centroid.y
    seed_point = (bbox_center_lon, bbox_center_lat)

    rhealpix_features = []
    seed_cell = rhealpix_dggs.cell_from_point(resolution, seed_point, plane=False)
    seed_cell_id = str(seed_cell)  # Unique identifier for the current cell
    seed_cell_polygon = rhealpix_cell_to_polygon(seed_cell)
    num_edges = 4
    
    if seed_cell_polygon.contains(bbox_polygon):
        rhealpix_feature = geodesic_dggs_to_feature('rhealpix',seed_cell_id,resolution,seed_cell_polygon,num_edges)
        rhealpix_features.append(rhealpix_feature)
        return {
            "type": "FeatureCollection",
            "features": rhealpix_features
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
            cell_polygon = rhealpix_cell_to_polygon(current_cell)

            # Skip cells that do not intersect the bounding box
            if not cell_polygon.intersects(bbox_polygon):
                continue

            # Get neighbors and add to queue
            neighbors = current_cell.neighbors(plane=False)
            for _, neighbor in neighbors.items():
                neighbor_id = str(neighbor)  # Unique identifier for the neighbor
                if neighbor_id not in covered_cells:
                    queue.append(neighbor)

        for cell_id in tqdm(covered_cells, desc="Processing Cells"):
            rhealpix_uids = (cell_id[0],) + tuple(map(int, cell_id[1:]))
            cell = rhealpix_dggs.cell(rhealpix_uids)    
            cell_polygon = rhealpix_cell_to_polygon(cell)           
            if cell_polygon.intersects(bbox_polygon):
                rhealpix_feature = geodesic_dggs_to_feature('rhealpix',seed_cell_id,resolution,seed_cell_polygon,num_edges)
                rhealpix_features.append(rhealpix_feature)
        
        return {
            "type": "FeatureCollection",
            "features": rhealpix_features
        }

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Generate RHEALPix grid.")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..15] of the grid")
    parser.add_argument(
        '-b', '--bbox', type=float, nargs=4, 
        help="Bounding box in the format: min_lon min_lat max_lon max_lat (default is the whole world)"
    )
    args = parser.parse_args()

    # Initialize RHEALPix DGGS
    
    resolution = args.resolution
    if resolution < 0 or resolution > 15:
            print(f"Please select a resolution in [0..15] range and try again ")
            return
    bbox = args.bbox if args.bbox else [-180, -90, 180, 90]
    if bbox == [-180, -90, 180, 90]:
        # Calculate the number of cells at the given resolution
        num_cells = rhealpix_dggs.num_cells(resolution)
        print(f"Resolution {resolution} will generate {num_cells} cells ")

        if num_cells > max_cells:
            print(
                f"which exceeds the limit of {max_cells}."
            )
            print("Please select a smaller resolution and try again.")
            return
        # Generate grid within the bounding box
        geojson_features = generate_grid(resolution)      
    else:
        # Generate grid within the bounding box
        geojson_features = generate_grid_within_bbox(resolution, bbox)
        # Define the GeoJSON file path
    
    geojson_path = f"rhealpix_grid_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump(geojson_features, f, indent=2)
    print(f"GeoJSON saved as {geojson_path}")
         
if __name__ == "__main__":
    main()
