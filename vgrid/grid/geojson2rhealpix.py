import argparse
import json
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from shapely.geometry import Polygon, box, Point, LineString, mapping
from tqdm import tqdm
from pyproj import Geod
import os
from vgrid.grid.rhealpixgrid import generate_grid_within_bbox

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


# Function to generate grid for Point
def point_to_grid(rhealpix_dggs, resolution, point):
    features = []
    # Convert point to the seed cell
    seed_cell = rhealpix_dggs.cell_from_point(resolution, (point.x, point.y), plane=False)
    seed_cell_id = str(seed_cell)  # Unique identifier for the current cell
    seed_cell_polygon = cell_to_polygon(seed_cell)
    geod = Geod(ellps="WGS84")
    
    # Get the bounds and area of the cell
    min_x, min_y, max_x, max_y = seed_cell_polygon.bounds
    center_lon = seed_cell_polygon.centroid.x
    center_lat = seed_cell_polygon.centroid.y
    cell_area = abs(geod.geometry_area_perimeter(seed_cell_polygon)[0])
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


# Function to generate grid for Polyline
def polyline_to_grid(rhealpix_dggs, resolution, geometry):
    features = []
    geod = Geod(ellps="WGS84")
    # Extract points from polyline
    if geometry.geom_type == 'LineString':
        # Handle single Polygon as before
        polylines = [geometry]
    elif geometry.geom_type == 'MultiLineString':
        # Handle MultiPolygon: process each polygon separately
        polylines = list(geometry)

    for polyline in polylines:
        minx, miny, maxx, maxy = polyline.bounds
        # Create a bounding box polygon
        bbox_polygon = box(minx, miny, maxx, maxy)

        bbox_center_lon = bbox_polygon.centroid.x
        bbox_center_lat = bbox_polygon.centroid.y
        seed_point = (bbox_center_lon, bbox_center_lat)

        seed_cell = rhealpix_dggs.cell_from_point(resolution, seed_point, plane=False)
        seed_cell_id = str(seed_cell)  # Unique identifier for the current cell
        seed_cell_polygon = cell_to_polygon(seed_cell)

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
                if cell_polygon.intersects(polyline):
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
        
# Function to generate grid for Polygon
def polygon_to_grid(rhealpix_dggs, resolution, geometry):
    features = []
    geod = Geod(ellps="WGS84")
    
    if geometry.geom_type == 'Polygon':
        # Handle single Polygon as before
        polygons = [geometry]
    elif geometry.geom_type == 'MultiPolygon':
        # Handle MultiPolygon: process each polygon separately
        polygons = list(geometry)

    for polygon in polygons:
        # Processing each polygon (either single Polygon or MultiPolygon)
        polygon_center_lon = polygon.representative_point().x
        polygon_center_lat = polygon.representative_point().y
        
        seed_point = (polygon_center_lon, polygon_center_lat)

        seed_cell = rhealpix_dggs.cell_from_point(resolution, seed_point, plane=False)
        seed_cell_id = str(seed_cell)  # Unique identifier for the current cell
        seed_cell_polygon = cell_to_polygon(seed_cell)

        if seed_cell_polygon.contains(polygon):
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
        else:
            # Process grid for non-contained cells
            covered_cells = set()  # Cells that have been processed
            queue = [seed_cell]  # Queue for BFS exploration
            while queue:
                current_cell = queue.pop()
                current_cell_id = str(current_cell)

                if current_cell_id in covered_cells:
                    continue

                covered_cells.add(current_cell_id)
                cell_polygon = cell_to_polygon(current_cell)

                if not cell_polygon.intersects(polygon):
                    continue

                neighbors = current_cell.neighbors(plane=False)
                for _, neighbor in neighbors.items():
                    neighbor_id = str(neighbor)
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
                if cell_polygon.intersects(polygon):
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

# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Generate RHEALPix grid for shapes in GeoJSON format")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON string with Point, Polyline or Polygon"
    )
    args = parser.parse_args()
    geojson = args.geojson
     # Initialize RHEALPix DGGS
    rhealpix_dggs = RHEALPixDGGS()
    resolution = args.resolution
    
    if resolution < 1 or resolution > 15:
        print(f"Please select a resolution in [1..15] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = []

    # Iterate over each feature in the GeoJSON
    for feature in geojson_data['features']:
        if feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Polygon':
                coordinates = feature['geometry']['coordinates'][0]  # Only access the first set of coordinates
                polygon = Polygon(coordinates)
                polygon_features = polygon_to_grid(rhealpix_dggs, resolution, polygon)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                for coords in feature['geometry']['coordinates']:
                    polygon = Polygon(coords[0])  # Handle the first set of coordinates for each multi-polygon
                    polygon_features = polygon_to_grid(rhealpix_dggs, resolution, polygon)
                    geojson_features.extend(polygon_features['features'])

        # Handle other geometry types (Point, Polyline) if needed
        elif feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                coordinates = feature['geometry']['coordinates'] 
                point = Point(coordinates)
                point_features = point_to_grid(rhealpix_dggs, resolution, point)
                geojson_features.extend(point_features['features'])
    
        # Process Polyline feature if needed
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                coordinates = feature['geometry']['coordinates'] 
                polyline = LineString(coordinates)
                polyline_features = polyline_to_grid(rhealpix_dggs, resolution, polyline)
                geojson_features.extend(polyline_features['features'])

            # elif feature['geometry']['type'] == 'MultiLineString':
            #     for coords in feature['geometry']['coordinates']:
            #         polyline = LineString(coords[0])  # Handle the first set of coordinates for each multi-polygon
            #         polyline_features = polyline_to_grid(rhealpix_dggs, resolution, polyline)
            #         geojson_features.extend(polyline_features['features'])

    # Save the results to GeoJSON
    geojson_path = f"geojson2rhealpix_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
