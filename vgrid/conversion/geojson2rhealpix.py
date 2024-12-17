import argparse
import json
from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from shapely.geometry import Polygon, box, Point, LineString, mapping
from pyproj import Geod
import os
from vgrid.generator.rhealpixgrid import fix_rhealpix_antimeridian_cells

# Function to convert cell vertices to a Shapely Polygon
def cell_to_polygon(cell):
    vertices = [tuple(my_round(coord, 14) for coord in vertex) for vertex in cell.vertices(plane=False)]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    vertices = fix_rhealpix_antimeridian_cells(vertices)
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
    center_lat = round(seed_cell_polygon.centroid.y,7)
    center_lon = round(seed_cell_polygon.centroid.x,7)
    cell_area = round(abs(geod.geometry_area_perimeter(seed_cell_polygon)[0]),2)
    cell_perimeter = abs(geod.geometry_area_perimeter(seed_cell_polygon)[1])
    avg_edge_len = round(cell_perimeter/4,2)
    if seed_cell.ellipsoidal_shape() == 'dart':
        avg_edge_len = round(cell_perimeter/3,2)  
    
    features.append({
                "type": "Feature",
                "geometry": mapping(seed_cell_polygon),
                "properties": {
                        "rhealpix": seed_cell_id,
                        "center_lat": center_lat,
                        "center_lon": center_lon,
                        "cell_area": cell_area,
                        "avg_edge_len": avg_edge_len,
                        "resolution": resolution
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
            center_lat = round(seed_cell_polygon.centroid.y,7)
            center_lon = round(seed_cell_polygon.centroid.x,7)
            cell_area = abs(geod.geometry_area_perimeter(seed_cell_polygon)[0])  # Area in square meters                
            cell_perimeter = abs(geod.geometry_area_perimeter(seed_cell_polygon)[1])
            avg_edge_len = round(cell_perimeter/4,2)
            if seed_cell.ellipsoidal_shape() == 'dart':
                avg_edge_len = round(cell_perimeter/3,2)  
            
            features.append({
                "type": "Feature",
                "geometry": mapping(seed_cell_polygon),
                "properties": {
                        "rhealpix": seed_cell_id,
                        "center_lat": center_lat,
                        "center_lon": center_lon,
                        "cell_area": cell_area,
                        "avg_edge_len": avg_edge_len,
                        "resolution": resolution
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
                center_lat = round(cell_polygon.centroid.y,7)
                center_lon = round(cell_polygon.centroid.x,7)
                cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters                
                cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
                avg_edge_len = round(cell_perimeter/4,2)
                if seed_cell.ellipsoidal_shape() == 'dart':
                    avg_edge_len = round(cell_perimeter/3,2)  
               
                if cell_polygon.intersects(polyline):
                    features.append({
                        "type": "Feature",
                        "geometry": mapping(cell_polygon),
                        "properties": {
                                "rhealpix": cell_id,
                                "center_lat": center_lat,
                                "center_lon": center_lon,
                                "cell_area": cell_area,
                                "avg_edge_len": avg_edge_len,
                                "resolution": resolution
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
            center_lat = round(seed_cell_polygon.centroid.y,7)
            center_lon = round(seed_cell_polygon.centroid.x,7)
            cell_area = round(abs(geod.geometry_area_perimeter(seed_cell_polygon)[0]),2)  # Area in square meters
            cell_perimeter = abs(geod.geometry_area_perimeter(seed_cell_polygon)[1])
            avg_edge_len = round(cell_perimeter/4,2)
            if seed_cell.ellipsoidal_shape() == 'dart':
                avg_edge_len = round(cell_perimeter/3,2)  
            features.append({
                "type": "Feature",
                "geometry": mapping(seed_cell_polygon),
                "properties": {
                        "rhealpix": seed_cell_id,
                        "center_lat": center_lat,
                        "center_lon": center_lon,
                        "cell_area": cell_area,
                        "avg_edge_len": avg_edge_len,
                        "resolution": resolution
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
                center_lat = round(cell_polygon.centroid.y,7)
                center_lon = round(cell_polygon.centroid.x,7)
                cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters
                cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
                avg_edge_len = round(cell_perimeter/4,2)
                if seed_cell.ellipsoidal_shape() == 'dart':
                    avg_edge_len = round(cell_perimeter/3,2)  
                if cell_polygon.intersects(polygon):
                    features.append({
                        "type": "Feature",
                        "geometry": mapping(cell_polygon),
                        "properties": {
                            "rhealpix": cell_id,
                            "center_lat": center_lat,
                            "center_lon": center_lon,
                            "cell_area": cell_area,
                            "avg_edge_len": avg_edge_len
                        },
                    })

    return {
        "type": "FeatureCollection",
        "features": features,
    }

# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Generate RHEALPix grid for shapes in GeoJSON format")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [0..15]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON string with Point, Polyline or Polygon"
    )
    args = parser.parse_args()
    geojson = args.geojson
     # Initialize RHEALPix DGGS
    rhealpix_dggs = RHEALPixDGGS()
    resolution = args.resolution
    
    if resolution < 0 or resolution > 15:
        print(f"Please select a resolution in [0..15] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = []

    for feature in geojson_data['features']:      
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(rhealpix_dggs, resolution, point)
                geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)  # Create Point for each coordinate set
                    point_features = point_to_grid(rhealpix_dggs, resolution, point)
                    geojson_features.extend(point_features['features'])
        
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                # Directly process LineString geometry
                polyline = LineString(coordinates)
                polyline_features = polyline_to_grid(rhealpix_dggs, resolution, polyline)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                # Iterate through each line in MultiLineString geometry
                for line_coords in coordinates:
                    polyline = LineString(line_coords)  # Use each part's coordinates
                    polyline_features = polyline_to_grid(rhealpix_dggs, resolution, polyline)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                # Create Polygon with exterior and interior rings
                exterior_ring = coordinates[0]  # The first coordinate set is the exterior ring
                interior_rings = coordinates[1:]  # Remaining coordinate sets are interior rings (holes)
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = polygon_to_grid(rhealpix_dggs, resolution, polygon)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = polygon_to_grid(rhealpix_dggs, resolution, polygon)
                    geojson_features.extend(polygon_features['features'])

                    
    # Save the results to GeoJSON
    geojson_path = f"geojson2rhealpix_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
