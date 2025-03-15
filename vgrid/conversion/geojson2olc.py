from vgrid.utils import olc
from shapely.geometry import Point, LineString, Polygon, mapping, box
import argparse
import json
from tqdm import tqdm
import os
from pyproj import Geod
geod = Geod(ellps="WGS84")
from vgrid.generator.olcgrid import generate_grid,refine_cell

# Function to generate grid for Point
def point_to_grid(resolution, point):    
    features = []
    olc_cellid = olc.encode(point.y, point.x, resolution)
    coord = olc.decode(olc_cellid)    
    if coord:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
        max_lat, max_lon = coord.latitudeHi, coord.longitudeHi

        center_lat = round(coord.latitudeCenter,7)
        center_lon = round(coord.longitudeCenter,7)
        resolution = coord.codeLength 

        # Define the polygon based on the bounding box
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        
        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),2)  # Area in square meters     
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])       
        avg_edge_len = round(cell_perimeter / 4,2)
        features.append({
                "type": "Feature",
                "geometry": mapping(cell_polygon),          
                "properties": {
                    "olc": olc_cellid,  
                    "resolution": resolution,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "avg_edge_len": avg_edge_len,
                    "cell_area": cell_area                
                }
            })

    geojson_features = {
        'type': 'FeatureCollection',
        'features': features
    }

    return geojson_features

# Function to generate grid for Polyline
def poly_to_grid(resolution, geometry):
    # Extract points from polyline
    if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon':
        # Handle single Polygon as before
        polys = [geometry]
    elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'Multipolygon':
        # Handle MultiPolygon: process each polygon separately
        polys = list(geometry)

    for poly in polys:  
        base_resolution = 2
        base_cells = generate_grid(base_resolution)

        # Step 2: Identify seed cells that intersect with the bounding box
        seed_cells = []
        for base_cell in base_cells["features"]:
            base_cell_poly = Polygon(base_cell["geometry"]["coordinates"][0])
            if poly.intersects(base_cell_poly):
                seed_cells.append(base_cell)

        refined_features = []

        # Step 3: Iterate over seed cells and refine to the output resolution
        for seed_cell in seed_cells:
            seed_cell_poly = Polygon(seed_cell["geometry"]["coordinates"][0])

            if seed_cell_poly.contains(poly) and resolution == base_resolution:
                # Append the seed cell directly if fully contained and resolution matches
                refined_features.append(seed_cell)
            else:
                # Refine the seed cell to the output resolution and add it to the output
                refined_features.extend(
                    refine_cell(seed_cell_poly.bounds, base_resolution, resolution, poly)
                )

            resolution_features = [
                feature for feature in refined_features if feature["properties"]["resolution"] == resolution
            ]

            final_features = []
            seen_olc_codes = set()  # Reset the set for final feature filtering

            for feature in resolution_features:
                olc_code = feature["properties"]["olc"]
                if olc_code not in seen_olc_codes:  # Check if OLC code is already in the set
                    final_features.append(feature)
                    seen_olc_codes.add(olc_code)

            return {
                "type": "FeatureCollection",
                "features": final_features
            }
    

# Main function to handle different GeoJSON shapes
def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to Open_Eaggr OLC grid")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution of the grid [2, 4, 6, 8, 10..15]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="Point, Polyline or Polygon in GeoJSOn format"
    )
    args = parser.parse_args()
    geojson = args.geojson
     # Initialize h3 DGGS
    resolution = args.resolution
    
    if resolution not in [2, 4, 6, 8, 10, 11, 12, 13, 14, 15]:
        print(f"Please select a resolution in [2, 4, 6, 8, 10..15] and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing GeoJSON features"):
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(resolution, point)
                geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)  # Create Point for each coordinate set
                    point_features = point_to_grid(resolution, point)
                    geojson_features.extend(point_features['features'])
        
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                # Directly process LineString geometry
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(resolution, polyline)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                # Iterate through each line in MultiLineString geometry
                for line_coords in coordinates:
                    polyline = LineString(line_coords)  # Use each part's coordinates
                    polyline_features = poly_to_grid(resolution, polyline)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                # Create Polygon with exterior and interior rings
                exterior_ring = coordinates[0]  # The first coordinate set is the exterior ring
                interior_rings = coordinates[1:]  # Remaining coordinate sets are interior rings (holes)
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(resolution, polygon)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(resolution, polygon)
                    geojson_features.extend(polygon_features['features'])

    # Save the results to GeoJSON
    geojson_path = f"geojson2olc_{resolution}.geojson"
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
