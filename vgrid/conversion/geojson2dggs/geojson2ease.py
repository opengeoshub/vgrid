import argparse, json, os
from tqdm import tqdm
from shapely.geometry import Polygon, box, Point, LineString
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.utils.easedggs.constants import grid_spec, ease_crs, geo_crs, levels_specs
from vgrid.utils.easedggs.dggs.grid_addressing import grid_ids_to_geos, geos_to_grid_ids, geo_polygon_to_grid_ids
from vgrid.conversion.dggscompact import ease_compact

# Function to generate grid for Point
def point_to_grid(resolution, point, feature_properties):
    ease_features = []
    latitude = point.y
    longitude = point.x
    ease_cell = geos_to_grid_ids([(longitude,latitude)],level = resolution)
    ease_id = ease_cell['result']['data'][0]

    level = int(ease_id[1])  # Get the level (e.g., 'L0' -> 0)
    # Get level specs
    level_spec = levels_specs[level]
    n_row = level_spec["n_row"]
    n_col = level_spec["n_col"]
        
    geo = grid_ids_to_geos([ease_id])
    center_lon, center_lat = geo['result']['data'][0] 

    cell_min_lat = center_lat - (180 / (2 * n_row))
    cell_max_lat = center_lat + (180 / (2 * n_row))
    cell_min_lon = center_lon - (360 / (2 * n_col))
    cell_max_lon = center_lon + (360 / (2 * n_col))

    cell_polygon = Polygon([
        [cell_min_lon, cell_min_lat],
        [cell_max_lon, cell_min_lat],
        [cell_max_lon, cell_max_lat],
        [cell_min_lon, cell_max_lat],
        [cell_min_lon, cell_min_lat]
     ])

    if cell_polygon:   
        resolution = level
        num_edges = 4
        ease_feature = geodesic_dggs_to_feature("ease",ease_id,resolution,cell_polygon,num_edges)  
        ease_feature["properties"].update(feature_properties) 
        ease_features.append(ease_feature)

    return {
        "type": "FeatureCollection",
        "features": ease_features
    }

def poly_to_grid(resolution, geometry, feature_properties,compact=None):    
    ease_features = []
    if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon':
        polys = [geometry]
    elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'MultiPolygon':
        polys = list(geometry)

    for poly in polys:
        poly_bbox = box(*poly.bounds)
        # Get all grid cells within the bounding box
        polygon_bbox_wkt = poly_bbox.wkt
        cells_bbox = geo_polygon_to_grid_ids(polygon_bbox_wkt, resolution, geo_crs, ease_crs, levels_specs, return_centroids = True, wkt_geom=True)
        cells = cells_bbox['result']['data']   
        
        if compact:
            cells = ease_compact(cells)
             
        for cell in cells:
            cell_resolution =  int(cell[1])  # Get the level (e.g., 'L0' -> 0)
            level_spec = levels_specs[cell_resolution]
            n_row = level_spec["n_row"]
            n_col = level_spec["n_col"]

            geo = grid_ids_to_geos([cell])
            center_lon, center_lat = geo['result']['data'][0] 

            cell_min_lat = center_lat - (180 / (2 * n_row))
            cell_max_lat = center_lat + (180 / (2 * n_row))
            cell_min_lon = center_lon - (360 / (2 * n_col))
            cell_max_lon = center_lon + (360 / (2 * n_col))

            cell_polygon = Polygon([
                [cell_min_lon, cell_min_lat],
                [cell_max_lon, cell_min_lat],
                [cell_max_lon, cell_max_lat],
                [cell_min_lon, cell_max_lat],
                [cell_min_lon, cell_min_lat]
            ])
            if cell_polygon.intersects(poly):
                num_edges = 4
                ease_feature = geodesic_dggs_to_feature('ease', str(cell), cell_resolution, cell_polygon, num_edges)
                ease_feature["properties"].update(feature_properties)
                ease_features.append(ease_feature)            
   
    return {
        "type": "FeatureCollection",
        "features": ease_features
    }

        
def main():
    parser = argparse.ArgumentParser(description="Convert GeoJSON to EASE-DGGS")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..6]")
    # actual resolution range: [0..40]

    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON file path (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable EASE compact mode")

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    compact = args.compact  

    if resolution < 0 or resolution > 6:
        print(f"Please select a resolution in [0..6] range and try again ")
        return
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing features"):
        feature_properties = feature['properties'] 
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(resolution, point,feature_properties)
                geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)  # Create Point for each coordinate set
                    point_features = point_to_grid(resolution, point,feature_properties)
                    geojson_features.extend(point_features['features'])
        
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                # Directly process LineString geometry
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(resolution, polyline,feature_properties)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                # Iterate through each line in MultiLineString geometry
                for line_coords in coordinates:
                    polyline = LineString(line_coords)  # Use each part's coordinates
                    polyline_features = poly_to_grid(resolution, polyline,feature_properties)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                # Create Polygon with exterior and interior rings
                exterior_ring = coordinates[0]  # The first coordinate set is the exterior ring
                interior_rings = coordinates[1:]  # Remaining coordinate sets are interior rings (holes)
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(resolution, polygon,feature_properties,compact)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                # Handle each sub-polygon in MultiPolygon geometry
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]  # The first coordinate set is the exterior ring
                    interior_rings = sub_polygon_coords[1:]  # Remaining coordinate sets are interior rings (holes)
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(resolution, polygon,feature_properties,compact)
                    geojson_features.extend(polygon_features['features'])

                    
    geojson_name = os.path.splitext(os.path.basename(geojson))[0]
    geojson_path = f"{geojson_name}2ease_{resolution}.geojson"
    if compact:        
        geojson_path = f"{geojson_name}2ease_{resolution}_compacted.geojson"
        
    with open(geojson_path, 'w') as f:
        json.dump({"type": "FeatureCollection", "features": geojson_features}, f, indent=2)

    print(f"GeoJSON saved as {geojson_path}")


if __name__ == "__main__":
    main()
