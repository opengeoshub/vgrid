from vgrid.utils import s2
from shapely.geometry import Point, LineString, Polygon
import argparse
import json
from tqdm import tqdm
import os
from vgrid.generator.s2grid import s2_cell_to_polygon
from vgrid.generator.settings import geodesic_dggs_to_feature

def point_to_grid(resolution, point, feature_properties):    
    s2_features = []
    # Convert point to the seed cell
    latitude = point.y
    longitude = point.x
    lat_lng = s2.LatLng.from_degrees(latitude, longitude)
    cell_id_max_res = s2.CellId.from_lat_lng(lat_lng)
    cell_id = cell_id_max_res.parent(resolution)
    s2_cell = s2.Cell(cell_id)
    cell_token = s2.CellId.to_token(s2_cell.id())    
    if s2_cell:
        cell_polygon = s2_cell_to_polygon(cell_id) # Fix antimeridian
        resolution = cell_id.level()
        num_edges = 4
        s2_feature = geodesic_dggs_to_feature("s2",cell_token,resolution,cell_polygon,num_edges)   
        s2_feature["properties"].update(feature_properties)
        s2_features.append(s2_feature)

    return {
        "type": "FeatureCollection",
        "features": s2_features
    }
        
def poly_to_grid(resolution, geometry, feature_properties, compact=None):
    s2_features = []
    
    if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon':
        polys = [geometry]
    elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'MultiPolygon':
        polys = list(geometry)

    for poly in polys:    
        min_lng, min_lat, max_lng, max_lat = poly.bounds
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
        cell_ids = covering  
        if compact:
            covering = s2.CellUnion(covering)
            covering.normalize()
            cell_ids = covering.cell_ids()  
            
        for cell_id in cell_ids:
            cell_polygon = s2_cell_to_polygon(cell_id)
          
            if cell_polygon.intersects(poly):
                cell_token = s2.CellId.to_token(cell_id)  
                cell_resolution = cell_id.level()
                num_edges = 4
                s2_feature = geodesic_dggs_to_feature("s2",cell_token,cell_resolution,cell_polygon,num_edges)   
                s2_feature["properties"].update(feature_properties)    
                s2_features.append(s2_feature)
                            
    return {
        "type": "FeatureCollection",
        "features": s2_features,
    }

def geojson2s2(geojson_data, resolution, compact=False):
    """
    Convert GeoJSON data to S2 DGGS format.
    
    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): S2 resolution level [0..30]
        compact (bool): Enable S2 compact mode for polygons
        
    Returns:
        dict: GeoJSON FeatureCollection with S2 cells
    """
    if resolution < 0 or resolution > 30:
        raise ValueError("Resolution must be in range [0..30]")
        
    geojson_features = []

    for feature in tqdm(geojson_data['features'], desc="Processing features"):
        feature_properties = feature['properties']
        if feature['geometry']['type'] in ['Point', 'MultiPoint']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'Point':
                point = Point(coordinates)
                point_features = point_to_grid(resolution, point, feature_properties)
                geojson_features.extend(point_features['features'])

            elif feature['geometry']['type'] == 'MultiPoint':
                for point_coords in coordinates:
                    point = Point(point_coords)
                    point_features = point_to_grid(resolution, point, feature_properties)
                    geojson_features.extend(point_features['features'])
        
        elif feature['geometry']['type'] in ['LineString', 'MultiLineString']:
            coordinates = feature['geometry']['coordinates']
            if feature['geometry']['type'] == 'LineString':
                polyline = LineString(coordinates)
                polyline_features = poly_to_grid(resolution, polyline, feature_properties)
                geojson_features.extend(polyline_features['features'])

            elif feature['geometry']['type'] == 'MultiLineString':
                for line_coords in coordinates:
                    polyline = LineString(line_coords)
                    polyline_features = poly_to_grid(resolution, polyline, feature_properties)
                    geojson_features.extend(polyline_features['features'])
            
        elif feature['geometry']['type'] in ['Polygon', 'MultiPolygon']:
            coordinates = feature['geometry']['coordinates']

            if feature['geometry']['type'] == 'Polygon':
                exterior_ring = coordinates[0]
                interior_rings = coordinates[1:]
                polygon = Polygon(exterior_ring, interior_rings)
                polygon_features = poly_to_grid(resolution, polygon, feature_properties, compact)
                geojson_features.extend(polygon_features['features'])

            elif feature['geometry']['type'] == 'MultiPolygon':
                for sub_polygon_coords in coordinates:
                    exterior_ring = sub_polygon_coords[0]
                    interior_rings = sub_polygon_coords[1:]
                    polygon = Polygon(exterior_ring, interior_rings)
                    polygon_features = poly_to_grid(resolution, polygon, feature_properties, compact)
                    geojson_features.extend(polygon_features['features'])

    return {
        "type": "FeatureCollection",
        "features": geojson_features
    }

def geojson2s2_cli():
    """Command line interface for converting GeoJSON to S2 DGGS format."""
    parser = argparse.ArgumentParser(description="Convert GeoJSON to S2 DGGS")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [0..30]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, help="GeoJSON file path (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable S2 compact mode - for polygon only")

    args = parser.parse_args()
    geojson = args.geojson
    resolution = args.resolution
    compact = args.compact  
    
    if not os.path.exists(geojson):
        print(f"Error: The file {geojson} does not exist.")
        return

    with open(geojson, 'r', encoding='utf-8') as f:
        geojson_data = json.load(f)
    
    try:
        result = geojson2s2(geojson_data, resolution, compact)
        
        # Save the results to GeoJSON
        geojson_name = os.path.splitext(os.path.basename(geojson))[0]
        geojson_path = f"{geojson_name}2s2_{resolution}.geojson"
        if compact:
            geojson_path = f"{geojson_name}2s2_{resolution}_compacted.geojson"
        
        with open(geojson_path, 'w') as f:
            json.dump(result, f, indent=2)

        print(f"GeoJSON saved as {geojson_path}")
    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    geojson2s2_cli()