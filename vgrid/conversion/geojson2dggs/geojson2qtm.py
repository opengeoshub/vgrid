from shapely.geometry import shape, Point, LineString, Polygon
import argparse
import os, json
from tqdm import tqdm
from vgrid.utils import qtm
from vgrid.generator.settings import geodesic_dggs_to_feature
from vgrid.conversion.dggscompact import qtmcompact
import requests
from urllib.parse import urlparse

p90_n180, p90_n90, p90_p0, p90_p90, p90_p180 = (90.0, -180.0), (90.0, -90.0), (90.0, 0.0), (90.0, 90.0), (90.0, 180.0)
p0_n180, p0_n90, p0_p0, p0_p90, p0_p180 = (0.0, -180.0), (0.0, -90.0), (0.0, 0.0), (0.0, 90.0), (0.0, 180.0)
n90_n180, n90_n90, n90_p0, n90_p90, n90_p180 = (-90.0, -180.0), (-90.0, -90.0), (-90.0, 0.0), (-90.0, 90.0), (-90.0, 180.0)

# Function to generate grid for Point
def point_to_grid(resolution, point, feature_properties):
    qtm_features = []
    # Convert point to the seed cell
    latitude = point.y
    longitude = point.x
    qtm_id = qtm.latlon_to_qtm_id(latitude, longitude, resolution) 
    facet = qtm.qtm_id_to_facet(qtm_id)
    cell_polygon = qtm.constructGeometry(facet) 
    if cell_polygon:
        num_edges = 3
        qtm_feature = geodesic_dggs_to_feature("qtm",qtm_id,resolution,cell_polygon,num_edges)   
        qtm_feature["properties"].update(feature_properties)
        qtm_features.append(qtm_feature)

    return {
        "type": "FeatureCollection",
        "features": qtm_features
    }


def poly_to_grid(resolution, geometry, feature_properties,compact=None):    
    if geometry.geom_type == 'LineString' or geometry.geom_type == 'Polygon':
        polys = [geometry]
    elif geometry.geom_type == 'MultiLineString' or geometry.geom_type == 'MultiPolygon':
        polys = list(geometry)

    for poly in polys:
        levelFacets = {}
        QTMID = {}
        qtm_features = []    
        for lvl in range(resolution):
            levelFacets[lvl] = []
            QTMID[lvl] = []

            if lvl == 0:
                initial_facets = [
                    [p0_n180, p0_n90, p90_n90, p90_n180, p0_n180, True],
                    [p0_n90, p0_p0, p90_p0, p90_n90, p0_n90, True],
                    [p0_p0, p0_p90, p90_p90, p90_p0, p0_p0, True],
                    [p0_p90, p0_p180, p90_p180, p90_p90, p0_p90, True],
                    [n90_n180, n90_n90, p0_n90, p0_n180, n90_n180, False],
                    [n90_n90, n90_p0, p0_p0, p0_n90, n90_n90, False],
                    [n90_p0, n90_p90, p0_p90, p0_p0, n90_p0, False],
                    [n90_p90, n90_p180, p0_p180, p0_p90, n90_p90, False],
                ]

                for i, facet in enumerate(initial_facets):
                    QTMID[0].append(str(i + 1))
                    levelFacets[0].append(facet)
                    facet_geom = qtm.constructGeometry(facet)                  

                    if shape(facet_geom).intersects(poly) and resolution == 1 :
                        qtm_id = QTMID[0][i]
                        num_edges = 3
                        qtm_feature = geodesic_dggs_to_feature("qtm",qtm_id,resolution,facet_geom,num_edges)   
                        qtm_feature["properties"].update(feature_properties)
                        qtm_features.append(qtm_feature)

                        return {
                            "type": "FeatureCollection",
                            "features": qtm_features
                        }                       
            else:
                for i, pf in enumerate(levelFacets[lvl - 1]):
                    subdivided_facets = qtm.divideFacet(pf)
                    for j, subfacet in enumerate(subdivided_facets):
                        subfacet_geom = qtm.constructGeometry(subfacet)
                        if shape(subfacet_geom).intersects(poly):  # Only keep intersecting facets
                            new_id = QTMID[lvl - 1][i] + str(j)
                            QTMID[lvl].append(new_id)
                            levelFacets[lvl].append(subfacet)
                            if lvl == resolution - 1:  # Only store final resolution in GeoJSON
                                num_edges = 3
                                qtm_feature = geodesic_dggs_to_feature("qtm",new_id,resolution,subfacet_geom,num_edges)   
                                qtm_feature["properties"].update(feature_properties)
                                qtm_features.append(qtm_feature)
    
    qtm_geosjon = {
        "type": "FeatureCollection",
        "features": qtm_features
    }

    if compact:
        return qtmcompact(qtm_geosjon)

    else: return qtm_geosjon
                          
    

def geojson2qtm(geojson_data, resolution, compact=False):
    """
    Convert GeoJSON data to QTM DGGS format.
    
    Args:
        geojson_data (dict): GeoJSON data as a dictionary
        resolution (int): Resolution level [1..24]
        compact (bool): Whether to use compact mode
        
    Returns:
        dict: GeoJSON data in QTM DGGS format
    """
    if resolution < 1 or resolution > 24:
        raise ValueError("Resolution must be in range [1..24]")
    
    geojson_features = []

    # Process GeoJSON features
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

def is_url(path):
    """Check if the given path is a URL."""
    try:
        result = urlparse(path)
        return all([result.scheme, result.netloc])
    except:
        return False

def read_geojson_file(geojson_path):
    """Read GeoJSON from either a local file or URL."""
    if is_url(geojson_path):
        try:
            response = requests.get(geojson_path)
            response.raise_for_status()
            return json.loads(response.text)
        except requests.RequestException as e:
            print(f"Error: Failed to download GeoJSON from URL {geojson_path}: {str(e)}")
            return None
    else:
        if not os.path.exists(geojson_path):
            print(f"Error: The file {geojson_path} does not exist.")
            return None
        try:
            with open(geojson_path, 'r', encoding='utf-8') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error reading GeoJSON file: {e}")
            return None

def geojson2qtm_cli():
    """Command line interface for converting GeoJSON to QTM DGGS format."""
    parser = argparse.ArgumentParser(description="Convert GeoJSON to QTM DGGS")
    parser.add_argument('-r', '--resolution', type=int, required=True, help="Resolution [1..24]")
    parser.add_argument(
        '-geojson', '--geojson', type=str, required=True, 
        help="GeoJSON file path or URL (Point, Polyline or Polygon)"
    )
    parser.add_argument('-compact', action='store_true', help="Enable Tilecode compact mode")

    args = parser.parse_args()
    
    # Read GeoJSON data from file or URL
    geojson_data = read_geojson_file(args.geojson)
    if geojson_data is None:
        return

    try:
        # Convert the GeoJSON data
        result = geojson2qtm(geojson_data, args.resolution, args.compact)
        
        # Save the result
        geojson_name = os.path.splitext(os.path.basename(args.geojson))[0]
        geojson_path = f"{geojson_name}2qtm_{args.resolution}.geojson"
        if args.compact:        
            geojson_path = f"{geojson_name}2qtm_{args.resolution}_compacted.geojson"
            
        with open(geojson_path, 'w') as f:
            json.dump(result, f)

        print(f"GeoJSON saved as {geojson_path}")
        
    except ValueError as e:
        print(f"Error: {str(e)}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")