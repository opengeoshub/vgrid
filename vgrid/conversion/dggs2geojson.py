from vgrid.utils import s2, olc, geohash, georef, mgrs, mercantile, maidenhead
from vgrid.utils.gars import garsgrid
from vgrid.utils.qtm import constructGeometry, qtm_id_to_facet
import h3

from vgrid.utils.rhealpixdggs.dggs import RHEALPixDGGS
from vgrid.utils.rhealpixdggs.utils import my_round
from vgrid.utils.rhealpixdggs.ellipsoids import WGS84_ELLIPSOID
import platform

if (platform.system() == 'Windows'):   
    from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
    from vgrid.utils.eaggr.eaggr import Eaggr
    from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
    from vgrid.utils.eaggr.enums.model import Model
    from vgrid.generator.isea4tgrid import fix_isea4t_wkt, fix_isea4t_antimeridian_cells


if (platform.system() == 'Linux'):
    from vgrid.utils.dggrid4py import DGGRIDv7, dggs_types
    from vgrid.utils.dggrid4py.dggrid_runner import input_address_types


from vgrid.utils.easedggs.constants import levels_specs
from vgrid.utils.easedggs.dggs.grid_addressing import grid_ids_to_geos

from shapely.wkt import loads
from shapely.geometry import shape, Point,Polygon,mapping

import json, re,os,argparse
from vgrid.generator.h3grid import fix_h3_antimeridian_cells
from vgrid.generator.rhealpixgrid import fix_rhealpix_antimeridian_cells

from vgrid.utils.antimeridian import fix_polygon

from vgrid.generator.settings import graticule_dggs_to_feature, geodesic_dggs_to_feature,isea3h_accuracy_res_dict


from pyproj import Geod
geod = Geod(ellps="WGS84")
E = WGS84_ELLIPSOID

def h32geojson(h3_id):
    cell_boundary = h3.cell_to_boundary(h3_id)   
    h3_features = [] 
    if cell_boundary:
        filtered_boundary = fix_h3_antimeridian_cells(cell_boundary)
        # Reverse lat/lon to lon/lat for GeoJSON compatibility
        reversed_boundary = [(lon, lat) for lat, lon in filtered_boundary]
        cell_polygon = Polygon(reversed_boundary)
        resolution = h3.get_resolution(h3_id) 
        num_edges = 6
        if (h3.is_pentagon(h3_id)):
            num_edges = 5
        h3_feature = geodesic_dggs_to_feature("h3",h3_id,resolution,cell_polygon,num_edges)   
        h3_features.append(h3_feature)

    return {
        "type": "FeatureCollection",
        "features": h3_features
    }
       
    
def h32geojson_cli():
    """
    Command-line interface for h32geojson.
    """
    parser = argparse.ArgumentParser(description="Convert H3 cell ID to GeoJSON")
    parser.add_argument("h3", help="Input H3 cell ID, e.g., h32geojson 8d65b56628e46bf")
    args = parser.parse_args()
    geojson_data = json.dumps(h32geojson(args.h3))
    print(geojson_data)
    

def s22geojson(s2_token):
    # Create an S2 cell from the given cell ID
    s2_features = [] 
    cell_id = s2.CellId.from_token(s2_token)
    cell = s2.Cell(cell_id)
    if cell:
        # Get the vertices of the cell (4 vertices for a rectangular cell)
        vertices = [cell.get_vertex(i) for i in range(4)]
        # Prepare vertices in (longitude, latitude) format for Shapely
        shapely_vertices = []
        for vertex in vertices:
            lat_lng = s2.LatLng.from_point(vertex)  # Convert Point to LatLng
            longitude = lat_lng.lng().degrees  # Access longitude in degrees
            latitude = lat_lng.lat().degrees   # Access latitude in degrees
            shapely_vertices.append((longitude, latitude))

        # Close the polygon by adding the first vertex again
        shapely_vertices.append(shapely_vertices[0])  # Closing the polygon
        # Create a Shapely Polygon
        cell_polygon = fix_polygon(Polygon(shapely_vertices)) # Fix antimeridian
        resolution = cell_id.level()
        num_edges = 4
        s2_feature = geodesic_dggs_to_feature("s2",s2_token,resolution,cell_polygon,num_edges)   
        s2_features.append(s2_feature)

    return {
        "type": "FeatureCollection",
        "features": s2_features
    }
 
       
def s22geojson_cli():
    """
    Command-line interface for s22geojson.
    """
    parser = argparse.ArgumentParser(description="Convert S2 cell token to GeoJSON")
    parser.add_argument("s2", help="Input S2 cell token, e.g., s22geojson 31752f45cc94")
    args = parser.parse_args()
    geojson_data = json.dumps(s22geojson(args.s2))
    print(geojson_data)


def rhealpix_cell_to_polygon(cell):
    vertices = [tuple(my_round(coord, 14) for coord in vertex) for vertex in cell.vertices(plane=False)]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    vertices = fix_rhealpix_antimeridian_cells(vertices)
    return Polygon(vertices)

def rhealpix2geojson(rhealpix_id):
    rhealpix_uids = (rhealpix_id[0],) + tuple(map(int, rhealpix_id[1:]))
    rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3) 
    rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
    rhealpix_features = []
    if rhealpix_cell:
        resolution = rhealpix_cell.resolution        
        cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
        num_edges = 4
        if rhealpix_cell.ellipsoidal_shape() == 'dart':
            num_edges = 3
        rhealpix_feature = geodesic_dggs_to_feature("rhealpix",rhealpix_id,resolution,cell_polygon,num_edges)   
        rhealpix_features.append(rhealpix_feature)

    return {
        "type": "FeatureCollection",
        "features": rhealpix_features
    }
 
def rhealpix2geojson_cli():
    """
    Command-line interface for rhealpix2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Rhealpix cell ID to GeoJSON")
    parser.add_argument("rhealpix", help="Input Rhealpix cell ID, e.g., rhealpix2geojson R31260335553825")
    args = parser.parse_args()
    geojson_data = json.dumps(rhealpix2geojson(args.rhealpix))
    print(geojson_data)
    

def isea4t2geojson(isea4t_id):
    if (platform.system() == 'Windows'): 
        isea4t_dggs = Eaggr(Model.ISEA4T)
        cell_to_shape = isea4t_dggs.convert_dggs_cell_outline_to_shape_string(DggsCell(isea4t_id),ShapeStringFormat.WKT)
        cell_to_shape_fixed = loads(fix_isea4t_wkt(cell_to_shape))
        if isea4t_id.startswith('00') or isea4t_id.startswith('09') or isea4t_id.startswith('14')\
            or isea4t_id.startswith('04') or isea4t_id.startswith('19'):
            cell_to_shape_fixed = fix_isea4t_antimeridian_cells(cell_to_shape_fixed)
        
        isea4t_features = []
        if cell_to_shape_fixed:
            resolution = len(isea4t_id)-2
            num_edges = 3
            cell_polygon = Polygon(list(cell_to_shape_fixed.exterior.coords))
            isea4t_feature = geodesic_dggs_to_feature("isea4t",isea4t_id,resolution,cell_polygon,num_edges)   
            isea4t_features.append(isea4t_feature)

        return {
            "type": "FeatureCollection",
            "features": isea4t_features
        }

def isea4t2geojson_cli():
    """
    Command-line interface for isea4t2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Open-Eaggr ISEA4T cell ID to GeoJSON")
    parser.add_argument("isea4t", help="Input isea4t code, e.g., isea4t2geojson 131023133313201333311333")
    args = parser.parse_args()
    geojson_data = json.dumps(isea4t2geojson(args.isea4t))
    print(geojson_data)


def isea3h_cell_to_polygon(isea3h_id):
    if (platform.system() == 'Windows'):
        isea3h_dggs = Eaggr(Model.ISEA3H)
        cell_to_shape = isea3h_dggs.convert_dggs_cell_outline_to_shape_string(DggsCell(isea3h_id),ShapeStringFormat.WKT)
        if cell_to_shape:
            coordinates_part = cell_to_shape.replace("POLYGON ((", "").replace("))", "")
            coordinates = []
            for coord_pair in coordinates_part.split(","):
                lon, lat = map(float, coord_pair.strip().split())
                coordinates.append([lon, lat])

            # Ensure the polygon is closed (first and last point must be the same)
            if coordinates[0] != coordinates[-1]:
                coordinates.append(coordinates[0])

        cell_polygon = Polygon(coordinates)
        fixed_polygon = fix_polygon(cell_polygon)    
        return fixed_polygon


def isea3h2geojson(isea3h_id):
    if (platform.system() == 'Windows'):
        isea3h_dggs = Eaggr(Model.ISEA3H)
        cell_polygon = isea3h_cell_to_polygon(isea3h_id)
    
        cell_centroid = cell_polygon.centroid
        center_lat =  round(cell_centroid.y, 7)
        center_lon = round(cell_centroid.x, 7)
        
        cell_area = round(abs(geod.geometry_area_perimeter(cell_polygon)[0]),3)
        cell_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])
        
        isea3h2point = isea3h_dggs.convert_dggs_cell_to_point(DggsCell(isea3h_id))      
        accuracy = isea3h2point._accuracy
            
        avg_edge_len = cell_perimeter / 6
        resolution  = isea3h_accuracy_res_dict.get(accuracy)
        
        if (resolution == 0): # icosahedron faces at resolution = 0
            avg_edge_len = cell_perimeter / 3
        
        if accuracy == 0.0:
            if round(avg_edge_len,2) == 0.06:
                resolution = 33
            elif round(avg_edge_len,2) == 0.03:
                resolution = 34
            elif round(avg_edge_len,2) == 0.02:
                resolution = 35
            elif round(avg_edge_len,2) == 0.01:
                resolution = 36
            
            elif round(avg_edge_len,3) == 0.007:
                resolution = 37
            elif round(avg_edge_len,3) == 0.004:
                resolution = 38
            elif round(avg_edge_len,3) == 0.002:
                resolution = 39
            elif round(avg_edge_len,3) <= 0.001:
                resolution = 40
                
        feature = {
            "type": "Feature",
            "geometry": mapping(cell_polygon),
            "properties": {
                    "isea3h": isea3h_id,
                    "resolution": resolution,
                    "center_lat": center_lat,
                    "center_lon": center_lon,
                    "avg_edge_len": round(avg_edge_len,3),
                    "cell_area": cell_area
                    }
        }

        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        return  feature_collection


def isea3h2geojson_cli():
    """
    Command-line interface for isea3h2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert ISEA3H cell ID to GeoJSON")
    parser.add_argument("isea3h", help="Input ISEA3H cell ID, e.g., isea3h2geojson 1327916769,-55086")
    args = parser.parse_args()
    geojson_data = json.dumps(isea3h2geojson(args.isea3h))
    print(geojson_data)


def dggrid2geojson(dggrid_id,dggs_type,resolution):
    if (platform.system() == 'Linux'):
        dggrid_instance = DGGRIDv7(executable='/usr/local/bin/dggrid', working_dir='.', capture_logs=False, silent=True, tmp_geo_out_legacy=False, debug=False)
        dggrid_cell =  dggrid_instance.grid_cell_polygons_from_cellids([dggrid_id],dggs_type,resolution,split_dateline=True)    
      
        gdf = dggrid_cell.set_geometry("geometry")  # Ensure the geometry column is set
        # Check and set CRS to EPSG:4326 if needed
        if gdf.crs is None:
            gdf.set_crs(epsg=4326, inplace=True)
        elif not gdf.crs.equals("EPSG:4326"):
            gdf = gdf.to_crs(epsg=4326)
        
        feature_collection = gdf.to_json()
        return feature_collection


def dggrid2geojson_cli():
    """
    Command-line interface for dggrid2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert DGGRID code to GeoJSON. \
                                     Usage: dggrid2geojson <SEQNUM> <dggs_type> <res>. \
                                     Ex: dggrid2geojson 783229476878 ISEA7H 13")
    parser.add_argument("dggrid", help="Input DGGRID code in SEQNUM format")
    parser.add_argument("type", choices=dggs_types, help="Select a DGGS type from the available options.")
    parser.add_argument("res", type=int, help="resolution")
    # parser.add_argument("address", choices=input_address_types, help="Address type")

    args = parser.parse_args()
    geojson_data = dggrid2geojson(args.dggrid,args.type, args.res)
    print(geojson_data)


def ease2geojson(ease_id):
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

    ease_features = []
    if cell_polygon:
        resolution = level
        num_edges = 4
        ease_feature = geodesic_dggs_to_feature("ease",ease_id,resolution,cell_polygon,num_edges)   
        ease_features.append(ease_feature)

    return {
        "type": "FeatureCollection",
        "features": ease_features
    }


def ease2geojson_cli():
    """
    Command-line interface for ease2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert EASE-DGGS code to GeoJSON")
    parser.add_argument("ease", help="Input ASE-DGGS code, e.g., ease2geojson L4.165767.02.02.20.71")
    args = parser.parse_args()
    geojson_data = json.dumps(ease2geojson(args.ease))
    print(geojson_data)

def qtm2geojson(qtm_id):
    facet = qtm_id_to_facet(qtm_id)
    cell_polygon = constructGeometry(facet)    
    resolution = len(qtm_id)
    num_edges = 3
    qtm_features = []
    qtm_feature = geodesic_dggs_to_feature("qtm",qtm_id,resolution,cell_polygon,num_edges)   
    qtm_features.append(qtm_feature)

    return {
        "type": "FeatureCollection",
        "features": qtm_features
    }

def qtm2geojson_cli():
    """
    Command-line interface for qtm2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert QTM cell ID to GeoJSON")
    parser.add_argument("qtm", help="Input QTM cell ID, e.g., qtm2geojson 42012323")
    args = parser.parse_args()
    geojson_data = json.dumps(qtm2geojson(args.qtm))
    print(geojson_data)

    
def olc2geojson(olc_id):
    # Decode the Open Location Code into a CodeArea object
    coord = olc.decode(olc_id)    
    if coord:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
        max_lat, max_lon = coord.latitudeHi, coord.longitudeHi        
        resolution = coord.codeLength 

        # Define the polygon based on the bounding box
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        olc_features = []
        olc_feature = graticule_dggs_to_feature("olc",olc_id,resolution,cell_polygon)   
        olc_features.append(olc_feature)

    return {
        "type": "FeatureCollection",
        "features": olc_features
    }

def olc2geojson_cli():
    """
    Command-line interface for olc2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert OLC/ Google Plus Codes to GeoJSON")
    parser.add_argument("olc", help="Input OLC, e.g., olc2geojson 7P28QPG4+4P7")
    args = parser.parse_args()
    geojson_data = json.dumps(olc2geojson(args.olc))
    print(geojson_data)


def geohash2geojson(geohash_id):
    # Decode the Open Location Code into a CodeArea object
    bbox =  geohash.bbox(geohash_id)
    geohash_features = []
    if bbox:
        min_lat, min_lon = bbox['s'], bbox['w']  # Southwest corner
        max_lat, max_lon = bbox['n'], bbox['e']  # Northeast corner        
        resolution =  len(geohash_id)

        # Define the polygon based on the bounding box
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        geohash_feature = graticule_dggs_to_feature("geohash",geohash_id,resolution,cell_polygon)   
        geohash_features.append(geohash_feature)

    return {
        "type": "FeatureCollection",
        "features": geohash_features
    }
       
    
def geohash2geojson_cli():
    """
    Command-line interface for geohash2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Geohash cell ID to GeoJSON")
    parser.add_argument("geohash", help="Input Geohash cell ID, e.g., geohash2geojson w3gvk1td8")
    args = parser.parse_args()
    geojson_data = json.dumps(geohash2geojson(args.geohash))
    print(geojson_data)


def mgrs2geojson_old(mgrs_id,lat=None,lon=None):
    origin_lat, origin_lon, min_lat, min_lon, max_lat, max_lon,resolution = mgrs.mgrscell(mgrs_id)
    if origin_lat:
        # Define the polygon based on the bounding box
        origin_lat = round(origin_lat,7)
        origin_lon = round(origin_lon,7)
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
                      
        feature = {
            "type": "Feature",
            "geometry": mapping(cell_polygon),
            "properties": {
                "mgrs": mgrs_id,
                "resolution": resolution,
                "origin_lat": origin_lat,
                "origin_lon": origin_lon,
                "avg_edge_len" : avg_edge_len,
                "cell_area": cell_area
                }
            }
        
        if lat is not None and lon is not None:
            # Load the GZD JSON file (treated as GeoJSON format) from the same folder
            gzd_json_path = os.path.join(os.path.dirname(__file__), 'gzd.geojson')
            with open(gzd_json_path) as f:
                gzd_json = json.load(f)

            # Convert GZD GeoJSON features to Shapely polygons
            gzd_polygons = [
                {"geometry": shape(feature["geometry"]), "properties": feature["properties"]}
                for feature in gzd_json["features"]
            ]

            # Perform the intersection with the MGRS polygon
            intersection_features = []
            for gzd_polygon_data in gzd_polygons:
                gzd_polygon = gzd_polygon_data["geometry"]
                if cell_polygon.intersects(gzd_polygon):
                    # Find the intersection polygon
                    intersection_polygon = cell_polygon.intersection(gzd_polygon)
                    # intersection_min_lon, intersection_min_lat, intersection_max_lon, intersection_max_lat = intersection_polygon.bounds  # Bounds of the polygon
                    interection_area = round(abs(geod.geometry_area_perimeter(intersection_polygon)[0]),2)  # Area in square meters     
                    intersection_perimeter = abs(geod.geometry_area_perimeter(intersection_polygon)[1])       
                    avg_edge_len = round(intersection_perimeter / 4,2)
                    # intersection_width = round(geod.line_length([intersection_min_lon, intersection_max_lon], [intersection_min_lat, intersection_min_lat]),3)
                    # intersection_height = round(geod.line_length([intersection_min_lon, intersection_min_lon], [intersection_min_lat, intersection_max_lat]),3)

                    # Convert lat/lon to a Shapely point
                    point = Point(lon, lat)

                    # Check if the point is inside the intersection polygon
                    if intersection_polygon.contains(point):
                        # Manually construct the intersection as a JSON-like structure
                        intersection_feature = {
                            "type": "Feature",
                            "geometry": {
                                "type": "Polygon",
                                "coordinates": [list(intersection_polygon.exterior.coords)]
                            },
                            "properties": {
                                "mgrs": mgrs_id,
                                "resolution": resolution,
                                "origin_lat": origin_lat,
                                "origin_lon": origin_lon,
                                "avg_edge_len": avg_edge_len,
                                "cell_area": interection_area
                                # **gzd_polygon_data["properties"],  # Include properties from GZD
                            }
                        }
                        intersection_features.append(intersection_feature)

            # If intersections are found, wrap them in a FeatureCollection
            if intersection_features:
                intersection_feature_collection = {
                    "type": "FeatureCollection",
                    "features": intersection_features
                }
                return intersection_feature_collection


        # If no intersection or point not contained, return the original MGRS GeoJSON
        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        
        return feature_collection

def mgrs2geojson(mgrs_id,lat=None,lon=None):
    origin_lat, origin_lon, min_lat, min_lon, max_lat, max_lon,resolution = mgrs.mgrscell(mgrs_id)
    if min_lat:
        # Define the polygon based on the bounding box
        # origin_lat = round(origin_lat,7)
        # origin_lon = round(origin_lon,7)
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        
        mgrs_feature = graticule_dggs_to_feature("mgrs",mgrs_id,resolution,cell_polygon)                      
      
        if lat is not None and lon is not None:
            # Load the GZD JSON file (treated as GeoJSON format) from the same folder
            gzd_json_path = os.path.join(os.path.dirname(__file__), 'gzd.geojson')
            with open(gzd_json_path) as f:
                gzd_json = json.load(f)

            # Convert GZD GeoJSON features to Shapely polygons
            gzd_polygons = [
                {"geometry": shape(gzd_feature["geometry"]), "properties": gzd_feature["properties"]}
                for gzd_feature in gzd_json["features"]
            ]

            # Perform the intersection with the MGRS polygon
            intersection_features = []
            for gzd_polygon_data in gzd_polygons:
                gzd_polygon = gzd_polygon_data["geometry"]
                if cell_polygon.intersects(gzd_polygon):
                    # Find the intersection polygon
                    intersection_polygon = cell_polygon.intersection(gzd_polygon)
                    min_lon, min_lat, max_lon, max_lat = intersection_polygon.bounds
                    interection_center_lat = round((min_lat + max_lat) / 2,7)
                    interection_center_lon = round((min_lon + max_lon) / 2,7)
                    interection_cell_width = round(geod.line_length([min_lon, max_lon], [min_lat, min_lat]),3)
                    interection_cell_height = round(geod.line_length([min_lon, min_lon], [min_lat, max_lat]),3)
                    interection_area = round(abs(geod.geometry_area_perimeter(intersection_polygon)[0]),3)  # Area in square meters     
                    point = Point(lon, lat)
                    # Check if the point is inside the intersection polygon
                    if intersection_polygon.contains(point):                        
                        # Manually construct the intersection as a JSON-like structure
                        intersection_feature = {
                            "type": "Feature",
                            "geometry": {
                                "type": "Polygon",
                                "coordinates": [list(intersection_polygon.exterior.coords)]
                            },
                            "properties": {
                                "mgrs": mgrs_id,
                                "resolution": resolution,
                                "center_lat": interection_center_lat,
                                "center_lon": interection_center_lon,
                                "cell_width": interection_cell_width,
                                "cell_height": interection_cell_height,
                                "cell_area": interection_area
                                # **gzd_polygon_data["properties"],  # Include properties from GZD
                            }
                        }
                        intersection_features.append(intersection_feature)

            # If intersections are found, wrap them in a FeatureCollection
            if intersection_features:
                intersection_feature_collection = {
                    "type": "FeatureCollection",
                    "features": intersection_features
                }
                return intersection_feature_collection


        # If no intersection or point not contained, return the original MGRS GeoJSON
        feature_collection = {
            "type": "FeatureCollection",
            "features": [mgrs_feature]
        }
        
        return feature_collection
    
def mgrs2geojson_cli():
    """
    Command-line interface for mgrs2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert MGRS cell ID to GeoJSON")
    parser.add_argument("mgrs", help="Input MGRS cell ID, e.g., mgrs2geojson 48PXS866916")
    args = parser.parse_args()
    geojson_data = json.dumps(mgrs2geojson(args.mgrs))
    print(geojson_data)


def georef2geojson(georef_id):
    center_lat, center_lon, min_lat, min_lon, max_lat, max_lon,resolution = georef.georefcell(georef_id)
    georef_features = []
    if center_lat:
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])

        georef_feature = graticule_dggs_to_feature("georef",georef_id,resolution,cell_polygon)   
        georef_features.append(georef_feature)

    return {
        "type": "FeatureCollection",
        "features": georef_features
    }
    
def georef2geojson_cli():
    """
    Command-line interface for georef2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert GEOREF code to GeoJSON")
    parser.add_argument("georef", help="Input GEOREF code, e.g., georef2geojson VGBL42404651")
    args = parser.parse_args()
    geojson_data = json.dumps(georef2geojson(args.georef))
    print(geojson_data)


def tilecode2geojson(tilecode_id):  
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode_id)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Get the bounds of the tile in (west, south, east, north)
    bounds = mercantile.bounds(x, y, z)    
    tilecode_features = []
    if bounds:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        
        resolution = z
        tilecode_feature = graticule_dggs_to_feature("tilecode",tilecode_id,resolution,cell_polygon)   
        tilecode_features.append(tilecode_feature)

    return {
        "type": "FeatureCollection",
        "features": tilecode_features
    }
       
def tilecode2geojson_cli():
    """
    Command-line interface for tilecode2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Tilecode to GeoJSON")
    parser.add_argument("tilecode", help="Input Tilecode, e.g. z0x0y0")
    args = parser.parse_args()

    # Generate the GeoJSON feature
    geojson_data = json.dumps(tilecode2geojson(args.tilecode))
    print(geojson_data)


def quadkey2geojson(quadkey_id):
    tile = mercantile.quadkey_to_tile(quadkey_id)    
    # Format as tilecode
    z = tile.z
    x = tile.x
    y = tile.y
    # Get the bounds of the tile in (west, south, east, north)
    bounds = mercantile.bounds(x, y, z)    
    quadkey_features = []
    
    if bounds:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east

        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        
        resolution = z
        quadkey_feature = graticule_dggs_to_feature("quadkey",quadkey_id,resolution,cell_polygon)   
        quadkey_features.append(quadkey_feature)
    
    return {
        "type": "FeatureCollection",
        "features": quadkey_features
    }
        
def quadkey2geojson_cli():
    """
    Command-line interface for quadkey2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Quadkey to GeoJSON")
    parser.add_argument("quadkey", help="Input Quadkey, e.g. 13223011131020220011133")
    args = parser.parse_args()

    # Generate the GeoJSON feature
    geojson_data = json.dumps(quadkey2geojson(args.quadkey))
    print(geojson_data)


def maidenhead2geojson(maidenhead_id):
    # Decode the Open Location Code into a CodeArea object
    _, _, min_lat, min_lon, max_lat, max_lon, _ = maidenhead.maidenGrid(maidenhead_id)
    maidenhead_features = []
    if min_lat:
        resolution = int(len(maidenhead_id)/2)   
        # Define the polygon based on the bounding box
        cell_polygon = Polygon([
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ])
        maidenhead_feature = graticule_dggs_to_feature("gars",maidenhead_id,resolution,cell_polygon)   
        maidenhead_features.append(maidenhead_feature)

    return {
        "type": "FeatureCollection",
        "features": maidenhead_features
    }

def maidenhead2geojson_cli():
    """
    Command-line interface for maidenhead2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Maidenhead cell ID to GeoJSON")
    parser.add_argument("maidenhead", help="Input Maidenhead cell ID, e.g., maidenhead2geojson OK30is46")
    args = parser.parse_args()
    geojson_data = json.dumps(maidenhead2geojson(args.maidenhead))
    print(geojson_data)


def gars2geojson(gars_id):
    gars_grid = garsgrid.GARSGrid(gars_id)
    wkt_polygon = gars_grid.polygon
    gars_features = []
    
    if wkt_polygon:
        # # Create the bounding box coordinates for the polygon
        resolution_minute = gars_grid.resolution        
        cell_polygon = Polygon(list(wkt_polygon.exterior.coords))
        gars_feature = graticule_dggs_to_feature("gars",gars_id,resolution_minute,cell_polygon)   
        gars_features.append(gars_feature)

    return {
        "type": "FeatureCollection",
        "features": gars_features
    }
    

def gars2geojson_cli():
    """
    Command-line interface for gars2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert GARS cell ID to GeoJSON")
    parser.add_argument("gars", help="Input GARS cell ID, e.g., gars2geojson 574JK1918")
    args = parser.parse_args()
    geojson_data = json.dumps(gars2geojson(args.gars))
    print(geojson_data)