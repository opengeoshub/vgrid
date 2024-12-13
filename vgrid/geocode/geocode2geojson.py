from vgrid.utils import geohash
import h3
from vgrid.utils import georef
from vgrid.utils import maidenhead
from vgrid.utils import mgrs
from vgrid.utils import olc
from vgrid.utils.s2 import LatLng, CellId

from vgrid.utils.gars.garsgrid import GARSGrid
from vgrid.utils import mercantile

from rhealpixdggs.dggs import RHEALPixDGGS
from rhealpixdggs.utils import my_round
from rhealpixdggs.ellipsoids import WGS84_ELLIPSOID

from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.dggs_shape import DggsShape
from vgrid.utils.eaggr.shapes.dggs_polygon import DggsPolygon
from vgrid.utils.eaggr.enums.dggs_shape_location import DggsShapeLocation
from vgrid.utils.eaggr.enums.model import Model
from shapely.wkt import loads
from shapely.geometry import Polygon,mapping
from pyproj import Geod
import math
import json, re
import argparse

from vgrid.utils import s2
from vgrid.utils import tilecode

def haversine(lat1, lon1, lat2, lon2):
    # Radius of the Earth in meters
    R = 6371000  

    # Convert latitude and longitude from degrees to radians
    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lon2 - lon1)

    # Haversine formula
    a = math.sin(delta_phi / 2) ** 2 + math.cos(phi1) * math.cos(phi2) * math.sin(delta_lambda / 2) ** 2
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))

    return R * c  # Distance in meters

def olc2geojson(olc_code):
    # Decode the Open Location Code into a CodeArea object
    coord = olc.decode(olc_code)
    
    if coord:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
        max_lat, max_lon = coord.latitudeHi, coord.longitudeHi

        center_lat, center_lon = coord.latitudeCenter, coord.longitudeCenter
        resolution = coord.codeLength 
       
        lat_len = haversine(min_lat, min_lon, max_lat, min_lon)
        lon_len = haversine(min_lat, min_lon, min_lat, max_lon)

        bbox_width =  f'{round(lon_len,1)} m'
        bbox_height =  f'{round(lat_len,1)} m'
        if lon_len >= 10000:
            bbox_width = f'{round(lon_len/1000,1)} km'
            bbox_height = f'{round(lat_len/1000,1)} km'

        # Define the polygon based on the bounding box
        polygon_coords = [
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ]
        
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [polygon_coords]
            },
            "properties": {
                "olc": olc_code,  # Include the OLC as a property
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "resolution": resolution
            }
        }

        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        
        return feature_collection

def olc2geojson_cli():
    """
    Command-line interface for olc2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert OLC/ Google Plus Codes to GeoJSON")
    parser.add_argument("olc", help="Input OLC, e.g., olc2geojson 7P28QPG4+4P7")
    args = parser.parse_args()
    geojson_data = json.dumps(olc2geojson(args.olc))
    print(geojson_data)


def geohash2geojson(geohash_code):
    # Decode the Open Location Code into a CodeArea object
    bbox =  geohash.bbox(geohash_code)
    if bbox:
        min_lat, min_lon = bbox['s'], bbox['w']  # Southwest corner
        max_lat, max_lon = bbox['n'], bbox['e']  # Northeast corner
        
        center_lat = (min_lat + max_lat) / 2
        center_lon = (min_lon + max_lon) / 2
        
        resolution =  len(geohash_code)

        lat_len = haversine(min_lat, min_lon, max_lat, min_lon)
        lon_len = haversine(min_lat, min_lon, min_lat, max_lon)
 
        bbox_width =  f'{round(lon_len,1)} m'
        bbox_height =  f'{round(lat_len,1)} m'
        if lon_len >= 10000:
            bbox_width = f'{round(lon_len/1000,1)} km'
            bbox_height = f'{round(lat_len/1000,1)} km'
            
        # Define the polygon based on the bounding box
        polygon_coords = [
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ]

        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [polygon_coords]  # Directly use the coordinates list
            },
            "properties": {
                "geohash": geohash_code,
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "resolution": resolution
                }
            }
        
        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        
        return feature_collection
    
def geohash2geojson_cli():
    """
    Command-line interface for geohash2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Geohash code to GeoJSON")
    parser.add_argument("geohash", help="Input Geohash code, e.g., geohash2geojson w3gvk1td8")
    args = parser.parse_args()
    geojson_data = json.dumps(geohash2geojson(args.geohash))
    print(geojson_data)


def mgrs2geojson(mgrs_code,lat=None,lon=None):
    origin_lat, origin_lon, min_lat, min_lon, max_lat, max_lon,resolution = mgrs.mgrscell(mgrs_code)

    lat_len = haversine(min_lat, min_lon, max_lat, min_lon)
    lon_len = haversine(min_lat, min_lon, min_lat, max_lon)
  
    bbox_width =  f'{round(lon_len,1)} m'
    bbox_height =  f'{round(lat_len,1)} m'
    
    if lon_len >= 10000:
        bbox_width = f'{round(lon_len/1000,1)} km'
        bbox_height = f'{round(lat_len/1000,1)} km'
        
    if origin_lat:
        # Define the polygon based on the bounding box
        polygon_coords = [
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ]

        mgrs_polygon = Polygon(polygon_coords)
        
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [polygon_coords]  # Directly use the coordinates list
            },
            "properties": {
                "mgrs": mgrs_code,
                "origin_lat": origin_lat,
                "origin_lon": origin_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "resolution": resolution
                }
            }
        
        # if lat is not None and lon is not None:
        #     # Load the GZD JSON file (treated as GeoJSON format) from the same folder
        #     gzd_json_path = os.path.join(os.path.dirname(__file__), 'gzd.geojson')
        #     with open(gzd_json_path) as f:
        #         gzd_json = json.load(f)

        #     # Convert the GZD JSON to a GeoDataFrame
        #     gzd_gdf = gpd.GeoDataFrame.from_features(gzd_json['features'], crs="EPSG:4326")

        #     # Convert the MGRS polygon to a GeoDataFrame for intersection
        #     mgrs_gdf = gpd.GeoDataFrame(geometry=[mgrs_polygon], crs="EPSG:4326")

        #     # Perform the intersection
        #     intersection_gdf = gpd.overlay(mgrs_gdf, gzd_gdf, how='intersection')

        #     # Check if the intersection result is empty
        #     if not intersection_gdf.empty:
        #         # Convert lat/lon to a Shapely point
        #         point = Point(lon, lat)

        #         # Check if the point is inside any of the intersection polygons
        #         for intersection_polygon in intersection_gdf.geometry:
        #             if intersection_polygon.contains(point):
        #                 # Manually construct the intersection as a JSON-like structure
        #                 intersection_feature = {
        #                     "type": "Feature",
        #                     "geometry": {
        #                         "type": "Polygon",
        #                         "coordinates": [list(intersection_polygon.exterior.coords)]
        #                     },
        #                     "properties": {
        #                         "mgrs": mgrs_code,
        #                         "origin_lat": origin_lat,
        #                         "origin_lon": origin_lon,
        #                         "bbox_height": bbox_height,
        #                         "bbox_width": bbox_width,
        #                         "resolution": resolution
        #                     }
        #                 }

        #                 # Wrap the feature in a FeatureCollection
        #                 intersection_feature_collection = {
        #                     "type": "FeatureCollection",
        #                     "features": [intersection_feature]
        #                 }

        #                 return intersection_feature_collection

        # If no intersection or point not contained, return the original MGRS GeoJSON
        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        
        return feature_collection
    
def mgrs2geojson_cli():
    """
    Command-line interface for mgrs2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert MGRS code to GeoJSON")
    parser.add_argument("mgrs", help="Input MGRS code, e.g., mgrs2geojson 34TGK56063228")
    args = parser.parse_args()
    geojson_data = json.dumps(mgrs2geojson(args.mgrs))
    print(geojson_data)


def georef2geojson(georef_code):
    center_lat, center_lon, min_lat, min_lon, max_lat, max_lon,resolution = georef.georefcell(georef_code)

    lat_len = haversine(min_lat, min_lon, max_lat, min_lon)
    lon_len = haversine(min_lat, min_lon, min_lat, max_lon)
  
    bbox_width =  f'{round(lon_len,1)} m'
    bbox_height =  f'{round(lat_len,1)} m'
    
    if lon_len >= 10000:
        bbox_width = f'{round(lon_len/1000,1)} km'
        bbox_height = f'{round(lat_len/1000,1)} km'
        
    if center_lat:
        # Define the polygon based on the bounding box
        polygon_coords = [
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ]

        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [polygon_coords]  # Directly use the coordinates list
            },
            "properties": {
                "georef": georef_code,
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "resolution": resolution
                }
            }
        
        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        
        return feature_collection

def georef2geojson_cli():
    """
    Command-line interface for georef2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert GEOREF code to GeoJSON")
    parser.add_argument("georef", help="Input GEOREF code, e.g., georef2geojson VGBL42404651")
    args = parser.parse_args()
    geojson_data = json.dumps(georef2geojson(args.georef))
    print(geojson_data)


def h32geojson(h3_code):
    # Get the boundary coordinates of the H3 cell
    boundary = h3.cell_to_boundary(h3_code)
    
    if boundary:
        # Get the center coordinates of the H3 cell
        center_lat, center_lon = h3.cell_to_latlng(h3_code)
        resolution = h3.get_resolution(h3_code)
        avg_edge_len = h3.average_hexagon_edge_length(resolution,unit='m')
        
        boundary = list(boundary)
        # Ensure the polygon boundary is closed
        if boundary[0] != boundary[-1]:
            boundary.append(boundary[0])
        
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
               "coordinates": [[
                [lon, lat] for lat, lon in boundary  # Convert boundary to the correct coordinate order
                ]]
            },
            "properties": {
                "h3": h3_code,
                "center_lat": center_lat,
                "center_lon": center_lon,
                "avg_edge_len": avg_edge_len,
                "resolution": resolution
            }
        }
        # Wrap the feature in a FeatureCollection
        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }

        # Convert the feature collection to JSON formatted string
        return feature_collection

def h32geojson_cli():
    """
    Command-line interface for h32geojson.
    """
    parser = argparse.ArgumentParser(description="Convert H3 code to GeoJSON")
    parser.add_argument("h3", help="Input H3 code, e.g., h32geojson 8d65b56628e46bf")
    args = parser.parse_args()
    geojson_data = json.dumps(h32geojson(args.h3))
    print(geojson_data)

def s22geojson(cell_id_token):
    # Create an S2 cell from the given cell ID
    cell_id = CellId.from_token(cell_id_token)
    cell = s2.Cell(cell_id)
    
    if cell:
        # Get the vertices of the cell (4 vertices for a rectangular cell)
        vertices = [cell.get_vertex(i) for i in range(4)]
        
        # Prepare vertices in [longitude, latitude] format
        json_vertices = []
        for vertex in vertices:
            lat_lng = LatLng.from_point(vertex)  # Convert Point to LatLng
            longitude = lat_lng.lng().degrees  # Access longitude in degrees
            latitude = lat_lng.lat().degrees    # Access latitude in degrees
            json_vertices.append([longitude, latitude])

        # Close the polygon by adding the first vertex again
        json_vertices.append(json_vertices[0])  # Closing the polygon

        # Create a JSON-compatible polygon structure
        json_polygon = {
            "type": "Polygon",
            "coordinates": [json_vertices]
        }

        # Get the center of the cell
        center = cell.get_center()
        center_lat_lng = LatLng.from_point(center)
        center_lat = center_lat_lng.lat().degrees
        center_lon = center_lat_lng.lng().degrees

        # Get rectangular bounds of the cell
        rect_bound = cell.get_rect_bound()
        min_lat = rect_bound.lat_lo().degrees
        max_lat = rect_bound.lat_hi().degrees
        min_lon = rect_bound.lng_lo().degrees
        max_lon = rect_bound.lng_hi().degrees
        
        # Calculate width and height of the bounding box
        lat_len = haversine(min_lat, min_lon, max_lat, min_lon)
        lon_len = haversine(min_lat, min_lon, min_lat, max_lon)
  
        bbox_width = f'{round(lon_len, 1)} m'
        bbox_height = f'{round(lat_len, 1)} m'
        cell_size= cell_id.get_size_ij()

        if lon_len >= 10000:
            bbox_width = f'{round(lon_len / 1000, 1)} km'
            bbox_height = f'{round(lat_len / 1000, 1)} km'

        # Create properties for the Feature
        properties = {
            # "s2": cell_id.id(),
            "s2": cell_id_token,
            "center_lat": center_lat,
            "center_lon": center_lon,
            "bbox_width": bbox_width,
            "bbox_height": bbox_height,
            "cell_size": cell_size,
            "level": cell_id.level()
        }

        # Manually create the Feature and FeatureCollection
        feature = {
            "type": "Feature",
            "geometry": json_polygon,
            "properties": properties
        }
        
        # Create the FeatureCollection
        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }

        # Convert to JSON format
        return feature_collection

def s22geojson_cli():
    """
    Command-line interface for s22geojson.
    """
    parser = argparse.ArgumentParser(description="Convert S2 cell token to GeoJSON")
    parser.add_argument("s2", help="Input S2 cell token, e.g., s22geojson 31752f45cc94")
    args = parser.parse_args()
    geojson_data = json.dumps(s22geojson(args.s2))
    print(geojson_data)

def tilecode2geojson(tilecode):
    """
    Converts a tilecode (e.g., 'z8x11y14') to a GeoJSON Feature with a Polygon geometry
    representing the tile's bounds and includes the original tilecode as a property.

    Args:
        tilecode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        dict: A GeoJSON Feature with a Polygon geometry and tilecode as a property.
    """
    # Extract z, x, y from the tilecode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', tilecode)
    if not match:
        raise ValueError("Invalid tilecode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Get the bounds of the tile in (west, south, east, north)
    bounds = mercantile.bounds(x, y, z)    

    if bounds:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = bounds.south, bounds.west
        max_lat, max_lon = bounds.north, bounds.east

        # tile = mercantile.Tile(x, y, z)
        # quadkey = mercantile.quadkey(tile)

        center_lat = (min_lat + max_lat) / 2
        center_lon = (min_lon + max_lon) / 2
               
        lat_len = haversine(min_lat, min_lon, max_lat, min_lon)
        lon_len = haversine(min_lat, min_lon, min_lat, max_lon)

        bbox_width =  f'{round(lon_len,1)} m'
        bbox_height =  f'{round(lat_len,1)} m'
        if lon_len >= 10000:
            bbox_width = f'{round(lon_len/1000,1)} km'
            bbox_height = f'{round(lat_len/1000,1)} km'

        # Define the polygon based on the bounding box
        polygon_coords = [
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ]
        
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [polygon_coords]
            },
            "properties": {
                "tilecode": tilecode,  # Include the OLC as a property
                # "quadkey": quadkey,
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "resolution": z  # Using the code length as precision
            }
        }

        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        
        return feature_collection

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


def maidenhead2geojson(maidenhead_code):
    # Decode the Open Location Code into a CodeArea object
    center_lat, center_lon, min_lat, min_lon, max_lat, max_lon, _ = maidenhead.maidenGrid(maidenhead_code)
    resolution = int(len(maidenhead_code)/2)
    
    lat_len = haversine(min_lat, min_lon, max_lat, min_lon)
    lon_len = haversine(min_lat, min_lon, min_lat, max_lon)

    bbox_width =  f'{round(lon_len,1)} m'
    bbox_height =  f'{round(lat_len,1)} m'
    
    if lon_len >= 10000:
        bbox_width = f'{round(lon_len/1000,1)} km'
        bbox_height = f'{round(lat_len/1000,1)} km'
        
    if center_lat:
        # Define the polygon based on the bounding box
        polygon_coords = [
            [min_lon, min_lat],  # Bottom-left corner
            [max_lon, min_lat],  # Bottom-right corner
            [max_lon, max_lat],  # Top-right corner
            [min_lon, max_lat],  # Top-left corner
            [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ]

        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [polygon_coords]
            },
            "properties": {
                "maidenhead": maidenhead_code,  # Include the OLC as a property
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "resolution": resolution
            }
        }

        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        
        return feature_collection

def maidenhead2geojson_cli():
    """
    Command-line interface for maidenhead2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Maidenhead code to GeoJSON")
    parser.add_argument("maidenhead", help="Input Maidenhead code, e.g., maidenhead2geojson OK30is46")
    args = parser.parse_args()
    geojson_data = json.dumps(maidenhead2geojson(args.maidenhead))
    print(geojson_data)

# SOS: Convert gars_code object to str first
def gars2geojson(gars_code):
    gars_grid = GARSGrid(gars_code)
    wkt_polygon = gars_grid.polygon
    if wkt_polygon:
        # # Create the bounding box coordinates for the polygon
        x, y = wkt_polygon.exterior.xy
        resolution_minute = gars_grid.resolution
        
        min_lon = min(x)
        max_lon = max(x)
        min_lat = min(y)
        max_lat = max(y)

        # Calculate center latitude and longitude
        center_lon = (min_lon + max_lon) / 2
        center_lat = (min_lat + max_lat) / 2

        # Calculate bounding box width and height
        lat_len = haversine(min_lat, min_lon, max_lat, min_lon)
        lon_len = haversine(min_lat, min_lon, min_lat, max_lon)
 
        bbox_width =  f'{round(lon_len,1)} m'
        bbox_height =  f'{round(lat_len,1)} m'

        if lon_len >= 10000:
            bbox_width = f'{round(lon_len/1000,1)} km'
            bbox_height = f'{round(lat_len/1000,1)} km'

        polygon_coords = list(wkt_polygon.exterior.coords)

        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [polygon_coords]  # Directly use the coordinates list
            },
            "properties": {
                "gars": gars_code,
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "resolution_minute": resolution_minute
                }
            }
        
        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        
        return feature_collection

def gars2geojson_cli():
    """
    Command-line interface for gars2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert GARS code to GeoJSON")
    parser.add_argument("gars", help="Input GARS code, e.g., gars2geojson 574JK1918")
    args = parser.parse_args()
    geojson_data = json.dumps(gars2geojson(args.gars))
    print(geojson_data)

def filter_antimeridian_cells(boundary, threshold=-128):
    if any(lon < threshold for lon, _ in boundary):
        return [(lon - 360 if lon > 0 else lon, lat) for lon, lat in boundary]
    return boundary

def rhealpix_cell_to_polygon(cell):
    vertices = [tuple(my_round(coord, 14) for coord in vertex) for vertex in cell.vertices(plane=False)]
    if vertices[0] != vertices[-1]:
        vertices.append(vertices[0])
    vertices = filter_antimeridian_cells(vertices)
    return Polygon(vertices)


def rhealpix2geojson(rhealpix_code):
    rhealpix_code = str(rhealpix_code)
    rhealpix_uids = (rhealpix_code[0],) + tuple(map(int, rhealpix_code[1:]))
    
    E = WGS84_ELLIPSOID
    rhealpix_dggs = RHEALPixDGGS(ellipsoid=E, north_square=1, south_square=3, N_side=3) 
    rhealpix_cell = rhealpix_dggs.cell(rhealpix_uids)
    resolution = rhealpix_cell.resolution
       
    cell_polygon = rhealpix_cell_to_polygon(rhealpix_cell)
    features = []
    
    if cell_polygon:
        geod = Geod(ellps="WGS84")
        center_lon = cell_polygon.centroid.x
        center_lat = cell_polygon.centroid.y
        
        min_x, min_y, max_x, max_y = cell_polygon.bounds
        cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters                
        _, _, cell_width = geod.inv(min_x, min_y, max_x, min_y)
        _, _, cell_height = geod.inv(min_x, min_y, min_x, max_y)
        
        cell_width_str =  f'{round(cell_width,2)} m'
        cell_height_str =  f'{round(cell_height,2)} m'
        cell_area_str=  f'{round(cell_area,2)} m2'

        if cell_area >= 1000_000:
            cell_width_str = f'{round(cell_width/1000,2)} km'
            cell_height_str = f'{round(cell_height/1000,2)} km'
            cell_area_str = f'{round(cell_area/(10**6),2)} km2'
        
        features. append({
                "type": "Feature",
                "geometry": mapping(cell_polygon),
                "properties": {
                        "rhealpix": str(rhealpix_cell),
                        "center_lat": center_lat,
                        "center_lon": center_lon,
                        "cell_width": cell_width_str,
                        "cell_height": cell_height_str,
                        "cell_area": cell_area_str,
                        "resolution": resolution
                    }
                })
        
        feature_collection =  {
            "type": "FeatureCollection",
            "features": features,
        }                
        return feature_collection

def rhealpix2geojson_cli():
    """
    Command-line interface for rhealpix2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert Rhealpix code to GeoJSON")
    parser.add_argument("rhealpix", help="Input Rhealpix code, e.g., rhealpix2geojson R31260335553825")
    args = parser.parse_args()
    geojson_data = json.dumps(rhealpix2geojson(args.rhealpix))
    print(geojson_data)


def fix_eaggr_wkt(eaggr_wkt):
        # Extract the coordinate section
        coords_section = eaggr_wkt[eaggr_wkt.index("((") + 2 : eaggr_wkt.index("))")]
        coords = coords_section.split(",")
        # Append the first point to the end if not already closed
        if coords[0] != coords[-1]:
            coords.append(coords[0])
        fixed_coords = ", ".join(coords)
        return f"POLYGON (({fixed_coords}))"

def eaggrisea4t2geojson(eaggrisea4t):
    eaggr_dggs = Eaggr(Model.ISEA4T)
    eaggr_cell_shape = DggsShape(DggsCell(eaggrisea4t), DggsShapeLocation.ONE_FACE)._shape
    cell_to_shp = eaggr_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell_shape,ShapeStringFormat.WKT)
    cell_to_shp_fixed = fix_eaggr_wkt(cell_to_shp)
    cell_polygon = loads(cell_to_shp_fixed)

    resolution = len(eaggrisea4t)-2
    # Compute centroid
    cell_centroid = cell_polygon.centroid
    center_lat, center_lon = round(cell_centroid.y,7), round(cell_centroid.x,7)
    # Compute area using PyProj Geod
    geod = Geod(ellps="WGS84")
    cell_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])  # Area in square meters
    # Compute perimeter using PyProj Geod
    edge_len = abs(geod.geometry_area_perimeter(cell_polygon)[1])/3  # Perimeter in meters/ 3
    
    edge_len_str =  f'{round(edge_len,2)} m'
    cell_area_str=  f'{round(cell_area,2)} m2'

    if cell_area >= 1000_000:
        edge_len_str = f'{round(edge_len/1000,2)} km'
        cell_area_str = f'{round(cell_area/(10**6),2)} km2'
    

    if cell_polygon:
        coordinates = list(cell_polygon.exterior.coords)
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [coordinates]  
            },
            "properties": {
                 "eaggr_isea4t": eaggrisea4t,
                 "center_lat": center_lat,
                 "center_lon": center_lon,
                 "cell_area": cell_area_str,
                 "edge_len": edge_len_str,
                 "resolution": resolution,
                    }
        }

        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        return  feature_collection

def eaggrisea4t2geojson_cli():
    """
    Command-line interface for eaggrisea4t2geojson.
    """
    parser = argparse.ArgumentParser(description="Convert EaggrISEA4T code to GeoJSON")
    parser.add_argument("eaggrisea4t", help="Input EaggrISEA4T code, e.g., eaggrisea4t2geojson 131023133313201333311333")
    args = parser.parse_args()
    geojson_data = json.dumps(eaggrisea4t2geojson(args.eaggrisea4t))
    print(geojson_data)

def eaggrisea3hgeojson(eaggrisea3h):
    eaggr_dggs = Eaggr(Model.ISEA3H)
    eaggr_cell_shape = DggsShape(DggsCell(eaggrisea3h), DggsShapeLocation.ONE_FACE)._shape
    cell_to_shp = eaggr_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell_shape,ShapeStringFormat.WKT)
    
    if cell_to_shp:
        coordinates_part = cell_to_shp.replace("POLYGON ((", "").replace("))", "")
        coordinates = []
        for coord_pair in coordinates_part.split(","):
            lon, lat = map(float, coord_pair.strip().split())
            coordinates.append([lon, lat])

        # Ensure the polygon is closed (first and last point must be the same)
        if coordinates[0] != coordinates[-1]:
            coordinates.append(coordinates[0])

        # Step 3: Construct the GeoJSON feature
        feature = {
            "type": "Feature",
            "geometry": {
                "type": "Polygon",
                "coordinates": [coordinates]  # Directly use the coordinates list
            },
            "properties": {
                    }
        }

        # Step 4: Construct the FeatureCollection
        feature_collection = {
            "type": "FeatureCollection",
            "features": [feature]
        }
        return  feature_collection


def eaggrisea3h2geojson_cli():
    """
    Command-line interface for eaggriseah32geojson.
    """
    parser = argparse.ArgumentParser(description="Convert EeaggrISEA3H code to GeoJSON")
    parser.add_argument("eaggrisea3h", help="Input EeaggrISEA3H code, e.g., eaggrisea3h2geojson '07024,0'")
    args = parser.parse_args()
    geojson_data = json.dumps(eaggrisea4t2geojson(args.eaggrisea3h))
    print(geojson_data)