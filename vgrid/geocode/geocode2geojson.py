from vgrid.geocode import mgrs, maidenhead, geohash, georef, olc
from geopy.distance import geodesic
import geojson, os
import geopandas as gpd
from shapely.geometry import Polygon, Point
import vgrid.utils.mercantile as mercantile
import re

def vcode2geojson(vcode):
    """
    Converts a vcode (e.g., 'z8x11y14') to a GeoJSON Feature with a Polygon geometry
    representing the tile's bounds and includes the original vcode as a property.

    Args:
        vcode (str): The tile code in the format 'zXxYyZ'.

    Returns:
        dict: A GeoJSON Feature with a Polygon geometry and vcode as a property.
    """
    # Extract z, x, y from the vcode using regex
    match = re.match(r'z(\d+)x(\d+)y(\d+)', vcode)
    if not match:
        raise ValueError("Invalid vcode format. Expected format: 'zXxYyZ'")

    # Convert matched groups to integers
    z = int(match.group(1))
    x = int(match.group(2))
    y = int(match.group(3))

    # Get the bounds of the tile in (west, south, east, north)
    bbox = mercantile.bounds(x, y, z)
    if bbox:
        min_lat, min_lon = bbox.south, bbox.west  # Southwest corner
        max_lat, max_lon = bbox.north, bbox.east  # Northeast corner
        
        center_lat = (min_lat + max_lat) / 2
        center_lon = (min_lon + max_lon) / 2
        
        precision = z

        lat_len = geodesic((min_lat, min_lon), (max_lat, min_lon)).meters 
        lon_len = geodesic((min_lat, min_lon), (min_lat, max_lon)).meters  
        
        bbox_width =  f'{round(lon_len,1)} m'
        bbox_height =  f'{round(lat_len,1)} m'
        if lon_len >= 10000:
            bbox_width = f'{round(lon_len/1000,1)} km'
            bbox_height = f'{round(lat_len/1000,1)} km'
            
        # Create the coordinates of the polygon
        polygon_coords = [
           [min_lon, min_lat],  # Bottom-left corner
           [max_lon, min_lat],  # Bottom-right corner
           [max_lon, max_lat],  # Top-right corner
           [min_lon, max_lat],  # Top-left corner
           [min_lon, min_lat]   # Closing the polygon (same as the first point)
        ]

        geojson_feature = geojson.Feature(
        geometry=geojson.Polygon([polygon_coords]),
        properties={
            "vcode": vcode,
            "center_lat": center_lat,
            "center_lon": center_lon,
            "bbox_height": bbox_height,
            "bbox_width": bbox_width,
            "precision": precision  # Using the code length as precisi
            }
        )

    return geojson_feature

def olc2geojson(olc_code):
    # Decode the Open Location Code into a CodeArea object
    coord = olc.decode(olc_code)
    
    if coord:
        # Create the bounding box coordinates for the polygon
        min_lat, min_lon = coord.latitudeLo, coord.longitudeLo
        max_lat, max_lon = coord.latitudeHi, coord.longitudeHi

        center_lat, center_lon = coord.latitudeCenter, coord.longitudeCenter
        precision = coord.codeLength 
       
        lat_len = geodesic((min_lat, min_lon), (max_lat, min_lon)).meters 
        lon_len = geodesic((min_lat, min_lon), (min_lat, max_lon)).meters  

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
        
        geojson_feature = geojson.Feature(
            geometry=geojson.Polygon([polygon_coords]),
            properties={
                "olc": olc_code,  # Include the OLC as a property
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "precision": precision  # Using the code length as precision
                }
            )

        return geojson_feature


def maidenhead2geojson(maidenhead_code):
    # Decode the Open Location Code into a CodeArea object
    center_lat, center_lon, min_lat, min_lon, max_lat, max_lon, _ = maidenhead.maidenGrid(maidenhead_code)
    precision = int(len(maidenhead_code)/2)
    
    lat_len = geodesic((min_lat, min_lon), (max_lat, min_lon)).meters 
    lon_len = geodesic((min_lat, min_lon), (min_lat, max_lon)).meters  
    
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

        # Create the GeoJSON structure
        geojson_feature = geojson.Feature(
            geometry=geojson.Polygon([polygon_coords]),
            properties={
                "maidenhead": maidenhead_code,  # Include the OLC as a property
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "precision": precision  # Using the code length as precision
                }
            )

        return geojson_feature

def geohash2geojson(geohash_code):
    # Decode the Open Location Code into a CodeArea object
    bbox =  geohash.bbox(geohash_code)
    if bbox:
        min_lat, min_lon = bbox['s'], bbox['w']  # Southwest corner
        max_lat, max_lon = bbox['n'], bbox['e']  # Northeast corner
        
        center_lat = (min_lat + max_lat) / 2
        center_lon = (min_lon + max_lon) / 2
        
        precision =  len(geohash_code)

        lat_len = geodesic((min_lat, min_lon), (max_lat, min_lon)).meters 
        lon_len = geodesic((min_lat, min_lon), (min_lat, max_lon)).meters  
        
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

        # Create the GeoJSON structure
        geojson_feature = geojson.Feature(
            geometry=geojson.Polygon([polygon_coords]),
            properties={
                "geohsah": geohash_code,  
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "precision": precision  # Using the code length as precision
                }
            )

        return geojson_feature

def mgrs2geojson(mgrs_code,lat,lon):
    origin_lat, origin_lon, min_lat, min_lon, max_lat, max_lon,precision = mgrs.mgrscell(mgrs_code)

    lat_len = geodesic((min_lat, min_lon), (max_lat, min_lon)).meters 
    lon_len = geodesic((min_lat, min_lon), (min_lat, max_lon)).meters  
    
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

        # # Create the GeoJSON structure
        # geojson_feature = geojson.Feature(
        #     geometry=geojson.Polygon([polygon_coords]),
        #     properties={
        #         "mgrs": mgrs_code,  # Include the OLC as a property
        #         "origin_lat": origin_lat,
        #         "origin_lon": origin_lon,
        #         "bbox_height": bbox_height,
        #         "bbox_width": bbox_width,
        #         "precision": precision  # Using the code length as precision
        #         }
        #     )

        # return geojson_feature
        # Convert to a Shapely polygon
        # Convert to a Shapely polygon
        mgrs_polygon = Polygon(polygon_coords)

        # Create the original GeoJSON feature
        geojson_feature = geojson.Feature(
            geometry=geojson.Polygon([polygon_coords]),
            properties={
                "mgrs": mgrs_code,
                "origin_lat": origin_lat,
                "origin_lon": origin_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "precision": precision
            }
        )

        # Load the GZD GeoJSON file from the same folder
        gzd_geojson_path = os.path.join(os.path.dirname(__file__), 'gzd.geojson')
        with open(gzd_geojson_path) as f:
            gzd_geojson = geojson.load(f)

        # Convert the GZD GeoJSON to a GeoDataFrame
        gzd_gdf = gpd.GeoDataFrame.from_features(gzd_geojson['features'], crs="EPSG:4326")

        # Convert the MGRS polygon to a GeoSeries for intersection
        mgrs_gdf = gpd.GeoDataFrame(geometry=[mgrs_polygon], crs="EPSG:4326")

        # Perform the intersection
        intersection_gdf = gpd.overlay(mgrs_gdf, gzd_gdf, how='intersection')

        # Check if the intersection result is empty
        if not intersection_gdf.empty:
            # Convert lat/lon to a Shapely point
            point = Point(lon, lat)

            # Check if the point is inside any of the intersection polygons
            for intersection_polygon in intersection_gdf.geometry:
                if intersection_polygon.contains(point):
                    # Return the intersection as GeoJSON if the point is inside
                    intersection_geojson = geojson.Feature(
                        geometry=geojson.Polygon([list(intersection_polygon.exterior.coords)]),
                        properties={
                            "mgrs": mgrs_code,
                            "origin_lat": origin_lat,
                            "origin_lon": origin_lon,
                            "bbox_height": bbox_height,
                            "bbox_width": bbox_width,
                            "precision": precision
                        }
                    )
                    return intersection_geojson

        # If no intersection or point not contained, return the original MGRS GeoJSON
        return geojson_feature


def georef2geojson(georef_code):
    center_lat, center_lon, min_lat, min_lon, max_lat, max_lon,precision = georef.georefcell(georef_code)

    lat_len = geodesic((min_lat, min_lon), (max_lat, min_lon)).meters 
    lon_len = geodesic((min_lat, min_lon), (min_lat, max_lon)).meters
    
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

        # Create the GeoJSON structure
        geojson_feature = geojson.Feature(
            geometry=geojson.Polygon([polygon_coords]),
            properties={
                "georef": georef_code,  # Include the OLC as a property
                "center_lat": center_lat,
                "center_lon": center_lon,
                "bbox_height": bbox_height,
                "bbox_width": bbox_width,
                "precision": precision  # Using the code length as precision
                }
            )

        return geojson_feature
