from vgrid.utils.easedggs.dggs.grid_addressing import _grid_xy_to_grid_id, geos_to_grid_ids, grid_ids_to_geos, grid_ids_to_ease,geo_polygon_to_grid_ids
from vgrid.utils.easedggs.grid_align import  easedggs_grid_bounds,grid_id_to_corner_coord
from vgrid.utils.easedggs.constants import grid_spec, levels_specs, ease_crs, geo_crs, cell_scale_factors
from vgrid.utils.easedggs.dggs.transforms import *
from shapely.wkt import loads
from shapely.geometry import Point, Polygon, mapping
import numpy as np
import json
import pyproj
from pyproj import Transformer
from vgrid.utils.easedggs.dggs.checks import check_gid_l0_index,check_coords_range


latitude, longitude = 10.775275567242561, 106.70679737574993# 
resolution = 6 #[0..6]
coords_ease = coords_lon_lat_to_coords_ease([(longitude,latitude)],source_crs=geo_crs,target_crs=ease_crs)
print(coords_ease[0].x)


# Extract the EASE and Geo min/max bounds
ease_min_x, ease_min_y = grid_spec['ease']['min_x'], grid_spec['ease']['min_y']
ease_max_x, ease_max_y = grid_spec['ease']['max_x'], grid_spec['ease']['max_y']

geo_min_x, geo_min_y = grid_spec['geo']['min_x'], grid_spec['geo']['min_y']
geo_max_x, geo_max_y = grid_spec['geo']['max_x'], grid_spec['geo']['max_y']

# Calculate scaling factors for EASE to Geo transformation
scale_x = (geo_max_x - geo_min_x) / (ease_max_x - ease_min_x)
scale_y = (geo_max_y - geo_min_y) / (ease_max_y - ease_min_y)

# Example EASE grid bounding box (min_x, min_y, max_x, max_y)
# ease_box = {
#     'min_x': -17367530.445161372,  # Example min_x of the bounding box
#     'min_y': -7314540.830638504,  # Example min_y of the bounding box
#     'max_x': 17367530.445161372,
#     'max_y': 7314540.830638504}

ease_box = {
    'min_x': coords_ease[0].x,  # Example min_x of the bounding box
    'min_y': coords_ease[0].y,  # Example min_y of the bounding box
    'max_x': coords_ease[0].x+14616000,
    'max_y': coords_ease[0].y+34704000
}
# Convert EASE bounding box to Geo (WGS84) coordinates
geo_min_x_converted = geo_min_x + (ease_box['min_x'] - ease_min_x) * scale_x
geo_min_y_converted = geo_min_y + (ease_box['min_y'] - ease_min_y) * scale_y
geo_max_x_converted = geo_min_x + (ease_box['max_x'] - ease_min_x) * scale_x
geo_max_y_converted = geo_min_y + (ease_box['max_y'] - ease_min_y) * scale_y

# Create a polygon using the converted Geo coordinates
polygon = Polygon([
    (geo_min_x_converted, geo_min_y_converted),  # Bottom-left
    (geo_max_x_converted, geo_min_y_converted),  # Bottom-right
    (geo_max_x_converted, geo_max_y_converted),  # Top-right
    (geo_min_x_converted, geo_max_y_converted),  # Top-left
    (geo_min_x_converted, geo_min_y_converted)   # Close the polygon
])

# Print the polygon
print(polygon)

# Convert the polygon to GeoJSON format (FeatureCollection)
geojson_feature_collection = {
    "type": "FeatureCollection",
    "features": [
        {
            "type": "Feature",
            "geometry": polygon.__geo_interface__,  # Get GeoJSON representation of the polygon
            "properties": {}  # You can add any properties here if needed
        }
    ]
}

# Convert the dictionary to a JSON string (GeoJSON format)
geojson_str = json.dumps(geojson_feature_collection)

# Print the GeoJSON output
print(geojson_str)



# # coords_grid = coords_ease_to_coords_grid(coords_ease)
# # print(coords_grid)
# # grid_ids = coords_grid.apply(lambda coord: _grid_xy_to_grid_id(coord, level = resolution))
# easedggs_cell = geos_to_grid_ids([(longitude,latitude)],level = resolution)
# easedggs_cell_id = easedggs_cell['result']['data'][0]
# print(easedggs_cell_id)
# geo = grid_ids_to_geos([easedggs_cell_id])
# print (geo)
# # polygon_wkt = 'POLYGON((106.7050103098 10.7834344861, 106.7125825211 10.7834344861, 106.7125825211 10.7781930751, 106.7050103098 10.7781930751, 106.7050103098 10.7834344861))'
# polygon_wkt = 'POLYGON((-180 -90, 180 -90, 180 90, -180 90, -180 -90))'

# print(check_coords_range([(longitude,latitude)]))
# print(geo_polygon_to_grid_ids(polygon_wkt,level=3))
# import geopandas as gpd

# def create_cell_polygon_from_point(levels_specs, level, point, origin=(0, 0)):
#     """
#     Creates a polygon for a grid cell based on the provided level's specifications and a Point geometry.
    
#     Parameters:
#     ----------
#     levels_specs : dict
#         The dictionary containing the grid specifications for each level.
#     level : int
#         The grid level for which the polygon needs to be generated.
#     point : shapely.geometry.Point
#         The Point geometry representing the center of the grid cell.
#     origin : tuple, optional
#         The origin (starting point) of the grid. Default is (0, 0).
        
#     Returns:
#     -------
#     Polygon
#         A Shapely Polygon representing the grid cell.
#     """
    
#     # Extract specifications for the given level
#     level_spec = levels_specs.get(level)
#     if not level_spec:
#         raise ValueError(f"Level {level} not found in levels_specs")
    
#     n_row = level_spec['n_row']
#     n_col = level_spec['n_col']
#     x_length = level_spec['x_length']
#     y_length = level_spec['y_length']
    
#     # Extract the x, y coordinates from the Point
#     x_center = point.x
#     y_center = point.y
    
#     # Calculate the bottom-left corner of the cell
#     x_min = x_center - x_length / 2
#     y_min = y_center - y_length / 2
    
#     # The top-right corner is just one step away in x and y directions
#     x_max = x_center + x_length / 2
#     y_max = y_center + y_length / 2
    
#     # Create a Polygon for the grid cell (using the corners)
#     cell_polygon = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max), (x_min, y_max)])
    
#     return cell_polygon

# polygon = create_cell_polygon_from_point(levels_specs,resolution,coords_grid)
# print(polygon)

# def convert_polygon_to_wgs84(polygon, source_crs=6933):
#     """
#     Convert a polygon from the given CRS to WGS84 (EPSG:4326).
    
#     Parameters:
#     ----------
#     polygon : shapely.geometry.Polygon
#         The polygon geometry to be converted.
#     source_crs : int, optional
#         The source CRS EPSG code. Default is 6933 (EASE Grid v2).
    
#     Returns:
#     -------
#     shapely.geometry.Polygon
#         The reprojected polygon in WGS84 (EPSG:4326).
#     """
#     # Create a GeoDataFrame with the polygon
#     gdf = gpd.GeoDataFrame({'geometry': [polygon]}, crs=f"EPSG:{source_crs}")
    
#     # Convert the GeoDataFrame to WGS84 (EPSG:4326)
#     gdf_wgs84 = gdf.to_crs(epsg=4326)
    
#     # Extract the reprojected polygon
#     return gdf_wgs84.geometry.iloc[0]


# polygon_wgs84 = convert_polygon_to_wgs84(polygon)

# # Print the WGS84 polygon coordinates
# print(f"Converted WGS84 polygon coordinates: {polygon_wgs84}")

# def polygon_to_geojson_featurecollection(polygon):
#     """
#     Convert a Shapely polygon to GeoJSON FeatureCollection format.
    
#     Parameters:
#     ----------
#     polygon : shapely.geometry.Polygon
#         The polygon to convert to GeoJSON.
        
#     Returns:
#     -------
#     str
#         The GeoJSON string representing a FeatureCollection containing the polygon.
#     """
#     # Extract the coordinates from the Polygon
#     geojson = {
#         "type": "FeatureCollection",
#         "features": [
#             {
#                 "type": "Feature",
#                 "geometry": {
#                     "type": "Polygon",
#                     "coordinates": [list(polygon.exterior.coords)]
#                 },
#                 "properties": {}
#             }
#         ]
#     }
    
#     # Convert the dictionary to a GeoJSON string
#     return json.dumps(geojson)

# geojson_featurecollection = polygon_to_geojson_featurecollection(polygon_wgs84)

# # Print the GeoJSON FeatureCollection
# print(geojson_featurecollection)

# cell2latlong = grid_ids_to_ease([easedggs_cell_id])
# easedggs_cell_id = 'L3.165767.02.02.22'
# level = 6
# cell2latlong = grid_ids_to_geos([easedggs_cell_id])
# print(cell2latlong['result']['data'][0])

# corner_coord = grid_id_to_corner_coord(easedggs_cell_id,0)
# corners = {
#         "upper_left": grid_id_to_corner_coord(easedggs_cell_id, level),
#         "upper_right": grid_id_to_corner_coord(easedggs_cell_id, level, shift=True),
#         "lower_right": grid_id_to_corner_coord(easedggs_cell_id, level, shift=True),
#         "lower_left": grid_id_to_corner_coord(easedggs_cell_id, level),
#     }

# print(corners)
# # # Convert cell_id to JSON Polygon
# # def ease_to_wgs84(x, y):
# #     return Transformer.from_crs(ease_crs, geo_crs, always_xy=True).transform
# # upper_left = ease_to_wgs84(corner_coord[0], corner_coord[1])
# # print (upper_left)

# # def cell_id_to_geojson(cell_id, level, levels_specs):
# #     """
# #     Converts an EASE-DGGS cell ID to a GeoJSON FeatureCollection in WGS84 coordinates.

# #     Parameters
# #     ----------
# #     cell_id : str
# #         The cell ID (e.g., 'L3.165767.02.02.22').
# #     level : int
# #         The level of the cell.
# #     levels_specs : dict
# #         Specifications for the grid levels, containing x_length and y_length.

# #     Returns
# #     ----------
# #     geojson : dict
# #         A GeoJSON FeatureCollection representing the cell as a polygon in WGS84.
# #     """
# #     corners = {
# #         "upper_left": grid_id_to_corner_coord(cell_id, level),
# #         "upper_right": grid_id_to_corner_coord(cell_id, level, shift=True),
# #         "lower_right": grid_id_to_corner_coord(cell_id, level, shift=True),
# #         "lower_left": grid_id_to_corner_coord(cell_id, level),
# #     }

# #     # Convert coordinates to WGS84 using pyproj (assuming EASE-DGGS is in some local CRS)
# #     proj_ease = pyproj.Proj(init="epsg:3395")  # Replace with correct CRS for EASE-DGGS
# #     proj_wgs84 = pyproj.Proj(init="epsg:4326")

# #     def ease_to_wgs84(x, y):
# #         return Transformer.from_crs(proj_ease, proj_wgs84, always_xy=True).transform

# #     # Convert corner coordinates to WGS84
# #     coordinates_wgs84 = [
# #         list(ease_to_wgs84(corners["upper_left"][0], corners["upper_left"][1])),
# #         list(ease_to_wgs84(corners["upper_right"][0], corners["upper_right"][1])),
# #         list(ease_to_wgs84(corners["lower_right"][0], corners["lower_right"][1])),
# #         list(ease_to_wgs84(corners["lower_left"][0], corners["lower_left"][1])),
# #         list(ease_to_wgs84(corners["upper_left"][0], corners["upper_left"][1]))  # Close the polygon
# #     ]

# #     # Create GeoJSON FeatureCollection in WGS84
# #     geojson = {
# #         "type": "FeatureCollection",
# #         "features": [
# #             {
# #                 "type": "Feature",
# #                 "geometry": {
# #                     "type": "Polygon",
# #                     "coordinates": [coordinates_wgs84]
# #                 },
# #                 "properties": {
# #                     "cell_id": cell_id,
# #                     "level": level
# #                 }
# #             }
# #         ]
# #     }

# #     return geojson

# # cell_id = 'L3.165767.02.02.22'
# # level = 3
# # polygon_json = cell_id_to_geojson(cell_id, level,levels_specs)

# # # Print JSON Polygon
# # print(json.dumps(polygon_json))
