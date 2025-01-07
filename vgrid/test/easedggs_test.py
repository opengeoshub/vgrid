from vgrid.utils.easedggs.dggs.grid_addressing import _grid_xy_to_grid_id, ease_polygon_to_grid_ids,geos_to_grid_ids, grid_ids_to_geos, grid_ids_to_ease,geo_polygon_to_grid_ids
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
import geopandas as gpd
from vgrid.utils.easedggs.constants import levels_specs
from vgrid.utils.easedggs.dggs.utils import gen_point_grid
# def get_cells_bbox(resolution, bbox):
#     min_lon, min_lat, max_lon, max_lat = bbox
#     bbox_coords = [
#             (min_lon, min_lat),
#             (min_lon, max_lat),
#             (max_lon, max_lat),
#             (max_lon, min_lat),
#             (min_lon, min_lat)
#         ]
#     bbox_polygon = Polygon(bbox_coords)
#     cells_bbox = geo_polygon_to_grid_ids(bbox_polygon.wkt, level=resolution, source_crs = geo_crs, target_crs = ease_crs, levels_specs = levels_specs, return_centroids = True, wkt_geom=True)
#     return cells_bbox

# resolution = 3
# bbox = [106.6990073571, 10.7628112647, 106.71767427, 10.7786496202]
# cells = get_cells_bbox(resolution, bbox)
# print(cells)
# print(grid_ids_to_geos(['L2.163763.03.21', 'L2.163763.03.11']))
# print(grid_ids_to_geos(['L3.165767.02.02.22']))
# print(grid_id_to_corner_coord('L2.163763.03.21',2))
from pyproj import Proj, transform
# Define the projection for UTM EPSG:6933
lon_min, lat_min, lon_max, lat_max = 106.6990073571, 10.7628112647, 106.71767427, 10.7786496202

print(gen_point_grid(lon_min, lat_min, lon_max, lat_max, 1, 1, target_crs = ease_crs))
    

# Coordinates in UTM (EPSG:6933)
# utm_x, utm_y = grid_id_to_corner_coord('L2.163763.03.11',2)
lon_min, lat_min, lon_max, lat_max = 106.6990073571, 10.7628112647, 106.71767427, 10.7786496202

# Create a transformer from WGS84 to UTM 6933
transformer = Transformer.from_crs("EPSG:4326", "EPSG:6933", always_xy=True)

# Transform the bounding box corners to UTM 6933
utm_min_x, utm_min_y = transformer.transform(lon_min, lat_min)
utm_max_x, utm_max_y = transformer.transform(lon_max, lat_max)

bound = easedggs_grid_bounds((utm_min_x, utm_min_y,utm_max_x, utm_max_y),ease_crs, 2)
# print (bound)


utm_min_x, utm_min_y, utm_max_x, utm_max_y = 10293204.420126759, 1363219.0218020855, 10299209.790266857, 1369224.3919421826

# Create a transformer from UTM 6933 to WGS84 (EPSG:4326)
transformer = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)

# Transform the bounding box corners back to WGS84
lon_min, lat_min = transformer.transform(utm_min_x, utm_min_y)
lon_max, lat_max = transformer.transform(utm_max_x, utm_max_y)

# print (lon_min, lat_min, lon_max, lat_max)


# transformer = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)

# # Perform the transformation
# lon, lat = transformer.transform(utm_x, utm_y)

# print(lon, lat)

# print(cells['result']['data'])
# for cell in cells['result']['data']:
#     print(cell)
#     geo = grid_ids_to_geos([cell])
#     print(geo)
#     # center_lon, center_lat = geo['result']['data'][0]
        
# def get_cells(res):
#     """
#     Generate a list of cell IDs based on the resolution, row, and column.
#     """
#     n_row = levels_specs[res]["n_row"]
#     n_col = levels_specs[res]["n_col"]

#     # Generate list of cell IDs
#     cell_ids = []

#     # Loop through all rows and columns at the specified resolution
#     for row in range(n_row):
#         for col in range(n_col):
#             # Generate base ID (e.g., L0.RRRCCC for res=0)
#             base_id = f"L{res}.{row:03d}{col:03d}"

#             # Add additional ".RC" for each higher resolution
#             cell_id = base_id
#             for i in range(1, res + 1):
#                 cell_id += f".{row:1d}{col:1d}"  # For res=1: L0.RRRCCC.RC, res=2: L0.RRRCCC.RC.RC, etc.

#             # Append the generated cell ID to the list
#             cell_ids.append(cell_id)

#     return cell_ids


# resolution = 0  # Change to desired resolution (0..6)
# cell_ids = get_cells(resolution)
# print(cell_ids)

# latitude, longitude = 10.775275567242561, 106.70679737574993# 
# resolution = 0 #[0..6]
# from pyproj import Transformer

# from pyproj import Geod

# # Define the geodetic calculator (WGS84 ellipsoid)
# geod = Geod(ellps="WGS84")
# # Center point (centroid)
# center_lon = 106.61825726141079
# center_lat = 10.651251665700356

# # # Dimensions in meters
# # width = 36032.22084058376
# # height = 36032.22084058376

# # # Calculate half dimensions
# half_width = width / 2
# half_height = height / 2

# Calculate the four corners of the bounding box
# Top-left (northwest)
# nw_lon, nw_lat, _ = geod.fwd(center_lon, center_lat, -90, half_width)
# nw_lon, nw_lat, _ = geod.fwd(nw_lon, nw_lat, 0, half_height)

# # Bottom-left (southwest)
# sw_lon, sw_lat, _ = geod.fwd(center_lon, center_lat, -90, half_width)
# sw_lon, sw_lat, _ = geod.fwd(sw_lon, sw_lat, 180, half_height)

# # Top-right (northeast)
# ne_lon, ne_lat, _ = geod.fwd(center_lon, center_lat, 90, half_width)
# ne_lon, ne_lat, _ = geod.fwd(ne_lon, ne_lat, 0, half_height)

# # Bottom-right (southeast)
# se_lon, se_lat, _ = geod.fwd(center_lon, center_lat, 90, half_width)
# se_lon, se_lat, _ = geod.fwd(se_lon, se_lat, 180, half_height)

# # Create GeoJSON
# geojson_features = {
#     "type": "FeatureCollection",
#     "features": [
#         {
#             "type": "Feature",
#             "geometry": {
#                 "type": "Polygon",
#                 "coordinates": [[
#                     [nw_lon, nw_lat],  # Top-left
#                     [ne_lon, ne_lat],  # Top-right
#                     [se_lon, se_lat],  # Bottom-right
#                     [sw_lon, sw_lat],  # Bottom-left
#                     [nw_lon, nw_lat],  # Close the polygon
#                 ]]
#             },
#             "properties": {
#                 "center": [center_lon, center_lat],
#                 "width_meters": width,
#                 "height_meters": height
#             }
#         }
#     ]
# }



# print(json.dumps(geojson_features))

# Define the bounding box in EPSG:6933
# bbox_6933 = {
#     'min_x': -17367530.445161372,
#     'min_y': -7314540.830638504,
#     'max_x': 17367530.445161372,
#     'max_y': 7314540.830638504
# }

# # Create a transformer from EPSG:6933 to EPSG:4326
# transformer = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)

# # Transform the corners of the bounding box
# min_lon, min_lat = transformer.transform(bbox_6933['min_x'], bbox_6933['min_y'])
# max_lon, max_lat = transformer.transform(bbox_6933['max_x'], bbox_6933['max_y'])

# # Output the converted bounding box
# bbox_4326 = {
#     'min_lon': min_lon,
#     'min_lat': min_lat,
#     'max_lon': max_lon,
#     'max_lat': max_lat
# }

# print("Bounding box in EPSG:4326:", bbox_4326)

# coords_ease = coords_lon_lat_to_coords_ease([(longitude,latitude)],source_crs=geo_crs,target_crs=ease_crs)
# for coord in coords_ease:
#     print(coord.x)
#     print(coord.y)

# coords_grid = coords_ease_to_coords_grid(coords_ease)
# for coord in coords_grid:
#     print(coord.x)
#     print(coord.y)
# print(latitude, longitude)

# easedggs_cell = geos_to_grid_ids([(longitude,latitude)],level = resolution)
# easedggs_cell_id = easedggs_cell['result']['data'][0]
# print (easedggs_cell_id)
# geo = grid_ids_to_geos([easedggs_cell_id])
# longitude, latitude = geo['result']['data'][0]  # Unpack tuple
# print (longitude, latitude)



# # corner_coord_1 = grid_id_to_corner_coord(easedggs_cell_id,False)
# # print(corner_coord_1)

# corners = {
#         "upper_left": grid_id_to_corner_coord(easedggs_cell_id, resolution),
#         "upper_right": grid_id_to_corner_coord(easedggs_cell_id, resolution, shift=True),
#         "lower_right": grid_id_to_corner_coord(easedggs_cell_id, resolution, shift=True),
#         "lower_left": grid_id_to_corner_coord(easedggs_cell_id, resolution),
#     }
# print(corners)
# x_values = [point[0] for point in corners.values()]
# y_values = [point[1] for point in corners.values()]

# # Compute bounding box
# min_x = min(x_values)
# max_x = max(x_values)
# min_y = min(y_values)
# max_y = max(y_values)

# # Bounding box
# easedggs_bounds = [min_x, min_y, max_x, max_y]
# transformer = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)

# # Extract the bounding box corners
# min_x, min_y, max_x, max_y = easedggs_bounds

# # Transform the corners
# min_lon, min_lat = transformer.transform(min_x, min_y)
# max_lon, max_lat = transformer.transform(max_x, max_y)

# # Construct the transformed bounding box
# bounding_box_epsg_4326 = [min_lon, min_lat, max_lon, max_lat]

# print("Bounding box in EPSG:4326:", bounding_box_epsg_4326)


# def convert_easedggs_bounds_to_wgs84(easedggs_bounds):
#     """
#     Convert EASE-Grid 2.0 (EPSG:6933) bounds to WGS84 (EPSG:4326).

#     Parameters
#     ----------
#     easedggs_bounds : list or tuple
#         Bounds in EPSG:6933 format (min_x, min_y, max_x, max_y).

#     Returns
#     -------
#     wgs84_bounds : tuple
#         Bounds in WGS84 (longitude, latitude) format (min_lon, min_lat, max_lon, max_lat).
#     """
#     # Define a transformer from EPSG:6933 to EPSG:4326
#     transformer = Transformer.from_crs("EPSG:6933", "EPSG:4326", always_xy=True)

#     # Unpack the bounds
#     min_x, min_y, max_x, max_y = easedggs_bounds

#     # Transform each corner of the bounds
#     min_lon, min_lat = transformer.transform(min_x, min_y)
#     max_lon, max_lat = transformer.transform(max_x, max_y)

#     # Return the bounds in WGS84
#     return min_lon, min_lat, max_lon, max_lat

# print (convert_easedggs_bounds_to_wgs84(corners))

# polygon_wkt = 'POLYGON((106.7050103098 10.7834344861, 106.7125825211 10.7834344861, 106.7125825211 10.7781930751, 106.7050103098 10.7781930751, 106.7050103098 10.7834344861))'
# # polygon_wkt = 'POLYGON((-180 -90, 180 -90, 180 90, -180 90, -180 -90))'


# bounds = (106.705, 10.783, 106.712, 10.778)  # (min_x, min_y, max_x, max_y)
# grid_bounds = easedggs_grid_bounds(bounds, geo_crs, 3)
# print("Ease Bounds:", grid_bounds)
# # Example usage
# # easedggs_bounds = [10295206.210173458, 1366221.7068721335, 10297208.000220157, 1368223.4969188329]
# wgs84_bounds = convert_easedggs_bounds_to_wgs84(grid_bounds)
# print("WGS84 Bounds:", wgs84_bounds)

# # print(ease_polygon_to_grid_ids(polygon_wkt,level=3))

# import geopandas as gpd

# def create_cell_polygon_from_point(levels_specs, level, point, origin=(0, 0)):
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


# # polygon = create_cell_polygon_from_point(levels_specs,resolution,coords_grid)
# # print(polygon)

# # def convert_polygon_to_wgs84(polygon, source_crs=6933):
# #     """
# #     Convert a polygon from the given CRS to WGS84 (EPSG:4326).
    
# #     Parameters:
# #     ----------
# #     polygon : shapely.geometry.Polygon
# #         The polygon geometry to be converted.
# #     source_crs : int, optional
# #         The source CRS EPSG code. Default is 6933 (EASE Grid v2).
    
# #     Returns:
# #     -------
# #     shapely.geometry.Polygon
# #         The reprojected polygon in WGS84 (EPSG:4326).
# #     """
# #     # Create a GeoDataFrame with the polygon
# #     gdf = gpd.GeoDataFrame({'geometry': [polygon]}, crs=f"EPSG:{source_crs}")
    
# #     # Convert the GeoDataFrame to WGS84 (EPSG:4326)
# #     gdf_wgs84 = gdf.to_crs(epsg=4326)
    
# #     # Extract the reprojected polygon
# #     return gdf_wgs84.geometry.iloc[0]


# # polygon_wgs84 = convert_polygon_to_wgs84(polygon)

# # Print the WGS84 polygon coordinates
# # print(f"Converted WGS84 polygon coordinates: {polygon_wgs84}")

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

# # geojson_featurecollection = polygon_to_geojson_featurecollection(polygon_wgs84)

# # # Print the GeoJSON FeatureCollection
# # print(geojson_featurecollection)

# # cell2latlong = grid_ids_to_ease([easedggs_cell_id])
# # easedggs_cell_id = 'L3.165767.02.02.22'
# # level = 6
# # cell2latlong = grid_ids_to_geos([easedggs_cell_id])
# # print(cell2latlong['result']['data'][0])

# # corner_coord = grid_id_to_corner_coord(easedggs_cell_id,0)
# # corners = {
# #         "upper_left": grid_id_to_corner_coord(easedggs_cell_id, level),
# #         "upper_right": grid_id_to_corner_coord(easedggs_cell_id, level, shift=True),
# #         "lower_right": grid_id_to_corner_coord(easedggs_cell_id, level, shift=True),
# #         "lower_left": grid_id_to_corner_coord(easedggs_cell_id, level),
# #     }

# # print(corners)
# # # # Convert cell_id to JSON Polygon
# # # def ease_to_wgs84(x, y):
# # #     return Transformer.from_crs(ease_crs, geo_crs, always_xy=True).transform
# # # upper_left = ease_to_wgs84(corner_coord[0], corner_coord[1])
# # # print (upper_left)

# def cell_id_to_geojson(cell_id, level, levels_specs):
#     """
#     Converts an EASE-DGGS cell ID to a GeoJSON FeatureCollection in WGS84 coordinates.

#     Parameters
#     ----------
#     cell_id : str
#         The cell ID (e.g., 'L3.165767.02.02.22').
#     level : int
#         The level of the cell.
#     levels_specs : dict
#         Specifications for the grid levels, containing x_length and y_length.

#     Returns
#     ----------
#     geojson : dict
#         A GeoJSON FeatureCollection representing the cell as a polygon in WGS84.
#     """
#     # Coordinates of the cell corners in some local CRS
#     corners = {
#         "upper_left": grid_id_to_corner_coord(cell_id, level),
#         "upper_right": grid_id_to_corner_coord(cell_id, level, shift=True),
#         "lower_right": grid_id_to_corner_coord(cell_id, level, shift=True),
#         "lower_left": grid_id_to_corner_coord(cell_id, level),
#     }

#     # Define projections using EPSG codes or PROJ strings
#     proj_ease = pyproj.CRS("EPSG:6933")  # Replace with the correct EPSG code or PROJ string for EASE-DGGS
#     proj_wgs84 = pyproj.CRS("EPSG:4326")  # WGS84 CRS

#     # Define transformer
#     transformer = Transformer.from_crs(proj_ease, proj_wgs84, always_xy=True)

#     # Convert corner coordinates to WGS84
#     coordinates_wgs84 = [
#         list(transformer.transform(corners["upper_left"][0], corners["upper_left"][1])),
#         list(transformer.transform(corners["upper_right"][0], corners["upper_right"][1])),
#         list(transformer.transform(corners["lower_right"][0], corners["lower_right"][1])),
#         list(transformer.transform(corners["lower_left"][0], corners["lower_left"][1])),
#         list(transformer.transform(corners["upper_left"][0], corners["upper_left"][1]))  # Close the polygon
#     ]

#     # Create GeoJSON FeatureCollection
#     geojson = {
#         "type": "FeatureCollection",
#         "features": [
#             {
#                 "type": "Feature",
#                 "geometry": {
#                     "type": "Polygon",
#                     "coordinates": [coordinates_wgs84]
#                 },
#                 "properties": {
#                     "cell_id": cell_id,
#                     "level": level
#                 }
#             }
#         ]
#     }

#     return geojson

# Example usage
# latitude, longitude = 10.775275567242561, 106.70679737574993# 
# resolution = 6 #[0..6]
# coords_ease = coords_lon_lat_to_coords_ease([(longitude,latitude)],source_crs=geo_crs,target_crs=ease_crs)
# ids = geos_to_grid_ids(coords_lon_lat=[(106.70679737574993, 10.775275567242561)])
# print(ids)
# # cell_id = 'L3.165767.02.02.22'
# cell_id = 'L0.165767'
# print(grid_ids_to_geos(grid_ids = ['L0.165767']))

# ease = grid_ids_to_ease([cell_id], cell_scale_factors  = cell_scale_factors, target_crs = ease_crs)
# geojson = ease.to_json()

# # Print GeoJSON
# print(geojson)

# level = 0
# polygon_json = cell_id_to_geojson(cell_id, level, levels_specs)

# # Print GeoJSON Polygon
# print(json.dumps(polygon_json))

# projected_crs = pyproj.CRS('EPSG:6933')  # Example: UTM Zone 33N
# wgs84_crs = pyproj.CRS('EPSG:4326')  # WGS84

# # Transformer to convert the coordinates
# transformer = pyproj.Transformer.from_crs(projected_crs, wgs84_crs, always_xy=True)

# # Input GeoJSON (example from your question)
# geojson = {
#     "type": "Point",
#     "coordinates": [10287199.049986664, 1351208.2815218912],
#     "bbox": [10287199.049986664, 1351208.2815218912, 10287199.049986664, 1351208.2815218912]
# }

# # Extract the coordinates and transform them
# x, y = geojson["coordinates"]
# longitude, latitude = transformer.transform(x, y)

# # Convert bounding box coordinates to WGS84 (also transformed)
# bbox = geojson["bbox"]
# x_min, y_min, x_max, y_max = bbox
# lon_min, lat_min = transformer.transform(x_min, y_min)
# lon_max, lat_max = transformer.transform(x_max, y_max)

# # Create the Point in WGS84
# point = Point(lon_min, lat_min)

# # Create the BBox Polygon in WGS84
# bbox_polygon = Polygon([(lon_min, lat_min), 
#                         (lon_max, lat_min), 
#                         (lon_max, lat_max), 
#                         (lon_min, lat_max), 
#                         (lon_min, lat_min)])

# # Create GeoJSON with Point and BBox Polygon
# multi_type_geojson = {
#     "type": "FeatureCollection",
#     "features": [
#         {
#             "type": "Feature",
#             "geometry": point.__geo_interface__,
#             "properties": {"type": "Point"}
#         },
#         {
#             "type": "Feature",
#             "geometry": bbox_polygon.__geo_interface__,
#             "properties": {"type": "Bounding Box Polygon"}
#         }
#     ]
# }

# # Print the result
# print(json.dumps(multi_type_geojson))