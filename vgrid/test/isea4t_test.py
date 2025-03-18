from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.enums.model import Model
from vgrid.utils.eaggr.enums.dggs_shape_type import DggsShapeType
from vgrid.utils.eaggr.enums.dggs_shape_location import DggsShapeLocation
from vgrid.utils.eaggr.enums.shape_string_format import ShapeStringFormat
from vgrid.utils.eaggr.exceptions.eaggr_exception import EaggrException
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from vgrid.utils.eaggr.shapes.lat_long_linestring import LatLongLinestring
from vgrid.utils.eaggr.shapes.lat_long_polygon import LatLongPolygon
from vgrid.utils.eaggr.shapes.dggs_cell import DggsCell
from vgrid.utils.eaggr.shapes.dggs_linestring import DggsLinestring
from vgrid.utils.eaggr.shapes.dggs_polygon import DggsPolygon
from vgrid.utils.eaggr.enums.dggs_analysis_type import DggsAnalysisType
from vgrid.utils.eaggr.shapes.dggs_shape import DggsShape
from shapely.wkt import loads
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from shapely.geometry import Polygon
from shapely.ops import transform
from pyproj import Geod
from vgrid.conversion.latlon2dggs import *
from vgrid.conversion.dggs2geojson import *
from vgrid.utils.eaggr.shapes.lat_long_linestring import LatLongLinestring
import random
# eaggr_dggs = Eaggr(Model.ISEA4T)
# accuracy = 10**-10 # len = 41
# accuracy = 5*10**-10 # len = 40
# accuracy = 10**-9 # len = 39
# accuracy = 10**-8 # len = 38
# accuracy = 5*10**-8 # len = 37
# accuracy = 10**-7 # len = 36
# accuracy = 5*10**-7 # len = 35
# accuracy = 10**-6 # len = 34
# accuracy = 5*10**-6 # len = 33
# accuracy = 5*10**-5 # len = 32
# accuracy = 10**-4 # len = 31
# accuracy = 5*10**-4 # len = 30
# accuracy = 9*10**-4 # len = 29
# accuracy = 5*10**-3 # len = 28
# accuracy = 2*10**-2 # len = 27
# accuracy = 5*10**-2 # len = 26
# accuracy = 5*10**-1 # len = 25
# accuracy = 1 # len = 24
# accuracy = 10 # len = 23
# accuracy = 5*10 # len = 22
# accuracy = 10**2 # len = 21
# accuracy = 5*10**2 # len = 20
# accuracy = 10**3 # len = 19
# accuracy = 5*10**3 # len = 18
# accuracy = 5*10**4 # len = 17
# accuracy = 10**5 # len = 16
# accuracy = 5*10**5 # len = 15
# accuracy = 10**6 # len = 14
# accuracy = 5*10**6 # len = 13
# accuracy = 5*10**7 # len = 12
# accuracy = 10**8 # len = 11
# accuracy = 5*10**8 # len = 10
# accuracy = 10**9 # len = 9
# accuracy = 10**10 # len = 8
# accuracy = 5*10**10 # len = 7
# accuracy = 10**11 # len = 6
# accuracy = 5*10**11 # len = 5
# accuracy = 10**12 # len = 4
# accuracy = 5*10**12 # len = 3
accuracy = 5*10**13 # len = 2

eaggr_dggs = Eaggr(Model.ISEA4T)
latitude, longitude = 10.775275567242561, 106.70679737574993# 
lat_long_point = LatLongPoint(latitude, longitude,accuracy)

dggs_cell = eaggr_dggs.convert_point_to_dggs_cell(lat_long_point)
print(dggs_cell.get_cell_id())
# print(len(dggs_cell.get_cell_id()))


# wkt_string = 'LINESTRING(106.1 10.234, 107.2 14.322)'
# wkt_string = 'POLYGON((106.58897968 10.76341669, 106.69972187 10.76341669, 106.69972187 10.82837501, 106.58897968 10.82837501, 106.58897968 10.76341669))'
# string_format = ShapeStringFormat.WKT
# # Invalid input types
# shapes = eaggr_dggs.convert_shape_string_to_dggs_shapes(wkt_string, string_format, accuracy)
# for shape in shapes:
#     cells = shape.get_shape().get_outer_ring().get_cells()
#     # bounding_cell = eaggr_dggs.get_bounding_dggs_cell(cells)
#     # print(bounding_cell.get_cell_id())
#     for cell in cells:
#         # print(cell.get_cell_id())
#         print (len(cell.get_cell_id()))
    
# print(eaggr_dggs.convert_dggs_cells_to_shape_string([DggsCell('01'),DggsCell('02')], string_format))

# cell_to_shp =  eaggr_dggs.convert_dggs_cell_outline_to_shape_string(DggsCell('06'), ShapeStringFormat.WKT)
# print(cell_to_shp)

# eaggr_cell_shape = DggsShape(DggsCell('06'), DggsShapeLocation.ONE_FACE)._shape
# cell_to_shp = eaggr_dggs.convert_dggs_cell_outline_to_shape_string(eaggr_cell_shape, ShapeStringFormat.WKT)
# print(cell_to_shp)
# latitude, longitude = 10.775275567242561, 106.70679737574993# 
# latitude, longitude = -21.35158998051384, 179.9999999994089
# lat_long_point = LatLongPoint(latitude, longitude,resolution)
# # Initialise the DGGS model
# # Convert the lat/long point
# dggs_cell = eaggr_dggs.convert_point_to_dggs_cell(lat_long_point)
# print(dggs_cell.get_cell_id())
# # dggs_cell_2 = DggsCell('14010000')
# # Initialise the DGGS model
# # Convert the DGGS cell
# lat_long_point = eaggr_dggs.convert_dggs_cell_to_point(dggs_cell)
# print(lat_long_point._latitude, lat_long_point._longitude)


# dggs_shape_cell = DggsShape(DggsCell("05"), DggsShapeLocation.ONE_FACE)
# wkt = eaggr_dggs.convert_dggs_cell_outline_to_shape_string(DggsCell("14"), ShapeStringFormat.WKT)
# print(wkt)
# dggs_shape_cell = DggsShape(DggsCell("14220230"), DggsShapeLocation.ONE_FACE)
# # Compare the cell with another cell
# another_dggs_shape_cell = DggsShape(DggsCell("14220231"), DggsShapeLocation.ONE_FACE)


# analysisResult = eaggr_dggs.compare_dggs_shapes(dggs_shape_cell, another_dggs_shape_cell,
# DggsAnalysisType.TOUCHES)
# print(analysisResult)
# bounding_cell = eaggr_dggs.get_bounding_dggs_cell([DggsCell("14220010"),DggsCell("14220011"),DggsCell("14220012"),\
#                                                 DggsCell("14220200"),
#                                                 DggsCell("14220201"),
#                                                 DggsCell("14220203"),
#                                                 DggsCell("14220230"),
#                                                 DggsCell("14220231"),
#                                                 DggsCell("14220232"),
#                                                 DggsCell("14220233"),
#                                                 DggsCell("14221230"),
#                                                 DggsCell("14221232"),
#                                                 DggsCell("14221233")
#                                                 ])
# print(bounding_cell.get_cell_id())

# get_dggs_cell_siblings
# isea4t2geojson 13102313331320133331133
# isea4t_cell_id =  latlon2isea4t (latitude,longitude, 21) # latlon2isea4t <lat> <lon> <res> [0..38]

# isea4t_cell_id =  latlon2isea4t (latitude,longitude, 21)

# # isea4t_cell = DggsCell('13102313331323203322')
# # isea4t_cell_parent = eaggr_dggs.get_dggs_cell_parents(isea4t_cell)
# # for i in isea4t_cell_parent:
# #     print(i.get_cell_id())

# isea4t_cell = DggsCell('14220200')
# isea4t_cell_siblings = eaggr_dggs.get_dggs_cell_siblings(isea4t_cell)
# for i in isea4t_cell_siblings:
#     print(i.get_cell_id())

# siblings_of_parent = (eaggr_dggs.get_dggs_cell_siblings(isea4t_cell_parent[0]))
 
# for sibling_parent in siblings_of_parent:
#     print (sibling_parent.get_cell_id())
#     sibling_parent_children = eaggr_dggs.get_dggs_cell_children(sibling_parent)
#     for sibling_parent_child in sibling_parent_children:
#         print (sibling_parent_child.get_cell_id())

# isea4t_siblings = eaggr_dggs.get_dggs_cell_siblings(isea4t_cell_parent)
# for sibling in isea4t_siblings:
#     print (sibling.get_cell_id())


# def convert_point_to_dggs_cell_in_thread(latitude, longitude):
#     # Create the lat/long points
#     lat_long_point = LatLongPoint(latitude, longitude, 10)
#     # Initialise the DGGS model
#     dggs = Eaggr(Model.ISEA4T)
#     # Convert the lat/long point
#     dggs_cell = dggs.convert_point_to_dggs_cell(lat_long_point)
#     # Convert back to a lat/long point
#     # converted_point = dggs.convert_dggs_cell_to_point(dggs_cell)
#     print(dggs_cell._cell_id)
# convert_point_to_dggs_cell_in_thread(latitude, longitude)
# for latitude in range(-16, 16):
#     for longitude in range(-34, 34):
#         dggsRunners = []
#         dggsRunners.append(threading.Thread(convert_point_to_dggs_cell_in_thread(5 * latitude, 5 * longitude)))
#     # Start the threads
#     for runner in dggsRunners:
#         runner.start()
#     # Join the threads to wait for completion
#     for runner in dggsRunners:
#         runner.join()

# res = -10
# accuracy_sq_meters = 10**(res)
# # 10^14 = maximum resolution - base cells; 10^(-10)= minimum resolution
#  # Create the lat/long points
# # accuracy_sq_meters = 

# lat_long_point = LatLongPoint(10.775275567242561, 106.70679737574993, accuracy_sq_meters)
# # Initialise the DGGS model
# # Convert the first lat/long point
# dggs_cell = dggs.convert_point_to_dggs_cell(lat_long_point)
# cell_id_len = 23
# dggs_cell = DggsCell(dggs_cell._cell_id[:cell_id_len])
# cell_to_shp = dggs.convert_dggs_cell_outline_to_shape_string(dggs_cell,ShapeStringFormat.WKT)
# def fix_wkt(wkt):
#     # Extract the coordinate section
#     coords_section = wkt[wkt.index("((") + 2 : wkt.index("))")]
#     coords = coords_section.split(",")
#     # Append the first point to the end if not already closed
#     if coords[0] != coords[-1]:
#         coords.append(coords[0])
#     fixed_coords = ", ".join(coords)
#     return f"POLYGON (({fixed_coords}))"

# fixed_wkt_string = fix_wkt(cell_to_shp)

# # Convert fixed WKT to Shapely Polygon
# cell_polygon = loads(fixed_wkt_string)
# centroid = cell_polygon.centroid
# center_lon, center_lat = centroid.x, centroid.y
# print(center_lat,center_lon)
# geod = Geod(ellps="WGS84")
# geod_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
# geod_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])  # Perimeter in meters

# print(geod_area)
# print(geod_perimeter/3)


# # coordinates_part = cell_to_shp.replace("POLYGON ((", "").replace("))", "")
# # print(cell_to_shp)
# # triangle_coords = [
# #     (103.36956187681459, 6.112677457571387),
# #     (98.66047431289643, 13.781719950199278),
# #     (108.0, 12.373115939394218)
# # ]
# from shapely.wkt import loads


# # Create a Shapely Polygon
# triangle = Polygon(coordinates_part)
# # print(triangle)

# print(dggs_cell._cell_id)
# # dggs_cell = DggsCell('000')
# print(dggs_cell._cell_id)
# cell_to_point = dggs.convert_dggs_cell_to_point(dggs_cell)
# print(dggs_cell._cell_id)
# print(len(dggs_cell._cell_id))
# cell_to_shp = dggs.convert_dggs_cell_outline_to_shape_string(dggs_cell,ShapeStringFormat.WKT)
# print(cell_to_shp)
# coordinates_part = cell_to_shp.replace("POLYGON ((", "").replace("))", "")
# coordinates = []

# for coord_pair in coordinates_part.split(","):
#     lon, lat = map(float, coord_pair.strip().split())
#     coordinates.append([lon, lat])

# if coordinates[0] != coordinates[-1]:
#     coordinates.append(coordinates[0])
# print (coordinates)
# x_coords = [v[0] for v in coordinates]
# y_coords = [v[1] for v in coordinates]
# # 1. Centroid
# center_lon = sum(x_coords[:-1]) / 3
# center_lat = sum(y_coords[:-1]) / 3

# # 2. Area using the shoelace formula
# area = 0.5 * abs(
#     sum(x_coords[i] * y_coords[i + 1] for i in range(3)) -
#     sum(y_coords[i] * x_coords[i + 1] for i in range(3))
# )

# # 3. Average edge length
# def edge_length(x1, y1, x2, y2):
#     return math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
# edges = [
#     edge_length(x_coords[i], y_coords[i], x_coords[i + 1], y_coords[i + 1])
#     for i in range(3)
# ]
# average_edge_length = sum(edges) / 3

# print(center_lat)
# print(center_lon)
# print(area)
# print(average_edge_length)


# coordinates_part = cell_to_shp.replace("POLYGON ((", "").replace("))", "")
# # Step 2: Split and parse the coordinates
# # Convert "lon lat" to [lon, lat]
# coordinates = []
# for coord_pair in coordinates_part.split(","):
#     lon, lat = map(float, coord_pair.strip().split())
#     coordinates.append([lon, lat])

# # Ensure the polygon is closed (first and last point must be the same)
# if coordinates[0] != coordinates[-1]:
#     coordinates.append(coordinates[0])

# # Step 3: Construct the GeoJSON feature
# feature = {
#     "type": "Feature",
#     "geometry": {
#         "type": "Polygon",
#         "coordinates": [coordinates]  # Directly use the coordinates list
#     },
#     "properties": {
#             }
# }

# # Step 4: Construct the FeatureCollection
# feature_collection = {
#     "type": "FeatureCollection",
#     "features": [feature]
# }

# print(feature_collection)