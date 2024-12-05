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
import math
from pyproj import Geod


# dggs = Eaggr(Model.ISEA3H)
dggs = Eaggr(Model.ISEA4T)

latitude, longitude = 10.775275567242561, 106.70679737574993# 

res = -10
accuracy_sq_meters = 10**(res)
# 10^14 = maximum resolution - base cells; 10^(-10)= minimum resolution
 # Create the lat/long points
# accuracy_sq_meters = 

lat_long_point = LatLongPoint(10.775275567242561, 106.70679737574993, accuracy_sq_meters)
# Initialise the DGGS model
# Convert the first lat/long point
dggs_cell = dggs.convert_point_to_dggs_cell(lat_long_point)
cell_id_len = 23
dggs_cell = DggsCell(dggs_cell._cell_id[:cell_id_len])
cell_to_shp = dggs.convert_dggs_cell_outline_to_shape_string(dggs_cell,ShapeStringFormat.WKT)
def fix_wkt(wkt):
    # Extract the coordinate section
    coords_section = wkt[wkt.index("((") + 2 : wkt.index("))")]
    coords = coords_section.split(",")
    # Append the first point to the end if not already closed
    if coords[0] != coords[-1]:
        coords.append(coords[0])
    fixed_coords = ", ".join(coords)
    return f"POLYGON (({fixed_coords}))"

fixed_wkt_string = fix_wkt(cell_to_shp)

# Convert fixed WKT to Shapely Polygon
cell_polygon = loads(fixed_wkt_string)
centroid = cell_polygon.centroid
center_lon, center_lat = centroid.x, centroid.y
print(center_lat,center_lon)
geod = Geod(ellps="WGS84")
geod_area = abs(geod.geometry_area_perimeter(cell_polygon)[0])
geod_perimeter = abs(geod.geometry_area_perimeter(cell_polygon)[1])  # Perimeter in meters

print(geod_area)
print(geod_perimeter/3)


# coordinates_part = cell_to_shp.replace("POLYGON ((", "").replace("))", "")
# print(cell_to_shp)
# triangle_coords = [
#     (103.36956187681459, 6.112677457571387),
#     (98.66047431289643, 13.781719950199278),
#     (108.0, 12.373115939394218)
# ]
from shapely.wkt import loads


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