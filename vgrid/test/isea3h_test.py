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
from vgrid.conversion.latlon2cell import *
from vgrid.conversion.cell2geojson import *
from vgrid.utils.eaggr.shapes.lat_long_linestring import LatLongLinestring

eaggr_dggs = Eaggr(Model.ISEA3H)
# eaggr_dggs = Eaggr(Model.ISEA4T)

res = 0
max_accuracy = 10**14 # 10**14 is min cell_id length (7), eg. : 03000,0
max_accuracy = 10**13 # 10**14 is min cell_id length (7), eg. : 03000,0
# max_accuracy = 10**-5 # 10**-5 is max cell_id length (26), eg: 0339-551800102,-1099037679

# latitude, longitude = 10.775275567242561, 106.70679737574993# 
# latitude, longitude = 11.88595373, 8.31375000
latitude, longitude = 28.33265016,37.42848519
lat_long_point = LatLongPoint(latitude, longitude,max_accuracy)

eaggr_cell_max_accuracy = eaggr_dggs.convert_point_to_dggs_cell(lat_long_point)
print(eaggr_cell_max_accuracy.get_cell_id())
cell_id_len = res+7
eaggr_cell = DggsCell(eaggr_cell_max_accuracy._cell_id[:cell_id_len])


# print(eaggr_cell.get_cell_id())

# center = DggsCell("07000,0")
# print(center.get_cell_id())
# siblings = eaggr_dggs.get_dggs_cell_children(center)
# for sibling in siblings:
#     print (sibling.get_cell_id())
# parent = eaggr_dggs.get_dggs_cell_parents(dggs_cell)
# print('parent: ')
# for p in parent: 
#     print(p.get_cell_id())

# print('siblings: ')
# siblings = eaggr_dggs.get_dggs_cell_siblings(dggs_cell)
# for sibling in siblings:
#     print(sibling.get_cell_id())

# print("Children:")
# children = eaggr_dggs.get_dggs_cell_children(dggs_cell)
# for child in children:
#     print(child.get_cell_id())
# # print(len(dggs_cell.get_cell_id()))
# base_cells = [
#     '00000,0', '01000,0', '02000,0', '03000,0', '04000,0', '05000,0', '06000,0', '07000,0', '08000,0', '09000,0',
#     '10000,0', '11000,0', '12000,0', '13000,0', '14000,0', '15000,0', '16000,0', '17000,0', '18000,0', '19000,0'
# ]
# print (base_cells[7])