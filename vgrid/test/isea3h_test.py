from vgrid.utils.eaggr.eaggr import Eaggr
from vgrid.utils.eaggr.enums.model import Model
from shapely.wkt import loads
from vgrid.utils.eaggr.shapes.lat_long_point import LatLongPoint
from vgrid.conversion.latlon2cell import *
from vgrid.conversion.cell2geojson import *

eaggr_dggs = Eaggr(Model.ISEA3H)

# max_accuracy = 10**14 # 10**14 is min cell_id length (7), eg. : 03000,0
# max_accuracy = 10**13 # 10**14 is min cell_id length (7), eg. : 03000,0
# max_accuracy = 10**-5 # 10**-5 is max cell_id length (26), eg: 0339-551800102,-1099037679
res_accuracy_dict = {
    18: 10**-5,   
    17: 10**-4,
    16: 10**-3,
    15: 10**-2,
    14: 10**-1,
    13: 10**0,
    12: 10**1,  
    11: 10**2,
    10: 10**3, 
    9: 10**4, 
    8: 5*10**4, 
    7: 10**5, 
    6: 10**7,
    5: 10**8,
    4: 10**9,
    3: 10**10,
    2: 10**11,
    1: 10**12,
    0: 10**14
}
res = 11 # res=11 is suitable for geocoding with avg_edg_len = 4.87m
accuracy = res_accuracy_dict.get(res)
latitude, longitude = 10.775275567242561, 106.70679737574993# 
# latitude, longitude = 11.88595373, 8.31375000
# latitude, longitude = 28.33265016,37.42848519
lat_long_point = LatLongPoint(latitude, longitude, accuracy)
isea3h_cell = eaggr_dggs.convert_point_to_dggs_cell(lat_long_point)
print(isea3h_cell.get_cell_id())
# print(len(isea3h_cell.get_cell_id()))
geojson_data = json.dumps(isea3h2geojson(isea3h_cell.get_cell_id()))
print(geojson_data)
point = eaggr_dggs.convert_dggs_cell_to_point(isea3h_cell)
print(point._latitude, point._longitude, point._accuracy)