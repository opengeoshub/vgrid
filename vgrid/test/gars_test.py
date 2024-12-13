
# from vgrid.utils.gars.garsgrid import GARSGrid
# from vgrid.geocode.geocode2geojson import *
# import json 

# latitude, longitude = 10.775275567242561, 106.70679737574993# GARS encoding
# gars_resolution = 1 # 1, 5, 15, 30 minutes
# gars_grid = GARSGrid.from_latlon(latitude, longitude, gars_resolution)
# gars_code = gars_grid.gars_id
# print(gars_code)

# data = gars2geojson(gars_code)
# print(data)
# output_file = f'gars_{gars_resolution}.geojson'
# with open(output_file, 'w') as f:
#     json.dump(data, f, indent=2)  
# print(f'GeoJSON written to {output_file}')

from vgrid.utils import s2, olc, geohash, georef, mgrs, tilecode, maidenhead, gars
import h3, json
from vgrid.geocode.geocode2geojson import *
from vgrid.geocode.latlon2geocode import *

latitude, longitude = 10.775275567242561, 106.70679737574993

from vgrid.utils import s2, olc, geohash, georef, mgrs, tilecode, maidenhead, gars
print('\nGARS:')
gars_resolution = 1 # [1, 5, 15, 30 minutes]
gars_grid = gars.garsgrid.GARSGrid.from_latlon(latitude, longitude, gars_resolution)
gars_code = gars_grid.gars_id
print(gars_code)

data = gars2geojson(gars_code)
output_file = f'gars_{gars_resolution}.geojson'
with open(output_file, 'w') as f:
    json.dump(data, f, indent=2)  
print(f'GeoJSON written to {output_file}')
